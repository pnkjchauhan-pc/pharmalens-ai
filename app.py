import streamlit as st
import pandas as pd
import numpy as np
from pathlib import Path
from io import BytesIO

from rdkit import Chem
from rdkit.Chem import (
    Descriptors,
    Crippen,
    Lipinski,
    rdMolDescriptors,
    AllChem,
    Draw,
    DataStructs,
    inchi
)

import plotly.express as px
import plotly.graph_objects as go

from reportlab.lib import colors
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.enums import TA_CENTER
from reportlab.lib.units import inch
from reportlab.platypus import (
    SimpleDocTemplate,
    Paragraph,
    Spacer,
    Table,
    TableStyle,
    Image as ReportLabImage,
    PageBreak
)


# ============================================================
# PAGE CONFIGURATION
# ============================================================

st.set_page_config(
    page_title="PharmaLens AI",
    page_icon=None,
    layout="wide",
    initial_sidebar_state="expanded"
)


# ============================================================
# GLOBAL CSS
# ============================================================

st.markdown(
    """
    <style>
    .stApp {
        background:
            radial-gradient(circle at top right, rgba(40, 80, 150, 0.18), transparent 30%),
            radial-gradient(circle at bottom left, rgba(70, 40, 120, 0.12), transparent 35%),
            #07111f;
    }

    [data-testid="stSidebar"] {
        background: #091523;
        border-right: 1px solid rgba(255,255,255,0.08);
    }

    [data-testid="stSidebar"] * {
        color: #e8eef7;
    }

    h1, h2, h3 {
        letter-spacing: -0.4px;
    }

    .main-title {
        font-size: 2.7rem;
        font-weight: 750;
        margin-bottom: 0.1rem;
    }

    .subtitle {
        color: #9eacbf;
        font-size: 1rem;
        margin-bottom: 1.5rem;
    }

    .section-title {
        font-size: 1.45rem;
        font-weight: 700;
        margin-top: 1.2rem;
        margin-bottom: 0.8rem;
    }

    .small-muted {
        color: #8d9bad;
        font-size: 0.88rem;
    }

    div[data-testid="stMetric"] {
        background: rgba(16, 29, 47, 0.78);
        border: 1px solid rgba(255,255,255,0.07);
        border-radius: 14px;
        padding: 16px;
    }

    div[data-testid="stMetricLabel"] {
        color: #9eacbf;
    }

    div[data-testid="stMetricValue"] {
        color: #f3f7fb;
    }

    .stButton > button {
        border-radius: 10px;
        border: 1px solid rgba(255,255,255,0.1);
        font-weight: 600;
    }

    .stDownloadButton > button {
        border-radius: 10px;
        font-weight: 600;
    }

    .footer {
        text-align: center;
        color: #748399;
        padding: 30px 0 10px 0;
        font-size: 0.82rem;
    }
    </style>
    """,
    unsafe_allow_html=True
)


# ============================================================
# PATHS
# ============================================================

BASE_DIR = Path(__file__).resolve().parent
DATABASE_PATH = BASE_DIR / "data" / "compounds.csv"


# ============================================================
# DATABASE
# ============================================================

@st.cache_data
def load_database():
    if not DATABASE_PATH.exists():
        return pd.DataFrame(columns=["Name", "SMILES", "Category"])

    df = pd.read_csv(DATABASE_PATH)

    required_columns = ["Name", "SMILES", "Category"]

    for col in required_columns:
        if col not in df.columns:
            df[col] = ""

    df = df[required_columns].copy()

    return df


database = load_database()


# ============================================================
# MOLECULE HELPERS
# ============================================================

def mol_from_smiles(smiles):
    if not smiles:
        return None

    try:
        mol = Chem.MolFromSmiles(str(smiles))
        return mol
    except Exception:
        return None


def get_properties(mol):
    if mol is None:
        return {}

    return {
        "Molecular Weight": Descriptors.MolWt(mol),
        "LogP": Crippen.MolLogP(mol),
        "TPSA": rdMolDescriptors.CalcTPSA(mol),
        "HBD": Lipinski.NumHDonors(mol),
        "HBA": Lipinski.NumHAcceptors(mol),
        "Rotatable Bonds": Lipinski.NumRotatableBonds(mol),
        "Aromatic Rings": rdMolDescriptors.CalcNumAromaticRings(mol),
        "Heavy Atoms": Lipinski.HeavyAtomCount(mol),
        "Fraction Csp3": rdMolDescriptors.CalcFractionCSP3(mol),
        "Formula": rdMolDescriptors.CalcMolFormula(mol)
    }


def fingerprint(mol):
    if mol is None:
        return None

    generator = AllChem.GetMorganGenerator(
        radius=2,
        fpSize=2048
    )

    return generator.GetFingerprint(mol)


def fingerprint_array(mol):
    fp = fingerprint(mol)

    if fp is None:
        return None

    arr = np.zeros((2048,), dtype=int)
    DataStructs.ConvertToNumpyArray(fp, arr)

    return arr


def molecule_image(mol):
    if mol is None:
        return None

    return Draw.MolToImage(
        mol,
        size=(500, 400)
    )


def lipinski_status(mol):
    if mol is None:
        return "Invalid"

    mw = Descriptors.MolWt(mol)
    logp = Crippen.MolLogP(mol)
    hbd = Lipinski.NumHDonors(mol)
    hba = Lipinski.NumHAcceptors(mol)

    violations = 0

    if mw > 500:
        violations += 1

    if logp > 5:
        violations += 1

    if hbd > 5:
        violations += 1

    if hba > 10:
        violations += 1

    return "PASS" if violations <= 0 else f"FAIL ({violations} violation{'s' if violations != 1 else ''})"


def lipinski_table(mol):
    if mol is None:
        return pd.DataFrame()

    mw = Descriptors.MolWt(mol)
    logp = Crippen.MolLogP(mol)
    hbd = Lipinski.NumHDonors(mol)
    hba = Lipinski.NumHAcceptors(mol)

    return pd.DataFrame({
        "Rule": [
            "Molecular Weight ≤ 500",
            "LogP ≤ 5",
            "H-Bond Donors ≤ 5",
            "H-Bond Acceptors ≤ 10"
        ],
        "Value": [
            round(mw, 2),
            round(logp, 2),
            hbd,
            hba
        ],
        "Status": [
            "PASS" if mw <= 500 else "FAIL",
            "PASS" if logp <= 5 else "FAIL",
            "PASS" if hbd <= 5 else "FAIL",
            "PASS" if hba <= 10 else "FAIL"
        ]
    })


def molecular_inchi_key(mol):
    if mol is None:
        return "Unavailable"

    try:
        return inchi.MolToInchiKey(mol)
    except Exception:
        return "Unavailable"


# ============================================================
# BUILD DATABASE RECORDS
# ============================================================

records = []

for _, row in database.iterrows():
    mol = mol_from_smiles(row["SMILES"])

    if mol is None:
        continue

    props = get_properties(mol)

    records.append({
        "Name": row["Name"],
        "SMILES": row["SMILES"],
        "Category": row["Category"],
        "Molecular Weight": props["Molecular Weight"],
        "LogP": props["LogP"],
        "TPSA": props["TPSA"],
        "HBD": props["HBD"],
        "HBA": props["HBA"],
        "Rotatable Bonds": props["Rotatable Bonds"],
        "Aromatic Rings": props["Aromatic Rings"],
        "Heavy Atoms": props["Heavy Atoms"],
        "Fraction Csp3": props["Fraction Csp3"],
        "Formula": props["Formula"],
        "Lipinski": lipinski_status(mol)
    })


property_df = pd.DataFrame(records)


# ============================================================
# ADMET ENGINE
# ============================================================

def absorption_assessment(p):
    flags = []

    score = 3

    if p["Molecular Weight"] > 500:
        flags.append("High molecular weight may reduce oral absorption.")
        score -= 1

    if p["LogP"] > 5:
        flags.append("High lipophilicity may reduce balanced permeability.")
        score -= 1

    if p["TPSA"] > 140:
        flags.append("High TPSA may reduce passive membrane permeability.")
        score -= 1

    if p["HBD"] > 5:
        flags.append("High hydrogen-bond donor count may affect permeability.")
        score -= 1

    score = max(1, score)

    if score == 3:
        status = "Favorable profile"
    elif score == 2:
        status = "Intermediate profile"
    else:
        status = "Requires review"

    return status, flags, score


def distribution_assessment(p):
    flags = []

    score = 3

    if p["Molecular Weight"] > 600:
        flags.append("High molecular size may affect tissue distribution.")
        score -= 1

    if p["LogP"] > 5:
        flags.append("High lipophilicity may increase nonspecific distribution.")
        score -= 1

    if p["TPSA"] > 140:
        flags.append("High polarity may limit membrane passage.")
        score -= 1

    score = max(1, score)

    if score == 3:
        status = "Balanced profile"
    elif score == 2:
        status = "Intermediate profile"
    else:
        status = "Requires review"

    return status, flags, score


def metabolism_assessment(mol, p):
    flags = []

    score = 3

    if p["LogP"] > 5:
        flags.append("High lipophilicity may increase metabolic liability.")
        score -= 1

    if p["Rotatable Bonds"] > 10:
        flags.append("High flexibility may introduce additional metabolic sites.")
        score -= 1

    if p["Aromatic Rings"] > 4:
        flags.append("Multiple aromatic regions may require metabolic assessment.")
        score -= 1

    score = max(1, score)

    if score == 3:
        status = "Lower descriptor-based concern"
    elif score == 2:
        status = "Intermediate"
    else:
        status = "Requires review"

    return status, flags, score


def excretion_assessment(p):
    flags = []

    score = 3

    if p["Molecular Weight"] > 500:
        flags.append("Higher molecular size may affect renal elimination.")
        score -= 1

    if p["LogP"] > 5:
        flags.append("High lipophilicity may favor metabolic/biliary handling.")
        score -= 1

    if p["TPSA"] < 20:
        flags.append("Very low polarity may reduce direct renal clearance.")
        score -= 1

    score = max(1, score)

    if score == 3:
        status = "Balanced descriptor profile"
    elif score == 2:
        status = "Intermediate"
    else:
        status = "Requires review"

    return status, flags, score


def toxicity_assessment(mol, p):
    flags = []

    score = 4

    if p["LogP"] > 5:
        flags.append("High lipophilicity can be associated with nonspecific toxicity risk.")
        score -= 1

    if p["Molecular Weight"] > 600:
        flags.append("High molecular size warrants additional toxicity review.")
        score -= 1

    if p["Aromatic Rings"] > 4:
        flags.append("High aromatic ring count warrants structural review.")
        score -= 1

    if p["HBD"] > 5:
        flags.append("High hydrogen-bond donor count may affect exposure.")
        score -= 1

    score = max(1, score)

    if score >= 4:
        status = "No major descriptor flags"
    elif score >= 2:
        status = "Some descriptor flags"
    else:
        status = "Requires review"

    return status, flags, score


def run_admet(mol):
    p = get_properties(mol)

    absorption = absorption_assessment(p)
    distribution = distribution_assessment(p)
    metabolism = metabolism_assessment(mol, p)
    excretion = excretion_assessment(p)
    toxicity = toxicity_assessment(mol, p)

    return {
        "Absorption": absorption,
        "Distribution": distribution,
        "Metabolism": metabolism,
        "Excretion": excretion,
        "Toxicity": toxicity
    }


def overall_admet_score(results):
    scores = [v[2] for v in results.values()]

    total = sum(scores)
    maximum = len(scores) * 3

    return total, maximum


def admet_radar_values(p):
    absorption = 100

    if p["Molecular Weight"] > 500:
        absorption -= 20

    if p["LogP"] > 5:
        absorption -= 20

    if p["TPSA"] > 140:
        absorption -= 20

    distribution = 100

    if p["Molecular Weight"] > 600:
        distribution -= 20

    if p["LogP"] > 5:
        distribution -= 20

    metabolism = 100

    if p["LogP"] > 5:
        metabolism -= 20

    if p["Rotatable Bonds"] > 10:
        metabolism -= 15

    if p["Aromatic Rings"] > 4:
        metabolism -= 15

    excretion = 100

    if p["Molecular Weight"] > 500:
        excretion -= 20

    if p["LogP"] > 5:
        excretion -= 20

    toxicity = 100

    if p["LogP"] > 5:
        toxicity -= 20

    if p["Molecular Weight"] > 600:
        toxicity -= 15

    if p["Aromatic Rings"] > 4:
        toxicity -= 15

    return {
        "Absorption": max(0, absorption),
        "Distribution": max(0, distribution),
        "Metabolism": max(0, metabolism),
        "Excretion": max(0, excretion),
        "Toxicity": max(0, toxicity)
    }


# ============================================================
# BATCH ANALYZER
# ============================================================

def analyze_batch_smiles(smiles_text):
    rows = []
    invalid = []

    lines = [
        line.strip()
        for line in smiles_text.splitlines()
        if line.strip()
    ]

    for i, smiles in enumerate(lines, start=1):
        mol = mol_from_smiles(smiles)

        if mol is None:
            invalid.append({
                "Line": i,
                "SMILES": smiles
            })
            continue

        p = get_properties(mol)

        rows.append({
            "Compound": f"Compound {i}",
            "SMILES": smiles,
            "Formula": p["Formula"],
            "MW": round(p["Molecular Weight"], 2),
            "LogP": round(p["LogP"], 2),
            "TPSA": round(p["TPSA"], 2),
            "HBD": p["HBD"],
            "HBA": p["HBA"],
            "Rotatable Bonds": p["Rotatable Bonds"],
            "Aromatic Rings": p["Aromatic Rings"],
            "Heavy Atoms": p["Heavy Atoms"],
            "Fraction Csp3": round(p["Fraction Csp3"], 3),
            "Lipinski": lipinski_status(mol)
        })

    return pd.DataFrame(rows), pd.DataFrame(invalid)


# ============================================================
# PDF REPORT GENERATOR
# ============================================================

def create_pdf_report(
    compound_name,
    category,
    smiles,
    mol,
    admet_results,
    admet_values
):
    buffer = BytesIO()

    doc = SimpleDocTemplate(
        buffer,
        pagesize=A4,
        rightMargin=40,
        leftMargin=40,
        topMargin=40,
        bottomMargin=40
    )

    styles = getSampleStyleSheet()

    title_style = ParagraphStyle(
        "CustomTitle",
        parent=styles["Title"],
        fontSize=22,
        leading=26,
        alignment=TA_CENTER,
        textColor=colors.HexColor("#152238")
    )

    heading_style = ParagraphStyle(
        "Heading",
        parent=styles["Heading2"],
        fontSize=14,
        leading=18,
        textColor=colors.HexColor("#152238"),
        spaceBefore=12,
        spaceAfter=8
    )

    normal_style = ParagraphStyle(
        "NormalCustom",
        parent=styles["BodyText"],
        fontSize=9,
        leading=13
    )

    story = []

    story.append(
        Paragraph(
            "PharmaLens AI",
            title_style
        )
    )

    story.append(
        Paragraph(
            "Professional Molecular Analysis Report",
            styles["Heading3"]
        )
    )

    story.append(Spacer(1, 12))

    story.append(
        Paragraph(
            f"<b>Compound:</b> {compound_name}<br/>"
            f"<b>Category:</b> {category}",
            normal_style
        )
    )

    story.append(Spacer(1, 15))

    image = molecule_image(mol)

    if image is not None:
        image_buffer = BytesIO()
        image.save(image_buffer, format="PNG")
        image_buffer.seek(0)

        story.append(
            ReportLabImage(
                image_buffer,
                width=3.8 * inch,
                height=3.0 * inch
            )
        )

    story.append(Spacer(1, 10))

    story.append(
        Paragraph(
            "Molecular Identifiers",
            heading_style
        )
    )

    p = get_properties(mol)

    identifier_data = [
        ["Parameter", "Value"],
        ["Compound", compound_name],
        ["Category", category],
        ["SMILES", smiles],
        ["Formula", p["Formula"]],
        ["InChI Key", molecular_inchi_key(mol)]
    ]

    identifier_table = Table(
        identifier_data,
        colWidths=[1.5 * inch, 5.4 * inch]
    )

    identifier_table.setStyle(
        TableStyle([
            ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#152238")),
            ("TEXTCOLOR", (0, 0), (-1, 0), colors.white),
            ("GRID", (0, 0), (-1, -1), 0.4, colors.grey),
            ("VALIGN", (0, 0), (-1, -1), "TOP"),
            ("FONTSIZE", (0, 0), (-1, -1), 8),
            ("BOTTOMPADDING", (0, 0), (-1, -1), 6),
            ("TOPPADDING", (0, 0), (-1, -1), 6)
        ])
    )

    story.append(identifier_table)

    story.append(
        Paragraph(
            "Molecular Descriptors",
            heading_style
        )
    )

    descriptor_data = [
        ["Descriptor", "Value"],
        ["Molecular Weight", f"{p['Molecular Weight']:.2f} Da"],
        ["LogP", f"{p['LogP']:.2f}"],
        ["TPSA", f"{p['TPSA']:.2f} Å²"],
        ["H-Bond Donors", str(p["HBD"])],
        ["H-Bond Acceptors", str(p["HBA"])],
        ["Rotatable Bonds", str(p["Rotatable Bonds"])],
        ["Aromatic Rings", str(p["Aromatic Rings"])],
        ["Heavy Atoms", str(p["Heavy Atoms"])],
        ["Fraction Csp3", f"{p['Fraction Csp3']:.3f}"]
    ]

    descriptor_table = Table(
        descriptor_data,
        colWidths=[3.0 * inch, 3.9 * inch]
    )

    descriptor_table.setStyle(
        TableStyle([
            ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#152238")),
            ("TEXTCOLOR", (0, 0), (-1, 0), colors.white),
            ("GRID", (0, 0), (-1, -1), 0.4, colors.grey),
            ("FONTSIZE", (0, 0), (-1, -1), 8),
            ("BOTTOMPADDING", (0, 0), (-1, -1), 5),
            ("TOPPADDING", (0, 0), (-1, -1), 5)
        ])
    )

    story.append(descriptor_table)

    story.append(
        Paragraph(
            "Lipinski Rule of Five",
            heading_style
        )
    )

    lipinski_df = lipinski_table(mol)

    lipinski_data = [
        ["Rule", "Value", "Status"]
    ]

    for _, row in lipinski_df.iterrows():
        lipinski_data.append([
            row["Rule"],
            str(row["Value"]),
            row["Status"]
        ])

    lipinski_table_pdf = Table(
        lipinski_data,
        colWidths=[3.7 * inch, 1.1 * inch, 1.4 * inch]
    )

    lipinski_table_pdf.setStyle(
        TableStyle([
            ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#152238")),
            ("TEXTCOLOR", (0, 0), (-1, 0), colors.white),
            ("GRID", (0, 0), (-1, -1), 0.4, colors.grey),
            ("FONTSIZE", (0, 0), (-1, -1), 8),
            ("ALIGN", (1, 1), (-1, -1), "CENTER")
        ])
    )

    story.append(lipinski_table_pdf)

    story.append(PageBreak())

    story.append(
        Paragraph(
            "ADMET Intelligence",
            heading_style
        )
    )

    admet_data = [
        ["Domain", "Assessment", "Profile Score"]
    ]

    for domain, values in admet_results.items():
        admet_data.append([
            domain,
            values[0],
            str(values[2])
        ])

    admet_table = Table(
        admet_data,
        colWidths=[1.5 * inch, 4.1 * inch, 1.3 * inch]
    )

    admet_table.setStyle(
        TableStyle([
            ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#152238")),
            ("TEXTCOLOR", (0, 0), (-1, 0), colors.white),
            ("GRID", (0, 0), (-1, -1), 0.4, colors.grey),
            ("FONTSIZE", (0, 0), (-1, -1), 8),
            ("ALIGN", (2, 1), (2, -1), "CENTER")
        ])
    )

    story.append(admet_table)

    total, maximum = overall_admet_score(admet_results)

    story.append(Spacer(1, 12))

    story.append(
        Paragraph(
            f"<b>Descriptor Screening Score:</b> {total}/{maximum}",
            normal_style
        )
    )

    story.append(
        Paragraph(
            "ADMET Profile Values",
            heading_style
        )
    )

    radar_data = [
        ["Domain", "Score"]
    ]

    for domain, value in admet_values.items():
        radar_data.append([
            domain,
            f"{value:.1f}%"
        ])

    radar_table = Table(
        radar_data,
        colWidths=[4.5 * inch, 2.4 * inch]
    )

    radar_table.setStyle(
        TableStyle([
            ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#152238")),
            ("TEXTCOLOR", (0, 0), (-1, 0), colors.white),
            ("GRID", (0, 0), (-1, -1), 0.4, colors.grey),
            ("FONTSIZE", (0, 0), (-1, -1), 8)
        ])
    )

    story.append(radar_table)

    story.append(
        Paragraph(
            "ADMET Flags",
            heading_style
        )
    )

    any_flags = False

    for domain, values in admet_results.items():
        flags = values[1]

        if flags:
            any_flags = True

            story.append(
                Paragraph(
                    f"<b>{domain}</b>",
                    normal_style
                )
            )

            for flag in flags:
                story.append(
                    Paragraph(
                        f"• {flag}",
                        normal_style
                    )
                )

            story.append(Spacer(1, 5))

    if not any_flags:
        story.append(
            Paragraph(
                "No descriptor-based flags were identified.",
                normal_style
            )
        )

    story.append(Spacer(1, 15))

    story.append(
        Paragraph(
            "<b>Scientific Disclaimer:</b> This report is intended for "
            "research and educational screening only. Descriptor-based "
            "ADMET assessments are not experimental measurements, clinical "
            "predictions, or regulatory conclusions.",
            normal_style
        )
    )

    story.append(Spacer(1, 20))

    story.append(
        Paragraph(
            "Generated by PharmaLens AI",
            ParagraphStyle(
                "Footer",
                parent=normal_style,
                alignment=TA_CENTER,
                textColor=colors.grey
            )
        )
    )

    doc.build(story)

    buffer.seek(0)

    return buffer.getvalue()


# ============================================================
# SIDEBAR NAVIGATION
# ============================================================

st.sidebar.title("PharmaLens AI")
st.sidebar.caption("Drug Discovery & Molecular Intelligence")

page = st.sidebar.radio(
    "Navigation",
    [
        "Dashboard",
        "Drug Discovery Analytics",
        "Molecular Analysis",
        "Molecular Insights",
        "Lead Discovery",
        "ADMET Intelligence",
        "ADMET Comparison",
        "Compound Database",
        "Compound Comparison",
        "Chemical Space",
        "Custom Molecule",
        "Similarity Search",
        "Similarity Explorer",
        "Drug-Likeness",
        "Batch SMILES Analyzer",
        "Compound Report Generator"
    ]
)

st.sidebar.divider()

st.sidebar.caption(
    f"Database: {len(property_df)} valid compounds"
)

st.sidebar.caption(
    f"Categories: {property_df['Category'].nunique() if not property_df.empty else 0}"
)


# ============================================================
# DASHBOARD
# ============================================================

if page == "Dashboard":

    st.markdown(
        '<div class="main-title">PharmaLens AI</div>',
        unsafe_allow_html=True
    )

    st.markdown(
        '<div class="subtitle">Molecular intelligence platform for drug discovery research</div>',
        unsafe_allow_html=True
    )

    total_compounds = len(property_df)
    total_categories = property_df["Category"].nunique() if not property_df.empty else 0
    avg_mw = property_df["Molecular Weight"].mean() if not property_df.empty else 0

    lipinski_pass = (
        property_df["Lipinski"].astype(str).str.startswith("PASS").mean() * 100
        if not property_df.empty
        else 0
    )

    c1, c2, c3, c4 = st.columns(4)

    c1.metric("Total Compounds", total_compounds)
    c2.metric("Compound Categories", total_categories)
    c3.metric("Average Molecular Weight", f"{avg_mw:.1f} Da")
    c4.metric("Lipinski Pass Rate", f"{lipinski_pass:.1f}%")

    st.markdown("### Compound Landscape")

    if not property_df.empty:

        col1, col2 = st.columns(2)

        with col1:
            category_counts = (
                property_df["Category"]
                .value_counts()
                .reset_index()
            )

            category_counts.columns = ["Category", "Count"]

            fig = px.bar(
                category_counts,
                x="Category",
                y="Count",
                title="Compounds by Category"
            )

            fig.update_layout(
                template="plotly_dark",
                paper_bgcolor="rgba(0,0,0,0)",
                plot_bgcolor="rgba(0,0,0,0)"
            )

            st.plotly_chart(
                fig,
                use_container_width=True
            )

        with col2:
            fig = px.scatter(
                property_df,
                x="Molecular Weight",
                y="LogP",
                color="Category",
                hover_name="Name",
                hover_data=[
                    "TPSA",
                    "HBD",
                    "HBA"
                ],
                title="Molecular Weight vs LogP"
            )

            fig.update_layout(
                template="plotly_dark",
                paper_bgcolor="rgba(0,0,0,0)",
                plot_bgcolor="rgba(0,0,0,0)"
            )

            st.plotly_chart(
                fig,
                use_container_width=True
            )


# ============================================================
# DRUG DISCOVERY ANALYTICS
# ============================================================

elif page == "Drug Discovery Analytics":

    st.markdown(
        '<div class="main-title">Drug Discovery Analytics</div>',
        unsafe_allow_html=True
    )

    st.markdown(
        '<div class="subtitle">Interactive analysis of molecular properties, chemical diversity and drug-likeness</div>',
        unsafe_allow_html=True
    )

    if property_df.empty:
        st.warning("No valid compound data available.")
        st.stop()

    # --------------------------------------------------------
    # FILTERS
    # --------------------------------------------------------

    st.markdown("### Analytics Filters")

    f1, f2, f3 = st.columns(3)

    with f1:
        categories = st.multiselect(
            "Compound Categories",
            sorted(property_df["Category"].dropna().unique()),
            default=sorted(property_df["Category"].dropna().unique())
        )

    with f2:
        lipinski_filter = st.selectbox(
            "Lipinski Status",
            ["All", "PASS", "FAIL"]
        )

    with f3:
        top_n = st.slider(
            "Top compounds to display",
            min_value=5,
            max_value=min(30, len(property_df)),
            value=min(10, len(property_df))
        )

    filtered = property_df[
        property_df["Category"].isin(categories)
    ].copy()

    if lipinski_filter == "PASS":
        filtered = filtered[
            filtered["Lipinski"].astype(str).str.startswith("PASS")
        ]

    elif lipinski_filter == "FAIL":
        filtered = filtered[
            ~filtered["Lipinski"].astype(str).str.startswith("PASS")
        ]

    st.markdown("### Dataset Overview")

    a1, a2, a3, a4, a5 = st.columns(5)

    a1.metric("Filtered Compounds", len(filtered))

    if not filtered.empty:
        a2.metric(
            "Average MW",
            f"{filtered['Molecular Weight'].mean():.1f} Da"
        )

        a3.metric(
            "Average LogP",
            f"{filtered['LogP'].mean():.2f}"
        )

        a4.metric(
            "Average TPSA",
            f"{filtered['TPSA'].mean():.1f}"
        )

        a5.metric(
            "Average HBA",
            f"{filtered['HBA'].mean():.1f}"
        )

    if filtered.empty:
        st.info("No compounds match the selected filters.")
        st.stop()

    # --------------------------------------------------------
    # PROPERTY DISTRIBUTIONS
    # --------------------------------------------------------

    st.markdown("### Molecular Property Distributions")

    p1, p2 = st.columns(2)

    with p1:

        fig = px.histogram(
            filtered,
            x="Molecular Weight",
            nbins=15,
            title="Molecular Weight Distribution"
        )

        fig.update_layout(
            template="plotly_dark",
            paper_bgcolor="rgba(0,0,0,0)",
            plot_bgcolor="rgba(0,0,0,0)"
        )

        st.plotly_chart(
            fig,
            use_container_width=True
        )

    with p2:

        fig = px.histogram(
            filtered,
            x="LogP",
            nbins=15,
            title="LogP Distribution"
        )

        fig.update_layout(
            template="plotly_dark",
            paper_bgcolor="rgba(0,0,0,0)",
            plot_bgcolor="rgba(0,0,0,0)"
        )

        st.plotly_chart(
            fig,
            use_container_width=True
        )

    p3, p4 = st.columns(2)

    with p3:

        fig = px.histogram(
            filtered,
            x="TPSA",
            nbins=15,
            title="TPSA Distribution"
        )

        fig.update_layout(
            template="plotly_dark",
            paper_bgcolor="rgba(0,0,0,0)",
            plot_bgcolor="rgba(0,0,0,0)"
        )

        st.plotly_chart(
            fig,
            use_container_width=True
        )

    with p4:

        fig = px.histogram(
            filtered,
            x="Fraction Csp3",
            nbins=15,
            title="Fraction Csp3 Distribution"
        )

        fig.update_layout(
            template="plotly_dark",
            paper_bgcolor="rgba(0,0,0,0)",
            plot_bgcolor="rgba(0,0,0,0)"
        )

        st.plotly_chart(
            fig,
            use_container_width=True
        )

    # --------------------------------------------------------
    # PROPERTY RELATIONSHIPS
    # --------------------------------------------------------

    st.markdown("### Molecular Property Relationships")

    relationship = st.selectbox(
        "Select relationship",
        [
            "Molecular Weight vs LogP",
            "Molecular Weight vs TPSA",
            "LogP vs TPSA",
            "TPSA vs HBD",
            "Molecular Weight vs Rotatable Bonds"
        ]
    )

    relationships = {
        "Molecular Weight vs LogP":
            ("Molecular Weight", "LogP"),

        "Molecular Weight vs TPSA":
            ("Molecular Weight", "TPSA"),

        "LogP vs TPSA":
            ("LogP", "TPSA"),

        "TPSA vs HBD":
            ("TPSA", "HBD"),

        "Molecular Weight vs Rotatable Bonds":
            ("Molecular Weight", "Rotatable Bonds")
    }

    x_col, y_col = relationships[relationship]

    fig = px.scatter(
        filtered,
        x=x_col,
        y=y_col,
        color="Category",
        hover_name="Name",
        hover_data=[
            "Molecular Weight",
            "LogP",
            "TPSA",
            "HBD",
            "HBA"
        ],
        title=relationship
    )

    fig.update_layout(
        template="plotly_dark",
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)"
    )

    st.plotly_chart(
        fig,
        use_container_width=True
    )

    # --------------------------------------------------------
    # CATEGORY ANALYTICS
    # --------------------------------------------------------

    st.markdown("### Category Analytics")

    category_summary = (
        filtered
        .groupby("Category")
        .agg(
            Compounds=("Name", "count"),
            Avg_MW=("Molecular Weight", "mean"),
            Avg_LogP=("LogP", "mean"),
            Avg_TPSA=("TPSA", "mean"),
            Avg_HBD=("HBD", "mean"),
            Avg_HBA=("HBA", "mean")
        )
        .reset_index()
    )

    category_summary["Avg_MW"] = category_summary["Avg_MW"].round(2)
    category_summary["Avg_LogP"] = category_summary["Avg_LogP"].round(2)
    category_summary["Avg_TPSA"] = category_summary["Avg_TPSA"].round(2)
    category_summary["Avg_HBD"] = category_summary["Avg_HBD"].round(2)
    category_summary["Avg_HBA"] = category_summary["Avg_HBA"].round(2)

    st.dataframe(
        category_summary,
        use_container_width=True,
        hide_index=True
    )

    # --------------------------------------------------------
    # TOP COMPOUNDS
    # --------------------------------------------------------

    st.markdown("### Molecular Property Leaders")

    leader_property = st.selectbox(
        "Rank compounds by",
        [
            "Molecular Weight",
            "LogP",
            "TPSA",
            "HBD",
            "HBA",
            "Rotatable Bonds",
            "Aromatic Rings",
            "Heavy Atoms",
            "Fraction Csp3"
        ]
    )

    leaders = (
        filtered
        .sort_values(
            leader_property,
            ascending=False
        )
        .head(top_n)
        [
            [
                "Name",
                "Category",
                leader_property,
                "Molecular Weight",
                "LogP",
                "TPSA",
                "Lipinski"
            ]
        ]
    )

    st.dataframe(
        leaders,
        use_container_width=True,
        hide_index=True
    )

    # --------------------------------------------------------
    # LIPINSKI ANALYTICS
    # --------------------------------------------------------

    st.markdown("### Lipinski Statistics")

    pass_count = int(
        filtered["Lipinski"]
        .astype(str)
        .str.startswith("PASS")
        .sum()
    )

    fail_count = len(filtered) - pass_count

    l1, l2 = st.columns(2)

    with l1:

        lipinski_data = pd.DataFrame({
            "Status": ["PASS", "FAIL"],
            "Count": [pass_count, fail_count]
        })

        fig = px.pie(
            lipinski_data,
            names="Status",
            values="Count",
            title="Lipinski Rule of Five Distribution",
            hole=0.55
        )

        fig.update_layout(
            template="plotly_dark",
            paper_bgcolor="rgba(0,0,0,0)",
            plot_bgcolor="rgba(0,0,0,0)"
        )

        st.plotly_chart(
            fig,
            use_container_width=True
        )

    with l2:

        pass_rate = (
            pass_count / len(filtered) * 100
            if len(filtered) > 0
            else 0
        )

        st.metric(
            "Lipinski Pass Rate",
            f"{pass_rate:.1f}%"
        )

        st.write(
            "The Lipinski pass rate indicates how many compounds "
            "satisfy the four standard Rule of Five descriptor thresholds."
        )

        st.write(
            "This is a drug-likeness screening metric and does not "
            "guarantee biological activity, efficacy or safety."
        )

    # --------------------------------------------------------
    # CORRELATION MATRIX
    # --------------------------------------------------------

    st.markdown("### Molecular Descriptor Correlation")

    correlation_columns = [
        "Molecular Weight",
        "LogP",
        "TPSA",
        "HBD",
        "HBA",
        "Rotatable Bonds",
        "Aromatic Rings",
        "Heavy Atoms",
        "Fraction Csp3"
    ]

    correlation = filtered[correlation_columns].corr()

    fig = px.imshow(
        correlation,
        text_auto=".2f",
        aspect="auto",
        title="Descriptor Correlation Matrix"
    )

    fig.update_layout(
        template="plotly_dark",
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)"
    )

    st.plotly_chart(
        fig,
        use_container_width=True
    )

    # --------------------------------------------------------
    # DOWNLOAD ANALYTICS
    # --------------------------------------------------------

    st.markdown("### Export Analytics")

    csv_data = filtered.to_csv(index=False).encode("utf-8")

    st.download_button(
        "Download Filtered Analytics CSV",
        data=csv_data,
        file_name="pharmalens_analytics.csv",
        mime="text/csv"
    )


# ============================================================
# MOLECULAR ANALYSIS
# ============================================================

elif page == "Molecular Analysis":

    st.markdown(
        '<div class="main-title">Molecular Analysis</div>',
        unsafe_allow_html=True
    )

    if property_df.empty:
        st.warning("No compounds available.")
        st.stop()

    selected = st.selectbox(
        "Select compound",
        property_df["Name"].tolist(),
        key="analysis_compound"
    )

    row = property_df[
        property_df["Name"] == selected
    ].iloc[0]

    mol = mol_from_smiles(row["SMILES"])
    props = get_properties(mol)

    col1, col2 = st.columns([1, 2])

    with col1:
        st.image(
            molecule_image(mol),
            caption=selected,
            use_container_width=True
        )

    with col2:

        st.markdown("### Molecular Properties")

        c1, c2, c3 = st.columns(3)

        c1.metric(
            "Molecular Weight",
            f"{props['Molecular Weight']:.2f} Da"
        )

        c2.metric(
            "LogP",
            f"{props['LogP']:.2f}"
        )

        c3.metric(
            "TPSA",
            f"{props['TPSA']:.2f}"
        )

        c4, c5, c6 = st.columns(3)

        c4.metric("HBD", props["HBD"])
        c5.metric("HBA", props["HBA"])
        c6.metric("Rotatable Bonds", props["Rotatable Bonds"])

        c7, c8, c9 = st.columns(3)

        c7.metric("Aromatic Rings", props["Aromatic Rings"])
        c8.metric("Heavy Atoms", props["Heavy Atoms"])
        c9.metric(
            "Fraction Csp3",
            f"{props['Fraction Csp3']:.3f}"
        )

    st.markdown("### Molecular Identifiers")

    st.write("SMILES:", row["SMILES"])
    st.write("Formula:", props["Formula"])
    st.write("InChI Key:", molecular_inchi_key(mol))

    st.markdown("### Lipinski Rule of Five")

    st.dataframe(
        lipinski_table(mol),
        use_container_width=True,
        hide_index=True
    )


# ============================================================
# MOLECULAR INSIGHTS
# ============================================================

elif page == "Molecular Insights":

    st.markdown(
        '<div class="main-title">Molecular Insights</div>',
        unsafe_allow_html=True
    )

    selected = st.selectbox(
        "Select compound",
        property_df["Name"].tolist(),
        key="insights_compound"
    )

    row = property_df[
        property_df["Name"] == selected
    ].iloc[0]

    mol = mol_from_smiles(row["SMILES"])
    p = get_properties(mol)

    st.image(
        molecule_image(mol),
        width=450
    )

    st.markdown("### Property Snapshot")

    c1, c2, c3, c4 = st.columns(4)

    c1.metric("Formula", p["Formula"])
    c2.metric("Rotatable Bonds", p["Rotatable Bonds"])
    c3.metric("Aromatic Rings", p["Aromatic Rings"])
    c4.metric("Fraction Csp3", f"{p['Fraction Csp3']:.3f}")

    st.markdown("### Descriptor Profile")

    profile = pd.DataFrame({
        "Descriptor": [
            "Molecular Weight",
            "LogP",
            "TPSA",
            "HBD",
            "HBA",
            "Rotatable Bonds",
            "Aromatic Rings",
            "Heavy Atoms",
            "Fraction Csp3"
        ],
        "Value": [
            p["Molecular Weight"],
            p["LogP"],
            p["TPSA"],
            p["HBD"],
            p["HBA"],
            p["Rotatable Bonds"],
            p["Aromatic Rings"],
            p["Heavy Atoms"],
            p["Fraction Csp3"]
        ]
    })

    st.dataframe(
        profile,
        use_container_width=True,
        hide_index=True
    )


# ============================================================
# LEAD DISCOVERY
# ============================================================

elif page == "Lead Discovery":

    st.markdown(
        '<div class="main-title">Lead Discovery</div>',
        unsafe_allow_html=True
    )

    if property_df.empty:
        st.warning("No data available.")
        st.stop()

    st.markdown("### Molecular Filters")

    c1, c2 = st.columns(2)

    with c1:
        mw_range = st.slider(
            "Molecular Weight",
            0.0,
            float(max(600, property_df["Molecular Weight"].max())),
            (
                0.0,
                float(max(600, property_df["Molecular Weight"].max()))
            )
        )

        logp_range = st.slider(
            "LogP",
            -5.0,
            float(max(10, property_df["LogP"].max())),
            (
                -5.0,
                float(max(10, property_df["LogP"].max()))
            )
        )

        tpsa_range = st.slider(
            "TPSA",
            0.0,
            float(max(250, property_df["TPSA"].max())),
            (
                0.0,
                float(max(250, property_df["TPSA"].max()))
            )
        )

    with c2:
        hbd_range = st.slider(
            "HBD",
            0,
            15,
            (0, 15)
        )

        hba_range = st.slider(
            "HBA",
            0,
            20,
            (0, 20)
        )

        rot_range = st.slider(
            "Rotatable Bonds",
            0,
            20,
            (0, 20)
        )

    categories = st.multiselect(
        "Categories",
        sorted(property_df["Category"].unique()),
        default=sorted(property_df["Category"].unique())
    )

    lipinski_filter = st.selectbox(
        "Lipinski",
        ["All", "PASS", "FAIL"]
    )

    filtered = property_df[
        (property_df["Molecular Weight"] >= mw_range[0]) &
        (property_df["Molecular Weight"] <= mw_range[1]) &
        (property_df["LogP"] >= logp_range[0]) &
        (property_df["LogP"] <= logp_range[1]) &
        (property_df["TPSA"] >= tpsa_range[0]) &
        (property_df["TPSA"] <= tpsa_range[1]) &
        (property_df["HBD"] >= hbd_range[0]) &
        (property_df["HBD"] <= hbd_range[1]) &
        (property_df["HBA"] >= hba_range[0]) &
        (property_df["HBA"] <= hba_range[1]) &
        (property_df["Rotatable Bonds"] >= rot_range[0]) &
        (property_df["Rotatable Bonds"] <= rot_range[1]) &
        property_df["Category"].isin(categories)
    ]

    if lipinski_filter == "PASS":
        filtered = filtered[
            filtered["Lipinski"].str.startswith("PASS")
        ]

    elif lipinski_filter == "FAIL":
        filtered = filtered[
            ~filtered["Lipinski"].str.startswith("PASS")
        ]

    st.metric(
        "Matching Compounds",
        len(filtered)
    )

    st.dataframe(
        filtered,
        use_container_width=True,
        hide_index=True
    )


# ============================================================
# ADMET INTELLIGENCE
# ============================================================

elif page == "ADMET Intelligence":

    st.markdown(
        '<div class="main-title">ADMET Intelligence</div>',
        unsafe_allow_html=True
    )

    st.warning(
        "This is a descriptor-based screening system for research and "
        "education. It is not an experimental or clinically validated "
        "ADMET prediction model."
    )

    selected = st.selectbox(
        "Select compound",
        property_df["Name"].tolist(),
        key="admet_compound"
    )

    row = property_df[
        property_df["Name"] == selected
    ].iloc[0]

    mol = mol_from_smiles(row["SMILES"])
    p = get_properties(mol)

    results = run_admet(mol)
    radar_values = admet_radar_values(p)

    st.markdown("### Compound Overview")

    c1, c2 = st.columns([1, 2])

    with c1:
        st.image(
            molecule_image(mol),
            caption=selected,
            use_container_width=True
        )

    with c2:

        c21, c22, c23 = st.columns(3)

        c21.metric(
            "Molecular Weight",
            f"{p['Molecular Weight']:.2f}"
        )

        c22.metric(
            "LogP",
            f"{p['LogP']:.2f}"
        )

        c23.metric(
            "TPSA",
            f"{p['TPSA']:.2f}"
        )

        c24, c25, c26 = st.columns(3)

        c24.metric("HBD", p["HBD"])
        c25.metric("HBA", p["HBA"])
        c26.metric("Rotatable Bonds", p["Rotatable Bonds"])

    st.markdown("### ADMET Domain Assessment")

    domains = list(results.keys())

    cols = st.columns(5)

    for col, domain in zip(cols, domains):
        with col:
            st.metric(
                domain,
                results[domain][0]
            )

    st.markdown("### ADMET Profile")

    radar = go.Figure()

    radar_categories = list(radar_values.keys())
    radar_scores = list(radar_values.values())

    radar.add_trace(
        go.Scatterpolar(
            r=radar_scores + [radar_scores[0]],
            theta=radar_categories + [radar_categories[0]],
            fill="toself",
            name=selected
        )
    )

    radar.update_layout(
        polar=dict(
            radialaxis=dict(
                visible=True,
                range=[0, 100]
            )
        ),
        template="plotly_dark",
        paper_bgcolor="rgba(0,0,0,0)",
        title="Descriptor-Based ADMET Profile"
    )

    st.plotly_chart(
        radar,
        use_container_width=True
    )

    total, maximum = overall_admet_score(results)

    st.metric(
        "Descriptor Screening Score",
        f"{total}/{maximum}"
    )

    for domain, values in results.items():

        with st.expander(domain):

            st.write(
                f"Assessment: **{values[0]}**"
            )

            if values[1]:
                for flag in values[1]:
                    st.warning(flag)
            else:
                st.success(
                    "No descriptor-based flags identified."
                )

    st.markdown("### ADMET Descriptor Table")

    admet_descriptor_table = pd.DataFrame({
        "Descriptor": [
            "Molecular Weight",
            "LogP",
            "TPSA",
            "HBD",
            "HBA",
            "Rotatable Bonds",
            "Aromatic Rings",
            "Fraction Csp3"
        ],
        "Value": [
            round(p["Molecular Weight"], 2),
            round(p["LogP"], 2),
            round(p["TPSA"], 2),
            p["HBD"],
            p["HBA"],
            p["Rotatable Bonds"],
            p["Aromatic Rings"],
            round(p["Fraction Csp3"], 3)
        ]
    })

    st.dataframe(
        admet_descriptor_table,
        use_container_width=True,
        hide_index=True
    )


# ============================================================
# ADMET COMPARISON
# ============================================================

elif page == "ADMET Comparison":

    st.markdown(
        '<div class="main-title">ADMET Comparison</div>',
        unsafe_allow_html=True
    )

    selected = st.multiselect(
        "Select 2–3 compounds",
        property_df["Name"].tolist(),
        max_selections=3,
        key="admet_compare"
    )

    if len(selected) < 2:
        st.info("Select at least two compounds.")
        st.stop()

    comparison_results = {}

    for name in selected:
        row = property_df[
            property_df["Name"] == name
        ].iloc[0]

        mol = mol_from_smiles(row["SMILES"])

        comparison_results[name] = {
            "mol": mol,
            "properties": get_properties(mol),
            "admet": run_admet(mol)
        }

    st.markdown("### ADMET Domain Comparison")

    comparison_rows = []

    for domain in [
        "Absorption",
        "Distribution",
        "Metabolism",
        "Excretion",
        "Toxicity"
    ]:

        row_data = {"Domain": domain}

        for name in selected:
            row_data[name] = comparison_results[name]["admet"][domain][0]

        comparison_rows.append(row_data)

    st.dataframe(
        pd.DataFrame(comparison_rows),
        use_container_width=True,
        hide_index=True
    )

    st.markdown("### Descriptor Screening Scores")

    score_rows = []

    for name in selected:
        results = comparison_results[name]["admet"]

        total, maximum = overall_admet_score(results)

        score_rows.append({
            "Compound": name,
            "Score": total,
            "Maximum": maximum
        })

    st.dataframe(
        pd.DataFrame(score_rows),
        use_container_width=True,
        hide_index=True
    )

    st.markdown("### ADMET Radar Comparison")

    radar = go.Figure()

    for name in selected:

        values = admet_radar_values(
            comparison_results[name]["properties"]
        )

        categories = list(values.keys())
        scores = list(values.values())

        radar.add_trace(
            go.Scatterpolar(
                r=scores + [scores[0]],
                theta=categories + [categories[0]],
                fill="toself",
                name=name
            )
        )

    radar.update_layout(
        polar=dict(
            radialaxis=dict(
                visible=True,
                range=[0, 100]
            )
        ),
        template="plotly_dark",
        paper_bgcolor="rgba(0,0,0,0)"
    )

    st.plotly_chart(
        radar,
        use_container_width=True
    )

    st.markdown("### Molecular Descriptor Comparison")

    descriptor_rows = []

    for name in selected:

        p = comparison_results[name]["properties"]

        descriptor_rows.append({
            "Compound": name,
            "MW": round(p["Molecular Weight"], 2),
            "LogP": round(p["LogP"], 2),
            "TPSA": round(p["TPSA"], 2),
            "HBD": p["HBD"],
            "HBA": p["HBA"],
            "Rotatable Bonds": p["Rotatable Bonds"],
            "Aromatic Rings": p["Aromatic Rings"],
            "Fraction Csp3": round(p["Fraction Csp3"], 3)
        })

    st.dataframe(
        pd.DataFrame(descriptor_rows),
        use_container_width=True,
        hide_index=True
    )

    st.markdown("### Molecular Structures")

    cols = st.columns(len(selected))

    for col, name in zip(cols, selected):

        with col:
            st.subheader(name)

            st.image(
                molecule_image(
                    comparison_results[name]["mol"]
                ),
                use_container_width=True
            )


# ============================================================
# COMPOUND DATABASE
# ============================================================

elif page == "Compound Database":

    st.markdown(
        '<div class="main-title">Compound Database</div>',
        unsafe_allow_html=True
    )

    st.dataframe(
        property_df,
        use_container_width=True,
        hide_index=True
    )

    csv = property_df.to_csv(index=False).encode("utf-8")

    st.download_button(
        "Download Compound Database",
        data=csv,
        file_name="compounds.csv",
        mime="text/csv"
    )


# ============================================================
# COMPOUND COMPARISON
# ============================================================

elif page == "Compound Comparison":

    st.markdown(
        '<div class="main-title">Compound Comparison</div>',
        unsafe_allow_html=True
    )

    selected = st.multiselect(
        "Select 2–3 compounds",
        property_df["Name"].tolist(),
        max_selections=3,
        key="compound_compare"
    )

    if len(selected) < 2:
        st.info("Select at least two compounds.")
        st.stop()

    rows = []

    cols = st.columns(len(selected))

    for col, name in zip(cols, selected):

        row = property_df[
            property_df["Name"] == name
        ].iloc[0]

        mol = mol_from_smiles(row["SMILES"])
        p = get_properties(mol)

        with col:

            st.subheader(name)

            st.image(
                molecule_image(mol),
                use_container_width=True
            )

            st.metric(
                "MW",
                f"{p['Molecular Weight']:.2f}"
            )

            st.metric(
                "LogP",
                f"{p['LogP']:.2f}"
            )

            st.metric(
                "TPSA",
                f"{p['TPSA']:.2f}"
            )

            st.metric(
                "Lipinski",
                lipinski_status(mol)
            )

        rows.append({
            "Compound": name,
            "MW": round(p["Molecular Weight"], 2),
            "LogP": round(p["LogP"], 2),
            "TPSA": round(p["TPSA"], 2),
            "HBD": p["HBD"],
            "HBA": p["HBA"],
            "Rotatable Bonds": p["Rotatable Bonds"],
            "Aromatic Rings": p["Aromatic Rings"],
            "Lipinski": lipinski_status(mol)
        })

    st.markdown("### Comparison Table")

    st.dataframe(
        pd.DataFrame(rows),
        use_container_width=True,
        hide_index=True
    )


# ============================================================
# CHEMICAL SPACE
# ============================================================

elif page == "Chemical Space":

    st.markdown(
        '<div class="main-title">Chemical Space</div>',
        unsafe_allow_html=True
    )

    if len(records) < 3:
        st.warning("Not enough compounds for chemical-space analysis.")
        st.stop()

    fp_matrix = []

    names = []
    categories = []
    metadata = []

    for record in records:

        mol = mol_from_smiles(record["SMILES"])

        arr = fingerprint_array(mol)

        if arr is not None:

            fp_matrix.append(arr)
            names.append(record["Name"])
            categories.append(record["Category"])
            metadata.append(record)

    fp_matrix = np.array(fp_matrix, dtype=float)

    centered = fp_matrix - fp_matrix.mean(axis=0)

    try:
        u, s, vt = np.linalg.svd(
            centered,
            full_matrices=False
        )

        coords = u[:, :2] * s[:2]

    except Exception:
        coords = np.zeros((len(fp_matrix), 2))

    space_df = pd.DataFrame({
        "PC1": coords[:, 0],
        "PC2": coords[:, 1],
        "Name": names,
        "Category": categories
    })

    for key in [
        "Molecular Weight",
        "LogP",
        "TPSA",
        "HBD",
        "HBA"
    ]:
        space_df[key] = [
            item[key]
            for item in metadata
        ]

    fig = px.scatter(
        space_df,
        x="PC1",
        y="PC2",
        color="Category",
        hover_name="Name",
        hover_data=[
            "Molecular Weight",
            "LogP",
            "TPSA",
            "HBD",
            "HBA"
        ],
        title="Chemical Space Map"
    )

    fig.update_layout(
        template="plotly_dark",
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)",
        xaxis_title="Structural Coordinate 1",
        yaxis_title="Structural Coordinate 2"
    )

    st.plotly_chart(
        fig,
        use_container_width=True
    )

    st.info(
        "The axes represent structural coordinates derived from molecular "
        "fingerprints. They are not direct physical or biological properties."
    )


# ============================================================
# CUSTOM MOLECULE
# ============================================================

elif page == "Custom Molecule":

    st.markdown(
        '<div class="main-title">Custom Molecule Analysis</div>',
        unsafe_allow_html=True
    )

    smiles = st.text_input(
        "Enter SMILES",
        value="C1=CC(=C(C=C1CCN)O)O",
        key="custom_smiles"
    )

    mol = mol_from_smiles(smiles)

    if mol is None:

        st.error(
            "Invalid SMILES. Please enter a valid molecular structure."
        )

    else:

        p = get_properties(mol)

        col1, col2 = st.columns([1, 2])

        with col1:
            st.image(
                molecule_image(mol),
                use_container_width=True
            )

        with col2:

            c1, c2, c3 = st.columns(3)

            c1.metric(
                "Molecular Weight",
                f"{p['Molecular Weight']:.2f}"
            )

            c2.metric(
                "LogP",
                f"{p['LogP']:.2f}"
            )

            c3.metric(
                "TPSA",
                f"{p['TPSA']:.2f}"
            )

            c4, c5, c6 = st.columns(3)

            c4.metric("HBD", p["HBD"])
            c5.metric("HBA", p["HBA"])
            c6.metric(
                "Rotatable Bonds",
                p["Rotatable Bonds"]
            )

            c7, c8, c9 = st.columns(3)

            c7.metric(
                "Aromatic Rings",
                p["Aromatic Rings"]
            )

            c8.metric(
                "Heavy Atoms",
                p["Heavy Atoms"]
            )

            c9.metric(
                "Fraction Csp3",
                f"{p['Fraction Csp3']:.3f}"
            )

        st.markdown("### Molecular Identifiers")

        st.write("Formula:", p["Formula"])
        st.write("InChI Key:", molecular_inchi_key(mol))

        st.markdown("### Lipinski Rule of Five")

        st.dataframe(
            lipinski_table(mol),
            use_container_width=True,
            hide_index=True
        )


# ============================================================
# SIMILARITY SEARCH
# ============================================================

elif page == "Similarity Search":

    st.markdown(
        '<div class="main-title">Similarity Search</div>',
        unsafe_allow_html=True
    )

    smiles = st.text_input(
        "Query SMILES",
        value="CC(=O)OC1=CC=CC=C1C(=O)O",
        key="similarity_smiles"
    )

    threshold = st.slider(
        "Similarity Threshold",
        0.0,
        1.0,
        0.15,
        0.01
    )

    max_results = st.slider(
        "Maximum Results",
        1,
        min(25, max(1, len(property_df))),
        min(10, max(1, len(property_df)))
    )

    mol = mol_from_smiles(smiles)

    if mol is None:

        st.error("Invalid SMILES.")

    else:

        query_fp = fingerprint(mol)

        results = []

        for record in records:

            target_mol = mol_from_smiles(record["SMILES"])
            target_fp = fingerprint(target_mol)

            similarity = DataStructs.TanimotoSimilarity(
                query_fp,
                target_fp
            )

            if similarity >= threshold:

                results.append({
                    "Name": record["Name"],
                    "Category": record["Category"],
                    "Similarity": similarity,
                    "Molecular Weight": record["Molecular Weight"],
                    "LogP": record["LogP"],
                    "TPSA": record["TPSA"]
                })

        results_df = pd.DataFrame(results)

        if not results_df.empty:

            results_df = (
                results_df
                .sort_values(
                    "Similarity",
                    ascending=False
                )
                .head(max_results)
            )

            results_df["Similarity"] = (
                results_df["Similarity"] * 100
            ).round(2)

            st.dataframe(
                results_df,
                use_container_width=True,
                hide_index=True
            )

        else:
            st.info(
                "No compounds passed the selected similarity threshold."
            )


# ============================================================
# SIMILARITY EXPLORER
# ============================================================

elif page == "Similarity Explorer":

    st.markdown(
        '<div class="main-title">Similarity Explorer</div>',
        unsafe_allow_html=True
    )

    selected = st.selectbox(
        "Select reference compound",
        property_df["Name"].tolist(),
        key="similarity_reference"
    )

    threshold = st.slider(
        "Similarity threshold",
        0.0,
        1.0,
        0.10,
        0.01
    )

    row = property_df[
        property_df["Name"] == selected
    ].iloc[0]

    query_mol = mol_from_smiles(row["SMILES"])
    query_fp = fingerprint(query_mol)

    results = []

    for record in records:

        if record["Name"] == selected:
            continue

        mol = mol_from_smiles(record["SMILES"])

        similarity = DataStructs.TanimotoSimilarity(
            query_fp,
            fingerprint(mol)
        )

        if similarity >= threshold:

            results.append({
                "Name": record["Name"],
                "Category": record["Category"],
                "Similarity": round(similarity, 4),
                "Molecular Weight": round(record["Molecular Weight"], 2),
                "LogP": round(record["LogP"], 2),
                "TPSA": round(record["TPSA"], 2)
            })

    results_df = pd.DataFrame(results)

    if results_df.empty:
        st.info("No similar compounds found.")
    else:

        results_df = results_df.sort_values(
            "Similarity",
            ascending=False
        )

        st.dataframe(
            results_df,
            use_container_width=True,
            hide_index=True
        )

        fig = px.bar(
            results_df.head(15),
            x="Similarity",
            y="Name",
            orientation="h",
            title=f"Similarity to {selected}"
        )

        fig.update_layout(
            template="plotly_dark",
            paper_bgcolor="rgba(0,0,0,0)",
            plot_bgcolor="rgba(0,0,0,0)"
        )

        st.plotly_chart(
            fig,
            use_container_width=True
        )


# ============================================================
# DRUG LIKENESS
# ============================================================

elif page == "Drug-Likeness":

    st.markdown(
        '<div class="main-title">Drug-Likeness</div>',
        unsafe_allow_html=True
    )

    selected = st.selectbox(
        "Select compound",
        property_df["Name"].tolist(),
        key="druglikeness_compound"
    )

    row = property_df[
        property_df["Name"] == selected
    ].iloc[0]

    mol = mol_from_smiles(row["SMILES"])

    p = get_properties(mol)

    st.markdown("### Lipinski Rule of Five")

    st.dataframe(
        lipinski_table(mol),
        use_container_width=True,
        hide_index=True
    )

    st.markdown("### Drug-Likeness Descriptors")

    c1, c2, c3, c4 = st.columns(4)

    c1.metric(
        "MW",
        f"{p['Molecular Weight']:.2f}"
    )

    c2.metric(
        "LogP",
        f"{p['LogP']:.2f}"
    )

    c3.metric(
        "HBD",
        p["HBD"]
    )

    c4.metric(
        "HBA",
        p["HBA"]
    )

    if lipinski_status(mol) == "PASS":
        st.success(
            "This compound satisfies the implemented Lipinski Rule of Five thresholds."
        )
    else:
        st.warning(
            "This compound violates one or more implemented Lipinski thresholds."
        )


# ============================================================
# BATCH SMILES ANALYZER
# ============================================================

elif page == "Batch SMILES Analyzer":

    st.markdown(
        '<div class="main-title">Batch SMILES Analyzer</div>',
        unsafe_allow_html=True
    )

    st.write(
        "Enter one SMILES string per line."
    )

    smiles_text = st.text_area(
        "SMILES input",
        height=220,
        placeholder="CCO\nCC(=O)OC1=CC=CC=C1C(=O)O\n..."
    )

    if st.button("Analyze Batch"):

        valid_df, invalid_df = analyze_batch_smiles(
            smiles_text
        )

        st.session_state["batch_valid"] = valid_df
        st.session_state["batch_invalid"] = invalid_df

    if "batch_valid" in st.session_state:

        valid_df = st.session_state["batch_valid"]
        invalid_df = st.session_state["batch_invalid"]

        c1, c2, c3, c4 = st.columns(4)

        c1.metric(
            "Valid",
            len(valid_df)
        )

        c2.metric(
            "Invalid",
            len(invalid_df)
        )

        if not valid_df.empty:

            c3.metric(
                "Average MW",
                f"{valid_df['MW'].mean():.2f}"
            )

            c4.metric(
                "Average LogP",
                f"{valid_df['LogP'].mean():.2f}"
            )

            st.markdown("### Valid Compounds")

            st.dataframe(
                valid_df,
                use_container_width=True,
                hide_index=True
            )

            pass_count = (
                valid_df["Lipinski"]
                .astype(str)
                .str.startswith("PASS")
                .sum()
            )

            st.metric(
                "Lipinski PASS Count",
                int(pass_count)
            )

            csv = valid_df.to_csv(
                index=False
            ).encode("utf-8")

            st.download_button(
                "Download Batch CSV",
                data=csv,
                file_name="batch_analysis.csv",
                mime="text/csv"
            )

        if not invalid_df.empty:

            st.markdown("### Invalid SMILES")

            st.dataframe(
                invalid_df,
                use_container_width=True,
                hide_index=True
            )


# ============================================================
# COMPOUND REPORT GENERATOR
# ============================================================

elif page == "Compound Report Generator":

    st.markdown(
        '<div class="main-title">Compound Report Generator</div>',
        unsafe_allow_html=True
    )

    st.markdown(
        '<div class="subtitle">Generate a detailed molecular research report</div>',
        unsafe_allow_html=True
    )

    selected = st.selectbox(
        "Select compound",
        property_df["Name"].tolist(),
        key="report_compound"
    )

    row = property_df[
        property_df["Name"] == selected
    ].iloc[0]

    mol = mol_from_smiles(row["SMILES"])

    p = get_properties(mol)

    admet_results = run_admet(mol)
    admet_values = admet_radar_values(p)

    st.markdown("### Compound Overview")

    c1, c2 = st.columns([1, 2])

    with c1:
        st.image(
            molecule_image(mol),
            use_container_width=True
        )

    with c2:

        st.write("Compound:", selected)
        st.write("Category:", row["Category"])
        st.write("SMILES:", row["SMILES"])
        st.write("Formula:", p["Formula"])
        st.write(
            "InChI Key:",
            molecular_inchi_key(mol)
        )

    st.markdown("### Molecular Descriptors")

    descriptor_table = pd.DataFrame({
        "Descriptor": [
            "Molecular Weight",
            "LogP",
            "TPSA",
            "HBD",
            "HBA",
            "Rotatable Bonds",
            "Aromatic Rings",
            "Heavy Atoms",
            "Fraction Csp3"
        ],
        "Value": [
            round(p["Molecular Weight"], 2),
            round(p["LogP"], 2),
            round(p["TPSA"], 2),
            p["HBD"],
            p["HBA"],
            p["Rotatable Bonds"],
            p["Aromatic Rings"],
            p["Heavy Atoms"],
            round(p["Fraction Csp3"], 3)
        ]
    })

    st.dataframe(
        descriptor_table,
        use_container_width=True,
        hide_index=True
    )

    st.markdown("### Lipinski Rule of Five")

    st.dataframe(
        lipinski_table(mol),
        use_container_width=True,
        hide_index=True
    )

    st.markdown("### ADMET Screening")

    admet_rows = []

    for domain, values in admet_results.items():

        admet_rows.append({
            "Domain": domain,
            "Assessment": values[0],
            "Score": values[2]
        })

    st.dataframe(
        pd.DataFrame(admet_rows),
        use_container_width=True,
        hide_index=True
    )

    total, maximum = overall_admet_score(
        admet_results
    )

    st.metric(
        "Descriptor Screening Score",
        f"{total}/{maximum}"
    )

    st.markdown("### ADMET Profile")

    radar = go.Figure()

    categories = list(admet_values.keys())
    scores = list(admet_values.values())

    radar.add_trace(
        go.Scatterpolar(
            r=scores + [scores[0]],
            theta=categories + [categories[0]],
            fill="toself",
            name=selected
        )
    )

    radar.update_layout(
        polar=dict(
            radialaxis=dict(
                visible=True,
                range=[0, 100]
            )
        ),
        template="plotly_dark",
        paper_bgcolor="rgba(0,0,0,0)"
    )

    st.plotly_chart(
        radar,
        use_container_width=True
    )

    st.markdown("### ADMET Flags")

    for domain, values in admet_results.items():

        if values[1]:

            with st.expander(domain):

                for flag in values[1]:
                    st.warning(flag)

    # --------------------------------------------------------
    # TXT REPORT
    # --------------------------------------------------------

    report_text = []

    report_text.append(
        "PHARMALENS AI - COMPOUND REPORT"
    )

    report_text.append("=" * 50)
    report_text.append("")

    report_text.append(
        f"Compound: {selected}"
    )

    report_text.append(
        f"Category: {row['Category']}"
    )

    report_text.append(
        f"SMILES: {row['SMILES']}"
    )

    report_text.append(
        f"Formula: {p['Formula']}"
    )

    report_text.append(
        f"InChI Key: {molecular_inchi_key(mol)}"
    )

    report_text.append("")

    report_text.append(
        "MOLECULAR DESCRIPTORS"
    )

    report_text.append("-" * 30)

    for key in [
        "Molecular Weight",
        "LogP",
        "TPSA",
        "HBD",
        "HBA",
        "Rotatable Bonds",
        "Aromatic Rings",
        "Heavy Atoms",
        "Fraction Csp3"
    ]:
        report_text.append(
            f"{key}: {p[key]}"
        )

    report_text.append("")

    report_text.append(
        f"Lipinski: {lipinski_status(mol)}"
    )

    report_text.append("")

    report_text.append(
        "ADMET SCREENING"
    )

    report_text.append("-" * 30)

    for domain, values in admet_results.items():

        report_text.append(
            f"{domain}: {values[0]} | Score: {values[2]}"
        )

        for flag in values[1]:
            report_text.append(
                f"  - {flag}"
            )

    report_text.append("")

    report_text.append(
        f"Descriptor Screening Score: {total}/{maximum}"
    )

    report_text.append("")

    report_text.append(
        "DISCLAIMER: Descriptor-based screening is intended "
        "for research and educational purposes and is not "
        "a validated clinical or experimental prediction."
    )

    report_text_content = "\n".join(report_text)

    st.download_button(
        "Download Text Report",
        data=report_text_content,
        file_name=f"{selected.replace(' ', '_')}_report.txt",
        mime="text/plain"
    )

    # --------------------------------------------------------
    # PDF REPORT
    # --------------------------------------------------------

    st.markdown("### Professional PDF Report")

    pdf_bytes = create_pdf_report(
        selected,
        row["Category"],
        row["SMILES"],
        mol,
        admet_results,
        admet_values
    )

    st.download_button(
        "Download Professional PDF Report",
        data=pdf_bytes,
        file_name=f"{selected.replace(' ', '_')}_PharmaLens_Report.pdf",
        mime="application/pdf"
    )


# ============================================================
# FOOTER
# ============================================================

st.markdown(
    """
    <div class="footer">
        PharmaLens AI | Molecular Intelligence & Drug Discovery Platform
    </div>
    """,
    unsafe_allow_html=True
)