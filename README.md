# PharmaLens AI

### AI-Powered Molecular Analysis & Drug Discovery Platform

PharmaLens AI is an interactive cheminformatics and drug-discovery platform built with **Python, RDKit, Pandas, NumPy, Plotly, and Streamlit**.

The platform provides a unified environment for molecular property analysis, drug-likeness evaluation, molecular similarity searching, chemical-space exploration, lead discovery, ADMET descriptor-based screening, batch analysis, and professional compound reporting.

> **Note:** PharmaLens AI is an educational and research-oriented screening platform. Its ADMET assessments are descriptor-based computational interpretations and should not be considered experimentally validated, clinical, or regulatory predictions.

---

## Project Overview

Drug discovery involves evaluating large numbers of chemical compounds using molecular descriptors, structural fingerprints, physicochemical properties, and drug-likeness rules.

PharmaLens AI brings several of these computational workflows into a single interactive application.

Users can:

* Analyze individual molecules from SMILES
* Calculate important molecular descriptors
* Visualize 2D molecular structures
* Evaluate Lipinski's Rule of Five
* Search compounds using molecular similarity
* Explore chemical space
* Compare multiple compounds
* Filter compounds for lead-discovery workflows
* Perform descriptor-based ADMET screening
* Analyze multiple SMILES simultaneously
* Generate professional compound reports
* Explore drug-discovery datasets through interactive analytics

---

# Key Features

## 1. Molecular Analysis

Analyze a molecule directly from its SMILES representation.

Calculated properties include:

* Molecular Weight
* LogP
* TPSA
* Hydrogen Bond Donors
* Hydrogen Bond Acceptors
* Rotatable Bonds
* Aromatic Rings
* Heavy Atoms
* Fraction Csp3
* Molecular Formula
* InChI Key

A 2D molecular structure is also generated using RDKit.

---

## 2. Drug-Likeness Analysis

PharmaLens evaluates compounds using **Lipinski's Rule of Five**.

The application evaluates:

* Molecular Weight
* LogP
* Hydrogen Bond Donors
* Hydrogen Bond Acceptors

The result is presented as an overall drug-likeness screening status along with the individual rule checks.

---

## 3. Molecular Similarity Search

PharmaLens uses **Morgan molecular fingerprints** and **Tanimoto similarity** to identify structurally similar compounds.

Current fingerprint configuration:

* Fingerprint type: Morgan
* Radius: 2
* Fingerprint size: 2048 bits
* Similarity metric: Tanimoto similarity

Users can enter a query SMILES and search the compound database for structurally similar molecules.

---

## 4. Compound Comparison

Compare multiple compounds side-by-side using:

* Molecular Weight
* LogP
* TPSA
* HBD
* HBA
* Rotatable Bonds
* Aromatic Rings
* Lipinski status
* Molecular structures

This provides a convenient way to inspect molecular-property differences between compounds.

---

## 5. Chemical Space Exploration

PharmaLens provides a 2D visualization of molecular chemical space using molecular fingerprint-derived structural features.

Interactive Plotly visualization allows users to:

* Zoom
* Pan
* Hover over compounds
* Inspect compound names
* View molecular categories
* Examine selected molecular descriptors

The plotted coordinates represent structural relationships and should not be interpreted as direct physical properties.

---

## 6. Lead Discovery

The Lead Discovery module allows compounds to be filtered using molecular-property criteria.

Available filters include:

* Molecular Weight
* LogP
* TPSA
* HBD
* HBA
* Rotatable Bonds
* Compound Category
* Lipinski status

The filtered compound set can then be inspected as a potential computational starting point for further investigation.

---

## 7. ADMET Intelligence

PharmaLens provides a descriptor-based ADMET screening layer covering:

### Absorption

Evaluates molecular descriptors associated with absorption-related behavior.

### Distribution

Screens descriptor patterns associated with molecular distribution.

### Metabolism

Evaluates structural and descriptor-based indicators related to metabolic concerns.

### Excretion

Provides a descriptor-based interpretation of molecular properties relevant to elimination.

### Toxicity

Flags selected molecular-property patterns associated with potential toxicity concerns.

The module provides:

* Domain-level status
* Descriptor flags
* Screening scores
* ADMET profile visualization
* Molecular descriptor summary

These results are intended for computational screening and educational interpretation, not clinical prediction.

---

## 8. ADMET Comparison

Multiple compounds can be compared across ADMET domains.

The comparison includes:

* Absorption
* Distribution
* Metabolism
* Excretion
* Toxicity
* Descriptor screening score
* Molecular descriptors
* Descriptor flags
* ADMET profile visualization

---

## 9. Batch SMILES Analyzer

Users can provide multiple SMILES strings for simultaneous analysis.

For each valid molecule, PharmaLens calculates:

* Molecular Formula
* Molecular Weight
* LogP
* TPSA
* HBD
* HBA
* Rotatable Bonds
* Aromatic Rings
* Heavy Atoms
* Fraction Csp3
* Lipinski status

The resulting dataset can be exported as a CSV file.

Invalid SMILES entries are separately identified.

---

## 10. Drug Discovery Analytics

The analytics dashboard provides interactive exploration of the compound database.

Available analyses include:

* Molecular Weight distribution
* LogP distribution
* TPSA distribution
* Fraction Csp3 distribution
* Property relationship analysis
* Category-level statistics
* Compound ranking
* Lipinski distribution
* Lipinski pass rate
* Descriptor correlation matrix
* Filtered dataset export

Supported property relationships include:

* Molecular Weight vs LogP
* Molecular Weight vs TPSA
* LogP vs TPSA
* TPSA vs HBD
* Molecular Weight vs Rotatable Bonds

---

## 11. Compound Report Generator

PharmaLens can generate detailed compound reports containing:

* Compound information
* Molecular structure
* SMILES
* Molecular Formula
* InChI Key
* Molecular descriptors
* Lipinski analysis
* ADMET screening
* Descriptor flags
* ADMET profile
* Scientific disclaimer

Reports can be exported as text files or generated as professional PDF reports.

---

# Technology Stack

| Technology | Purpose                                |
| ---------- | -------------------------------------- |
| Python     | Core programming language              |
| Streamlit  | Interactive web application            |
| RDKit      | Cheminformatics and molecular analysis |
| Pandas     | Data processing                        |
| NumPy      | Numerical computation                  |
| Plotly     | Interactive visualization              |
| ReportLab  | PDF report generation                  |
| Pytest     | Automated testing                      |
| Git        | Version control                        |
| GitHub     | Project hosting                        |

---

# Molecular Analysis Methodology

PharmaLens uses RDKit to convert SMILES strings into molecular representations.

The general computational workflow is:

```text
SMILES Input
     |
     v
RDKit Molecular Representation
     |
     +--------------------+
     |                    |
     v                    v
Molecular Descriptors   Molecular Fingerprint
     |                    |
     |                    v
     |              Similarity Search
     |
     +----------+
     |          |
     v          v
Lipinski    ADMET Screening
     |
     v
Drug-Likeness Assessment
```

---

# Molecular Descriptors

The platform calculates several commonly used molecular descriptors.

### Molecular Weight

Represents the molecular mass of a compound in Daltons.

### LogP

Provides an estimate of lipophilicity based on the compound's partitioning behavior.

### TPSA

Topological Polar Surface Area provides an estimate of molecular polar surface characteristics.

### HBD

Hydrogen Bond Donors.

### HBA

Hydrogen Bond Acceptors.

### Rotatable Bonds

Provides an indication of molecular flexibility.

### Aromatic Rings

Counts aromatic ring systems identified by RDKit.

### Heavy Atoms

Counts non-hydrogen atoms in the molecule.

### Fraction Csp3

Represents the fraction of carbon atoms with sp3 hybridization.

---

# Molecular Similarity Methodology

PharmaLens generates Morgan fingerprints using:

```text
Radius = 2
Fingerprint Size = 2048 bits
```

Similarity between molecular fingerprints is calculated using the **Tanimoto coefficient**.

Conceptually:

```text
Similarity =
Common fingerprint features
---------------------------
Total fingerprint features
```

Higher similarity values indicate greater fingerprint overlap.

---

# Lipinski Rule of Five

PharmaLens uses four commonly evaluated Lipinski descriptors:

```text
Molecular Weight
LogP
Hydrogen Bond Donors
Hydrogen Bond Acceptors
```

The module reports whether the compound satisfies the implemented screening criteria.

Lipinski's rules are useful for early-stage drug-likeness screening but do not guarantee biological activity, efficacy, safety, or clinical success.

---

# ADMET Screening

The ADMET module uses molecular descriptors and structural information to produce rule-based computational screening interpretations.

The five domains evaluated are:

```text
Absorption
Distribution
Metabolism
Excretion
Toxicity
```

The output is intended to help users identify compounds that may require further investigation.

These results should be treated as **hypothesis-generating computational assessments**, not experimental measurements.

---

# Dataset

The project contains a curated compound database used for demonstration and computational analysis.

The current application includes:

* Compound names
* SMILES representations
* Compound categories
* Molecular descriptors generated through RDKit

The dataset is intended for educational and research-oriented computational workflows.

---

# Project Structure

```text
pharmalens-ai/
│
├── app.py
├── create_database.py
├── main.py
├── README.md
├── requirements.txt
├── .gitignore
│
├── data/
│   ├── compounds.csv
│   └── drugs.json
│
├── reports/
│   ├── aspirin_structure.png
│   └── caffeine_structure.png
│
├── src/
│   ├── drug_database.py
│   ├── drug_input.py
│   ├── drug_likeness.py
│   ├── molecular_fingerprint.py
│   ├── molecular_similarity.py
│   ├── molecule_analyzer.py
│   ├── molecule_visualizer.py
│   └── similarity_search.py
│
└── tests/
    └── test_pharmalens.py
```

---

# Installation

Clone the repository:

```bash
git clone https://github.com/YOUR-USERNAME/pharmalens-ai.git
```

Move into the project directory:

```bash
cd pharmalens-ai
```

Create a virtual environment:

```bash
python -m venv .venv
```

Activate it on Windows:

```powershell
.venv\Scripts\Activate.ps1
```

Install the dependencies:

```bash
pip install -r requirements.txt
```

---

# Running the Application

Start the Streamlit application:

```bash
streamlit run app.py
```

The application will open locally at:

```text
http://localhost:8501
```

---

# Testing

The project includes automated tests using Pytest.

Run:

```bash
pytest
```

The tests cover important PharmaLens molecular-analysis functionality.

---

# Example Workflow

A typical PharmaLens workflow can be:

```text
1. Enter a molecule using SMILES
          ↓
2. Generate molecular structure
          ↓
3. Calculate molecular descriptors
          ↓
4. Evaluate Lipinski drug-likeness
          ↓
5. Search molecular similarity
          ↓
6. Explore chemical space
          ↓
7. Evaluate descriptor-based ADMET profile
          ↓
8. Compare candidate compounds
          ↓
9. Generate a compound report
```

---

# Scientific Limitations

PharmaLens AI is designed as a computational screening and educational platform.

The results should not be interpreted as:

* Clinical predictions
* Experimental measurements
* Regulatory assessments
* Confirmed pharmacokinetic properties
* Confirmed toxicity results
* Proof of therapeutic efficacy

Computational screening can help prioritize molecules for further investigation, but experimental validation remains essential.

---

# Future Scope

Potential future development areas include:

* Larger compound libraries
* Additional molecular descriptors
* Advanced molecular visualization
* External chemical databases
* Molecular docking integration
* QSAR modelling
* More extensive ADMET datasets
* Machine-learning-based property prediction
* Structure-based drug discovery workflows
* Cloud deployment
* API-based molecular analysis

---

# Project Status

**Current status: Functional prototype / portfolio project**

The current version focuses on:

* Cheminformatics
* Molecular descriptor analysis
* Molecular similarity
* Drug-likeness screening
* Chemical-space visualization
* Lead discovery workflows
* Descriptor-based ADMET screening
* Data analytics
* Automated testing
* Professional reporting

---

# Disclaimer

PharmaLens AI is an educational and research-oriented software project.

Computational results are dependent on the underlying molecular representations, descriptors, rules, datasets, and algorithms used by the application. They should not be used as a substitute for experimental studies, professional scientific judgment, clinical evaluation, or regulatory assessment.

---

# Author

**Pankaj Chauhan**

B.Tech Bioinformatics

Interests:

* Bioinformatics
* Computational Biology
* Cheminformatics
* Drug Discovery
* Artificial Intelligence
* Machine Learning

---

## PharmaLens AI

**From molecular structure to computational drug-discovery insight.**
