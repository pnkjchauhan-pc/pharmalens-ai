from src.molecule_visualizer import generate_molecule_image
from src.drug_input import validate_smiles
from src.molecule_analyzer import analyze_molecule
from src.drug_likeness import calculate_drug_likeness
from src.drug_database import get_drug, list_drugs


def run_pharmalens(drug_name, smiles):

    print("\n======================================")
    print("          PHARMALENS AI")
    print("     Molecular Intelligence Engine")
    print("======================================\n")

    print(f"Drug: {drug_name}")
    print(f"SMILES: {smiles}\n")

    # 1. Validate molecule
    print("🔍 Validating molecule...")

    if not validate_smiles(smiles):
        print("❌ Invalid SMILES. Analysis stopped.")
        return

    print("✅ Molecule validated successfully.")

    # 2. Molecular analysis
    print("\n🧬 Analyzing molecular properties...\n")

    molecular_profile = analyze_molecule(smiles)

    print("---------- MOLECULAR PROFILE ----------")

    for property_name, value in molecular_profile.items():
        print(f"{property_name}: {value}")

    # 3. Drug-likeness
    print("\n💊 Evaluating drug-likeness...\n")

    drug_likeness = calculate_drug_likeness(smiles)

    print("---------- DRUG-LIKENESS ----------")

    for property_name, details in drug_likeness["checks"].items():

        status = "PASS" if details["passes"] else "FAIL"

        print(
            f"{property_name}: "
            f"{details['value']} → {status}"
        )

    print(
        f"\nRule of Five violations: "
        f"{drug_likeness['violations']}"
    )

    if drug_likeness["rule_of_five_pass"]:
        print("Overall Rule of Five screen: PASS")
    else:
        print("Overall Rule of Five screen: FAIL")

    print("\n======================================")
        # ----------------------------------
    # 4. Molecular visualization
    # ----------------------------------

    print("\n🖼️ Generating molecular structure...")

    image_path = generate_molecule_image(
        smiles,
        drug_name
    )

    print(f"✅ Molecular structure saved to: {image_path}")
    print("✅ PharmaLens analysis completed.")
    print("======================================\n")


def main():

    print("\n======================================")
    print("          PHARMALENS AI")
    print("======================================")

    print("\nAvailable drugs:\n")

    drugs = list_drugs()

    for index, drug_id in enumerate(drugs, start=1):
        drug = get_drug(drug_id)
        print(f"{index}. {drug['name']}")

    print("\nEnter the drug ID from the database.")
    print("Example: aspirin\n")

    selected_drug = input("Drug: ").strip().lower()

    drug = get_drug(selected_drug)

    if drug is None:
        print("\n❌ Drug not found in PharmaLens database.")
        return

    print(f"\n✅ Selected: {drug['name']}")

    run_pharmalens(
        drug["name"],
        drug["smiles"]
    )


if __name__ == "__main__":
    main()