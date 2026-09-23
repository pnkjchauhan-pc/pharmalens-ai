from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski


def calculate_drug_likeness(smiles):
    """
    Evaluate basic drug-likeness properties using
    Lipinski Rule of Five criteria.
    """

    molecule = Chem.MolFromSmiles(smiles)

    if molecule is None:
        return {"error": "Invalid SMILES string"}

    molecular_weight = Descriptors.MolWt(molecule)
    log_p = Descriptors.MolLogP(molecule)
    h_bond_donors = Lipinski.NumHDonors(molecule)
    h_bond_acceptors = Lipinski.NumHAcceptors(molecule)

    checks = {
        "molecular_weight": {
            "value": round(molecular_weight, 2),
            "limit": 500,
            "passes": molecular_weight <= 500
        },
        "logP": {
            "value": round(log_p, 2),
            "limit": 5,
            "passes": log_p <= 5
        },
        "h_bond_donors": {
            "value": h_bond_donors,
            "limit": 5,
            "passes": h_bond_donors <= 5
        },
        "h_bond_acceptors": {
            "value": h_bond_acceptors,
            "limit": 10,
            "passes": h_bond_acceptors <= 10
        }
    }

    violations = sum(
        1 for check in checks.values()
        if not check["passes"]
    )

    return {
        "checks": checks,
        "violations": violations,
        "rule_of_five_pass": violations <= 1
    }


if __name__ == "__main__":

    aspirin_smiles = "CC(=O)OC1=CC=CC=C1C(=O)O"

    result = calculate_drug_likeness(aspirin_smiles)

    print("\n=== PharmaLens AI | Drug-Likeness Analysis ===\n")

    for property_name, details in result["checks"].items():
        status = "PASS" if details["passes"] else "FAIL"

        print(
            f"{property_name}: "
            f"{details['value']} "
            f"→ {status}"
        )

    print(f"\nRule of Five violations: {result['violations']}")

    if result["rule_of_five_pass"]:
        print("Overall Rule of Five screen: PASS")
    else:
        print("Overall Rule of Five screen: FAIL")