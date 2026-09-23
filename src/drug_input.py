from rdkit import Chem


def validate_smiles(smiles):
    """
    Check whether the provided SMILES string
    represents a valid molecule.
    """

    molecule = Chem.MolFromSmiles(smiles)

    if molecule is None:
        return False

    return True


def get_molecule_from_smiles(smiles):
    """
    Convert a valid SMILES string into an RDKit molecule.
    """

    molecule = Chem.MolFromSmiles(smiles)

    if molecule is None:
        raise ValueError("Invalid SMILES string.")

    return molecule


if __name__ == "__main__":

    test_smiles = {
        "Aspirin": "CC(=O)OC1=CC=CC=C1C(=O)O",
        "Caffeine": "Cn1c(=O)c2c(ncn2C)n(C)c1=O",
        "Paracetamol": "CC(=O)NC1=CC=C(C=C1)O"
    }

    print("\n=== PharmaLens AI | Drug Input Validation ===\n")

    for drug_name, smiles in test_smiles.items():

        if validate_smiles(smiles):
            print(f"✅ {drug_name}: Valid molecule")
        else:
            print(f"❌ {drug_name}: Invalid molecule")