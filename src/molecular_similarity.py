from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator


def calculate_similarity(smiles_1, smiles_2):
    """
    Calculate Tanimoto similarity between two molecules.
    """

    molecule_1 = Chem.MolFromSmiles(smiles_1)
    molecule_2 = Chem.MolFromSmiles(smiles_2)

    if molecule_1 is None:
        raise ValueError("First SMILES string is invalid.")

    if molecule_2 is None:
        raise ValueError("Second SMILES string is invalid.")

    fingerprint_generator = rdFingerprintGenerator.GetMorganGenerator(
        radius=2,
        fpSize=2048
    )

    fingerprint_1 = fingerprint_generator.GetFingerprint(molecule_1)
    fingerprint_2 = fingerprint_generator.GetFingerprint(molecule_2)

    similarity = DataStructs.TanimotoSimilarity(
        fingerprint_1,
        fingerprint_2
    )

    return round(similarity, 4)


if __name__ == "__main__":

    aspirin = "CC(=O)OC1=CC=CC=C1C(=O)O"
    ibuprofen = "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O"

    similarity = calculate_similarity(
        aspirin,
        ibuprofen
    )

    print("\n=== PharmaLens AI | Molecular Similarity ===\n")

    print("Molecule 1: Aspirin")
    print("Molecule 2: Ibuprofen")
    print(f"Tanimoto similarity: {similarity}")

    print("\n✅ Molecular similarity calculated successfully.")