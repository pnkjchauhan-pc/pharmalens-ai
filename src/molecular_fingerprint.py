from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator


def generate_fingerprint(smiles, radius=2, fp_size=2048):
    """
    Generate a Morgan molecular fingerprint from a SMILES string.
    """

    molecule = Chem.MolFromSmiles(smiles)

    if molecule is None:
        raise ValueError("Invalid SMILES string.")

    fingerprint_generator = rdFingerprintGenerator.GetMorganGenerator(
        radius=radius,
        fpSize=fp_size
    )

    fingerprint = fingerprint_generator.GetFingerprint(molecule)

    return fingerprint


if __name__ == "__main__":

    aspirin_smiles = "CC(=O)OC1=CC=CC=C1C(=O)O"

    fingerprint = generate_fingerprint(aspirin_smiles)

    print("\n=== PharmaLens AI | Molecular Fingerprint ===\n")

    print("Drug: Aspirin")
    print(f"Fingerprint size: {fingerprint.GetNumBits()}")
    print(f"Active bits: {fingerprint.GetNumOnBits()}")

    print("\n✅ Molecular fingerprint generated successfully.")