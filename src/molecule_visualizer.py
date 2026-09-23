from pathlib import Path

from rdkit import Chem
from rdkit.Chem import Draw


# Project root directory
PROJECT_ROOT = Path(__file__).resolve().parent.parent

# Folder where molecule images will be saved
REPORTS_DIR = PROJECT_ROOT / "reports"

# Create reports folder if it doesn't exist
REPORTS_DIR.mkdir(exist_ok=True)


def generate_molecule_image(smiles, drug_name):
    """
    Generate a 2D molecular structure image from SMILES.
    """

    molecule = Chem.MolFromSmiles(smiles)

    if molecule is None:
        raise ValueError("Invalid SMILES string.")

    # Create a safe filename
    filename = drug_name.lower().replace(" ", "_") + "_structure.png"

    output_path = REPORTS_DIR / filename

    # Generate and save the molecular image
    image = Draw.MolToImage(
        molecule,
        size=(800, 600)
    )

    image.save(output_path)

    return output_path


if __name__ == "__main__":

    drug_name = "Aspirin"

    smiles = "CC(=O)OC1=CC=CC=C1C(=O)O"

    output = generate_molecule_image(
        smiles,
        drug_name
    )

    print("\n=== PharmaLens AI | Molecular Visualization ===\n")

    print("✅ Molecular structure generated!")
    print(f"Saved to: {output}")