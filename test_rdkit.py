from rdkit import Chem
from rdkit.Chem import Descriptors

# Example molecule: Aspirin
smiles = "CC(=O)OC1=CC=CC=C1C(=O)O"

# Convert SMILES into a molecular structure
molecule = Chem.MolFromSmiles(smiles)

if molecule is None:
    print("❌ Invalid molecule")
else:
    print("✅ Molecule successfully loaded!")
    print("Molecular formula:", Chem.rdMolDescriptors.CalcMolFormula(molecule))
    print("Molecular weight:", round(Descriptors.MolWt(molecule), 2))
    print("Number of atoms:", molecule.GetNumAtoms())
    print("Number of bonds:", molecule.GetNumBonds())