from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski
from rdkit.Chem.rdMolDescriptors import CalcMolFormula


def analyze_molecule(smiles):
    """
    Analyze a molecule from its SMILES representation.
    Returns important molecular and drug-likeness properties.
    """

    molecule = Chem.MolFromSmiles(smiles)

    if molecule is None:
        return {"error": "Invalid SMILES string"}

    properties = {
        "molecular_formula": CalcMolFormula(molecule),
        "molecular_weight": round(Descriptors.MolWt(molecule), 2),
        "logP": round(Descriptors.MolLogP(molecule), 2),
        "tpsa": round(Descriptors.TPSA(molecule), 2),
        "h_bond_donors": Lipinski.NumHDonors(molecule),
        "h_bond_acceptors": Lipinski.NumHAcceptors(molecule),
        "rotatable_bonds": Lipinski.NumRotatableBonds(molecule),
        "heavy_atoms": molecule.GetNumHeavyAtoms(),
        "rings": Lipinski.RingCount(molecule),
    }

    return properties


if __name__ == "__main__":
    aspirin_smiles = "CC(=O)OC1=CC=CC=C1C(=O)O"

    result = analyze_molecule(aspirin_smiles)

    print("\n=== PharmaLens AI | Molecular Analysis ===\n")

    for property_name, value in result.items():
        print(f"{property_name}: {value}")