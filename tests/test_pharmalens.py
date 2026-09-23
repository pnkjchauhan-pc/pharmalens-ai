from src.drug_input import validate_smiles
from src.molecule_analyzer import analyze_molecule
from src.drug_likeness import calculate_drug_likeness
from src.drug_database import get_drug, list_drugs
from src.molecular_fingerprint import generate_fingerprint
from src.molecular_similarity import calculate_similarity


def test_valid_smiles():
    """Test that a valid SMILES is accepted."""

    smiles = "CC(=O)OC1=CC=CC=C1C(=O)O"

    assert validate_smiles(smiles) is True


def test_invalid_smiles():
    """Test that an invalid SMILES is rejected."""

    smiles = "THIS_IS_NOT_A_MOLECULE"

    assert validate_smiles(smiles) is False


def test_molecule_analysis():
    """Test molecular property calculation."""

    aspirin = get_drug("aspirin")

    result = analyze_molecule(aspirin["smiles"])

    assert result["molecular_formula"] == "C9H8O4"
    assert result["molecular_weight"] == 180.16


def test_drug_likeness():
    """Test Lipinski-based drug-likeness analysis."""

    aspirin = get_drug("aspirin")

    result = calculate_drug_likeness(aspirin["smiles"])

    assert result["violations"] == 0
    assert result["rule_of_five_pass"] is True


def test_drug_database():
    """Test that the drug database contains expected drugs."""

    drugs = list_drugs()

    assert "aspirin" in drugs
    assert "caffeine" in drugs
    assert "paracetamol" in drugs


def test_molecular_fingerprint():
    """Test molecular fingerprint generation."""

    aspirin = get_drug("aspirin")

    fingerprint = generate_fingerprint(
        aspirin["smiles"]
    )

    assert fingerprint.GetNumBits() == 2048
    assert fingerprint.GetNumOnBits() > 0


def test_molecular_similarity():
    """Test molecular similarity calculation."""

    aspirin = get_drug("aspirin")
    ibuprofen = get_drug("ibuprofen")

    similarity = calculate_similarity(
        aspirin["smiles"],
        ibuprofen["smiles"]
    )

    assert 0 <= similarity <= 1
    assert similarity > 0