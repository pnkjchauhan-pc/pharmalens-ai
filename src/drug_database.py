import json
from pathlib import Path


# Find the project root directory
PROJECT_ROOT = Path(__file__).resolve().parent.parent

# Location of our drug database
DATABASE_PATH = PROJECT_ROOT / "data" / "drugs.json"


def load_drug_database():
    """
    Load the drug database from drugs.json.
    """

    with open(DATABASE_PATH, "r", encoding="utf-8") as file:
        database = json.load(file)

    return database


def get_drug(drug_id):
    """
    Retrieve a drug from the database using its ID.
    """

    database = load_drug_database()

    drug = database.get(drug_id.lower())

    if drug is None:
        return None

    return drug


def list_drugs():
    """
    Return all available drug IDs.
    """

    database = load_drug_database()

    return list(database.keys())


if __name__ == "__main__":

    print("\n=== PharmaLens AI | Drug Database ===\n")

    drugs = list_drugs()

    print("Available drugs:")

    for drug in drugs:
        print(f"- {drug}")

    print("\nTesting database lookup...")

    aspirin = get_drug("aspirin")

    if aspirin:
        print("\n✅ Aspirin found!")
        print(f"Name: {aspirin['name']}")
        print(f"SMILES: {aspirin['smiles']}")
    else:
        print("\n❌ Aspirin not found.")