from drug_database import load_drug_database
from molecular_similarity import calculate_similarity


def search_similar_drugs(
    target_smiles,
    top_n=3,
    similarity_threshold=0.15
):
    """
    Compare a target molecule against drugs in the database.

    The target molecule itself is excluded from the results.

    Only molecules meeting the similarity threshold
    are returned.
    """

    database = load_drug_database()

    results = []

    for drug_id, drug in database.items():

        similarity = calculate_similarity(
            target_smiles,
            drug["smiles"]
        )

        # Exclude identical molecules
        if similarity >= 1.0:
            continue

        # Apply similarity threshold
        if similarity >= similarity_threshold:

            results.append({
                "drug_id": drug_id,
                "drug_name": drug["name"],
                "similarity": similarity
            })

    results.sort(
        key=lambda item: item["similarity"],
        reverse=True
    )

    return results[:top_n]


if __name__ == "__main__":

    aspirin_smiles = "CC(=O)OC1=CC=CC=C1C(=O)O"

    results = search_similar_drugs(
        aspirin_smiles,
        top_n=3,
        similarity_threshold=0.15
    )

    print("\n=== PharmaLens AI | Similarity Search ===\n")

    print("Target molecule: Aspirin")
    print("Similarity threshold: 0.15\n")

    print("Most structurally similar molecules:\n")

    if not results:
        print("No similar molecules found.")

    else:
        for index, result in enumerate(results, start=1):

            print(
                f"{index}. {result['drug_name']} "
                f"→ Tanimoto similarity: "
                f"{result['similarity']}"
            )

    print("\n✅ Similarity search completed successfully.")