import csv
from pathlib import Path

DATA_DIR = Path("data")
DATA_DIR.mkdir(exist_ok=True)

compounds = [
("Aspirin", "CC(=O)OC1=CC=CC=C1C(=O)O", "Analgesic"),
("Paracetamol", "CC(=O)NC1=CC=C(O)C=C1", "Analgesic"),
("Ibuprofen", "CC(C)CC1=CC=C(C=C1)[C@@H](C)C(=O)O", "Analgesic"),
("Naproxen", "COC1=CC2=C(C=C1)C=C(C=C2)[C@@H](C)C(=O)O", "Analgesic"),
("Diclofenac", "O=C(O)C1=C(NC2=CC=CC=C2Cl)C=C(Cl)C=C1", "Analgesic"),
("Ketoprofen", "CC(C(=O)C1=CC=CC=C1)C(=O)O", "Analgesic"),
("Idomethacin", "COC1=C(OC)C=C(C=C1)[C@H]1C2=CC=CC=C2NC1=O", "Analgesic"),


("Caffeine", "CN1C=NC2=C1C(=O)N(C(=O)N2C)C", "CNS"),
("Theophylline", "CN1C=NC2=C1C(=O)NC(=O)N2C", "CNS"),
("Theobromine", "CN1C=NC2=C1C(=O)NC(=O)N2", "CNS"),

("Adenine", "NC1=NC=NC2=C1N=CN2", "Biological"),
("Guanine", "NC1=NC2=C(N1)C(=O)NC(=N2)N", "Biological"),
("Cytosine", "NC1=NC=CC(=O)N1", "Biological"),
("Thymine", "CC1=CN(C(=O)NC1=O)C", "Biological"),
("Uracil", "O=C1NC=CC(=O)N1", "Biological"),
("Histamine", "NCCN1C=NC=C1", "Biological"),

("Ampicillin", "CC1(C(N2C(S1)C(C2=O)N)C(=O)O)N", "Antibiotic"),
("Penicillin G", "CC1(C(N2C(S1)C(C2=O)N)C(=O)O)C(=O)N1C2=CC=CC=C2", "Antibiotic"),
("Ciprofloxacin", "O=C(O)C1=CN(C2CC2)C(=O)C2=C(N1)C=C(F)C=C2N", "Antibiotic"),
("Metronidazole", "CC1=NC=C(N1CCO)C[N+](=O)[O-]", "Antibiotic"),
("Doxycycline", "CN(C)C1C(O)C2C(C(=O)C3=C(O)C=CC=C3O)=C(O)C(=O)C2(O)C1O", "Antibiotic"),

("Atorvastatin", "CC(C)C1=CC=C(C=C1)C(C)C(=O)NC(C)C(O)C(O)CC", "Cardiovascular"),
("Simvastatin", "CC(C)C1CCC(C)C1C(=O)OC(C)C", "Cardiovascular"),
("Warfarin", "CC(C(=O)C1=CC=CC=C1)C2=C(O)C=CC=C2C(=O)O", "Cardiovascular"),
("Amlodipine", "CCOC(=O)C1=C(NC(=C1)C(=O)OC)C2=CC=CC=C2Cl", "Cardiovascular"),
("Losartan", "CCCCC1=NC(C2=CC=CC=C2Cl)=NO1", "Cardiovascular"),
("Atenolol", "CC(C)NCC(COC1=CC=CC=C1)O", "Cardiovascular"),
("Propranolol", "CC(C)NCC(COC1=CC=CC=C1)O", "Cardiovascular"),
("Metoprolol", "COCCOC1=CC=C(C=C1)CC(C)N", "Cardiovascular"),

("Loratadine", "CCOC(=O)N1CCC[C@H]1C2=CC=C(C=C2)Cl", "Antihistamine"),
("Cetirizine", "O=C(O)COCCN1CCC(C1)C2=CC=C(C=C2)Cl", "Antihistamine"),
("Diphenhydramine", "CN(C)CCOC(C1=CC=CC=C1)C2=CC=CC=C2", "Antihistamine"),
("Chlorpheniramine", "CN(C)CCC(C1=CC=CC=C1)C2=CC=CC=C2Cl", "Antihistamine"),

("Omeprazole", "COC1=NC=NC2=C1N=CN2C3=CC=C(C=C3)S(=O)CC", "GI"),
("Pantoprazole", "COC1=NC=NC2=C1N=CN2C3=CC=C(C=C3)S(=O)C4=CC=NC=C4", "GI"),
("Lansoprazole", "CC1=NC=NC2=C1N=CN2C3=CC=C(C=C3)S(=O)CC", "GI"),

("Metformin", "CN(C)C(=N)NC(=N)N", "Antidiabetic"),
("Glipizide", "CC1=CN(C(=O)N1)C2=CC=C(C=C2)S(=O)NCC3=CC=CC=C3", "Antidiabetic"),
("Pioglitazone", "CC1=CC=C(C=C1)CC2=NC(=CS2)C3=CC=CC=C3", "Antidiabetic"),

("Fluoxetine", "CC(C)NCCC(C1=CC=C(C=C1)F)OC2=CC=CC=C2", "Neuroactive"),
("Sertraline", "CN(C)C1CCC(C2=CC=C(C=C2)Cl)C1", "Neuroactive"),
("Diazepam", "CN1C(=O)CN=C(C2=CC=CC=C2Cl)C2=CC=CC=C12", "Neuroactive"),
("Alprazolam", "CN1C=NC2=C(C1=O)C(=C(C=N2)Cl)C3=CC=CC=C3", "Neuroactive"),
("Nicotine", "CN1CCC[C@H]1C2=CN=CC=C2", "Neuroactive"),
("Cotinine", "CN1CCC[C@H]1C2=CN=CC(=O)C2", "Neuroactive"),
("Dopamine", "C1=CC(=C(C=C1CCN)O)O", "Neuroactive"),
("Epinephrine", "CNCC(C1=CC(=C(C=C1)O)O)O", "Neuroactive"),
("Serotonin", "C1=CC2=C(C=C1O)NC=C2CCN", "Neuroactive"),
("Melatonin", "COC1=CC2=C(C=C1)NC(=O)CC2", "Neuroactive"),

("Vitamin C", "C(C(C1C(=O)C(C(=O)O1)O)O)O", "Nutrient"),
("Vitamin B3", "CC1=CC=CC=C1C(=O)O", "Nutrient"),
("Vitamin B6", "CC1=NC=C(CO)C(O)=C1", "Nutrient"),
("Folic Acid", "NC1=NC=CC(=N1)C2=CC(=C(C=C2)C(=O)NCC3=CN=C(N=C3N)N)N", "Nutrient"),
("Riboflavin", "CC1=C(C)C2=C(NC3=C2C(=O)C(=O)N3)N(C)C1=O", "Nutrient"),

("Resveratrol", "C1=CC(=CC=C1C2=CC=C(C=C2)O)O", "Polyphenol"),
("Curcumin", "COC1=CC(=CC(=C1O)OC)C=CC(=O)CC(=O)C=CC2=CC(=C(C=C2)O)OC", "Polyphenol"),
("Quercetin", "C1=CC(=C(C=C1C2=C(C(=O)C3=C(C=C(C=C3O2)O)O)O)O)O", "Polyphenol"),


]

output_file = DATA_DIR / "compounds.csv"

rows = [["Name", "SMILES", "Category"]] + compounds

file = open(output_file, "w", newline="", encoding="utf-8")
writer = csv.writer(file)
writer.writerows(rows)
file.close()

print(f"Database created successfully: {output_file}")
print(f"Total compounds: {len(compounds)}")
print("CSV generation completed successfully.")

