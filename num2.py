# Define the input amino acid sequence
amino_acid_sequence = "NAN"

# Create a dictionary to map codons to their corresponding amino acids
codon_to_amino_acid = {
    "UUU": "F", "UUC": "F", "UUA": "L", "UUG": "L",
    "CUU": "L", "CUC": "L", "CUA": "L", "CUG": "L",
    "AUU": "I", "AUC": "I", "AUA": "I", "AUG": "M",
    "GUU": "V", "GUC": "V", "GUA": "V", "GUG": "V",
    "UCU": "S", "UCC": "S", "UCA": "S", "UCG": "S",
    "CCU": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACU": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCU": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "UAU": "Y", "UAC": "Y", "UAA": "STOP", "UAG": "STOP",
    "CAU": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAU": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAU": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "UGU": "C", "UGC": "C", "UGA": "STOP", "UGG": "W",
    "CGU": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGU": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGU": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

# Function to convert an amino acid sequence to mRNA
def amino_acid_to_mrna(amino_acid):
    mrna = ""
    for aa in amino_acid:
        for codon, aa_code in codon_to_amino_acid.items():
            if aa == aa_code:
                mrna += codon
                break
    return mrna

# Count the occurrences of specific codons
def count_codons(mrna, codon):
    count = 0
    start = 0
    while start < len(mrna):
        index = mrna.find(codon, start)
        if index == -1:
            break
        count += 1
        start = index + 1
    return count

# Convert amino acid sequence to mRNA
mrna_sequence = amino_acid_to_mrna(amino_acid_sequence)

# Count the occurrences of specific codons
codon_AAU_count = count_codons(mrna_sequence, "AAU")
codon_GCU_count = count_codons(mrna_sequence, "GCU")

# Print the results
print(f"Input Aminoacid = {amino_acid_sequence}")
print(f"mRNA = {mrna_sequence}")
print(f"AAU = {codon_AAU_count}")
print(f"GCU = {codon_GCU_count}")


