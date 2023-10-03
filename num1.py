def complement_DNA(dna_sequence):
    complement_dict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    complement_sequence = ''.join(complement_dict[base] for base in dna_sequence)
    return complement_sequence

def transcribe_to_mRNA(dna_sequence):
    return dna_sequence.replace('T', 'U')

def translate_to_amino_acids(mrna_sequence):
    codon_table = {
        'UUU': 'Phe (F)', 'UUC': 'Phe (F)', 'UUA': 'Leu (L)', 'UUG': 'Leu (L)',
        'CUU': 'Leu (L)', 'CUC': 'Leu (L)', 'CUA': 'Leu (L)', 'CUG': 'Leu (L)',
        'AUU': 'Ile (I)', 'AUC': 'Ile (I)', 'AUA': 'Ile (I)', 'AUG': 'Met (M)',
        'GUU': 'Val (V)', 'GUC': 'Val (V)', 'GUA': 'Val (V)', 'GUG': 'Val (V)',
        'UCU': 'Ser (S)', 'UCC': 'Ser (S)', 'UCA': 'Ser (S)', 'UCG': 'Ser (S)',
        'CCU': 'Pro (P)', 'CCC': 'Pro (P)', 'CCA': 'Pro (P)', 'CCG': 'Pro (P)',
        'ACU': 'Thr (T)', 'ACC': 'Thr (T)', 'ACA': 'Thr (T)', 'ACG': 'Thr (T)',
        'GCU': 'Ala (A)', 'GCC': 'Ala (A)', 'GCA': 'Ala (A)', 'GCG': 'Ala (A)',
        'UAU': 'Tyr (Y)', 'UAC': 'Tyr (Y)', 'UAA': 'STOP (*)', 'UAG': 'STOP (*)',
        'CAU': 'His (H)', 'CAC': 'His (H)', 'CAA': 'Gln (Q)', 'CAG': 'Gln (Q)',
        'AAU': 'Asn (N)', 'AAC': 'Asn (N)', 'AAA': 'Lys (K)', 'AAG': 'Lys (K)',
        'GAU': 'Asp (D)', 'GAC': 'Asp (D)', 'GAA': 'Glu (E)', 'GAG': 'Glu (E)',
        'UGU': 'Cys (C)', 'UGC': 'Cys (C)', 'UGA': 'STOP (*)', 'UGG': 'Trp (W)',
        'CGU': 'Arg (R)', 'CGC': 'Arg (R)', 'CGA': 'Arg (R)', 'CGG': 'Arg (R)',
        'AGU': 'Ser (S)', 'AGC': 'Ser (S)', 'AGA': 'Arg (R)', 'AGG': 'Arg (R)',
        'GGU': 'Gly (G)', 'GGC': 'Gly (G)', 'GGA': 'Gly (G)', 'GGG': 'Gly (G)'
    }

    amino_acids = []
    for i in range(0, len(mrna_sequence), 3):
        codon = mrna_sequence[i:i+3]
        amino_acid = codon_table.get(codon, 'Unknown')
        amino_acids.append(amino_acid)

    return amino_acids

# Input DNA sequence
input_dna = "TTACGA"

# Calculate complement
complement = complement_DNA(input_dna)

# Transcribe to mRNA
mRNA = transcribe_to_mRNA(complement)

# Translate to amino acids
amino_acids = translate_to_amino_acids(mRNA)

# Print the results
print("Input DNA =", input_dna)
print("Complement =", complement)
print("mRNA =", mRNA)
print("Aminoacid =", ' – '.join(amino_acids))
