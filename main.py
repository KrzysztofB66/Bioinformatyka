import re

#1
def validate_dna(dna_sequence):
    dna_sequence = dna_sequence.upper()
    if len(dna_sequence) < 10:
        return False, "Liczba nukleotydów < 10"

    if not set(dna_sequence).issubset({'A', 'T', 'C', 'G'}):
        return False, "Sekwencja zawiera nieprawidłowe znaki"
    return True


def transcription_dna_rna(dna_sequence):
    dna_sequence = dna_sequence.upper()
    
    if validate_dna(dna_sequence) == False:
        return False
    transcription_map = {'A': 'U', 'C': 'G', 'G': 'C', 'T': 'A'}
    rna_sequence = ""
    for nuc in dna_sequence:
        rna_sequence += transcription_map[nuc]

    return rna_sequence


print(transcription_dna_rna("TACTAGAGCATT"))


#2
def find_matching_bases(seq1, seq2):
    matching_bases = []
    for i in range(len(seq1)):
        if seq1[i] == seq2[i]:
            matching_bases.append((i, seq1[i]))
    return matching_bases


def calculate_match_percentage(seq1, seq2):
    match_count = 0
    for i in range(len(seq1)):
        if seq1[i] == seq2[i]:
            match_count += 1

    return (match_count / len(seq1)) * 100


def combine_sequences(seq1, seq2):
    """Łączy dwie sekwencje DNA."""
    return f"{seq1} \n {seq2}"


def compare_dna_sequences(seq1, seq2):
    """Główna funkcja zarządzająca porównaniem DNA."""
    valid1= validate_dna(seq1)
    valid2= validate_dna(seq2)

    if not valid1 or not valid2:
        return "Błąd: Jedna lub obie sekwencje są niepoprawne."

    if len(seq1) != len(seq2):
        return "Błąd: Sekwencje muszą mieć identyczną długość."

    matching_nucleotides = find_matching_bases(seq1, seq2)
    match_percentage = calculate_match_percentage(seq1, seq2)
    combined_sequences = combine_sequences(seq1, seq2)
    print(f"matching_nucleotides: {matching_nucleotides}")
    print(f"match_percentage: {round(match_percentage, 2)}%")
    print(f"combined_sequences: \n {combined_sequences}")


seq1 = "TTAAGTGTAGCCTTGTGTGACATGTATTTTTAT"
seq2 = "TTTCTAGGTAGTTGTGGTGAGTTTAGTTGATAT"

compare_dna_sequences(seq1, seq2)


#3
def determine_amino_acid(codon):

    codon_table = {
        'U': {
            'U': {'U': 'Phe (F)', 'C': 'Phe (F)', 'A': 'Leu (L)', 'G': 'Leu (L)'},
            'C': {'U': 'Ser (S)', 'C': 'Ser (S)', 'A': 'Ser (S)', 'G': 'Ser (S)'},
            'A': {'U': 'Tyr (Y)', 'C': 'Tyr (Y)', 'A': None, 'G': None},
            'G': {'U': 'Cys (C)', 'C': 'Cys (C)', 'A': None, 'G': 'Trp (W)'}
        },
        'C': {
            'U': {'U': 'Leu (L)', 'C': 'Leu (L)', 'A': 'Leu (L)', 'G': 'Leu (L)'},
            'C': {'U': 'Pro (P)', 'C': 'Pro (P)', 'A': 'Pro (P)', 'G': 'Pro (P)'},
            'A': {'U': 'His (H)', 'C': 'His (H)', 'A': 'Gln (Q)', 'G': 'Gln (Q)'},
            'G': {'U': 'Arg (R)', 'C': 'Arg (R)', 'A': 'Arg (R)', 'G': 'Arg (R)'}
        },
        'A': {
            'U': {'U': 'Ile (I)', 'C': 'Ile (I)', 'A': 'Ile (I)', 'G': 'Met (M)'},
            'C': {'U': 'Thr (T)', 'C': 'Thr (T)', 'A': 'Thr (T)', 'G': 'Thr (T)'},
            'A': {'U': 'Asn (N)', 'C': 'Asn (N)', 'A': 'Lys (K)', 'G': 'Lys (K)'},
            'G': {'U': 'Ser (S)', 'C': 'Ser (S)', 'A': 'Arg (R)', 'G': 'Arg (R)'}
        },
        'G': {
            'U': {'U': 'Val (V)', 'C': 'Val (V)', 'A': 'Val (V)', 'G': 'Val (V)'},
            'C': {'U': 'Ala (A)', 'C': 'Ala (A)', 'A': 'Ala (A)', 'G': 'Ala (A)'},
            'A': {'U': 'Asp (D)', 'C': 'Asp (D)', 'A': 'Glu (E)', 'G': 'Glu (E)'},
            'G': {'U': 'Gly (G)', 'C': 'Gly (G)', 'A': 'Gly (G)', 'G': 'Gly (G)'}
        }
    }

    if len(codon) != 3:
        return "Error: Codon must be 3 nucleotides long."

    codon = codon.upper()

    for nucleotide in codon:
        if nucleotide not in 'UCAG':
            return "Error: Codon contains invalid characters."

    amino_acid = codon_table[codon[0]][codon[1]][codon[2]]

    if amino_acid is None:
        return "STOP Codon"
    else:
        return amino_acid


codon = input("Enter a codon")
amino_acid = determine_amino_acid(codon)
print(f"Amino Acid: {amino_acid}")


# 4
def translate_mRNA(mRNA):
    """
    Tłumaczy sekwencję mRNA na sekwencję aminokwasów.

    Args:
        mRNA (str): Sekwencja mRNA.

    Returns:
        str: Sekwencja aminokwasów lub komunikat o błędzie.
    """

    try:
        mRNA = mRNA.upper()
        if len(mRNA) % 3 != 0:
            raise ValueError("Długość mRNA musi być wielokrotnością 3.")
        if not all(nucleotide in 'UCAG' for nucleotide in mRNA):
            raise ValueError("mRNA zawiera nieprawidłowe nukleotydy.")

        amino_acids = []
        for i in range(0, len(mRNA), 3):
            codon = mRNA[i:i + 3]
            amino_acid = determine_amino_acid(codon)
            if amino_acid == "STOP":
                break  # Przerywa translację po napotkaniu kodonu STOP
            amino_acids.append(amino_acid)

        return "-".join(amino_acids)

    except ValueError as e:
        return f"Błąd: {e}"