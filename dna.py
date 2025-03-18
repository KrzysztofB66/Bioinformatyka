import os


#6
class Aminokwasy:
    """
    Klasa do obliczania masy cząsteczkowej białek na podstawie sekwencji aminokwasów.
    """
    masa_aminokwasow = {'A':89,'V':117,'L':131,'I':131,'P':115,'F':165,'W':204,'M':149,'G':75,'S':105,'C':121,'T':119,
                        'Y':181,'N':132,'Q':146,'D':133,'E':147,'K':146,'R':174,'H':155}

    def __init__(self):
        pass

    def oblicz_mase_bialka_z_aminokwasow(self, sekwencja_aminokwasow):
        """
        Oblicza masę białka na podstawie sekwencji aminokwasów.

        Args:
        sekwencja_aminokwasow (str): Sekwencja aminokwasów.

        Returns:
            float: Masa cząsteczkowa białka.
        """

        masa = 0
        for aminokwas in sekwencja_aminokwasow.upper():
            if aminokwas in self.masa_aminokwasow:
                masa += self.masa_aminokwasow[aminokwas]
            else:
                return f"Błąd: Nieprawidłowy aminokwas w sekwencji."

        return round(masa, 2)

    def oblicz_mase_bialka_z_dna(self, dna_sequence):
        """
        Translacja sekwencji DNA na białko i oblicza jego masę.

        Args:
            dna_sequence (str): Sekwencja DNA.

        Returns:
            float: Masa białka.
            str: Komunikat o błędzie w przypadku problemu przy translacji mRNA na białko.
        """
        mRNA, valid = transcription_dna_rna(dna_sequence)
        if not valid:
            return "Błąd z sekwencją dna."
        else:
            bialko = translate_mRNA(mRNA)

        return self.oblicz_mase_bialka_z_aminokwasow(bialko)

# 1
def validate_dna(dna_sequence):
    """
    Waliduje, czy podana sekwencja DNA ma co najmniej 10 nukleotydów i zawiera tylko prawidłowe znaki (A, T, C, G).

    Args:
        dna_sequence (str): Sekwencja DNA do walidacji.

    Returns:
        tuple: (bool, str), gdzie bool wskazuje, czy sekwencja jest prawidłowa, a str zawiera komunikat o błędzie, jeśli nie jest prawidłowa.
    """
    dna_sequence = dna_sequence.upper()
    if len(dna_sequence) < 10:
        return False, "Liczba nukleotydów < 10"

    if not set(dna_sequence).issubset({'A', 'T', 'C', 'G'}):
        return False, "Sekwencja zawiera nieprawidłowe znaki"
    return True, "Sekwencja jest poprawna"


def transcription_dna_rna(dna_sequence):
    """
    Transkrybuje sekwencję DNA do sekwencji mRNA, zastępując nukleotydy zgodnie z mapą transkrypcyjną.

    Args:
        dna_sequence (str): Sekwencja DNA do transkrypcji.

    Returns:
        str: Transkrybowana sekwencja mRNA.
    """
    dna_sequence = dna_sequence.upper()
    valid, message = validate_dna(dna_sequence)

    # If the sequence is invalid, return the error message
    if not valid:
        return message, valid

    transcription_map = {'A': 'U', 'C': 'G', 'G': 'C', 'T': 'A'}
    rna_sequence = "".join([transcription_map[nuc] for nuc in dna_sequence])

    return rna_sequence , valid



# 2
def find_matching_bases(seq1, seq2):
    """
    Znajduje pozycje pasujących zasad pomiędzy dwoma sekwencjami DNA.

    Args:
        seq1 (str): Pierwsza sekwencja DNA.
        seq2 (str): Druga sekwencja DNA.

    Returns:
        list: Lista krotek, w której zasady pasują.
        """
    matching_bases = []
    for i in range(len(seq1)):
        if seq1[i] == seq2[i]:
            matching_bases.append((i, seq1[i]))
    return matching_bases


def calculate_match_percentage(seq1, seq2):
    """
    Oblicza procent pasujących zasad pomiędzy dwoma sekwencjami.

    Args:
        seq1 (str): Pierwsza sekwencja DNA.
        seq2 (str): Druga sekwencja DNA.

    Returns:
        float: Procent pasujących zasad.
       """
    match_count = 0
    for i in range(len(seq1)):
        if seq1[i] == seq2[i]:
            match_count += 1

    return (match_count / len(seq1)) * 100


def combine_sequences(seq1, seq2):
    """
    Łączy dwie sekwencje DNA w jeden ciąg.

    Args:
        seq1 (str): Pierwsza sekwencja DNA.
        seq2 (str): Druga sekwencja DNA.

    Returns:
        str: Połączony ciąg DNA.
        """

    return f"{seq1}{seq2}"


def compare_dna_sequences(seq1, seq2):
    """
    Porównuje dwie sekwencje DNA, sprawdzając ich poprawność, długość i dopasowanie zasad.

    Args:
        seq1 (str): Pierwsza sekwencja DNA.
        seq2 (str): Druga sekwencja DNA.

    Returns:
        tuple: Krotka zawierająca:
            - listę dopasowanych zasad (pozycja, zasada),
            - procent dopasowania zasad,
            - połączoną sekwencję.
        str: W przypadku błędów, zwraca komunikat o błędzie.
        """
    valid1 = validate_dna(seq1)
    valid2 = validate_dna(seq2)

    if not valid1 or not valid2:
        return "Błąd: Jedna lub obie sekwencje są niepoprawne."

    if len(seq1) != len(seq2):
        return "Błąd: Sekwencje muszą mieć identyczną długość."

    matching_nucleotides = find_matching_bases(seq1, seq2)
    match_percentage = calculate_match_percentage(seq1, seq2)
    combined_sequences = combine_sequences(seq1, seq2)

    return matching_nucleotides, round(match_percentage, 2), combined_sequences

# 3
def determine_amino_acid(codon):
    """
    Określa aminokwas na podstawie trójnukleotydowego kodonu mRNA.

    Args:
        codon (str): Kodon mRNA (ciąg trzech nukleotydów: U, C, A, G).

    Returns:
        str: Kodowany aminokwas, lub 'stop' w przypadku kodonu stop.
        str: Komunikat o błędzie, jeśli kodon jest nieprawidłowy.
    """
    codon_table = {
        'U': {
            'U': {'U': 'F', 'C': 'F', 'A': 'L', 'G': 'L'},
            'C': {'U': 'S', 'C': 'S', 'A': 'S', 'G': 'S'},
            'A': {'U': 'Y', 'C': 'Y', 'A': None, 'G': None},
            'G': {'U': 'C', 'C': 'C', 'A': None, 'G': 'W'}
        },
        'C': {
            'U': {'U': 'L', 'C': 'L', 'A': 'L', 'G': 'L'},
            'C': {'U': 'P', 'C': 'P', 'A': 'P', 'G': 'P'},
            'A': {'U': 'H', 'C': 'H', 'A': 'Q', 'G': 'Q'},
            'G': {'U': 'R', 'C': 'R', 'A': 'R', 'G': 'R'}
        },
        'A': {
            'U': {'U': 'I', 'C': 'I', 'A': 'I', 'G': 'M'},
            'C': {'U': 'T', 'C': 'T', 'A': 'T', 'G': 'T'},
            'A': {'U': 'N', 'C': 'N', 'A': 'K', 'G': 'K'},
            'G': {'U': 'S', 'C': 'S', 'A': 'R', 'G': 'R'}
        },
        'G': {
            'U': {'U': 'V', 'C': 'V', 'A': 'V', 'G': 'V'},
            'C': {'U': 'A', 'C': 'A', 'A': 'A', 'G': 'A'},
            'A': {'U': 'D', 'C': 'D', 'A': 'E', 'G': 'E'},
            'G': {'U': 'G', 'C': 'G', 'A': 'G', 'G': 'G'}
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
        return "stop"
    else:
        return amino_acid

# 4
def translate_mRNA(mRNA):
    """
    Translacja sekwencji mRNA do sekwencji aminokwasów (białka).

    Args:
        mRNA (str): Sekwencja mRNA (ciąg nukleotydów: U, C, A, G).

    Returns:
        str: Przetłumaczona sekwencja białka (ciąg aminokwasów).
        str: Komunikat o błędzie w przypadku nieprawidłowej sekwencji mRNA.
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
            if amino_acid == "stop":
                break
            amino_acids.append(amino_acid)

        return "".join(amino_acids)

    except ValueError as e:
        return f"Błąd: {e}"

#7
def wczytaj_fasta(sciezka):
    """
    Odczytuje plik .FASTA i zwraca sekwencję DNA jako ciąg znaków.

    Args:
        sciezka (str): Ścieżka do pliku.

    Returns:
        str: Sekwencja DNA.

    Raises:
        FileNotFoundError: Jeśli plik nie istnieje.
        ValueError: Jeśli plik jest pusty.
    """

    if not os.path.exists(sciezka):
        raise FileNotFoundError(f"Błąd: Plik {sciezka} nie istnieje.")

    plik = open(sciezka, "r").readlines()

    if not plik or len(plik) < 2:
        raise ValueError(f"Błąd: Plik {sciezka} jest pusty lub nie zawiera sekwencji.")

    sekwencja = ""
    for linia in plik[1:]:
        sekwencja += linia.strip()

    return sekwencja






