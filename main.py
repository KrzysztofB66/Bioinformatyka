import dna

#1 Test jednostkowy
print(dna.transcription_dna_rna('TACTAGAGCATT'))

#2
seq1 = "TTAAGTGTAGCCTTGTGTGACATGTATTTTTAT"
seq2 = "TTTCTAGGTAGTTGTGGTGAGTTTAGTTGATAT"

compare = dna.compare_dna_sequences(seq1, seq2)
print(f"zasady o identycznym położeniu: {compare[0]}")
print(f"procent zgodności obu sekwencji DNA: {compare[1]}%")
print(f"połączenie dwóch sekwencji: {compare[2]}")
#3
codon = "CUU"
print(f"Amino Acid: {dna.determine_amino_acid(codon)}")

#4
print(dna.translate_mRNA("CUUCUUCUU"))

#5
seqDNA = "ATCGAATGGCGCAAAACCTTTCGCGGTATGGCATGATAGCGCCCGGAAGAGAGTCAATTCAGG"

#Wyznacz ilość par zasad w łańcuchu
dna_length = len(seqDNA)
print(f"Ilość par zasad: {dna_length}")

#Wyznacz nazwy zasad: pierwszą, środkową oraz ostatnią
first_base = seqDNA[0]
if len(seqDNA) % 2 == 0:
    middle_base = seqDNA[dna_length // 2 - 1], seqDNA[dna_length // 2 + 1]

else:
    middle_base = seqDNA[dna_length // 2]

last_base = seqDNA[-1]

print(f"Pierwsza zasada: {first_base}")
print(f"Środkowa zasada: {middle_base}")
print(f"Ostatnia zasada: {last_base}")

#Oblicz liczebność poszczególnych zasad w łańcuchu DNA


A, T, C, G = 0, 0, 0, 0
for base in seqDNA:
    if base == 'A':
        A += 1
    elif base == 'T':
        T += 1
    elif base == 'C':
        C += 1
    elif base == 'G':
        G += 1
    else:
        break
        print("Błędna baza")

print(f"A:{A}, T:{T}, C:{C}, G:{G}")

#Oblicz zawartość % adeniny (A) i cytozyny (C)


print(f"Zawartość A: {round((A / dna_length )*100, 2)}%, C: {round((C / dna_length)*100, 2)}%")


dna_no_g = seqDNA.replace("G", "")
print(f"DNA po usunięciu G: {dna_no_g}")

#Odwróć listę zasad DNA (rewers znaków)
print(f"Odwrócona sekwencja DNA: {seqDNA[::-1]}")

#Wykonaj transkrypcję na mRNA
mRNA, valid = dna.transcription_dna_rna(seqDNA)
print(f"Transkrypcja mRNA: {mRNA}")

#Odwróć sekwencję mRNA

print(f"Odwrócona sekwencja mRNA: {mRNA[::-1]}")

#Wykonaj translację mRNA
print(f"Translacja mRNA: {dna.translate_mRNA(mRNA)}")


#7

Aminokwasy = dna.Aminokwasy()
fasta = dna.wczytaj_fasta('seqProtein1.fasta')
print(Aminokwasy.oblicz_mase_bialka_z_aminokwasow(fasta))
print(Aminokwasy.oblicz_mase_bialka_z_dna('TTAAGTGTATTATTA'))


