import unittest
from dna import Aminokwasy, transcription_dna_rna, translate_mRNA  # Zmień 'twoj_modul' na nazwę pliku


class TestAminokwasy(unittest.TestCase):

    def setUp(self):
        self.amino = Aminokwasy()

    def test_oblicz_mase_bialka_z_aminokwasow(self):
        self.assertEqual(self.amino.oblicz_mase_bialka_z_aminokwasow("A"), 89)
        self.assertEqual(self.amino.oblicz_mase_bialka_z_aminokwasow("AV"), 89 + 117)
        self.assertEqual(self.amino.oblicz_mase_bialka_z_aminokwasow("AVF"), 89 + 117 + 165)
        self.assertEqual(self.amino.oblicz_mase_bialka_z_aminokwasow(""), 0)
        self.assertEqual(self.amino.oblicz_mase_bialka_z_aminokwasow("AX"),
                         "Błąd: Nieprawidłowy aminokwas w sekwencji.")

    def test_oblicz_mase_bialka_z_dna(self):
        # Test dla poprawnej sekwencji DNA
        dna_seq = "TACTACTACTACTAC"  # Powinno być przekształcone na "M" (Metionina)
        self.assertEqual(self.amino.oblicz_mase_bialka_z_dna(dna_seq), 745)  # Masa Metioniny (M)

        # Test dla niepoprawnej sekwencji DNA (niepoprawne znaki)
        dna_invalid = "ATGXYZ"
        self.assertEqual(self.amino.oblicz_mase_bialka_z_dna(dna_invalid),
                         "Błąd z sekwencją dna.")

        # Test dla zbyt krótkiej sekwencji DNA
        dna_short = "ATG"
        self.assertEqual(self.amino.oblicz_mase_bialka_z_dna(dna_short),
                         "Błąd z sekwencją dna.")


if __name__ == "__main__":
    unittest.main()
