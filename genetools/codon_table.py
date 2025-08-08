import json

# class to load our codon table and compare codons
class CodonTable:

    def __init__(self,file_path="genetools/codon_table.json"):

        try:
            with open(file_path,"r") as file:
                self.codon_dict = json.load(file)

        except FileNotFoundError:
            print(f"Error: {file_path} not found.")
            self.codon_dict = {}

        except json.JSONDecodeError:
            print("Error: Invalid JSON format.")
            self.codon_dict = {}

    # returns a list of possible codons for a given amino acid 1 letter abbreviation
    def get_codons(self, amino_acid):
        return self.codon_dict.get(amino_acid,[])
    
    # finds the first instance of a codon in a given nucleotide sequence
    def find_codon(self, amino_acid, nuc_seq):

        possible_codons = self.get_codons(amino_acid)

        for idx,_ in enumerate(nuc_seq):

            # clear codon
            codon = None

            # take this nucleotide and the two after as a codon
            if idx+3 <= len(nuc_seq):
                codon = nuc_seq[idx:idx+3]

            # if the codon is in possible codons then return the index position of the start of the next codon
            if codon in possible_codons:
                return idx+3
            
            # if you reach the end of the sequence without finding codon return None
            elif codon == None:
                return None
            