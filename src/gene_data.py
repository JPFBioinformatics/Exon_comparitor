import requests
import pandas as pd
from codon_table import CodonTable
from Bio.Align import PairwiseAligner, substitution_matrices

class GeneData:

    def __init__(self, gene_name = None, species = 'Human'):
        self.BASE_URL = "https://rest.ensembl.org"
        self.gene_name = gene_name
        self.species = species
        self.aa_seq = None
        self.gene = None
        self.transcript_id = None                   # the most recent canonical transcript (where the external_name is gene_name-201)
        self.dna_sequence = None    
        self.cdna_sequence = None
        self.exons = []
        self.gene_start = None                      # the starting position of the gene on the genome
        self._fetch_gene_info()

    # Retreives gene_id, gene_name, and best transcript for an Ensembl gene search
    def _fetch_gene_info(self):

        url = f"{self.BASE_URL}/lookup/symbol/{self.species}/{self.gene_name}?content-type=application/json"
        response = requests.get(url)
        
        if response.status_code == 200:

            g_data = response.json()

            # store gene info
            self.gene = g_data

            primary_transcript = g_data.get('canonical_transcript')
            self.transcript_id = primary_transcript.split(".")[0]

        return None

    # gets the full DNA sequence, with both introns and exons
    def get_dna_sequence(self):
        
        url = f"{self.BASE_URL}/sequence/id/{self.gene['id']}?content-type=application/json"

        response = requests.get(url)

        if response.status_code == 200:
            sequence_info = response.json()
        else:
            print(f"{response.status_code}")

        if sequence_info.get('seq') is not None:
            self.dna_sequence = sequence_info.get('seq')
        else:
            print("No sequence data found")

    # gets the CDNA for each exon
    def get_cdna_sequence(self):
        
        url = f"{self.BASE_URL}/sequence/id/{self.transcript_id}?content-type=application/json;type=cdna"

        response = requests.get(url)

        if response.status_code == 200:
            sequence_info = response.json()
        else:
            print(f"{response.status_code}")

        if sequence_info.get('seq') is not None:
            self.cdna_sequence = sequence_info.get('seq')
        else:
            print("No sequence data found")

        return None

    # gets a list of exons for the gene in order of appearance in the CDNA
    def get_exons(self):

        url = f"{self.BASE_URL}/overlap/id/{self.transcript_id}/?feature=exon;content-type=application/json;type=cdna"
        response = requests.get(url)

        if response.status_code == 200:
            all_exons = response.json()

        for exon in all_exons:
            if exon.get('Parent') == self.transcript_id:
                self.exons.append(exon)

        self.gene_start = self.exons[0].get('start')

        exons = self.exons

        for exon in exons:
            
            exon_url = f"{self.BASE_URL}/sequence/id/{exon['id']}?content-type=application/json;type=cdna"
            response = requests.get(exon_url)

            if response.status_code == 200:
                data = response.json()
                exon['cdna_seq'] = data['seq']

            else:
                print("No exon cdna found")

    # gets the transcript for each exon for you gene
    def get_aa_sequence(self):
           
        url = f"{self.BASE_URL}/sequence/id/{self.transcript_id}?content-type=application/json;type=protein"
        response = requests.get(url)

        if response.status_code == 200:

            data = response.json()
            print(data)
            self.aa_seqeuence = data['seq']

        else:
            print(f"No aa sequence found for transcript {self.transcript_id}")

        return None

    # takes aa_sequence and the exon list and finds what group of aa a given exon codes
    def aa_map_exons(self):

        # create a codon table for comparison
        codon_table = CodonTable()

        # get list of exons for mapping
        exons = self.exons
        # find the sequence of amino acids for mapping
        seq = self.aa_seq
        # start a counter for where in the sequence we are
        seq_idx = 0

        for exon in exons:

            exon_aas = []
            exon_seq = exon['cdna_seq']
            start_idx = 0

            for _ in exon_seq:
                
                # take the section of the exon cdna for query
                section = exon_seq[start_idx:]

                # locate the rseidue we are searching for
                residue = seq[seq_idx]

                # find the first instance of a codon for that residue in the query section
                next_idx = codon_table.find_codon(residue,section)

                # if no matching codon found in exon print error message
                if next_idx == None:
                    print(f"Could not find codon for amino acid {residue} in exon {exon['rank']}")
                    return None
                
                else:
                    # increment sequence index and start index
                    seq_idx += 1
                    start_idx = next_idx
                    # add residue to list of resiudes for this amnino acid
                    exon_aas.append(residue)

            # save the exon amino acid sequence to exon dictionary
            exon['aa_seq'] = ''.join(exon_aas)

    # preforms a Smith_Waterman local alignment of two given exon amino acid sequences
    @staticmethod
    def algin_exon_aa(exon1, exon2, sub_matrix = "BLOSUM62", score_mode = "blastp"):
        
        # load the BLOSUM62 substituion matrix, other options are 45,50,80,90 corrosponding to expected % similarity in sequence
        blosum = substitution_matrices.load(sub_matrix)

        # create new PairwiseAligner object with default parameters and score_mode scoring
        aligner = PairwiseAligner(scoring=score_mode)

        # assign substituion matrix to aligner
        aligner.substitution_matrix = blosum

        # set mode to global, will align the entire sequence as best it can, swich it to local to align based on most similar "chunk" of exons
        aligner.mode = "global"

        # align input exons
        alginments = aligner.align(exon1['aa_seq'],exon2['aa_seq'])
        
        for alignment in alginments:
            print(alignment)
            print("")

    # preforms a Smith_Waterman local alignment of two given exon amino acid sequences
    @staticmethod
    def align_exon_cdna(exon1, exon2, score_mode = "blastn"):

        # create new PairwiseAligner object with default parameters and score_mode scoring
        aligner = PairwiseAligner(scoring=score_mode)

        # set mode to global, will align the entire sequence as best it can, swich it to local to align based on most similar "chunk" of exons
        aligner.mode = "global"

        # align input exons
        alignments = aligner.align(exon1['cdna_seq'],exon2['cdna_seq'])

        for alignment in alignments:
            print(alignment)
            print("")