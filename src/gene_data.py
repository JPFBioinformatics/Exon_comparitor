import requests
import pandas as pd

class GeneData:
    BASE_URL = "https://rest.ensembl.org"

    def __init__(self, gene_name, species = 'Human'):
        self.gene_name = gene_name
        self.species = species
        self.gene_id = None
        self.transscript = None                 # the most recent canonical transcript (index = -1)
        self.dna_sequence = None
        self.cdna_sequence = None
        self.exons = []
        self.gene_start = None                  # the starting position of the gene on the genome
        self._fetch_gene_info()

    def _fetch_gene_info(self):
        # Retreive gene information and store the Ensembl ID

        url = f"{self.BASE_URL}/lookup/symbol/{self.species}/{self.gene_name}?content-type=application/json"

        response = requests.get(url)

        if response.status_code == 200:
            data = response.json()
            self.gene_id = data.get('id')
        else:
            print(f"Error: Gene '{self.gene_name}' not found.")

    def get_transcripts(self):
        # Fetch and store gene trnascript data

        if not self.gene_id:
            print("No gene ID found, cannot fetch transcripts")
            return
        
        url = f"{self.BASE_URL}/overlap/id/{self.gene_id}/?feature=transcript;content-type=application/json"

        response = requests.get(url)
        
        temp_transcripts = []

        if response.status_code == 200:
            temp_transcripts = response.json()
        else:
            print(f"{response.status_code}")

        # Create lists of cannonical and non-cannonical transcripts, they are lsited from oldest to newest (newest is at index -1)
        # you will get more from this than a manual lookup, this function returns 58 total transcripts for human CFTR but only 38 shown on manual search 
        canonical_transcripts = []
        
        for transcript in temp_transcripts:
            if transcript.get('is_canonical') == 1:
                canonical_transcripts.append(transcript)


        # sets the primary transcript to the most recent canonical transcript (end of list)
        self.transcript = canonical_transcripts[-1]

    def get_dna_sequence(self):
        
        url = f"{self.BASE_URL}/sequence/id/{self.gene_id}?content-type=application/json"

        response = requests.get(url)

        if response.status_code == 200:
            sequence_info = response.json()
        else:
            print(f"{response.status_code}")

        if sequence_info.get('seq') is not None:
            self.dna_sequence = sequence_info.get('seq')
        else:
            print("No sequence data found")

    def get_cdna_sequence(self):
        
        url = f"{self.BASE_URL}/sequence/id/{self.transcript.get('id')}?content-type=application/json"

        response = requests.get(url)

        if response.status_code == 200:
            sequence_info = response.json()
        else:
            print(f"{response.status_code}")

        if sequence_info.get('seq') is not None:
            self.cdna_sequence = sequence_info.get('seq')
        else:
            print("No sequence data found")

    def get_exons(self):

        url = f"{self.BASE_URL}/overlap/id/{self.transcript.get('id')}/?feature=exon;content-type=application/json"
    
        response = requests.get(url)

        if response.status_code == 200:
            all_exons = response.json()

        for exon in all_exons:
            if exon.get('Parent') == self.transcript.get('id'):
                self.exons.append(exon)

        self.gene_start = self.exons[0].get('start')