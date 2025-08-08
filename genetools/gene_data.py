import requests
from genetools.codon_table import CodonTable
from Bio.Align import PairwiseAligner, substitution_matrices
from maxentpy.maxent import score5, score3
import plotly.graph_objs as go


class GeneData:

    def __init__(self, gene_name = None, species = 'Human'):
        self.BASE_URL = "https://rest.ensembl.org"
        self.gene_name = gene_name
        self.species = species
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

            # store canonical transacript ID seperate for easy use
            primary_transcript = g_data.get('canonical_transcript')
            self.transcript_id = primary_transcript.split(".")[0]
            
            # starting location of gene, not the same as start of fist exon
            self.gene_start = g_data['start']

        # get dna sequence and exons
        self.get_dna_sequence()
        self.get_exons()
        self.aa_map_exons()

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

    # gets the transcript for each exon for you gene
    def get_aa_sequence(self):
           
        url = f"{self.BASE_URL}/sequence/id/{self.transcript_id}?content-type=application/json;type=protein"
        response = requests.get(url)

        if response.status_code == 200:

            data = response.json()
            aa_seq = data['seq']

        else:
            print(f"No aa sequence found for transcript {self.transcript_id}")

        return aa_seq

    # takes aa_sequence and the exon list and finds what group of aa a given exon codes
    def aa_map_exons(self):
        
        # list to hold all aa sequences for confirmation
        list =[]

        # grab AA sequence for entire transcript to map
        seq = self.get_aa_sequence()

        # create a codon table for comparison
        codon_table = CodonTable()

        # get list of exons for mapping
        exons = self.exons

        # start a counter for where in the sequence we are
        seq_idx = 0

        for idx,exon in enumerate(exons):

            exon_aas = []
            exon_seq = self.get_exon_cdna(exons[idx])
            if idx+1 < len(exons):
                next_exon = self.get_exon_cdna(exons[idx+1])
            start_idx = 0
            

            # get start and end phases of exon
            start_phase = exon['ensembl_phase']
            end_phase = exon['ensembl_end_phase']

            # adjust sequence to start phase
            if start_phase == 1:
                exon_seq = exon_seq[2:]
            if start_phase == 2:
                exon_seq = exon_seq[1:]
            
            # adjust section to end phase
            if end_phase == 1:
                chunks = []
                chunks.append(exon_seq)
                chunks.append(next_exon[:2])
                exon_seq = ''.join(chunks)
            if end_phase == 2:
                chunks = []
                chunks.append(exon_seq)
                chunks.append(next_exon[:1])
                exon_seq = ''.join(chunks)
            
            section = exon_seq

            # loop over every nucleotide in sequence and map amino acids
            for _ in exon_seq:
                
                # take the section of the exon cdna for query
                section = section[start_idx:]                 

                # locate the rseidue we are searching for
                if seq_idx >= len(seq):
                    break
                residue = seq[seq_idx]

                # find the first instance of a codon for that residue in the query section
                next_idx = codon_table.find_codon(residue,section)

                # if no matching codon found in exon print error message
                if next_idx == None:
                    #print(f"Could not find codon for amino acid {residue} in exon {exon['rank']}")
                    break

                else:
                    # increment sequence index and start index
                    seq_idx += 1
                    start_idx = next_idx
                    # add residue to list of resiudes for this amnino acid
                    exon_aas.append(residue)

            # save the exon amino acid sequence to exon dictionary
            aas = ''.join(exon_aas)
            exon['aa_seq'] = aas
            list.append(aas)

        # check to make sure the mapped exons sum to be the same as the origonal sequence
        generated_seq = ''.join(list)
        if generated_seq == seq:
            print(f"Amino Acids successfully mapped to exons")

    # preforms a Smith_Waterman local alignment of two given exon amino acid sequences
    @staticmethod
    def align_aa(seq1, seq2, sub_matrix = "BLOSUM62", score_mode = "blastp", align_mode = "global"):
        
        # load the BLOSUM62 substituion matrix, other options are 45,50,80,90 corrosponding to expected % similarity in sequence
        blosum = substitution_matrices.load(sub_matrix)

        # create new PairwiseAligner object with default parameters and score_mode scoring
        aligner = PairwiseAligner(scoring=score_mode)

        # assign substituion matrix to aligner
        aligner.substitution_matrix = blosum

        # set mode to global, will align the entire sequence as best it can, swich it to "local" to align based on most similar "chunk" of exons
        aligner.mode = align_mode

        # align input exons
        alignments = aligner.align(seq1,seq2)
        
        # return best alignment
        return alignments[0]

    # preforms a Smith_Waterman local alignment of two given exon amino acid sequences
    @staticmethod
    def align_cdna(seq1, seq2, score_mode = "blastn", align_mode = "global"):

        # create new PairwiseAligner object with default parameters and score_mode scoring
        aligner = PairwiseAligner(scoring=score_mode)

        # set mode to global, will align the entire sequence as best it can, swich it to local to align based on most similar "chunk" of exons
        aligner.mode = align_mode

        # align input exons
        alignments = aligner.align(seq1,seq2)
        
        # return best alignment
        return alignments[0]
        
    # gets the cdna sequence of a given exon
    def get_exon_cdna(self,exon):

        # normalized start and end indices
        start = exon['start'] - self.gene['start']
        end = exon['end'] - self.gene['start'] + 1

        # grab exon's nuc sequence from dna
        cdna = self.dna_sequence[start:end]

        return cdna
    
    # makes sure all values in cdna are valid for splice site analysis
    @staticmethod
    def is_valid(seq):
        return all(base in "ACGT" for base in seq)

    # preforms a splice site analysis using maxentpy libray (from MaxEntScan)
    # <0 means non-splice site, 0-5 is ok splice site, 6-12 is good splice site (higher number better)
    def splice_site_analysis(self,idx):

        # location of first and last nuclotides in exon on the dna sequence
        normal_start = int(self.exons[idx]['start']) - int(self.gene['start'])
        # +1 to make sure 3' is in phase (GT consensus at +1,+2 base 1 gene reading)
        normal_end = int(self.exons[idx]['end']) - int(self.gene['start']) + 1

        # find normalized location of bounds for 3'ss analysis on full dna sequence
        splice3_start = normal_start - 20
        splice3_end = normal_start + 3
        
        # generate nucleotide sequence for 3'ss analysis
        try:
            splice3_seq = self.dna_sequence[splice3_start:splice3_end]
        except IndexError:
            print(f"3' splice site out of index")

        # find normalized location of bounds for 5'ss analysis on full dna sequence
        splice5_start = normal_end - 3
        splice5_end = normal_end + 6

        # generate nucleotide sequence for 5'ss analysis
        try:
            splice5_seq = self.dna_sequence[splice5_start:splice5_end]
        except IndexError:
            print(f"5' splice site out of index")

        # calculate scores if lengths are correct
        if len(splice3_seq) == 23:
            score_3ss = score3(splice3_seq)
        else:
            print(f"3'ss sequence wrong length: {len(splice3_seq)}")
            score_3ss = 'Error'

        if len(splice5_seq) == 9:
            score_5ss = score5(splice5_seq)
        else:
            print(f"5'ss sequence wrong length: {len(splice5_seq)}")
            score_5ss = 'Error'


        # replace 3' of exon 0 and 5' of last exon with 'N/A'
        if idx == 0:
            score_3ss = 'N/A'
        if idx == len(self.exons)-1:
            score_5ss = 'N/A'

        # return 3' and 5' splice site scores
        return score_3ss,score_5ss

    # swaps an exon (rank = num) out of sequence and generates new, short chunk for SSA
    def swap_exons(self,num,new_seq):

        # create list to hold the "chunks" of nucleotides to build the new sequence
        nuc_chunks = []

        # grab the exon to replace
        exon = self.exons[num]

        # idx of 1st nucleotide in exon
        og_start = int(exon['start']) - int(self.gene['start'])
        # idx of last nucleotide in exon
        og_end = int(exon['end']) - int(self.gene['start']) + 1

        # get the header (20 nucleotides before spice site from OG dna seq)
        header = self.dna_sequence[og_start-20:og_start]
        # get the footer (6 nucleotides after splice site from OG dna seq)
        footer = self.dna_sequence[og_end:og_end+6]

        # build new string
        nuc_chunks.append(header)
        nuc_chunks.append(new_seq)
        nuc_chunks.append(footer)

        # build string
        result = ''.join(nuc_chunks)

        # change sequence to 'N/A' if it is 3'ssa for exon 0 or 5'ssa for last exon
        if num == 0:
            return 1, result
        if num == len(self.exons)-1:
            return -1, result

        # idx = 20 is first nuc of exon, idx = -7 is last nuc of exon
        # result is a new sequence including 20 nucleotide header and 6 nucleotide footer for ssa
        return 0, result
    
    # runs SSA on the swapped chunk of sequence
    @staticmethod
    def swapped_ssa(check, sequence):
        
        # start/end indices are always the same on swapped sequence
        start_idx = 20
        end_idx = -7

        # get sequencse for SSAs
        try:
            splice3_seq = sequence[:start_idx+3]    # add 3 so we get the -20 to +3 for ssa
        except IndexError:
            print("3' ss out of index")
        try:
            splice5_seq = sequence[end_idx-2:]      # subtract 2 so we get the -3 to +6 for ssa
        except IndexError:
            print("5' ss out of idex")
    
        # check to make sure sequences are correct length
        if len(splice3_seq) != 23:
            print(f"3'ss sequence wrong length: {len(splice3_seq)}")
            return None, None
        if len(splice5_seq) != 9:
            print(f"5'ss sequence wrong length: {len(splice5_seq)}")
            return None, None
        
        # calculate scores
        score_3ss = score3(splice3_seq)
        score_5ss = score5(splice5_seq)

        # check to make sure fist exon 3'ssa and last exon 5'ssa are both returned as N/A
        if check == 1:
            score_3ss = 'N/A'
        if check == -1:
            score_5ss = 'N/A'

        # return ssa for hybrid gene
        return score_3ss,score_5ss
