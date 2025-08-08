from Bio.Align import PairwiseAligner
aligner = PairwiseAligner()
print(hasattr(aligner, "maximum_number_of_alignments"))  # should be True
