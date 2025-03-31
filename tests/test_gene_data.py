from gene_data import GeneData

human_gene = GeneData("CFTR", species = 'Human')
mouse_gene = GeneData("CFTR", species = "Mouse")

human_gene.get_cdna_sequence()
human_gene.get_aa_sequence()
human_gene.get_exons()