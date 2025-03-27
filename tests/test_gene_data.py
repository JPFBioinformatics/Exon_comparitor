from gene_data import GeneData

human_gene = GeneData("CFTR", species = 'Human')
mouse_gene = GeneData("CFTR", species = "Mouse")


human_gene.get_transcripts()
mouse_gene.get_transcripts()

human_gene.get_cdna_sequence()
mouse_gene.get_cdna_sequence()

human_gene.get_exons()
mouse_gene.get_exons()

print(f"Human cdna sequence lenth: {len(human_gene.cdna_sequence)}")
print(f"Mouse cdna sequence length: {len(mouse_gene.cdna_sequence)}")