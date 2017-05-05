from genoracle import enrichment
import numpy as np

gene_list = open('phosphosite_plus_kinome_homo_sapiens.tsv').readlines()[1].strip().split('\t')
n_genes = len(gene_list)
matrices = [np.random.rand(n_genes, n_genes)]
enrichments = enrichment.network_enrich(matrices, gene_list, id_type='symbol',
                                        resource_list=['gene_sets/phosphosite_plus_kinome_homo_sapiens.gmt'])
print enrichments
