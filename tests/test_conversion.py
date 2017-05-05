from genoracle import conversion
genes = [line.strip().split('\t') for line in open('phosphosite_plus_kinome_homo_sapiens.tsv')]
converted_genes = conversion.convert(genes, 'symbol', 'entrez', ignore_null=True)
for i, gene_list in enumerate(converted_genes):
    print '%s\thttp://www.phosphosite.org\t%s' % (genes[i][0], '\t'.join(gene_list[1:]))
