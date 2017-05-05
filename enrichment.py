import statsmodels.stats.multitest as smm
import numpy as np
import random
from scipy.stats import fisher_exact

import settings
import resources
import conversion

def _calculate_contingency_table(gene_list, gene_set, background_gene_list):
    table = np.zeros([2, 2])
    genes_in_set = np.array([gene in gene_set for gene in gene_list])
    bg_in_set = np.array([gene in gene_set for gene in background_gene_list])
    table[0, 0] = genes_in_set.sum()
    table[0, 1] = (1 - genes_in_set).sum()
    table[1, 0] = bg_in_set.sum()
    table[1, 1] = (1 - bg_in_set).sum()
    return table

def _convert_lists_to_entrez(gene_list, id_type):
    if id_type != 'entrez':
        return conversion.convert(gene_list, id_type, 'entrez')
    else:
        return gene_list

def correct_enrichment_pvalues(enrichments, method, sig_cutoff):
    corrected_enrichments = []
    for enrichment in enrichments:
        pvalues = enrichment.values()
        gene_set_names = enrichment.keys()
        if method == 'none' or method is None:
            corrected_pvalues = pvalues
            reject = pvalues > sig_cutoff
        else:
            reject, corrected_pvalues, _, _ = smm.multipletests(pvalues,
                                                        alpha=sig_cutoff,
                                                        method=method)
        accepted_indices = np.where(reject)[0]
        accepted_pvalues = dict([(gene_set_names[i], corrected_pvalues[i]) 
                                    for i in accepted_indices])
        corrected_enrichments.append(accepted_pvalues)
    return corrected_enrichments


def intramodular_connections(matrix, subset_indices):
    return matrix[subset_indices, :][:, subset_indices].sum()

#allowed id types are 'symbol', 'entrez', 'ensembl', 'ensembl_protein'
def network_enrich(adjacency_matrices, gene_list, id_type='entrez',
                   network_stat=intramodular_connections,
                   resource_list=settings.DEFAULT_GENE_SETS,
                   n_random_sets=1000,
                   sig_cutoff=0.05,
                   multitest_method='fdr_bh'):
    gene_list = np.array(gene_list)
    gene_sets_to_test = resources.read_resource_list(resource_list)
    processed_gene_list = _convert_lists_to_entrez(gene_list, id_type)
    enrichments = [{}]*len(adjacency_matrices)

    n_genes = len(processed_gene_list)
    all_indices = range(n_genes)
    for name, gene_set in gene_sets_to_test.iteritems():
        gene_set_indices = [i for i, g in enumerate(processed_gene_list) if g in gene_set]

        n_genes_in_set = len(gene_set_indices)
        if n_genes_in_set > n_genes:
            print 'Cannot test gene set %s, it is larger than gene list size!' % name
            continue
        if n_genes_in_set == 0:
            continue
        print 'Doing %s' % name
        random_sets = [random.sample(all_indices, n_genes_in_set) for _ in xrange(n_random_sets)]
        for matrix_idx, adj_matrix in enumerate(adjacency_matrices):
            pval = 0.
            gene_set_stat = network_stat(adj_matrix, gene_set_indices)
            for random_indices in random_sets:
                random_stat = network_stat(adj_matrix, random_indices)
                if gene_set_stat <= random_stat:
                    pval += 1
            pval /= n_random_sets
            enrichments[matrix_idx][name] = pval

    corrected_enrichments = correct_enrichment_pvalues(enrichments, 
                                                       multitest_method, 
                                                       sig_cutoff)
    return corrected_enrichments



def enrich(gene_lists, background_gene_list, id_type='entrez',
           resource_list=settings.DEFAULT_GENE_SETS, 
           sig_cutoff=0.05,
           test_fun=fisher_exact, multitest_method='fdr_bh'):
    gene_sets_to_test = resources.read_resource_list(resource_list)
    processed_gene_lists = _convert_lists_to_entrez(gene_lists, id_type)
    processed_bg_list = _convert_lists_to_entrez(background_gene_list, id_type)

    enrichments = []
    for gene_list in processed_gene_lists:
        contingency_tables = []
        for name, gene_set in gene_sets_to_test.iteritems(): 
            contingency_tables[name] = _calculate_contingency_table(gene_list, 
                                                                    gene_set, 
                                                                    processed_bg_list)
        pvalues = dict([(name, test_fun(ct)) for name, ct in contingency_tables.iteritems()])
        enrichments.append(pvalues)
    corrected_enrichments = correct_enrichment_pvalues(enrichments, 
                                                       multitest_method, 
                                                       sig_cutoff)
    return corrected_enrichments


