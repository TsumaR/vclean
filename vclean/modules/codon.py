from Bio import SeqIO, SeqUtils
import pandas as pd
import argparse
import itertools
from sklearn.decomposition import PCA
import numpy as np

def codon_calculate(gene_file):
    DNA = ['A', 'C', 'G', 'T']
    codon_list = []
    codon_tuple_list = list(itertools.product(DNA,DNA,DNA))
    for codon in codon_tuple_list:
        codon_list.append("".join(codon))

    codon_table = pd.DataFrame(index=[],columns=[codon_list])
    for cds in SeqIO.parse(gene_file, 'fasta'):
        id_part = cds.id
        desc_part = cds.description
        seq = cds.seq
        contig_id = id_part.rsplit("_", 1)[0]

        if contig_id not in codon_table.index:
            codon_table.loc[contig_id] = 0
        for c_number in range(0,len(seq)//3,3):
            codon = str(seq[c_number:c_number+3])
            if "N" not in codon:
                codon_table.loc[contig_id, codon] = int(codon_table.loc[contig_id, codon]) + 1
    codon_table = codon_table.apply(lambda x:x/sum(x),axis=1)

    return codon_table

def calculate_pca_variance(target_table):
    if len(target_table) >= 2:
        pca = PCA()
        pca.fit(target_table)
        feature = pca.transform(target_table)
        pca_variance = np.var(feature[:,0], ddof=1)
        max_min_pca = max(feature[:,0]) - min(feature[:,0])
    else:
        pca_variance = 0
        max_min_pca = 0

    return pca_variance, max_min_pca


def cos_sim_matrix(target_table):
    matrix = np.array(target_table)
    d = matrix @ matrix.T
    norm = (matrix * matrix).sum(axis=1, keepdims=True) ** .5
    res_cos_matrix = d / norm / norm.T
    return res_cos_matrix.min()
