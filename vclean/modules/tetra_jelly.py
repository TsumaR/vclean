from Bio import SeqIO, SeqUtils
import pandas as pd
import argparse
import itertools
import subprocess
import glob
import os
from sklearn.decomposition import PCA
import numpy as np

def calculate_pca_variance(input_tsv):
    target_table = input_tsv.T
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

def cos_sim_matrix(input_tsv):
    matrix = np.array(input_tsv).T
    d = matrix @ matrix.T
    norm = (matrix * matrix).sum(axis=1, keepdims=True) ** .5
    res_cos_matrix = d / norm / norm.T
    return res_cos_matrix.min()


def cal_4mer_freq(fasta_path, tmpdir):
    DNA = ['A', 'C', 'G', 'T']
    kmer_tuple_list = list(itertools.product(DNA,DNA,DNA,DNA))
    kmer_list=[]
    for codon in kmer_tuple_list:
        kmer_list.append("".join(codon))
    kmer_table = pd.DataFrame(index=kmer_list,columns=[])
    
    for contig in SeqIO.parse(fasta_path, 'fasta'):
        tmp_out = tmpdir + "/tmp.js"
        tmp_count = tmpdir + "/tmp_count.txt"
        tmp_fasta = tmpdir + "/tmp_contig.fasta"
        SeqIO.write(contig, tmp_fasta, "fasta")
        jelly_count_cmd = ["jellyfish", "count", "-m 4", "-s 10000000", "-o", tmp_out, tmp_fasta]
        jelly_dump_cmd = ["jellyfish", "dump", "-c", tmp_out, "-o", tmp_count]
        subprocess.run(jelly_count_cmd)
        subprocess.run(jelly_dump_cmd)
        tmp_res = pd.read_table(tmp_count, sep=" ", header=None, index_col=0).rename(columns={1: contig.description })
        kmer_table = pd.merge(kmer_table, tmp_res, left_index=True, right_index=True, how="outer").fillna(0)
        kmer_table = kmer_table/kmer_table.sum()

    return kmer_table
