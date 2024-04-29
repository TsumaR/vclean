import pandas as pd
import subprocess as sp

def count_redundant_protein(protein_fasta, tmp_dir):
    # To calculate protein redundancy, all proteins within a bin are clustered using Mmseqs2
    # Any proteins clustered within a sample, excluding those along the same scaffold, are considered redundant.
    mmseqs_out_prefix = tmp_dir + '/mmseqs'
    tmp_log = tmp_dir + '/mmseqs_log.txt'
    with open(tmp_log, "w") as f:
        sp.run(['mmseqs', 'easy-linclust', protein_fasta, mmseqs_out_prefix, tmp_dir, '--min-seq-id', '0.5', '-c', '0.8', '-e', '0.01', '--min-aln-len', '50', '--cluster-mode', '0', '--seq-id-mode', '0', '--alignment-mode', '3', '--cov-mode', '5', '--kmer-per-seq', '75'], stdout=f)

    mmseq_out = mmseqs_out_prefix + '_cluster.tsv'

    num_redundant_protein = 0
    with open(mmseq_out) as f:
        for line in f:
            each_cluster = line.split()
            cluster_member_contigs = set([x.rsplit('_', 1)[0] for x in each_cluster])
            num_redundant_protein += (len(cluster_member_contigs) - 1)

    print('++++++++++++++++++++++++ print num_redundant_protein ++++++++++++++++++++++++++++++')
    print(num_redundant_protein)
    print('++++++++++++++++++++++++ printed num_redundant_protein ++++++++++++++++++++++++++++++')
    return num_redundant_protein
