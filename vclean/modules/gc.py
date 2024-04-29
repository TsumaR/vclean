from Bio import SeqIO, SeqUtils
import pandas as pd
import re
import numpy as np

def cal_gc_values(fasta):
    id_lists = []
    gc_lists = []
    cpg_lists = []
    gc_skew_lists = []
    for record in SeqIO.parse(fasta, 'fasta'):
        id_part = record.id
        desc_part = record.description
        seq = record.seq
        seq_len = len(seq)

        id_lists.append(id_part)
        # GC
        gc_contents = 100 * SeqUtils.gc_fraction(seq, ambiguous="ignore")
        gc_lists.append(gc_contents)
        # CpG
        cpg_counts = 0
        for ss in range(seq_len):
            if seq[ss:ss+2]=="CG":
                cpg_counts += 1
        cpg_lists.append(cpg_counts/seq_len)
        # GC-skew
        gc_skew = SeqUtils.GC_skew(seq, seq_len)
        gc_skew_lists.append(gc_skew[0])

    gc_table = pd.DataFrame(data={"Sequence Id":id_lists, "GC_contents":gc_lists, "CPG_contents":cpg_lists, "GC_skew":gc_skew_lists},
                            columns = ["Sequence Id", "GC_contents", "CPG_contents", "GC_skew"])
    gc_table = gc_table.set_index("Sequence Id")

    return gc_table

def gc_to_features(gc_table):
    print(len(gc_table))
    if len(gc_table) >= 2:
        gc_var = np.var(gc_table["GC_contents"], ddof=1)
        gc_max_min = gc_table["GC_contents"].max() - gc_table["GC_contents"].min()
        cpg_var = np.var(gc_table["CPG_contents"], ddof=1)
        cpg_max_min = gc_table["CPG_contents"].max() - gc_table["CPG_contents"].min()
        skew_var = np.var(gc_table["GC_skew"], ddof=1)
        skew_max_min = gc_table["GC_skew"].max() - gc_table["GC_skew"].min()
        results = [gc_var, gc_max_min, cpg_var, cpg_max_min, skew_var, skew_max_min]
    else:
        results = [0,0,0,0,0,0]

    return results
