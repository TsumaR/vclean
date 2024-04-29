import subprocess as sp
import pandas as pd

def run_checkv(fasta, tmp_dir, ref_dir, cpu):
    checkv_command = sp.run(["checkv", "completeness", fasta, tmp_dir, "-d", ref_dir, "-t", cpu, "--quiet"])
    res = pd.read_table(tmp_dir + "/completeness.tsv")
    sum_comp = res['aai_completeness'].sum()
    return [int(sum_comp)]
    # print(res)
