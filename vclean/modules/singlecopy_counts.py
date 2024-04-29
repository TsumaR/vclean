import pandas as pd
import itertools
import subprocess as sp
import glob
import os

def run_interpro_pfam(protein_out, protein_fasta, basename, interproscan, cpu):
    asterix_removed_faa = protein_out + "/" + basename + "_remove_asterisx.faa"
    with open(asterix_removed_faa, 'w') as fp:
        rasterisx = sp.run(['sed', '-e', 's/*//g', protein_fasta], stdout=fp)

    interpro_out = protein_out + '/' + basename + '_interpro'
    sp.run([interproscan, '-f', 'tsv','-b', interpro_out, '--cpu', cpu, '-i', asterix_removed_faa, '-appl', 'Pfam'])
    interpro_res = pd.read_csv(interpro_out+".tsv", header=None, sep="\t", usecols=[0,1,2,3,4,8])

    # extract best hit
    interpro_res = interpro_res.reset_index()
    df_idx_max = interpro_res.loc[:,[0,8]].groupby(0).idxmax()
    df_idx_max.columns = ['max_idx']
    interpro_tophit = pd.merge(interpro_res, df_idx_max, left_on='index', right_on='max_idx')
    interpro_tophit = interpro_tophit.loc[:, [0,4]].rename(columns={0: "tmpname", 4: "pfam_id"})

    return interpro_tophit

def run_hmm_pfam(protein_out, protein_fasta, basename, pfamdb, cpu):
    pfam_out = protein_out + '/' + basename + '_pfam.txt'
    log_out = protein_out + '/' + 'log.txt'
    with open(log_out, "w") as f:
        sp.run(['hmmsearch', '--tblout', pfam_out, '--cpu', cpu, '-E', '1e-3', '--noali', pfamdb, protein_fasta], stdout=f)
    pfam_raw_result = pd.read_csv(pfam_out, comment="#", header=None,delim_whitespace=True).reset_index()
    pfam_tophit = pfam_raw_result.loc[pfam_raw_result.groupby(0)[4].idxmin()].loc[:,[0,3]].rename(columns={0:"contig_id", 3:"pfam_id"})
    print(pfam_tophit)
    print(pfam_tophit[pfam_tophit['pfam_id']=='PF03796.18'])
    return pfam_tophit


def simple_single_copy(pfam_tophit, single_like_pfam):
    single_pfam_list = list(pd.read_csv(single_like_pfam, header=None)[0])

    interpro_res2 = pd.concat([pfam_tophit, pfam_tophit["tmpname"].str.rsplit('_', 1, expand=True)], axis=1).rename(columns={0: "contig_id", 1: "gene_id"}).drop("tmpname", axis=1)
    interpro_res3 = interpro_res2[interpro_res2["pfam_id"].isin(single_pfam_list)]
    print("======== Detected Single Copy like genes ========")
    print(interpro_res3)
    if len(interpro_res3) == 0:
        num_single_like_overlaps = 0
    else:
        remove_dupli_in_same_contig = interpro_res3[~interpro_res3.duplicated(subset=['contig_id', 'pfam_id'])]
        df_dupli_pfam = remove_dupli_in_same_contig[remove_dupli_in_same_contig.duplicated(subset='pfam_id')]
        num_single_like_overlaps = len(df_dupli_pfam)

    return num_single_like_overlaps

def contraposition_single_copy(pfam_tophit, uniq_single_like_pfam):
    Terminase_N = ['PF03237.18', 'PF03354.18', 'PF05876.15']
    Terminase_C = ['PF17289.5', 'PF17288.5', 'PF20441.1', 'PF20454.1']
    minor_capsid = ['PF04233.17', 'PF11114.11']
    peptidase_S = ['PF14550.9', 'PF03420.16']
    holin = ['PF16945.8', 'PF04531.16']
    cont_single_list = [Terminase_N, Terminase_C, minor_capsid, peptidase_S, holin]
    num_count_single_like_overlaps = 0
    for t_pfams in cont_single_list:
        if len(pfam_tophit[pfam_tophit["pfam_id"].isin(t_pfams)]) >= 2:
            print(pfam_tophit[pfam_tophit["pfam_id"].isin(t_pfams)])
            num_count_single_like_overlaps += 1

    other_pfams = list(pd.read_csv(uniq_single_like_pfam, header=None)[0])
    for t_pfam in other_pfams:
        if len(pfam_tophit[pfam_tophit["pfam_id"] == t_pfam]) >= 2:
            print(pfam_tophit[pfam_tophit["pfam_id"] == t_pfam])
            num_count_single_like_overlaps += 1

    return num_count_single_like_overlaps
        
