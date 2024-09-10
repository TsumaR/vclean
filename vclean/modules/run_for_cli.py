#!/usr/bin/env python

import argparse
import subprocess as sp
import sys
import random, string
import os
import time
import pandas as pd
import numpy as np
import pickle
import glob
import concurrent.futures
import lightgbm as lgb
from Bio import SeqIO
import shutil

# original scripts
# import seqkit as seqkit
from vclean.modules import seqkit as seqkit
from vclean.modules import checkv as checkv
from vclean.modules import gc as gc
from vclean.modules import tetra_jelly as tetra
from vclean.modules import penta_jelly as penta
from vclean.modules import codon as codon
from vclean.modules import estimate_simu_contami as contami
from vclean.modules import singlecopy_counts as scopy
from vclean.modules import count_redundant_protein as red_pr
from vclean.modules import lgb_model as lgb_model

#################################################################################
############### Process for each file ###########################################
#################################################################################
def process_file(input_fasta, outdir, run_mode, tmpout, cpu, checkv_db_dir, pfam_db,
                input_protein, input_gene, trans_table, single_like_pfam):
    gene_set_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    single_like_pfam = os.path.join(gene_set_dir, 'db', 'single_pfam_099.txt')

    fasta_file_name = os.path.basename(input_fasta)
    basename = os.path.splitext(fasta_file_name)[0]

    sp.run(['mkdir', '-p', tmpout])
    random_val = ''.join(random.choices(string.ascii_letters + string.digits, k=15))
    tmp_dir = tmpout + "/" + random_val
    sp.run(['mkdir', '-p', tmp_dir])

    # --------------------------------------------------------------------------
    #                              Contami calculation
    # --------------------------------------------------------------------------
    contami_res = outdir + '/' + "Contami_value.txt"
    if run_mode == 'T':
        contami_val = contami.cal_contamination(input_fasta)
        print(contami_val)
        print(str(contami_val))
        print(basename + ',' + str(contami_val))
        with open(contami_res, 'a') as c_txt:
            print(basename + ',' + str(contami_val), file=c_txt)
        sys.exit()

    # --------------------------------------------------------------------------
    #                                run seqkit
    # --------------------------------------------------------------------------
    seqkit_list = seqkit.run_seqkit(input_fasta)

    # --------------------------------------------------------------------------
    #                                run checkV
    # --------------------------------------------------------------------------
    sum_comp = checkv.run_checkv(input_fasta, tmp_dir, checkv_db_dir, cpu)
    print("finish checkv")

    # --------------------------------------------------------------------------
    #                    gc-contents, gc-skew, cpg-contents 
    # --------------------------------------------------------------------------
    gc_table = gc.cal_gc_values(input_fasta)
    gc_res_6items = gc.gc_to_features(gc_table)

    # --------------------------------------------------------------------------
    #                            Tetra, codon, penta
    # --------------------------------------------------------------------------
    ## tetra
    tetra_table = tetra.cal_4mer_freq(input_fasta, tmp_dir)
    tetra_pca, mm_tetra_pca = tetra.calculate_pca_variance(tetra_table)
    tetra_cos_mm = tetra.cos_sim_matrix(tetra_table)

    tetra_res = [tetra_pca, mm_tetra_pca, tetra_cos_mm]

    ## penta
    penta_table = penta.cal_5mer_freq(input_fasta, tmp_dir)
    penta_pca, mm_penta_pca = penta.calculate_pca_variance(penta_table)
    penta_cos_mm = penta.cos_sim_matrix(penta_table)
    penta_res = [penta_pca, mm_penta_pca, penta_cos_mm]

    print("finish tetra, codon, penta")
    # --------------------------------------------------------------------------
    #                            Protein prediction
    # --------------------------------------------------------------------------
    if input_protein is None:
        protein_out = outdir + '/' + basename
        sp.run(['mkdir', '-p', protein_out])
        protein_fasta = protein_out + "/" +  basename + ".faa"
        gene_fasta = protein_out + "/" + basename + ".fna"
        gff_out = protein_out + "/" + basename + ".gff"

        sp.run(['prodigal', '-i', input_fasta, '-p', 'meta', '-a', protein_fasta,  '-d', gene_fasta, '-o', gff_out, '-f', 'gff', '-g', trans_table, '-q'])
    else:
        protein_fasta = input_protein
        gene_fasta = input_gene

    print("finish protein prediction")
    # --------------------------------------------------------------------------
    #                          Single-copy-like genes
    # --------------------------------------------------------------------------
    hmm_tophit = scopy.run_hmm_pfam(protein_out, protein_fasta, basename, pfam_db, cpu)
    # num_single_like_overlaps = scopy.contraposition_single_copy(hmm_tophit, uniq_single_like_pfam)
    num_single_like_overlaps = scopy.contraposition_single_copy(hmm_tophit, single_like_pfam)

    # ----------------------protein redundancy -------------------------
    redundant_protein = red_pr.count_redundant_protein(protein_fasta, tmp_dir)

    # ------------------------ codon usage -----------------------------
    codon_table = codon.codon_calculate(gene_fasta)
    codon_pca, mm_codon_pca = codon.calculate_pca_variance(codon_table)
    codon_cos_mm = codon.cos_sim_matrix(codon_table)

    codon_res = [codon_pca, mm_codon_pca, codon_cos_mm]

    # --------------------------------------------------------------------------
    #                            Summary
    # --------------------------------------------------------------------------
    result_list = seqkit_list + sum_comp + gc_res_6items + tetra_res + penta_res + [num_single_like_overlaps] + codon_res + [redundant_protein]
    result_columns = ['file', 'num_seqs', 'sum_len', 'min_len', 'avg_len', 'sum_completeness', 'gc_var', 'gc_max_min', 'cpg_var', 'cpg_max_min', 'skew_var', 'skew_max_min', 'tetra_pca', 'mm_tetra_pca', 'tetra_cos_mm', 'penta_pca', 'mm_penta_pca', 'penta_cos_mm', 'overlap_singl_like_gene', 'codon_pca', 'mm_codon_pca', 'codon_cos_mm', 'redundant_protein']
    result_df = pd.DataFrame(index=[], columns=result_columns)
    result_list = [basename] + result_list
    add_series = pd.Series(result_list, index=result_df.columns, name=basename)
    # result_df = result_df.append(add_series)
    add_series = add_series.to_frame().T
    result_df = pd.concat([result_df, add_series], ignore_index=True)
    result_path = outdir + '/' + basename + '_result.tsv'
    result_df.to_csv(result_path, sep='\t', index=False)
    # --------------------------------------------------------------------------
    #                           Cleaning
    # --------------------------------------------------------------------------
    sp.run(['rm', '-r', tmp_dir])

#################################################################################
############### parallel works  #################################################
#################################################################################
def process_files_in_folder(input_dir, outdir, run_mode, tmpout, cpu, checkv_db_dir, pfam_db, input_protein, input_gene, trans_table, single_like_pfam):
    # 出力ディレクトリを作成します
    os.makedirs(outdir, exist_ok=True)

    with concurrent.futures.ThreadPoolExecutor() as executor:
        # 各ファイルについて並列に処理を行います
        futures = [executor.submit(process_file, entry.path, outdir, run_mode, tmpout, cpu, checkv_db_dir, pfam_db, input_protein, input_gene, trans_table, single_like_pfam) for entry in os.scandir(input_dir) if entry.is_file()]
        for future in concurrent.futures.as_completed(futures):
            try:
                future.result()
            except Exception as e:
                print(f"Exception occurred while processing file: {e}")

#################################################################################
############### Main  ###########################################################
#################################################################################
def fetch_arguments(parser):
    parser.set_defaults(func=main)
    parser.set_defaults(program="run")
    
    parser.add_argument('input', type=str, help='Put the input fasta directory')
    parser.add_argument('output', type=str, help='Put the output directory')
    parser.add_argument('-d', '--db', action='store', dest='db', help='Set the database directory path. By default the VCLEANDB environment vairable is used')
    parser.add_argument('-tmp', action='store', dest='tmp', default='./tmp', help='Set the path of temporary file directory')
    parser.add_argument('-t', '--threads', action='store', dest='threads', default='1', help='Put the number of CPU to use. default=1')
    parser.add_argument('-p', '--protein', action='store', dest='protein', default=None,
                    help='Not nessesary. you can input protein fasta file if you have. In default, vDeteCon predict CDS using prodigal from nucleotide fasta file.')
    parser.add_argument('-n', '--nucleotide', action='store', dest='gene', default=None)
    parser.add_argument('--translate_table', action='store', dest='t_table', default=11, help='put the translate table, default=11')
    parser.add_argument('-m', '--mode', action='store', dest='mode', default='F', help='True or False. If you want to calculate contamination value of simulation data, set this value True')
    parser.add_argument('--skip_feature_table', action='store', dest='skip_feature_table', default=False, help='If you set True, skip features prediction step.')
    parser.add_argument('--skip_lgb_step', action='store', dest='skip_lgb_step', default=False, help='If you set True, skip contamination prediction step.')
    parser.add_argument('--f_table', action='store', dest='f_table', help="You want to only run lgb model, you have to input the features table path.")
    parser.add_argument('-pr', '--threshold', action='store', dest='threshold', default=0.6, help='Put the threshold for the Contamination probability rate value. default=0.6. if the contamination probability value is over the set score, the input fasta are assigned as CONTAMINATION.')

    # --------------------------------------------------------------------------
    # interproscan = "/home/mako43/packages/interproscan/interproscan-5.47-82.0/interproscan.sh"
    # single_like_pfam = "/home/wagatsuma/analysis/virus/Tool/tmp/vDeteCon/db/single_pfam.txt"
    # uniq_single_like_pfam = "/home/wagatsuma/analysis/virus/Tool/tmp/vDeteCon/db/unique_single_pfam.txt"
    # --------------------------------------------------------------------------

def check_db(dbdir):
    "check existence of database"
    if dbdir is None:
        if "VCLEANDB" not in os.environ:
            msg="Error: database dir not specified\nUse -d or set VCLEANDB environmental variable"
            sys.exit(msg)
        else:
            print("use VCLEANDB environmental variable")
            dbdir = os.environ["VCLEANDB"]
    dbdir = os.path.abspath(dbdir)
    if not os.path.exists(dbdir):
        msg = f"Error: database dir not found '{dbdir}'"
        sys.exit(msg)

    return dbdir

def path_to_db(dbdir):
    all_dbdir = os.listdir(dbdir)
    # checkvdb
    checkv_db_dir = None
    for d in all_dbdir:
        if d.startswith('checkv-db') and os.path.isdir(os.path.join(dbdir, d)):
            checkv_db_dir = os.path.join(dbdir, d)
            break
    if checkv_db_dir:
        print("Found checkv-db directory:", checkv_db_dir)
    else:
        msg = f"Error: No checkv-db directory found '{dbdir}'"
        sys.exit(msg)
    # pfamdb 
    pfam_db = None
    for d in all_dbdir:
        if d.startswith('Pfam') and os.path.isdir(os.path.join(dbdir, d)):
            pfam_db = os.path.join(dbdir, d, "Pfam-A.hmm")
            break
    if pfam_db:
        print("Found pfam-db:", pfam_db)
    else:
        msg = f"Error: No pfam-db '{dbdir}'"
        sys.exit(msg)

    return checkv_db_dir, pfam_db

def main(args):
    program_start = time.time()
    input_fasta_dir = args["input"]
    tmpout = args["tmp"]
    cpu = args["threads"]
    input_protein = args["protein"]
    input_gene = args["gene"]
    
    final_out_dir = args["output"]
    sp.run(['mkdir', '-p', final_out_dir])
    outdir = final_out_dir + "/vClean_tmp"
    out_fasta_dir = final_out_dir + "/purified_fasta"
    os.makedirs(outdir, exist_ok=True)
    os.makedirs(out_fasta_dir, exist_ok=True)
    
    trans_table = str(args["t_table"])
    program_mode = args["mode"]
    contami_threshold = args["threshold"]
    f_table = args["f_table"]
    db_dir = check_db(args["db"])
    checkv_db_dir, pfam_db = path_to_db(db_dir)

    gene_set_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    single_like_pfam = os.path.join(gene_set_dir, 'db', 'single_pfam_099.txt')

    # --------------------------------------------------------------------------
    # sp.run(['mkdir', '-p', tmpout])
    # random_val = ''.join(random.choices(string.ascii_letters + string.digits, k=15))
    # tmp_dir = tmpout + "/" + random_val
    # sp.run(['mkdir', '-p', tmp_dir])

    # ------------------------- produce feature tables -------------------------------------
    
    if args["skip_feature_table"] == False:
        process_files_in_folder(input_fasta_dir, outdir, program_mode, tmpout, cpu, checkv_db_dir, pfam_db, input_protein, input_gene, trans_table, single_like_pfam)
        tsv_files = glob.glob(os.path.join(outdir, '*.tsv'))
        result_dataframes = [pd.read_csv(file, sep='\t') for file in tsv_files]
        features_table = pd.concat(result_dataframes, ignore_index=True)
        features_table_out_name = final_out_dir + "/features_result.tsv"
        features_table.to_csv(features_table_out_name, sep="\t", index=False)
    else:
        features_table_out_name = final_out_dir + "/features_result.tsv"
        feature_table = pd.read_csv(features_table_out_name, sep="\t")
    
    # --------------------------------------------------------------------------
    #                          Run lightgbm
    # --------------------------------------------------------------------------
    
    if args["skip_lgb_step"] == False:
        result_for_input = lgb_model.table_reformat(features_table)
    
        res_drop_columns =  ["penta_pca", "mm_codon_pca", "gc_max_min", "mm_penta_pca", "mm_tetra_pca", "tetra_cos_mm", "penta_cos_mm", "codon_cos_mm"]
        result_for_input = result_for_input.drop(res_drop_columns, axis=1)

        script_dir = os.path.dirname(os.path.abspath(__file__))
        model_dir = os.path.dirname(os.path.dirname(script_dir))
        model_path = os.path.join(model_dir, 'models', 'best_lgb_model_v1.3.2.pickle')

        y_pred_proba_list = []
        for i in range(5):
            model_name = f'best_lgb_model_fold{i}.txt'
            model_path = os.path.join(model_dir, 'models', model_name)
            loaded_model = lgb.Booster(model_file=model_path)
            y_pred_proba = loaded_model.predict(result_for_input)
            y_pred_proba_list.append(y_pred_proba)
    
        y_pred_proba_avg = np.array(y_pred_proba_list).mean(axis=0)
        y_pred = np.where(y_pred_proba_avg >= float(contami_threshold), 1, 0)
    
        contamination_probability_result = pd.DataFrame({'file' : features_table['file'],
                                'y_pred1':y_pred_proba_list[0],
                                'y_pred2':y_pred_proba_list[1],
                                'y_pred3':y_pred_proba_list[2],
                                'y_pred4':y_pred_proba_list[3],
                                'y_pred5':y_pred_proba_list[4],
                                'result_prediction': y_pred,
                                'y_ave':y_pred_proba_avg})
        # contamination_probability_result = pd.merge(contamination_probability_result_over2contigs, features_table['file'], on='file', how='right')
        # column_to_fill = ['y_pred1', 'y_pred2', 'y_pred3', 'y_pred4', 'y_pred5', 'result_prediction', 'y_ave']
        # contamination_probability_result[columns_to_fill] = contamination_probability_result[columns_to_fill].fillna(0)
        print(contamination_probability_result)
        # contamination_probability_result = contamination_probability_result.assign(file=contamination_probability_result.file.str)
        contami_prob_res_name = final_out_dir + "/Contamination_probability.tsv"
        contamination_probability_result.to_csv(contami_prob_res_name, sep="\t", index=False)
    else:
        contami_prob_res_name = final_out_dir + "/Contamination_probability.tsv"
        contamination_probability_result = pd.read_csv(contami_prob_res_name, sep="\t")

    # --------------------------------------------------------------------------
    #                        Purifying vMAGs and vSAGs
    # --------------------------------------------------------------------------
    purification_target_list = contamination_probability_result[contamination_probability_result['result_prediction']==1]['file'].tolist()
    for fasta_file in os.listdir(input_fasta_dir):
        file_path = os.path.join(input_fasta_dir, fasta_file)

        # parse the fasta file
        records = list(SeqIO.parse(file_path, 'fasta'))
        basename = os.path.splitext(fasta_file)[0]
        """
        if '_edit' in tmp_basename:
            basename = tmp_basename.replace('_edit', '')
        elif '_2000' in tmp_basename:
            basename = tmp_basename.replace('_2000', '')
            print(basename)
        elif '_5000' in tmp_basename:
            basename = tmp_basename.replace('_5000', '')
        else:
            basename = tmp_basename
            print(basename)
        """

        if basename in purification_target_list:
            longest_seq = max(records, key=len)
            # create a new file in the output directory
            output_file = os.path.join(out_fasta_dir, basename + "_contami.fasta")
            with open(output_file, 'w') as f:
                SeqIO.write(longest_seq, f, 'fasta')
        else:
            output_file = os.path.join(out_fasta_dir, basename + "_noncontami.fasta")
            shutil.copy(file_path, output_file)

if __name__ == "__main__":
    status = main()
    sys.exit(status)
