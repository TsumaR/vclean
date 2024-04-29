import pandas as pd
import matplotlib.pyplot as plt
import sklearn #機械学習のライブラリ
from sklearn.decomposition import PCA #主成分分析器
import numpy as np
import seaborn as sns
import re
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score
import warnings

#--------- require
## ./seqkit_stats_Caudo_doublet.txt
## ./quality_summary.tsv
## ./single_copy_features.tsv
## ./ALL_Caudo_sep_10000_summary/ALL_Caudo_sep_10000_virus_genes.tsv
## ./gc_table.tsv
## ./codon_usage.tsv
## ./tetra_calculates.txt



def read_data():
    seqkit = pd.read_table("./seqkit_stats_Caudo_doublet.txt", delim_whitespace=True)
    seqkit = seqkit.replace(',', '', regex=True)
    features_table = seqkit.copy()
    features_table["gene_count"]=0
    features_table["viral_genes"]=0
    features_table["host_genes"]=0
    features_table["sum_completeness"]=0
    features_table["sum_contamination"]=0
    features_table["count_medium"]=0
    features_table["count_high"]=0
    features_table["Status"] = np.where(features_table["file"].str.contains("doublet"), "Doublet", "Single")
    return feature_table

def checkv(input):
    #checkV
    checkv = pd.read_table("./quality_summary.tsv")
    checkv[["NC", "genome_id", "contig_no"]] = checkv.contig_id.str.split("_", expand=True)
    checkv_impo = checkv[["genome_id", "gene_count", "viral_genes", "host_genes", "completeness", "contamination"]].groupby("genome_id").sum()
    for index, row in checkv_impo.iterrows():
        features_table["gene_count"] = np.where(features_table["file"].str.contains(index), features_table["gene_count"]+row[0], features_table["gene_count"])
        features_table["viral_genes"] = np.where(features_table["file"].str.contains(index), features_table["viral_genes"]+row[1], features_table["viral_genes"])
        features_table["host_genes"] = np.where(features_table["file"].str.contains(index), features_table["host_genes"]+row[2], features_table["host_genes"])
        features_table["sum_completeness"] = np.where(features_table["file"].str.contains(index), features_table["sum_completeness"]+row[3], features_table["sum_completeness"])
        features_table["sum_contamination"] = np.where(features_table["file"].str.contains(index), features_table["sum_contamination"]+row[4], features_table["sum_contamination"])
    return feature table

def single_copy_list():
     single_copy = pd.read_table("./single_copy_features.tsv", index_col = 0)
     single_copy_list = list(single_copy.index.values)

     gene_table = pd.read_table("./ALL_Caudo_sep_10000_summary/ALL_Caudo_sep_10000_virus_genes.tsv").dropna(subset=["annotation_description"])
     gene_table = gene_table[gene_table["annotation_description"].isin(single_copy_list)]
     gene_table[["NC", "genome_id", "contig_no", "gene_no"]] = gene_table.gene.str.split("_", expand=True)
     each_geno_count = gene_table[["gene","genome_id", "annotation_description"]].groupby(["genome_id", "annotation_description"]).count()

     for s_gene in single_copy_list:
         features_table[s_gene] = 0

     for index, row in each_geno_count.iterrows():
         geno_id = index[0]
         gene_name = index[1]
         num_g = row["gene"]
         features_table[gene_name] = np.where(features_table["file"].str.contains(geno_id) , features_table[gene_name]+num_g, features_table[gene_name])
     return feature_table

def gc_values():
    gc_table = pd.read_table("./gc_table.tsv")
    gc_table[["NC", "genome_id", "contig_no"]] = gc_table["Sequence Id"].str.split("_", expand=True)
    gc_table = gc_table.drop(["NC", "contig_no"], axis=1)
    # gc_table_var = gc_table.groupby("genome_id").var()

    for gc in ["GC_contents", "CPG_contents", "GC_skew"]:
        features_table[gc] = 0
    
    for vSAG in list(features_table.file):
        ids = re.split("_|\.",vSAG)
        if len(ids) >= 8:
            id_1 =  ids[2]
            id_2 = ids[5]
            target_gc = gc_table[(gc_table["genome_id"].str.contains(id_1)) | (gc_table["genome_id"].str.contains(id_2))]
            features_table.loc[features_table["file"]==vSAG,"GC_contents"] = target_gc["GC_contents"].var()
            features_table.loc[features_table["file"]==vSAG,"CPG_contents"] = target_gc["CPG_contents"].var()
            features_table.loc[features_table["file"]==vSAG, "GC_skew"] = target_gc["GC_skew"].var()

        else:
            id_1 =  ids[2]
            target_gc = gc_table[(gc_table["genome_id"].str.contains(id_1))]
            features_table.loc[features_table["file"]==vSAG,"GC_contents"] = target_gc["GC_contents"].var()
            features_table.loc[features_table["file"]==vSAG,"CPG_contents"] = target_gc["CPG_contents"].var()
            features_table.loc[features_table["file"]==vSAG, "GC_skew"] = target_gc["GC_skew"].var()
    return feature_table


def codon_usage():
    codon_table = pd.read_table("./codon_usage.tsv")
    features_table["codon_var"] = 0
    for vSAG in list(features_table.file):
    ids = re.split("_|\.",vSAG)
    if len(ids) >= 8:
        id_1 =  ids[2]
        id_2 = ids[5]
        target_codon = codon_table[(codon_table["Unnamed: 0"].str.contains(id_1)) | (codon_table["Unnamed: 0"].str.contains(id_2))]
        pca = PCA()
        pca.fit(target_codon.iloc[:,1:])
        feature = pca.transform(target_codon.iloc[:,1:])
        features_table.loc[features_table["file"]==vSAG,"codon_var"] = feature[:, 0].var()

    else:
        id_1 =  ids[2]
        target_codon = codon_table[(codon_table["Unnamed: 0"].str.contains(id_1))]
        if len(target_codon) >= 2:
            pca = PCA()
            pca.fit(target_codon.iloc[:,1:])
            feature = pca.transform(target_codon.iloc[:,1:])
            features_table.loc[features_table["file"]==vSAG,"codon_var"] = feature[:, 0].var()
        else:
            features_table.loc[features_table["file"]==vSAG,"codon_var"] = 0
    return feature_table

def tetra_freq():
    tetra_table = pd.read_table("./tetra_calculates.txt")
    tetra_table["tetra_var"] = 0
    for vSAG in list(features_table.file):
    ids = re.split("_|\.",vSAG)
    if len(ids) >= 8:
        id_1 =  ids[2]
        id_2 = ids[5]
        target_codon = tetra_table[(tetra_table["Sequence Id"].str.contains(id_1)) | (tetra_table["Sequence Id"].str.contains(id_2))]
        pca = PCA()
        pca.fit(target_codon.iloc[:,1:])
        feature = pca.transform(target_codon.iloc[:,1:])
        features_table.loc[features_table["file"]==vSAG,"tetra_var"] = feature[:, 0].var()

    else:
        id_1 =  ids[2]
        target_codon = tetra_table[(tetra_table["Sequence Id"].str.contains(id_1))]
        if len(target_codon) >= 2:
            pca = PCA()
            pca.fit(target_codon.iloc[:,1:])
            feature = pca.transform(target_codon.iloc[:,1:])
            features_table.loc[features_table["file"]==vSAG,"tetra_var"] = feature[:, 0].var()
        else:
            features_table.loc[features_table["file"]==vSAG,"tetra_var"] = 0

    return feature_table

def main():
    



