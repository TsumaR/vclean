import pandas as pd

def func_cate(x):
    if  x == 0:
        return 0
    elif x >0 and x <= 1:
        return 1
    elif x > 1 and x <= 2:
        return 2
    elif x >= 3:
        return 3

def table_reformat(tbl):
    result_tbl = tbl.copy()
    # float64の列をfloat8に変換
    float64_cols = result_tbl.select_dtypes(include='float64').columns
    result_tbl[float64_cols] = result_tbl[float64_cols].astype('float16')
    
    # result_tbl = result_tbl[result_tbl["num_seqs"]==1]
    result_tbl["min_contig_rate"] = 100 * (result_tbl["min_len"] / result_tbl["sum_len"])
    result_tbl = result_tbl.drop(['file', 'avg_len'], axis=1)
    
    # 1. overlap single copy like geneを0 or 1に
    result_tbl["overlap_sp_one_or_no"] = 0
    result_tbl["overlap_sp_one_or_no"].where(result_tbl["overlap_singl_like_gene"] == 0, 1 ,inplace=True)
    
    # redundancy
    result_tbl['redundant_cat'] = result_tbl['redundant_protein'].apply(func_cate)
    
    return result_tbl
