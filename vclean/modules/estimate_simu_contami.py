from Bio import SeqIO, SeqUtils

def cal_contamination(fasta):
    seq_len_dict = {}
    for record in SeqIO.parse(fasta, 'fasta'):
        # id_part = record.id
        desc_part = record.description
        desc_list = desc_part.split('_')[:-1]
        seq_name = '_'.join(desc_list)
        seq = record.seq
        seq_len = len(seq)
        if seq_name in seq_len_dict:
            seq_len_dict[seq_name] += seq_len
        else:
            seq_len_dict[seq_name] = seq_len

    print(seq_len_dict)
    if len(seq_len_dict) >= 2:
        result = 100 * ((sum(seq_len_dict.values()) -  max(seq_len_dict.values())) / sum(seq_len_dict.values()))
    else:
        result = 0

    return result
