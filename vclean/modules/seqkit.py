import subprocess as sp

def run_seqkit(fasta):
    seqkit_stats_command = sp.run(["seqkit", "stats", fasta], capture_output=True)
    sed_seqkit_command = sp.run(['sed', '-e', 's/,//g'], input=seqkit_stats_command.stdout, capture_output=True)
    res_seqkit = sp.run(['awk', 'NR>1{print $4,$5,$6,$7}'], input=sed_seqkit_command.stdout, capture_output=True).stdout.decode().split()
    return res_seqkit
