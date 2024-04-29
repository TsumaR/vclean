#!/usr/bin/env python

import subprocess as sp

def download_database(db_dir):
    checkv_download_command = sp.run(["checkv", "download_database", db_dir])
    pfam_download_command = sp.run(["wget", "-P", db_dir, "https://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam35.0/Pfam-A.hmm.gz"])
    pfam_db_dir = db_dir + "/Pfam"
    mkdir_pfam = sp.run(["mkdir", "-p", pfam_db_dir])
    pfam_db = pfam_db_dir + "/Pfam-A.hmm.gz"
    unzipped_pfam_db = pfam_db_dir + "/Pfam-A.hmm"
    pfam_unzip = sp.run(["gunzip", pfam_db])
    hmmpress_command = sp.run(["hmmpress", unzipped_pfam_db])

if __name__ == "__main__":
    download_database(db_dir) 
