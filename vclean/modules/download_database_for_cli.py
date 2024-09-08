#!/usr/bin/env python

import subprocess as sp
import time
import os

def fetch_arguments(parser):
    parser.set_defaults(func=download_database)
    parser.set_defaults(program="download_database")
    parser.add_argument(
        "destination",
        type=str,
        help="Directory where the database will be downloaded to.",
    )

def download_database(args):
    program_start = time.time()
    db_dir = args["destination"]
    if not os.path.exists(db_dir):
        os.makedirs(db_dir)

    checkv_download_command = sp.run(["checkv", "download_database", db_dir])
    pfam_db_dir = db_dir + "/Pfam"
    if not os.path.exists(pfam_db_dir):
        os.makedirs(pfam_db_dir)
    pfam_download_command = sp.run(["wget", "-P", pfam_db_dir, "https://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam35.0/Pfam-A.hmm.gz"])
    pfam_db = pfam_db_dir + "/Pfam-A.hmm.gz"
    unzipped_pfam_db = pfam_db_dir + "/Pfam-A.hmm"
    pfam_unzip = sp.run(["gunzip", pfam_db])
    hmmpress_command = sp.run(["hmmpress", unzipped_pfam_db])

if __name__ == "__main__":
    download_database(args)
