# vClean

vClean is a fully automated command-line pipeline for assessing the multiple virus-derived contamination risk of environmental viral genomes.

# Quickstart
Run the following commands to install vClean, download the database, and run the program.
Replace `</path/to/database>`, `<input_fasta_dir>` and `<output_dir>` with the correct paths.
```bash
conda install -c conda-forge -c bioconda vclean
vclean download_database </path/to/database>
export VCLEANDB=</path/to/database>
vclean run <input_fasta_dir> <output_dir> [options]
```

# Installation
You can install vClean as follows:
```bash
conda install -c conda-forge -c bioconda vclean python=3.9
```

*Specify python version to be 3.9*

# Database installation
You have to download the databases.
Please replace `</path/to/database>` with the desired path for downloading the database:
```bash
vclean download_database </path/to/database>
```

You'll need to use the `-d` flag or update the `VCLEANDB` environment variable to specify the database location:
```bash
export VCLEANDB=<path/to/database>
```

# Run
You can simply run vClean as follows:
```bash
vclean run <input_fasta_dir> <output_dir> [OPTIONS]
```

# Output files
In the output directory, two files are generated: `Contamination_probability.tsv` and `feature_table.tsv`. For sequences confirmed to be contaminated, new FASTA files containing only the longest contigs are stored in the `purified_fasta` directory within the output directory.
