# Turvey and Cavalcanti, 2024

Repository for the paper: 

> Turvey, AK and Cavalcanti ARO (2024) Human disease-causing missense genetic variants are enriched in the evolutionarily ancient, conserved domains of the aminoacyl-tRNA synthetase proteins.

**REQUIREMENTS:**

You need python 3 installed with the following python packages:

- pandas
- plotly
- plotly-kaleido
- scipy
- numpy
- statsmodel
- scikit-bio

Additionally, you need the executables for Blast, MAFFT and DSSP in the path. You can install these from bioconda.

**FILES:**

The scripts necessary to replicate the analyses of the paper are `RunAllAARS.py` and `Summarize.py`. These scripts use the data in the `aaRS_domains.xlsx` file to load the domains from Guo et *al.*, and the files in the `hgmd` directory to load the pathogenic variant data from the HGMD database.

The Results folder has html files with the results of the analyses using a coverage of 5% (as in the paper) and 1%. All the results used in the paper are in these html files.

**INSTRUCTIONS:**

To repeat the calculations of the paper first run `RunAllAARS.py`, this script will download the data and perform all the calculations for each aaRS. It creates one directory for each aaRS. In each of this directories, among other files, it creates an html file with the results for that aaRS.

Then, run `Summarize.py`. This will create a `summary.html` file with a compilation of the results for all aaRSs.


