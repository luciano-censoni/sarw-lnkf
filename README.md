# sarw-lnkf
USAGE: python read-folding-rate.py

REQUIRES: bash utilities (ls, grep, xargs, rm, wget, gunzip) plus numpy, scipy, matplotlib


The script will create a "structures" folder under its working directory and will download over a hundred structures directly from the pdb repository. Average information will be calculated with and without redundant contacts for each structure, and bootstrapped correlations will be plotted. The entire process may take about 15-20 minutes.

Complete lists of proteins included in and excluded from the analysis will be generated automatically and written to the included_proteins.log and excluded_proteins.log files.



TODO: include source code for mypdb.py instead of compiled ".pyc" file.
