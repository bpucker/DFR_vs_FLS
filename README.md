# DFR vs FLS
This repository harbours all scripts used to study DFR and FLS and their contribution to the competing branches of the flavonoid biosynthesis.

## Coexp_plot
'''
Usage:
python3 Coexp_plot.py --count_table_path <FILE> --gene_x <STR> --gene_y <STR> --gene_x_name <STR> --gene_y_name <STR>

Mandatory:
  --count_table_path  FILE  Expression file
  --gene_x            STR   Gene ID(s) to be plotted on X-axis
  --gene_y            STR   Gene ID(s) to be plotted on Y-axis 
  --gene_x_name       STR   Name of the gene to be plotted on X-axis 
  --gene_y_name       STR   Name of the gene to be plotted on Y-axis
  '''
  python3 Coexp_plot.py
                        --count_table_path <full path to count table>
                        --gene_x <geneID to be plotted in X-axis enclosed in quotes, if isoforms exist, multiple comma-separated geneIDs could also be provided > 
                        --gene_y <geneID to be plotted in Y-axis enclosed in quotes, if isoforms exist, multiple comma-separated geneIDs could also be provided > 
                        --gene_x_name <name of gene to be plotted in X-axis>
                        --gene_y_name <name of gene to be plotted in Y-axis>
## Reference
Choudhary N. & Pucker B. (2023). Conserved amino acid residues and gene expression patterns associated with the substrate preferences of the competing enzymes FLS and DFR. bioRxiv 2023.11.05.565693; doi: [10.1101/2023.11.05.565693](https://doi.org/10.1101/2023.11.05.565693).
