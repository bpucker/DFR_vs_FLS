# DFR vs FLS
This repository harbours all scripts used to study DFR and FLS and their contribution to the competing branches of the flavonoid biosynthesis.

## Create a Co-expression plot between any 2 genes of interest based on expression data
```
Usage:
python3 Coexp_plot.py --count_table_path <FILE> --gene_x <STR> --gene_y <STR> --gene_x_name <STR> --gene_y_name <STR>

Mandatory:
  --count_table_path  FILE  Expression file
  --gene_x            STR   Gene ID(s) to be plotted on X-axis
  --gene_y            STR   Gene ID(s) to be plotted on Y-axis 
  --gene_x_name       STR   Name of the gene to be plotted on X-axis 
  --gene_y_name       STR   Name of the gene to be plotted on Y-axis
  ```

``` --count_table_path``` specifies the full path to the expression data file
```--gene_x``` specifies the geneID to be plotted in X-axis. The geneID should be enclosed in quotes, e.g.,"GeneA". If isoforms of the genes exist, multiple comma-separated geneIDs could also be provided, e.g., "GeneA1,GeneA2,GeneA3.." 
```--gene_y``` specifies the geneID to be plotted in Y-axis. The geneID should be enclosed in quotes, e.g.,"GeneB". If isoforms of the genes exist, multiple comma-separated geneIDs could also be provided, e.g., "GeneB1,GeneB2,GeneB3.."
```--gene_x_name```specifies the name of the gene to be plotted in X-axis. 
```--gene_y_name```specifies the name of the gene to be plotted in Y-axis.

## Reference
Choudhary N. & Pucker B. (2023). Conserved amino acid residues and gene expression patterns associated with the substrate preferences of the competing enzymes FLS and DFR. bioRxiv 2023.11.05.565693; doi: [10.1101/2023.11.05.565693](https://doi.org/10.1101/2023.11.05.565693).
