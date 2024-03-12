#Python code to create 2D density heatmap with marginal histograms using any gene IDs of interest and expression data (from count table) 
__usage__="""
                        python3 Coexp_plot.py
                        --count_table_path <full path to count table>
                        --gene_x <geneID to be plotted in X-axis enclosed in quotes, if isoforms exist, multiple comma-separated geneIDs could also be provided > 
                        --gene_y <geneID to be plotted in Y-axis enclosed in quotes, if isoforms exist, multiple comma-separated geneIDs could also be provided > 
                        --gene_x_name <name of gene to be plotted in X-axis>
                        --gene_y_name <name of gene to be plotted in Y-axis>
                        """

import os, sys
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.colors as mcolors

def main(arguments):
    count_table_path = arguments[arguments.index('--count_table_path') + 1]
    count_table = pd.read_csv(count_table_path, sep='\t', index_col=0)
    gene_x = arguments[arguments.index('--gene_x') + 1]
    gene_y = arguments[arguments.index('--gene_y') + 1]
    gene_x_name = arguments[arguments.index('--gene_x_name') + 1]
    gene_y_name = arguments[arguments.index('--gene_y_name') + 1]
    
    if ',' in gene_x:            
        gene_x = gene_x.split(',')
    else:
        gene_x = [gene_x]

    if ',' in gene_y:            
        gene_y = gene_y.split(',')
    else:
        gene_y = [gene_y]        

    gene_x_values = count_table.loc[gene_x].sum()
    gene_y_values = count_table.loc[gene_y].sum()

    marghist(x_values=gene_x_values, y_values=gene_y_values, gene_x_name=gene_x_name, gene_y_name=gene_y_name)
    
def marghist(x_values=None, y_values=None, cmap='gnuplot2_r', gene_x_name=None, gene_y_name=None, bins=25, figsize=(5, 4), snsbins=50):
    ax1 = sns.jointplot(x=x_values,y= y_values, color='#cb181d', marginal_kws=dict(bins=snsbins))
    ax1.fig.set_size_inches(figsize[0], figsize[1])
    ax1.ax_joint.cla()
    plt.sca(ax1.ax_joint)
    plt.hist2d(x_values, y_values, norm=mcolors.LogNorm(), cmap=cmap, bins=bins)
    sample_size_text = f'n = {len(x_values)}'
    plt.text(0.5, 0.95, sample_size_text, transform=plt.gca().transAxes, ha='center', va='center', fontsize=10,
             bbox=dict(facecolor='white', alpha=0.8))
    plt.xlabel(r'$\it{' + gene_x_name + '}$ Expression (in TPMs)', fontsize=12)  
    plt.ylabel(r'$\it{' + gene_y_name + '}$ Expression (in TPMs)', fontsize=12)
    cbar_ax = ax1.fig.add_axes([1, 0.1, .03, .7])
    cb = plt.colorbar(cax=cbar_ax)
    cb.set_label(r'$\log_{10}$ density of points', fontsize=13)
    plt.suptitle(r'$\it{' + gene_x_name + '}$ vs $\it{' + gene_y_name + '}$ Expression')

if __name__ == '__main__':
    if '--count_table_path' in sys.argv and '--gene_x' in sys.argv and '--gene_y' in sys.argv and '--gene_x_name' in sys.argv and '--gene_y_name' in sys.argv:
        main(sys.argv)
        plt.savefig("Co-expression_plot.png", bbox_inches='tight', dpi=600)
    else:
        sys.exit(__usage__)

