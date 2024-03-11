import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import matplotlib.colors as mcolors

# Replace 'your_count_table.csv' with the actual file path of your count table
count_table_path = 'F_N_merged_count_table.txt'

# Read the count table
count_table = pd.read_csv(count_table_path, sep='\t', index_col=0)

gene_x = 'FLS_F'
# Replace 'GeneX' with the actual gene name you want to plot on the x-axis
gene_x_values = count_table.loc[gene_x]

gene_y = 'DFR_N'
# Replace 'GeneY' with the actual gene name you want to plot on the y-axis
gene_y_values = count_table.loc[gene_y]

def makesweetgraph(x=None, y=None, cmap='gnuplot2_r', ylab=None, xlab=None, bins=25, figsize=(5, 4), snsbins=50):
    ax1 = sns.jointplot(x=gene_x_values, y=gene_y_values, color='#cb181d', marginal_kws=dict(bins=snsbins))
    ax1.fig.set_size_inches(figsize[0], figsize[1])
    ax1.ax_joint.cla()
    plt.sca(ax1.ax_joint)
    # Create a 2D histogram
    plt.hist2d(gene_x_values, gene_y_values, norm=mcolors.LogNorm(), cmap=cmap, bins=bins)
    # Add sample size (N) as text
    sample_size_text = f'n = {len(gene_x_values)}'
    plt.text(0.5, 0.95, sample_size_text, transform=plt.gca().transAxes, ha='center', va='center', fontsize=10,
             bbox=dict(facecolor='white', alpha=0.8))
    plt.xlabel('${FLS_\mathrm{F}}$ expression (in TPMs)', fontsize=12)
    plt.ylabel('${DFR_\mathrm{N}}$ Expression (in TPMs)', fontsize=12)
    # Add a colorbar
    cbar_ax = ax1.fig.add_axes([0.2, 0.02, 0.6, 0.02])
    cb = plt.colorbar(cax=cbar_ax, orientation='horizontal'
    cb.set_label(r'$\log_{10}$ density of points', fontsize=13)"""
    # Add label (a) at top left

# Call the function without sets
makesweetgraph()

plt.savefig("FFvsDN", bbox_inches='tight', dpi=1000)
plt.show()
