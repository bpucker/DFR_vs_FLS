#Python code to create 2D density heatmap with marginal histograms using any gene IDs of interest and expression data (from count table) 

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.colors as mcolors

count_table_path = '20211016_Vitis_vinifera_phytozome_tpms_clean.txt'
count_table = pd.read_csv(count_table_path, sep='\t', index_col=0)

gene_x = 'VIT_218s0001g03490.1' #replace with the actual gene ID you want to plot on the x-axis
gene_x_values = count_table.loc[genex]

gene_y = 'VIT_218s0001g12800.1' #Replace with the actual gene ID you want to plot on the y-axis
gene_y_values = count_table.loc[geney]

def marghist(x=None, y=None, cmap='gnuplot2_r', ylab=None, xlab=None, bins=25, figsize=(5, 4), snsbins=50):
    ax1 = sns.jointplot(x=gene_x_values, y=gene_y_values, color='#cb181d', marginal_kws=dict(bins=snsbins))
    ax1.fig.set_size_inches(figsize[0], figsize[1])
    ax1.ax_joint.cla()
    plt.sca(ax1.ax_joint)
    plt.hist2d(gene_x_values, gene_y_values, norm=mcolors.LogNorm(), cmap=cmap, bins=bins)
    sample_size_text = f'n = {len(gene_x_values)}'
    plt.text(0.5, 0.95, sample_size_text, transform=plt.gca().transAxes, ha='center', va='center', fontsize=10,
             bbox=dict(facecolor='white', alpha=0.8))
    plt.xlabel(r'$\it{' + gene_x + '}$ Expression (in TPMs)', fontsize=12)  
    plt.ylabel(r'$\it{' + gene_y + '}$ Expression (in TPMs)', fontsize=12)
    cbar_ax = ax1.fig.add_axes([1, 0.1, .03, .7])
    cb = plt.colorbar(cax=cbar_ax)
    cb.set_label(r'$\log_{10}$ density of points', fontsize=13)
    #plt.suptitle(r'$\it{' + gene_x + '}$ vs $\it{' + gene_y + '}$ Expression')

marghist()

plt.show()
