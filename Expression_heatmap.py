#Python code for creating expression heatmap
#Customized code used for creating AdditionalFile7 (https://doi.org/10.1101/2023.11.05.565693).

import numpy as np
import csv
import matplotlib
from matplotlib import pyplot as plt
import seaborn as sns
import pandas as pd
sns.set()
csv_reader=pd.read_csv('HEATmap_edit.csv')
asp = csv_reader.shape[0]/float(csv_reader.shape[1])
figw = 14
figh = figw*4

sp1_1=csv_reader.loc[:9]
sp1_2=csv_reader.loc[10:12]
sp2_1=csv_reader.loc[13:52]
sp2_2=csv_reader.loc[53:64]
sp3_1=csv_reader.loc[65:84]
sp3_2=csv_reader.loc[85:90]
sp4_1=csv_reader.loc[91:100]
sp4_2=csv_reader.loc[101:103]
sp5_1=csv_reader.loc[104:123]
sp5_2=csv_reader.loc[124:129]
sp6_1=csv_reader.loc[130:189]
sp6_2=csv_reader.loc[190:207]
sp7_1=csv_reader.loc[208:247]
sp7_2=csv_reader.loc[248:259]
sp8_1=csv_reader.loc[260:319]
sp8_2=csv_reader.loc[320:337]
sp9_1=csv_reader.loc[338:407]
sp9_2=csv_reader.loc[408:428]
sp10_1=csv_reader.loc[429:438]
sp10_2=csv_reader.loc[439:441]
sp11_1=csv_reader.loc[442:461]
sp11_2=csv_reader.loc[462:467]
sp12_1=csv_reader.loc[468:527]
sp12_2=csv_reader.loc[528:545]
sp13_1=csv_reader.loc[546:585]
sp13_2=csv_reader.loc[586:597]
sp14_1=csv_reader.loc[598:637]
sp14_2=csv_reader.loc[638:649]
sp15_1=csv_reader.loc[650:699]
sp15_2=csv_reader.loc[700:714]
sp16_1=csv_reader.loc[715:744]
sp16_2=csv_reader.loc[745:753]
sp17_1=csv_reader.loc[754:793]
sp17_2=csv_reader.loc[794:805]
sp18_1=csv_reader.loc[806:845]
sp18_2=csv_reader.loc[846:857]
sp19_1=csv_reader.loc[858:897]
sp19_2=csv_reader.loc[898:909]
sp20_1=csv_reader.loc[910:959]
sp20_2=csv_reader.loc[960:974]
sp21_1=csv_reader.loc[975:1034]
sp21_2=csv_reader.loc[1035:1052]
sp22_1=csv_reader.loc[1053:1112]
sp22_2=csv_reader.loc[1113:1130]
sp23_1=csv_reader.loc[1131:1160]
sp23_2=csv_reader.loc[1161:1169]
sp24_1=csv_reader.loc[1170:1189]
sp24_2=csv_reader.loc[1190:1195]
sp25_1=csv_reader.loc[1196:1235]
sp25_2=csv_reader.loc[1236:1247]
sp26_1=csv_reader.loc[1248:1257]
sp26_2=csv_reader.loc[1258:1260]
sp27_1=csv_reader.loc[1261:1290]
sp27_2=csv_reader.loc[1291:1299]
sp28_1=csv_reader.loc[1300:1329]
sp28_2=csv_reader.loc[1330:1338]
sp29_1=csv_reader.loc[1339:1388]
sp29_2=csv_reader.loc[1389:1403]
sp30_1=csv_reader.loc[1404:1443]
sp30_2=csv_reader.loc[1444:1455]
sp31_1=csv_reader.loc[1456:1465]
sp31_2=csv_reader.loc[1466:1468]
sp32_1=csv_reader.loc[1469:1528]
sp32_2=csv_reader.loc[1529:1546]
sp33_1=csv_reader.loc[1547:1576]
sp33_2=csv_reader.loc[1577:1585]
sp34_1=csv_reader.loc[1586:1635]
sp34_2=csv_reader.loc[1636:1650]
sp35_1=csv_reader.loc[1651:1690]
sp35_2=csv_reader.loc[1691:1702]
sp36_1=csv_reader.loc[1703:1742]
sp36_2=csv_reader.loc[1743:1754]
sp37_1=csv_reader.loc[1755:1774]
sp37_2=csv_reader.loc[1775:1780]
sp38_1=csv_reader.loc[1781:1820]
sp38_2=csv_reader.loc[1821:1832]
sp39_1=csv_reader.loc[1833:1862]
sp39_2=csv_reader.loc[1863:1871]
sp40_1=csv_reader.loc[1872:1901]
sp40_2=csv_reader.loc[1902:1910]
sp41_1=csv_reader.loc[1911:1960]
sp41_2=csv_reader.loc[1961:1975]
sp42_1=csv_reader.loc[1976:2005]
sp42_2=csv_reader.loc[2006:2014]
sp43_1=csv_reader.loc[2015:2044]
sp43_2=csv_reader.loc[2045:2053]

sp1_1 = sp1_1.pivot_table(index='Tissue', columns='Gene', values='Log2Value', sort=False)
sp1_2 = sp1_2.pivot_table(index='Tissue', columns='Gene', values='Log2Value', sort=False)
sp2_1 = sp2_1.pivot_table(index='Tissue', columns='Gene', values='Log2Value', sort=False)
sp2_2 = sp2_2.pivot_table(index='Tissue', columns='Gene', values='Log2Value', sort=False)
sp3_1 = sp3_1.pivot_table(index='Tissue', columns='Gene', values='Log2Value', sort=False)
sp3_2 = sp3_2.pivot_table(index='Tissue', columns='Gene', values='Log2Value', sort=False)
sp4_1 = sp4_1.pivot_table(index='Tissue', columns='Gene', values='Log2Value', sort=False)
sp4_2 = sp4_2.pivot_table(index='Tissue', columns='Gene', values='Log2Value', sort=False)
sp5_1 = sp5_1.pivot_table(index='Tissue', columns='Gene', values='Log2Value', sort=False)
sp5_2 = sp5_2.pivot_table(index='Tissue', columns='Gene', values='Log2Value', sort=False)
sp6_1 = sp6_1.pivot_table(index='Tissue', columns='Gene', values='Log2Value', sort=False)
sp6_2 = sp6_2.pivot_table(index='Tissue', columns='Gene', values='Log2Value', sort=False)
sp7_1 = sp7_1.pivot_table(index='Tissue', columns='Gene', values='Log2Value', sort=False)
sp7_2 = sp7_2.pivot_table(index='Tissue', columns='Gene', values='Log2Value', sort=False)
sp8_1 = sp8_1.pivot_table(index='Tissue', columns='Gene', values='Log2Value', sort=False)
sp8_2 = sp8_2.pivot_table(index='Tissue', columns='Gene', values='Log2Value', sort=False)
sp9_1 = sp9_1.pivot_table(index='Tissue', columns='Gene', values='Log2Value', sort=False)
sp9_2 = sp9_2.pivot_table(index='Tissue', columns='Gene', values='Log2Value', sort=False)
sp10_1 = sp10_1.pivot_table(index='Tissue', columns='Gene', values='Log2Value', sort=False)
sp10_2 = sp10_2.pivot_table(index='Tissue', columns='Gene', values='Log2Value', sort=False)
sp11_1 = sp11_1.pivot_table(index='Tissue', columns='Gene', values='Log2Value', sort=False)
sp11_2 = sp11_2.pivot_table(index='Tissue', columns='Gene', values='Log2Value', sort=False)
sp12_1 = sp12_1.pivot_table(index='Tissue', columns='Gene', values='Log2Value', sort=False)
sp12_2 = sp12_2.pivot_table(index='Tissue', columns='Gene', values='Log2Value', sort=False)
sp13_1 = sp13_1.pivot_table(index='Tissue', columns='Gene', values='Log2Value', sort=False)
sp13_2 = sp13_2.pivot_table(index='Tissue', columns='Gene', values='Log2Value', sort=False)
sp14_1 = sp14_1.pivot_table(index='Tissue', columns='Gene', values='Log2Value', sort=False)
sp14_2 = sp14_2.pivot_table(index='Tissue', columns='Gene', values='Log2Value', sort=False)
sp15_1 = sp15_1.pivot_table(index='Tissue', columns='Gene', values='Log2Value', sort=False)
sp15_2 = sp15_2.pivot_table(index='Tissue', columns='Gene', values='Log2Value', sort=False)
sp16_1 = sp16_1.pivot_table(index='Tissue', columns='Gene', values='Log2Value', sort=False)
sp16_2 = sp16_2.pivot_table(index='Tissue', columns='Gene', values='Log2Value', sort=False)
sp17_1 = sp17_1.pivot_table(index='Tissue', columns='Gene', values='Log2Value', sort=False)
sp17_2 = sp17_2.pivot_table(index='Tissue', columns='Gene', values='Log2Value', sort=False)
sp18_1 = sp18_1.pivot_table(index='Tissue', columns='Gene', values='Log2Value', sort=False)
sp18_2 = sp18_2.pivot_table(index='Tissue', columns='Gene', values='Log2Value', sort=False)
sp19_1 = sp19_1.pivot_table(index='Tissue', columns='Gene', values='Log2Value', sort=False)
sp19_2 = sp19_2.pivot_table(index='Tissue', columns='Gene', values='Log2Value', sort=False)
sp20_1 = sp20_1.pivot_table(index='Tissue', columns='Gene', values='Log2Value', sort=False)
sp20_2 = sp20_2.pivot_table(index='Tissue', columns='Gene', values='Log2Value', sort=False)
sp21_1 = sp21_1.pivot_table(index='Tissue', columns='Gene', values='Log2Value', sort=False)
sp21_2 = sp21_2.pivot_table(index='Tissue', columns='Gene', values='Log2Value', sort=False)
sp22_1 = sp22_1.pivot_table(index='Tissue', columns='Gene', values='Log2Value', sort=False)
sp22_2 = sp22_2.pivot_table(index='Tissue', columns='Gene', values='Log2Value', sort=False)
sp23_1 = sp23_1.pivot_table(index='Tissue', columns='Gene', values='Log2Value', sort=False)
sp23_2 = sp23_2.pivot_table(index='Tissue', columns='Gene', values='Log2Value', sort=False)
sp24_1 = sp24_1.pivot_table(index='Tissue', columns='Gene', values='Log2Value', sort=False)
sp24_2 = sp24_2.pivot_table(index='Tissue', columns='Gene', values='Log2Value', sort=False)
sp25_1 = sp25_1.pivot_table(index='Tissue', columns='Gene', values='Log2Value', sort=False)
sp25_2 = sp25_2.pivot_table(index='Tissue', columns='Gene', values='Log2Value', sort=False)
sp26_1 = sp26_1.pivot_table(index='Tissue', columns='Gene', values='Log2Value', sort=False)
sp26_2 = sp26_2.pivot_table(index='Tissue', columns='Gene', values='Log2Value', sort=False)
sp27_1 = sp27_1.pivot_table(index='Tissue', columns='Gene', values='Log2Value', sort=False)
sp27_2 = sp27_2.pivot_table(index='Tissue', columns='Gene', values='Log2Value', sort=False)
sp28_1 = sp28_1.pivot_table(index='Tissue', columns='Gene', values='Log2Value', sort=False)
sp28_2 = sp28_2.pivot_table(index='Tissue', columns='Gene', values='Log2Value', sort=False)
sp29_1 = sp29_1.pivot_table(index='Tissue', columns='Gene', values='Log2Value', sort=False)
sp29_2 = sp29_2.pivot_table(index='Tissue', columns='Gene', values='Log2Value', sort=False)
sp30_1 = sp30_1.pivot_table(index='Tissue', columns='Gene', values='Log2Value', sort=False)
sp30_2 = sp30_2.pivot_table(index='Tissue', columns='Gene', values='Log2Value', sort=False)
sp31_1 = sp31_1.pivot_table(index='Tissue', columns='Gene', values='Log2Value', sort=False)
sp31_2 = sp31_2.pivot_table(index='Tissue', columns='Gene', values='Log2Value', sort=False)
sp32_1 = sp32_1.pivot_table(index='Tissue', columns='Gene', values='Log2Value', sort=False)
sp32_2 = sp32_2.pivot_table(index='Tissue', columns='Gene', values='Log2Value', sort=False)
sp33_1 = sp33_1.pivot_table(index='Tissue', columns='Gene', values='Log2Value', sort=False)
sp33_2 = sp33_2.pivot_table(index='Tissue', columns='Gene', values='Log2Value', sort=False)
sp34_1 = sp34_1.pivot_table(index='Tissue', columns='Gene', values='Log2Value', sort=False)
sp34_2 = sp34_2.pivot_table(index='Tissue', columns='Gene', values='Log2Value', sort=False)
sp35_1 = sp35_1.pivot_table(index='Tissue', columns='Gene', values='Log2Value', sort=False)
sp35_2 = sp35_2.pivot_table(index='Tissue', columns='Gene', values='Log2Value', sort=False)
sp36_1 = sp36_1.pivot_table(index='Tissue', columns='Gene', values='Log2Value', sort=False)
sp36_2 = sp36_2.pivot_table(index='Tissue', columns='Gene', values='Log2Value', sort=False)
sp37_1 = sp37_1.pivot_table(index='Tissue', columns='Gene', values='Log2Value', sort=False)
sp37_2 = sp37_2.pivot_table(index='Tissue', columns='Gene', values='Log2Value', sort=False)
sp38_1 = sp38_1.pivot_table(index='Tissue', columns='Gene', values='Log2Value', sort=False)
sp38_2 = sp38_2.pivot_table(index='Tissue', columns='Gene', values='Log2Value', sort=False)
sp39_1 = sp39_1.pivot_table(index='Tissue', columns='Gene', values='Log2Value', sort=False)
sp39_2 = sp39_2.pivot_table(index='Tissue', columns='Gene', values='Log2Value', sort=False)
sp40_1 = sp40_1.pivot_table(index='Tissue', columns='Gene', values='Log2Value', sort=False)
sp40_2 = sp40_2.pivot_table(index='Tissue', columns='Gene', values='Log2Value', sort=False)
sp41_1 = sp41_1.pivot_table(index='Tissue', columns='Gene', values='Log2Value', sort=False)
sp41_2 = sp41_2.pivot_table(index='Tissue', columns='Gene', values='Log2Value', sort=False)
sp42_1 = sp42_1.pivot_table(index='Tissue', columns='Gene', values='Log2Value', sort=False)
sp42_2 = sp42_2.pivot_table(index='Tissue', columns='Gene', values='Log2Value', sort=False)
sp43_1 = sp43_1.pivot_table(index='Tissue', columns='Gene', values='Log2Value', sort=False)
sp43_2 = sp43_2.pivot_table(index='Tissue', columns='Gene', values='Log2Value', sort=False)

my_cmap1 = plt.get_cmap('YlGnBu').copy()
my_cmap1.set_over('whitesmoke')
my_cmap2 =plt.get_cmap('rocket_r').copy()
my_cmap2.set_over('whitesmoke')
norm1 = matplotlib.colors.Normalize(vmin= -7.212843061, vmax=10.01308325)
norm2 = matplotlib.colors.Normalize(vmin= -15.1397545, vmax=14.97349558)
                     
gridspec_kw = {"height_ratios":[1,4,2,1,2,6,4,6,7,1,2,6,4,4,5,3,4,4,4,5,6,6,3,2,4,1,3,3,5,4,1,6,3,5,4,4,2,4,3,3,5,3,3], "width_ratios":[10,3]}
heatmapkws1 = dict(annot_kws = {'size':8}, fmt =".2f", square=False, cbar=False, cmap = my_cmap1, clip_on=False, linewidths=0.5,vmin= -7.212843061, vmax=10.01308325) 
heatmapkws2 = dict(annot_kws = {'size':8}, fmt =".2f", square=False, cbar=False, cmap = my_cmap2, clip_on=False, linewidths=0.5,vmin= -15.1397545, vmax=14.97349558) 
tickskw =  dict(xticklabels=False, yticklabels=False)
Genes1 = ["$\\it{F3H}$", "$\\it{F3'H}$", "$\\it{F3'5'H}$", "$\\it{FLS_{H}}$", "$\\it{FLS_{F}}$", "$\\it{FLS_{Y}}$", "$\\it{DFR_{N}}$", "$\\it{DFR_{D}}$", "$\\it{DFR_{A}}$","$\\it{DFR*}$" ]                     
Genes2 = ["$\\it{F3'H/F3H}$","$\\it{F3'5'H/F3H}$","$\\it{DFR/FLS}$"]  
left = 0.07; right=0.87
bottom = 0.1; top = 0.9
fig, axes = plt.subplots(ncols=2, nrows=43, figsize=(figw, figh), gridspec_kw=gridspec_kw)
plt.subplots_adjust(left=left, right=right,bottom=bottom, top=top, wspace=0.02, hspace=0.1*2 )
labels1_1=  np.array([[0.732778,50.3783495,1.703475,90000,15.2797,90000,22.1472,90000,90000,90000]])
sp1_1=sns.heatmap(sp1_1, annot = labels1_1, ax=axes[0,0], xticklabels=Genes1, yticklabels=False, **heatmapkws1)
sp1_1.set(xlabel=None)
sp1_1.set_ylabel('$\\it{Nymphae\ colorata}$', rotation =0, ha='right', va='center')
sp1_1.xaxis.tick_top()
sp1_1.patch.set_linewidth(2)
sp1_1.patch.set_edgecolor('k')
for t in sp1_1.texts:
    if float(t.get_text())!=90000:
        t.set_text(t.get_text()) 
    else:
        t.set_text("")
labels1_2=  np.array([[68.74981,2.324680872,1.449452542]])
sp1_2=sns.heatmap(sp1_2, annot = labels1_2, ax=axes[0,1], xticklabels=Genes2, yticklabels=True, **heatmapkws2)
sp1_2.set(xlabel=None)
sp1_2.set(ylabel=None)
sp1_2.xaxis.tick_top()
sp1_2.set_yticklabels(sp1_2.get_yticklabels(), rotation = 0, fontsize = 9)
sp1_2.yaxis.tick_right()
sp1_2.patch.set_linewidth(2)
sp1_2.patch.set_edgecolor('k')
labels2_1=  np.array([[0.0896253,2.4102478,0.4284246,90000,0.0357204,90000,0.12393,90000,90000,90000],
                      [0.0290858,2.5374015,3.46656,90000,0.2664485,90000,0.01,90000,90000,90000],
                      [0.01,4.18242115,5.3017041,90000,0.09820815,90000,0.7589445,90000,90000,90000],
                      [0.0519369,0.0324031,1.0195912,90000,0.0324031,90000,0.795422,90000,90000,90000]])
sp2_1=sns.heatmap(sp2_1, annot=labels2_1, ax=axes[1,0], xticklabels=False, yticklabels=False, **heatmapkws1)
sp2_1.set(xlabel=None)
sp2_1.set_ylabel('$\\it{Thinopyrum\ intermedium}$', rotation =0, ha='right', va='center')
sp2_1.patch.set_linewidth(2)
sp2_1.patch.set_edgecolor('k')
for t in sp2_1.texts:
    if float(t.get_text())!=90000:
        t.set_text(t.get_text()) 
    else:
        t.set_text("")
labels2_2=  np.array([[26.89249353,4.780174794,3.469558012],
                      [87.23849782,119.1839317,0.037530705],
                      [418.242115,530.17041,7.727917693],
                      [0.623893609,19.63134496,24.54771303]])
sp2_2=sns.heatmap(sp2_2, annot=labels2_2, ax=axes[1,1], xticklabels=False, yticklabels=True, **heatmapkws2)
sp2_2.set(xlabel=None)
sp2_2.set(ylabel=None)
sp2_2.set_yticklabels(sp2_2.get_yticklabels(), rotation = 0, fontsize = 9)
sp2_2.yaxis.tick_right()
sp2_2.patch.set_linewidth(2)
sp2_2.patch.set_edgecolor('k')
labels3_1=  np.array([[1.418053,41.012829,17.848,90000,7.01433,90000,0.8881035,90000,90000,90000],
                      [0.2906849,48.455583,15.2617,90000,5.42678,90000,0.416928,90000,90000,90000]])
sp3_1=sns.heatmap(sp3_1, annot=labels3_1, ax=axes[2,0], xticklabels=False, yticklabels=False, **heatmapkws1)
sp3_1.set(xlabel=None)
sp3_1.set_ylabel('$\\it{Miscanthus\ sinensis}$', rotation =0, ha='right', va='center')
sp3_1.patch.set_linewidth(2)
sp3_1.patch.set_edgecolor('k')
for t in sp3_1.texts:
    if float(t.get_text())!=90000:
        t.set_text(t.get_text()) 
    else:
        t.set_text("")
labels3_2=  np.array([[28.92192958,12.58627146,0.126612734],[166.6945307,52.50255517,0.07682788]])
sp3_2=sns.heatmap(sp3_2, annot=labels3_2, ax=axes[2,1], xticklabels=False, yticklabels=True, **heatmapkws2)
sp3_2.set(xlabel=None)
sp3_2.set(ylabel=None)
sp3_2.set_yticklabels(sp3_2.get_yticklabels(), rotation = 0, fontsize = 9)
sp3_2.yaxis.tick_right()
sp3_2.patch.set_linewidth(2)
sp3_2.patch.set_edgecolor('k')
labels4_1=  np.array([[0.278629,18.6872617,0.0847113,90000,0.01,90000,1.749275,90000,90000,90000]])
sp4_1=sns.heatmap(sp4_1, annot=labels4_1, ax=axes[3,0], xticklabels=False, yticklabels=False, **heatmapkws1)
sp4_1.set(xlabel=None)
sp4_1.set_ylabel('$\\it{Panicum\ hallii}$', rotation =0, ha='right', va='center')
sp4_1.patch.set_linewidth(2)
sp4_1.patch.set_edgecolor('k')
for t in sp4_1.texts:
    if float(t.get_text())!=90000:
        t.set_text(t.get_text()) 
    else:
        t.set_text("")
labels4_2=  np.array([[67.06861705,0.304029013,174.9275]])
sp4_2=sns.heatmap(sp4_2, annot=labels4_2, ax=axes[3,1], xticklabels=False, yticklabels=True, **heatmapkws2)
sp4_2.set(xlabel=None)
sp4_2.set(ylabel=None)
sp4_2.set_yticklabels(sp4_2.get_yticklabels(), rotation = 0, fontsize = 9)
sp4_2.yaxis.tick_right()
sp4_2.patch.set_linewidth(2)
sp4_2.patch.set_edgecolor('k')
labels5_1=  np.array([[0.571719,5.399874,1.17180145,90000,4.575075,90000,0.1815775,90000,90000,90000],
[1.990584,20.9840134,2.371787,90000,3.933433,90000,1.66488,90000,90000,90000]])
sp5_1=sns.heatmap(sp5_1, annot=labels5_1, ax=axes[4,0], xticklabels=False, yticklabels=False, **heatmapkws1)
sp5_1.set(xlabel=None)
sp5_1.set_ylabel('$\\it{Panicum\ virgatum}$', rotation =0, ha='right', va='center')
sp5_1.patch.set_linewidth(2)
sp5_1.patch.set_edgecolor('k')
for t in sp5_1.texts:
    if float(t.get_text())!=90000:
        t.set_text(t.get_text()) 
    else:
        t.set_text("")
labels5_2=  np.array([[9.444979089,2.049610823,0.03968842],
[10.54163673,0.113028283,0.423263851]])
sp5_2=sns.heatmap(sp5_2, annot=labels5_2, ax=axes[4,1], xticklabels=False, yticklabels=True, **heatmapkws2)
sp5_2.set(xlabel=None)
sp5_2.set(ylabel=None)
sp5_2.set_yticklabels(sp5_2.get_yticklabels(), rotation = 0, fontsize = 9)
sp5_2.yaxis.tick_right()
sp5_2.patch.set_linewidth(2)
sp5_2.patch.set_edgecolor('k')
labels6_1=  np.array([[0.01,7.492232,0.509888,90000,0.42166,90000,0.05584305,90000,90000,90000],
[0.01,0.26995695,0.05690385,90000,0.01,90000,0.01,90000,90000,90000],
[0.01,33.776539,11.15655,90000,1.01752,90000,0.1279165,90000,90000,90000],
[0.00956205,32.3254447,6.313875,90000,0.241437,90000,0.08442785,90000,90000,90000],
[0.503351,0.1797557,0.066041,90000,0.793338,90000,0.01,90000,90000,90000],
[12.74151,1.9840603,0.6288275,90000,1.38635,90000,3.935844395,90000,90000,90000]])
sp6_1=sns.heatmap(sp6_1, annot=labels6_1, ax=axes[5,0], xticklabels=False, yticklabels=False, **heatmapkws1)
sp6_1.set(xlabel=None)
sp6_1.set_ylabel('$\\it{Sorghum\ bicolor}$', rotation =0, ha='right', va='center')
sp6_1.patch.set_linewidth(2)
sp6_1.patch.set_edgecolor('k')
for t in sp6_1.texts:
    if float(t.get_text())!=90000:
        t.set_text(t.get_text()) 
    else:
        t.set_text("")
labels6_2=  np.array([[749.2232,50.9888,0.132436205],
[26.995695,5.690385,1],
[3377.6539,1115.655,0.125713991],
[3380.597748,660.305583,0.349688946],
[0.357117995,0.13120268,0.01],
[0.155716261,0.049352667,2.838997652]])
sp6_2=sns.heatmap(sp6_2, annot=labels6_2, ax=axes[5,1], xticklabels=False, yticklabels=True, **heatmapkws2)
sp6_2.set(xlabel=None)
sp6_2.set(ylabel=None)
sp6_2.set_yticklabels(sp6_2.get_yticklabels(), rotation = 0, fontsize = 9)
sp6_2.yaxis.tick_right()
sp6_2.patch.set_linewidth(2)
sp6_2.patch.set_edgecolor('k')
labels7_1=  np.array([[66.42925,0.2621075,17.5243635,90000,49.85990535,0.01,16.3295921,90000,90000,0.01],
[68.7863,0.232368,28.4176875,90000,13.157935,0.01,20.74989,90000,90000,0.18149],
[77.566745,2.19211,13.2336277,90000,0.9508385,0.01,39.74412,90000,90000,0.01],
[18.32675,2.65393,3.18053,90000,0.691102,0.01,5.25418,90000,90000,0.01]])
sp7_1=sns.heatmap(sp7_1, annot=labels7_1, ax=axes[6,0], xticklabels=False, yticklabels=False, **heatmapkws1)
sp7_1.set(xlabel=None)
sp7_1.set_ylabel('$\\it{Musa\ acuminata}$', rotation =0, ha='right', va='center')
sp7_1.patch.set_linewidth(2)
sp7_1.patch.set_edgecolor('k')
for t in sp7_1.texts:
    if float(t.get_text())!=90000:
        t.set_text(t.get_text()) 
    else:
        t.set_text("")
labels7_2=  np.array([[0.003945664,0.263804928,0.327509488543383],
[0.003378115,0.413130049,1.57698681442035],
[0.028260951,0.170609553,41.7990226521118],
[0.144811819,0.173545773,7.60261148137323]])
sp7_2=sns.heatmap(sp7_2, annot=labels7_2, ax=axes[6,1], xticklabels=False, yticklabels=True, **heatmapkws2)
sp7_2.set(xlabel=None)
sp7_2.set(ylabel=None)
sp7_2.set_yticklabels(sp7_2.get_yticklabels(), rotation = 0, fontsize = 9)
sp7_2.yaxis.tick_right()
sp7_2.patch.set_linewidth(2)
sp7_2.patch.set_edgecolor('k')
labels8_1=  np.array([[226.5593,14.4211045,1.98503,63.84413,90000,0.01,10.350351,1.885122,90000,5.72332],
[547.2476,518.877622,22.1176,152.757202,90000,0.0487237,43.7443097,47.294473,90000,3.01721],
[43.516534,35.2312091,3.47304,10.93851,90000,0.326098,19.02512,2.9529583,90000,34.8814],
[1021.396355,60.7064915,6.67513,389.18105,90000,0.0237557,11.778815,10.90394,90000,17.02645],
[40.55431,42.258802,0.0782259,82.98398,90000,0.545894,13.05184,0.235006,90000,46.9824],
[629.4126,54.5288463,106.903,74.44845,90000,0.01,178.24828,110.55305,90000,6.31519]])
sp8_1=sns.heatmap(sp8_1, annot=labels8_1, ax=axes[7,0], xticklabels=False, yticklabels=False, **heatmapkws1)
sp8_1.set(xlabel=None)
sp8_1.set_ylabel('$\\it{Lotus\ japonicus}$', rotation =0, ha='right', va='center')
sp8_1.patch.set_linewidth(2)
sp8_1.patch.set_edgecolor('k')
for t in sp8_1.texts:
    if float(t.get_text())!=90000:
        t.set_text(t.get_text()) 
    else:
        t.set_text("")
labels8_2=  np.array([[0.063652671,0.008761635,0.283142601],
                      [0.94815879,0.040416075,0.606868708],
                      [0.809605128,0.079809665,6.336549285],
                      [0.059434803,0.006535298,0.113041448],
                      [1.042029861,0.001928917,0.719043803],
                      [0.086634501,0.169845662,3.902178084]])
sp8_2=sns.heatmap(sp8_2, annot=labels8_2, ax=axes[7,1], xticklabels=False, yticklabels=True, **heatmapkws2)
sp8_2.set(xlabel=None)
sp8_2.set(ylabel=None)
sp8_2.set_yticklabels(sp8_2.get_yticklabels(), rotation = 0, fontsize = 9)
sp8_2.yaxis.tick_right()
sp8_2.patch.set_linewidth(2)
sp8_2.patch.set_edgecolor('k')
labels9_1=  np.array([[38.6536,28.3259906,0.0434196,8.0633,90000,90000,0.01,0.530377,90000,90000],
                      [67.7722,24.1328,74.3371,31.2533,90000,90000,3.8472,5.21114,90000,90000],
                      [4.08456,1.3793398,2.13063,0.01,90000,90000,0.9160325,2.56425,90000,90000],
                      [47.4506,9.072005,2.74186,10.1912,90000,90000,1.410354,2.31907,90000,90000],
                      [11.42165,13.80390965,2.66433,3.326665,90000,90000,1.3505865,2.028045,90000,90000],
                      [1.88887,14.8589,0.133573,9.43833,90000,90000,0.0709902,2.028045,90000,90000],
                      [33.84095,21.15542035,0.329265,9.40496,90000,90000,26.40028,11.24386,90000,90000]])
sp9_1=sns.heatmap(sp9_1, annot=labels9_1, ax=axes[8,0], xticklabels=False, yticklabels=False, **heatmapkws1)
sp9_1.set(xlabel=None)
sp9_1.set_ylabel('$\\it{Vigna\ unguiculata}$', rotation =0, ha='right', va='center')
sp9_1.patch.set_linewidth(2)
sp9_1.patch.set_edgecolor('k')
for t in sp9_1.texts:
    if float(t.get_text())!=90000:
        t.set_text(t.get_text()) 
    else:
        t.set_text("")
labels9_2=  np.array([[0.732816364,0.0011233,0.081426835],
                      [0.356087009,1.096867152,0.289836273],
                      [0.337696055,0.521630237,344.6408],
                      [0.191188415,0.057783463,0.35952616],
                      [1.208574037,0.233270149,1.0156212],
                      [7.866555136,0.070715825,0.068754324],
                      [0.625142626,0.00972978,4.097379468]])
sp9_2=sns.heatmap(sp9_2, annot=labels9_2, ax=axes[8,1], xticklabels=False, yticklabels=True, **heatmapkws2)
sp9_2.set(xlabel=None)
sp9_2.set(ylabel=None)
sp9_2.set_yticklabels(sp9_2.get_yticklabels(), rotation = 0, fontsize = 9)
sp9_2.yaxis.tick_right()
sp9_2.patch.set_linewidth(2)
sp9_2.patch.set_edgecolor('k')
labels10_1=  np.array([[11.0235,33.48837,1.43833,38.3136,90000,90000,0.496657,90000,90000,90000]])
sp10_1=sns.heatmap(sp10_1, annot=labels10_1, ax=axes[9,0], xticklabels=False, yticklabels=False, **heatmapkws1)
sp10_1.set(xlabel=None)
sp10_1.set_ylabel('$\\it{Phaseolus\ acutifolius}$', rotation =0, ha='right', va='center')
sp10_1.patch.set_linewidth(2)
sp10_1.patch.set_edgecolor('k')
for t in sp10_1.texts:
    if float(t.get_text())!=90000:
        t.set_text(t.get_text()) 
    else:
        t.set_text("")
labels10_2=  np.array([[3.037907198,0.130478523,0.012962943]])
sp10_2=sns.heatmap(sp10_2, annot=labels10_2, ax=axes[9,1], xticklabels=False, yticklabels=True, **heatmapkws2)
sp10_2.set(xlabel=None)
sp10_2.set(ylabel=None)
sp10_2.set_yticklabels(sp10_2.get_yticklabels(), rotation = 0, fontsize = 9)
sp10_2.yaxis.tick_right()
sp10_2.patch.set_linewidth(2)
sp10_2.patch.set_edgecolor('k')
labels11_1=  np.array([[124.877117,34.662838,0.01,0.01,90000,90000,0.01,0.03089265,90000,90000],
				[34.204756,14.33875,0.01,1.41444,90000,90000,11.2506,6.46569,90000,90000]])
sp11_1=sns.heatmap(sp11_1, annot=labels11_1, ax=axes[10,0], xticklabels=False, yticklabels=False, **heatmapkws1)
sp11_1.set(xlabel=None)
sp11_1.set_ylabel('$\\it{Medicago\ truncatula}$', rotation =0, ha='right', va='center')
sp11_1.patch.set_linewidth(2)
sp11_1.patch.set_edgecolor('k')
for t in sp11_1.texts:
    if float(t.get_text())!=90000:
        t.set_text(t.get_text()) 
    else:
        t.set_text("")
labels11_2=  np.array([[0.277575579,0.0000801,21.2086],
[0.4192034,0.000292357,11.0272546]])
sp11_2=sns.heatmap(sp11_2, annot=labels11_2, ax=axes[10,1], xticklabels=False, yticklabels=True, **heatmapkws2)
sp11_2.set(xlabel=None)
sp11_2.set(ylabel=None)
sp11_2.set_yticklabels(sp11_2.get_yticklabels(), rotation = 0, fontsize = 9)
sp11_2.yaxis.tick_right()
sp11_2.patch.set_linewidth(2)
sp11_2.patch.set_edgecolor('k')
labels12_1=  np.array([[30.91066645,15.98782105,0.01,35.83025,90000,90000,90000,0.196118501,90000,90000],
[524.659013,18.3534731,0.088883,172.13893,90000,90000,90000,29.371941,90000,90000],
[17.9691,17.6109128,0.01,0.0512367,90000,90000,90000,0.01,90000,90000],
[115.7282202,37.3727897,0.01,27.0957,90000,90000,90000,33.913765,90000,90000],
[5.6655774,7.462398,0.01,2.55894,90000,90000,90000,0.205014,90000,90000],
[116.0678988,44.68872115,0.01,22.840943,90000,90000,90000,1.66606827,90000,90000]])
sp12_1=sns.heatmap(sp12_1, annot=labels12_1, ax=axes[11,0], xticklabels=False, yticklabels=False, **heatmapkws1)
sp12_1.set(xlabel=None)
sp12_1.set_ylabel('$\\it{Glycine\ max}$', rotation =0, ha='right', va='center')
sp12_1.patch.set_linewidth(2)
sp12_1.patch.set_edgecolor('k')
for t in sp12_1.texts:
    if float(t.get_text())!=90000:
        t.set_text(t.get_text()) 
    else:
        t.set_text("")
labels12_2=  np.array([[0.517226669,0.000323513,0.005473545],
[0.034981717,0.000169411,0.170629276],
[0.980066492,0.000556511,0.195172601],
[0.322935837,0.0000864,1.251629041],
[1.317146951,0.001765045,0.080116767],
[0.385022229,0.0000862,0.072942184]])
sp12_2=sns.heatmap(sp12_2, annot=labels12_2, ax=axes[11,1], xticklabels=False, yticklabels=True, **heatmapkws2)
sp12_2.set(xlabel=None)
sp12_2.set(ylabel=None)
sp12_2.set_yticklabels(sp12_2.get_yticklabels(), rotation = 0, fontsize = 9)
sp12_2.yaxis.tick_right()
sp12_2.patch.set_linewidth(2)
sp12_2.patch.set_edgecolor('k')
labels13_1=  np.array([[19.442487,13.9977839,90000,133.164,90000,90000,0.01,0.304542,90000,90000],
[0.3388683,0.9717825,90000,0.880537,90000,90000,0.0208702,0.0591913,90000,90000],
[16.9584094,9.92222995,90000,100.87495,90000,90000,0.04213515,0.265615,90000,90000],
[121.65171,0.1279924,90000,21.3006,90000,90000,0.01,0.0476442,90000,90000]])
sp13_1=sns.heatmap(sp13_1, annot=labels13_1, ax=axes[12,0], xticklabels=False, yticklabels=False, **heatmapkws1)
sp13_1.set(xlabel=None)
sp13_1.set_ylabel('$\\it{Cicer\ arietinum}$', rotation =0, ha='right', va='center')
sp13_1.patch.set_linewidth(2)
sp13_1.patch.set_edgecolor('k')
for t in sp13_1.texts:
    if float(t.get_text())!=90000:
        t.set_text(t.get_text()) 
    else:
        t.set_text("")
labels13_2=  np.array([[0.719958506,90000,0.003514133],
[2.867729144,90000,0.126391168],
[0.585092016,90000,0.002889873],
[0.001052122,90000,0.00247286]])
sp13_2=sns.heatmap(sp13_2, annot=labels13_2, ax=axes[12,1], xticklabels=False, yticklabels=True, **heatmapkws2)
sp13_2.set(xlabel=None)
sp13_2.set(ylabel=None)
sp13_2.set_yticklabels(sp13_2.get_yticklabels(), rotation = 0, fontsize = 9)
sp13_2.yaxis.tick_right()
sp13_2.patch.set_linewidth(2)
sp13_2.patch.set_edgecolor('k')
for t in sp13_2.texts:
    if float(t.get_text())!=90000:
        t.set_text(t.get_text()) 
    else:
        t.set_text("")
labels14_1=  np.array([[165.514,91.69378815,90000,13.7035315,90000,90000,42.0104,90000,90000,90000],
[296.912,94.239663,90000,0.01,90000,90000,321.715,90000,90000,90000],
[12.5581,43.869782,90000,0.300926,90000,90000,19.254,90000,90000,90000],
[27.29185,177.7760985,90000,0.01,90000,90000,10.8067,90000,90000,90000]])
sp14_1=sns.heatmap(sp14_1, annot=labels14_1, ax=axes[13,0], xticklabels=False, yticklabels=False, **heatmapkws1)
sp14_1.set(xlabel=None)
sp14_1.set_ylabel('$\\it{Prunus\ persica}$', rotation =0, ha='right', va='center')
sp14_1.patch.set_linewidth(2)
sp14_1.patch.set_edgecolor('k')
for t in sp14_1.texts:
    if float(t.get_text())!=90000:
        t.set_text(t.get_text()) 
    else:
        t.set_text("")
labels14_2=  np.array([[0.553994152,90000,3.065662308],
[0.317399307,90000,32171.5],
[3.49334549,90000,63.98250733],
[6.513889623,90000,1080.67]])
sp14_2=sns.heatmap(sp14_2, annot=labels14_2, ax=axes[13,1], xticklabels=False, yticklabels=True, **heatmapkws2)
sp14_2.set(xlabel=None)
sp14_2.set(ylabel=None)
sp14_2.set_yticklabels(sp14_2.get_yticklabels(), rotation = 0, fontsize = 9)
sp14_2.yaxis.tick_right()
sp14_2.patch.set_linewidth(2)
sp14_2.patch.set_edgecolor('k')
for t in sp14_2.texts:
    if float(t.get_text())!=90000:
        t.set_text(t.get_text()) 
    else:
        t.set_text("")
labels15_1=  np.array([[47.1647,24.272927,90000,38.6557,90000,90000,32.073,90000,23.6061,90000],
[14.2639,9.842902,90000,0.0686805,90000,90000,6.91401,90000,0.238206,90000],
[338.4395,143.274662,90000,0.637255,90000,90000,183.89705,90000,51.361286,90000],
[41.94825,10.35023735,90000,0.03593415,90000,90000,7.460615,90000,7.639399,90000],
[113.401,3.7734375,90000,0.05038705,90000,90000,77.84895,90000,19.985591,90000]])
sp15_1=sns.heatmap(sp15_1, annot=labels15_1, ax=axes[14,0], xticklabels=False, yticklabels=False, **heatmapkws1)
sp15_1.set(xlabel=None)
sp15_1.set_ylabel('$\\it{Fragaria\ vesca}$', rotation =0, ha='right', va='center')
sp15_1.patch.set_linewidth(2)
sp15_1.patch.set_edgecolor('k')
for t in sp15_1.texts:
    if float(t.get_text())!=90000:
        t.set_text(t.get_text()) 
    else:
        t.set_text("")
labels15_2=  np.array([[0.514641819,90000,1.562098733],
[0.690056857,90000,103.2078654],
[0.423339067,90000,388.0794807],
[0.246738239,90000,420.6277872],
[0.03327517,90000,1952.267835]])
sp15_2=sns.heatmap(sp15_2, annot=labels15_2, ax=axes[14,1], xticklabels=False, yticklabels=True, **heatmapkws2)
sp15_2.set(xlabel=None)
sp15_2.set(ylabel=None)
sp15_2.set_yticklabels(sp15_2.get_yticklabels(), rotation = 0, fontsize = 9)
sp15_2.yaxis.tick_right()
sp15_2.patch.set_linewidth(2)
sp15_2.patch.set_edgecolor('k')
for t in sp15_2.texts:
    if float(t.get_text())!=90000:
        t.set_text(t.get_text()) 
    else:
        t.set_text("")
labels16_1=  np.array([[119.1996,87.248558,90000,12.19211,90000,90000,66.3186,90000,3.754957,90000],
[99.7133,56.973234,90000,134.94175,90000,90000,39.96585,90000,6.804285,90000],
[342.4582,11.75670195,90000,0.897624,90000,90000,135.42978,90000,183.029,90000]])
sp16_1=sns.heatmap(sp16_1, annot=labels16_1, ax=axes[15,0], xticklabels=False, yticklabels=False, **heatmapkws1)
sp16_1.set(xlabel=None)
sp16_1.set_ylabel('$\\it{Fragaria\ x\ ananassa}$', rotation =0, ha='right', va='center')
sp16_1.patch.set_linewidth(2)
sp16_1.patch.set_edgecolor('k')
for t in sp16_1.texts:
    if float(t.get_text())!=90000:
        t.set_text(t.get_text()) 
    else:
        t.set_text("")
labels16_2=  np.array([[0.731953446,90000,5.852461715],
[0.571370459,90000,0.363613633],
[0.034330327,90000,348.9487413]])
sp16_2=sns.heatmap(sp16_2, annot=labels16_2, ax=axes[15,1], xticklabels=False, yticklabels=True, **heatmapkws2)
sp16_2.set(xlabel=None)
sp16_2.set(ylabel=None)
sp16_2.set_yticklabels(sp16_2.get_yticklabels(), rotation = 0, fontsize = 9)
sp16_2.yaxis.tick_right()
sp16_2.patch.set_linewidth(2)
sp16_2.patch.set_edgecolor('k')
for t in sp16_2.texts:
    if float(t.get_text())!=90000:
        t.set_text(t.get_text()) 
    else:
        t.set_text("")
labels17_1=  np.array([[110.3643,91.1788,3.23231,6.673394,90000,90000,107.213,3.67119,90000,90000],
[256.62321,49.956875,11.1402,110.707969,90000,90000,138.6815,3.46391,90000,90000],
[48.2015,32.173077,1.13857,0.01,90000,90000,107.156,0.233158,90000,90000],
[226.5593,14.4211045,1.98503,63.84413,90000,90000,20.2949,0.214822,90000,90000]])
sp17_1=sns.heatmap(sp17_1, annot=labels17_1, ax=axes[16,0], xticklabels=False, yticklabels=False, **heatmapkws1)
sp17_1.set(xlabel=None)
sp17_1.set_ylabel('$\\it{Carya\ illinoinensis}$', rotation =0, ha='right', va='center')
sp17_1.patch.set_linewidth(2)
sp17_1.patch.set_edgecolor('k')
for t in sp17_1.texts:
    if float(t.get_text())!=90000:
        t.set_text(t.get_text()) 
    else:
        t.set_text("")
labels17_2=  np.array([[0.826162083,0.029287641,16.72492288],
[0.194670135,0.043410727,1.277075998],
[0.667470452,0.023621049,10736.4683],
[0.063652671,0.008761635,0.283142601]])
sp17_2=sns.heatmap(sp17_2, annot=labels17_2, ax=axes[16,1], xticklabels=False, yticklabels=True, **heatmapkws2)
sp17_2.set(xlabel=None)
sp17_2.set(ylabel=None)
sp17_2.set_yticklabels(sp17_2.get_yticklabels(), rotation = 0, fontsize = 9)
sp17_2.yaxis.tick_right()
sp17_2.patch.set_linewidth(2)
sp17_2.patch.set_edgecolor('k')
labels18_1=  np.array([[53.812443,44.24694807,1.4295,24.0105,90000,90000,90000,8.29622,90000,90000],
[175.21075,78.0510173,9.45901,59.5349,90000,90000,90000,36.378,90000,90000],
[243.5609,34.0190976,0.511097,0.0268892,90000,90000,90000,27.0697,90000,90000],
[11.161247,25.180356,0.01,1.14284,90000,90000,90000,4.21319,90000,90000]])
sp18_1=sns.heatmap(sp18_1, annot=labels18_1, ax=axes[17,0], xticklabels=False, yticklabels=False, **heatmapkws1)
sp18_1.set(xlabel=None)
sp18_1.set_ylabel('$\\it{Populus\ trichocarpa}$', rotation =0, ha='right', va='center')
sp18_1.patch.set_linewidth(2)
sp18_1.patch.set_edgecolor('k')
for t in sp18_1.texts:
    if float(t.get_text())!=90000:
        t.set_text(t.get_text()) 
    else:
        t.set_text("")
labels18_2=  np.array([[0.822243808,0.026564488,0.345524666],
[0.445469341,0.053986471,0.611036552],
[0.139673887,0.002098436,1006.712732],
[2.256052214,0.000895957,3.686596549]])
sp18_2=sns.heatmap(sp18_2, annot=labels18_2, ax=axes[17,1], xticklabels=False, yticklabels=True, **heatmapkws2)
sp18_2.set(xlabel=None)
sp18_2.set(ylabel=None)
sp18_2.set_yticklabels(sp18_2.get_yticklabels(), rotation = 0, fontsize = 9)
sp18_2.yaxis.tick_right()
sp18_2.patch.set_linewidth(2)
sp18_2.patch.set_edgecolor('k')
labels19_1=  np.array([[12.0728,93.73441,90000,0.267805,90000,90000,15.4294,90000,90000,90000],
[111.8425,88.48765,90000,11.017875,90000,90000,71.98795,90000,90000,90000],
[212.6645,132.8535251,90000,0.05996195,90000,90000,194.293,90000,90000,90000],
[22.6686,49.0357409,90000,1.360463,90000,90000,54.922,90000,90000,90000]])
sp19_1=sns.heatmap(sp19_1, annot=labels19_1, ax=axes[18,0], xticklabels=False, yticklabels=False, **heatmapkws1)
sp19_1.set(xlabel=None)
sp19_1.set_ylabel('$\\it{Theobroma\ cacao}$', rotation =0, ha='right', va='center')
sp19_1.patch.set_linewidth(2)
sp19_1.patch.set_edgecolor('k')
for t in sp19_1.texts:
    if float(t.get_text())!=90000:
        t.set_text(t.get_text()) 
    else:
        t.set_text("")
labels19_2=  np.array([[7.764098635,90000,57.61430892],
[0.791180902,90000,6.533741761],
[0.624709461,90000,3240.271539],
[2.163157006,90000,40.37007989]])
sp19_2=sns.heatmap(sp19_2, annot=labels19_2, ax=axes[18,1], xticklabels=False, yticklabels=True, **heatmapkws2)
sp19_2.set(xlabel=None)
sp19_2.set(ylabel=None)
sp19_2.set_yticklabels(sp19_2.get_yticklabels(), rotation = 0, fontsize = 9)
sp19_2.yaxis.tick_right()
sp19_2.patch.set_linewidth(2)
sp19_2.patch.set_edgecolor('k')
for t in sp19_2.texts:
    if float(t.get_text())!=90000:
        t.set_text(t.get_text()) 
    else:
        t.set_text("")
labels20_1=  np.array([[471.29943,97.7724555,150.184,27.677225,90000,0.21961,31.90189,124.674735,90000,90000],
[728.18807,179.3020784,178.515,242.65727,90000,0.0386759,200.70079,91.949786,90000,90000],
[98.247072,53.31358,47.7779,0.109112,90000,0.01,20.49162,47.741965,90000,90000],
[236.40175,62.97446,22.8004,1.328338,90000,0.0628952,50.0575,113.001803,90000,90000],
[70.750657,15.0005525,9.343355,48.51065,90000,0.01,3.20285,7.4681,90000,90000]])
sp20_1=sns.heatmap(sp20_1, annot=labels20_1, ax=axes[19,0], xticklabels=False, yticklabels=False, **heatmapkws1)
sp20_1.set(xlabel=None)
sp20_1.set_ylabel('$\\it{Gossypium\ raimondii}$', rotation =0, ha='right', va='center')
sp20_1.patch.set_linewidth(2)
sp20_1.patch.set_edgecolor('k')
for t in sp20_1.texts:
    if float(t.get_text())!=90000:
        t.set_text(t.get_text()) 
    else:
        t.set_text("")
labels20_2=  np.array([[0.207452947,0.318659414,6.161765871],
[0.246230453,0.245149581,1.431935131],
[0.542648029,0.486303551,598.6564264],
[0.266387453,0.096447679,116.4849773],
[0.212019975,0.132060328,0.219971285]])
sp20_2=sns.heatmap(sp20_2, annot=labels20_2, ax=axes[19,1], xticklabels=False, yticklabels=True, **heatmapkws2)
sp20_2.set(xlabel=None)
sp20_2.set(ylabel=None)
sp20_2.set_yticklabels(sp20_2.get_yticklabels(), rotation = 0, fontsize = 9)
sp20_2.yaxis.tick_right()
sp20_2.patch.set_linewidth(2)
sp20_2.patch.set_edgecolor('k')
labels21_1=  np.array([[647.7803,131.0632471,213.8257,271.3614055,90000,0.313979,158.7232,232.0917,90000,90000],
[489.9819,74.47815,107.029,42.59482,90000,0.0189445,271.4573,114.823,90000,90000],
[124.71767,26.116367,71.4525,0.3031056,90000,0.01,37.07221,59.6819,90000,90000],
[91.53736,18.4505118,44.7054,2.971001,90000,0.0514728,44.3939,34.048,90000,90000],
[64.3521,12.886368,10.9509,2.3820139,90000,0.01,6.872474,13.75719,90000,90000],
[758.9096,12.503154,152.4148,4.522884,90000,0.0376942,186.2063,342.747,90000,90000]])
sp21_1=sns.heatmap(sp21_1, annot=labels21_1, ax=axes[20,0], xticklabels=False, yticklabels=False, **heatmapkws1)
sp21_1.set(xlabel=None)
sp21_1.set_ylabel('$\\it{Gossypium\ hirsutum}$', rotation =0, ha='right', va='center')
sp21_1.patch.set_linewidth(2)
sp21_1.patch.set_edgecolor('k')
for t in sp21_1.texts:
    if float(t.get_text())!=90000:
        t.set_text(t.get_text()) 
    else:
        t.set_text("")
labels21_2=  np.array([[0.202326695,0.330089847,1.539525565],
[0.152001839,0.218434599,9.371216383],
[0.209403904,0.572914006,319.2092459],
[0.201562638,0.488384196,22.862961],
[0.200247824,0.170171603,8.660597656],
[0.016475156,0.200833933,102.2211175]])
sp21_2=sns.heatmap(sp21_2, annot=labels21_2, ax=axes[20,1], xticklabels=False, yticklabels=True, **heatmapkws2)
sp21_2.set(xlabel=None)
sp21_2.set(ylabel=None)
sp21_2.set_yticklabels(sp21_2.get_yticklabels(), rotation = 0, fontsize = 9)
sp21_2.yaxis.tick_right()
sp21_2.patch.set_linewidth(2)
sp21_2.patch.set_edgecolor('k')
labels22_1=  np.array([[186.473525,72.71678525,71.969465,12.8514441,90000,0.01,8.42859,42.087705,90000,90000],
[508.42645,82.2322534,95.27086,105.11375,90000,0.01,32.35445,60.74034,90000,90000],
[349.20724,45.982865,251.0184,0.163662025,90000,0.01,54.53349,134.28425,90000,90000],
[282.0322,43.25406665,114.37992,56.714932,90000,0.2863745,45.05775,80.83385,90000,90000],
[127.858285,38.6908815,104.10104,3.180923155,90000,0.01,19.4355,32.79146,90000,90000],
[199.8227,8.685,69.307857,0.2049353,90000,0.01,2.16193,211.5684,90000,90000]])
sp22_1=sns.heatmap(sp22_1, annot=labels22_1, ax=axes[21,0], xticklabels=False, yticklabels=False, **heatmapkws1)
sp22_1.set(xlabel=None)
sp22_1.set_ylabel('$\\it{Gossypium\ barbadense}$', rotation =0, ha='right', va='center')
sp22_1.patch.set_linewidth(2)
sp22_1.patch.set_edgecolor('k')
for t in sp22_1.texts:
    if float(t.get_text())!=90000:
        t.set_text(t.get_text()) 
    else:
        t.set_text("")
labels22_2=  np.array([[0.389957691,0.385950043,270.0065094],
[0.16173874,0.18738376,0.90241662],
[0.131677868,0.718823585,1153.705265],
[0.153365703,0.405556245,2.191496775],
[0.302607543,0.814190805,16.27976693],
[0.04346353,0.346846765,1035.990486]])
sp22_2=sns.heatmap(sp22_2, annot=labels22_2, ax=axes[21,1], xticklabels=False, yticklabels=True, **heatmapkws2)
sp22_2.set(xlabel=None)
sp22_2.set(ylabel=None)
sp22_2.set_yticklabels(sp22_2.get_yticklabels(), rotation = 0, fontsize = 9)
sp22_2.yaxis.tick_right()
sp22_2.patch.set_linewidth(2)
sp22_2.patch.set_edgecolor('k')
labels23_1=  np.array([[26.7567,18.923759,90000,95.6988,90000,90000,0.17471,90000,90000,90000],
[57.511834,10.443759,90000,110.818,90000,90000,7.494966,90000,90000,90000],
[0.426504,15.6872184,90000,0.496939,90000,90000,0.01,90000,90000,90000]])
sp23_1=sns.heatmap(sp23_1, annot=labels23_1, ax=axes[22,0], xticklabels=False, yticklabels=False, **heatmapkws1)
sp23_1.set(xlabel=None)
sp23_1.set_ylabel('$\\it{Eruca\ vesicaria}$', rotation =0, ha='right', va='center')
sp23_1.patch.set_linewidth(2)
sp23_1.patch.set_edgecolor('k')
for t in sp23_1.texts:
    if float(t.get_text())!=90000:
        t.set_text(t.get_text()) 
    else:
        t.set_text("")
labels23_2=  np.array([[0.707253099,90000,0.001825624],
[0.181593218,90000,0.06763311],
[36.78094086,90000,0.020123194]])
sp23_2=sns.heatmap(sp23_2, annot=labels23_2, ax=axes[22,1], xticklabels=False, yticklabels=True, **heatmapkws2)
sp23_2.set(xlabel=None)
sp23_2.set(ylabel=None)
sp23_2.set_yticklabels(sp23_2.get_yticklabels(), rotation = 0, fontsize = 9)
sp23_2.yaxis.tick_right()
sp23_2.patch.set_linewidth(2)
sp23_2.patch.set_edgecolor('k')
for t in sp23_2.texts:
    if float(t.get_text())!=90000:
        t.set_text(t.get_text()) 
    else:
        t.set_text("")
labels24_1=  np.array([[0.0983734,14.289006,90000,1.44595,90000,90000,0.01,90000,90000,90000],
[0.436485,3.5811871,90000,3.37085,90000,90000,0.3541495,90000,90000,90000]])
sp24_1=sns.heatmap(sp24_1, annot=labels24_1, ax=axes[23,0], xticklabels=False, yticklabels=False, **heatmapkws1)
sp24_1.set(xlabel=None)
sp24_1.set_ylabel('$\\it{Isatis\ tinctoria}$', rotation =0, ha='right', va='center')
sp24_1.patch.set_linewidth(2)
sp24_1.patch.set_edgecolor('k')
for t in sp24_1.texts:
    if float(t.get_text())!=90000:
        t.set_text(t.get_text()) 
    else:
        t.set_text("")
labels24_2=  np.array([[145.2527411,90000,0.006915868],
[8.204605198,90000,0.105062373]])
sp24_2=sns.heatmap(sp24_2, annot=labels24_2, ax=axes[23,1], xticklabels=False, yticklabels=True, **heatmapkws2)
sp24_2.set(xlabel=None)
sp24_2.set(ylabel=None)
sp24_2.set_yticklabels(sp24_2.get_yticklabels(), rotation = 0, fontsize = 9)
sp24_2.yaxis.tick_right()
sp24_2.patch.set_linewidth(2)
sp24_2.patch.set_edgecolor('k')
for t in sp24_2.texts:
    if float(t.get_text())!=90000:
        t.set_text(t.get_text()) 
    else:
        t.set_text("")
labels25_1=  np.array([[0.01,12.069514,90000,8.91884,90000,90000,0.138588,90000,90000,90000],
[369.422,24.58312815,90000,278.577,90000,90000,0.3225625,90000,90000,90000],
[0.0138071,0.25051775,90000,0.2388655,90000,90000,0.01,90000,90000,90000],
[63.10445,14.1977205,90000,12.56925,90000,90000,109.531,90000,90000,90000]])
sp25_1=sns.heatmap(sp25_1, annot=labels25_1, ax=axes[24,0], xticklabels=False, yticklabels=False, **heatmapkws1)
sp25_1.set(xlabel=None)
sp25_1.set_ylabel('$\\it{Capsella\ grandiflora}$', rotation =0, ha='right', va='center')
sp25_1.patch.set_linewidth(2)
sp25_1.patch.set_edgecolor('k')
for t in sp25_1.texts:
    if float(t.get_text())!=90000:
        t.set_text(t.get_text()) 
    else:
        t.set_text("")
labels25_2=  np.array([[1206.9514,90000,0.015538792],
[0.066544841,90000,0.001157894],
[18.14412512,90000,0.041864564],
[0.224987628,90000,8.714203314]])
sp25_2=sns.heatmap(sp25_2, annot=labels25_2, ax=axes[24,1], xticklabels=False, yticklabels=True, **heatmapkws2)
sp25_2.set(xlabel=None)
sp25_2.set(ylabel=None)
sp25_2.set_yticklabels(sp25_2.get_yticklabels(), rotation = 0, fontsize = 9)
sp25_2.yaxis.tick_right()
sp25_2.patch.set_linewidth(2)
sp25_2.patch.set_edgecolor('k')
for t in sp25_2.texts:
    if float(t.get_text())!=90000:
        t.set_text(t.get_text()) 
    else:
        t.set_text("")
labels26_1=  np.array([[0.8231885,18.3430235,90000,10.50565,90000,90000,2.0691345,90000,90000,90000]])
sp26_1=sns.heatmap(sp26_1, annot=labels26_1, ax=axes[25,0], xticklabels=False, yticklabels=False, **heatmapkws1)
sp26_1.set(xlabel=None)
sp26_1.set_ylabel('$\\it{Eutrema\ salsugineum}$', rotation =0, ha='right', va='center')
sp26_1.patch.set_linewidth(2)
sp26_1.patch.set_edgecolor('k')
for t in sp26_1.texts:
    if float(t.get_text())!=90000:
        t.set_text(t.get_text()) 
    else:
        t.set_text("")
labels26_2=  np.array([[22.28289572,90000,0.196954448]])
sp26_2=sns.heatmap(sp26_2, annot=labels26_2, ax=axes[25,1], xticklabels=False, yticklabels=True, **heatmapkws2)
sp26_2.set(xlabel=None)
sp26_2.set(ylabel=None)
sp26_2.set_yticklabels(sp26_2.get_yticklabels(), rotation = 0, fontsize = 9)
sp26_2.yaxis.tick_right()
sp26_2.patch.set_linewidth(2)
sp26_2.patch.set_edgecolor('k')
for t in sp26_2.texts:
    if float(t.get_text())!=90000:
        t.set_text(t.get_text()) 
    else:
        t.set_text("")
labels27_1=  np.array([[65.9166,2.7484677,90000,169.028,90000,90000,0.6362,90000,90000,90000],
[14.2209298,38.0993351,90000,42.63515,90000,90000,0.3301958,90000,90000,90000],
[233.724848,84.91224041,90000,261.502,90000,90000,266.084465,90000,90000,90000]])
sp27_1=sns.heatmap(sp27_1, annot=labels27_1, ax=axes[26,0], xticklabels=False, yticklabels=False, **heatmapkws1)
sp27_1.set(xlabel=None)
sp27_1.set_ylabel('$\\it{Sinapis\ alba}$', rotation =0, ha='right', va='center')
sp27_1.patch.set_linewidth(2)
sp27_1.patch.set_edgecolor('k')
for t in sp27_1.texts:
    if float(t.get_text())!=90000:
        t.set_text(t.get_text()) 
    else:
        t.set_text("")
labels27_2=  np.array([[0.041696139,90000,0.003763873],
[2.679102959,90000,0.007744685],
[0.363300013,90000,1.017523633]])
sp27_2=sns.heatmap(sp27_2, annot=labels27_2, ax=axes[26,1], xticklabels=False, yticklabels=True, **heatmapkws2)
sp27_2.set(xlabel=None)
sp27_2.set(ylabel=None)
sp27_2.set_yticklabels(sp27_2.get_yticklabels(), rotation = 0, fontsize = 9)
sp27_2.yaxis.tick_right()
sp27_2.patch.set_linewidth(2)
sp27_2.patch.set_edgecolor('k')
for t in sp27_2.texts:
    if float(t.get_text())!=90000:
        t.set_text(t.get_text()) 
    else:
        t.set_text("")
labels28_1=  np.array([[55.138845,10.1343785,90000,37.8135,90000,90000,4.29505,90000,90000,90000],
[538.29819,55.4371752,90000,259.8,90000,90000,0.426784,90000,90000,90000],
[146.3924,2.725906,90000,79.60975,90000,90000,2.4951446,90000,90000,90000]])
sp28_1=sns.heatmap(sp28_1, annot=labels28_1, ax=axes[27,0], xticklabels=False, yticklabels=False, **heatmapkws1)
sp28_1.set(xlabel=None)
sp28_1.set_ylabel('$\\it{Arabidopsis\ lyrata}$', rotation =0, ha='right', va='center')
sp28_1.patch.set_linewidth(2)
sp28_1.patch.set_edgecolor('k')
for t in sp28_1.texts:
    if float(t.get_text())!=90000:
        t.set_text(t.get_text()) 
    else:
        t.set_text("")
labels28_2=  np.array([[0.183797439,90000,0.113585095],
[0.102985996,90000,0.001642741],
[0.018620543,90000,0.031342199]])
sp28_2=sns.heatmap(sp28_2, annot=labels28_2, ax=axes[27,1], xticklabels=False, yticklabels=True, **heatmapkws2)
sp28_2.set(xlabel=None)
sp28_2.set(ylabel=None)
sp28_2.set_yticklabels(sp28_2.get_yticklabels(), rotation = 0, fontsize = 9)
sp28_2.yaxis.tick_right()
sp28_2.patch.set_linewidth(2)
sp28_2.patch.set_edgecolor('k')
for t in sp28_2.texts:
    if float(t.get_text())!=90000:
        t.set_text(t.get_text()) 
    else:
        t.set_text("")
labels29_1=  np.array([[4.50061,0.245951,90000,6.77431,90000,90000,0.310677,90000,90000,90000],
[74.89375,20.84475,90000,155.5415,90000,90000,0.084128,90000,90000,90000],
[51.4126,3.08119,90000,21.8948,90000,90000,0.01,90000,90000,90000],
[18.07205,1.85115,90000,3.877115,90000,90000,7.189785,90000,90000,90000],
[4.50061,0.245951,90000,6.77431,90000,90000,0.310677,90000,90000,90000]])
sp29_1=sns.heatmap(sp29_1, annot=labels29_1, ax=axes[28,0], xticklabels=False, yticklabels=False, **heatmapkws1)
sp29_1.set(xlabel=None)
sp29_1.set_ylabel('$\\it{Arabidopsis\ thaliana}$', rotation =0, ha='right', va='center')
sp29_1.patch.set_linewidth(2)
sp29_1.patch.set_edgecolor('k')
for t in sp29_1.texts:
    if float(t.get_text())!=90000:
        t.set_text(t.get_text()) 
    else:
        t.set_text("")
labels29_2=  np.array([[0.05464837,90000,0.045861054],
[0.278324293,90000,0.000540872],
[0.05993064,90000,0.000456729],
[0.102431656,90000,1.854416235],
[0.05464837,90000,0.045861054]])
sp29_2=sns.heatmap(sp29_2, annot=labels29_2, ax=axes[28,1], xticklabels=False, yticklabels=True, **heatmapkws2)
sp29_2.set(xlabel=None)
sp29_2.set(ylabel=None)
sp29_2.set_yticklabels(sp29_2.get_yticklabels(), rotation = 0, fontsize = 9)
sp29_2.yaxis.tick_right()
sp29_2.patch.set_linewidth(2)
sp29_2.patch.set_edgecolor('k')
for t in sp29_2.texts:
    if float(t.get_text())!=90000:
        t.set_text(t.get_text()) 
    else:
        t.set_text("")
labels30_1=  np.array([[32.5872,11.13877205,90000,46.60255,90000,90000,90000,90000,90000,0.04540545],
[10.6272,10.279686,90000,13.1792,90000,90000,90000,90000,90000,0.2426],
[8.65518,45.86158,90000,0.0922218,90000,90000,90000,90000,90000,4.528545],
[39.4019,3.4159,90000,0.01,90000,90000,90000,90000,90000,0.01]])
sp30_1=sns.heatmap(sp30_1, annot=labels30_1, ax=axes[29,0], xticklabels=False, yticklabels=False, **heatmapkws1)
sp30_1.set(xlabel=None)
sp30_1.set_ylabel('$\\it{Carica\ papaya}$', rotation =0, ha='right', va='center')
sp30_1.patch.set_linewidth(2)
sp30_1.patch.set_edgecolor('k')
for t in sp30_1.texts:
    if float(t.get_text())!=90000:
        t.set_text(t.get_text()) 
    else:
        t.set_text("")
labels30_2=  np.array([[0.341814334,90000,0.000974313],
[0.967299571,90000,0.018407794],
[5.298743643,90000,49.10492964],
[0.086693789,90000,1]])
sp30_2=sns.heatmap(sp30_2, annot=labels30_2, ax=axes[29,1], xticklabels=False, yticklabels=True, **heatmapkws2)
sp30_2.set(xlabel=None)
sp30_2.set(ylabel=None)
sp30_2.set_yticklabels(sp30_2.get_yticklabels(), rotation = 0, fontsize = 9)
sp30_2.yaxis.tick_right()
sp30_2.patch.set_linewidth(2)
sp30_2.patch.set_edgecolor('k')
for t in sp30_2.texts:
    if float(t.get_text())!=90000:
        t.set_text(t.get_text()) 
    else:
        t.set_text("")
labels31_1=  np.array([[147.836,3.23613,90000,0.64579495,90000,90000,90000,928.768775,90000,90000]])
sp31_1=sns.heatmap(sp31_1, annot=labels31_1, ax=axes[30,0], xticklabels=False, yticklabels=False, **heatmapkws1)
sp31_1.set(xlabel=None)
sp31_1.set_ylabel('$\\it{Anacardium\ occidentale}$', rotation =0, ha='right', va='center')
sp31_1.patch.set_linewidth(2)
sp31_1.patch.set_edgecolor('k')
for t in sp31_1.texts:
    if float(t.get_text())!=90000:
        t.set_text(t.get_text()) 
    else:
        t.set_text("")
labels31_2=  np.array([[0.02189,90000,1438.179061]])
sp31_2=sns.heatmap(sp31_2, annot=labels31_2, ax=axes[30,1], xticklabels=False, yticklabels=True, **heatmapkws2)
sp31_2.set(xlabel=None)
sp31_2.set(ylabel=None)
sp31_2.set_yticklabels(sp31_2.get_yticklabels(), rotation = 0, fontsize = 9)
sp31_2.yaxis.tick_right()
sp31_2.patch.set_linewidth(2)
sp31_2.patch.set_edgecolor('k')
for t in sp31_2.texts:
    if float(t.get_text())!=90000:
        t.set_text(t.get_text()) 
    else:
        t.set_text("")
labels32_1 =  np.array([[56.03105,30.17494787,0.01,134.8035,90000,90000,90000,0.01,90000,90000],
[20.2508,69.00359165,0.01152605,50.46955,90000,90000,90000,4.03498,90000,90000],
[1.42086,13.613025,0.01,46.36958,90000,90000,90000,0.01,90000,90000],
[92.80495,66.77729802,0.01,241.6814943,90000,90000,90000,3.871535,90000,90000],
[84.0518,20.686488,0.01,92.182,90000,90000,90000,0.0593258,90000,90000],
[899.3365,439.6446655,0.0502064,243.0625719,90000,90000,90000,0.0067409,90000,90000]])
sp32_1=sns.heatmap(sp32_1, annot=labels32_1, ax=axes[31,0], xticklabels=False, yticklabels=False, **heatmapkws1)
sp32_1.set(xlabel=None)
sp32_1.set_ylabel('$\\it{Citrus\ trifoliata}$', rotation =0, ha='right', va='center')
sp32_1.patch.set_linewidth(2)
sp32_1.patch.set_edgecolor('k')
for t in sp32_1.texts:
    if float(t.get_text())!=90000:
        t.set_text(t.get_text()) 
    else:
        t.set_text("")
labels32_2=  np.array([[0.538539754,0.000178472,7.42E-05],
[3.407450158,0.000569165,0.079948801],
[9.580834847,0.007037991,0.000215659],
[0.719544572,0.000107753,0.016019162],
[0.246115943,0.000118974,0.000643572],
[0.488854467,0.0000558,0.0000277]])
sp32_2=sns.heatmap(sp32_2, annot=labels32_2, ax=axes[31,1], xticklabels=False, yticklabels=True, **heatmapkws2)
sp32_2.set(xlabel=None)
sp32_2.set(ylabel=None)
sp32_2.set_yticklabels(sp32_2.get_yticklabels(), rotation = 0, fontsize = 9)
sp32_2.yaxis.tick_right()
sp32_2.patch.set_linewidth(2)
sp32_2.patch.set_edgecolor('k')
labels33_1=  np.array([[2.177682,8.05147,0.01,85.8263,90000,90000,90000,0.01,90000,90000],
[0.705257,40.730129,0.01,193.736,90000,90000,90000,0.01,90000,90000],
[34.5180295,32.064825,0.01,220.05,90000,90000,90000,0.01,90000,90000]])
sp33_1=sns.heatmap(sp33_1, annot=labels33_1, ax=axes[32,0], xticklabels=False, yticklabels=False, **heatmapkws1)
sp33_1.set(xlabel=None)
sp33_1.set_ylabel('$\\it{Citrus\ sinensis}$', rotation =0, ha='right', va='center')
sp33_1.patch.set_linewidth(2)
sp33_1.patch.set_edgecolor('k')
for t in sp33_1.texts:
    if float(t.get_text())!=90000:
        t.set_text(t.get_text()) 
    else:
        t.set_text("")
labels33_2=  np.array([[3.697266176,0.004592039,0.000116514],
[57.7521797,0.014179228,0.0000516],
[0.928929764,0.000289704,0.0000454]])
sp33_2=sns.heatmap(sp33_2, annot=labels33_2, ax=axes[32,1], xticklabels=False, yticklabels=True, **heatmapkws2)
sp33_2.set(xlabel=None)
sp33_2.set(ylabel=None)
sp33_2.set_yticklabels(sp33_2.get_yticklabels(), rotation = 0, fontsize = 9)
sp33_2.yaxis.tick_right()
sp33_2.patch.set_linewidth(2)
sp33_2.patch.set_edgecolor('k')
labels34_1=  np.array([[468.566,177.338248,98.5465,25.0447379,90000,90000,110.462,90000,90000,90000],
[232.53905,134.3905407,4.7946695,59.83429305,90000,90000,67.71245,90000,90000,90000],
[60.50854,48.812918,9.2353605,0.1197975,90000,90000,35.10215,90000,90000,90000],
[557.477,186.0530034,2.6119228,0.3944377,90000,90000,170.695,90000,90000,90000],
[222.8518,134.3905407,1.4816674,7.806057,90000,90000,67.1763,90000,90000,90000]])
sp34_1=sns.heatmap(sp34_1, annot=labels34_1, ax=axes[33,0], xticklabels=False, yticklabels=False, **heatmapkws1)
sp34_1.set(xlabel=None)
sp34_1.set_ylabel('$\\it{Vitis\ vinifera}$', rotation =0, ha='right', va='center')
sp34_1.patch.set_linewidth(2)
sp34_1.patch.set_edgecolor('k')
for t in sp34_1.texts:
    if float(t.get_text())!=90000:
        t.set_text(t.get_text()) 
    else:
        t.set_text("")
labels34_2=  np.array([[0.378470158,0.210315089,4.410587184],
[0.577926764,0.020618771,1.131666249],
[0.806711218,0.152629042,293.012375],
[0.087560416,0.016566353,432.7552868],
[0.603048935,0.006648667,8.605663525]])
sp34_2=sns.heatmap(sp34_2, annot=labels34_2, ax=axes[33,1], xticklabels=False, yticklabels=True, **heatmapkws2)
sp34_2.set(xlabel=None)
sp34_2.set(ylabel=None)
sp34_2.set_yticklabels(sp34_2.get_yticklabels(), rotation = 0, fontsize = 9)
sp34_2.yaxis.tick_right()
sp34_2.patch.set_linewidth(2)
sp34_2.patch.set_edgecolor('k')
labels35_1=  np.array([[28.04254,0.03661115,8.754656,10.5915345,90000,90000,90000,1.212045,90000,90000],
[471.175,0.927505,66.6084,31.77083,90000,90000,90000,47.8659,90000,90000],
[823.494,1.19574,329.495,12.63388,90000,90000,90000,1.53237,90000,90000],
[372.4613,1.65234,48.9721,5.948756,90000,90000,90000,1.61369,90000,90000]])
sp35_1=sns.heatmap(sp35_1, annot=labels35_1, ax=axes[34,0], xticklabels=False, yticklabels=False, **heatmapkws1)
sp35_1.set(xlabel=None)
sp35_1.set_ylabel('$\\it{Kalanchoe\ laxiflora}$', rotation =0, ha='right', va='center')
sp35_1.patch.set_linewidth(2)
sp35_1.patch.set_edgecolor('k')
for t in sp35_1.texts:
    if float(t.get_text())!=90000:
        t.set_text(t.get_text()) 
    else:
        t.set_text("")
labels35_2=  np.array([[0.001305558,0.312191977,0.11443526],
[0.001968494,0.141366584,1.506598978],
[0.001452032,0.400118277,0.12129053],
[0.004436273,0.131482385,0.271265118]])
sp35_2=sns.heatmap(sp35_2, annot=labels35_2, ax=axes[34,1], xticklabels=False, yticklabels=True, **heatmapkws2)
sp35_2.set(xlabel=None)
sp35_2.set(ylabel=None)
sp35_2.set_yticklabels(sp35_2.get_yticklabels(), rotation = 0, fontsize = 9)
sp35_2.yaxis.tick_right()
sp35_2.patch.set_linewidth(2)
sp35_2.patch.set_edgecolor('k')
labels36_1=  np.array([[81.736341,16.176699,90000,20.0873,90000,90000,1.808913,90000,90000,90000],
                      [224.430338,18.07722,90000,37.47141,90000,90000,226.4927,90000,90000,90000],
                      [263.353058,6.4867447,90000,39.633935,90000,90000,299.8325185,90000,90000,90000],
                      [4.19184885,12.59345425,90000,1.69929,90000,90000,0.030696,90000,90000,90000]])
sp36_1=sns.heatmap(sp36_1, annot=labels36_1, ax=axes[35,0], xticklabels=False, yticklabels=False, **heatmapkws1)
sp36_1.set(xlabel=None)
sp36_1.set_ylabel('$\\it{Dianthus\ caryophyllus}$', rotation =0, ha='right', va='center')
sp36_1.patch.set_linewidth(2)
sp36_1.patch.set_edgecolor('k')
for t in sp36_1.texts:
    if float(t.get_text())!=90000:
        t.set_text(t.get_text()) 
    else:
        t.set_text("")
labels36_2=  np.array([[0.197913177933913,90000,0.090052570529638],
                       [0.080547131734035,90000,6.0444135942576],
                       [0.024631362738913,90000,7.56504542130374],
                       [3.004272029,90000,0.018064015]])
sp36_2=sns.heatmap(sp36_2, annot=labels36_2, ax=axes[35,1], xticklabels=False, yticklabels=True, **heatmapkws2)
sp36_2.set(xlabel=None)
sp36_2.set(ylabel=None)
sp36_2.set_yticklabels(sp36_2.get_yticklabels(), rotation = 0, fontsize = 9)
sp36_2.yaxis.tick_right()
sp36_2.patch.set_linewidth(2)
sp36_2.patch.set_edgecolor('k')
for t in sp36_2.texts:
    if float(t.get_text())!=90000:
        t.set_text(t.get_text()) 
    else:
        t.set_text("")
labels37_1=  np.array([[15.21934325,1.76283,90000,4.72757,90000,90000,0.01,90000,90000,90000],
[26.94864,8.86334,90000,0.05127865,90000,90000,0.01,90000,90000,90000]])
sp37_1=sns.heatmap(sp37_1, annot=labels37_1, ax=axes[36,0], xticklabels=False, yticklabels=False, **heatmapkws1)
sp37_1.set(xlabel=None)
sp37_1.set_ylabel('$\\it{Beta\ vulgaris}$', rotation =0, ha='right', va='center')
sp37_1.patch.set_linewidth(2)
sp37_1.patch.set_edgecolor('k')
for t in sp37_1.texts:
    if float(t.get_text())!=90000:
        t.set_text(t.get_text()) 
    else:
        t.set_text("")
labels37_2=  np.array([[0.115828257,90000,0.002115252],
[0.328897488,90000,0.195012934]])
sp37_2=sns.heatmap(sp37_2, annot=labels37_2, ax=axes[36,1], xticklabels=False, yticklabels=True, **heatmapkws2)
sp37_2.set(xlabel=None)
sp37_2.set(ylabel=None)
sp37_2.set_yticklabels(sp37_2.get_yticklabels(), rotation = 0, fontsize = 9)
sp37_2.yaxis.tick_right()
sp37_2.patch.set_linewidth(2)
sp37_2.patch.set_edgecolor('k')
for t in sp37_2.texts:
    if float(t.get_text())!=90000:
        t.set_text(t.get_text()) 
    else:
        t.set_text("")
labels38_1=  np.array([[45.83726,45.633452,90000,47.02709,90000,90000,1.58158,90000,90000,90000],
[292.7279,142.3930819,90000,106.2388,90000,90000,361.429,90000,90000,90000],
[1.96225725,3.75595725,90000,0.536495,90000,90000,0.01177905,90000,90000,90000],
[10.416285,13.559449,90000,16.84775,90000,90000,1.033586,90000,90000,90000]])
sp38_1=sns.heatmap(sp38_1, annot=labels38_1, ax=axes[37,0], xticklabels=False, yticklabels=False, **heatmapkws1)
sp38_1.set(xlabel=None)
sp38_1.set_ylabel('$\\it{Lactuca\ sativa}$', rotation =0, ha='right', va='center')
sp38_1.patch.set_linewidth(2)
sp38_1.patch.set_edgecolor('k')
for t in sp38_1.texts:
    if float(t.get_text())!=90000:
        t.set_text(t.get_text()) 
    else:
        t.set_text("")
labels38_2=  np.array([[0.995553661,90000,0.033631254],
[0.486434952,90000,3.402043321],
[1.91410033,90000,0.021955563],
[1.3017548,90000,0.061348607]])
sp38_2=sns.heatmap(sp38_2, annot=labels38_2, ax=axes[37,1], xticklabels=False, yticklabels=True, **heatmapkws2)
sp38_2.set(xlabel=None)
sp38_2.set(ylabel=None)
sp38_2.set_yticklabels(sp38_2.get_yticklabels(), rotation = 0, fontsize = 9)
sp38_2.yaxis.tick_right()
sp38_2.patch.set_linewidth(2)
sp38_2.patch.set_edgecolor('k')
for t in sp38_2.texts:
    if float(t.get_text())!=90000:
        t.set_text(t.get_text()) 
    else:
        t.set_text("")
labels39_1=  np.array([[1.2825025,0.9923595,90000,8.66905,90000,90000,0.01,90000,90000,90000],
[2.1824955,34.408758,90000,19.32085,90000,90000,0.01,90000,90000,90000],
[8.0652825,17.934608,90000,1.10398,90000,90000,75.93463,90000,90000,90000]])
sp39_1=sns.heatmap(sp39_1, annot=labels39_1, ax=axes[38,0], xticklabels=False, yticklabels=False, **heatmapkws1)
sp39_1.set(xlabel=None)
sp39_1.set_ylabel('$\\it{Helianthus\ annuus}$', rotation =0, ha='right', va='center')
sp39_1.patch.set_linewidth(2)
sp39_1.patch.set_edgecolor('k')
for t in sp39_1.texts:
    if float(t.get_text())!=90000:
        t.set_text(t.get_text()) 
    else:
        t.set_text("")
labels39_2=  np.array([[0.773768082,90000,0.001153529],
[15.76578646,90000,0.000517576],
[2.223680076,90000,68.78261382]])
sp39_2=sns.heatmap(sp39_2, annot=labels39_2, ax=axes[38,1], xticklabels=False, yticklabels=True, **heatmapkws2)
sp39_2.set(xlabel=None)
sp39_2.set(ylabel=None)
sp39_2.set_yticklabels(sp39_2.get_yticklabels(), rotation = 0, fontsize = 9)
sp39_2.yaxis.tick_right()
sp39_2.patch.set_linewidth(2)
sp39_2.patch.set_edgecolor('k')
for t in sp39_2.texts:
    if float(t.get_text())!=90000:
        t.set_text(t.get_text()) 
    else:
        t.set_text("")
labels40_1=  np.array([[14.55995,31.26833,90000,26.3293,90000,90000,90000,2.908662,90000,90000],
[471.0275,303.364188,90000,1033.3285,90000,90000,90000,50.848215,90000,90000],
[14.1242,22.296024,90000,1.8309,90000,90000,90000,17.19085,90000,90000]])
sp40_1=sns.heatmap(sp40_1, annot=labels40_1, ax=axes[39,0], xticklabels=False, yticklabels=False, **heatmapkws1)
sp40_1.set(xlabel=None)
sp40_1.set_ylabel('$\\it{Solanum\ tuberosum}$', rotation =0, ha='right', va='center')
sp40_1.patch.set_linewidth(2)
sp40_1.patch.set_edgecolor('k')
for t in sp40_1.texts:
    if float(t.get_text())!=90000:
        t.set_text(t.get_text()) 
    else:
        t.set_text("")
labels40_2=  np.array([[2.147557512,90000,0.110472439],
[0.644047721,90000,0.04920818],
[1.578568981,90000,9.389289421]])
sp40_2=sns.heatmap(sp40_2, annot=labels40_2, ax=axes[39,1], xticklabels=False, yticklabels=True, **heatmapkws2)
sp40_2.set(xlabel=None)
sp40_2.set(ylabel=None)
sp40_2.set_yticklabels(sp40_2.get_yticklabels(), rotation = 0, fontsize = 9)
sp40_2.yaxis.tick_right()
sp40_2.patch.set_linewidth(2)
sp40_2.patch.set_edgecolor('k')
for t in sp40_2.texts:
    if float(t.get_text())!=90000:
        t.set_text(t.get_text()) 
    else:
        t.set_text("")
labels41_1=  np.array([[3.26271,9.1894765,0.2756865,12.88905,90000,90000,90000,0.4339,90000,90000],
[154.4675,11.8492,0.410126,127.7016,90000,90000,90000,0.353537,90000,90000],
[1.359085,15.4101225,0.01,1.48288,90000,90000,90000,0.04947685,90000,90000],
[10.73544,52.8789065,0.193818,49.7925,90000,90000,90000,0.65723035,90000,90000],
[48.07315,5.066896,0.0175803,112.33765,90000,90000,90000,0.192665,90000,90000]])
sp41_1=sns.heatmap(sp41_1, annot=labels41_1, ax=axes[40,0], xticklabels=False, yticklabels=False, **heatmapkws1)
sp41_1.set(xlabel=None)
sp41_1.set_ylabel('$\\it{Solanum\ lycopersicum}$', rotation =0, ha='right', va='center')
sp41_1.patch.set_linewidth(2)
sp41_1.patch.set_edgecolor('k')
for t in sp41_1.texts:
    if float(t.get_text())!=90000:
        t.set_text(t.get_text()) 
    else:
        t.set_text("")
labels41_2=  np.array([[2.816516485,0.08449617,0.033664234],
[0.076709988,0.002655096,0.002768462],
[11.33860097,0.007357892,0.033365377],
[4.925639424,0.018054034,0.013199384],
[0.105399709,0.000365699,0.001715053]])
sp41_2=sns.heatmap(sp41_2, annot=labels41_2, ax=axes[40,1], xticklabels=False, yticklabels=True, **heatmapkws2)
sp41_2.set(xlabel=None)
sp41_2.set(ylabel=None)
sp41_2.set_yticklabels(sp41_2.get_yticklabels(), rotation = 0, fontsize = 9)
sp41_2.yaxis.tick_right()
sp41_2.patch.set_linewidth(2)
sp41_2.patch.set_edgecolor('k')
labels42_1=  np.array([[0.0637741,577.0334284,90000,19.3355,90000,90000,0.11798,90000,90000,90000],
[0.0486327,213.6589459,90000,0.01,90000,90000,0.0464564,90000,90000,90000],
[0.019547,116.689786,90000,1.24283,90000,90000,0.0737128,90000,90000,90000]])
sp42_1=sns.heatmap(sp42_1, annot=labels42_1, ax=axes[41,0], xticklabels=False, yticklabels=False, **heatmapkws1)
sp42_1.set(xlabel=None)
sp42_1.set_ylabel('$\\it{Sesamum\ indicum}$', rotation =0, ha='right', va='center')
sp42_1.patch.set_linewidth(2)
sp42_1.patch.set_edgecolor('k')
for t in sp42_1.texts:
    if float(t.get_text())!=90000:
        t.set_text(t.get_text()) 
    else:
        t.set_text("")
labels42_2=  np.array([[9048.084229,90000,0.00610173],
[4393.318609,90000,4.64564],
[5969.703075,90000,0.059310445]])
sp42_2=sns.heatmap(sp42_2, annot=labels42_2, ax=axes[41,1], xticklabels=False, yticklabels=True, **heatmapkws2)
sp42_2.set(xlabel=None)
sp42_2.set(ylabel=None)
sp42_2.set_yticklabels(sp42_2.get_yticklabels(), rotation = 0, fontsize = 9)
sp42_2.yaxis.tick_right()
sp42_2.patch.set_linewidth(2)
sp42_2.patch.set_edgecolor('k')
for t in sp42_2.texts:
    if float(t.get_text())!=90000:
        t.set_text(t.get_text()) 
    else:
        t.set_text("")
labels43_1=  np.array([[26.5252812,227.403968,90000,41.2094,90000,90000,144.294,90000,90000,90000],
[31.476462,73.8105827,90000,128.585,90000,90000,0.36654,90000,90000,90000],
[10.20165,28.66948615,90000,8.64008,90000,90000,37.4065,90000,90000,90000]])
sp43_1=sns.heatmap(sp43_1, annot=labels43_1, ax=axes[42,0], xticklabels=Genes1, yticklabels=False, **heatmapkws1)
sp43_1.set(xlabel=None)
sp43_1.set_ylabel('$\\it{Erythranthe\ guttata}$', rotation =0, ha='right', va='center')
sp43_1.patch.set_linewidth(2)
sp43_1.patch.set_edgecolor('k')
for t in sp43_1.texts:
    if float(t.get_text())!=90000:
        t.set_text(t.get_text()) 
    else:
        t.set_text("")
labels43_2=  np.array([[8.573103006,90000,3.501482671],
[2.344945334,90000,0.002850566],
[2.810279332,90000,4.329415931]])
sp43_2=sns.heatmap(sp43_2, annot=labels43_2, ax=axes[42,1], xticklabels=Genes2, yticklabels=True, **heatmapkws2)
sp43_2.set(xlabel=None)
sp43_2.set(ylabel=None)
sp43_2.set_yticklabels(sp43_2.get_yticklabels(), rotation = 0, fontsize = 9)
sp43_2.yaxis.tick_right()
sp43_2.patch.set_linewidth(2)
sp43_2.patch.set_edgecolor('k')
for t in sp43_2.texts:
    if float(t.get_text())!=90000:
        t.set_text(t.get_text()) 
    else:
        t.set_text("")

sm1 = matplotlib.cm.ScalarMappable(cmap=my_cmap1,norm=norm1)
sm1.set_array([])
clb1_1=fig.colorbar(sm1, cax=fig.add_axes([0.06,0.93,0.6, 0.01]),  orientation='horizontal',ticks=[-6,2,10])
clb1_1.ax.set_title('Gene Expression (TPMs)')
clb1_1.ax.set_xticklabels(['0','4','1033'])
clb1_2=fig.colorbar(sm1, cax=fig.add_axes([0.06,0.07,0.6, 0.01]),  orientation='horizontal', label='Gene Expression (TPMs)',ticks=[-6,2,10])
clb1_2.ax.set_xticklabels(['0','4','1033'])
sm2 = matplotlib.cm.ScalarMappable(cmap=my_cmap2,norm=norm2)
sm2.set_array([])
clb2_2=fig.colorbar(sm2, cax=fig.add_axes([0.68,0.07,0.19, 0.01]),  orientation='horizontal',label='Ratio of Gene Expression (TPMs)',ticks=[-10,90000,14.9])
clb2_2.ax.set_xticklabels(['0','1','32000'])
clb2_1=fig.colorbar(sm2, cax=fig.add_axes([0.68,0.93,0.19, 0.01]),  orientation='horizontal',ticks=[-10,90000,14.9])
clb2_1.ax.set_xticklabels(['0','1','32000'])
clb2_1.ax.set_title('Ratio of Gene Expression (TPMs)')

plt.savefig("AdditionalFile7.pdf", dpi=800, bbox_inches='tight')
