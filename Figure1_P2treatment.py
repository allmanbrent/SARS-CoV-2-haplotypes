#!/usr/bin/env python


import pandas as pd
#import PyVCF
import sys
import numpy as np
import os
import matplotlib.pyplot as plt
from datetime import date as dt
import seaborn as sns 
import itertools
from matplotlib.lines import Line2D
import matplotlib.artist as artists
from functools import reduce
from matplotlib.ticker import MaxNLocator
import matplotlib.patches as patches
from brokenaxes import brokenaxes
from matplotlib.gridspec import GridSpec


fig = plt.figure(figsize=(3, 5), constrained_layout=True)


directory = '../Data/20220404_'


numcols=['#bdbde6', '#6abdcd', '#de3131', '#fc7f03', '#5a4a73'] #butterfree pallette
#dpiNum=[[30], [28,39,53,75], [51,56,79,91,101], [77, 100, 109], [192]]
dpiNum=[[28,39,53,75], [51,56,79,91,101], [77, 100, 109]]
shapes = ['o', 's', 'v', '^', 'P']

alph = 0.6

################
## PATIENT 2 ###
################
samples=['705E', '798T', '962B', '1064Z']
date='20220210'
dpi2=['28','39','53','75']
# First positive sample taken 20210104
# Symptom onset 20210104
# mAb d4
# Sample  Patient Date
# GA-EHC-705E IC_patient2 20210201
# GA-EHC-798T IC_patient2 20210212
# GA-EHC-962B IC_patient2 20210226
# GA-EHC-1064Z    IC_patient2 20210320

#POSITION and TIME iSNVs. HIGH FREQUENCY
dpiNum=[[28,39,53,75], [51,56,79,91,101], [77, 100, 109]]
marker = itertools.cycle(('o', 's', 'v', '^', 'P', 'X'))
markerpos = itertools.cycle(('o', 's', 'v', '^'))
pal = itertools.cycle(sns.color_palette("colorblind", n_colors=16))
lins = itertools.cycle(('-', '--', '-.', ':'))

alph = 0.6

P2data = pd.read_csv('P2-HIGH-FREQ-TABLE-10_2022-04-28.tsv', sep = '\t')  
P2data = pd.read_csv('P2-HIGH-FREQ-TABLE-10_2022-05-06-WITH-WUHAN.txt', sep = '\t')  

samples=['705E', '798T', '962B', '1064Z']
days = dpiNum[0]
start_date = '2020-12-28'#'2021-01-25'
end_date = '2021-03-20'
loci = [2419, 7086, 21557, 22363, 22999, 23027, 23029, 23562]
P2705E = ['T', 'C', 'T', 'T', 'G', 'G', 'T', 'C']
var = ['C', 'T', 'G', 'C', 'A', 'A', 'C', 'T']
Wuhan = ['G', 'A', 'T']
seq_file = 'GA-P2-period-aligned-705E.fasta'#'nextclade_P2_timeframe.fasta'
Wuhan_idx = True
 
counter = 0


#timePlot.legend(title = "iSNV position")
#legs = [Line2D(range(1), range(1), color="white", marker=i, markersize=10, markerfacecolor="w", markeredgecolor='k') for i in shapes]
#[axis.legend(leg, 'd'+dpiString, title=patientNum) for axis,dpiString,patientNum,leg in zip([axp2, axp3, axp6], [dpi2, dpi3, dpi6], ['P2','P3','P4'], [legs[0:3], legs[0,4], legs[0,2]])]

# current_handles, current_labels = posPlot.get_legend_handles_labels()
# print(posPlot.get_legend_handles_labels())
# posPlot.legend(handles = current_handles, labels =current_labels,  framealpha = 0.5, loc = 'upper left', bbox_to_anchor=(0.03, 1), fontsize=10)
# leg = posPlot.get_legend()
# leg = ax.get_legend()
# hl_dict = {handle.get_label(): handle for handle in leg.legendHandles}
import matplotlib.lines as mlines

cols = ['#ffff00', '#00d48f', '#008eac', '#3f475d']


# fig= plt.subplots(nrows=2,ncols=1, sharex='col',
#                     gridspec_kw={
#                            'width_ratios': [1],
#                            'height_ratios': [1,1],
#                        'wspace': 0.16,
#                        'hspace': 0.1})

fig = plt.figure(constrained_layout=True)
#axs = [P2]
alph = 0.3

#fig = plt.figure(constrained_layout=False)
# # create left/right subfigs
# (Ct_fig, lib_pos) = fig.subfigures(2, 1, hspace=0.1, wspace=0.2, width_ratios=[1], height_ratios=[1.2,1])
# gs = lib_pos.add_gridspec(1, 4)
# # HEIGHT by WIDTH
# axt1 = lib_pos.add_subplot(gs[0, 0])
# axt2 = lib_pos.add_subplot(gs[0, 1], sharey=axt1)
# axt3 = lib_pos.add_subplot(gs[0, 2], sharey=axt1)
# axt4 = lib_pos.add_subplot(gs[0, 3], sharey=axt1)
samples=['705E', '798T', '962B', '1064Z']
date='20220210'
directory = '../Data/20220404_'
dpi2=['28','39','53','75']

new_cols = ['yellow', 'red', 'blue', 'purple', 'white']
new_cols = ['#FFFF00', '#FF0000', '#00FFFF', '#BA55D3', '#FFFFFF']

directory = '20220121_'
# gs2 = Ct_fig.add_gridspec(1,1)
# Ct_subplot = Ct_fig.add_subplot(gs2[0, 0])
gs = plt.GridSpec(1, 1, figure = fig)
Ct_subplot = fig.add_subplot(gs[0,:])
# axs = [Ct_fig.gca()]
axs = [Ct_subplot]

numcols=['#bdbde6', '#6abdcd', '#de3131', '#fc7f03', '#5a4a73'] #butterfree pallette
#dpiNum=[[30], [28,39,53,75], [51,56,79,91,101], [77, 100, 109], [192]]
dpiNum=[[28,39,53,75]]

#CtP1 = [26] #25.6
CtP2 = [18, 15, 16, 16]
CtP3 = [23, 19, 23, 21, 23]
CtP6 = [23, 23, 24]
CtAll = [CtP2, CtP3, CtP6]
remdesP2 = [[28,32], [40,45]]
remdesP3 = [52,58]
remdesP6 = [[0,5], [101,110]]
labs = ['P2', 'P3', 'P4']
#remdesivir_trtmt = [[30,35], [[28,32], [40,45]], [52,58], [[0,5], [101,110]], [196,200]]

i = 0
for remdesivir_trtmt in [remdesP2]:
	#print(remdesivir_trtmt)
	axs[i].plot(dpiNum[i], CtAll[i], label=labs[i], marker='o', alpha=alph, color='black')
	#[axs[i].text(dpi-3, Ct-3, Ct) for dpi, Ct in zip(dpiNum[i], CtAll[i])]
	if isinstance(remdesivir_trtmt[0],list):
		#print(remdesivir_trtmt)
		#print(1)	
		for remwindow in remdesivir_trtmt:
			#rect = patches.Rectangle((remwindow[0], .2), (remwindow[1]-remwindow[0]), 0.1, linewidth=1, edgecolor='none', facecolor=cols[i],alpha=0.7)
			#rectGen = patches.Rectangle((remwindow[0], -.001), (remwindow[1]-remwindow[0]), 0.001, linewidth=1, edgecolor='none', facecolor=cols[i],alpha=0.7)
			#axs[i].add_patch(rectGen)
			rectGen = patches.Rectangle((remwindow[0], -.001), (remwindow[1]-remwindow[0]), 30, linewidth=1, edgecolor='none', facecolor="grey",alpha=0.3)
			axs[i].add_patch(rectGen)
	else:
		#print(2)
		rectGen = patches.Rectangle((remdesivir_trtmt[0], -.001), (remdesivir_trtmt[1]-remdesivir_trtmt[0]), 30, linewidth=1, edgecolor='none', facecolor="grey",alpha=0.3)
		axs[i].add_patch(rectGen)
	
	#mAb Treatment
	if labs[i] == 'P2':
		#axs[i].scatter(x=[4],y=[3],marker="1",s=100,color='green')
		axs[i].axvline(x=[4], color='green', linestyle='--', lw=1, alpha = 0.8)
		axs[i].scatter(x=[86],y=[10],marker="x",s=80,color='k')
	elif labs[i] == 'P3':
		P3.scatter(x=[8],y=[3],marker="1",s=100)
	elif labs[i] == 'P4':
		P6.scatter(x=[0,104],y=[3,3],marker="1",s=100)
	
	axs[i].set_ylim(0, 30)
	i+=1

Ct_subplot.set_ylim(9, 21)
plt.yticks(ticks=[10, 15, 20], labels = ['10', '15', '20'], fontsize=10)
#plt.yticks(ticks=[5, 10, 15, 25, 30], labels = ['5', '', '15', '25', '30'], fontsize=10)
# Ct_fig.gca().set_ylim(0, 20)
# days = [4,28,39,53,75]
# Ct_fig.gca().set_xlim(min(days)-5, max(days)+7)
# Ct_fig.gca().set_xticks(days) 
# Ct_fig.gca().set_xticklabels(days, fontsize=10)

# #P2.text(0, 25, 'Patient 2', fontsize=12)


# axt1.set_title('B)', loc='left', fontsize=12)
# Ct_fig.gca().set_title('A)', loc='left', fontsize=12)

# Ct_fig.supylabel('Ct')
# Ct_fig.supxlabel('day', y = -0.06, fontsize=12)

days = [0,28,39,53,75,86]
Ct_subplot.set_xlim(min(days)-5, max(days)+7)
Ct_subplot.set_xticks(days) 
Ct_subplot.set_xticklabels(days, fontsize=10)

#P2.text(0, 25, 'Patient 2', fontsize=12)


#axt1.set_title('B)', loc='left', fontsize=12)
#Ct_subplot.set_title('A)', loc='left', fontsize=12)

Ct_subplot.set_ylabel('Ct', fontsize=12)
Ct_subplot.set_xlabel('Days after first positive test', y = -0.06, fontsize=12)
fig.align_labels()

fig.set_size_inches(3, 2.5)

# plt.tight_layout()
del date
from datetime import date

fig.savefig('Figure-1-Treatment-'+date.today().isoformat()+'.pdf',bbox_inches='tight')
fig.savefig('Figure-1-Treatment-'+date.today().isoformat()+'.png',bbox_inches='tight')














