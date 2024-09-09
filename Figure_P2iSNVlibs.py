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


#fig = plt.figure(constrained_layout=False)
# create left/right subfigs
#(posPlot,timePlot) = fig.subplots(1, 2, hspace=0.05, width_ratios=[1, 1])

fig = plt.figure(figsize=(3, 5), constrained_layout=True)
# posPlot = fig.add_subplot(grid[0, 0])
# timePlot = fig.add_subplot(grid[0, 1])
posPlot1, timePlot1 = plt.GridSpec(1, 2, hspace=0.25, wspace=0.2, width_ratios=[1,1.4])

posPlot = brokenaxes(xlims=((0, 10000), (20000, 26000)), subplot_spec=posPlot1)
timePlot = brokenaxes(xlims=((-1, 5), (26, 75)), subplot_spec=timePlot1)

#posPlot = freqPlots.add_subplot(gs2[0, 0])
cols = ['#ffff00', '#00d48f', '#008eac', '#3f475d']
#directory = '20220121_'
#directory = '20220210_'
directory = '../Data/20220404_'


numcols=['#bdbde6', '#6abdcd', '#de3131', '#fc7f03', '#5a4a73'] #butterfree pallette
#dpiNum=[[30], [28,39,53,75], [51,56,79,91,101], [77, 100, 109], [192]]
dpiNum=[[28,39,53,75], [51,56,79,91,101], [77, 100, 109]]
shapes = ['o', 's', 'v', '^', 'P']

alph = 0.6

def compareLibs(samples = [], date=[], dpiString=[], dpiNum=[], complot= [], colors=[], patient = 'P2',alph = 0.5):
	L1List=pd.read_csv(directory+samples+'/'+samples+'_L1_lofreq-output.txt', sep='\t')
	L2List=pd.read_csv(directory+samples+'/'+samples+'_L2_lofreq-output.txt', sep='\t')
	combined_libs = pd.merge(L1List, L2List, on = 'POSITION', how = 'outer', 
	suffixes = ['_L1','_L2'])[['POSITION', 'ALLELE-FREQUENCY_L1', 'ALLELE-FREQUENCY_L2']].fillna(0)

	condition = []
	count = 0
	# new_cols = ['dusty', 'red', 'blue', 'grey']
	for i in combined_libs.index:
		L1_val = combined_libs.iloc[i]['ALLELE-FREQUENCY_L1']
		L2_val = combined_libs.iloc[i]['ALLELE-FREQUENCY_L2']
		if (L1_val > 0.02) and (L2_val > 0.02):
			count+= 1
		# if (L1_val < 0.02) and (L2_val <0.02):
		# 	condition.append('0')
		# elif (L1_val > 0.1) and (L2_val > 0.1):
		# 	# both are above 10%
		# 	condition.append('1')
		# elif ((L1_val < 0.1) and (L2_val < 0.1)) and ((L1_val > 0.02) and (L2_val > 0.02)):
		# 	# both are above 2% but are both below 10%
		# 	condition.append('2')
		# elif ((L1_val > 0.1) ^ (L2_val > 0.1)) and ((L1_val > 0.02) and (L2_val > 0.02)):
		# 	#one is above 10%, and both are above 2%
		# 	condition.append('3')
		# else:
		# 	#one is above 10%, and the other is below 2%
		# 	#print(combined_libs.iloc[i])
		# 	condition.append('4')
		if combined_libs.iloc[i]['POSITION'] not in [2419, 7068,20467,21557,21626,22363,22999,23027,23029,23050,23562,24782]:
			condition.append('4')
		else:
			condition.append('1')

	combined_libs['condition'] = condition
	combined_libs = combined_libs.groupby('condition')
	
	for name, group in combined_libs:
		if name == '0':
			c = colors[0]
		elif name == '1':
			c = colors[1]
		elif name == '2':
			c = colors[2]
		elif name == '3':
			c = colors[3]
		elif name == '4':
			c = colors[4]
		complot.scatter(group['ALLELE-FREQUENCY_L1'], group['ALLELE-FREQUENCY_L2'], 
			color=c, edgecolor='#333333', alpha=alph, s = 8, linewidth = 0.2)
		
	complot.set_ylim(-0.05, 1.05)
	complot.set_xlim(-0.05, 1.05)
	#complot.set_xlabel('L1 frequency')
	#complot.set_ylabel('L2 frequency')
	complot.axhline(y=0.02, color='red', linestyle='--', lw=1, alpha = 0.7)
	complot.axvline(x=0.02, color='red', linestyle='--', lw=1, alpha = 0.7)
	complot.axhline(y=0.1, color='red', linestyle=':', lw=1, alpha = 0.7)
	complot.axvline(x=0.1, color='red', linestyle=':', lw=1, alpha = 0.7)

	#complot.set_title(patient + ' d_'+ dpiString + ' sample ' + samples, loc = 'left')
	#complot.set_title(patient + ' d'+ dpiString, loc = 'left')
	complot.text(0.2, 0.9, 'd'+dpiString, size=10, color='k')
	xpoints = ypoints = ax.get_xlim()
	complot.plot(xpoints, ypoints, linestyle='--', color='k', lw=1, scalex=False, scaley=False, alpha=0.4)
	print(dpiString)
	
	#print(combined_libs)
	print(count)
	return(complot)

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
for i_pos in P2data['POSITION']:
	if i_pos != 7613:
		j = 0
		c = next(pal)
		m1 = next(marker)
		if counter ==0:
			m2 = next(markerpos)
		l = next(lins)
		for d in days:
			d = str(d)
			if counter < 1:
				if i_pos < 10000:
					if i_pos in [23027, 22363]:
						print("d"+d+" "+str(i_pos) + " " + str(float(1-P2data[P2data['POSITION'] == i_pos]['AF_d'+d])))
						posPlot.axs[0].plot(i_pos,float(1-P2data[P2data['POSITION'] == i_pos]['AF_d'+d]),
							color=c, markeredgecolor=c, alpha=0.4, label = "d"+d, marker = m2, markersize=5, linestyle='')
					else:
						posPlot.axs[0].plot(i_pos,float(P2data[P2data['POSITION'] == i_pos]['AF_d'+d]),
							color=c, markeredgecolor=c, alpha=0.4, label = "d"+d, marker = m2, markersize=5, linestyle='')
					m2 = next(markerpos)
					j = j+1
			else:
				if i_pos < 10000:
					t = 0
				else:
					t = 1
				if i_pos in [23027, 22363]:
					print("d"+d+" "+str(i_pos) + " " + str(float(1-P2data[P2data['POSITION'] == i_pos]['AF_d'+d])))
					posPlot.axs[t].plot(i_pos,float(1-P2data[P2data['POSITION'] == i_pos]['AF_d'+d]),
						color=c, markeredgecolor=c, alpha=0.4, marker = m2, markersize=5, linestyle='')
				else:
					posPlot.axs[t].plot(i_pos,float(P2data[P2data['POSITION'] == i_pos]['AF_d'+d]),
						color=c, markeredgecolor=c, alpha=0.4, marker = m2, markersize=5, linestyle='')
				m2 = next(markerpos)
				j = j+1
		counter+=1
		if Wuhan_idx:
			ds = [0,28,39,53,75]
			ref = P2data[P2data['POSITION'] == i_pos]['REF'].reset_index(drop=True)[0]
			var = P2data[P2data['POSITION'] == i_pos]['VAR'].reset_index(drop=True)[0]
			if i_pos in [7068, 23562]:
				aa = 'S'
			else:
				aa = 'NS'
			if i_pos in [23027, 22363]: #we plot these in relation to the inferred d0 haplotype
				curr_lab = str(var)+str(int(P2data[P2data['POSITION'] == i_pos]['WUHAN']))+str(ref)+' ('+ aa + ')'
				#print("d"+d+" "+str(i_pos) + " " + str(float(1-P2data[P2data['POSITION'] == i_pos]['AF_d'+d])))
				dat=[]
				q = 0
				for d in ds:
					if d == 0:
						dat.append(0)
						q+=1
					else:
						dat.append( 1 - float(P2data[P2data['POSITION'] == i_pos]['AF_d'+str(d)]))
						q+=1
				#timePlot.plot(days, [1 - float(P2data[P2data['POSITION'] == i_pos]['AF_d'+str(d)]) for d in days], 
				timePlot.plot(ds, dat, 
					label=curr_lab, marker=m1, markersize=5,linestyle=l, alpha=0.4, color = c)
			else:
				dat=[]
				q = 0
				for d in ds:
					if d == 0:
						dat.append(0)
						q+=1
					else:
						dat.append(float(P2data[P2data['POSITION'] == i_pos]['AF_d'+str(d)]))
						q+=1
				curr_lab = str(ref)+str(int(P2data[P2data['POSITION'] == i_pos]['WUHAN']))+str(var)+' ('+ aa + ')'
				#timePlot.plot(days, [float(P2data[P2data['POSITION'] == i_pos]['AF_d'+str(d)]) for d in days], 
				timePlot.plot(ds, dat, 
					label=curr_lab, marker=m1, markersize=5,linestyle=l, alpha=0.4, color = c)
		else:
			ref = P2data[P2data['POSITION'] == i_pos]['REF'].reset_index(drop=True)[0]
			var = P2data[P2data['POSITION'] == i_pos]['VAR'].reset_index(drop=True)[0]
			if i_pos in [7068, 23562]:
				aa = 'S'
			else:
				aa = 'NS'
			curr_lab = str(ref)+str(i_pos)+str(var)+' ('+ aa + ')'
			timePlot.plot(days, [float(P2data[P2data['POSITION'] == i_pos]['AF_d'+str(d)]) for d in days], 
				label=curr_lab, marker=m1, linestyle=l, alpha=0.4, color = c, markersize = 5)

fig.set_size_inches(6,3)

timePlot.axs[0].set_xticks([])
timePlot.axs[1].set_xticks([])
timePlot.axs[0].set_yticks([])
posPlot.axs[0].set_xticks([])
posPlot.axs[0].set_yticks([])

#timePlot.big_ax.set_ylim(-0.05, 1.05)
timePlot.axs[1].set_xlim(min(days)-3, max(days)+2)
timePlot.axs[1].set_xticks(days) 
timePlot.axs[1].set_xticklabels(days, fontsize=10)
timePlot.axs[0].set_xticks([0])
timePlot.axs[0].set_xticklabels([0], fontsize=10)
timePlot.axs[0].set_yticks([0, 0.25, 0.5, 0.75, 1.0]) 
timePlot.axs[0].set_yticklabels(['', '', '', '', '' ])
posPlot.axs[0].set_xticks([1, 7000]) 
posPlot.axs[0].set_xticklabels([1, 7000], fontsize=10)
posPlot.axs[1].set_xticks([20000, 25000]) 
posPlot.axs[1].set_xticklabels([20000, 25000], fontsize=10)
posPlot.axs[0].set_yticks([0, 0.25, 0.5, 0.75, 1.0]) 
posPlot.axs[0].set_yticklabels([0, 0.25, 0.5, 0.75, 1.0])

# posPlot.big_ax.set_ylim(-0.05, 1.05)
# timePlot.big_ax.set_ylim(-0.05, 1.05)
posPlot.axs[0].set_ylim(-0.05, 1.05)
timePlot.axs[0].set_ylim(-0.05, 1.05)
posPlot.axs[1].set_ylim(-0.05, 1.05)
timePlot.axs[1].set_ylim(-0.05, 1.05)

posPlot.axs[1].set_xlabel('position', fontsize=12)
timePlot.axs[1].set_xlabel('day', fontsize=12, y = -.1)

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
markerpos = itertools.cycle(('o', 's', 'v', '^'))

d1 = mlines.Line2D([], [], color='k', marker=next(markerpos), ls='', label='d28', alpha = 0.4)
d2 = mlines.Line2D([], [], color='k', marker=next(markerpos), ls='', label='d39', alpha = 0.4)
d3 = mlines.Line2D([], [], color='k', marker=next(markerpos), ls='', label='d53', alpha = 0.4)
d4 = mlines.Line2D([], [], color='k', marker=next(markerpos), ls='', label='d75', alpha = 0.4)
posPlot.legend(handles=[d1, d2, d3, d4], framealpha = 0.5, loc = 'upper left', bbox_to_anchor=(0.01, 1), fontsize=7)
# print(current_handles)
# [h.set_color('k') for h in current_handles]
# leg = posPlot.get_legend()
# print(leg)
# [leg.legendHandles[i].set_color('k') for i in range(len(days))]
posPlot.axs[0].set_ylabel('iSNV frequency', fontsize=12)

timePlot.legend(title="iSNV position", loc='center left', bbox_to_anchor=(1.05, 0.5))

posPlot.set_title('A)', loc='left')
timePlot.set_title('B)', loc='left')
fig.savefig('Figure-2-P2-high-freq-loci-'+dt.today().isoformat()+'.pdf',bbox_inches='tight')



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


gs = plt.GridSpec(2, 4, figure = fig)#, wspace=-1, hspace=-1, top=0.95, bottom=0.05, left=0.17, right=0.845)
#gs = lib_pos.add_gridspec(1, 4)
# HEIGHT by WIDTH
axt1 = fig.add_subplot(gs[1, 0])
axt2 = fig.add_subplot(gs[1, 1], sharey=axt1)
axt3 = fig.add_subplot(gs[1, 2], sharey=axt1)
axt4 = fig.add_subplot(gs[1, 3], sharey=axt1)

axs=axt1,axt2,axt3,axt4
for ax,dpi,samp in zip(axs, dpi2, samples):
	ax = compareLibs(samples = samp, date=date, dpiString=dpi, complot= ax, colors=new_cols, patient = 'P2',alph = 0.5)
#fig.set_size_inches(3,3)
[ax.tick_params(labelleft=False) for ax in [axt3, axt2, axt4]]
#[ax.tick_params(labelbottom=False) for ax in [axt1, axt2]]

[ax.tick_params(axis = 'x', labelsize=10) for ax in [axt1, axt2, axt3, axt4]]
[ax.tick_params(axis = 'y', labelsize=10) for ax in [axt1, axt3]]
#plt.tight_layout()
# lib_pos.supylabel('L2 frequency', x=-0.05, fontsize=12)
# lib_pos.supxlabel('L1 frequency', y = -0.06, fontsize=12)
axt1.set_ylabel('L2 frequency', x=-0.05, fontsize=12)
#axt1.set_xlabel('L1 frequency', y = -0.06, fontsize=12)
fig.text(0.5, -.02, 'L1 frequency', ha='center', fontsize=12)

directory = '20220121_'
# gs2 = Ct_fig.add_gridspec(1,1)
# Ct_subplot = Ct_fig.add_subplot(gs2[0, 0])
Ct_subplot = fig.add_subplot(gs[0,:])
# axs = [Ct_fig.gca()]
axs = [Ct_subplot]
[ax.set_ylim(-0.05, 1.05) for ax in [axt1, axt2, axt3, axt4]]
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
	[axs[i].text(dpi-3, Ct-3, Ct) for dpi, Ct in zip(dpiNum[i], CtAll[i])]
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
		axs[i].scatter(x=[4],y=[3],marker="1",s=100,color='green')
		axs[i].scatter(x=[86],y=[3],marker="x",s=80,color='k')
	elif labs[i] == 'P3':
		P3.scatter(x=[8],y=[3],marker="1",s=100)
	elif labs[i] == 'P4':
		P6.scatter(x=[0,104],y=[3,3],marker="1",s=100)
	
	axs[i].set_ylim(0, 30)
	i+=1

plt.yticks(ticks=[5, 10, 15, 25, 30], labels = ['5', '', '15', '25', '30'], fontsize=10)
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
Ct_subplot.set_ylim(10, 20)
days = [0,4,28,39,53,75,86]
Ct_subplot.set_xlim(min(days)-5, max(days)+7)
Ct_subplot.set_xticks(days) 
Ct_subplot.set_xticklabels(days, fontsize=10)

#P2.text(0, 25, 'Patient 2', fontsize=12)


axt1.set_title('B)', loc='left', fontsize=12)
Ct_subplot.set_title('A)', loc='left', fontsize=12)

Ct_subplot.set_ylabel('Ct', fontsize=12)
Ct_subplot.set_xlabel('day', y = -0.06, fontsize=12)
fig.align_labels()

fig.set_size_inches(4, 3)

# plt.tight_layout()
del date
from datetime import date

fig.savefig('Figure-1-Treatment-Libs-'+date.today().isoformat()+'.pdf',bbox_inches='tight')













