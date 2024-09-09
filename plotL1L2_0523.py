#!/usr/bin/env python

import pandas as pd
import sys
import numpy as np
import os
from datetime import date as dt
import seaborn as sns 
import itertools
from matplotlib.lines import Line2D
import matplotlib.artist as artists
import matplotlib.pyplot as plt
from functools import reduce

cols = ['#ffff00', '#00d48f', '#008eac', '#3f475d']
#directory = '20220121_'
#directory = '20220210_'
directory = '../Data/20220404_'

# fig, ((axp1, axt1),(axp2,axt2),(axp3, axt3),(axp6,axt6), (axp7, axt7)) = plt.subplots(nrows=5,ncols=2, sharex='col',
# 					gridspec_kw={
#                            'width_ratios': [3.5, 2.5],
#                            'height_ratios': [2, 2,2,2,2],
#                        'wspace': 0.16,
#                        'hspace': 0.1})
# axs=(axp1, axt1),(axp2,axt2),(axp3, axt3),(axp6,axt6), (axp7, axt7)

fig, ((axp2,axt2),(axp3, axt3),(axp6,axt6)) = plt.subplots(nrows=3,ncols=2, #sharex='col',
                    gridspec_kw={
                           'width_ratios': [3.5, 2.5],
                           'height_ratios': [2,2,2],
                       'wspace': 0.16,
                       'hspace': 0.1})
axs=(axp2,axt2),(axp3, axt3),(axp6,axt6)


fig.suptitle('')
plt.rc('legend',fontsize=8, loc='upper right')

#THIS WONT WORK IN GENERAL YET
# def freqVtimeDF(subs_pos, freqs, dpi, samples,patient_str):
# 	data=pd.DataFrame()

# 	for i_pos in subs_pos.POS:
# 		row = pd.DataFrame({'POS': [i_pos], dpi[0]: freqs[(samples[0], i_pos)]})
# 		data=data.append(row, ignore_index = True)
# 	data.to_csv('table_freq_over_time_'+patient_str+'.csv', sep='\t', index=True, header=True)

figNum, axnum = plt.subplots(constrained_layout=True)
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

  complot.scatter(combined_libs['ALLELE-FREQUENCY_L1'], combined_libs['ALLELE-FREQUENCY_L2'], 
    color=colors, edgecolor='#333333', alpha=alph, s = 10, linewidth = 0.6)
  complot.set_ylim(-0.05, 1.05)
  complot.set_xlim(-0.05, 1.05)
  #complot.set_xlabel('L1 frequency')
  #complot.set_ylabel('L2 frequency')
  complot.axhline(y=0.02, color='red', linestyle='--', lw=1)
  complot.axvline(x=0.02, color='red', linestyle='--', lw=1)
  #complot.set_title(patient + ' d_'+ dpiString + ' sample ' + samples, loc = 'left')
  #complot.set_title(patient + ' d'+ dpiString, loc = 'left')
  complot.text(0.1, 0.9, 'd'+dpiString, size=12, color='k')
  xpoints = ypoints = ax.get_xlim()
  complot.plot(xpoints, ypoints, linestyle='--', color='k', lw=1, scalex=False, scaley=False, alpha=0.4)
  return(complot)
    

################
## PATIENT 1 ###
################
# samples=['706F']
# date='20211104'
# i=0
# dpi=['30']
# First positive sample taken 20210101
# Symptom onset 20201229
# Sample     Patient     Date
# 706F    IC_patient1 20210131

# [axp1, axt1] = freqVtimePlot(samples = ['706F'], date = '20211104', dpiString = ['30'], 
#     dpiNum = [30], patientNum = 'P1', posPlot = axp1, timePlot = axt1, colors = cols)

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

fig, ((axt1,axt2),(axt3,axt4)) = plt.subplots(nrows=2,ncols=2, sharex='col',
                    gridspec_kw={
                           'width_ratios': [1, 1],
                           'height_ratios': [1,1],
                       'wspace': 0.16,
                       'hspace': 0.1})
axs=axt1,axt2,axt3,axt4
for ax,dpi,samp,c in zip(axs, dpi2, samples, numcols[0:4]):
	ax = compareLibs(samples = samp, date=date, dpiString=dpi, complot= ax, colors=c, patient = 'P2',alph = 0.5)
fig.set_size_inches(3,3)
[ax.tick_params(labelleft=False) for ax in [axt2, axt4]]
[ax.tick_params(labelbottom=False) for ax in [axt1, axt2]]


#plt.tight_layout()
del date
from datetime import date
fig.supylabel('L2 frequency', x=-0.05)
fig.supxlabel('L1 frequency', y = -0.03)
#plt.setp(axt1, ylabel='L2 frequency')
fig.savefig('P2-L1-L2-'+date.today().isoformat()+'.pdf',bbox_inches='tight')

################
## PATIENT 3 ###
################

cols = ['#ffff00', '#00d48f', '#008eac', '#3f475d', '#000000']
samples=['760H', '769Q', '997K', '1034V', '1100J']
date='20220210'
dpi3=['51','56','79','91','101']
# First positive sample either 20201215 or 20210110 (unclear from 111021 covid timeline table)
# Symptom onset unsure. 20210110?
# mAb d8
 # Sample     Patient     Date
 # 760H   IC_patient3 20210204
 # 769Q   IC_patient3 20210209
 # 997K   IC_patient3 20210304
 # 1034V  IC_patient3 20210316
 # 1100J  IC_patient3 20210326

# fig, (axt1,axt2,axt3,axt4, axt5) = plt.subplots(nrows=1,ncols=5, #sharex='col',
#                     gridspec_kw={
#                            'width_ratios': [1, 1, 1, 1, 1],
#                            'height_ratios': [1],
#                        'wspace': 0.16,
#                        'hspace': 0.1})
# axs=axt1,axt2,axt3,axt4,axt5
# for ax,dpi,samp,c in zip(axs, dpi3, samples, numcols[0:5]):
# 	ax = compareLibs(samples = samp, date=date, dpiString=dpi, complot= ax, colors=c, patient = 'P3',alph = 0.5)
# [ax.tick_params(labelleft=False) for ax in [axt2, axt3, axt4, axt5]]

# plt.setp(axt1, ylabel='L2 frequency')
# fig.set_size_inches(10, 10/5)

# #plt.tight_layout()
# del date
# from datetime import date
#fig.savefig('P3-L1-L2-'+date.today().isoformat()+'.pdf',bbox_inches='tight')

################
## PATIENT 6 ###
################
samples=['1063Y', '1215U', '1253G']
dates=['20220210', '20220210', '20220303'] 
dpi6=['77', '100', '109']
# First positive sample 20210101
# Symptom onset 20210101
# CP d0 and d104
# Sample     Patient     Date
# 1063Y   IC_patient6 20210319
# 1215U   IC_patient6 20210411
# 1253G   IC_patient6 20210420

fig, (axt1,axt2,axt3) = plt.subplots(nrows=1,ncols=3, #sharex='col',
                    gridspec_kw={
                           'width_ratios': [1, 1, 1],
                           'height_ratios': [1],
                       'wspace': 0.16,
                       'hspace': 0.1})
# axs=axt1,axt2,axt3,axt4
# for ax,dpi,samp,c,date in zip(axs, dpi6, samples, numcols[0:3], dates):
# 	ax = compareLibs(samples = samp, date=date, dpiString=dpi, complot= ax, colors=c, patient = 'P6',alph = 0.5)
# fig.set_size_inches(4.5, 4.5/3)
# plt.setp(axt1, ylabel='L2 frequency')
# [ax.tick_params(labelleft=False) for ax in [axt2, axt3]]
# #plt.tight_layout()
# del date
# from datetime import date
# fig.savefig('P6-P4-L1-L2-'+date.today().isoformat()+'.pdf',bbox_inches='tight')

################
## PATIENT 7 ###
################
# samples=['1183O']
# date='20211104'
# dpi=['192']
# First positive sample 20200924
# Symptom onset 20200901
# Sample     Patient     Date
# 1183O   IC_patient7 20210421

# [axp7, axt7] = freqVtimePlot(samples = samples, date = date, dpiString = dpi, 
#     dpiNum = [192], patientNum = 'P5', posPlot = axp7, timePlot = axt7, colors = cols)


# #axt1.tick_params(labelleft=False)
# axt2.tick_params(labelleft=False)
# axt3.tick_params(labelleft=False)
# axt6.tick_params(labelleft=False)
# #axt7.tick_params(labelleft=False)


# plt.setp(axp6, xlabel='nt position')

# plt.setp(axt6, xlabel='day')
# axt2.set_xlim(28 - 5, 75 +5)
# axt3.set_xlim(51 - 5, 101 +5)
# axt6.set_xlim(77 - 5, 109 +5)

# fig.supylabel('iSNV frequency')
