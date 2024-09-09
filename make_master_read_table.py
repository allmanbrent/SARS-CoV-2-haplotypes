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

directory = '../Data/20220404_'
date = '20220404'

def countISNVdepth(samfile = '',  df = [], sample = '', library = 'L1'):
	# Note that currently this code only looks for iSNVs that all occur together when there are more than 2 iSNVs. It does not look for all combos of concurrent iSNVs
	#header_list = ['READNAME', 'FLAG', 'REFERENCE', 'STARTPOS',  'MAPQ', 'CIGAR', 'RNEXT', 'PNEXT', 'TLEN', 'SEQ', 'QUAL' ]
	print(samfile)
	my_reads = pd.read_csv(samfile, sep='\t', comment='@', header = None, usecols=[0,1,2,3,4,5,6,7,8,9])
	sequences = my_reads[9]
	seqStarts = my_reads[3] # note that these are indexed from "1"
	i = 0
	indic = 0

	pos_depth = {k: 0 for k in df['POSITION']}
	var_depth = {k: 0 for k in df['POSITION']}
	for read, start_pos in zip(sequences, seqStarts):
		#determine if the read spans the loci of interest
		read_len = len(read)
		for k in df['POSITION']:
			if (start_pos <= k) and (read_len + start_pos > k):
				#spans_all_sites += 1
				pos_depth[k] += 1
				#the read is indexed from 0 so we need to reposition the iSNV loci as such
				site_reindexed = k - start_pos
				
				curr_nt = read[site_reindexed] #itemgetter(*sites_reindexed)(read)
				#print(df[df['POSITION'] == k]['VAR-NT'].item())
				if df[df['POSITION'] == k]['VAR-NT'].shape[0] == 2:
					#print('position ' + str(k) + ' has more than one variant')
					#print(df[df['POSITION'] == k])
					if indic == 0:
						indic = 1
						var_depth[k] = [0,0]
					first_nt = df[df['POSITION'] == k]['VAR-NT'].head(1).item()
					secnd_nt = df[df['POSITION'] == k]['VAR-NT'].tail(1).item()
					if first_nt == curr_nt:
						var_depth[k][0] += 1 
					elif secnd_nt == curr_nt:
						var_depth[k][1] += 1 
				elif curr_nt == df[df['POSITION'] == k]['VAR-NT'].item():
					var_depth[k] += 1
	return([pos_depth, var_depth])

def countISNVdepthTriAllelic(samfile = '',  df = [], sample = '', library = 'L1'):
	# Note that currently this code only looks for iSNVs that all occur together when there are more than 2 iSNVs. It does not look for all combos of concurrent iSNVs
	#header_list = ['READNAME', 'FLAG', 'REFERENCE', 'STARTPOS',  'MAPQ', 'CIGAR', 'RNEXT', 'PNEXT', 'TLEN', 'SEQ', 'QUAL' ]
	print(samfile)
	my_reads = pd.read_csv(samfile, sep='\t', comment='@', header = None, usecols=[0,1,2,3,4,5,6,7,8,9])
	sequences = my_reads[9]
	seqStarts = my_reads[3] # note that these are indexed from "1"
	i = 0
	indic = 0

	pos_depth = {(k,v): 0 for k,v in zip(df['POSITION'], df['VAR-NT'])}
	var_depth = {(k,v): 0 for k,v in zip(df['POSITION'], df['VAR-NT'])}
	for read, start_pos in zip(sequences, seqStarts):
		#determine if the read spans the loci of interest
		read_len = len(read)
		for k,v in zip(df['POSITION'], df['VAR-NT']):
			if (start_pos <= k) and (read_len + start_pos > k):
				#spans_all_sites += 1
				pos_depth[(k,v)] += 1
				#the read is indexed from 0 so we need to reposition the iSNV loci as such
				site_reindexed = k - start_pos
				
				curr_nt = read[site_reindexed] #itemgetter(*sites_reindexed)(read)
				if curr_nt == v:
					var_depth[(k,v)] += 1
	return([pos_depth, var_depth])


	#print(s)
	#pd.DataFrame.from_dict(hap_counts).to_csv(sample+'sites-'+s+'-Haplotypes-'+dt.today().isoformat()+'.tsv', sep='\t', index=True, header=True)

	# with open(sample+'-'+L+'-sites-'+s+'-Haplotypes-'+dt.today().isoformat()+'.tsv', 'w') as csvfile:
	# 	writer = csv.writer(csvfile, delimiter = '\t')
	# 	for key in hap_counts.keys():
	# 		#writer.writeheader()
	# 	 	writer.writerow([key, hap_counts[key]])

	# 	csvfile.close()

def compareLibs(samples = [], dpiString=[], dpiNum=[]):
	L1List=pd.read_csv(directory+samples+'/'+samples+'_L1_lofreq-output.txt', sep='\t')
	L2List=pd.read_csv(directory+samples+'/'+samples+'_L2_lofreq-output.txt', sep='\t')
	L1List['VAR-DEPTH'] = L1List['FWD-VAR'] + L1List['REV-VAR']
	L2List['VAR-DEPTH'] = L2List['FWD-VAR'] + L2List['REV-VAR']

	L1List = L1List.rename(columns={"ALLELE-FREQUENCY": "AF"})
	L2List = L2List.rename(columns={"ALLELE-FREQUENCY": "AF"})
	#print(L1List.columns)
	to_keep = ['POSITION', 'REF-NT', 'VAR-NT', 'AF', 'DEPTH', 'VAR-DEPTH']
	for col in ['DEPTH', 'VAR-DEPTH']:
		L1List[col] = L1List[col].astype('int64')
		L2List[col] = L2List[col].astype('int64')

	combined_libs = pd.merge(L1List[to_keep], L2List[to_keep], on = ['POSITION', 'REF-NT', 'VAR-NT'], how = 'outer', 
		suffixes = ['_d'+dpiString+'_L1','_d'+dpiString+'_L2']).fillna(0)
	for col in ['DEPTH_d'+dpiString+'_L1', 'DEPTH_d'+dpiString+'_L2','VAR-DEPTH_d'+dpiString+'_L1', 'VAR-DEPTH_d'+dpiString+'_L2']:
		combined_libs[col] = combined_libs[col].astype('int64')


	return(combined_libs)

samples=['705E', '798T', '962B', '1064Z']
date='20220404'
dpi2=['28','39','53','75']
dpiNum=[28,39,53,75]
combLibs = []
i = 0

for samp in samples:
	combLibs.append(compareLibs(samples = samp, dpiString = dpi2[i], dpiNum = dpiNum[i]))
	i = i+1

print([c['POSITION'] for c in combLibs])
df_final = reduce(lambda left, right: pd.merge_ordered(left, right, on=['POSITION', 'REF-NT', 'VAR-NT'], how='outer'), combLibs).fillna(0)
print(df_final.columns)
print(df_final['POSITION'])
fname = 'preprocessed_iSNV_calls_without_supporting_reads_d28_coord_'+dt.today().isoformat()+'.tsv'
df_final.to_csv(fname, sep="\t", index = False)
i = 0

for s in samples:
	day = dpi2[i]
	i = i+1
	for L in ['L1', 'L2']:
		my_path = '../Data/'+date+'_'+s+'/'+s+'_'+L
		pos_no_depth = df_final[df_final['DEPTH_d'+day+'_'+L] < 1]['POSITION']
		curr_no_depth = df_final[df_final['POSITION'].isin(pos_no_depth)]
		[pos_depth, var_depth] = countISNVdepthTriAllelic(samfile = my_path+'_bbamp_preened.sam', 
			df = curr_no_depth, sample = s, library = L)
		print(pos_depth)
		print(var_depth)
		for r,p in curr_no_depth['POSITION'].iteritems():
			#print(p)
			c = 'DEPTH_d'+day+'_'+L
			v = 'VAR-DEPTH_d'+day+'_'+L
			var_nt = df_final[df_final['POSITION'] == p]['VAR-NT']
			#print(df_final[df_final['POSITION'] == p][['POSITION', 'VAR-NT']])
			#print(df_final[df_final['POSITION'] == p])
			# df_final[(df_final['POSITION'] == p) and (df_final[c] < 1)][c] = pos_depth[p]
			# df_final[(df_final['POSITION'] == p) and (df_final[v] < 1)][v] = var_depth[p]
			# df_final[(df_final['POSITION'] == p)][c] = pos_depth[p]
			# df_final[(df_final['POSITION'] == p)][v] = var_depth[p]
			#print(df_final.loc[(df_final['POSITION'] == p), c])
			#if not(isinstance(var_depth[(p, var_nt)], int)): #if there is more than one variant at this position
			if var_nt.size > 1:
				print(df_final[df_final['POSITION'] == p][['POSITION', 'VAR-NT']])
				#print(df_final.loc[(df_final['POSITION'] == p), c])
				# df_final.loc[(df_final['POSITION'] == p), [c]][0] = var_depth[(p, var_nt)][0] #var_depth[p][0]
				# df_final.loc[(df_final['POSITION'] == p), [c]][1] = var_depth[(p, var_nt)][1]#var_depth[p][1]
				df_final.loc[df_final['POSITION'] == p, [c]][0] = pos_depth[(p,var_nt.iloc[0])]
				df_final.loc[df_final['POSITION'] == p, [c]][1] = pos_depth[(p,var_nt.iloc[1])]
				df_final.loc[df_final['POSITION'] == p, [v]][0] = var_depth[(p,var_nt.iloc[0])]
				df_final.loc[df_final['POSITION'] == p, [v]][1] = var_depth[(p,var_nt.iloc[1])]
			else:
				# df_final.loc[(df_final['POSITION'] == p), c] = pos_depth[p]
				# df_final.loc[(df_final['POSITION'] == p), v] = var_depth[p]
				df_final.loc[df_final['POSITION'] == p, [c]] = pos_depth[(p,var_nt.iloc[0])]
				df_final.loc[df_final['POSITION'] == p, [v]] = var_depth[(p,var_nt.iloc[0])]

print(df_final)
fname = 'iSNV_calls_with_supporting_reads_d28_coord_'+dt.today().isoformat()+'.tsv'
df_final.to_csv(fname, sep="\t", index = False)


