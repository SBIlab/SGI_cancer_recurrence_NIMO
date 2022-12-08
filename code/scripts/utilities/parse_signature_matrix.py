## parse immune cell signature matrix
import pandas as pd
from collections import defaultdict
import time, os
import scipy.stats as stat
import numpy as np

def import_signature_genes(sig_type):
	'''
	sig_type = 'LM22', 'LM22_FANTOM', 'LM22_FANTOM_NET'
	'''
	LM22_dir = '/data/user/junghokong/co_work/SGI_CRC_seq/data'
	FANTOM_LM22_dir = '/data/user/junghokong/co_work/SGI_CRC_seq/result/3_transcriptome_methylome_signature_comparison'

	if sig_type == 'LM22':
		df = pd.read_csv('%s/LM22.txt'%LM22_dir, sep='\t')
	else:
		fi_dir = FANTOM_LM22_dir
		if sig_type == 'LM22_FANTOM_NET':
			df = pd.read_csv('%s/FANTOM_LM22_string_network_expansion.txt'%fi_dir, sep='\t')
		if sig_type == 'LM22_FANTOM':
			df = pd.read_csv('%s/FANTOM_min_distance_LM22.txt'%fi_dir, sep='\t')
	return list(df['Gene symbol'])



def import_signature_matrix(sig_type):
	'''
	sig_type = 'LM22', 'LM22_FANTOM', 'LM22_FANTOM_NET'
	'''
	LM22_dir = '/data/user/junghokong/co_work/SGI_CRC_seq/data'
	FANTOM_LM22_dir = '/data/user/junghokong/co_work/SGI_CRC_seq/result/3_transcriptome_methylome_signature_comparison'

	if sig_type == 'LM22':
		df = pd.read_csv('%s/LM22.txt'%LM22_dir, sep='\t')
	else:
		fi_dir = FANTOM_LM22_dir
		if sig_type == 'LM22_FANTOM_NET':
			df = pd.read_csv('%s/FANTOM_LM22_string_network_expansion.txt'%fi_dir, sep='\t')
		if sig_type == 'LM22_FANTOM':
			df = pd.read_csv('%s/FANTOM_min_distance_LM22.txt'%fi_dir, sep='\t')
	return df


def import_LM22_DEG_markers():
	'''
	return LM22 DEG dataset as pandas dataframe
	'''
	fo_dir = '/data/user/junghokong/co_work/SGI_CRC_seq/data'
	df = pd.read_csv('%s/LM22_cell_type_DEGs.txt'%fo_dir, skiprows=3, sep='\t')
	return df

def parse_specific_LM22_DEG_markers():
	'''
	return LM22 DEG dataset as pandas dataframe
	add column with DEG marker specificity
	'''
	df = import_LM22_DEG_markers()
	columns = df.columns
	#col.append('marker_specificity')
	ms_list = []
	for i in range(len(df)):
		count = 0
		for col in df.columns[1:]:
			m_stat = int(df[col][i])
			count += m_stat
		ms_list.append(count)
	output = pd.DataFrame(data=df, columns=columns)
	output['marker_specificity']=ms_list
	return output
