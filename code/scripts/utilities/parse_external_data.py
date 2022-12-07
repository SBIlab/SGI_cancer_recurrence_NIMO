## parse external datasets
import pandas as pd
import numpy as np
import scipy.stats as stat
import os, time
from collections import defaultdict

# FACS + Bulk RNA seq
def FACS_plus_BulkRNAseq(dataset):
	'''
	dataset : 'Linsley', 'Zimmermann'
	---
	input is case insensitive.
	---
	Returns
	output = FACS dataframe
	output2 = expression dataframe
	'''
	fi_dir = '/data/user/junghokong/co_work/SGI_CRC_seq/data/external_dataset'
	if dataset.upper() == 'Linsley'.upper():
		d_name = 'Linsley'
		exp_fiName = 'Linsley_et_al_whole_blood_TMM_normalized_exp.txt'
		FACS_fiName = 'Linsley_cell_proportions_from_Racle_et_al_eLife_2017.txt'
	
	if dataset.upper() == 'Zimmermann'.upper():
		d_name = 'Zimmermann'
		exp_fiName = 'Zimmermann_et_al_tmm_normalized_RNA_seq_prevaccination_samples.txt'
		FACS_fiName = 'Zimmermann_cell_proportions_from_Racle_et_al_eLife_2017.txt'
	
	FACS_df = pd.read_csv('%s/%s_et_al/%s'%(fi_dir, d_name, FACS_fiName), sep='\t')
	exp_df = pd.read_csv('%s/%s_et_al/%s'%(fi_dir, d_name, exp_fiName), sep='\t')
	return FACS_df, exp_df


# IHC + Bulk RNA seq
def IHC_plus_BulkRNAseq(cancerType):
	'''
	cancerType : 'CRC', 'lung', 'melanoma'
	dataset : quanTIseq tumor dataset
	---
	Input is case insensitive.
	---
	output = IHC dataframe
	output2 = expression dataframe
	'''
	fi_dir = '/data/user/junghokong/co_work/SGI_CRC_seq/data/external_dataset/quanTIseq_data'	
	df1 = pd.read_csv('%s/IHC.txt'%fi_dir, sep='\t')
	df2 = pd.read_csv('%s/Leiden_CRC_IF.txt'%fi_dir, sep='\t')
	df = pd.concat([df1, df2], sort=False)
	cancerType2 = cancerType
	hospital = 'Vanderbilt'
	if cancerType.upper() == 'CRC':
		cancerType = 'CRC'
		cancerType2 = 'Colorectal cancer'
		hospital = 'Leiden'
	if cancerType.lower() == 'lung':
		cancerType = 'lung_cancer'
		cancerType2 = 'Lung cancer'
	if cancerType.lower() == 'melanoma':
		cancerType = 'melanoma'
		cancerType2 = 'Melanoma'
	
	output = df.loc[df['Cancer.type']==cancerType2,:]
	output2 = pd.read_csv('%s/%s_%s_TPM.txt'%(fi_dir, hospital, cancerType), sep='\t')
	return output, output2
