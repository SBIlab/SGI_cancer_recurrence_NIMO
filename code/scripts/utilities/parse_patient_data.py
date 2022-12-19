### parse patient data ###

## CAVEAT no. 1
## !!! WARNING !!!
## sample recur/non-recur info has changed, thus updated code was written to account for changes


import numpy as np
import pandas as pd
import time, os
from collections import defaultdict
from datetime import date
execfile('pathway_utilities.py', globals())
ensg2gene, gene2uniprot, uniprot2gene = ensembl2geneID(), geneID2uniprot(), uniprot2geneID()

# -------------------------------------------
## UPDATED CODE // using new cohort info
# -------------------------------------------
def parse_updated_retrospective_cohort_recurrence():
	'''
	output = { sample : recurrence status }
	-------
	sample = '10001625'
	'''
	output = {}
	df = pd.read_csv('../../data/clinical_SMC_cohort.txt', sep='\t')
	for i in range(len(df)):
		sample, recurrence_status = str(df['sampleID'][i]), df['Recur'][i]
		if not sample in output:
			if recurrence_status == 'R':
				recurrence = 'TRUE'
			elif recurrence_status == 'N':
				recurrence = 'FALSE'
			output[sample] = recurrence
	return output



def parse_updated_GSE_cohort_recurrence():
	'''
	output = { sample : recurrence status }
	-------
	sample = '10001625'
	'''
	output = {}
	df = pd.read_csv('../../data/clinical_GSE107422.txt', sep='\t')
	for i in range(len(df)):
		sample, recurrence_status = str(df['sample_id'][i]), df['Recurrence'][i]
		if not sample in output:
			if recurrence_status == 'R':
				recurrence = 'TRUE'
			elif recurrence_status == 'NR':
				recurrence = 'FALSE'
			else: continue
			output[sample] = recurrence
	return output


