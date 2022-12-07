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
	df = pd.read_csv('/home/junghokong/PROJECT/co_work/SGI_CRC_seq/clinical_final_20200121_DHK.txt', sep='\t')
	for i in range(len(df)):
		sample, recurrence_status = str(df['sampleID'][i]), df['Recur'][i]
		if not sample in output:
			if recurrence_status == 'R':
				recurrence = 'TRUE'
			elif recurrence_status == 'N':
				recurrence = 'FALSE'
			output[sample] = recurrence
	return output

def parse_updated_TCGA_cohort_recurrence():
	'''
	output = { sample : recurrence status }
	-------
	sample = '10001625'
	'''
	output = {}
	df = pd.read_csv('/home/kimdh4962/coPROJECT/for_jhkong/SGI_CRC/data/TCGA/clinical_patient_COAD_from_cBioportal.tsv', sep='\t')
	for i in range(len(df)):
		sample, recurrence_status = str(df['Sample ID'][i]), df['Disease Free Status'][i]
		if not sample in output:
			if recurrence_status == '1:Recurred/Progressed':
				recurrence = 'TRUE'
			elif recurrence_status == '0:DiseaseFree':
				recurrence = 'FALSE'
			else: continue
			output[sample] = recurrence
	return output

def parse_updated_GSE_cohort_recurrence():
	'''
	output = { sample : recurrence status }
	-------
	sample = '10001625'
	'''
	output = {}
	df = pd.read_csv('/home/kwang/PROJECT/samsung_gene/nimo/data/GSE107422_clinical.txt', sep='\t')
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

def parse_updated_GSE12032_cohort_recurrence():
	'''
	output = { sample : recurrence status }
	-------
	sample = '10001625'
	'''
	output = {}
	df = pd.read_csv('/home/kwang/PROJECT/samsung_gene/nimo/data/GSE12032/GSE12032_clinical.txt', sep='\t')
	for i in range(len(df)):
		sample, recurrence_status = str(df['geo_accession'][i]), df['recurrence'][i]
		if not sample in output:
			if recurrence_status == 'R':
				recurrence = 'TRUE'
			elif recurrence_status == 'NR':
				recurrence = 'FALSE'
			else: continue
			output[sample] = recurrence
	return output



def parse_updated_recur_norecur_DEG_analysis_stats():
	'''
	output = { gene : { 'pvalue', 'adj_pvalue', 'log2FC' } }
	output2 = { uniprot : { 'pvalue', 'adj_pvalue', 'log2FC' } } 
	'''
	output, output2 = {}, {}
	f = open('/data/user/junghokong/co_work/SGI_CRC_seq/result/2_transcriptome_analysis/DESeq2_DEG_analysis_updated_cohort_KDH.txt', 'r')
	for line in f.xreadlines():
		line = line.strip().split('\t')
		if (not '#' in line[0]) and (not 'GeneID' in line[0]):
			ensg, log2FC, padj, pvalue = line[0], line[2], line[6], line[5]
			if (not log2FC == 'NA') and (not padj == 'NA') and (not pvalue == 'NA'):
				log2FC, padj = float(log2FC), float(padj)
				if ensg in ensg2gene:
					gene = ensg2gene[ensg]
					if not gene in output:
						output[gene] = {}
					for value, key in zip([pvalue, padj, log2FC], ['pvalue', 'adj_pvalue', 'log2FC']):
						output[gene][key] = value
					# uniprot ID
					if gene in gene2uniprot:
						uniprot = gene2uniprot[gene]
						if not uniprot in output2:
							output2[uniprot] = {}
						for value, key in zip([pvalue, padj, log2FC], ['pvalue', 'adj_pvalue', 'log2FC']):
							output2[uniprot][key] = value
	f.close()
	return output, output2						


def parse_survival_status_days():
	'''
	output = pandas dataframe
	output columns : sample, OpDate, Last follow up date, Death date, days, status at the time of last follow up
	'''
	output = defaultdict(list)
	output_col = ['sample', 'OpDate', 'Last_followup_date', 'Death_date', 'days', 'status']

	df = pd.read_csv('/home/junghokong/PROJECT/co_work/SGI_CRC_seq/200213_clinical_data_survival_JHK.reordered.txt', sep='\t')
	samples = df['sampleID'].tolist()
	OP_dates = df['OpDate'].tolist()
	followup_dates = df['Last follow up date'].tolist()
	death_dates = df['Death date'].tolist()
	status_list = df['status at the time of last follow up'].tolist()
	for sample, op_date, followup_date, death_date, status in zip(samples, OP_dates, followup_dates, death_dates, status_list):
		# start date
		start_year, start_month, start_date = op_date.split('-')
		# end date
		if not 'death' in status.lower():
			end_year, end_month, end_date = followup_date.split('-')
		if 'death' in status.lower():
			if not 'nan' == str(death_date):
				end_year, end_month, end_date = death_date.split('-')
			else:
				print 'No survival data for sample %s'%sample
				print 'raising error, %s, '%time.ctime()
				raise ValueError
		# days
		start_d = date(int(start_year), int(start_month), int(start_date))
		end_d = date(int(end_year), int(end_month), int(end_date))
		delta = end_d - start_d
		days = delta.days
		# output
		output['sample'].append(sample)
		output['OpDate'].append(op_date)
		output['Last_followup_date'].append(followup_date)
		output['Death_date'].append(death_date)
		output['days'].append(days)
		output['status'].append(status)
	output = pd.DataFrame(data=output, columns=output_col)
	return output









# -------------------------------------------
## OUTDATED CODE // using old cohort info
# -------------------------------------------
 
def parse_retrospective_cohort_recurrence():
	'''
	output = { sample : recurrence status }
	-------
	sample = '10001625' 
	'''
	output = {}
	df = pd.read_csv('/home/junghokong/PROJECT/co_work/SGI_CRC_seq/retrospective_sample_info.txt', sep='\t')
	for i in range(len(df)):
		sample, recurrence_status = str(df['ID'][i]).split('_')[0], df['Tumor_recurrence_status'][i]
		if not sample in output:
			if recurrence_status == 'yes':
				recurrence = 'TRUE'
			elif recurrence_status == 'no':
				recurrence = 'FALSE'
			output[sample] = recurrence
	return output


def parse_recur_norecur_DEGs(padj_cutoff, abs_log2FC_cutoff):
	'''
	output = { 'recur_up', 'recur_down' : [ genes ] }
	output2 = { 'recur_up', 'recur_down' : [ genes in uniprot ID ] }
	----------------------------------------
	padj_cutoff = adjusted pvalue cutoff
	abs_log2FC_cutoff = log2 fold change cutoff
	'''
	output, output2 = defaultdict(list), defaultdict(list) # { 'recur_up', 'recur_down' : [ genes ] }
	f = open('/data/user/junghokong/co_work/SGI_CRC_seq/result/2_transcriptome_analysis/DESeq2_DEG_analysis.txt', 'r')
	for line in f.xreadlines():
		line = line.strip().split('\t')
		if (not '#' in line[0]) and (not 'GeneID' in line[0]):
			ensg, log2FC, padj = line[0], line[2], line[6]
			if (not log2FC == 'NA') and (not padj == 'NA'):
				log2FC, padj = float(log2FC), float(padj)
				# DEGs
				if (np.abs(log2FC) >= abs_log2FC_cutoff) and (padj <= padj_cutoff):
					if ensg in ensg2gene:
						# gene ID
						gene = ensg2gene[ensg]
						if log2FC > 0:
							output['recur_up'].append(gene)
						elif log2FC < 0:
							output['recur_down'].append(gene)
						# uniprot ID
						if gene in gene2uniprot:
							uniprot = gene2uniprot[gene]
							if log2FC > 0:
								output2['recur_up'].append(uniprot)
							elif log2FC < 0:
								output2['recur_down'].append(uniprot)
	f.close()
	return output, output2


def parse_recur_norecur_DEG_analysis_stats():
	'''
	output = { gene : { 'pvalue', 'adj_pvalue', 'log2FC' } }
	output2 = { uniprot : { 'pvalue', 'adj_pvalue', 'log2FC' } } 
	'''
	output, output2 = {}, {}
	f = open('/data/user/junghokong/co_work/SGI_CRC_seq/result/2_transcriptome_analysis/DESeq2_DEG_analysis.txt', 'r')
	for line in f.xreadlines():
		line = line.strip().split('\t')
		if (not '#' in line[0]) and (not 'GeneID' in line[0]):
			ensg, log2FC, padj, pvalue = line[0], line[2], line[6], line[5]
			if (not log2FC == 'NA') and (not padj == 'NA') and (not pvalue == 'NA'):
				log2FC, padj = float(log2FC), float(padj)
				if ensg in ensg2gene:
					gene = ensg2gene[ensg]
					if not gene in output:
						output[gene] = {}
					for value, key in zip([pvalue, padj, log2FC], ['pvalue', 'adj_pvalue', 'log2FC']):
						output[gene][key] = value
					# uniprot ID
					if gene in gene2uniprot:
						uniprot = gene2uniprot[gene]
						if not uniprot in output2:
							output2[uniprot] = {}
						for value, key in zip([pvalue, padj, log2FC], ['pvalue', 'adj_pvalue', 'log2FC']):
							output2[uniprot][key] = value
	f.close()
	return output, output2						
