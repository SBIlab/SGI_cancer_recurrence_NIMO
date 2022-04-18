## parse genomic data including : gene chromosome sites, gene promoter and enhancer regions ##
## promoter, enhancer or TSS/TTS region from : FANTOM
## signature methylation sites from 'MethylCIBERSORT' reference datasets
## signature transcriptome genes from 'CIBERSORT' reference datasets

import numpy as np
import scipy.stats as stat
from collections import defaultdict
import time, os
import pandas as pd
execfile('transcriptome_utilities.py', globals())
execfile('pathway_utilities.py', globals())
gene2ensg, ensg2gene = annotation()
gene2uniprot, uniprot2gene = geneID2uniprot(), uniprot2geneID()

def illumina_Chr_info():
	'''
	output = { methylation site : ( chromosome, position, strand ) }
	----------------------------------------
	Returns
	illumina CpG site methylation sites info
	'''
	output = {}
	f = open('/home/webpeace/workbench/SGI_Methyl_deconv/data/illumChrInfo.txt','r')
	for line in f.xreadlines():
		line = line.strip().split('\t')
		if not 'chr' in line[0]:
			output[line[0]] = (line[1], float(line[2]), line[3])
	f.close()
	return output


def methylcibersort_reference_signature_methylation_sites():
	'''
	output = [ list of signature methylation cites ]
	output2 = { cell type : { methylation site : methylation score } }
	'''
	output, output2 = [], {}
	df = pd.read_csv('/data/user/junghokong/co_work/SGI_CRC_seq/data/SGI_COAD_methylation/Total_Signature_Int.txt', sep='\t')
	output = list(df['NAME'])
	for i in range(len(df)):
		m_site = df['NAME'][i]
		for cell_type in df.columns:
			if not 'NAME' in str(cell_type):
				m_score = float(df[cell_type][i])
				if not cell_type in output2:
					output2[cell_type] = {}
				output2[cell_type][m_site] = m_score
	return output, output2


def FANTOM_CpG_mapping():
	'''
	output = { chr : [ ( 'promoter'/'enhancer', position_start, position_end, strand, gene, uniprot ) ] }
	'''
	output = defaultdict(list)
	fi_dir = '/home/webpeace/workbench/SGI_Methyl_deconv/data'
	for reg_element in ['promoter', 'enhancer']:
		f = open('%s/%s---transcript.prec90.txt'%(fi_dir, reg_element),'r')
		for line in f.xreadlines():
			line = line.strip().split('\t')
			if reg_element in line[0]:
				tmp = line
				for clmn_index, clmn in enumerate(line):
					if 'gene' in clmn:
						gene_index = clmn_index
			if not reg_element in line[0]:
				tmp = line[0].split(',')
				if len(tmp) == 1:
					strand = 'na'
				else:
					strand = tmp[1]
				Chr = tmp[0].split('@')[1].split(':')[0]
				position_start = int(tmp[0].split('@')[1].split(':')[1].split('..')[0])
				position_end = int(tmp[0].split('@')[1].split(':')[1].split('..')[1])
				gene = line[gene_index]
				if gene in gene2uniprot:
					uniprot = gene2uniprot[gene]
				else:
					uniprot = 'na'
				
				output[Chr].append( (reg_element, position_start, position_end, strand, gene, uniprot) )
		f.close()
	return output




def LM22_signature_genes():
	'''
	output = [ LM22 signature genes ]
	output2 = [ LM 22 signature genes in uniprot ] 
	output3 = { cell type : { gene : expression } }
	output4 = { cell type : { uniprot : expression } } 
	'''
	output, output2, output3, output4 = [], [], {}, {}
	df = pd.read_csv('/data/user/junghokong/co_work/SGI_CRC_seq/data/LM22.txt', sep='\t')
	output = list(df['Gene symbol'])
	for i in range(len(df)):
		gene = df['Gene symbol'][i]
		for col in df.columns:
			if not str(col) == 'Gene symbol':
				cell_type = col
				exp = float(df[col][i])
				for tmpDic in [output3, output4]:
					if not cell_type in tmpDic:
						tmpDic[cell_type] = {}
				if gene in gene2uniprot:
					uniprot = gene2uniprot[gene]
					output2.append(uniprot)
					output4[cell_type][uniprot] = exp
				output3[cell_type][gene] = exp
	return output, list(set(output2)), output3, output4


def LM22_cell_type_DEGs():
	'''
	output = [ LM 22 signature genes ]
	output2 = [ LM 22 signature genes in uniprot ]
	output3 = { cell type : { gene (cell type specific DEG) : expression } }
	output4 = { cell type : { uniprot (cell type specific DEG) : expression } }
	'''
	output3, output4 = {}, {}
	output, output2, tmp_output3, tmp_output4 = LM22_signature_genes()
	f = open('/data/user/junghokong/co_work/SGI_CRC_seq/data/LM22_cell_type_DEGs.txt', 'r')
	for line in f.xreadlines():
		line = line.strip().split('\t')
		if not '#' in line[0]:
			if 'leukocyte' in line[0]:
				cell_types = line[1:]
			else:
				gene = line[0]
				DEGs = line[1:]
				for index, DEG in enumerate(DEGs):
					cell_type = cell_types[index]
					if not cell_type in output3:
						output3[cell_type] = {}
					if not cell_type in output4:
						output4[cell_type] = {}
					if int(DEG) == 1:
						output3[cell_type][gene] = tmp_output3[cell_type][gene]
						if gene in gene2uniprot:
							uniprot = gene2uniprot[gene]
							output4[cell_type][uniprot] = tmp_output4[cell_type][uniprot]
	f.close()
	return output, output2, output3, output4

def parse_multiomics_signature_genes(coordinate_source='FANTOM'):
	'''
	Returns signature genes from transcriptome + methylome CIBERSORT/MethylCIBERSORT datasets
	output = [ list of signature genes ]
	output2 = [ list of signature genes in uniprot ID ]
	-----------------------
	Input
	coordinate_source = 'FANTOM'
	'''
	output, output2 = [], []
	fi_dir = '/data/user/junghokong/co_work/SGI_CRC_seq/result/3_transcriptome_methylome_signature_comparison'
	fiList = os.listdir(fi_dir)
	for fi in fiList:
		if '%s_LM22_string_network_expansion.txt'%coordinate_source == fi:
			f = open('%s/%s'%(fi_dir, fi), 'r')
			for line in f.xreadlines():
				line = line.strip().split('\t')
				if not 'Gene symbol' in line[0]:
					gene = line[0]
					output.append(gene)
					if gene in gene2uniprot:
						output2.append(gene2uniprot[gene])
			f.close()
			break
	return list(set(output)), list(set(output2))
