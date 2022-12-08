### transcriptome utilities ###
import pandas as pd
from collections import defaultdict
import numpy as np
import time, os

def ensg_hugo_annotation():
	'''
	output = { ensembl gene id : [HUGO IDs] }
	output2 = { HUGO ID : [ ensembl gene ids ] }
	'''
	output, output2 = defaultdict(set), defaultdict(set)
	df = pd.read_csv('/home/junghokong/PROJECT/co_work/SGI_CRC_seq/data/ensg_HUGO.txt', sep='\t')
	for i in range(len(df)):
		ensg, HGNC = df['Gene stable ID'][i], df['HGNC symbol'][i]
		if (pd.isnull(ensg)==False) and (pd.isnull(HGNC)==False):
			output[ensg].add(HGNC)
			output2[HGNC].add(ensg)
	for tmp_dic in [output, output2]:
		for key in tmp_dic:
			tmp_dic[key] = list(tmp_dic[key])
	return output, output2



def import_transcriptome_data( tumor_only=True, expression_unit='TPM' ):
	'''
	output = { sample : { ensg : exp } } 
	output2 = { sample : { HUGO gene symbol : exp } }
	-----------------
	Input
	tumor_only = True (tumor samples only), False (tumor and normal samples)
	expression_unit = 'TPM', 'FPKM', 'TMM'
	'''
	output, output2 = {}, {}
	ensg2hugo, hugo2ensg = ensg_hugo_annotation()

	exp_dir = '/data/user/junghokong/co_work/SGI_CRC_seq/data/UCSC/GRCh38/RSEM/calculate_expression'

	# FPKM, TPM normalized counts
	if expression_unit.upper() in ['FPKM', 'TPM']:
		for fldr in os.listdir(exp_dir):
			if ('TT' in fldr) or ('TN' in fldr):
				test_continue = True
				if tumor_only == True:
					if 'TN' in fldr:
						test_continue = False
				if test_continue == True:
					sampleID = fldr
					if '%s.genes.results'%sampleID in os.listdir('%s/%s'%(exp_dir, fldr)):
						statinfo = os.stat('%s/%s/%s.genes.results'%(exp_dir,fldr,sampleID))
						if statinfo.st_size > 1000000: # file size greater than 1000000 bytes
							output[sampleID], output2[sampleID] = {}, {}
							f = open('%s/%s/%s.genes.results'%(exp_dir,fldr,sampleID), 'r')
							for line in f.xreadlines():
								line = line.strip().split('\t')
								if (not 'gene_id' in line[0]):
									ensg, TPM, FPKM = line[0], float(line[5]), float(line[6])
									if expression_unit == 'TPM':
										exp = TPM
									elif expression_unit == 'FPKM':
										exp = FPKM

									output[sampleID][ensg] = exp
									if ensg in ensg2hugo:
										for hugo in ensg2hugo[ensg]:
											output2[sampleID][hugo] = exp
							f.close()

	# TMM normalized counts
	if expression_unit.upper() == 'TMM':
		f = open('%s/TMM_normalized_expected_count.txt'%exp_dir, 'r')
		for line in f.xreadlines():
			line = line.strip().split('\t')
			if 'genes' in line[0]:
				samples = line[1:]
			else:
				ensg = line[0]
				expList = map(float, line[1:])
				for exp, sample in zip(expList, samples):
					if not sample in output:
						output[sample] = {}
					output[sample][ensg] = exp
					if ensg in ensg2hugo:
						if not sample in output2:
							output2[sample]= {}
						for hugo in ensg2hugo[ensg]:
							output2[sample][hugo] = exp
		f.close()
	return output, output2
