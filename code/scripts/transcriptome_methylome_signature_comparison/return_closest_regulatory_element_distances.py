## return genomic distance to the most proximal regulatory element ##
import numpy as np
import math
import scipy.stats as stat
from collections import defaultdict
import time, os
execfile('parse_genomic_data.py', globals())

## INITIALIZE
fo_dir = '/data/user/junghokong/co_work/SGI_CRC_seq/result/3_transcriptome_methylome_signature_comparison'
LM22_signature, LM22_uniprot, LM22_exp, LM22_exp_uniprot = LM22_signature_genes()
coordinate_sources = ['FANTOM'] # 'FANTOM'


for coordinate_source in coordinate_sources:
	print '\n---> testing %s , '%coordinate_source, time.ctime()

	if coordinate_source == 'FANTOM':
		cpg_map = FANTOM_CpG_mapping() # { chr : [ ( 'promoter'/'enhancer', position_start, position_end, strand, gene, uniprot ) ] }

	
	fo_list = os.listdir(fo_dir)


	# calculate distance - return all distances
	if not '%s_MethylCIBERSORT_Signature_Distance_To_RegulatoryElements.txt'%coordinate_source in fo_list:
		print 'calculating all distances, ', time.ctime()
		m_signature, _ = methylcibersort_reference_signature_methylation_sites() # [ list of signature methylation cites ]
		chr_info = illumina_Chr_info() # { methylation site : ( chromosome, position, strand ) }

		fo = open('%s/%s_MethylCIBERSORT_Signature_Distance_To_RegulatoryElements.txt'%(fo_dir, coordinate_source), 'w')
		print >> fo, '\t'.join(['signature_site', 'signature_chr', 'signature_pos', 'signature_strand','regulatory_element', 'gene_pos_start', 'gene_pos_end', 'gene', 'uniprot', 'distance', 'within_regulatory_element'])

		for m_site in m_signature:
			if m_site in chr_info:
				sig_chr, sig_pos, sig_strand = chr_info[m_site][0], chr_info[m_site][1], chr_info[m_site][2]
				if sig_chr in cpg_map:
					for gene_loc in cpg_map[sig_chr]:
						gene = 'na'
						if coordinate_source == 'FANTOM':
							reg_element, pos_start, pos_end, strand, gene, uniprot = gene_loc[0], gene_loc[1], gene_loc[2], gene_loc[3], gene_loc[4], gene_loc[5]
						# measure distance
						if not gene == 'na':
							distance_list, within_reg_element, equal_strand = [], False, False
							if (not strand == 'na') and (strand == sig_strand):
								equal_strand = True
							if strand == 'na':
								equal_strand = True

							# matching strand
							if equal_strand == True:
								if (sig_pos >= pos_start) and (sig_pos <= pos_end):
									within_reg_element = True
								for pos in [pos_start, pos_end]:
									distance = np.abs( pos - sig_pos )
									distance_list.append(distance)
								min_distance = np.min(distance_list)
								# return output
								output = [ m_site, sig_chr, sig_pos, sig_strand, reg_element, pos_start, pos_end, gene, uniprot, min_distance, within_reg_element ]
								print >> fo, '\t'.join(map(str, output))
		fo.close()
	else:
		print 'we have precalculated distances, ', time.ctime()


	# return minimum distance
	print 'calculating proximal distances, ', time.ctime()
	fo = open('%s/%s_min_distance.txt'%(fo_dir, coordinate_source), 'w')
	print >> fo, '\t'.join(['signature_site', 'signature_chr', 'signature_pos', 'signature_strand', 'regulatory_element', 'gene_pos_start', 'gene_pos_end', 'gene', 'uniprot', 'distance', 'within_regulatory_element'])

	proxiDic = {} # { signature site : ( 'signature_chr', 'signature_pos', 'signature_strand', 'regulatory_element', 'gene_pos_start', 'gene_pos_end', 'gene', 'uniprot', 'distance', 'within_regulatory_element' ) }
	f = open('%s/%s_MethylCIBERSORT_Signature_Distance_To_RegulatoryElements.txt'%(fo_dir, coordinate_source), 'r')
	for line in f.xreadlines():
		line = line.strip().split('\t')
		if not 'signature_site' in line[0]:
			sig_site = line[0]
			tmpList = ( line[1], line[2], line[3], line[4], line[5], line[6], line[7], line[8], line[9], line[10] )
			if sig_site in proxiDic:
				oldDist, newDist = float(proxiDic[sig_site][8]), float(tmpList[8])
				if newDist < oldDist:
					proxiDic[sig_site] = tmpList
			if not sig_site in proxiDic:
				proxiDic[sig_site] = tmpList
	f.close()
	
	for sig_site in proxiDic:
		print >> fo, str(sig_site) + '\t' + '\t'.join(map(str, proxiDic[sig_site]))
	fo.close()


	# methylCIBERSORT signature genes
	m_signature_uniprot_list = []
	for sig_site in proxiDic:
		uniprot = proxiDic[sig_site][7]
		m_signature_uniprot_list.append(uniprot)
	m_signature_uniprot_list = list(set(m_signature_uniprot_list))


	# overlap with LM22 signature gene sets
	print 'overlap with LM22 signature genes & new LM22 reference datasets, ', time.ctime()
	fo = open('%s/%s_min_distance_LM22.txt'%(fo_dir, coordinate_source), 'w')
	cell_types = LM22_exp.keys()
	cell_types.sort()
	print >> fo, 'Gene symbol' +  '\t' + '\t'.join(map(str, cell_types))
	for uniprot in LM22_uniprot:
		gene = uniprot2gene[uniprot]
		if uniprot in m_signature_uniprot_list:
			tmpExpList = []
			for cell_type in cell_types:
				exp = LM22_exp_uniprot[cell_type][uniprot]
				tmpExpList.append(exp)
			print >> fo, gene + '\t' + '\t'.join(map(str, tmpExpList))
	fo.close()

	print 'process completed, ', time.ctime()
