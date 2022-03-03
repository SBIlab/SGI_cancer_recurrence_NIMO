## network expansion of genes for cell type deconvolution ##
## among LM22 signature genes, first neighbors of seed genes are taken as signature genes
## seed genes are generated from comparison of methylCIBERSORT signature genes and CIBERSORT signature genes
import numpy as np
import scipy.stats as stat
from collections import defaultdict
import os, time
execfile('network_utilities.py', globals())
execfile('pathway_utilities.py', globals())
execfile('parse_genomic_data.py', globals())
gene2uniprot, uniprot2gene = geneID2uniprot(), uniprot2geneID()

def parse_multiomics_signature_genes(coordinate_source):
	'''
	output = [ list of signature genes ]
	output2 = [ list of signature genes in uniprot ID ]
	---------------------
	coordinate_source = 'FANTOM'
	'''
	output, output2 = [], [] 
	fi_dir = '../../result' 
	fiList = os.listdir(fo_dir)
	for fi in fiList:
		if '%s_min_distance_LM22.txt'%coordinate_source == fi:
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



## INITIALIZE
coordinate_sources = ['FANTOM'] # 'FANTOM'
fo_dir = '../../result'

print 'importing network, ', time.ctime()
network = 'string'
if network == 'string':
	network_parameter_dic = {'score_cutoff' : 700 }
	networkDic_tmp, G_tmp, networkName = import_protein_interaction_network_uniprot( network, network_parameter_dic ) # NETWORK IN UNIPROT ID
	network_nodeList, networkDic, G = return_LCC( G_tmp, networkDic_tmp.keys() )



print 'importing LM22 signatures, ', time.ctime()
LM22_signature, LM22_uniprot, LM22_exp, LM22_exp_uniprot = LM22_signature_genes()



## EXPANSION
for coordinate_source in coordinate_sources:
	print '--- expanding genes for %s, '%coordinate_source, time.ctime()
	
	sig_glist, sig_ulist = [], []
	seed_glist, seed_ulist = parse_multiomics_signature_genes(coordinate_source)
	#sig_glist, sig_ulist = seed_glist, seed_ulist
	for seed_u in seed_ulist:
		seed_g = uniprot2gene[seed_u]
		sig_glist.append(seed_g)
		sig_ulist.append(seed_u)
		if seed_u in network_nodeList:
			seed_neighbors = G[seed_u].keys()
			if len( set(seed_neighbors) & set(LM22_uniprot) ) > 0:
				for seed_neighbor in list( set(seed_neighbors) & set(LM22_uniprot) ):
					sig_ulist.append(seed_neighbor)
					sig_glist.append(uniprot2gene[seed_neighbor])
	sig_glist, sig_ulist = list(set(sig_glist)), list(set(sig_ulist))
						
	# 
	fo = open('%s/%s_LM22_%s_network_expansion.txt'%(fo_dir, coordinate_source, network), 'w')
	cell_types = LM22_exp.keys()
	cell_types.sort()
	print >> fo, 'Gene symbol' +  '\t' + '\t'.join(map(str, cell_types))
	for uniprot in sig_ulist:
		gene = uniprot2gene[uniprot]
		tmp = [ gene ]
		for cell_type in cell_types:
			exp = LM22_exp_uniprot[cell_type][uniprot]
			tmp.append(exp)
		print >> fo, '\t'.join(map(str, tmp))
	fo.close()
	print 'expansion complete, ', time.ctime()
