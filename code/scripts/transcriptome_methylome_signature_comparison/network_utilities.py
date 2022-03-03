#########################################################################
# Network Utility
#########################################################################
import networkx as nx
from collections import defaultdict, OrderedDict
import random, os, time, scipy, json
import numpy as np
import scipy.stats as stat
import pandas as pd
import numpy, scipy, pandas, networkx
import csv

## import protein interaction network
def import_protein_interaction_network( network, network_parameter_dic ):
    """
    Return protein-protein interaction network in { source_node : [ list of neighbor nodes ] } and G // G : networkx format
    network : 'string', 'regulatory', 'combined' // 'combined' combines string and tissue-specific regulatory network
    network_parameter_dic : { "score_cutoff" : cutoff score } or { "tissue" : target tissue, "confidence_score" : network confidence score } or { "score_cutoff" : cutoff score, "tissue" : target tissue, "confidence_score" : network confidence score }
    Returns { source_node : [ list of neighbor nodes ] }, G, networkName
    """
    output = defaultdict(list)
    G = nx.Graph()
    if network.upper() == 'STRING':
        if 'score_cutoff' in network_parameter_dic:
            score_cutoff = network_parameter_dic['score_cutoff']
            output = import_STRING_network_geneID( score_cutoff )
            networkName = 'STRING_%s' %(score_cutoff)
        else:
            print 'input network_parameter_dic in { "score_cutoff" : cutoff score } format'
            raise KeyError
    elif network.upper() == 'REGULATORY':
        if ( 'tissue' in network_parameter_dic ) and ( 'confidence_score' in network_parameter_dic ):
            tissue, confidence_score = network_parameter_dic['tissue'], network_parameter_dic['confidence_score']
            output = import_regulatory_network_geneID( tissue, confidence_score )
            networkName = 'REGULATORY_%s_%s' %(tissue, confidence_score)
        else:
            print 'input network_parameter_dic in { "tissue" : target tissue, "confidence_score" : network confidence score } format'
            raise KeyError
    elif network.upper() == 'COMBINED':
        if ( 'score_cutoff' in network_parameter_dic ) and ( 'tissue' in network_parameter_dic ) and ( 'confidence_score' in network_parameter_dic ):
            score_cutoff = network_parameter_dic['score_cutoff']
            tissue, confidence_score = network_parameter_dic['tissue'], network_parameter_dic['confidence_score']
            networkName = 'STRING_%s_REGULATORY_%s_%s_COMBINED' %(score_cutoff, tissue, confidence_score)
            stringnet, regnet = import_STRING_network_geneID( score_cutoff ), import_regulatory_network_geneID( tissue, confidence_score )
            output = stringnet
            for p1 in regnet:
                for p2 in regnet[p1]:
                    if p1 in output:
                        if not p2 in output[p1]:
                            output[p1].append(p2)
                    else:
                        output[p1] = []
                        output[p1].append(p2)
        else:
            print 'input network_parameter_dic in { "score_cutoff" : cutoff score, "tissue" : target tissue, "confidence_score" : network confidence score } format'
            raise KeyError
    for g1 in output:
        for g2 in output[g1]:
            G.add_edge(g1,g2)
            
    return output, G, networkName

def import_protein_interaction_network_uniprot( network, network_parameter_dic ):
    """
    Return protein-protein interaction network in { source_node : [ list of neighbor nodes ] } and G // G : networkx format
    network : 'string', 'regulatory', 'combined', 'barabasi_drug_combination' // 'combined' combines string and tissue-specific regulatory network
    network_parameter_dic : { "score_cutoff" : cutoff score } or { "tissue" : target tissue, "confidence_score" : network confidence score } or { "score_cutoff" : cutoff score, "tissue" : target tissue, "confidence_score" : network confidence score }
    Returns { source_node : [ list of neighbor nodes ] }, G, networkName
    """
    output = defaultdict(list)
    G = nx.Graph()
    if network.upper() == 'STRING':
        if 'score_cutoff' in network_parameter_dic:
            score_cutoff = network_parameter_dic['score_cutoff']
            output = import_STRING_network_uniprot( score_cutoff )
            networkName = 'STRING_%s' %(score_cutoff)
        else:
            print 'input network_parameter_dic in { "score_cutoff" : cutoff score } format'
            raise KeyError
    if network.upper() == 'BARABASI_DRUG_COMBINATION':
        output = import_barabasi_drug_combination_network_uniprotID()
        networkName = network.upper()
        
            
    elif network.upper() == 'REGULATORY':
        if ( 'tissue' in network_parameter_dic ) and ( 'confidence_score' in network_parameter_dic ):
            tissue, confidence_score = network_parameter_dic['tissue'], network_parameter_dic['confidence_score']
            output = import_regulatory_network_uniprot( tissue, confidence_score )
            networkName = 'REGULATORY_%s_%s' %(tissue, confidence_score)
        else:
            print 'input network_parameter_dic in { "tissue" : target tissue, "confidence_score" : network confidence score } format'
            raise KeyError
    elif network.upper() == 'COMBINED':
        if ( 'score_cutoff' in network_parameter_dic ) and ( 'tissue' in network_parameter_dic ) and ( 'confidence_score' in network_parameter_dic ):
            score_cutoff = network_parameter_dic['score_cutoff']
            tissue, confidence_score = network_parameter_dic['tissue'], network_parameter_dic['confidence_score']
            networkName = 'STRING_%s_REGULATORY_%s_%s_COMBINED' %(score_cutoff, tissue, confidence_score)
            stringnet, regnet = import_STRING_network_uniprot( score_cutoff ), import_regulatory_network_uniprot( tissue, confidence_score )
            output = stringnet
            for p1 in regnet:
                for p2 in regnet[p1]:
                    if p1 in output:
                        if not p2 in output[p1]:
                            output[p1].append(p2)
                    else:
                        output[p1] = []
                        output[p1].append(p2)
        else:
            print 'input network_parameter_dic in { "score_cutoff" : cutoff score, "tissue" : target tissue, "confidence_score" : network confidence score } format'
            raise KeyError
    for g1 in output:
        for g2 in output[g1]:
            G.add_edge(g1,g2)
            
    return output, G, networkName

# import STRING network
def import_STRING_network_uniprot( score_cutoff = 900 ): # { uniprot ID : [ list of first neighbors in uniprot ] }
    output = defaultdict(set)
    ensembl_geneID = ensembl2geneID()
    geneID_uniprot = geneID2uniprot()
    fi_directory = '/home/junghokong/PROJECT/bladder_cancer/code/9_RegulatoryNetwork/data/9606.protein.links.v10.5.txt'
    f = open('%s/string_gene_links_v10.5_20180806.txt' %fi_directory, 'r')
    for line in f.xreadlines():
        if not 'protein1' in line:
            line = line.strip().split('\t')
            ens_g1, ens_g2, score = line[0], line[1], float(line[2])
            if score >= score_cutoff:
                if (ens_g1 in ensembl_geneID) and (ens_g2 in ensembl_geneID):
                    g1, g2 = ensembl_geneID[ens_g1], ensembl_geneID[ens_g2]
                    if (g1 in geneID_uniprot) and (g2 in geneID_uniprot):
                        uniprot1, uniprot2 = geneID_uniprot[g1], geneID_uniprot[g2]
                        output[uniprot1].add(uniprot2)
                        output[uniprot2].add(uniprot1)
    f.close()
    for g in output:
        output[g] = list(output[g])
    return output

def import_STRING_network_geneID( score_cutoff = 900 ): # { uniprot ID : [ list of first neighbors in uniprot ] }
    output = defaultdict(set)
    stringnet = import_STRING_network_uniprot( score_cutoff = 900 ) # { uniprot ID : [ list of first neighbors in uniprot ] }
    uniprot2gene = uniprot2geneID()
    for key in stringnet:
        for value in stringnet[key]:
            if ( key in uniprot2gene ) and ( value in uniprot2gene ):
                g1, g2 = uniprot2gene[key], uniprot2gene[value]
                output[g1].add(g2)
    for gene in output:
        output[gene] = list(output[gene])
    return output

        

# RETURN LARGEST CONNECTED COMPONENT (LCC)
def return_LCC( G, geneList ):
	"""
	INPUT
	G : network
	geneList : list of genes
	------------------------------------------------------------------------------------------------
	OUTPUT
	Returns a list of genes from geneList that form a largest connected component (LCC) in network G
	Returns node-edge information as dictionary for LCC genes
	Returns GC (Giant Component)
	"""
	edgeDic = defaultdict(list)
	g = nx.Graph()
	if len(geneList) == 0:
		GC = []
	
	else:
		for i, g1 in enumerate(geneList):
			for j, g2 in enumerate(geneList):
				if i<j:
					if G.has_edge(g1,g2) == True:
						g.add_edge(g1,g2)
		if len(g) > 0:
			GC = max(nx.connected_component_subgraphs(g), key=len)
			for key in GC:
				edgeDic[key] = GC[key].keys()
		else:
			GC = g
	return edgeDic.keys(), edgeDic, GC

    
    

## annotation
def ensembl2geneID():
    output = {} # { ensembl : geneID }
    fi_directory = '/home/junghokong/PROJECT/bladder_cancer/code/8_hESC/data'
    df = pd.read_table('%s/2017_07_31_biomart_protein_coding_genes.txt' %fi_directory)
    for i in range(len(df)):
        ensembl, gene = df['Gene stable ID'][i], df['Gene name'][i]
        output[ensembl] = gene
    return output

def geneID2uniprot():
    output = {} # { gene ID : uniprot ID }
    fi_directory = '/home/junghokong/PROJECT/bladder_cancer/code/8_hESC/data'
    df = pd.read_table('%s/uniprot_homoSapiens_multipleGeneName_20180802.tab' %fi_directory)
    for i in range(len(df)):
        uniprot, geneList = df['Entry'][i], df['Gene names'][i]
        if pd.isnull(geneList) == False:
            geneList = geneList.split()
            for gene in geneList:
                output[gene] = uniprot
    return output

def uniprot2geneID():
    output = {} # { uniprot ID : gene ID }
    fi_directory = '/home/junghokong/PROJECT/bladder_cancer/code/8_hESC/data'
    df = pd.read_table('%s/uniprot_homoSapiens_multipleGeneName_20180802.tab' %fi_directory)
    for i in range(len(df)):
        uniprot, geneList = df['Entry'][i], df['Gene names'][i]
        if pd.isnull(geneList) == False:
            geneList = geneList.split()
            gene = geneList[0]
            output[uniprot] = gene
    return output

def entrez2geneID(): # NCBI ID
    output = {} # { entrez : gene ID }
    fi_directory = '/home/junghokong/PROJECT/bladder_cancer/code/15_network_based_analysis/data'
    df = pd.read_table('%s/biomart_ENSG_Uniprot_NCBI_geneName_20190322.txt' %fi_directory)
    for i in range(len(df)):
        entrez, geneName = str(df['NCBI gene ID'][i]).split('.')[0], df['Gene name'][i]
        if (pd.isnull(entrez)==False) and (pd.isnull(geneName)==False):
            output[entrez] = geneName
    return output
