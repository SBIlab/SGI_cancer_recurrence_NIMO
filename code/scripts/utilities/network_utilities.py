#########################################################################
# Network Utility
#
# This code includes
# - create a network
# - filter a network based on degree
# - find paths and (connected) components of given network
# - randomize a network
# - network proximity
#
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

# import regulatory network
def import_regulatory_network_uniprot( tissue = 'bladder_adult', confidence_score = 0.05 ): # provide in uniprot ID
    output = defaultdict(set) # { tf : [ list of target genes ] }
    geneID_uniprot = geneID2uniprot()
    fi_directory = '/data/user/junghokong/data/394_individual_networks'
    fiList = os.listdir(fi_directory)
    for fi in fiList:
        if tissue == fi.split('.')[0]:
            print tissue, fi, time.ctime() # print regulatory network file name
            df = pd.read_csv('%s/%s' %(fi_directory, fi), compression='gzip', sep ='\t', names=['tf','target','score'])
            for i in range(len(df)):
                tf, target, score = df['tf'][i], df['target'][i], df['score'][i]
                if score >= confidence_score:
                    if (tf in geneID_uniprot) and (target in geneID_uniprot):
                        tf_uniprot, target_uniprot = geneID_uniprot[tf], geneID_uniprot[target]
                        output[tf_uniprot].add(target_uniprot)
    for tf in output:
        output[tf] = list(output[tf])
    print 'regulatory network generated for ', tissue, time.ctime()
    return output

def import_regulatory_network_geneID( tissue = 'bladder_adult', confidence_score = 0.05 ): # provide in gene ID
    output = defaultdict(set)
    regnet = import_regulatory_network_uniprot( tissue = 'bladder_adult', confidence_score = 0.05 ) # provide in uniprot ID
    uniprot2gene = uniprot2geneID()
    for key in regnet:
        for value in regnet[key]:
            if ( key in uniprot2gene ) and ( value in uniprot2gene ):
                g1, g2 = uniprot2gene[key], uniprot2gene[value]
                output[g1].add(g2)
    for gene in output:
        output[gene] = list(output[gene])
    return output

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

# MENCHE SCIENCE COMBINED NETWORK
def import_Menche_SCIENCE_network_entrezID():
    output = defaultdict(set)
    fi_directory = '/data/user/junghokong/data/Incomplete_interactome_SCIENCE/data'
    f = open('%s/DataS1_interactome.tsv' %fi_directory, 'r')
    for line in f.xreadlines():
        line = line.strip().split('\t')
        if not '#' in line[0]:
            gene1, gene2 = line[0], line[1]
            output[gene1].add(gene2)
            output[gene2].add(gene1)
    f.close()
    for g in output:
        output[g] = list(output[g])
    return output

# BARABASI 2019 DRUG COMBINATION NETWORK
def import_barabasi_drug_combination_network_entrezID():
    output = defaultdict(list)
    fi_directory = '/home/junghokong/PROJECT/bladder_cancer/code/15_network_based_analysis/data/barabasi_2019_drug_combination'
    f = open('%s/network.txt' %fi_directory, 'r')
    for line in f.xreadlines():
        line = line.strip().split('\t')
        if not 'Protein_A_Entrez_ID' in line[0]:
            p1, p2 = line[0], line[1]
            output[p1].append(p2)
            output[p2].append(p1)
    f.close()
    for p in output:
        output[p] = list(output[p])
    return output

def import_barabasi_drug_combination_network_geneID():
    output = defaultdict(list)
    entrez_network = import_barabasi_drug_combination_network_entrezID()
    entrez2gene = entrez2geneID()
    for e1 in entrez_network:
        for e2 in entrez_network[e1]:
            if (e1 in entrez2gene) and (e2 in entrez2gene):
                g1, g2 = entrez2gene[e1], entrez2gene[e2]
                output[g1].append(g2)
                output[g2].append(g1)
    for g in output:
        output[g] = list(set(output[g]))
    return output

def import_barabasi_drug_combination_network_uniprotID():
    output = defaultdict(list)
    gene_network = import_barabasi_drug_combination_network_geneID()
    gene2uniprot = geneID2uniprot()
    for g1 in gene_network:
        for g2 in gene_network[g1]:
            if (g1 in gene2uniprot) and (g2 in gene2uniprot):
                u1, u2 = gene2uniprot[g1], gene2uniprot[g2]
                output[u1].append(u2)
                output[u2].append(u1)
    for g in output:
        output[g] = list(set(output[g]))
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

## Network separation scores
# =============================================================================
def get_pathlengths_for_single_set(G,given_gene_set):
    
    """
    calculate the shortest paths of a given set of genes in a
    given network. The results are stored in a dictionary of
    dictionaries:
    all_path_lenghts[gene1][gene2] = l
    with gene1 < gene2, so each pair is stored only once!

    PARAMETERS:
    -----------
        - G: network
        - gene_set: gene set for which paths should be computed

    RETURNS:
    --------
        - all_path_lenghts[gene1][gene2] = l for all pairs of genes
          with gene1 < gene2

    """ 

    # remove all nodes that are not in the network
    all_genes_in_network = set(G.nodes())
    gene_set = given_gene_set & all_genes_in_network

    all_path_lenghts = {}
    
    # calculate the distance of all possible pairs
    for gene1 in gene_set:
        if not all_path_lenghts.has_key(gene1):
            all_path_lenghts[gene1] = {}
        for gene2 in gene_set:
            if gene1 < gene2:
                try:
                    l = nx.shortest_path_length(G, source=gene1, target=gene2)
                    all_path_lenghts[gene1][gene2] = l
                except:
                    continue

    return all_path_lenghts



# =============================================================================
def get_pathlengths_for_two_sets(G,given_gene_set1,given_gene_set2):
    
    """
    calculate the shortest paths between two given set of genes in a
    given network. The results are stored in a dictionary of
    dictionaries: all_path_lenghts[gene1][gene2] = l with gene1 <
    gene2, so each pair is stored only once!

    PARAMETERS:
    -----------
        - G: network
        - gene_set1/2: gene sets for which paths should be computed

    RETURNS:
    --------
        - all_path_lenghts[gene1][gene2] = l for all pairs of genes
          with gene1 < gene2

    """ 

    # remove all nodes that are not in the network
    all_genes_in_network = set(G.nodes())
    gene_set1 = given_gene_set1 & all_genes_in_network
    gene_set2 = given_gene_set2 & all_genes_in_network

    all_path_lenghts = {}
    
    # calculate the distance of all possible pairs
    for gene1 in gene_set1:
        if not all_path_lenghts.has_key(gene1):
            all_path_lenghts[gene1] = {}
        for gene2 in gene_set2:
            if gene1 != gene2:
                try:
                    l = nx.shortest_path_length(G, source=gene1, target=gene2)
                    if gene1 < gene2:
                        all_path_lenghts[gene1][gene2] = l
                    else:
                        if not all_path_lenghts.has_key(gene2):
                            all_path_lenghts[gene2] = {}
                        all_path_lenghts[gene2][gene1] = l
                except:
                    continue

    return all_path_lenghts

def calc_single_set_distance(G,given_gene_set):

    """
    Calculates the mean shortest distance for a set of genes on a
    given network    
    

    PARAMETERS:
    -----------
        - G: network
        - gene_set: gene set for which distance will be computed 

    RETURNS:
    --------
         - mean shortest distance 

    """


    # remove all nodes that are not in the network, just to be safe
    all_genes_in_network = set(G.nodes())
    gene_set = given_gene_set & all_genes_in_network

    # get the network distances for all gene pairs:
    all_path_lenghts = get_pathlengths_for_single_set(G,gene_set)

    all_distances = []

    # going through all gene pairs
    for geneA in gene_set:

        all_distances_A = []
        for geneB in gene_set:

            # I have to check which gene is 'smaller' in order to know
            # where to look up the distance of that pair in the
            # all_path_lengths dict
            if geneA < geneB:
                if all_path_lenghts[geneA].has_key(geneB):
                    all_distances_A.append(all_path_lenghts[geneA][geneB])
            else:
                if all_path_lenghts[geneB].has_key(geneA):
                    all_distances_A.append(all_path_lenghts[geneB][geneA])

        if len(all_distances_A) > 0:
            l_min = min(all_distances_A)
            all_distances.append(l_min)

    # calculate mean shortest distance
    mean_shortest_distance = np.mean(all_distances)

    return mean_shortest_distance


# =============================================================================
def calc_set_pair_distances(G,given_gene_set1,given_gene_set2):

    """
    Calculates the mean shortest distance between two sets of genes on
    a given network
    
    PARAMETERS:
    -----------
        - G: network
        - gene_set1/2: gene sets for which distance will be computed 

    RETURNS:
    --------
         - mean shortest distance 

    """

    # remove all nodes that are not in the network
    all_genes_in_network = set(G.nodes())
    gene_set1 = given_gene_set1 & all_genes_in_network
    gene_set2 = given_gene_set2 & all_genes_in_network

    # get the network distances for all gene pairs:
    all_path_lenghts = get_pathlengths_for_two_sets(G,gene_set1,gene_set2)

    all_distances = []

    # going through all pairs starting from set 1 
    for geneA in gene_set1:

        all_distances_A = []
        for geneB in gene_set2:

            # the genes are the same, so their distance is 0
            if geneA == geneB:
                all_distances_A.append(0)
                
            # I have to check which gene is 'smaller' in order to know
            # where to look up the distance of that pair in the
            # all_path_lengths dict
            else:
                if geneA < geneB:
                    try:
                        all_distances_A.append(all_path_lenghts[geneA][geneB])
                    except:
                        pass

                else:
                    try:
                        all_distances_A.append(all_path_lenghts[geneB][geneA])
                    except:
                        pass


        if len(all_distances_A) > 0:
            l_min = min(all_distances_A)
            all_distances.append(l_min)

    # going through all pairs starting from disease B
    for geneA in gene_set2:

        all_distances_A = []
        for geneB in gene_set1:

            # the genes are the same, so their distance is 0
            if geneA == geneB:
                all_distances_A.append(0)

            # I have to check which gene is 'smaller' in order to know
            # where to look up the distance of that pair in the
            # all_path_lengths dict
            else:
                if geneA < geneB:
                    try:
                        all_distances_A.append(all_path_lenghts[geneA][geneB])
                    except:
                        pass
                        
                else:
                    try:
                        all_distances_A.append(all_path_lenghts[geneB][geneA])
                    except:
                        pass

        if len(all_distances_A) > 0:
            l_min = min(all_distances_A)
            all_distances.append(l_min)


    # calculate mean shortest distance
    mean_shortest_distance = np.mean(all_distances)

    return mean_shortest_distance


def get_separation(network, nodes1, nodes2):
    """
    s_AB = d_AB - (d_AA + d_BB)/2
    """
    d_AA = calc_single_set_distance(network,set(nodes1))
    d_BB = calc_single_set_distance(network,set(nodes2))
    d_AB = calc_set_pair_distances(network, set(nodes1), set(nodes2))
    s_AB = d_AB - (d_AA + d_BB)/2.
    return s_AB


## Network distances and proximity calculation
def get_shortest_path_lengths(G):
    return networkx.shortest_path_length(G)

def get_shortest_paths(G):
    return networkx.shortest_path(G)

def get_shortest_path_between(G, source_id, target_id):
    return networkx.shortest_path(G, source_id, target_id)

def get_shortest_path_length_between(G, source_id, target_id):
    return networkx.shortest_path_length(G, source_id, target_id)

def get_shortest_path_from_source_to_gene_group(G, source_id, geneGroup):
    output = []
    tmp, tmpLength = [], []
    if not type(geneGroup) == list:
        geneGroup = [geneGroup]
    for targetGene in geneGroup:
        try:
            tmp.append(get_shortest_path_between(G, source_id, targetGene))
            tmpLength.append(get_shortest_path_length_between(G, source_id, targetGene))
        except:
            pass
    if len(tmpLength) > 0:
        minLength = np.min(tmpLength)
        min_distance = minLength
        for pathList, pathLength in zip(tmp, tmpLength):
            if pathLength == minLength:
                output.append(pathList)
    else:
        output, min_distance = None, None
    return output, min_distance
        
        

def calculate_closest_distance(network, nodes_from, nodes_to, lengths=None):
    """
    Average shortest path length of nodes_from and nodes_to
    """
    
    values_outer = []
    if lengths is None:
        for node_from in nodes_from:
            values = []
            for node_to in nodes_to:
                val = get_shortest_path_length_between(network, node_from, node_to)
                values.append(val)
            d = min(values)
            #print d,
            values_outer.append(d)
    else:
        for node_from in nodes_from:
            values = []
            vals = lengths[node_from]
            for node_to in nodes_to:
                val = vals[node_to]
                values.append(val)
            d = min(values)
            values_outer.append(d)
    d = numpy.mean(values_outer)
    #print d
    return d


def calculate_shortest_distance(network, nodes_from, nodes_to, lengths = None):
    """
    Shortest distance to all disease genes
    """

    values = []
    if lengths is None:
        for node_from in nodes_from:
            num_nodes_to = len(nodes_to)
            value = 0
            for node_to in nodes_to:
                d = get_shortest_path_length_between(network, node_from, node_to)
                value += d
            values.append( value/float(num_nodes_to) )
    d = np.sum(values)/float(len(nodes_from))
    return d

def calculate_kernel_distance(network, nodes_from, nodes_to, lengths = None):
    """
    Kernel distance
    exponential penalty to further distances
    """
    LN_values = []
    if lengths is None:
        for node_from in nodes_from:
            num_nodes_to = len(nodes_to)
            value = 0
            for node_to in nodes_to:
                d = get_shortest_path_length_between(network, node_from, node_to)
                value += ( np.power(np.e, (-1*(d+1))) )/float(num_nodes_to)
            LN_value = np.log(value)
            LN_values.append(LN_value)
        kernel_distance = -1*(np.sum(LN_values))/float(len(nodes_from))
    return kernel_distance


def get_shortest_path_length_between(G, source_id, target_id):
    return nx.shortest_path_length(G, source_id, target_id)

def get_degree_binning(g, bin_size, lengths=None):
    degree_to_nodes = {}
    for node, degree in g.degree().iteritems(): #.iteritems(): # iterator in networkx 2.0
        if lengths is not None and node not in lengths:
            continue
        degree_to_nodes.setdefault(degree, []).append(node)
    values = degree_to_nodes.keys()
    values.sort()
    bins = []
    i = 0
    while i < len(values):
        low = values[i]
        val = degree_to_nodes[values[i]]
        while len(val) < bin_size:
            i += 1
            if i == len(values):
                break
            val.extend(degree_to_nodes[values[i]])
        if i == len(values):
            i -= 1
        high = values[i]
        i += 1 
        #print i, low, high, len(val) 
        if len(val) < bin_size:
            low_, high_, val_ = bins[-1]
            bins[-1] = (low_, high, val_ + val)
        else:
            bins.append((low, high, val))
    return bins

def get_degree_equivalents(seeds, bins, g):
    seed_to_nodes = {}
    for seed in seeds:
        d = g.degree(seed)
        for l, h, nodes in bins:
            if l <= d and h >= d:
                mod_nodes = list(nodes)
                mod_nodes.remove(seed)
                seed_to_nodes[seed] = mod_nodes
                break
    return seed_to_nodes
       

def pick_random_nodes_matching_selected(network, bins, nodes_selected, n_random, degree_aware=True, connected=False, seed=None):
    """
    Use get_degree_binning to get bins
    """
    if seed is not None:
        random.seed(seed)
    values = []
    nodes = network.nodes()
    for i in xrange(n_random):
        if degree_aware:
            if connected:
                raise ValueError("Not implemented!")
            nodes_random = set()
            node_to_equivalent_nodes = get_degree_equivalents(nodes_selected, bins, network)
            for node, equivalent_nodes in node_to_equivalent_nodes.iteritems():
                #nodes_random.append(random.choice(equivalent_nodes))
                chosen = random.choice(equivalent_nodes)
                for k in xrange(20): # Try to find a distinct node (at most 20 times)
                    if chosen in nodes_random:
                        chosen = random.choice(equivalent_nodes)
                nodes_random.add(chosen)
            nodes_random = list(nodes_random)
        else:
            if connected:
                nodes_random = [ random.choice(nodes) ]
                k = 1
                while True:
                    if k == len(nodes_selected):
                        break
                    node_random = random.choice(nodes_random)
                    node_selected = random.choice(network.neighbors(node_random))
                    if node_selected in nodes_random:
                        continue
                    nodes_random.append(node_selected)
                    k += 1
            else:
                nodes_random = random.sample(nodes, len(nodes_selected))
        values.append(nodes_random)
    return values

def get_random_nodes(nodes, network, bins=None, n_random=1000, min_bin_size=100, degree_aware=True, seed=None):
    if bins is None:
        # Get degree bins of the network
        bins = get_degree_binning(network, min_bin_size) 
    nodes_random = pick_random_nodes_matching_selected(network, bins, nodes, n_random, degree_aware, seed=seed) 
    return nodes_random

#============================
##### Proximity related #####

def calculate_proximity(network, nodes_from, nodes_to, nodes_from_random=None, nodes_to_random=None, bins=None, n_random=1000, min_bin_size=100, seed=452456, lengths=None): # closest & closest_inverted
    """
    Calculate proximity from nodes_from to nodes_to
    If degree binning or random nodes are not given, they are generated
    lengths: precalculated shortest path length dictionary
    """
    #distance = "closest"
    #lengths = get_shortest_path_lengths(network, "../data/toy.sif.pcl")
    #d = get_separation(network, lengths, nodes_from, nodes_to, distance, parameters = {})
    nodes_network = set(network.nodes())
    nodes_from = set(nodes_from) & nodes_network 
    nodes_to = set(nodes_to) & nodes_network
    if len(nodes_from) == 0 or len(nodes_to) == 0:
        return None # At least one of the node group not in network
    d = calculate_closest_distance(network, nodes_from, nodes_to, lengths)
    if bins is None and (nodes_from_random is None or nodes_to_random is None):
        bins = get_degree_binning(network, min_bin_size, lengths) # if lengths is given, it will only use those nodes
    if nodes_from_random is None:
        nodes_from_random = get_random_nodes(nodes_from, network, bins = bins, n_random = n_random, min_bin_size = min_bin_size, seed = seed)
    if nodes_to_random is None:
        nodes_to_random = get_random_nodes(nodes_to, network, bins = bins, n_random = n_random, min_bin_size = min_bin_size, seed = seed)
    random_values_list = zip(nodes_from_random, nodes_to_random)
    values = numpy.empty(len(nodes_from_random)) #n_random
    for i, values_random in enumerate(random_values_list):
        nodes_from, nodes_to = values_random
        #values[i] = get_separation(network, lengths, nodes_from, nodes_to, distance, parameters = {})
        values[i] = calculate_closest_distance(network, nodes_from, nodes_to, lengths)
    #pval = float(sum(values <= d)) / len(values) # needs high number of n_random
    m, s = numpy.mean(values), numpy.std(values)
    if s == 0:
        z = 0.0
    else:
        z = (d - m) / s
    return d, z, m, s, values #(z, pval)


def calculate_shortest_proximity(network, nodes_from, nodes_to, nodes_from_random=None, nodes_to_random=None, bins=None, n_random=1000, min_bin_size=100, seed=452456, lengths=None):
    """
    Calculate shortest proximity from nodes_from to nodes_to
    nodes_from = drug targets // nodes_to = disease genes (subtyping genes)
    """
    nodes_network = set(network.nodes())
    nodes_from = set(nodes_from) & nodes_network
    nodes_to = set(nodes_to) & nodes_network
    if len(nodes_from) == 0 or len(nodes_to) == 0:
        return None # At least one of the node group not in network
    # calculate observed shortest distance proximity
    d = calculate_shortest_distance(network, nodes_from, nodes_to, lengths)

    # random shortest distance proximity
    if bins is None and (nodes_from_random is None or nodes_to_random is None):
        bins = get_degree_binning(network, min_bin_size, lengths) # if lengths is given, it will only use those nodes
    if nodes_from_random is None:
        nodes_from_random = get_random_nodes(nodes_from, network, bins = bins, n_random = n_random, min_bin_size = min_bin_size, seed = seed)
    if nodes_to_random is None:
        nodes_to_random = get_random_nodes(nodes_to, network, bins = bins, n_random = n_random, min_bin_size = min_bin_size, seed = seed)
    random_values_list = zip(nodes_from_random, nodes_to_random)
    values = numpy.empty(len(nodes_from_random)) #n_random
    for i, values_random in enumerate(random_values_list):
        nodes_from, nodes_to = values_random
        #values[i] = get_separation(network, lengths, nodes_from, nodes_to, distance, parameters = {})
        values[i] = calculate_shortest_distance(network, nodes_from, nodes_to, lengths)
        
    #pval = float(sum(values <= d)) / len(values) # needs high number of n_random
    m, s = numpy.mean(values), numpy.std(values)
    if s == 0:
        z = 0.0
    else:
        z = (d - m) / s
    return d, z, m, s, values #(z, pval)


def calculate_kernel_proximity(network, nodes_from, nodes_to, nodes_from_random=None, nodes_to_random=None, bins=None, n_random=1000, min_bin_size=100, seed=452456, lengths=None):
    """
    calculate kernel proximity
    """
    nodes_network = set(network.nodes())
    nodes_from = set(nodes_from) & nodes_network
    nodes_to = set(nodes_to) & nodes_network
    if len(nodes_from) == 0 or len(nodes_to) == 0:
        return None # At least one of the node group not in network

    # kernel distance
    d = calculate_kernel_distance(network, nodes_from, nodes_to, lengths)

    # random kernel distance proximity
    if bins is None and (nodes_from_random is None or nodes_to_random is None):
        bins = get_degree_binning(network, min_bin_size, lengths) # if lengths is given, it will only use those nodes
    if nodes_from_random is None:
        nodes_from_random = get_random_nodes(nodes_from, network, bins = bins, n_random = n_random, min_bin_size = min_bin_size, seed = seed)
    if nodes_to_random is None:
        nodes_to_random = get_random_nodes(nodes_to, network, bins = bins, n_random = n_random, min_bin_size = min_bin_size, seed = seed)
    random_values_list = zip(nodes_from_random, nodes_to_random)
    values = numpy.empty(len(nodes_from_random)) #n_random
    for i, values_random in enumerate(random_values_list):
        nodes_from, nodes_to = values_random
        values[i] = calculate_kernel_distance(network, nodes_from, nodes_to, lengths)

    m, s = np.mean(values), np.std(values)
    if s == 0:
        z = 0.0
    else:
        z = (d - m) / s
    return d, z, m, s, values

    
    
    

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
