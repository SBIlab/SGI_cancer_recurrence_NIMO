################## provide pathway gene list ##################
import pandas as pd
import numpy as np
import scipy.stats as stat
from collections import defaultdict
import os, time



## cancer geneset  // PROVIDING IN GENE IDs
# MutSigDB Hallmark pathway genes
def hallmark_pathway():
    output = defaultdict(list)
    fi_directory = '/home/junghokong/PROJECT/bladder_cancer/code/1_SubtypeSimilarity/data/MSigDB_gene_set'
    f = open('%s/h.all.v6.1.symbols_2017_12_14.txt' %fi_directory, 'r')
    for line in f.readlines():
        line = line.strip().split('\t')
        pathway, geneList = line[0], line[2:]
        output[pathway] = geneList
    f.close()
    return output

def hallmark_pathway_uniprot():
    # CONVERT GENE ID TO UNIPROT FOR HALLMARK GENE SETS
    hallmark = hallmark_pathway()
    hallmark_uniprot = defaultdict(set)
    gene2uniprot = geneID2uniprot()
    for pathway in hallmark:
        for gene in hallmark[pathway]:
            if gene in gene2uniprot:
                uniprot = gene2uniprot[gene]
                hallmark_uniprot[pathway].add(uniprot)
    for pathway in hallmark_uniprot:
        hallmark_uniprot[pathway] = list(hallmark_uniprot[pathway])
    return hallmark_uniprot

def hallmark_pathway_total_geneList():
    output = set()
    hallmark = hallmark_pathway()
    for pathway in hallmark:
        for gene in hallmark[pathway]:
            output.add(gene)
    return list(output)
    
# CGC genes
def CGC_genes():
    fi_directory = '/home/junghokong/PROJECT/bladder_cancer/code/1_SubtypeSimilarity/data/'
    df = pd.read_csv('%s/Cancer_Genome_Census_allFri Mar 30 05_12_48 2018.tsv' %fi_directory, sep='\t')
    geneList = list(set(df['Gene Symbol']))
    return geneList

# REACTOME genes
def reactome_genes(): # provide in a dictionary
    output = defaultdict(list)
    output_list = []
    fi_directory = '/home/junghokong/PROJECT/bladder_cancer/co_work/propagate_and_NMF_cluster/data'
    f = open('%s/MSigDB_50_hallmark_gene_set/msigdb.v6.1.symbols.gmt.txt' %fi_directory,'r')
    for line in f.xreadlines():
        line = line.strip().split('\t')
        if 'REACTOME' in line[0]:
            reactome = line[0]
            output_list.append(reactome)
            for i in range(2, len(line)):
                gene = line[i]
                output[reactome].append(gene)
    f.close()
    return output

def reactome_genes_uniprot():
    output, reactome = defaultdict(list), reactome_genes()
    gene2uniprot = geneID2uniprot()
    for pathway in reactome:
        for gene in reactome[pathway]:
            if gene in gene2uniprot:
                uniprot = gene2uniprot[gene]
                if not uniprot in output[pathway]:
                    output[pathway].append(uniprot)
    return output

# KEGG genes
def kegg_genes(): # provide in a dictionary
    output = defaultdict(list)
    output_list = []
    fi_directory = '/home/junghokong/PROJECT/bladder_cancer/co_work/propagate_and_NMF_cluster/data'
    f = open('%s/MSigDB_50_hallmark_gene_set/msigdb.v6.1.symbols.gmt.txt' %fi_directory,'r')
    for line in f.xreadlines():
        line = line.strip().split('\t')
        if 'KEGG' in line[0]:
            kegg = line[0]
            output_list.append(kegg)
            for i in range(2, len(line)):
                gene = line[i]
                output[kegg].append(gene)
    f.close()
    return output

def kegg_genes_uniprot():
    output, kegg = defaultdict(list), kegg_genes()
    gene2uniprot = geneID2uniprot()
    for pathway in kegg:
        for gene in kegg[pathway]:
            if gene in gene2uniprot:
                uniprot = gene2uniprot[gene]
                if not uniprot in output[pathway]:
                    output[pathway].append(uniprot)
    return output
   
# ------------------------------------------------------------------------------------------------ 
## gene annotation conversion utilities

def convert_geneList_to_uniprotList( input_geneList ):
    output = []
    for gene in input_geneList:
        if gene in gene2uniprot:
            output.append(gene2uniprot[gene])
    return list(set(output))

def convert_uniprotList_to_geneList( input_uniprotList ):
    output = []
    for uniprot in input_uniprotList:
        if uniprot in uniprot2gene:
            output.append(uniprot2gene[uniprot])
    return list(set(output))
    

## gene annotation    
# ensembl gene annotation
def annotation():
    geneID2ensembl, ensembl2geneID = defaultdict(set), {}
    fi_directory = '/home/junghokong/PROJECT/bladder_cancer/code/1_SubtypeSimilarity/data'
    df = pd.read_csv('%s/2017_07_31_biomart_protein_coding_genes.txt' %fi_directory, sep='\t')
    for i in range(len(df)):
        geneID, ensembl = df['Gene name'][i], df['Gene stable ID'][i]
        geneID2ensembl[ geneID ].add( ensembl )
        ensembl2geneID[ ensembl ] = geneID

    for geneID in geneID2ensembl:
        geneID2ensembl[geneID] = list(geneID2ensembl[geneID])
    return geneID2ensembl, ensembl2geneID

def ensembl2geneID():
    output = {} # { ensembl : geneID }
    fi_directory = '/home/junghokong/PROJECT/bladder_cancer/code/8_hESC/data'
    df = pd.read_csv('%s/2017_07_31_biomart_protein_coding_genes.txt' %fi_directory, sep='\t')
    for i in range(len(df)):
        ensembl, gene = df['Gene stable ID'][i], df['Gene name'][i]
        output[ensembl] = gene
    return output

def geneID2uniprot():
    output = {} # { gene ID : uniprot ID }
    fi_directory = '/home/junghokong/PROJECT/bladder_cancer/code/8_hESC/data'
    df = pd.read_csv('%s/uniprot_homoSapiens_multipleGeneName_20180802.tab' %fi_directory, sep='\t')
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
    df = pd.read_csv('%s/uniprot_homoSapiens_multipleGeneName_20180802.tab' %fi_directory, sep='\t')
    for i in range(len(df)):
        uniprot, geneList = df['Entry'][i], df['Gene names'][i]
        if pd.isnull(geneList) == False:
            geneList = geneList.split()
            gene = geneList[0]
            output[uniprot] = gene
    return output
gene2uniprot, uniprot2gene = geneID2uniprot(), uniprot2geneID()
