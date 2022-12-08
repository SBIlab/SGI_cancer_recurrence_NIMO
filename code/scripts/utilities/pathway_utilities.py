################## provide pathway gene list ##################
import pandas as pd
import numpy as np
import scipy.stats as stat
from collections import defaultdict
import os, time


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
    fi_directory = '../../data'
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
    fi_directory = '../../data'
    df = pd.read_csv('%s/2017_07_31_biomart_protein_coding_genes.txt' %fi_directory, sep='\t')
    for i in range(len(df)):
        ensembl, gene = df['Gene stable ID'][i], df['Gene name'][i]
        output[ensembl] = gene
    return output

def geneID2uniprot():
    output = {} # { gene ID : uniprot ID }
    fi_directory = '../../data'
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
    fi_directory = '../../data'
    df = pd.read_csv('%s/uniprot_homoSapiens_multipleGeneName_20180802.tab' %fi_directory, sep='\t')
    for i in range(len(df)):
        uniprot, geneList = df['Entry'][i], df['Gene names'][i]
        if pd.isnull(geneList) == False:
            geneList = geneList.split()
            gene = geneList[0]
            output[uniprot] = gene
    return output
gene2uniprot, uniprot2gene = geneID2uniprot(), uniprot2geneID()
