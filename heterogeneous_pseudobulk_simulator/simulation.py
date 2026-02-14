import argparse
import anndata
import numpy as np
import scanpy as sc
import pandas as pd
from tqdm import tqdm
import random
from numpy.random import choice
import warnings
from tools import generate_simulated_data, main_gene_selection, normRank, generate_simulated_data_per_target, pseudobulk_norm



# Setup argparser
parser = argparse.ArgumentParser(description='Basic tool for simulating pseudobulk data from single-cell data.')

parser.add_argument('--sc_path', type=str, default='./', help='Input single-cell data path.')
parser.add_argument('--sc_layer', type=str, default='unspecified', help='Input layer of single-cell data to use for normalization.')
parser.add_argument('--samplenum', type=int, default=5000, help='Number of pseudobulk samples to be simulated.')
parser.add_argument('--props', type=str, default=None, help='If desired, specify path to a csv file that provides proportions for the desired amount of samples. Must contain all cell types given in sc data.')
parser.add_argument('--target', type=str, default='no_target', help='Specify whether simulation should be specific to individuals or condition as source of single cells.')
parser.add_argument('--target_name', type=str, default=None, help='Specify whether a specific individual or condition should be used for simulation. Only used if target is specified.')
parser.add_argument('--ncells', type=int, default=1000, help='Number of pseudobulk samples to be simulated.')
parser.add_argument('--rare', type=float, default=0.3, help='Probability for rare cell types.')
parser.add_argument('--norm', type=str, default='CPM', choices=['rank', 'CPM', 'raw'], help='Normalization strategy for pseudobulk after aggregation of single cells.')
parser.add_argument('--filter_genes', type=str, default='mRNA', choices=['all', 'mRNA', '3k'], help='Selection of set of genes from pseudobulks: mRNA, all genes as in single-cels, 3k most variable genes.')
parser.add_argument('--sparse', type=float, default=0.2, help='Probability for sparse cell types (e.g. cell type is not present).')
parser.add_argument('--threads', type=int, default=1, help='Number of threads to use for computation.')
parser.add_argument('--outname', type=str, default='NoDate', help='Ideally, specifiy a date or tissue ID.')

args = parser.parse_args()

# Extract arguments
ncells = args.ncells
samplenum = args.samplenum
sc_path = args.sc_path
outputname = args.outname
sparse = args.sparse
rare = args.rare
propPath = args.props
sc_layer = args.sc_layer
norm = args.norm
filter_genes = args.filter_genes
target = args.target
target_name = args.target_name
threads = args.threads


# Step 1: Load single-cell data, select relevant pre-normalited layer; scData should already be QC'ed beforehand; load proportions if supplied

# Load single-cell data
inData = sc.read_h5ad(sc_path)

if sc_layer not in inData.layers.keys():
    if sc_layer == "unspecified":
        sc_layer = 'X'
    else:
        print('Specified sc_layer ' + sc_layer + ' not in scRNA-seq data object. Using default layer.')
        sc_layer = 'X'
        

# Extract relevant layer from sc data
if sc_layer == 'X':
    print('Using default layer for pseudobulk simulation.')
    scData = inData.to_df() # get X to dataframe
    scData['CellType'] = inData.obs['cell_type'] # annotate a CellType column from cell_type column in obs
elif sc_layer != 'X':
    print('Using layer ' + sc_layer + ' from scRNA-seq data for pseudobulk simulation.')
    scData = inData.to_df(layer = sc_layer) # get X to dataframe
    scData['CellType'] = inData.obs['cell_type'] # annotate a CellType column from cell_type column in obs

if 'individual' in inData.obs.columns.tolist():
    scData['individual'] = inData.obs['individual']
else:
    scData['individual'] = 'unspecified'

if 'condition' in inData.obs.columns.tolist():
    scData['condition'] = inData.obs['condition']
else:
    scData['condition'] = 'unspecified'


# Load proportions, if supplied

if propPath == None:
    props = None
elif propPath.endswith('.csv'):
    props = pd.read_csv(propPath, index_col=0)
    if props.shape[0] != samplenum:
        print('Number of samples in proportions is not matching specified sample number. Sample number is now adjusted to ' + str(props.shape[0]) +'.')
        samplenum = props.shape[0]
    elif props.shape[1] != len(scData['CellType'].value_counts()):
        raise ValueError('Number of cell types in proportions is not matching number of cell types in scRNA-seq data.')
    else:
        print('Props match specified parameters.')
else:
    raise ValueError('Proportions not in csv format. Please supply as csv file of samples x cell types.')


# Step 2: Create simulated data

if target not in ['condition', 'individual']:
    target = None

if target == None:
    pseudobulks = generate_simulated_data(scData,
                                      n = ncells,
                                      samplenum = samplenum,
                                      props=props,
                                      sparse=True,
                                      sparse_prob=sparse,
                                      rare=True,
                                      rare_percentage=rare,
                                      n_jobs=threads)

elif target != None:
    pseudobulks = generate_simulated_data_per_target(scData, target=target, target_name=target_name,
                                                     n = ncells,
                                                     samplenum = samplenum,
                                                     props=props,
                                                     sparse=True,
                                                     sparse_prob=sparse,
                                                     rare=True,
                                                     rare_percentage=rare,
                                                     n_jobs=threads)


# Step 3: Normalize simulated data and filter data in different scenarios

proportionsDF, pseudobulkDF = pseudobulk_norm(pseudobulks, norm, filter_genes)


#Step 4: Export as csv
print('Writing pseudobulks to output.')

proportionsDF.to_csv(outputname+'_pseudobulk_proprotions.csv')
pseudobulkDF.to_csv(outputname+'_pseudobulks.csv')