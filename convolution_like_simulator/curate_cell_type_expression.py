import pandas as pd
import numpy as np
import scanpy as sc
from anndata import AnnData
import argparse


# Argparser
parser = argparse.ArgumentParser(description='Convert cell type expression in csv format to h5ad with CPM and log1p layers.')

parser.add_argument('--path', type=str, default='./', help='Input cell type expression data path in csv format. Contains individual column, which is written into observations.')

args = parser.parse_args()

# Get dataframe from csv
df = pd.read_csv(args.path, index_col=0)

# Extract patient info
obs = pd.DataFrame({'individual' : list(df['individual']), 'cell_type' : df.index.tolist()})
print(f"Shape observations: {obs.shape}")
obs.index = [f"{obs.iloc[i,0]}_{obs.iloc[i,1]}" for i in range(df.shape[0])]

df.index = obs.index.tolist()
# Drop patient info from gene expression
df.drop(['individual'], axis=1, inplace=True)

# Create anndata object
adata = AnnData(X = df)
adata.obs = obs

# Create CPM and log1p layers
adata.layers["CPM"] = adata.X
sc.pp.normalize_total(adata, layer='CPM', target_sum=1e6)
adata.layers["log1p"] = np.log(adata.layers["CPM"]+1)

# Write to h5ad file
adata.write_h5ad('Data_for_bulk_mixing.h5ad')