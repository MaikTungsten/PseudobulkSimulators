# Pseudobulk simulators for convolution and deconvolution assessment

In this repo, we provide two Python-based pseudobulk simulators.

## Heterogeneous pseudobulk simulator

The heterogeneous pseudobulk simulator was created as an extension of the simulator used in Chen, 2022 (TAPE and Python-implementation of Scaden). Specifically, we introduced multi-threading to accelerate simulation (especially from large single-cell references) and included the option to account for the origin of a single cell (e.g. specific individual or condition), which represents the heterogeneous simulation option. Please note that the scRNA-seq data should be provide as an h5ad file (AnnData object) and must contain at least a column **cell_type** in the observations. Additional **condition** and **individual** columns can be included for heterogeneous pseudobulk simulation. We recommend to supply raw counts rather than CP10K- or CPM-normalized single-cell count data. Further, predefined cell type proportions can be provided for simulation or randomly simulated. 

For exact usage parameters, refer to ``python simulation.py --help``:

```bash
usage: simulation.py [-h] [--sc_path SC_PATH] [--sc_layer SC_LAYER] [--samplenum SAMPLENUM] [--props PROPS] [--target TARGET] [--target_name TARGET_NAME] [--ncells NCELLS] [--rare RARE] [--norm {rank,CPM,raw}] [--filter_genes {all,mRNA,3k}] [--sparse SPARSE] [--threads THREADS] [--outname OUTNAME]

Basic tool for simulating pseudobulk data from single-cell data.

options:
  -h, --help            show this help message and exit
  --sc_path SC_PATH     Input single-cell data path.
  --sc_layer SC_LAYER   Input layer of single-cell data to use for normalization.
  --samplenum SAMPLENUM
                        Number of pseudobulk samples to be simulated.
  --props PROPS         If desired, specify path to a csv file that provides proportions for the desired amount of samples. Must contain all cell types given in sc data.
  --target TARGET       Specify whether simulation should be specific to individuals or condition as source of single cells.
  --target_name TARGET_NAME
                        Specify whether a specific individual or condition should be used for simulation. Only used if target is specified.
  --ncells NCELLS       Number of pseudobulk samples to be simulated.
  --rare RARE           Probability for rare cell types.
  --norm {rank,CPM,raw}
                        Normalization strategy for pseudobulk after aggregation of single cells.
  --filter_genes {all,mRNA,3k}
                        Selection of set of genes from pseudobulks: mRNA, all genes as in single-cels, 3k most variable genes.
  --sparse SPARSE       Probability for sparse cell types (e.g. cell type is not present).
  --threads THREADS     Number of threads to use for computation.
  --outname OUTNAME     Ideally, specifiy a date or tissue ID.
```


## Convolution-like pseudobulk simulator

The convolution-like simulator aims to simulate pseudobulks in a manner reflected in the holistic transcriptome model, where cell type expression and cell type proportions are convoluted to a bulk expression profile. In our simulator, single-cell expression is summarized to cell type-specific expression, for which different modi are available (median, mean or also scValue-based metric). On the basis of the cell type expression matrix, random or user-defined proportions are used to compute pseudobulks.

For exact usage parameters, refer to ``python conv_simulator.py --help``:

```bash
usage: conv_simulator.py [-h] --sc_path SC_PATH [--sc_layer SC_LAYER] [--cg_layer {mean,median,tukeys-biweight,scV-thr,scV-avg}] [--props PROPS] [--outname OUTNAME] [--additional_samples ADDITIONAL_SAMPLES]
                         [--norm {raw,CPM,log,rank,z-score,min-max}] [--from_celltypist] [--gene_filtering {none,mRNA,zeroes} [{none,mRNA,zeroes} ...]]

Tool for simulating pseudobulk data from single-cell data.

options:
  -h, --help            show this help message and exit
  --sc_path SC_PATH     Input single cell matrix path (h5ad format, containing layers of CPM normalized data and log1p-normalized data!). If using celltype annotations, that obs has to be named cell_type.
  --sc_layer SC_LAYER   Input layer of single-cell data to use for calculations
  --cg_layer {mean,median,tukeys-biweight,scV-thr,scV-avg}
                        The layer of the celltype gene matrix that will be used to construct the pseudobulk in the last step.
  --props PROPS         If desired, specify path to a csv file that provides proportions for the desired amount of samples in the format samples (y) * celltypes (x).
  --outname OUTNAME     Ideally, specifiy a date or tissue ID. The resulting files will be saved to: [OUTNAME]_[sc_layer]_[celltype annotation strategy]_[cg_layer]_[normalization
                        strategy]_[proportions/pseudobulks].csv
  --additional_samples ADDITIONAL_SAMPLES
                        Number of additional rows (samples) to simulate and add to the proportions
  --norm {raw,CPM,log,rank,z-score,min-max}
                        How the resulting pseudobulk shall be normalized.
  --from_celltypist     set flag to indicate that cell_types were annotated using CellTypist (not recommended to use) (the celltype obs here are called majority_vote and predicted_labels)
  --gene_filtering {none,mRNA,zeroes} [{none,mRNA,zeroes} ...]
                        Filter genes in resulting pseudobulks for certain gene groups? (e.g. mRNA genes only, an/or remove near-zero expressed genes)
```

## Environment

Refer to the ``env`` subdirectory for conda or pip installation of required packages.

## Acknowledgments

We want to acknowledge instrumental work from other teams, on which our code is based and by which our work is inspired by.

Chen and colleagues created a Python implementation of a standard pseudobulk simulator, which we adapted for our heterogeneous pseudobulk simulator:
[GitHub](https://github.com/poseidonchan/TAPE) and [Publication](https://doi.org/10.1038/s41467-022-34550-9)

The state-of-the-art, standalone simulator SimBu created by Dietrich and colleagues:
[Publication](https://doi.org/10.1093/bioinformatics/btac499)
