import argparse
import scanpy as sc
import numpy as np
import pandas as pd
# from scvalue import SCValue
# from scvalue.scvalue import SCValue
import scipy.sparse as sp
import os

import matplotlib.pyplot as plt


parser = argparse.ArgumentParser(description='Tool for simulating pseudobulk data from single-cell data.')

parser.add_argument('--sc_path', type=str, required=True, help='Input single cell matrix path (h5ad format, containing layers of CPM normalized data and log1p-normalized data!). If using celltype annotations, that obs has to be named cell_type.')
parser.add_argument('--sc_layer', type=str, default='X', help='Input layer of single-cell data to use for calculations')
parser.add_argument('--cg_layer', type=str, default='mean', choices=['mean', 'median', 'tukeys-biweight', 'scV-thr', 'scV-avg'], help='The layer of the celltype gene matrix that will be used to construct the pseudobulk in the last step.')
parser.add_argument('--props', type=str, default=None, help='If desired, specify path to a csv file that provides proportions for the desired amount of samples in the format samples (y) * celltypes (x).')
parser.add_argument('--outname', type=str, default='NoDate', help='Ideally, specifiy a date or tissue ID. The resulting files will be saved to: [OUTNAME]_[sc_layer]_[celltype annotation strategy]_[cg_layer]_[normalization strategy]_[proportions/pseudobulks].csv')
parser.add_argument('--additional_samples', type=int, default=1, help='Number of additional rows (samples) to simulate and add to the proportions')
parser.add_argument('--norm', type=str, default="CPM", choices=['raw', 'CPM', 'log', 'rank', 'z-score', 'min-max'], help="How the resulting pseudobulk shall be normalized.")
parser.add_argument('--from_celltypist', action='store_true', help="set flag to indicate that cell_types were annotated using CellTypist (not recommended to use) (the celltype obs here are called majority_vote and predicted_labels)")
parser.add_argument('--gene_filtering', type=str, nargs='+', default=["none"], choices=['none', 'mRNA', 'zeroes'], help='Filter genes in resulting pseudobulks for certain gene groups? (e.g. mRNA genes only, an/or remove near-zero expressed genes)')

args = parser.parse_args()

sc_path = args.sc_path
sc_layer = args.sc_layer
cg_layer = args.cg_layer
propPath = args.props
outname = args.outname
additional_samples = args.additional_samples
norm = args.norm
from_celltypist = args.from_celltypist
gene_filtering = args.gene_filtering
step_counter = 0

# setting additional_samples to 0 if proportions were supplied and additional_samples parameter has not been set by the user
if args.props and args.additional_samples is None:
    additional_samples = 0
else:
    additional_samples = args.additional_samples

# setting the base directory to correctly handle the file called within
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
# handles a future warning that appears when filling nans after mRNA filtering
pd.set_option('future.no_silent_downcasting', True)

# EXPLANATION OF CG_LAYER OPTIONS:
# mean      -> aggregated mean, mean of all cells of a celltype combined into the celltype gene matrix
# scV-thr   -> using SCValues RFValue to assign quality score to each cell, using the threshold to subset the cells, then apply aggregated mean
# scV-avg   -> using SCValues RFValue to assign quality score to each cell, then weight cells by associated quality scores!


"""
WORKFLOW

0. extracting cg matrices from sc input => sc_to_celltypegene
0.1 adding all the returned cg matrices (1 if not from celltpyist, else 2) to a list
=> this list will later be used to continue at step 1

1. reading the celltype subsets
Extracts the relevant layer(s) based on the annotation, typically "cell_type", but the user can also set the flag to indicate if their celltypes were annotated using CellTypists automized model 
If the CellTypists flag is set, then there will be two layers (one for predicted label and one for majority vote)
Then, if supplied, the proportions are loaded, else they are simulated randomly using a dirichlet distribution based on the celltypes present in the single cell input.

2. weighted average 
If different cell types are in the propotions and the single cell data, awaits user input to confirm subsetting the proportions to the intersection of celltypes.
This includes removing celltypes, the user is warned.
The proportions are scaled back to 1 before continuing with the weighted average in any case.

3. normalizing the data
-> CPM, log, rank, min-max, z-score

4. filtering the data
-> mRNA genes
-> near-zero expressed genes

5. saving pseudobulks and proportions to file
"""

# main method
def sc_to_celltypexgene(adata, celltypist, layer):
    """ converts a single cell adata to a celltypexgene adata
    - parameters:
        adata: the adata object that will be converted
        layer: the layer that will be used for the aggregation
    - returns:
        concatenated adata objects: one for each celltype annotation - aggregation combination
    """   

    def aggregate_by_label(annotation_key, layer):
        # annotation key referes to either predicted labels or majority voting
        
        def calculate_tukey_biweight(data, c=6.0):
            median = np.median(data, axis=0)
            mad = np.median(np.abs(data - median), axis=0)
            if np.any(mad == 0):
                mad = mad + 0.001
            scaled_distance = (data - median) / (c * mad)
            weights = np.where(np.abs(scaled_distance) < 1, (1 - scaled_distance**2)**2, 0)  # Weights = 0 if |u| > 1
            weighted_mean = np.sum(weights * data) / (np.sum(weights) + 1)
            return weighted_mean  

        
        # extracting the celltypes
        celltypes = adata.obs[annotation_key].cat.categories

        new_adata = None

        # mean aggregation # TODO: separate the new_adata creation from mean calculation to speed up things
        if cg_layer == "mean":
            print("Mean")
            means = []
            for celltype in celltypes:
                bdata = adata[adata.obs[annotation_key] == celltype]
                if bdata.n_obs == 0:
                    print(f"Warning: No cells found for celltype '{celltype}'")
                    continue  # Skip to next celltype
                if sp.issparse(bdata.layers[layer]):
                    data = bdata.layers[layer].toarray()
                else:
                    data = bdata.layers[layer]
                gene_mean = np.asarray(data.mean(axis=0)).flatten()
                means.append(gene_mean)
            # Stack to create matrix of shape (n_celltypes, n_genes)
            mean_matrix = np.vstack(means)
            new_adata = sc.AnnData(X=mean_matrix)
            new_adata.obs = pd.DataFrame({'CellType': celltypes}, index=celltypes)
            new_adata.var = adata.var.copy()
            new_adata.layers['mean'] = mean_matrix

        # Median
        if cg_layer == "median":
            print("Median")
            medians = []
            for celltype in celltypes:
                bdata = adata[adata.obs[annotation_key] == celltype]
                if bdata.n_obs == 0:
                    print(f"Warning: No cells found for celltype '{celltype}'")
                    continue
                #print(celltype, bdata.shape)
                if sp.issparse(bdata.layers[layer]):
                    data = bdata.layers[layer].toarray()
                else:
                    data = bdata.layers[layer]
                gene_median = np.median(data, axis=0).flatten()  # Take median along cells (rows)
                medians.append(gene_median)
            median_matrix = np.vstack(medians)
            new_adata = sc.AnnData(X=median_matrix)
            new_adata.obs = pd.DataFrame({'CellType': celltypes}, index=celltypes)
            new_adata.var = adata.var.copy()
            new_adata.layers["median"] = median_matrix

        # Tukeys Biweight
        if cg_layer == "tukeys-biweight":
            print("Tukeys Biweight")
            tukeys = []
            for celltype in celltypes:
                bdata = adata[adata.obs[annotation_key] == celltype]
                if bdata.n_obs == 0:
                    print(f"Warning: No cells found for celltype '{celltype}'")
                    continue
                if sp.issparse(bdata.layers[layer]):
                    data = bdata.layers[layer].toarray()
                else:
                    data = bdata.layers[layer]
                #tukey_biweight = calculate_tukey_biweight(data)
                tukey_biweight_per_gene = np.apply_along_axis(calculate_tukey_biweight, axis=0, arr=data)
                tukeys.append(tukey_biweight_per_gene)
            # Stack to create matrix of shape (n_celltypes, n_genes)
            tukey_matrix = np.vstack(tukeys)
            new_adata = sc.AnnData(X=tukey_matrix)
            new_adata.obs = pd.DataFrame({'CellType': celltypes}, index=celltypes)
            new_adata.var = adata.var.copy()
            new_adata.layers['tukeys-biweight'] = tukey_matrix

        # # scV-methods
        # if (cg_layer == "scV-thr") or (cg_layer == "scV-avg"):
        #     print("You have chosen one of the scValue methods for cell-filtering pre-simulation of pseudobulks")

        #     # NEW
        #     # 1. create a copy and convert to dense if sparse
        #     adata_copy = adata.copy()
        #     if sp.issparse(adata_copy.layers[layer]):
        #         data = adata_copy.layers[layer].toarray()
        #     else:
        #         data = adata_copy.layers[layer]
        #     # 2. preprocessing -> with CPM and log1p layers, checked in the beginning
        #     if "highly_variable" not in adata.var.columns:
        #         sc.pp.highly_variable_genes(adata_copy, layer="log1p") # fixed! now uses the correct layer to prevent overflows...
        #     if "pca" not in adata.uns:
        #         sc.pp.pca(adata_copy)
        #     # 3. scvalue class creation (default parameters for now)
        #     scv = SCValue(adata=adata_copy, 
        #         sketch_size=0.1, 
        #         use_rep="X_pca",
        #         cell_type_key="cell_type", 
        #         n_trees=100,
        #         strategy="FB", # OPT: enable user to use the different binning methods? (FB, MTB, TP)
        #         prop_sampling=False,
        #         write_dv=False, # i dont need to save them, as they are kept in the adata obs as "dv"
        #         seed=42)
        #     scv.value() # Ich zitiere: "perform value computation and value-based subsampling" (https://github.com/LHBCB/scvalue/tree/main)
        #     adata_sub = scv.adata_sub # saving the result to adata
        #     #print(adata_sub)
        #     # 4. celltype x gene matrices (für scv thr und scv avg)
        #     if adata_sub.obs["dv"].sum() > 0:  # can only attempt if there are non-zero dv scores...
        #         if (cg_layer == "scV-thr"):
        #             print("using a threshhold to select the cells per celltypes based on the scValue quality scores (dv)")
        #             sc_thr_data = []
        #             kept_celltypes = []
        #             for celltype in celltypes:
        #                 bdata = adata_sub[adata_sub.obs[annotation_key] == celltype]
        #                 celltype_dv_scores = bdata.obs['dv'].values
        #                 dv_threshold = 0.2 # OPT: allow user input to modify the quality threshold
        #                 selected_cells = celltype_dv_scores >= dv_threshold
        #                 selected_expr = bdata[selected_cells, :].layers[layer]
        #                 if selected_expr.shape[0] > 0: # only attempt if any cell of this celltype reached the quality threshold, else would throw an error :o
        #                     sc_thr = np.array(selected_expr.mean(axis=0)).flatten()
        #                     sc_thr_data.append(sc_thr)
        #                     kept_celltypes.append(celltype)
        #             scV_thr_matrix = np.vstack(sc_thr_data)
        #             #new_adata.layers["scV-thr"] = scV_thr_matrix
        #             if sc_thr_data:
        #                 scV_thr_matrix = np.vstack(sc_thr_data)
        #                 new_adata = sc.AnnData(X=scV_thr_matrix)
        #                 new_adata.obs = pd.DataFrame({'CellType': kept_celltypes}, index=kept_celltypes)
        #                 new_adata.var = adata.var.copy()
        #                 new_adata.layers["scV-thr"] = scV_thr_matrix
        #             else:
        #                 raise ValueError("No celltypes passed the scValue threshold. Cannot perform scValue-based thresholding!")
        #         if (cg_layer == "scV-avg"):
        #             print("using the scValue quality scores of cells (dv) for a weighted aggregation to celltype x gene matrix")
        #             sc_avg_data = []
        #             kept_celltypes = []
        #             for celltype in celltypes:
        #                 bdata = adata_sub[adata_sub.obs[annotation_key] == celltype]
        #                 if bdata.n_obs == 0:
        #                     print(f"No cells were found for celltype {celltype}. Excluding {celltype} from the scValue quality weighted aggregation.")
        #                     continue
        #                 celltype_dv_scores = bdata.obs['dv'].values
        #                 if np.sum(celltype_dv_scores) == 0:
        #                     print(f"The quality scores of the cells of type '{celltype}' are all zero. Excluding {celltype} from the aggregation")
        #                     continue
        #                 selected_expr = bdata.layers[layer]
        #                 if sp.issparse(selected_expr):
        #                     selected_expr = selected_expr.toarray()
        #                 try:
        #                     sc_avg = np.average(selected_expr, axis=0, weights=celltype_dv_scores)
        #                     sc_avg_data.append(sc_avg)
        #                     kept_celltypes.append(celltype)
        #                 except ZeroDivisionError:
        #                     print(f"Encountered zero division error while attempting weight normalization using cell quality scores of {celltype}")
        #                     continue
        #             if sc_avg_data:
        #                 scV_avg_matrix = np.vstack(sc_avg_data)
        #                 new_adata = sc.AnnData(X=scV_avg_matrix)
        #                 new_adata.obs = pd.DataFrame({'CellType': kept_celltypes}, index=kept_celltypes)
        #                 new_adata.var = adata.var.copy()
        #                 new_adata.layers["scV-avg"] = scV_avg_matrix
        #             else:
        #                 raise ValueError("All celltypes had non-valid cell quality scores > 0. Cannot create scV-avg layer.") # this is very very unlikely
        #     else:
        #         print("None of the cells have a quality score greater than 0, using scValue scores for aggregation not possible.")
        if new_adata is None:
            raise RuntimeError("new_adata was not created. The aggregation method you chose as cg_layer parameter might not exist?")
        return new_adata    
    

    
    ## accounting for different cell type annotations of CellTypist
    # predicted_labels vs majority_vote
    if celltypist == True:
        #print(celltypist)
        #print("Using celltypists cell type annotation")
        Celltypist_pl = aggregate_by_label("predicted_labels", layer)
        Celltypist_mv = aggregate_by_label("majority_voting", layer)
        return Celltypist_pl, Celltypist_mv
    else: # other cell type annotation, e.g. pregiven in data:
        # need to ensure that this was renamed to cell_type !
        #print("using cell_type annotation")
        initial_celltypes = aggregate_by_label('cell_type', layer) # CellType / initial_cell_type / cell_type
        return initial_celltypes


### STEP 0 ###
print(f"\n\n=== STEP {step_counter}: reading sc matrix, converting to cg matric ===")
step_counter += 1

# handling user-error: neither proportions nor additional_samples
if propPath == None:
    if additional_samples == 0:
        raise ValueError('Neither proportions were provided, nor additional samples shall be randomly simulated. Please check your input parameters.')

scData = sc.read_h5ad(sc_path)
total_sc_cells = scData.n_obs

# remove the None- genes (if present)
try:
    print("None gene symbols present :(")
    none_genes = scData.var['gene_symbols'].astype(str).str.startswith('None')
    genes_to_keep = ~(none_genes)
    scData = scData[:, genes_to_keep]
except:
    print()

# TODO: handling user errors: CellType not in obs
if "cell_type" not in scData.obs:
    raise ValueError("single cell input needs cell type annotation as the obs: cell_type")

# check for CPM and log1p layer being present, else abort
if "CPM" not in scData.layers:
    raise ValueError("CPM normalized data needs to be stored in CPM-layer of the input single cell data!")
if "log1p" not in scData.layers:
    raise ValueError("log1p normalized data needs to be stored in log1p-layer of the input single cell data!")

# handle missing sc_layers
if sc_layer not in scData.layers.keys(): # in this case, use default adata.X!
    sc_layer = 'X' # TODO: when trying to call .layers on scData, check if sc_layer = "X" -> then use scData.X
    print(f"Specified sc_layer '{args.sc_layer}' not in scRNA-seq data. Using default layer 'X'.")
    scData = scData.copy()
    scData.layers["X"] = scData.X.copy()

# constructing the cg matrices from sc input (with the specified sc input layer (sc_layer))
cg_matrices = []
annotation_strategy = []
if from_celltypist == True:
    print("Using Celltypist-annotated celltypes of input sc_data.")
    pl, mv = sc_to_celltypexgene(scData, celltypist=True, layer=sc_layer)
    cg_matrices.append(pl)
    annotation_strategy.append("Celltypist_pl")
    cg_matrices.append(mv)
    annotation_strategy.append("Celltypist_mv")
else:
    print("Using annotated inital celltypes of input sc_data.")
    # not from celltypist, only returning one cg matrix
    ct = sc_to_celltypexgene(scData, celltypist = False, layer=sc_layer)
    cg_matrices.append(ct)
    annotation_strategy.append("celltype") # using the celltypes saved under cell_type obs


### STEP 1 ###
print(f"\n\n=== STEP {step_counter}: reading cg matrix ===")
step_counter += 1

for i in range(len(cg_matrices)): 
    # this loop is accounting for when we have a celltypist annotation and two instead of only one aggregated celltype x gene matrix!
    # using cell_type => only one cg_matrix
    # using CellTypists celltypes (not recommended bc quality issues) => two cg matrices: predicted_labels (pl) and majority_vote (mv)
    print(annotation_strategy[i]) # <- and this is the cg matrix used, printed just for the convenience of allowing the user to check
    inData = cg_matrices[i]

    #print(inData)

    if cg_layer not in inData.layers.keys():
        cg_layer = 'mean'
        print('Specified cg_layer ' + cg_layer + ' not in cell type gene matrix. Using default layer mean.')

    cgData = inData.to_df(layer = cg_layer)

    #extracting the present cell types from input sc-h5ad
    present_celltypes = inData.obs['CellType'].tolist()

    print("Collecting proportions...")
    # loading proportions
    #props_larger = False
    if propPath == None:
        props = None
        # if we have no props as input, we have to randomly simulate based on the present celltypes -> later 
    elif propPath.endswith('.csv'):
        props_longer = False
        props = pd.read_csv(propPath, index_col=0)
        # stripping and converting to lowerspace to make celltypes comparable
        celltypes_in_props = [celltype.strip().lower() for celltype in props.columns]
        celltypes_in_sc = [celltype.strip().lower() for celltype in cgData.index]
        # catch common genes
        common_celltypes = list(set(celltypes_in_props).intersection(celltypes_in_sc))
        if not common_celltypes:
            raise ValueError("No overlapping celltypes between proportions and single cell input.")
        only_in_props = list(set(celltypes_in_props) - set(celltypes_in_sc))
        only_in_sc = list(set(celltypes_in_sc) - set(celltypes_in_props))
        if len(celltypes_in_sc) != len(celltypes_in_props) or only_in_props or only_in_sc:
            print("\n- ATTENTION: MISMATCH IN CELLTYPES ----------------------------------------------------")
            print("-")
            print("- There are celltypes that are either only in the proportions or the single cell input:")
            print("- only in proportions:\t", only_in_props)
            print("- only in sc input:\t", only_in_sc)
            print("-")
            print("-----------------------------------------------------------------------------------------")
            valid_celltypes = [col for col in props.columns if col.strip().lower() in common_celltypes]
            response = input(f"\nThe celltypes given in the proportions and the single cell input do not match.\nDo you want to continue and remove any unmatched celltypes?\nPLEASE NOTE: This will keep only the following celltypes: {valid_celltypes} (yes/no): \n>> ").strip().lower()
            if response == 'yes':
                print("removing unmatched celltypes...")
                valid_celltypes = [col for col in props.columns if col.strip().lower() in common_celltypes]
                print("\nKeeping only the following celltypes:")
                print(valid_celltypes)
                # filtering the cgData
                rows_to_keep = [row for row in cgData.index if row in valid_celltypes]
                cgData = cgData.loc[rows_to_keep]
                # filtering the sc data
                cells_to_keep = inData.obs[inData.obs['CellType'].isin(valid_celltypes)].index
                inData = inData[cells_to_keep].copy()
                # resetting the present celltypes
                present_celltypes = valid_celltypes
                # filter and rescale the proportions
                props = props[valid_celltypes]
                props = props.div(props.sum(axis=1), axis=0)
            else:
                raise ValueError("Operation aborted by user due to mismatched celltypes.")
    else:
        raise ValueError('Proportions not in csv format. Please supply as csv file of samples x cell types.')

    ## creating additional random sample proportions if specified by the user via "additional_samples"
    if additional_samples > 0:
        simulated_props = []
        rng = np.random.default_rng(seed=42) # seed to increase reproducibility :)
        for _ in range(additional_samples):
            skewness = 1 
            # OPTIONAL: The skewness parameter could be optimized dependent on what tissue is fed into the pipeline.
            # A low skewness/alpha leads to skewed distributions, while higher values smooth things out.
            # I am using skewenss=1 as this is also used by default in the simulator built into the Scaden deconvolution tool, as this enhances comparability.
            random_props = rng.dirichlet(np.full(len(present_celltypes), skewness), size=1)
            simulated_props.append(random_props.flatten())
        simulated_props_df = pd.DataFrame(simulated_props, columns=present_celltypes)
        if props is not None:
            props = pd.concat([props, simulated_props_df], ignore_index=True)
        else:
            props = simulated_props_df    
    props = props.fillna(0)

    ### STEP 2 ###
    print(f"\n\n=== STEP 2: creating simulated data (weighted averages) ===")
    step_counter += 1

    cgData_indexed = cgData.copy()
    cgData_indexed['CellType'] = inData.obs['CellType'].values
    cgData_indexed = cgData_indexed.set_index('CellType')
    cgData_ordered = cgData_indexed.loc[props.columns]
    cgData_ordered = cgData_ordered.fillna(0)

    # simulating pseudobulks (weighted average)
    pseudobulks = []
    for idx, row in props.iterrows():
        weights = row.values
        pseudo = np.sum(cgData_ordered.values * weights[:, np.newaxis], axis=0)
        pseudo_scaled = pseudo * total_sc_cells
        pseudobulks.append(pseudo_scaled)

    # creating the pseudobulks dataframe with genes as columns
    pseudobulk_df = pd.DataFrame(pseudobulks, index=props.index, columns=cgData_ordered.columns)

    print(f"Simulated pseudobulk matrix of shape {pseudobulk_df.shape}:")
    print(pseudobulk_df.head(3))

    ### STEP 3 ###
    if "zeroes" in gene_filtering or "mRNA" in gene_filtering:
        print(f"\n\n=== STEP {step_counter}: gene filtering ===")
        step_counter += 1
        # Filtering out of near-zero expressed genes handles the observed sparsity of single cell data
        # The threshold is set to 0.05 as default, meaning that genes that have a expression of 0 in less than 5% of cells will be excluded
        # Because bulks cannot be filtered in a similar and comparable way, bulks that are compared to pseudobulks that were filtered using this method shall be filtered to include only non-zero expressed genes

    final_df = pseudobulk_df.copy()

    # 1. filter for mRNA genes (if gene_filtering = mRNA)
    if "mRNA" in gene_filtering:

        # rebuild of the mRNA filtering:
        # 1. create new df with the mRNA indices
        # 2. if the index is a prefix of the input indices, take the associated value and add it to the new df
        # 3. fill NaN (.fillna(0))
        # 4. save the temporary df back to final_df

        print("Filtering for mRNA genes only")
        gene_list_path = os.path.join(BASE_DIR, "mRNA_gene_list.tsv")
        gene_list_df = pd.read_csv(gene_list_path, header=0, delimiter="\t")
        gene_list = list(gene_list_df['gene_name'])
        new_df = pd.DataFrame(index=final_df.T.columns, columns=gene_list)
        
        try:
            # extract the indices that are both in the mRNA annotation list and the input = filtered_indices
            #filtered_indices = [idx for idx in final_df.T.index if any(idx.startswith(prefix) for prefix in gene_list)] # edge case: SOX21 in mRNA annotation and SOX21-AS1 in the sandbox dataset (G01, pancreatic islets)
            final_df = final_df.T # it is easier to handle the filtering in the transposed orientation

            for prefix in gene_list:
                matching_cols = [col for col in final_df.T.columns if col.startswith(prefix)]
                if matching_cols:
                    new_df[prefix] = final_df.T[matching_cols].mean(axis=1)
                else:
                    # No match: keep NaNs for now, will fill with 0 later
                    pass
            final_df = new_df.fillna(0.0)

            print(f"After mRNA filtering, df is of shape: {final_df.shape}:")
            print(final_df.head(3))
        except KeyError:
            gene_filtering = "none"
            print("Filtering failed: Please check the gene symbol annotation in file provided as sc_path")
            print("Continuing with pseudobulk that was NOT filtered for mRNA genes...")
    # 2. filter out non-expressed genes to adress sparsity
    if "zeroes" in gene_filtering:
        print("Filtering out near-zero expressed genes") 
        # because different pseudobulks might have different 0-expressed genes, just replace these occurences with NaN
        final_df.replace(0, np.nan, inplace=True)
        mean_nan_count = final_df.isna().sum(axis=1).mean()
        num_genes = final_df.shape[1]
        nan_percentage = (mean_nan_count / num_genes) * 100
        print(f"The average number of non-expressed genes in all pseudobulks is {mean_nan_count:.2f} ({nan_percentage:.2f}%).")
        print("Non-expressed genes (0.0) were filled with NaNs.")
        print(f"After zeroes filtering, df is of shape: {final_df.shape}:")
        print(final_df.head(3))


    ### STEP 4 ###
    print(f"\n\n=== STEP {step_counter}: normalization of pseudobulks ===")
    step_counter += 1

    if norm == "raw":
        print("Returning raw aggregated pseudobulks.")
        final_df = final_df.copy()
    elif norm == 'CPM':
        print('Scaling pseudobulks to CPM.')
        final_df = final_df.div(final_df.sum(axis=1), axis=0) * 1e6
    elif norm == "log": 
        print("Scaling to log-norm.")
        cpm_df = final_df.div(final_df.sum(axis=1), axis=0) * 1e4  
        final_df = np.log1p(cpm_df)
    elif norm == "rank": 
        print("Normalizing to rank-order (min).")
        df = final_df.copy()
        final_df = df.rank(axis=1, method="min", na_option="keep", ascending=True)
        final_df = final_df.div(final_df.max(axis=1), axis=0)
        # min rank normalization plus min-max scaling to make result comparable to Scadens' simulator
    elif norm == "z-score":
        print("performing z-score normalization") # sample-wise
        final_df = final_df.sub(final_df.mean(axis=1), axis=0) # x - mü
        final_df = final_df.div(final_df.std(axis=1), axis=0) # divided by standard deviation
    elif norm == "min-max":
        print("performing min-max scaling") # sample-wise
        min_vals = final_df.min(axis=1)
        max_vals = final_df.max(axis=1)
        final_df = final_df.sub(min_vals, axis=0)
        final_df = final_df.div(max_vals - min_vals, axis=0)
        final_df = final_df.fillna(0.0) # in case the minimal and maximal value are identical, which is very very unlikely for transcriptomics data
    else:
        print("This is never reached :)")


    print(f"Final pseudobulk matrix of shape {final_df.shape}:")
    print(final_df.head(3))


    ### STEP 5 ###
    # saving the generated pseudobulk to csv
    print(f"\n\n=== STEP {step_counter}: saving to {outname} ===")
    filter_tag = "none" if gene_filtering == ["none"] or gene_filtering == "none" else "_".join(sorted(gene_filtering))


    props.to_csv(f"{outname}_pseudobulk_proportions.csv")
    print("Saved proportions to: " + f"{outname}_pseudobulk_proportions.csv")
    final_df.to_csv(f"{outname}_pseudobulks.csv")
    print("Saved pseudobulks to: " + f"{outname}_pseudobulks.csv")

# final outname should follow this:
# [[ID]_[sample/patient]]_[sc_layer]_[celltype annotation strategy]_[cg_layer]

# [[ID]_[sample/patient]] : this has to be specified as the outname :)
# celltype annotation strategy: (Celltypist_pl, Celltypist_mv, initial_celltypes)
# cg_layer : (mean, mas, scV-thr, scV-avg)
