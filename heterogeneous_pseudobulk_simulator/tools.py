import pathlib
import anndata
import numpy as np
import scanpy as sc
import pandas as pd
from tqdm import tqdm
import random
from numpy.random import choice
import warnings
from scipy.stats import rankdata
from joblib import Parallel, delayed

# Function for simulation of a single pseudobulk

def simulate_single_pseudobulk(sample_prop, scData, celltype_groups, allcellname):

    """
    Function for simulating a single pseudobulk given:
    sample_prop: proportions for a single sample
    scData: single-cell data in pd.DataFrame
    celltype_groups
    allcellname
    """
    
    single_pseudobulk = np.zeros(scData.shape[1])
    for j, cellname in enumerate(allcellname):
        if sample_prop[j] > 0:
            select_index = choice(celltype_groups[cellname], size=int(sample_prop[j]), replace=True)
            single_pseudobulk += scData.loc[select_index, :].sum(axis=0).values
    return single_pseudobulk

# Simulator for proportions

def simulate_proportions(props, random_state, d_prior, num_celltype, samplenum, sparse, sparse_prob, rare, rare_percentage):

    """
    Simulates proportions (or returns already given proportions)
    Parameters:
    num_celltype: number of cell types in reference
    samplenum: number of pseudobulks to simulate
    sparsity parameters: sparse, sparse_prob
    rare cell types: rare, rare_percentage
    random_state: can be set for reproducibility
    """

    if not isinstance(props, pd.DataFrame):
        # generate random cell type proportions

        if random_state is not None and isinstance(random_state, int):
            print('Random state specified. This will improve the reproducibility.')

        if d_prior is None:
            print('Generating cell fractions using Dirichlet distribution without prior info (actually random).')
            if isinstance(random_state, int):
                np.random.seed(random_state)
            prop = np.random.dirichlet(np.ones(num_celltype), samplenum) # randomly alter proportion of 1 per cell type via dirichlet distribution for samplenum, e.g. 5000 times
            print('Random cell fractions are generated.')
        elif d_prior is not None:
            print('Using prior info to generate cell fractions in Dirichlet distribution')
            assert len(d_prior) == num_celltype, 'dirichlet prior is a vector, its length should equals ' \
                                                'to the number of cell types'
            if isinstance(random_state, int):
                np.random.seed(random_state)
            prop = np.random.dirichlet(d_prior, samplenum) # d_prior would set an initial cell type distribution to randomly and slightly alter for samplenum
            print('Dirichlet cell fractions are generated.')

        prop = prop / np.sum(prop, axis=1).reshape(-1, 1) # scale to have percentage proportions for each cell type summing up to 100 %; each row is one simulated sample and its proportions
        
        
        # sparse cell fractions
        if sparse:
            print("You set sparse as True, some cell's fraction will be zero, the probability is", sparse_prob)
            ## Only partial simulated data is composed of sparse celltype distribution
            for i in range(int(prop.shape[0] * sparse_prob)):
                indices = np.random.choice(np.arange(prop.shape[1]), replace=False, size=int(prop.shape[1] * sparse_prob)) # chose a random cell type index from 0 to no. of cell types to drop out per sample to be simulated
                prop[i, indices] = 0

            prop = prop / np.sum(prop, axis=1).reshape(-1, 1) # resize proportions

        if rare:
            print(
                'Selected rare, thus some cell type fractions are very small (<3%), '
                'this celltype is randomly chosen by percentage you set before.')
            ## choose celltype to be rare
            np.random.seed(0)
            indices = np.random.choice(np.arange(prop.shape[1]), replace=False, size=int(prop.shape[1] * rare_percentage))
            prop = prop / np.sum(prop, axis=1).reshape(-1, 1)

            for i in range(int(0.5 * prop.shape[0]) + int(int(rare_percentage * 0.5 * prop.shape[0]))):
                prop[i, indices] = np.random.uniform(0, 0.03, len(indices)) # rare between 0 % and 3 %
                buf = prop[i, indices].copy()
                prop[i, indices] = 0
                prop[i] = (1 - np.sum(buf)) * prop[i] / np.sum(prop[i])
                prop[i, indices] = buf
    
    elif isinstance(props, pd.DataFrame):
        prop = np.array(props)
    
    return prop


# Simulate data

def generate_simulated_data(scData, d_prior=None,
                            n=1000, samplenum=5000, props = None,
                            random_state=None, sparse=True, sparse_prob=0.2,
                            rare=False, rare_percentage=0.3, n_jobs = 8):
    
    """
    Pseudobulk simulation with respect only to "cell_type" column in scData
    """
    
    # scData is a pd.DataFrame (a cells x genes matrix)
    
    genes = list(scData.iloc[:,0:scData.shape[1]-3].columns) # save gene names for later

    # Get number and groups of cell types present
      
    num_celltype = len(scData['CellType'].value_counts())
    celltype_groups = scData.groupby('CellType').groups # dictionary of cell types
    scData = scData.drop(columns=['CellType', 'individual', 'condition'])
    scData = scData.astype(np.float32)
    
    # Simulate proportions
    prop = simulate_proportions(props, random_state, d_prior, num_celltype, samplenum,
                                sparse, sparse_prob, rare, rare_percentage)

    # precise number for each celltype
    cell_num = np.floor(n * prop) # from proportions get the number of cells per cell type to be sampled given n (e.g. n=1000 cells)
    
    # precise proportion based on cell_num
    prop = cell_num / np.sum(cell_num, axis=1).reshape(-1, 1) # then, obtain accurate proportions that fit the given number n

    # start sampling
    allcellname = celltype_groups.keys() # cell type names
    print('Sampling cells to compose pseudo-bulk data...')
    
    results = Parallel(n_jobs=n_jobs)(
        delayed(simulate_single_pseudobulk)(
            cell_num[i],
            scData,
            celltype_groups,
            allcellname
        )
        for i in tqdm(range(cell_num.shape[0]))
    )

    sample = np.stack(results)

    # Create final dataframes
    if not isinstance(props, pd.DataFrame):
        sampleDF = pd.DataFrame(sample, index=['Sample_'+str(i) for i in range(1,prop.shape[0]+1)], columns = genes[0:sample.shape[1]])
        prop = pd.DataFrame(prop, index = ['Sample_'+str(i) for i in range(1,prop.shape[0]+1)], columns=celltype_groups.keys())

    elif isinstance(props, pd.DataFrame):
        sampleDF = pd.DataFrame(sample, index = props.index.tolist(), columns = genes[0:sample.shape[1]])
        prop = pd.DataFrame(prop, index = props.index.tolist(), columns = props.columns.tolist())
    
    simudata = anndata.AnnData(X=sampleDF,
                               obs=prop,
                               #var=genes[0:sample.shape[1]]
                               )

    print('Sampling is done.')
    
    return simudata # anndata object containing ground-truth proportions in observations and having geneIDs as numbers for hidden layers


# Simulate data with respect to specific conditions and or individuals/patients

def generate_simulated_data_per_target(scData,
                            target = None, target_name = None, d_prior=None,
                            n=1000, samplenum=5000, props = None,
                            random_state=None, sparse=True, sparse_prob=0.2,
                            rare=False, rare_percentage=0.3, n_jobs = 8):
    """
    Each patient/individual should have the same cell types.

    Pseudobulk simulation with respect to "cell_type" and "target", which can be individual or condition column in scData

    Target-based simulation can also be directed to a single target only
    """

    if target == 'individual':

        # scData is a pd.DataFrame (cells x genes matrix)
        genes = list(scData.iloc[:,0:scData.shape[1]-3].columns) # save gene names for later

        # Get number and groups of cell types present
        num_celltype = len(scData['CellType'].value_counts())
        
        if target_name == None:
            print('Simulating cell type proportions and pseudobulks for all individuals.')
            
            # Define list to store proportions and pseudobulks per individual
            prop_list = []
            pseudobulk_list = []

            all_celltypes = sorted(scData['CellType'].unique().tolist()) # all cell types from scData, might vary per individual

            for patient in list(set(scData['individual'].tolist())):

                print(f"Simulating pseudobulks for individual {patient}...")

                # Subset scData to individual selected
                df_work = scData.loc[scData['individual'] == patient]
                df_work['CellType'] = df_work['CellType'].cat.remove_unused_categories()
                celltype_groups = df_work.groupby('CellType').groups # dictionary of cell types
                df_work = df_work.drop(columns=['CellType', 'individual', 'condition'])

                # Sample number per individual
                samplenum_per_individual = int(np.floor(samplenum/len(list(set(scData['individual'].tolist())))))

                present_ct = list(celltype_groups.keys())
                missing_cell_types = list(set(all_celltypes) - set(present_ct))

                # number of CTs to simulate = number actually present in this patient
                num_celltype_patient = len(present_ct)

                if len(missing_cell_types) > 0:
                    print(f"{len(missing_cell_types)} cell types are missing for individual {patient}.",
                    f"Fractions of missing cell types will be automatically set to 0.")
                else:
                    missing_cell_types = None
                
                # Simulate proportions
                prop = simulate_proportions(props, random_state, d_prior, num_celltype_patient, samplenum_per_individual,
                                            sparse, sparse_prob, rare, rare_percentage)
                
                # Fill proportions with zeros for missing cell types at appropriate position
                if missing_cell_types != None:
                    propDF = pd.DataFrame(0, index=np.arange(prop.shape[0]), columns=list(all_celltypes))
                    propDF[present_ct] = prop
                    prop = propDF.to_numpy()

                # precise number for each celltype
                n_adapt = int(np.floor(df_work.shape[0]/2)) # maximum number of cells, if cells per patient too few
                if n_adapt > n:
                    n_adapt = n
                    print(f"Using {str(n_adapt)} cells for sampling.")
                else: 
                    print(f"Using {str(n_adapt)} cells for sampling.")
                cell_num = np.floor(n_adapt * prop) # from proportions get the number of cells per cell type to be sampled given n (e.g. n=1000 cells)
 
                # precise proportion based on cell_num
                prop = cell_num / np.sum(cell_num, axis=1).reshape(-1, 1) # then, obtain accurate proportions that fit the given number n

                # start sampling
                allcellname = all_celltypes # cell type names
                print('Sampling cells to compose pseudo-bulk data...')
  
                results = Parallel(n_jobs=n_jobs)(
                    delayed(simulate_single_pseudobulk)(
                        cell_num[i],
                        df_work,
                        celltype_groups,
                        allcellname
                    )
                    for i in tqdm(range(cell_num.shape[0]))
                )

                sample = np.stack(results)
     
                sampleDF = pd.DataFrame(sample, index=[f"{patient}_sample_{i}" for i in range(1,prop.shape[0]+1)], columns = genes[0:sample.shape[1]])
                prop = pd.DataFrame(prop, index = [f"{patient}_sample_{i}" for i in range(1,prop.shape[0]+1)], columns=all_celltypes)
                if missing_cell_types == None:
                    print(f"All cell types present for individual")
                else:
                    prop[missing_cell_types] = 0

                # Append props and bulks from current iteration to list
                prop_list.append(prop)
                pseudobulk_list.append(sampleDF)

            # Concatenate all simulated data in a single df
            pseudobulks = pd.concat(pseudobulk_list, axis=0)
            proportions = pd.concat(prop_list, axis=0)

        elif target_name != None:
            print('Simulating cell type proportions and pseudobulks for specified individual.')
            
            df_work = scData.loc[scData['individual'] == target_name]
            df_work['CellType'] = df_work['CellType'].cat.remove_unused_categories()
            num_celltype = len(df_work['CellType'].value_counts())
            celltype_groups = df_work.groupby('CellType').groups # dictionary of cell types
            df_work = df_work.drop(columns=['CellType', 'individual', 'condition'])
                
            # Simulate proportions
            prop = simulate_proportions(props, random_state, d_prior, num_celltype, samplenum,
                                        sparse, sparse_prob, rare, rare_percentage)
                
            # precise number for each celltype
            n_adapt = int(np.floor(df_work.shape[0]/2)) # maximum number of cells, if cells per patient too few
            if n_adapt > n:
                n_adapt = n
                print(f"Using {str(n_adapt)} cells for sampling.")
            else: 
                print(f"Using {str(n_adapt)} cells for sampling.")
            cell_num = np.floor(n_adapt * prop) # from proportions get the number of cells per cell type to be sampled given n (e.g. n=1000 cells)

            # precise proportion based on cell_num
            prop = cell_num / np.sum(cell_num, axis=1).reshape(-1, 1) # then, obtain accurate proportions that fit the given number n

            # start sampling
            allcellname = celltype_groups.keys() # cell type names
            print(f"Allcellname: {allcellname}")
            print('Sampling cells to compose pseudo-bulk data...')

            results = Parallel(n_jobs=n_jobs)(
                delayed(simulate_single_pseudobulk)(
                    cell_num[i],
                    df_work,
                    celltype_groups,
                    allcellname
                )
                for i in tqdm(range(cell_num.shape[0]))
            )

            sample = np.stack(results)

            pseudobulks = pd.DataFrame(sample, index=[f"{target_name}_sample_{i}" for i in range(1,prop.shape[0]+1)], columns = genes[0:sample.shape[1]])
            proportions = pd.DataFrame(prop, index=[f"{target_name}_sample_{i}" for i in range(1,prop.shape[0]+1)], columns=celltype_groups.keys())

        # Finalize simulated data in AnnData object
        simudata = anndata.AnnData(X=pseudobulks, obs=proportions)
    
    if target == 'condition':

        # scData is a pd.DataFrame (cells x genes matrix)
        genes = list(scData.iloc[:,0:scData.shape[1]-3].columns) # save gene names for later

        # Get number and groups of cell types present
        num_celltype = len(scData['CellType'].value_counts())

        if target_name == None:
            print('Simulating cell type proportions and pseudobulks for all conditions.' \
            'Assuming that all cell types are present in all conditions.')
            
            # Define list to store proportions and pseudobulks per individual
            prop_list = []
            pseudobulk_list = []

            all_celltypes = scData.groupby('CellType').groups.keys() # all cell types from scData, might vary per condition

            for condition in list(set(scData['condition'].tolist())):
                
                print(f"Simulating pseudobulks for condition {condition}...")

                # Subset scData to individual selected
                df_work = scData.loc[scData['condition'] == condition]
                df_work['CellType'] = df_work['CellType'].cat.remove_unused_categories()
                celltype_groups = df_work.groupby('CellType').groups # dictionary of cell types
                df_work = df_work.drop(columns=['CellType', 'individual', 'condition'])

                # Sample number per individual
                samplenum_per_condition = int(np.floor(samplenum/len(list(set(scData['condition'].tolist())))))

                missing = list(set(all_celltypes) - set(celltype_groups.keys()))

                if len(missing) > 0:
                    print(f"{len(missing)} cell types are missing for condition {condition}.",
                    f"Fractions of missing cell types will be automatically set to 0.")
                    num_celltype_condition = len(celltype_groups.keys())
                    missing_cell_types = missing
                else:
                    num_celltype_condition = num_celltype
                    missing_cell_types = None
                
                # Simulate proportions
                prop = simulate_proportions(props, random_state, d_prior, num_celltype_condition, samplenum_per_condition,
                                            sparse, sparse_prob, rare, rare_percentage)
                
                # Fill proportions with zeros for missing cell types at appropriate position
                if missing_cell_types != None:
                    propDF = pd.DataFrame(0, index=np.arange(prop.shape[0]), columns=list(all_celltypes))
                    present_ct = list(celltype_groups.keys())
                    propDF[present_ct] = prop
                    prop = propDF.to_numpy()

                # precise number for each celltype
                n_adapt = int(np.floor(df_work.shape[0]/2)) # maximum number of cells, if cells per patient too few
                if n_adapt > n:
                    n_adapt = n
                    print(f"Using {str(n_adapt)} cells for sampling.")
                else:
                    print(f"Using {str(n_adapt)} cells for sampling.")
                cell_num = np.floor(n_adapt * prop) # from proportions get the number of cells per cell type to be sampled given n (e.g. n=1000 cells)

                # precise proportion based on cell_num
                prop = cell_num / np.sum(cell_num, axis=1).reshape(-1, 1) # then, obtain accurate proportions that fit the given number n

                # start sampling
                allcellname = all_celltypes # cell type names
                print('Sampling cells to compose pseudo-bulk data...')

                results = Parallel(n_jobs=n_jobs)(
                    delayed(simulate_single_pseudobulk)(
                        cell_num[i],
                        df_work,
                        celltype_groups,
                        allcellname
                    )
                    for i in tqdm(range(cell_num.shape[0]))
                )

                sample = np.stack(results)

                sampleDF = pd.DataFrame(sample, index=[f"{condition}_sample_{i}" for i in range(1,prop.shape[0]+1)], columns = genes[0:sample.shape[1]])
                prop = pd.DataFrame(prop, index = [f"{condition}_sample_{i}" for i in range(1,prop.shape[0]+1)], columns=all_celltypes)
                
                if missing_cell_types == None:
                    print(f"All cell types present for individual")
                else:
                    prop[missing_cell_types] = 0

                # Append props and bulks from current iteration to list
                prop_list.append(prop)
                pseudobulk_list.append(sampleDF)
            
            # Concatenate all simulated data in a single df
            pseudobulks = pd.concat(pseudobulk_list, axis=0)
            proportions = pd.concat(prop_list, axis=0)
            
        elif target_name != None:
            print('Simulating cell type proportions and pseudobulks for specified condition.')
            
            df_work = scData.loc[scData['condition'] == target_name]
            df_work['CellType'] = df_work['CellType'].cat.remove_unused_categories()
            num_celltype = len(df_work['CellType'].value_counts())
            celltype_groups = df_work.groupby('CellType').groups # dictionary of cell types
            df_work = df_work.drop(columns=['CellType', 'individual', 'condition'])
                
            # Simulate proportions
            prop = simulate_proportions(props, random_state, d_prior, num_celltype, samplenum,
                                        sparse, sparse_prob, rare, rare_percentage)
                
            # precise number for each celltype
            n_adapt = int(np.floor(df_work.shape[0]/2)) # maximum number of cells, if cells per patient too few
            if n_adapt > n:
                n_adapt = n
                print(f"Using {str(n_adapt)} cells for sampling.")
            else: 
                print(f"Using {str(n_adapt)} cells for sampling.")
            cell_num = np.floor(n_adapt * prop) # from proportions get the number of cells per cell type to be sampled given n (e.g. n=1000 cells)

            # precise proportion based on cell_num
            prop = cell_num / np.sum(cell_num, axis=1).reshape(-1, 1) # then, obtain accurate proportions that fit the given number n

            # start sampling
            allcellname = celltype_groups.keys() # cell type names
            print('Sampling cells to compose pseudo-bulk data...')

            results = Parallel(n_jobs=n_jobs)(
                delayed(simulate_single_pseudobulk)(
                    cell_num[i],
                    df_work,
                    celltype_groups,
                    allcellname
                )
                for i in tqdm(range(cell_num.shape[0]))
            )

            sample = np.stack(results)

            pseudobulks = pd.DataFrame(sample, index=[f"{target_name}_sample_{i}" for i in range(1,prop.shape[0]+1)], columns = genes[0:sample.shape[1]])
            proportions = pd.DataFrame(prop, index = [f"{target_name}_sample_{i}" for i in range(1,prop.shape[0]+1)], columns=celltype_groups.keys())

        # Finalize simulated data in AnnData object
        simudata = anndata.AnnData(X=pseudobulks,obs=proportions)

    print('Sampling is done.')
    
    return simudata # anndata object containing ground-truth proportions in observations and having geneIDs as numbers for hidden layers


# Select mRNA genes

def main_gene_selection(X_df, gene_list):
    """
    Describe:
        rebuild the input adata to select target genes that encode proteins 
    Parameters:
        adata->`~anndata.AnnData` object: adata with var index_name by gene symbol
        gene_list->list: wanted target gene 
    Returns:
        adata_new->`~anndata.AnnData` object
        to_fill_columns->list: zero padding gene
    """
    X_df = X_df.fillna(0)
    to_fill_columns = list(set(gene_list) - set(X_df.columns))

    padding_df = pd.DataFrame(np.zeros((X_df.shape[0], len(to_fill_columns))), 
                              columns=to_fill_columns, 
                              index=X_df.index)
    
    X_df = pd.DataFrame(np.concatenate([df.values for df in [X_df, padding_df]], axis=1), 
                        index=X_df.index, 
                        columns=list(X_df.columns) + list(padding_df.columns))
    
    X_df = X_df[gene_list]
    
    var = pd.DataFrame(index=X_df.columns)
    var['mask'] = [1 if i in to_fill_columns else 0 for i in list(var.index)]
    
    return X_df, to_fill_columns, var


# Function for rank normalization

def normRank(pseudobulk):

    """
    Rank normalize gene expression per pseudobulk sample
    """

    df = np.array(pseudobulk.to_df())

    ranked_array = rankdata(df, axis = 1, method = 'min')
    ranked_array = ranked_array/ranked_array.shape[1]

    ranked_df = pd.DataFrame(ranked_array)
    ranked_df.index = pseudobulk.obs_names.to_list()
    ranked_df.columns = pseudobulk.var_names.to_list()

    pseudobulk.layers['ranked'] = ranked_df
    return pseudobulk


# Function for bulk normalization post gene selection

def bulk_norm(bulk, norm):

    """
    Function for bulk normalization
    """

    if norm == "CPM":
        
        # Convert to anndata, scale to CPM and convert back to df
        bulks = anndata.AnnData(X=bulk)
        sc.pp.normalize_total(bulks, target_sum=1e6)
        normbulk = pd.DataFrame(bulks.X, index=bulks.obs_names, columns=bulks.var_names)
    
    elif norm == "rank":
        
        # Convert to anndata, scale to CPM and convert back to df
        bulks = anndata.AnnData(X=bulk)
        bulks = normRank(bulks)
        normbulk = pd.DataFrame(bulks.layers['ranked'], index=bulks.obs_names, columns=bulks.var_names)

    elif norm == "log":

        # Convert to anndata, scale to CPM and convert back to df
        bulks = anndata.AnnData(X=bulk)
        sc.pp.log1p(bulks)
        normbulk = pd.DataFrame(bulks.X, index=bulks.obs_names, columns=bulks.var_names)

    elif norm == "none":
        normbulk = bulk

    return normbulk

# Function for gene selection and normalization of pseudobulks

def pseudobulk_norm(pseudobulks, norm, filter_genes):

    """
    Pseudobulk subsetting (filtering) and normalization as used in simulator
    """

    # mRNA_file_path 
    script_dir = pathlib.Path(__file__).resolve().parent

    # path to the text file
    file_path = script_dir / "mRNA_annotation.tsv"
    
    # CPM normalization option
    if norm == 'CPM':
        print('Scaling pseudobulks to CPM.')
        proportionsDF = pd.DataFrame(pseudobulks.obs) # does not need normalization --> proportions

        ##### All genes #####
        if filter_genes == "all":
            sc.pp.normalize_total(pseudobulks, target_sum=1e6)
            pseudobulkDF = pd.DataFrame(pseudobulks.X, index=pseudobulks.obs_names, columns=pseudobulks.var_names) 
            
        ##### Only mRNA genes #####
        elif filter_genes == "mRNA":
            pseudobulkDF = pd.DataFrame(pseudobulks.X, index=pseudobulks.obs_names, columns=pseudobulks.var_names) 
            proportionsDF = pd.DataFrame(pseudobulks.obs) # does not need normalization --> proportions

            # Import gene list for filtering
            gene_list_df = pd.read_csv(file_path, header=0, delimiter='\t')
            gene_list = list(gene_list_df['gene_name'])

            # Select (and add) genes as necessary
            pseudobulkDF, _, _ = main_gene_selection(pseudobulkDF, gene_list)
            
            # Convert to anndata, scale to CPM and convert back to df
            pseudobulks = anndata.AnnData(X=pseudobulkDF)
            sc.pp.normalize_total(pseudobulks, target_sum=1e6)
            pseudobulkDF = pd.DataFrame(pseudobulks.X, index=pseudobulks.obs_names, columns=pseudobulks.var_names)

        ##### remove zero variance genes #####
        elif filter_genes == "non_zero":
            df = pseudobulks.to_df()
            df = remove_zero_variance(df)
            subset_pseudobulks = anndata.AnnData(X=df.values, obs=pseudobulks.obs, var=pd.DataFrame(index=df.columns))
            sc.pp.normalize_total(subset_pseudobulks, target_sum=1e6)
            pseudobulkDF = pd.DataFrame(subset_pseudobulks.X, index=subset_pseudobulks.obs_names, columns=subset_pseudobulks.var_names)
            proportionsDF = pd.DataFrame(subset_pseudobulks.obs) # does not need normalization --> proportions

    elif norm == 'log':
        print('Log normalizing pseudobulks.')
        proportionsDF = pd.DataFrame(pseudobulks.obs) # does not need normalization --> proportions

        ##### All genes #####
        if filter_genes == "all":
            sc.pp.log1p(pseudobulks, target_sum=1e6)
            pseudobulkDF = pd.DataFrame(pseudobulks.X, index=pseudobulks.obs_names, columns=pseudobulks.var_names) 
            
        ##### Only mRNA genes #####
        elif filter_genes == "mRNA":
            pseudobulkDF = pd.DataFrame(pseudobulks.X, index=pseudobulks.obs_names, columns=pseudobulks.var_names) 
            proportionsDF = pd.DataFrame(pseudobulks.obs) # does not need normalization --> proportions

            # Import gene list for filtering
            gene_list_df = pd.read_csv(file_path, header=0, delimiter='\t')
            gene_list = list(gene_list_df['gene_name'])

            # Select (and add) genes as necessary
            pseudobulkDF, _, _ = main_gene_selection(pseudobulkDF, gene_list)
            
            # Convert to anndata, scale to CPM and convert back to df
            pseudobulks = anndata.AnnData(X=pseudobulkDF)
            sc.pp.log1p(pseudobulks, target_sum=1e6)
            pseudobulkDF = pd.DataFrame(pseudobulks.X, index=pseudobulks.obs_names, columns=pseudobulks.var_names)

        ##### remove zero variance genes #####
        elif filter_genes == "non_zero":
            df = pseudobulks.to_df()
            df = remove_zero_variance(df)
            subset_pseudobulks = anndata.AnnData(X=df.values, obs=pseudobulks.obs, var=pd.DataFrame(index=df.columns))
            sc.pp.log1p(subset_pseudobulks)
            pseudobulkDF = pd.DataFrame(subset_pseudobulks.X, index=subset_pseudobulks.obs_names, columns=subset_pseudobulks.var_names)
            proportionsDF = pd.DataFrame(subset_pseudobulks.obs) # does not need normalization --> proportions

    # Rank normalization option
    elif norm == 'rank':
        print('Ranking genes in pseudobulks.')
        proportionsDF = pd.DataFrame(pseudobulks.obs) # does not need normalization --> proportions
        
        ##### All genes #####
        if filter_genes == "all":
            pseudobulks = normRank(pseudobulks)
            pseudobulkDF = pd.DataFrame(pseudobulks.layers['ranked'], index=pseudobulks.obs_names, columns=pseudobulks.var_names) 


        ##### Only mRNA genes #####
        if filter_genes == "mRNA":
            # Import gene list for filtering
            gene_list_df = pd.read_csv(file_path, header=0, delimiter='\t')
            gene_list = list(gene_list_df['gene_name'])

            # Select (and add) genes as necessary
            pseudobulkDF, _, _ = main_gene_selection(pseudobulkDF, gene_list)

            # Convert to anndata, rank norm and convert back to df
            pseudobulks = anndata.AnnData(X=pseudobulkDF)
            pseudobulks = normRank(pseudobulks)
            pseudobulkDF = pd.DataFrame(pseudobulks.layers['ranked'], index=pseudobulks.obs_names, columns=pseudobulks.var_names)

        ##### remove zero variance genes #####
        elif filter_genes == "non_zero":
            df = pseudobulks.to_df()
            pseudobulkDF = remove_zero_variance(df)

            # Convert to anndata, rank norm and convert back to df
            pseudobulks = anndata.AnnData(X=pseudobulkDF)
            pseudobulks = normRank(pseudobulks)
            pseudobulkDF = pd.DataFrame(pseudobulks.layers['ranked'], index=pseudobulks.obs_names, columns=pseudobulks.var_names)
                
    # Raw counts option - no normalization
    elif norm == 'none':
        print('Returning raw summed counts in pseudobulks.')
        pseudobulkDF = pd.DataFrame(pseudobulks.X, index=pseudobulks.obs_names, columns=pseudobulks.var_names) 
        proportionsDF = pd.DataFrame(pseudobulks.obs) # does not need normalization --> proportions

        ##### Only mRNA genes #####
        if filter_genes == "mRNA":
            # Import gene list for filtering
            gene_list_df = pd.read_csv('mRNA_annotation.tsv', header=0, delimiter='\t')
            gene_list = list(gene_list_df['gene_name'])

            # Select (and add) genes as necessary
            pseudobulkDF, _, _ = main_gene_selection(pseudobulkDF, gene_list)

        ##### 3000 highly variable genes; not scaled #####
        elif filter_genes == "non_zero":
            df = pseudobulks.to_df()
            df = remove_zero_variance(df)
            pseudobulkDF = df
            proportionsDF = pd.DataFrame(subset_pseudobulks.obs) # does not need normalization --> proportions

    return proportionsDF, pseudobulkDF


def create_train_dir(path):

    """
    Create training and model directories 
    """

    char = "/"
    while path[-1] != char:
        path = path[:-1]
        
    pathlib.Path(path).mkdir(parents=True, exist_ok=True)
    pathlib.Path(path + 'train/model').mkdir(parents=True, exist_ok=True)

    train_path = path + 'train/'
    model_path = path + 'train/model/'

    return train_path, model_path


def create_pseudobulk_dir(path):

    """
    Create pseudobulk directory, for storing simulated pseudobulks
    """

    char = "/"
    while path[-1] != char:
        path = path[:-1]
        
    pathlib.Path(path).mkdir(parents=True, exist_ok=True)
    pathlib.Path(path + 'pseudobulk').mkdir(parents=True, exist_ok=True)

    pseudobulk_path = path + 'pseudobulk/'

    return pseudobulk_path


def create_pred_dir(path):

    """
    Create directory for storing inference/prediction results
    """

    char = "/"
    while path[-1] != char:
        path = path[:-1]
        
    pathlib.Path(path + 'Prediction').mkdir(parents=True, exist_ok=True)

    outpath = path + 'Prediction/'

    return outpath

def create_inference_dir(path):

    """
    Create directory for storing inference/prediction results
    """

    char = "/"
    while path[-1] != char:
        path = path[:-1]
        
    pathlib.Path(path + 'inference').mkdir(parents=True, exist_ok=True)

    outpath = path + 'inference/'

    return outpath


def remove_zero_variance(df: pd.DataFrame) -> pd.DataFrame:
    """
    Drop columns (genes) with zero variance
    
    Args:
        df (pd.DataFrame): DataFrame with samples as rows and genes as columns
        
    Returns:
        pd.DataFrame: DataFrame with zero-variance columns removed
    """
    variances = df.var(axis=0)
    keep = variances > 0
    return df.loc[:, keep]
