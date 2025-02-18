from src.general import *
import numpy as np
import pandas as pd
import sidetable
import matplotlib.pyplot as plt
import seaborn as sns


def read_data(file_path, index_col = None):
    assert file_path.exists()
    return pd.read_csv(file_path, index_col = index_col)


def insert_datetime(df, dayfirst=False):
    df.loc[df['Time'].str.startswith('0,'), 'Time'] = np.nan
    df.insert(0, 'datetime', df['Date'] + ' ' + df['Time'])
    df['datetime'] = pd.to_datetime(df.datetime, dayfirst=dayfirst)
    return df

def interpolate_wrong_times(df):
    df['datetime'] = df.datetime.values.astype('int64')
    df.loc[df.datetime < 0, 'datetime'] = np.nan
    df['datetime'] = pd.to_datetime(df.datetime.interpolate())
    return df

def clean_names(patches):
    patches.columns = patches.columns.str.capitalize()
    patches.columns = patches.columns.str.replace('_', ' ')

    return patches

def use_proxy_species(phylo: pd.DataFrame, species_dictionary: dict):

    for missing_species, proxy_species in species_dictionary.items():
        phylo.loc[missing_species] = phylo.loc[proxy_species]
        phylo.loc[:, missing_species] = phylo.loc[:, proxy_species]

    return phylo

def binarize(patches):
    return patches.gt(0)

def sort_species(patches):
    species = patches.columns
    species_by_abundance = patches[species].mean().sort_values(ascending=False).index
    return patches[species_by_abundance]

def remove_boring(patches, cutoff = 1.5):
    return patches[patches.sum(axis='columns') >= cutoff].reset_index(drop=True)


def keep_species(patches, n_keep: int = 16):
    patches = patches.iloc[:, :n_keep]
    return remove_boring(patches, cutoff=0.5)

def crop_clean(patches, n_keep: int = 16, cutoff: float = 1.5):
    return patches.pipe(keep_species, n_keep).pipe(remove_boring, cutoff).reset_index(drop=True)

def export_patches(patches: pd.DataFrame, paths: dict, file_name: str='', index : str=False):
    path = paths['outputs'] / file_name
    patches.to_csv(str(path), index=index)
    return path


def compute_common_patches(patches1, patches2):

    species = patches1.columns.tolist()
    patches2.columns = species

    s = patches1.stb.freq(species).set_index(species)['percent']
    s2 = patches2.stb.freq(species).set_index(species)['percent']
    
    common = s.to_frame().join(s2.to_frame(), how='inner', lsuffix='_l', rsuffix='_r').min(axis=1).sum()
    
    return common/100


def generate_synthetic_patches(patches: pd.DataFrame, n_patches : int=None):

    if n_patches == None:
        n_patches = patches.shape[0]
    
    synthetic_patches = pd.DataFrame()
    for i, species in enumerate(patches.columns):
        synthetic_patches[str(i+1)] = np.random.default_rng().binomial(1, patches[species].mean(), n_patches)
    
    return synthetic_patches


def get_patches(config: dict,
                data : str = 'testing_data'
                ):
    
    csvs = [p for p in (get_git_root() / 'outputs' / config[data]).glob('*.csv')]
    csv = csvs[-1]
    df = pd.read_csv(csv)
    # display(df)
    patches = df#.pipe(binarize)
    if config['species_order'] != 'phylogeny':
        patches = (patches
                        .pipe(sort_species)
                        .pipe(keep_species, config['n_species'])
                        .pipe(remove_boring, config['min_species'])
                        .pipe(sort_species)
                  )

        if config['species_order'] == 'reverse_abundance':
            patches = patches[reversed(list(patches.columns))]

        if config['species_order'] == 'random':
            patches = patches.sample(frac=1, axis=1)
            
        else:
            pass
    
    return patches


def get_species(location : 'str'):
    
    paths = set_paths(location)
    csvs = [p for p in paths['outputs'].glob('*.csv')]
    csv = csvs[0]
    if len(csvs) > 1:
        print(f'WARNING: List of CSVs found in  directory:\n{paths["outputs"]}\nUsing {csv}\nPlease make sure this is correct')
    return read_data(csv).pipe(sort_species).pipe(clean_names).columns


def subset_phylo(phylo : pd.DataFrame,
                 species_origin : list,
                 species_target : list,
                 from_species : int = 0,
                 to_species : int = 16):
    phylo_subset = phylo.loc[:, species_origin[from_species:to_species]].reindex(species_target).iloc[from_species:to_species, :]
    # phylo_subset.plot(kind='imshow', height = 800).show()
    return phylo_subset


def match_and_remove(phylo_subset : pd.DataFrame,
                     sp1 : str,
                     sp2 : str,
                     dictionary : dict):
    dictionary[sp2] = sp1
    phylo_subset = phylo_subset.drop(index=sp1).drop(columns=sp2)
    return phylo_subset, dictionary


def build_dictionary(config : dict,
                     phylo : pd.DataFrame,
                     species_origin : list,
                     species_target : list,
                     chunk_size : int = 16,
                     multi_chunks : bool = False):

    dictionary = {}
    
    if multi_chunks != True:
        
        for i in range(chunk_size):
            dictionary[species_origin[i]] = species_target[i] 
        
        phylo_subset = subset_phylo(phylo, species_origin, species_target, from_species = chunk_size, to_species = None)
        
        # # first loop to find exact matches
        for sp in phylo_subset.columns:
            try:
                phylo_subset, dictionary = match_and_remove(phylo_subset, sp, sp, dictionary)
                # print(f'{sp} : {sp}')
            except:
                continue

        # #second loop to find closest matches
        for sp in phylo_subset.columns:
            closest_species = phylo_subset.index[phylo_subset[sp].argmin()]
            phylo_subset, dictionary = match_and_remove(phylo_subset, closest_species, sp, dictionary)
            # print(f'{closest_species} : {sp}')

        
    else:
        n_chunks = config['n_species'] // chunk_size + 1
        remaining_species = config['n_species'] % chunk_size


        for n_chunk in range(0, n_chunks):

            # display(phylo)
            phylo_subset = subset_phylo(phylo, species_origin, species_target, from_species = n_chunk * chunk_size, to_species = (n_chunk + 1) * chunk_size)

            if n_chunk == n_chunks - 1:
                chunk_size = remaining_species
                if remaining_species == 0:
                    continue
            else:
                chunk_size = chunk_size

            # print(f'\nMatching chunk number {n_chunk} with {chunk_size} species\n')

            # # first loop to find exact matches
            for sp in phylo_subset.columns:
                try:
                    phylo_subset, dictionary = match_and_remove(phylo_subset, sp, sp, dictionary)
                    # print(f'{sp} : {sp}')
                except:
                    continue

            # #second loop to find closest matches
            for sp in phylo_subset.columns:
                closest_species = phylo_subset.index[phylo_subset[sp].argmin()]
                phylo_subset, dictionary = match_and_remove(phylo_subset, closest_species, sp, dictionary)
                # print(f'{closest_species} : {sp}')


    return dictionary