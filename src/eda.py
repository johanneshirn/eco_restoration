from src.general import *
import numpy as np
import pandas as pd
import sidetable
import matplotlib.pyplot as plt
# import seaborn as sns


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


def get_patches(config: dict,
                data : str = 'testing_data'
                ):
    
    csvs = [p for p in (get_git_root() / 'outputs' / config[data]).glob('*.csv')]
    csv = csvs[-1]
    df = pd.read_csv(csv)
    # display(df)
    patches = df#.pipe(binarize)
    
    return patches


def get_species(location : 'str'):
    
    paths = set_paths(location)
    csvs = [p for p in paths['outputs'].glob('*.csv')]
    csv = csvs[0]
    if len(csvs) > 1:
        print(f'WARNING: List of CSVs found in  directory:\n{paths["outputs"]}\nUsing {csv}\nPlease make sure this is correct')
    return read_data(csv).pipe(sort_species).pipe(clean_names).columns