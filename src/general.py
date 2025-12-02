from pathlib import Path
import subprocess
from tqdm import tqdm
import pandas as pd

def configure_defaults():
    config = {}
    config['n_species'] = 8
    config['min_species'] = 1
    config['dataset'] = 'Escuzar'
    config['n_bootstrap'] = 20
    config['other_species_handling'] = 'gather_species'
        
    return config

def get_git_root():
    return Path(subprocess.Popen(['git', 'rev-parse', '--show-toplevel'], stdout=subprocess.PIPE).communicate()[0].rstrip().decode('utf-8'))

def create_dir(path):
    path.mkdir(parents=True, exist_ok=True)

def set_paths(location: str) -> dict:
    root = get_git_root()
    paths = {'data' : root / 'data' / location,
            'outputs' : root / 'outputs' / location,
            'figures' : root / 'figures' / location,
            'models' : root / 'models' / location}
    create_dir(paths['outputs'])
    # [create_dir(p) for p in paths.values()]
    return paths 


def show(config: dict):
    
    return pd.DataFrame([config]).T