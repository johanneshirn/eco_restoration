from pathlib import Path
import subprocess
from tqdm import tqdm
import pandas as pd

def configure_defaults():
    config = {}
    config['n_species'] = 16
    config['min_species'] = 1
    config['species_order'] = 'abundance'
    config['training_data'] = 'petrer_limestone'
    config['testing_data'] = 'petrer_limestone'
    
    # config['training_start'] = 'tune_petrer_model'
    config['training_start'] = 'train_from_scratch'

    config['test_split'] = 0.33
    config['val_split'] = 0.5
    config['reconstruction_loss'] = 'bce'
    # config['reconstruction_loss'] = 'mse'
    config['depth'] = 32
    config['n_latent'] = 16
    config['beta'] = 1
    config['monitor'] = 'val_loss'
    config['mode'] = 'min'

    config['learning_rate'] = 1e-3 #1e-3
    config['batch_size'] = 64
    config['epochs'] = 200
    config['stride'] = 2

    config['n_bootstrap'] = 20

    config['n_tile'] = 8
    config['n_colors'] = 1
    config['fully_connected'] = False
    config['testing'] = False
    
    config['other_species_handling'] = 'drop_species'
    
    config['chunk_size'] = 1
    
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