# -*- coding: utf-8 -*-
"""
Generate settings and wiring JSON files to be used in the Megazord pipeline.

Settings JSON file
------------------------
SPIKE SORTER
    Which spike sorter SpikeInterface should use, default is Kilosort 2.5, for options see:
    https://spikeinterface.readthedocs.io/en/latest/modules/sorters.html#supported-spike-sorters
IDENTIFIER 
    An optional string which can be used to identify spike sorting runs with different settings
    Will be appended to the name of the output directory
DATA_FOLDER 
    Path to the top-level data folder (make sure to use \\ when on Windows machines)
USE_DOCKER
    Whether to use Docker while running the spike sorting or using a local installation of Kilosort
COMPRESS_RAW_DATA
    Whether to compress the raw bin files
N_CORES
    How many CPU cores to use (default is -1 which is all of them)
    
NIDQ wiring JSON file
------------------------
SYSTEM
    Should not be changed from 3B (otherwise known as 1.0). Currently Megazord only supports 
    Neuropixel 1.0
SYNC_WIRING_DIGITAL
    Dictionary with all the digital synchronization channels on the BNC breakout box that should be
    extracted and timestamped. One of these channels should be the square wave sync pulse from 
    the PXI chassis, in this example this is channel number 3.
SYNC_WIRING_ANALOG
    Dictionary with the wiring of any analog inputs that should be synchronized with the recording.
    
Probe wiring JSON file
------------------------
Should not be changed

Created on Mon Dec 4 2023 by Guido Meijer
"""

import json
from spikeinterface.sorters import get_default_sorter_params
from pathlib import Path
 
def main():
    
    # Get paths to where to save the configuration files (in the repository)
    project_root = Path(__file__).parent.parent.parent
    settings_file = project_root / 'config' / 'settings.json'
    wiring_dir = project_root / 'config' /'wiring'
    sorting_dir = project_root / 'config' / 'sorter_params'
    
    # Settings
    if settings_file.is_file():
        print(f'\nConfiguration file already exists at {settings_file}')
    else:
        
        # Generate example settings JSON file
        settings_dict = {
            "SPIKE_SORTER": "kilosort4",  
            "IDENTIFIER": "",
            "DATA_FOLDER": "C:\\path\\to\\data",
            "USE_DOCKER": False,
            "COMPRESS_RAW_DATA": True,
            "N_CORES": -1
        }
        with open(settings_file, 'w') as outfile:
            outfile.write(json.dumps(settings_dict, indent=4))
        
        print(f'\nExample configuration file generated at {settings_file}')
            
    # Wiring files
    if wiring_dir.is_dir():
        print(f'\nDirectory with wiring files already exists at {wiring_dir}')
    else:
        
        # NIDQ wiring JSON file
        nidq_wiring_dict = {
            "SYSTEM": "3B",
            "SYNC_WIRING_DIGITAL": {
                "P0.0": "imec_sync",
                "P0.1": "lick_detector",
                "P0.2": "camera"
            },
            "SYNC_WIRING_ANALOG": {
                "AI0": "breathing_sensor"
            }
        }
        
        # Neuropixel probe wiring JSON file
        probe_wiring_dict = {
            "SYSTEM": "3B",
            "SYNC_WIRING_DIGITAL": {
                "P0.6": "imec_sync"
            }
        }
        
        # Save example wiring configuration files to repo dir
        wiring_dir.mkdir()
        with open(wiring_dir / 'nidq.wiring.json', 'w') as outfile:
            outfile.write(json.dumps(nidq_wiring_dict, indent=4))
        with open(wiring_dir / '3B.wiring.json', 'w') as outfile:
            outfile.write(json.dumps(probe_wiring_dict, indent=4))
            
        print(f'\nExample wiring files generated in {wiring_dir}')
        
    # Spike sorter parameter files
    if sorting_dir.is_dir():
        print(f'\nSpike sorter parameter files already exist in {sorting_dir}')
    else:
        
        # Get default sorter params
        sorting_dir.mkdir()
        for sorter in ['kilosort2_5', 'kilosort3', 'kilosort4', 'pykilosort']:
            
            # Save to disk
            sorter_params = get_default_sorter_params(sorter)
            with open(sorting_dir / f'{sorter}_params.json', 'w') as outfile:
                outfile.write(json.dumps(sorter_params, indent=4))
        print(f'\nDefault settings for spike sorters generated in {sorting_dir}')
        
    # Example notch filter file
    if (project_root / 'config' / 'notch_filter.json').is_file():
        print('\nExample notch filter configuration file already exists.')
    else:
        notch_filter = {
            "FREQ": [4000, 12000],
            "Q": [20, 10]
        }
        with open(project_root / 'config' / 'notch_filter.json', 'w') as outfile:
            outfile.write(json.dumps(notch_filter, indent=4))
        print(f'\nExample notch filter file generated in {project_root / "config"}')
        
if __name__ == "__main__":
    main()
        
        
