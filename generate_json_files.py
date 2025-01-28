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
from os import mkdir
from os.path import join, dirname, realpath, isdir
 
# Settings JSON file
settings_dict = {
    "SPIKE_SORTER": "kilosort2_5",  
    "IDENTIFIER": "",
    "DATA_FOLDER": "C:\\path\\to\\data",
    "USE_DOCKER": True,
    "COMPRESS_RAW_DATA": True,
    "N_CORES": -1
}

# NIDQ wiring JSON file
nidq_wiring_dict = {
    "SYSTEM": "3B",
    "SYNC_WIRING_DIGITAL": {
        "P0.0": "example_sync_channel",
        "P0.3": "imec_sync"
    },
    "SYNC_WIRING_ANALOG": {
        "AI0": "example_sync_channel",
        "AI1": "example_sync_channel"
    }
}

# Neuropixel 1.0 probe wiring JSON file
probe_3B_wiring_dict = {
    "SYSTEM": "3B",
    "SYNC_WIRING_DIGITAL": {
        "P0.6": "imec_sync"
    }
}

# Save JSON files to repository directory
if not isdir(join(dirname(realpath(__file__)), 'wiring_files')):
    mkdir(join(dirname(realpath(__file__)), 'wiring_files'))
with open(join(dirname(realpath(__file__)), 'settings.json'), 'w') as outfile:
    outfile.write(json.dumps(settings_dict, indent=4))
with open(join(dirname(realpath(__file__)), 'wiring_files', 'nidq.wiring.json'), 'w') as outfile:
    outfile.write(json.dumps(nidq_wiring_dict, indent=4))
with open(join(dirname(realpath(__file__)), 'wiring_files', '3B.wiring.json'), 'w') as outfile:
    outfile.write(json.dumps(probe_3B_wiring_dict, indent=4))
