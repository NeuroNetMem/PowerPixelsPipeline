# -*- coding: utf-8 -*-
"""
Created on Wed Apr  5 14:02:41 2023 

by Guido Meijer
"""

from powerpixels import Pipeline

import os
from os.path import join, split, isdir
import numpy as np
from datetime import datetime
from glob import glob
from pathlib import Path

# Initialize power pixels pipeline
pp = Pipeline()
    
# Search for process_me.flag
print('Looking for process_me.flag..')
for root, directory, files in os.walk(pp.settings['DATA_FOLDER']):
    if 'process_me.flag' in files:
        print(f'\nStarting pipeline in {root} at {datetime.now().strftime("%H:%M")}\n')
        
        # Set session path
        pp.session_path = Path(root)
        
        # Restructure file and folders
        pp.restructure_files()
               
        # Create nidq synchronization files
        pp.nidq_synchronization()
        
        # Loop over multiple probes 
        probes = glob(join(root, 'raw_ephys_data', 'probe*'))
        probe_done = np.zeros(len(probes)).astype(bool)
        for i, probe_path in enumerate(probes):
            print(f'\nStarting preprocessing of {split(probe_path)[-1]}')
            
            # Set probe paths
            pp.set_probe_paths(Path(probe_path))
            
            # Check if probe is already processed
            if isdir(join(pp.session_path, pp.this_probe + pp.settings['IDENTIFIER'])):
                print('Probe already processed, moving on')
                probe_done[i] = True
                continue
            
            # Preprocessing
            rec = pp.preprocessing()
            
            # Spike sorting
            print(f'\nStarting {split(probe_path)[-1]} spike sorting at {datetime.now().strftime("%H:%M")}')
            sort = pp.spikesorting(rec, probe_path)   
            print(f'Detected {sort.get_num_units()} units\n')      
                                   
            # Create sorting analyzer for manual curation in SpikeInterface and save to disk
            pp.neuron_metrics(sort, rec)
            
            # Calculate raw ephys QC metrics
            pp.raw_ephys_qc()
            
            # Convert Kilosort output to ALF file format and move to results folder
            pp.convert_to_alf()
            asd
            # Add indication if neurons are good from several sources to the quality metrics
            pp.automatic_curation()
            
            # Synchronize spike sorting to the nidq clock
            pp.probe_synchronization()
            
            # Compress raw data (still hase issues, don't run)
            #pp.compress_raw_data()            
                        
            probe_done[i] = True
            print(f'Done! At {datetime.now().strftime("%H:%M")}')
        
        # Delete process_me.flag if all probes are processed
        if np.sum(probe_done) == len(probes):
            os.remove(os.path.join(root, 'process_me.flag'))
