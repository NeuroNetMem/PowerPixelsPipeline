# -*- coding: utf-8 -*-
"""

Written by Guido Meijer

"""

from powerpixels import Pipeline

import os
import numpy as np
from datetime import datetime
import spikeinterface.full as si
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
        if hasattr(pp, 'nidq_file'):
            pp.nidq_synchronization()
        
        # Loop over multiple probes 
        stream_names, _ = si.read_openephys(pp.session_path, stream_id='0').get_streams(pp.session_path)
        probes = np.unique([i[i.find('Probe'):i.find('Probe')+6] for i in stream_names])
        probe_done = np.zeros(len(probes)).astype(bool)
        for i, this_probe in enumerate(probes):
            print(f'\nStarting preprocessing of {this_probe}')
            
            # Set probe paths
            pp.set_probe_paths(this_probe)
            
            # Check if probe is already processed
            if pp.results_path.is_dir():
                print('Probe already processed, moving on')
                probe_done[i] = True
                continue
                        
            # Preprocessing
            rec = pp.preprocessing()
            
            # Spike sorting
            print(f'\nStarting {this_probe} spike sorting at {datetime.now().strftime("%H:%M")}')
            sort = pp.spikesorting(rec)   
            if sort is None:
                print('Spike sorting failed!')
                continue
            print(f'Detected {sort.get_num_units()} units\n')      
                                   
            # Create sorting analyzer for manual curation in SpikeInterface and save to disk
            pp.neuron_metrics(sort, rec)
                        
            # Export sorting results and LFP metrics
            pp.export_data(rec)
            
            # Add indication if neurons are good from several sources to the quality metrics
            pp.automatic_curation()
            
            # Synchronize spike sorting to the nidq clock
            if hasattr(pp, 'nidq_file'):
                pp.probe_synchronization()
            
            # Compress raw data 
            pp.compress_raw_data()            
                        
            probe_done[i] = True
            print(f'Done! At {datetime.now().strftime("%H:%M")}')
        
        # Delete process_me.flag if all probes are processed
        if np.sum(probe_done) == len(probes):
            os.remove(os.path.join(root, 'process_me.flag'))
   
