# -*- coding: utf-8 -*-
"""
Created on Tue Feb 20 15:42:00 2024

By Guido Meijer
"""

import numpy as np
import pandas as pd
import json
from os.path import join, isfile


def load_neural_data(session_path, probe, histology=True, only_good=True):
    """
    Helper function to read in the spike sorting output from the Power Pixels pipeline.

    Parameters
    ----------
    session_path : str
        Full path to the top-level folder of the session.
    probe : str
        Name of the probe to load in.
    histology : bool, optional
        Whether to load the channel location and brain regions from the output of the alignment GUI.
        If False, no brain regions will be provided. The default is True.
    only_good : bool, optional
        Whether to only load in neurons that have been manually labelled in Phy.
        The default is True.

    Returns
    -------
    spikes : dict
        A dictionary containing data per spike
    clusters : dict
        A dictionary containing data per cluster (i.e. neuron)
    channels : dict
        A dictionary containing data per channel 
    """
    
    # Load in spiking data
    spikes = dict()
    spikes['times'] = np.load(join(session_path, probe, 'spikes.times.npy'))
    spikes['clusters'] = np.load(join(session_path, probe, 'spikes.clusters.npy'))
    spikes['amps'] = np.load(join(session_path, probe, 'spikes.amps.npy'))
    spikes['depths'] = np.load(join(session_path, probe, 'spikes.depths.npy'))
    
    # Load in cluster data
    clusters = dict()
    clusters['channels'] = np.load(join(session_path, probe, 'clusters.channels.npy'))
    clusters['depths'] = np.load(join(session_path, probe, 'clusters.depths.npy'))
    clusters['amps'] = np.load(join(session_path, probe, 'clusters.amps.npy'))
    clusters['cluster_id'] = np.arange(clusters['channels'].shape[0])
    
    # Add cluster qc metrics
    if isfile(join(session_path, probe, 'clusters.bcUnitType.npy')):
        clusters['bc_label'] = np.load(join(session_path, probe, 'clusters.bcUnitType.npy'),
                                       allow_pickle=True)
    clusters['ks_label'] = pd.read_csv(join(session_path, probe, 'cluster_KSLabel.tsv'),
                                       sep='\t')['KSLabel']
    if isfile(join(session_path, probe, 'clusters.iblLabel.tsv')):
        clusters['ibl_label'] = pd.read_csv(join(session_path, probe, 'cluster_IBLLabel.tsv'),
                                            sep='\t')['ibl_label']
    if isfile(join(session_path, probe, 'cluster_group.tsv')):
        clusters['manual_label'] = pd.read_csv(join(session_path, probe, 'cluster_group.tsv'),
                                               sep='\t')['group']
    # Load in channel data
    channels = dict()
    if histology:
        if not isfile(join(session_path, probe, 'channel_locations.json')):
            raise Exception('No aligned channel locations found! Set histology to False to load data without brain regions.')
        
        # Load in alignment GUI output
        f = open(join(session_path, probe, 'channel_locations.json'))
        channel_locations = json.load(f)
        f.close()
        
        # Add channel information to channel dict        
        brain_region, brain_region_id, x, y, z = [], [], [], [], []
        for i, this_ch in enumerate(channel_locations.keys()):
            if this_ch[:7] != 'channel':
                continue
            brain_region.append(channel_locations[this_ch]['brain_region'])
            brain_region_id.append(channel_locations[this_ch]['brain_region_id'])
            x.append(channel_locations[this_ch]['x'])
            y.append(channel_locations[this_ch]['y'])
            z.append(channel_locations[this_ch]['z'])
        channels['acronym'] = np.array(brain_region)
        channels['atlas_id'] = np.array(brain_region_id)
        channels['x'] = np.array(x)
        channels['y'] = np.array(y)
        channels['z'] = np.array(z)
        
        # Use the channel location to infer the brain regions of the clusters
        clusters['acronym'] = channels['acronym'][clusters['channels']]
            
    # Load in the local coordinates of the probe
    local_coordinates = np.load(join(session_path, probe, 'channels.localCoordinates.npy'))  
    channels['lateral_um'] = local_coordinates[:, 0]
    channels['axial_um'] = local_coordinates[:, 1]
        
    # Only keep the neurons that are manually labeled as good
    if only_good:
        if 'manual_label' not in clusters.keys():
            raise Exception('No manual cluster labels found! Set only_good to False to load all neurons.')
        good_units = np.where(clusters['manual_label'] == 'good')[0]
        spikes['times'] = spikes['times'][np.isin(spikes['clusters'], good_units)]
        spikes['amps'] = spikes['amps'][np.isin(spikes['clusters'], good_units)]
        spikes['depths'] = spikes['depths'][np.isin(spikes['clusters'], good_units)]
        spikes['clusters'] = spikes['clusters'][np.isin(spikes['clusters'], good_units)]
        clusters['acronym'] = clusters['acronym'][good_units]
        clusters['depths'] = clusters['depths'][good_units]
        clusters['amps'] = clusters['amps'][good_units]
        clusters['cluster_id'] = clusters['cluster_id'][good_units]
    
    return spikes, clusters, channels
