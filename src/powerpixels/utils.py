#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul  4 11:42:07 2025

By Guido Meijer
"""

import numpy as np
import pandas as pd
import json
from pathlib import Path
from spikeinterface import load_sorting_analyzer
from spikeinterface.widgets import plot_sorting_summary

def manual_curation(results_path):
    """
    Launch the manual curation GUI

    Parameters
    ----------
    results_path : Path
        Path to where the final results of the spike sorting is saved.

    Returns
    -------
    None.

    """
    
    # Transfrom path into Pathlib if necessary
    if type(results_path) == str:
        results_path = Path(results_path)
    
    # Load in sorting analyzer from disk
    sorting_analyzer = load_sorting_analyzer(results_path / 'sorting')
    
    # Launch manual curation GUI       
    unit_properties = ['Bombcell', 'UnitRefine', 'IBL', 'Kilosort', 'firing_rate', 'rp_violations',
                       'snr', 'amplitude_median', 'presence_ratio']
    _ = plot_sorting_summary(sorting_analyzer=sorting_analyzer, curation=True,
                             displayed_unit_properties=unit_properties,
                             backend='spikeinterface_gui')
    
    # Extract manual curation labels and save in results folder
    if (results_path / 'sorting' / 'spikeinterface_gui' / 'curation_data.json').is_file():
        with open(results_path / 'sorting' / 'spikeinterface_gui' / 'curation_data.json') as f:
            label_dict = json.load(f)
        if (results_path / 'clusters.manualLabels.npy').is_file():
            manual_labels = np.load(results_path / 'clusters.manualLabels.npy')
        else:
            manual_labels = np.array(['no label'] * sorting_analyzer.unit_ids.shape[0])
        for this_unit in label_dict['manual_labels']:
            manual_labels[sorting_analyzer.unit_ids == this_unit['unit_id']] = this_unit['quality']
        np.save(results_path / 'clusters.manualLabels.npy', manual_labels)
    
    return       
            

def load_neural_data(session_path, probe, histology=False, keep_units='all'):
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
        If False, no brain regions will be provided. The default is False.
    keep_units : str, optional
        Which units to keep
        'all' = keep all units (default)
        'bombcell' = keep units classified as good by Bombcell
        'unitrefine' = keep units classified as good by UnitRefine
        'ibl' = keep units classified as good by IBL metrics
        'kilosort' = keep Kilosort good neurons
        'manual' = keep units manually annotated as good in the GUI
    Returns
    -------
    spikes : dict
        A dictionary containing data per spike
    clusters : dict
        A dictionary containing data per cluster (i.e. neuron)
    channels : dict
        A dictionary containing data per channel 
    """
    
    # Convert path to Pathlib if necessary
    if type(session_path) == str:
        session_path = Path(session_path)
    
    # Load in spiking data
    spikes = dict()
    spikes['times'] = np.load(session_path / probe / 'spikes.times.npy')
    spikes['clusters'] = np.load(session_path / probe / 'spikes.clusters.npy')
    spikes['amps'] = np.load(session_path /probe / 'spikes.amps.npy')
    spikes['amps'] = np.clip(spikes['amps'], 1e-5, None)  # remove negative spike amplitudes
    spikes['depths'] = np.load(session_path / probe / 'spikes.depths.npy')
    
    # Load in cluster data
    clusters = dict()
    clusters['channels'] = np.load(session_path / probe / 'clusters.channels.npy')
    clusters['cluster_id'] = np.arange(clusters['channels'].shape[0])
    
    # Add cluster qc metrics
    if (session_path / probe / 'clusters.bombcellLabels.npy').is_file():
        clusters['bombcell_label'] = np.load(session_path / probe / 'clusters.bombcellLabels.npy')
    if (session_path / probe / 'clusters.unitrefineLabels.npy').is_file():
        clusters['unitrefine_label'] = np.load(session_path / probe / 'clusters.unitrefineLabels.npy')
    elif (session_path / probe / 'clusters.MLLabel.npy').is_file():  # legacy
        clusters['unitrefine_label'] = np.load(session_path / probe / 'clusters.MLLabel.npy')
    if (session_path / probe / 'clusters.iblLabels.npy').is_file():
        clusters['ibl_label'] = np.load(session_path / probe / 'clusters.iblLabels.npy')
    elif (session_path / probe / 'clusters.IBLLabel.npy').is_file():  # legacy
        clusters['ibl_label'] = np.load(session_path / probe / 'clusters.IBLLabel.npy')
    if (session_path / probe / 'clusters.kilosortLabels.npy').is_file():
        clusters['kilosort_label'] = np.load(session_path / probe / 'clusters.kilosortLabels.npy')
    elif (session_path / probe / 'clusters.KSLabel.npy').is_file():  # legacy
        clusters['kilosort_label'] = np.load(session_path / probe / 'clusters.KSLabel.npy')
    if (session_path / probe / 'clusters.manualLabels.npy').is_file():
        clusters['manual_label'] = np.load(session_path / probe / 'clusters.manualLabels.npy')
        
    # Load in channel data
    channels = dict()
    if histology:
        if not (session_path / probe / 'channel_locations.json').is_file():
            raise Exception('No aligned channel locations found! Set histology to False to load data without brain regions.')
        
        # Load in alignment GUI output
        f = open(session_path / probe / 'channel_locations.json')
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
    local_coordinates = np.load(session_path / probe / 'channels.localCoordinates.npy')
    channels['lateral_um'] = local_coordinates[:, 0]
    channels['axial_um'] = local_coordinates[:, 1]
        
    # Only keep the neurons that are manually labeled as good
    if keep_units == 'all':
        return spikes, clusters, channels
    if keep_units == 'bombcell':
        if 'bombcell_label' not in clusters.keys():
            raise Exception('No Bombcell labels found! Set keep_units to "all" to load all neurons.')
        good_units = np.where(clusters['bombcell_label'] == 1)[0]
    elif keep_units == 'unitrefine':
        if 'unitrefine_label' not in clusters.keys():
            raise Exception('No UnitRefine labels found! Set keep_units to "all" to load all neurons.')
        good_units = np.where(clusters['unitrefine_label'] == 1)[0]
    elif keep_units == 'ibl':
        if 'ibl_label' not in clusters.keys():
            raise Exception('No IBL labels found! Set keep_units to "all" to load all neurons.')
        good_units = np.where(clusters['ibl_label'] == 1)[0]
    elif keep_units == 'kilosort':
        if 'kilosort_label' not in clusters.keys():
            raise Exception('No Kilosort labels found! Set keep_units to "all" to load all neurons.')
        good_units = np.where(clusters['kilosort_label'] == 1)[0]
    elif keep_units == 'manual':
        if 'manual_label' not in clusters.keys():
            raise Exception('No manual cluster labels found! Set keep_units to "all" to load all neurons.')
        good_units = np.where(clusters['manual_label'] == 'good')[0]
    else:
        raise Exception('keep_units shoud be all, bombcell, unitrefine, ibl or manual')
      
    spikes['times'] = spikes['times'][np.isin(spikes['clusters'], good_units)]
    spikes['amps'] = spikes['amps'][np.isin(spikes['clusters'], good_units)]
    spikes['depths'] = spikes['depths'][np.isin(spikes['clusters'], good_units)]
    spikes['clusters'] = spikes['clusters'][np.isin(spikes['clusters'], good_units)]
    clusters['cluster_id'] = clusters['cluster_id'][good_units]
    if histology:
        clusters['acronym'] = clusters['acronym'][good_units]
        
    return spikes, clusters, channels
        
        