#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Written by Guido Meijer

"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import os
from os.path import join, isfile, split, isdir, dirname, realpath
from pathlib import Path
import shutil
from glob import glob
import json

import spikeinterface.full as si

from one.api import ONE
from neuropixel import NP2Converter
import mtscomp
from atlaselectrophysiology.extract_files import extract_rmsmap
from brainbox.metrics.single_units import spike_sorting_metrics
from ibllib.ephys import ephysqc
from ibllib.ephys.spikes import ks2_to_alf, sync_spike_sorting
from ibllib.pipes.ephys_tasks import (EphysCompressNP1, EphysCompressNP21, EphysSyncPulses,
                                      EphysSyncRegisterRaw, EphysPulses)


class Pipeline:
    
    def __init__(self):
        
        # Load in setting files
        with open(join(dirname(realpath(__file__)), 'settings.json'), 'r') as openfile:
            self.settings = json.load(openfile)
        with open(join(dirname(realpath(__file__)), 'wiring_files', 'nidq.wiring.json'), 'r') as openfile:
            self.nidq_sync = json.load(openfile)
        with open(join(dirname(realpath(__file__)), 'wiring_files',
                       f'{self.nidq_sync["SYSTEM"]}.wiring.json'), 'r') as openfile:
            self.probe_sync = json.load(openfile)

        # Initialize spikeinterface parallel processing
        si.set_global_job_kwargs(n_jobs=self.settings['N_CORES'], progress_bar=True)
        
        # Load in spike sorting parameters
        if isfile(join(dirname(realpath(__file__)), 'spikesorter_param_files',
                       f'{self.settings["SPIKE_SORTER"]}_params.json')):
            with open(join(dirname(realpath(__file__)), 'spikesorter_param_files',
                           f'{self.settings["SPIKE_SORTER"]}_params.json'), 'r') as openfile:
                self.sorter_params = json.load(openfile)
        else:
            self.sorter_params = si.get_default_sorter_params(self.settings['SPIKE_SORTER'])
            
        # Initialize ONE connection (needed for some IBL steps for some reason)
        ONE.setup(base_url='https://openalyx.internationalbrainlab.org', silent=True)
        ONE(password='international')
            
        
    def set_probe_paths(self, probe_path):
        
        self.probe_path = probe_path
        self.sorter_out_path = Path(join(
            probe_path,
            self.settings['SPIKE_SORTER'] + self.settings['IDENTIFIER'],
            'sorter_output')
            )
        self.this_probe = split(probe_path)[1]
        self.results_path = Path(join(self.session_path,
                                      self.this_probe + self.settings['IDENTIFIER']))
        self.ap_file = Path(glob(join(self.probe_path, '*ap.*bin'))[0])
            
        return
    
    
    def restructure_files(self):
        """
        Restructure the raw data files from SpikeGLX (OpenEphys not supported)

        """
        
        # Restructure file and folders
        if len([i for i in os.listdir(join(self.session_path, 'raw_ephys_data')) if i[:5] == 'probe']) == 0:
            if len(os.listdir(join(self.session_path, 'raw_ephys_data'))) == 0:
                print('No ephys data found')
                return
            elif len(os.listdir(join(self.session_path, 'raw_ephys_data'))) > 1:
                print('More than one run found, not supported')
                return 
            orig_dir = os.listdir(join(self.session_path, 'raw_ephys_data'))[0]
            if orig_dir[-2] != 'g':
                print('Recording is not in SpikeGLX format, skipping file restructuring')
                return
            for i, this_dir in enumerate(os.listdir(join(self.session_path, 'raw_ephys_data', orig_dir))):
                shutil.move(join(self.session_path, 'raw_ephys_data', orig_dir, this_dir),
                            join(self.session_path, 'raw_ephys_data'))
            os.rmdir(join(self.session_path, 'raw_ephys_data', orig_dir))
            for i, this_path in enumerate(glob(join(self.session_path, 'raw_ephys_data', '*imec*'))):
                os.rename(this_path, join(self.session_path, 'raw_ephys_data', 'probe0' + this_path[-1]))
        return
    
    
    def nidq_synchronization(self):
        """
        Create synchronization file for the nidq

        """
        
        # Create synchronization file
        nidq_file = next(self.session_path.joinpath('raw_ephys_data').glob('*.nidq.*bin'))
        with open(nidq_file.with_suffix('.wiring.json'), 'w') as fp:
            json.dump(self.nidq_sync, fp, indent=1)
        
        for ap_file in self.session_path.joinpath('raw_ephys_data').rglob('*.ap.cbin'):
            with open(ap_file.with_suffix('.wiring.json'), 'w') as fp:
                json.dump(self.probe_sync, fp, indent=1)
        
        # Create nidq sync file        
        EphysSyncRegisterRaw(session_path=self.session_path, sync_collection='raw_ephys_data').run()
        
        return
                
    
    def preprocessing(self):
        """
        Run all the preprocessing steps before spike sorting.
        
        1. High pass filter
        2. Correct for the inter-sample shift in acquisition 
        3. Detect noisy channels and channels outside of the brain
        4. Remove channels out of the brain and interpolate over noisy channels
        5. - When single shank: perform destriping
           - When 4 shank: do common average referencing
        6. Apply notch filters if requested 

        Returns
        -------
        rec : SpikeInterface recording object
            The final preprocessed recording as a SpikeInterface object.

        """
        
        """
        # Decompress recording if it was compressed by a previous run
        if len(glob(join(self.probe_path, '*ap.cbin'))) > 0:

            # Recording is compressed by a previous run, decompress it before spike sorting
            cbin_path = glob(join(self.probe_path, '*ap.cbin'))[0]
            ch_path = glob(join(self.probe_path, '*ch'))[0]
            r = mtscomp.Reader(chunk_duration=1.)
            r.open(cbin_path, ch_path)
            r.tofile(cbin_path[:-4] + 'bin')
            r.close()
            
            # Remove compressed bin file after decompression
            if ((len(glob(join(self.session_path, 'raw_ephys_data', self.this_probe, '*ap.cbin'))) == 1)
                and (len(glob(join(self.session_path, 'raw_ephys_data', self.this_probe, '*ap.bin'))) == 1)):
                os.remove(glob(join(self.session_path, 'raw_ephys_data', self.this_probe, '*ap.cbin'))[0])
                self.ap_file = glob(join(self.session_path, 'raw_ephys_data', self.this_probe, '*ap.bin'))[0]
        """
        
        if len(glob(join(self.probe_path, '*ap.cbin'))) > 0:
            # Recording is already compressed by a previous run, loading in compressed data
            rec = si.read_cbin_ibl(self.probe_path)
        else: 
            rec = si.read_spikeglx(self.probe_path, stream_id=si.get_neo_streams('spikeglx', self.probe_path)[0][0])
                    
        # Apply high-pass filter
        print('\nApplying high-pass filter.. ')
        rec_filtered = si.highpass_filter(rec, ftype='bessel', dtype='float32')
                    
        # Correct for inter-sample phase shift
        print('Correcting for phase shift.. ')
        rec_shifted = si.phase_shift(rec_filtered)
        
        # Detect and interpolate over bad channels
        print('Detecting and interpolating over bad channels.. ')
        
        # Do common average referencing before detecting bad channels
        rec_comref = si.common_reference(rec_filtered)
        
        # Detect dead channels
        bad_channel_ids, all_channels = si.detect_bad_channels(rec_filtered, seed=42)
        prec_dead_ch = np.sum(all_channels == 'dead') / all_channels.shape[0]
        print(f'{np.sum(all_channels == "dead")} ({prec_dead_ch*100:.0f}%) dead channels')
        dead_channel_ids = rec_filtered.get_channel_ids()[all_channels == 'dead']
        prec_out_ch = np.sum(all_channels == 'out') / all_channels.shape[0]
        print(f'{np.sum(all_channels == "out")} ({prec_out_ch*100:.0f}%) channels outside of the brain')
        out_channel_ids = rec_filtered.get_channel_ids()[all_channels == 'out']
        
        # Detect noisy channels
        bad_channel_ids, all_channels = si.detect_bad_channels(rec_comref, method='mad',
                                                               std_mad_threshold=3, seed=42)
        prec_noise_ch = np.sum(all_channels == 'noise') / all_channels.shape[0]
        print(f'{np.sum(all_channels == "noise")} ({prec_noise_ch*100:.0f}%) noise channels')
        noisy_channel_ids = rec_comref.get_channel_ids()[all_channels == 'noise']
        
        # Remove channels that are outside of the brain
        rec_no_out = rec_shifted.remove_channels(remove_channel_ids=out_channel_ids)
        
        # Interpolate over bad channels          
        rec_interpolated = si.interpolate_bad_channels(rec_no_out, np.concatenate((
            dead_channel_ids, noisy_channel_ids)))
        
        # Destripe when there is one shank, CAR when there are four shanks
        if np.unique(rec_interpolated.get_property('group')).shape[0] == 1:
            print('Single shank recording; destriping')
            rec_processed = si.highpass_spatial_filter(rec_interpolated)
        else:
            print('Multi shank recording; common average reference')
            rec_processed = si.common_reference(rec_interpolated)
       
        # Plot spectral density
        print('Calculating power spectral density')
        data_chunk = si.get_random_data_chunks(rec_processed, num_chunks_per_segment=1,
                                               chunk_size=30000, seed=42)
        fig, ax = plt.subplots(figsize=(10, 7))
        for tr in data_chunk.T:
            p, f = ax.psd(tr, Fs=rec_processed.sampling_frequency, color="b")
        plt.savefig(join(self.probe_path, 'power spectral density.jpg'), dpi=600)
        
        # Apply notch filter 
        if isfile(join(self.probe_path, 'notch_filter.json')):
            
            # Load in notch filter settings
            with open(join(self.probe_path, 'notch_filter.json'), 'r') as openfile:
                notch_filter = json.load(openfile)
                
            # Apply filters
            rec_notch = rec_processed
            for freq, q in zip(notch_filter['FREQ'], notch_filter['Q']):
                print(f'Applying notch filter at {freq} Hz..')
                rec_notch = si.notch_filter(rec_notch, freq=freq, q=q)
                
            # Plot spectral density
            print('Calculating power spectral density')
            data_chunk = si.get_random_data_chunks(rec_notch, num_chunks_per_segment=1,
                                                   chunk_size=30000, seed=42)
            fig, ax = plt.subplots(figsize=(10, 7))
            for tr in data_chunk.T:
                p, f = ax.psd(tr, Fs=rec_processed.sampling_frequency, color="b")
            plt.savefig(join(self.probe_path, 'power spectral density after notch filter.jpg'), dpi=600)
            
            rec_final = rec_notch
        else:
            rec_final = rec_processed
            
        return rec_final
    
    
    def spikesorting(self, rec, probe_path):
        """
        Run spike sorting using SpikeInterface
        
        Returns
        -------
        sort : SpikeInterface object
            The spike sorted output 
            
        """
        
        # Run spike sorting
        try:
            sort = si.run_sorter(
                self.settings['SPIKE_SORTER'],
                rec,
                folder=join(probe_path, self.settings['SPIKE_SORTER']),
                verbose=True,
                docker_image=self.settings['USE_DOCKER'],
                **self.sorter_params)
        except Exception as err:
            
            # Log error to disk
            print(err)
            logf = open(os.path.join(probe_path, 'error_log.txt'), 'w')
            logf.write(str(err))
            logf.close()
            
            # Delete empty sorting directory
            if isdir(join(probe_path, self.settings['SPIKE_SORTER'] + self.settings['IDENTIFIER'])):
                shutil.rmtree(join(probe_path, self.settings['SPIKE_SORTER'] + self.settings['IDENTIFIER']))
            
            return None
        
        return sort
        
    
    def neuron_metrics(self, sort, rec):
        """
        Create sorting analyzer for manual curation in SpikeInterface and save to disk

        Parameters
        ----------
        sort : SpikeInterface object
            Result of the spike sorting
        rec : SpikeInterface object
            The preprocessed recording

        """
        
        if isdir(join(self.results_path, 'sorting')):
            return
        
        # Create a sorting analyzer and save to disk as folder
        sorting_analyzer = si.create_sorting_analyzer(
            sorting=sort,
            recording=rec,
            format='binary_folder',
            folder=join(self.results_path, 'sorting'),
            overwrite=True
            )           
        
        # Compute a bunch of attributes of the units
        sorting_analyzer.compute([
            'noise_levels',
            'correlograms',
            'isi_histograms',
            'random_spikes',
            'waveforms',
            #'principal_components',
            'templates',
            'template_similarity',
            'unit_locations',
            #'spike_locations',
            'spike_amplitudes',
            ])
        
        # Compute amplitude CV metrics
        #_ = si.compute_amplitude_cv_metrics(sorting_analyzer=sorting_analyzer)
        
        # Compute quality metrics
        _ = sorting_analyzer.compute('quality_metrics', metric_names=si.get_quality_metric_list())
        #_ = sorting_analyzer.compute('quality_metrics', metric_names=si.get_quality_pca_metric_list())        
                
        # Compute template metrics
        _ = si.compute_template_metrics(sorting_analyzer, include_multi_channel_metrics=True)
        
        # Compute drift metrics
        #_, _, _, = si.misc_metrics.compute_drift_metrics(sorting_analyzer)
                        
        return
        
        
    def raw_ephys_qc(self):
        """
        Calculate raw ephys QC metrics such as AP band RMS and LFP power per channel
        
        """
        
        # If there is no LF file (NP2 probes), generate it
        if len(glob(join(self.probe_path, '*lf.*bin'))) == 0:
            print('Generating LFP bin file (can take a while)')
            conv = NP2Converter(self.ap_file, compress=False)
            conv._process_NP21(assert_shanks=False)
            NP2_probe = True
        else:
            NP2_probe = False
                                    
        # Compute raw ephys QC metrics
        if not isfile(join(self.probe_path, '_iblqc_ephysSpectralDensityAP.power.npy')):
            task = ephysqc.EphysQC('', session_path=self.session_path, use_alyx=False)
            task.probe_path = self.probe_path
            task.run()                
            extract_rmsmap(self.ap_file, out_folder=self.probe_path, spectra=False)
        
        # If an LF bin file was generated, delete it (results in errors down the line)
        if NP2_probe and len(glob(join(self.probe_path, '*lf.*bin'))) == 1:
            os.remove(glob(join(self.probe_path, '*lf.*bin'))[0])
            os.remove(glob(join(self.probe_path, '*lf.*meta'))[0])
                
        return
    
    
    def convert_to_alf(self):
        """
        Convert Kilosort output to ALF files which are readable by the Open Neurophysiology 
        Environment (ONE), see more about this file format here: 
        https://int-brain-lab.github.io/iblenv/docs_external/alf_intro.html

        """        
        
        # Set the dat_file path correctly in params.py before conversion         
        with open(join(self.sorter_out_path, 'params.py'), 'r') as file:
            lines = file.readlines()
        lines[-1] = f"dat_path = '{self.ap_file}'\n"
        with open(join(self.sorter_out_path, 'params.py'), 'w') as file:
            file.writelines(lines)
            
        # Export as ALF files
        if not isdir(self.results_path):
            os.mkdir(self.results_path)
        ks2_to_alf(self.sorter_out_path, self.probe_path, self.results_path, bin_file=self.ap_file)
        
        # Delete phy waveforms (we won't use Phy)
        for phy_file in glob(join(self.results_path, '_phy_*')):
            os.remove(phy_file)
        
        # Move LFP power etc. to the alf folder
        qc_files = glob(join(self.probe_path, '_iblqc_*'))
        for ii, this_file in enumerate(qc_files):
            shutil.move(this_file, join(self.results_path, split(this_file)[1]))
        
        return
            
            
    def automatic_curation(self):
        """
        Add unit level QC from Kilosort and IBL to the quality metrics so that they show up
        during manual curation

        Parameters
        ----------
        results_path : Path
            Path to where the final results of the spike sorting will be saved.
        sorter_out_path : Path
            Path to the output of Kilosort.

        Returns
        -------
        None.

        """
        
        # Get kilosort good indication 
        ks_metric = pd.read_csv(join(self.sorter_out_path, 'cluster_KSLabel.tsv'), sep='\t')
        
        if hasattr(si, 'auto_label_units'):
            
            # Load in recording
            sorting_analyzer = si.load_sorting_analyzer(join(self.results_path, 'sorting'))
                     
            # Apply the sua/mua model
            ml_labels = si.auto_label_units(
                sorting_analyzer = sorting_analyzer,
                repo_id = 'AnoushkaJain3/sua_mua_classifier_lightweight',
                #repo_id = 'AnoushkaJain3/sua_mua_classifier',
                trust_model=True,
            )
                        
        else:
            print('Could not run machine learning model for unit curation\n'
                  'Update SpikeInterface to version >= 0.102.0\n'
                  'You might have to install SpikeInterface from source')
            ml_labels = np.zeros(ks_metric.shape[0])
        
            
        # Calculate IBL neuron level QC
        print('Calculating IBL neuron-level quality metrics..')
        spikes, clusters, channels = load_neural_data(self.session_path,
                                                      self.this_probe,
                                                      histology=False, only_good=False)
        df_units, rec_qc = spike_sorting_metrics(spikes['times'], spikes['clusters'],
                                                 spikes['amps'], spikes['depths'])
        
        # Add to quality metrics
        qc_metrics = pd.read_csv(join(self.results_path, 'sorting', 'extensions', 'quality_metrics',
                                      'metrics.csv'), index_col=0)
        qc_metrics['KS_label'] = (ks_metric['KSLabel'] == 'good').astype(int)
        qc_metrics.insert(0, 'KS_label', qc_metrics.pop('KS_label'))
        qc_metrics['IBL_label'] = df_units['label']
        qc_metrics.insert(0, 'IBL_label', qc_metrics.pop('IBL_label'))
        qc_metrics['ML_label'] = (ml_labels['prediction'] == 'sua').astype(int)
        qc_metrics.insert(0, 'ML_label', qc_metrics.pop('ML_label'))
        
        # Save to disk
        qc_metrics.to_csv(join(
            self.results_path, 'sorting', 'extensions', 'quality_metrics', 'metrics.csv'))
        np.save(join(self.results_path, 'clusters.IBLLabel.npy'), qc_metrics['IBL_label'])
        np.save(join(self.results_path, 'clusters.KSLabel.npy'), qc_metrics['KS_label'])
        np.save(join(self.results_path, 'clusters.MLLabel.npy'), qc_metrics['ML_label'])
        if isfile(join(self.results_path, 'cluster_KSLabel.tsv')):
            os.remove(join(self.results_path, 'cluster_KSLabel.tsv'))
        
        # Copy quality metrics to output folder
        shutil.copy(join(self.results_path, 'sorting', 'extensions', 'quality_metrics', 'metrics.csv'),
                    join(self.results_path, 'clusters.metrics.csv'))
        
        return
        
        
    def probe_synchronization(self):
        """
        Synchronize spikes of this probe to the nidq base station

        """
       
        # Create probe sync file
        task = EphysSyncPulses(session_path=self.session_path, sync='nidq', pname=self.this_probe,
                               sync_ext='bin', sync_namespace='spikeglx',
                               sync_collection='raw_ephys_data',
                               device_collection='raw_ephys_data')
        task.run()
        task = EphysPulses(session_path=self.session_path, pname=self.this_probe,
                           sync_collection='raw_ephys_data',
                           device_collection='raw_ephys_data')
        task.run()
        
        # Synchronize spike sorting to nidq clock
        sync_spike_sorting(self.ap_file, self.results_path)
        
        # Extract digital sync timestamps
        sync_times = np.load(join(self.session_path, 'raw_ephys_data', '_spikeglx_sync.times.npy'))
        sync_polarities = np.load(join(self.session_path, 'raw_ephys_data', '_spikeglx_sync.polarities.npy'))
        sync_channels = np.load(join(self.session_path, 'raw_ephys_data', '_spikeglx_sync.channels.npy'))
        for ii, ch_name in enumerate(self.nidq_sync['SYNC_WIRING_DIGITAL'].keys()):
            if ch_name == 'imec_sync':
                continue
            nidq_pulses = sync_times[(sync_channels == int(ch_name[-1])) & (sync_polarities == 1)]
            np.save(join(self.session_path, self.nidq_sync['SYNC_WIRING_DIGITAL'][ch_name] + '.times.npy'),
                    nidq_pulses)
        return
        
    
    def compress_raw_data(self):
        """
        Compress the raw bin file using mtscomp compression developed by IBL        
        After compression the file will be a .cbin file instead of .bin
        
        """
        
        # Load in recording 
        if len(glob(join(self.probe_path, '*.cbin'))) > 0:
            # Recording is already compressed by a previous run
            return
        else:
            rec = si.read_spikeglx(self.probe_path, stream_id=si.get_neo_streams('spikeglx', self.probe_path)[0][0])
        
        if self.settings['COMPRESS_RAW_DATA']:
            if len(glob(join(self.session_path, 'raw_ephys_data', self.this_probe, '*ap.cbin'))) == 0:
                print('Compressing raw binary file')
                mtscomp.compress(self.ap_file, str(self.ap_file)[:-3] + 'cbin', str(self.ap_file)[:-3] + 'ch',
                                 sample_rate=rec.get_sampling_frequency(),
                                 n_channels=rec.get_num_channels() + 1,
                                 dtype=rec.get_dtype())
                
            # Delete original raw data
            if ((len(glob(join(self.session_path, 'raw_ephys_data', self.this_probe, '*ap.cbin'))) == 1)
                and (len(glob(join(self.session_path, 'raw_ephys_data', self.this_probe, '*ap.bin'))) == 1)):
                try:
                    os.remove(glob(join(self.session_path, 'raw_ephys_data', self.this_probe, '*ap.bin'))[0])
                except:
                    print('Could not remove uncompressed ap bin file, delete manually')
                    return
        return
        
    
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
    
    # Load in sorting analyzer from disk
    sorting_analyzer = si.load_sorting_analyzer(join(results_path, 'sorting'))
    
    # Launch manual curation GUI            
    _ = si.plot_sorting_summary(sorting_analyzer=sorting_analyzer, curation=True,
                                backend='spikeinterface_gui')
    
    # Extract manual curation labels and save in results folder
    if isfile(join(results_path, 'sorting', 'spikeinterface_gui', 'curation_data.json')):
        with open(join(results_path, 'sorting', 'spikeinterface_gui', 'curation_data.json')) as f:
            label_dict = json.load(f)
        if isfile(join(results_path, 'clusters.manualLabels.npy')):
            manual_labels = np.load(join(results_path, 'clusters.manualLabels.npy'))
        else:
            manual_labels = np.array(['no label'] * sorting_analyzer.unit_ids.shape[0])
        for this_unit in label_dict['manual_labels']:
            manual_labels[sorting_analyzer.unit_ids == this_unit['unit_id']] = this_unit['quality']
        np.save(join(results_path, 'clusters.manualLabels.npy'), manual_labels)
    
    return       
            

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
