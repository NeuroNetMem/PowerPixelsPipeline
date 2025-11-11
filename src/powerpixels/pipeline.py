#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Written by Guido Meijer

"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
from os.path import join, isfile, isdir
from pathlib import Path
import shutil
from glob import glob
import json
from scipy.signal import welch, find_peaks
from .utils import load_neural_data
import bombcell as bc
import spikeinterface.full as si
import mtscomp
from brainbox.metrics.single_units import spike_sorting_metrics
from brainbox.metrics.single_units import METRICS_PARAMS as ibl_qc_default_params
from ibllib.ephys.spikes import sync_spike_sorting
from ibllib.pipes.ephys_tasks import EphysSyncPulses, EphysSyncRegisterRaw, EphysPulses


class Pipeline:
    
    def __init__(self):
        
        project_root = Path(__file__).parent.parent.parent
        config_dir = project_root / 'config'
        settings_file = config_dir / 'settings.json'
        bombcell_file = config_dir / 'bombcell_params.json'
        ibl_qc_file = config_dir / 'ibl_qc_params.json'

        # Check if the config file exists
        if not settings_file.is_file():
            raise FileNotFoundError(
                f'Configuration file not found at {settings_file}\n'
                'Please run "powerpixels-setup" to create the default files.'
            )
        
        # Load in config files
        with open(settings_file, 'r') as openfile:
            self.settings = json.load(openfile)
        if self.settings['SINGLE_SHANK'] not in ['car_local', 'car_global', 'destripe']:
            raise ValueError('SINGLE_SHANK should be set to "car_global", "car_local" or "destripe"')
        if self.settings['MULTI_SHANK'] not in ['car_local', 'car_global', 'destripe']:
            raise ValueError('MULTI_SHANK should be set to "car_global", "car_local" or "destripe"')
        if self.settings['USE_NIDAQ']:
            with open(config_dir / 'wiring' / 'nidq.wiring.json', 'r') as openfile:
                self.nidq_sync = json.load(openfile)
        with open(config_dir /'wiring' / '3B.wiring.json', 'r') as openfile:
            self.probe_sync = json.load(openfile)
        with open(bombcell_file, 'r') as openfile:
            self.bombcell_params = json.load(openfile)
        with open(ibl_qc_file, 'r') as openfile:
            self.ibl_qc_params = json.load(openfile)
            
        
        # Check spike sorter
        if self.settings['SPIKE_SORTER'][:8] != 'kilosort':
            print('\n --- WARNING: use Kilosort (any version) for full functionality of the pipeline --- \n')
        
        # Initialize spikeinterface parallel processing
        si.set_global_job_kwargs(n_jobs=self.settings['N_CORES'], progress_bar=True)
        
        # Load in spike sorting parameters
        if (config_dir / 'sorter_params' / f'{self.settings["SPIKE_SORTER"]}_params.json').is_file():
            with open(config_dir / 'sorter_params'
                      / f'{self.settings["SPIKE_SORTER"]}_params.json', 'r') as openfile:
                self.sorter_params = json.load(openfile)
        else:
            print('Did not find spike sorter parameter file, loading defaults..')
            self.sorter_params = si.get_default_sorter_params(self.settings['SPIKE_SORTER'])
            

    def detect_data_format(self):
                
        if len(list((self.session_path / 'raw_ephys_data').rglob('continuous.dat'))) > 0:
            self.data_format = 'openephys'
            print('OpenEphys recording detected')
        elif len(list((self.session_path / 'raw_ephys_data').rglob('*ap.meta'))) > 0:
            self.data_format = 'spikeglx'
            print('SpikeGLX recording detected')
        else:
            print(f'Could not detect recording in {self.session_path / "raw_ephys_data"}')
            return


    def restructure_files(self):
        """
        Restructure the raw data files
    
        """
        
        # Detect data format if necessary
        if not hasattr(self, 'data_format'):
            self.detect_data_format()
        
        # Restructure SpikeGLX recording folder structure
        if ((self.data_format == 'spikeglx')
            and len([i for i in os.listdir(self.session_path / 'raw_ephys_data') if i[:5] == 'probe']) == 0):
            orig_dir = os.listdir(self.session_path / 'raw_ephys_data')[0]
            for i, this_dir in enumerate(os.listdir(self.session_path / 'raw_ephys_data' / orig_dir)):
                shutil.move(self.session_path / 'raw_ephys_data' / orig_dir / this_dir,
                            self.session_path / 'raw_ephys_data')
            os.rmdir(self.session_path / 'raw_ephys_data' / orig_dir)
            for i, this_path in enumerate((self.session_path / 'raw_ephys_data').glob('*imec*')):
                this_path.rename(self.session_path / 'raw_ephys_data' / ('probe0' + str(this_path)[-1]))
                
        return

        
    def extract_sync_pulses(self):
        """
        Extract the synchronization pulses coming into the BNC breakout box of the NIDAQ
        This also prepares for the synchronization of the spike times to the NIDAQ clock later on
        
        """
        
        
        # Detect data format if necessary
        if not hasattr(self, 'data_format'):
            self.detect_data_format()

        # Extract sync pulses from spikeGLX recordings
        if self.data_format == 'spikeglx':
            # Set path to nidq file
            if len(list((self.session_path / 'raw_ephys_data').glob('*.nidq.*bin'))) == 1:
                self.nidq_file = list((self.session_path / 'raw_ephys_data').glob('*.nidq.*bin'))[0]        
    
            # Create synchronization file
            with open(self.nidq_file.with_suffix('.wiring.json'), 'w') as fp:
                json.dump(self.nidq_sync, fp, indent=1)
            
            # Create nidq sync file        
            EphysSyncRegisterRaw(session_path=self.session_path, sync_collection='raw_ephys_data').run()            
        
            # Extract sync pulses    
            task = EphysSyncPulses(session_path=self.session_path, sync='nidq',
                                   sync_ext='bin', sync_namespace='spikeglx',
                                   sync_collection='raw_ephys_data',
                                   device_collection='raw_ephys_data')
            task.run()
        
            # Extract digital sync timestamps
            sync_times = np.load(join(self.session_path, 'raw_ephys_data', '_spikeglx_sync.times.npy'))
            sync_polarities = np.load(join(self.session_path, 'raw_ephys_data', '_spikeglx_sync.polarities.npy'))
            sync_channels = np.load(join(self.session_path, 'raw_ephys_data', '_spikeglx_sync.channels.npy'))
            for ii, ch_name in enumerate(self.nidq_sync['SYNC_WIRING_DIGITAL'].keys()):
                if self.nidq_sync['SYNC_WIRING_DIGITAL'][ch_name] == 'imec_sync':
                    continue
                nidq_pulses = sync_times[(sync_channels == int(ch_name[-1])) & (sync_polarities == 1)]
                np.save(join(self.session_path, self.nidq_sync['SYNC_WIRING_DIGITAL'][ch_name] + '.times.npy'),
                        nidq_pulses)
                
        elif self.data_format == 'openephys':
            # TO DO
            print('TO DO')


    def set_probe_paths(self, probe):
        
        # Check if session path is defined
        if not hasattr(self, 'session_path'):
            raise ValueError('Define pp.session_path first!')
            
        # Detect data format if necessary
        if not hasattr(self, 'data_format'):
            self.detect_data_format()
        
        # Set the current probe and results directory
        self.this_probe = probe
        self.results_path = self.session_path / (self.this_probe + self.settings['IDENTIFIER'])
        self.sorter_path = (self.session_path
                            / (self.settings['SPIKE_SORTER'] + self.settings['IDENTIFIER'])
                            / self.this_probe)
        
        # Set SpikeGLX specific paths
        if self.data_format == 'spikeglx':
            self.probe_path = self.session_path / 'raw_ephys_data' / probe
            self.ap_file = list(self.probe_path.glob('*ap.*bin'))[0]
            if len(list((self.session_path / 'raw_ephys_data' / probe).glob('*ap.meta'))) == 1:
                self.meta_file = list(self.probe_path.glob('*ap.meta'))[0]
        
        # Set OpenEphys specific paths
        elif self.data_format == 'openephys':
            for ap_path in (self.session_path / 'raw_ephys_data').rglob('continuous*'):
                if ((ap_path.parent.parts[-1].split('.')[-1] == f'{self.this_probe}-AP')
                    or (ap_path.parent.parts[-1].split('.')[-1] == f'{self.this_probe}')):
                    self.ap_file = ap_path
            self.meta_file = list((self.session_path / 'raw_ephys_data').rglob('structure.oebin'))[0]

        return
    
    
    def decompress(self):
        """
        Decompress data before running the pipeline, some elements like Bombcell need to have
        uncompressed raw data as input
        
        """

        # If data is not compressed stop this process
        if (self.ap_file.suffix == '.bin') or (self.ap_file.suffix == '.dat'):
            return
        
        if self.ap_file.suffix == '.cbin':
            
            # Recording is compressed by a previous run, decompress it before spike sorting
            ch_path = list(self.probe_path.glob('*ch'))[0]
            r = mtscomp.Reader(chunk_duration=1.)
            r.open(self.ap_file, ch_path)
            r.tofile(self.ap_file.parent / (self.ap_file.stem + '.bin'))
            r.close()
            
            # Remove compressed bin file after decompression
            if self.ap_file.is_file() and Path(self.ap_file.stem + '.bin').is_file():
                os.remove(self.ap_file)
            self.ap_file = self.ap_file.parent / (self.ap_file.stem + '.bin')
            
        elif self.ap_file.suffix == '.zarr':
            
            # Decompress zarr file
            comp_rec = si.load_extractor(self.ap_file)
            si.write_binary_recording(
                comp_rec, file_paths=[self.ap_file.parent / (self.ap_file.stem + '.dat')])

            # Remove compressed bin file after decompression
            if self.ap_file.is_dir() and Path(self.ap_file.stem + '.dat').is_file():
                shutil.rmtree(self.ap_file)
            self.ap_file = self.ap_file.parent / (self.ap_file.stem + '.dat')
            
        return
            
    
    def load_raw_binary(self):
        
        # Detect data format if necessary
        if not hasattr(self, 'data_format'):
            self.detect_data_format()
            
        # Load in raw data
        if self.data_format == 'spikeglx':
            if len(glob(join(self.probe_path, '*ap.cbin'))) > 0:
                rec = si.read_cbin_ibl(self.probe_path)
            else: 
                rec = si.read_spikeglx(
                    self.probe_path,
                    stream_id=si.get_neo_streams('spikeglx', self.probe_path)[0][0])
        
        elif self.data_format == 'openephys':
            stream_names, _ = si.read_openephys(
                self.session_path, stream_id='0').get_streams(self.session_path)
            these_streams = [i for i in stream_names if self.this_probe in i]
            if len(these_streams) == 1:  # NP2 recording
                rec_stream = these_streams[0]
            elif len(these_streams) == 2:  # NP1 recording
                rec_stream = [i for i in stream_names if self.this_probe + '-AP' in i][0]
            rec = si.read_openephys(self.session_path, stream_name=rec_stream)
            
        return rec
    
    
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
        
        # Load in raw data
        rec = self.load_raw_binary()
                    
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
        
        # Perform spatial filtering
        if np.unique(rec_interpolated.get_property('group')).shape[0] == 1:
            print(f'Single shank recording detected, performing: {self.settings["SINGLE_SHANK"]}')
            if self.settings['SINGLE_SHANK'] == 'car_global':
                rec_processed = si.common_reference(rec_interpolated)
            elif self.settings['SINGLE_SHANK'] == 'car_local':
                rec_processed = si.common_reference(rec_interpolated, reference='local')
            elif self.settings['SINGLE_SHANK'] == 'destripe':
                rec_processed = si.highpass_spatial_filter(rec_interpolated)
        else:
            print(f'Multi shank recording detected, performing: {self.settings["MULTI_SHANK"]}')
            if self.settings['MULTI_SHANK'] == 'car_global':
                rec_processed = si.common_reference(rec_interpolated)
            elif self.settings['MULTI_SHANK'] == 'car_local':
                rec_processed = si.common_reference(rec_interpolated, reference='local')
            elif self.settings['MULTI_SHANK'] == 'destripe':
                rec_processed = si.highpass_spatial_filter(rec_interpolated)
        
        
        # Calculate power spectral density
        data_chunk = si.get_random_data_chunks(rec_processed, num_chunks_per_segment=1,
                                               chunk_size=30000, seed=42)
        all_power = []
        for tr in data_chunk.T:
            f, p = welch(tr, fs=rec_processed.sampling_frequency)
            all_power.append(p)
        mean_power = np.mean(np.vstack(all_power), axis=0)
        
        # Detect peaks
        peak_inds, peak_props = find_peaks(mean_power, threshold=self.settings['PEAK_THRESHOLD'])
        peak_freqs = f[peak_inds]
        peak_freqs = peak_freqs[peak_freqs > 2000]  # only select high frequency peaks
        print(f'Detected {peak_inds.shape[0]} peaks in the power spectrum')
        
        # Plot
        fig, ax = plt.subplots(figsize=(10, 7))
        ax.plot(f, mean_power, zorder=0)
        for peak_ind in peak_inds:
            ax.scatter(f[peak_ind], mean_power[peak_ind], marker='x', color='r', zorder=1)
        ax.set(ylabel='Power spectral density', xlabel='Frequency (Hz)')
        plt.tight_layout()
        plt.savefig(self.session_path / 'raw_ephys_data'
                    / f'{self.this_probe} power spectral density.jpg', dpi=600)
        
        # Apply notch filter 
        if peak_freqs.shape[0] > 0:
            rec_notch = rec_processed
            for freq in peak_freqs:
                print(f'Applying notch filter at {int(freq)} Hz..')
                rec_notch = si.notch_filter(rec_notch, freq=freq, q=10)
                
            # Calculate power spectral density
            data_chunk = si.get_random_data_chunks(rec_notch, num_chunks_per_segment=1,
                                                   chunk_size=30000, seed=42)
            all_power = []
            for tr in data_chunk.T:
                f, p = welch(tr, fs=rec_notch.sampling_frequency)
                all_power.append(p)
            mean_power = np.mean(np.vstack(all_power), axis=0)
            
            # Plot power spectrum after notch filters
            fig, ax = plt.subplots(figsize=(10, 7))
            ax.plot(f, mean_power)
            ax.set(ylabel='Power spectral density', xlabel='Frequency (Hz)')
            plt.tight_layout()
            plt.savefig(self.session_path / 'raw_ephys_data'
                        / f'{self.this_probe} power spectral density filtered.jpg', dpi=600)
            
            rec_final = rec_notch
        else:
            rec_final = rec_processed
        
        return rec_final
    
    
    def spikesorting(self, rec):
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
                folder=self.sorter_path,
                verbose=True,
                docker_image=self.settings['USE_DOCKER'],
                **self.sorter_params)
        except Exception as err:
            
            # Log error to disk
            print(err)
            logf = open(self.session_path / f'{self.this_probe}_error_log.txt', 'w')
            logf.write(str(err))
            logf.close()
            
            # Delete empty sorting directory
            if self.sorter_path.is_dir():
                shutil.rmtree(self.sorter_path)
            
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
            'templates',
            'template_similarity',
            'unit_locations',
            'spike_amplitudes',
            ])
                
        # Compute quality metrics
        _ = sorting_analyzer.compute('quality_metrics', metric_names=si.get_quality_metric_list())    
                
        # Compute template metrics
        _ = si.compute_template_metrics(sorting_analyzer, include_multi_channel_metrics=True)
                                
        return
        

    def export_data(self, rec):
        """
        Calculate raw ephys QC metrics such as AP band RMS and LFP power per channel
        Export spike sorted data in a format the alignment GUI can read in
        
        """
        
        # Load in sorting output
        sorting_analyzer = si.load_sorting_analyzer(self.results_path / 'sorting')
        
        # Load in raw LFP 
        rec_lfp = si.bandpass_filter(rec, freq_min=1, freq_max=300)
        
        # Export data to temporary folder
        si.export_to_ibl_gui(
            sorting_analyzer=sorting_analyzer,
            output_folder=self.results_path / 'exported_data',
            lfp_recording=rec_lfp,
            n_jobs=-1
        )
        
        # Copy the extracted data to the parent folder
        for file_path in (self.results_path / 'exported_data').iterdir():
            shutil.move(file_path, self.results_path)
        (self.results_path / 'exported_data').rmdir()
        
        # Export to NWB
        if self.settings['NWB_EXPORT']:
            
            print('/nExporting to NWB format/n')
            from neuroconv.datainterfaces import SpikeGLXRecordingInterface
            from neuroconv.datainterfaces import OpenEphysRecordingInterface
            from neuroconv.tools.spikeinterface import write_sorting_analyzer_to_nwbfile
            
            if self.data_format == 'openephys':
                interface = OpenEphysRecordingInterface(folder_path=self.session_path,
                                                        stream_name=rec.stream_name)   
            elif self.data_format == 'spikeglx':
                interface = SpikeGLXRecordingInterface(folder_path=self.session_path,
                                                        stream_name=rec.stream_name)   
            metadata = interface.get_metadata()
            nwbfile = interface.create_nwbfile(metadata=metadata)
            if not (self.results_path / 'NWB').is_dir():
                (self.results_path / 'NWB').mkdir()            
            write_sorting_analyzer_to_nwbfile(
                sorting_analyzer=sorting_analyzer,
                nwbfile_path=self.results_path / 'NWB' / 'sorting_output.nwb',
                nwbfile=nwbfile,
                metadata=metadata,
                overwrite=True,
                recording=rec
                )

    
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
        ks_metric = pd.read_csv(join(self.sorter_path / 'sorter_output', 'cluster_KSLabel.tsv'),
                                sep='\t')
        
        # Run Bombcell
        print('\nRunning Bombcell..', end=' ')
        if self.settings['SPIKE_SORTER'][-3:] == '2_5':
            kilosort_version = 2
        else:
            kilosort_version = int(self.settings['SPIKE_SORTER'][-1])
        param = bc.get_default_parameters(self.sorter_path / 'sorter_output', 
                                          raw_file=self.ap_file,
                                          meta_file=self.meta_file,
                                          kilosort_version=kilosort_version)
        param['plotGlobal'] = False
        param['verbose'] = False
        param.update(self.bombcell_params)
        quality_metrics, param, unit_type, unit_type_string = bc.run_bombcell(
            self.sorter_path / 'sorter_output', self.results_path / 'bombcell', param)
        print('Done')
        
        # Run UnitRefine
        print('\nRunning UnitRefine model..', end=' ')
       
        # Load in recording
        sorting_analyzer = si.load_sorting_analyzer(self.results_path / 'sorting')
                 
        # Apply the sua/mua model
        ml_labels = si.auto_label_units(
            sorting_analyzer = sorting_analyzer,
            repo_id = 'AnoushkaJain3/sua_mua_classifier_lightweight',
            trust_model=True,
        )
        print('Done')
                
        # Calculate IBL neuron level QC
        print('\nCalculating IBL neuron-level quality metrics..', end=' ')
        spikes, clusters, channels = load_neural_data(self.session_path, self.this_probe)
        ibl_qc_default_params.update(self.ibl_qc_params)
        df_units, rec_qc = spike_sorting_metrics(spikes['times'], spikes['clusters'],
                                                 spikes['amps'], spikes['depths'],
                                                 params=ibl_qc_default_params)
        print('Done')
        
        # Print results
        n_units = unit_type_string.shape[0]
        bc_perc = np.round((np.sum(unit_type_string == "GOOD") / n_units) * 100, 1)
        ur_perc = np.round((np.sum(ml_labels['prediction'] == 'sua') / n_units) * 100, 1)
        ibl_perc = np.round((np.sum(df_units['label'] == 1) / n_units) * 100, 1)
        print('\n---------------------------------------------------------\n',
              'Automatic curation results',
              '\n---------------------------------------------------------',
              f'\nBombcell: {np.sum(unit_type_string == "GOOD")} of {n_units} units classified as good ({bc_perc}%)',
              f'\nUnitRefine: {np.sum(ml_labels["prediction"] == "sua")} of {n_units} units classified as good ({ur_perc}%)',
              f'\nIBL: {np.sum(df_units["label"] == 1)} of {n_units} units classified as good ({ibl_perc}%)\n')
        
        
        # Add to quality metrics
        qc_metrics = pd.read_csv(join(self.results_path, 'sorting', 'extensions', 'quality_metrics',
                                      'metrics.csv'), index_col=0)
        qc_metrics['Kilosort'] = (ks_metric['KSLabel'] == 'good').astype(int)
        qc_metrics.insert(0, 'Kilosort', qc_metrics.pop('Kilosort'))
        qc_metrics['IBL'] = df_units['label']
        qc_metrics.insert(0, 'IBL', qc_metrics.pop('IBL'))
        qc_metrics['UnitRefine'] = (ml_labels['prediction'] == 'sua').astype(int)
        qc_metrics.insert(0, 'UnitRefine', qc_metrics.pop('UnitRefine'))
        qc_metrics['Bombcell'] = unit_type.astype(int)
        qc_metrics.insert(0, 'Bombcell', qc_metrics.pop('Bombcell'))
        
        # Save to disk
        qc_metrics.to_csv(join(
            self.results_path, 'sorting', 'extensions', 'quality_metrics', 'metrics.csv'))
        np.save(join(self.results_path, 'clusters.iblLabels.npy'), qc_metrics['IBL'])
        np.save(join(self.results_path, 'clusters.kilosortLabels.npy'), qc_metrics['Kilosort'])
        np.save(join(self.results_path, 'clusters.unitrefineLabels.npy'), qc_metrics['UnitRefine'])
        np.save(join(self.results_path, 'clusters.bombcellLabels.npy'), qc_metrics['Bombcell'])
        if isfile(join(self.results_path, 'cluster_KSLabel.tsv')):
            os.remove(join(self.results_path, 'cluster_KSLabel.tsv'))
        
        # Copy quality metrics to output folder
        shutil.copy(join(self.results_path, 'sorting', 'extensions', 'quality_metrics', 'metrics.csv'),
                    join(self.results_path, 'clusters.metrics.csv'))
        
        return
    
        
    def probe_synchronization(self):
        """
        Synchronize spikes of this probe to the NIDAQ base station clock

        """
       
        # Create probe wiring file
        with open(self.ap_file.with_suffix('.wiring.json'), 'w') as fp:
            json.dump(self.probe_sync, fp, indent=1)
       
        # Create probe sync file
        task = EphysPulses(session_path=self.session_path, pname=self.this_probe,
                           sync_collection='raw_ephys_data',
                           device_collection='raw_ephys_data')
        task.run()
        
        # Synchronize spike sorting to nidq clock
        sync_spike_sorting(self.ap_file, self.results_path)
        

        return
        
    
    def compress_raw_data(self):
        """
        Compress raw data using either zarr or mtscomp compression
        
        """
        
        # Load in raw binary
        rec = self.load_raw_binary()
        
        if self.settings['COMPRESSION'] == 'zarr':
            
            if self.ap_file.suffix == '.zarr':
                # Recording is already compressed by a previous run
                return
            
            # Compress raw data to zarr folder
            rec.save(folder=self.ap_file.parent.parent / 'continuous.zarr', format='zarr')
            
            # Move zarr folder to right place
            shutil.move(self.ap_file.parent.parent / 'continuous.zarr',
                        self.ap_file.parent / 'continuous.zarr')
            
            # Delete original .dat file
            os.remove(self.ap_file)
            
            # Update reference to ap file
            self.ap_file = self.ap_file.parent / (str(self.ap_file.stem) + '.zarr')
            
        elif self.settings['COMPRESSION'] == 'mtscomp':
        
            if self.ap_file.suffix == '.cbin':
                # Recording is already compressed by a previous run
                return
            
            # Compress
            mtscomp.compress(self.ap_file, str(self.ap_file)[:-3] + 'cbin',
                             str(self.ap_file)[:-3] + 'ch',
                             sample_rate=rec.get_sampling_frequency(),
                             n_channels=rec.get_num_channels() + 1,
                             dtype=rec.get_dtype())
            
            # Delete original bin file
            os.remove(self.ap_file)
            
            # Update reference to ap file
            self.ap_file = self.ap_file.parent / (str(self.ap_file.stem) + '.cbin')

        return
        
