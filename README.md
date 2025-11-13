![GitHub License](https://img.shields.io/github/license/NeuroNetMem/PowerPixelsPipeline)
# Power Pixels: A turnkey pipeline for processing of Neuropixel recordings ⚡ 
<img src="https://github.com/user-attachments/assets/37d9003a-788a-43d2-a5a6-ff4d7a29f780" alt="PowerPixels logo" width="30%" align="right" vspace="20"/>


📄 Please cite the [bioRxiv preprint](https://doi.org/10.1101/2025.06.27.661890) if you use the pipeline 📄 ⭐ And star this repository! ⭐

The Power Pixels pipeline combines several packages and workflows into one end-to-end pipeline. It supports Neuropixel 1.0 and 2.0 probes recorded on a National Instruments system using SpikeGLX or OpenEphys. 

This pipeline is all about combining existing modules and pipelines into one, which is especially useful for people who are just starting out doing Neuropixel recordings and maybe have heard of all these tools but need some help getting them all integrated. Power Pixels relies on these amazing open-source projects:
- [SpikeInterface](https://spikeinterface.readthedocs.io)
- [ibllib](https://github.com/int-brain-lab/ibllib)
- [Kilosort](https://github.com/MouseLand/Kilosort)
- [AP_histology](https://github.com/petersaj/AP_histology)
- [Universal Probe Finder](https://github.com/JorritMontijn/UniversalProbeFinder)
- [Bombcell](https://github.com/Julie-Fabre/bombcell)
- [UnitRefine](https://huggingface.co/SpikeInterface/UnitRefine_sua_mua_classifier)
- [neuroconv](https://neuroconv.readthedocs.io/en/stable/)

## Description of the pipeline elements

![pipeline process](https://github.com/user-attachments/assets/1a6b70e7-6f5f-4c3f-83d8-1de4c1d5ccce)

The pipeline contains the following elements:
- **High-frequency noise**: high-frequency noise in specific frequency bands is automatically filtered out using notch filters targeted to detected peaks in the power spectrum. 
- **Phase shift correction**: channels on a Neuropixel probe are not recorded simultaneously, there is a small delay in the order of microseconds between the acquisition of a block of channels. Correcting for this small delay greatly improves artifact removal.
- **Remove bad channels**: bad channels are detected by looking at both coherence with other channels and PSD power in the high-frequency range, then they are interpolated using neighboring channels. Channels outside of the brain are removed.
- **Artifact removal**: the user can decide whether to apply common average referencing, local average referencing, or destriping to remove electrical artifacts and noise.
- **Spike sorting**: a spike sorting algorithm is used to detect spikes and sort them into units. SpikeInterface supports many [spike sorters](https://spikeinterface.readthedocs.io/en/latest/modules/sorters.html#supported-spike-sorters) out of the box (recommended: Kilosort).
- **Automatic classification of single neurons**: The pipeline runs three algorithms for automatic classification of good single units: Bombcell, UnitRefine and the IBL quality criteria.
- **Synchronization**: each Neuropixel probe and the BNC breakout box has their own clock. This means one has to synchronize the spike times between the probes (if you use more than one) and the synchronization channels which carry timestamps of events (for example: behavioral events or pulses from a camera).
- **Compression**: the raw binary file is compressed using *zarr* or *mtscomp* compression which results in a 2-3x reduction in file size.
- **Histological tracing**: the fluorescent tracks of the probes are traced using AP_histology or Universal Probe Finder.
- **Ephys-histology alignment**: the brain regions along the probe, inferred from the tracing, are aligned to electrophysiological features.

## Installation

It is recommended to install Power Pixels in an Anaconda or Miniforge environment.
1. Install [Anaconda](https://www.anaconda.com/) or [Miniforge](https://github.com/conda-forge/miniforge) - Miniforge is the recommended option
2. Open the Anaconda or Miniforge prompt
3. Create a new environment by typing `conda create -n powerpixels python=3.10 git` (use `mamba` instead of `conda` when using Miniforge)
4. Navigate to the location on your computer you want the repository to be and clone the repository by typing `git clone https://github.com/NeuroNetMem/PowerPixelsPipeline`
5. Navigating to the repository directory you just cloned in your console (`cd PowerPixelPipeline`) and install PowerPixels with the command `pip install -e .`
6. You can now activate the environment by typing `conda activate powerpixels` (or `mamba`)
7. Install `iblapps` by cloning the repository `git clone https://github.com/int-brain-lab/iblapps` and installing it with the command `pip install -e iblapps`

### Spike sorting
To install a spike sorter there are two options: (1) directly install Kilosort4 in the python environment, or (2) use Docker to run the spike sorter in a container. Note: if you want to use a MATLAB-based spike sorter (like Kilosort 2.5) you will have to pick option 2. 

_Option 1: local installation of Kilosort4_

Kilosort4 is already installed with PowerPixels, what you'll have to do now is make sure it uses the GPU.
1. In a terminal window activate the `spikeinterface` environment
2. Remove the CPU version of PyTorch by typing `pip uninstall torch`
3. Install the GPU version of PyTorch (for CUDA 11.8) with `pip3 install torch --index-url https://download.pytorch.org/whl/cu118`.
4. Check whether the GPU is used by opening a python terminal with the command `ipython` and typing `import torch; torch.cuda.is_available()`. If the result is True you are good to go. If it's False you need to make sure PyTorch can use CUDA, good luck! (ask AI to help you out)  

_Option 2: run spike sorter in Docker_
1. Install [Docker Desktop](https://www.docker.com/products/docker-desktop/)
2. Create an account on [Docker Hub](https://hub.docker.com/)
3. Install WSL2
4. Open a PowerShell terminal and type `wsl --install`

### Probe tracing
The tracing of fluorescent probe trajectories in 2D histological brian slices can be done with [AP_Histology](https://github.com/petersaj/AP_histology) or [Universal Probe Finder]((https://github.com/JorritMontijn/UniversalProbeFinder). Both of these are MATLAB packages so unfortunatly you will need a MATLAB license for this. As far as I know there are no python-based packages for the tracing of probe trajectories in 2D brain slices.

## First time use
After installing all the necessary components you can set up your pipeline for use.
1. Open your Anaconda or Miniforge prompt.
2. Activate your environment by typing `conda activate powerpixels` (or `mamba`)
3. Set up the configuration files with the command `powerpixels-setup`
4. Navigate to where you cloned the repository.
5. Open config/settings.json and fill in your settings, see below for detailed explanation of all settings.
6. Open config/wiring/nidq.wiring.json and fill in the synchronization channels you have in use. If you do not use a BNC breakout box you can skip this step. If you do, make sure you wire the 1Hz synchronization pulse from the PXI chassis to the digital channel 0 of the breakout box (see the paper for more details).
7. (Optional) Adjust the spike sorter parameters to your liking by editing the parameter file in the config/sorter_params folder.
8. (Optional) Adjust the parameters of Bombcell and the IBL neuron-level QC metrics in the config folder.

## Settings
<img width="490" height="419" alt="Screenshot 2025-11-13 140701" src="https://github.com/user-attachments/assets/b91e811f-19c8-4453-ab9c-85f4fcbdc762" />

- SPIKE_SORTER: which spike sorter to use, accepts all spike sorters supported by SpikeInterface.
- IDENTIFIER: this text is appended to the final data folder to destinguish multiple spike sorting runs.
- DATA_FOLDER: path to the top level folder where your data lives.
- SINGLE_SHANK: artifact removal method used for single shank probes. Options: "car_global"; global median reference, "car_local"; local median reference, "destripe"; spatial filtering.
- MULTI_SHANK: artifact removal method for probes with multiple shanks.
- LOCAL_RADIUS: *only for car_local* the radius in um around each channel to select channels to subtract from it.
- PEAK_THRESHOLD: the threshold for peak detection in the power spectrum which will be filtered out to reduce high-frequency noise.
- USE_NIDAQ: whether you use a BNC-breakout box with synchronization channels.
- USE_DOCKER: whether to run the spike sorting in a Docker container.
- COMPRESS_RAW_DATA: whether to compress the raw data.
- COMPRESSION: compression method, options: zarr or mtscomp.
- NWB_EXPORT: whether to export the spike sorting results as an NWB file.
- N_CORES: how many CPU's to use for preprocessing (-1 is all).

## Folder structure
The pipeline is in principle agnostic to how your data is organized at a high level. The session folder can have any name and can be located anywhere in your top level data directory. However, each session folder does need to abide by some formatting requirements. Inside each session folder there needs to be a raw_ephys_data folder in which should be the output folder of SpikeGLX or OpenEphys. For the pipeline to find which session folder to process you need to create a process_me.flag file and put it in the session folder.
```
├── session_folder
|   ├── raw_ephys_data
└── process_me.flag
```
To facilitate the process you can run the helper function `python scripts\prepare_sessions.py` which creates the folders and flags for you (optional).

## Data output format
The data that comes out of the Power Pixels pipeline is (1) raw spike sorter output files (e.g. Kilosort files), (2) data files in [ALF filenaming convention](https://int-brain-lab.github.io/ONE/alf_intro.html), and (3) an NWB file. A helper function is included to load in your neural data `from powerpixels import load_neural_data`.

## Usage workflow

1. Before starting a recording prepare the folder structure, either manually or by running `python scripts\prepare_sessions.py`. 
2. Perform your Neuropixel recording and make sure the output folder of SpikeGLX or the OpenEphys GUI is the "raw_ephys_data" folder.
3. Start the pipeline by running the command `python scripts\run_pipeline_spikeglx.py` or `python scripts\run_pipeline_openephys.py`, this will search your top-level data folder for any sessions that have the process_me.flag. The pipeline will take a long time to run so best to do it overnight. After the pipeline has run there will be new probe folders for each of the probes in the top-level of the session folder which contain the spike sorted data and other quality metrics.
4. After you've done your histology, launch AP_Histology or Universal Probe Finder in MATLAB and trace the fluorescent probe tracts.
5. To transform the tracks to something the alignment GUI can read run the `scripts\convert_AP_Histology_probes.m` or `scripts\convert_UniversalProbeFinder_probes.m` script in MATLAB. This will save .json files for all the probe tracks in a format the alignment GUI can read.
6. Match these tracks to the recorded probes and move the .json files to the corresponsing probe folders that were created by the pipeline. Once it's in the correct probe folder, rename the .json file to `xyz_picks.json`.
7. Launch the alignment gui by typing `python iblapps\atlaselectrophysiology\ephys_atlas_gui.py -o True`, see instructions on how to use the GUI [here](https://github.com/int-brain-lab/iblapps/wiki/2.-Usage-instructions).
8. After the alignment is done click Upload and the final channel locations and brain regions will be saved in the `channel_locations.json` file.
9. You can do manual curation of the spike sorting output by running in a python terminal:
    ```
    from powerpixels import manual_curation
    manual_curation("path\to\sorting\results")
    ```
    The SpikeInterface manual curation GUI will launch which will include the automatic classification output from Bombcell, UnitRefine and the IBL. You can use the GUI to manually annote units as good, or you can use it to see which one of the automated classification metrics you like and use one of those. They are loaded in together with the neural data so you can easily use them to filter units to use.
10. You can load in the neural data of your recording with a supplied helper function like this:
    ```
    from powerpixels import load_neural_data
    spikes, clusters, channels = load_neural_data(session_path, probe)
    ```
    For extensive documentation as to what each dataset type in `spikes`, `clusters`, and `channels` means see the documentation [here](https://docs.google.com/document/d/1OqIqqakPakHXRAwceYLwFY9gOrm8_P62XIfCTnHwstg/).
    You can filter your neurons using one of the automatic curation criteria by adding the flag `keep_units` with the options `'bombcell'`, `'unitrefine'`, `'ibl'` and `kilosort`. Alternativly you can do it yourself by using the following dataset types:

    `clusters['bombcell_label']`: 0 = noise, 1 = good single neuron, 2 = multi-unit activity, 3 = non-somatic

    `clusters['unitrefine_label']`: 0 = multi-unit activity or noise, 1 = good single neuron

    `clusters['ibl_label']`: 0 = noise, 0.33-0.66 = multi-unit activity, 1 = good single neuron
    
    `clusters['kilosort_label']`: 0 = multi-unit activity or noise, 1 = good single neuron
   
That's it, enjoy your beautiful data!

*If you like this pipeline, you can star this repository and/or give me a shoutout on Bluesky ([@guidomeijer.bsky.social](https://bsky.app/profile/guidomeijer.bsky.social)) or X ([@guido_meijer](https://x.com/guido_meijer)).*





