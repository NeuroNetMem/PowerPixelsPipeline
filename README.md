# Pixelzord Pipeline for pre/post-processing and spike sorting of Neuropixel recordings

The Pixelzord pipeline combines the best parts of several pipelines (and even MATLAB) into one Python-based pipeline. It supports Neuropixel 1.0 probes (2.0 support coming soon) recorded on a National Instruments system (tested on NI PIXe-1071 with a BNC-2110 breakout board for synchronization channels) using SpikeGLX. OpenEphys GUI support can be implemented if there is interest. 

This pipeline is nothing new! It's all about combining existing modules and pipelines into one, which is especially useful for people who are just starting out doing Neuropixel recordings and maybe have heard of all these tools but need some help getting them all integrated. Pixelzord relies on these amazing open-source projects:
- SpikeInterface (https://spikeinterface.readthedocs.io)
- ibllib (https://github.com/int-brain-lab/ibllib)
- Kilosort (https://github.com/MouseLand/Kilosort)
- Bombcell (https://github.com/Julie-Fabre/bombcell)
- Phy (https://github.com/cortex-lab/phy)

## Description of the pipeline elements

![image](https://github.com/NeuroNetMem/PixelzordPipeline/assets/19360723/7d043851-de5e-4605-a5fa-48c036f68988)

The pipeline goes through the following steps:
- **Phase shift correction**: channels on a Neuropixel probe are not recorded simultaneously, there is a small (30 microsecond) delay between the acquisition of a block of channels. Correcting for this small delay greatly improves artifact removal at the "Destriping" step.
- **Remove bad channels**: bad channels are detected by looking at both coherence with other channels and PSD power in the high-frequency range, then they are interpolated using neighboring channels.
- **Destriping**: removes electrical artifacts by applying a high pass spatial filter over the depth of the probe.
- **Spike sorting**: a spike sorting algorithm is used to detect spikes and sort them into units. SpikeInterface supports many spike sorters out of the box (https://spikeinterface.readthedocs.io/en/latest/modules/sorters.html#supported-spike-sorters)
- **Bombcell**: a MATLAB package that calculates neuron-level QC metrics.
- **Synchronization**: each Neuropixel probe and the BNC breakout box has their own clock. This means one has to synchronize the spike times between the probes (if you use more than one) and the synchronization channels which carry timestamps of events (for example: behavioral events or pulses from a camera).
- **Compress raw data**: the raw binary file is compressed using *mtscomp* which results in a 2-3x reduction in file size.

## Installation

The best is to make a conda or mamba environment for the Pixelzord pipeline. Clone this repository and navigate to it and do ```conda env create -f environment.yml``` or ```mamba env create -f environment.yml``` depending on whether you use anaconda or mamba forge.

### Docker
SpikeInterface uses Docker to launch spike sorters in a docker container, this is great because it means that you don't need to tinker with grapic card drivers or have MATLAB installed. Instructions to set up Docker on Windows
1. Install Docker Desktop (https://www.docker.com/products/docker-desktop/)
2. Create an account on Docker Hub (https://hub.docker.com/)
3. Install WSL2
4. Open a PowerShell terminal and type ```wsl --install```

### MATLAB 
If you want to use the Bombcell (can be turned off in settings) you need to have MATLAB installed and you need to set up the MATLAB python engine so that python can run the toolbox. Follow these instructions to set up the engine (tested with MATLAB 20223b): https://nl.mathworks.com/help/matlab/matlab_external/install-the-matlab-engine-for-python.html

### Phy
For manual curation of the spike sorting output you need to install Phy, follow the instructions here (recommended to install in it's own environment): https://github.com/cortex-lab/phy

## First time use

After installing all the necessary components you can set up your pipeline for use.
1. Activate your environment ```conda activate pixelzord```
2. Navigate to the cloned repository
3. Generate setting JSON files ```python generate_setting_files.py```
4. Open settings.json and fill in your settings (explanations of each item can be found in generate_setting_files.py)
5. Open nidq.wiring.json and fill in the synchronization channels you have in use

## Usage workflow

1. Before starting a recording prepare the folder structure by doing ```python prepare_sessions.py```, this will create a folder for your animal, the recording day, and several folders for raw data. It also creates a spikesort_me.flag that will be used by the pipeline to find recordings that need to be spike sorted.
2. Perform your Neuropixel recording and make sure the output folder of SpikeGLX is the created raw_ephys_data folder.
3. Start the pipeline by doing ```python run_pipeline.py```, this will search your data folder for sessions that have the spikesort_me.flag and process them. The pipeline will take a long time to run so best to do it overnight. 



