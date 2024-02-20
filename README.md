# Power Pixels pipeline for processing of Neuropixel recordings

The Power Pixels pipeline combines the best parts of several pipelines (and even MATLAB) into one Python-based pipeline. It supports Neuropixel 1.0 probes (2.0 support coming soon) recorded on a National Instruments system (tested on NI PIXe-1071 with a BNC-2110 breakout board for synchronization channels) using SpikeGLX. OpenEphys GUI support can be implemented if there is interest. 

This pipeline is nothing new! It's all about combining existing modules and pipelines into one, which is especially useful for people who are just starting out doing Neuropixel recordings and maybe have heard of all these tools but need some help getting them all integrated. Power Pixels relies on these amazing open-source projects:
- [SpikeInterface](https://spikeinterface.readthedocs.io)
- [ibllib](https://github.com/int-brain-lab/ibllib)
- [Kilosort](https://github.com/MouseLand/Kilosort)
- [Bombcell](https://github.com/Julie-Fabre/bombcell)
- [Universal Probe Finder](https://github.com/JorritMontijn/UniversalProbeFinder)

## Description of the pipeline elements

![image](https://github.com/NeuroNetMem/PowerPixelsPipeline/assets/19360723/70ebe00d-ef57-4550-99d9-5d00a3d818f0)

The pipeline goes through the following steps:
- **Phase shift correction**: channels on a Neuropixel probe are not recorded simultaneously, there is a small (30 microsecond) delay between the acquisition of a block of channels. Correcting for this small delay greatly improves artifact removal at the "Destriping" step.
- **Remove bad channels**: bad channels are detected by looking at both coherence with other channels and PSD power in the high-frequency range, then they are interpolated using neighboring channels.
- **Destriping**: removes electrical artifacts by applying a high pass spatial filter over the depth of the probe.
- **Spike sorting**: a spike sorting algorithm is used to detect spikes and sort them into units. SpikeInterface supports many [spike sorters](https://spikeinterface.readthedocs.io/en/latest/modules/sorters.html#supported-spike-sorters) out of the box 
- **Bombcell**: a MATLAB package that calculates neuron-level QC metrics (optional).
- **Synchronization**: each Neuropixel probe and the BNC breakout box has their own clock. This means one has to synchronize the spike times between the probes (if you use more than one) and the synchronization channels which carry timestamps of events (for example: behavioral events or pulses from a camera).
- **Compression**: the raw binary file is compressed using *mtscomp* which results in a 2-3x reduction in file size.
- **Histological tracing**: the fluorescent tracks of the probes are traced using the Universal Probe Finder package.
- **Ephys-histology alignment**: the brain regions along the probe, inferred from the tracing, are aligned to electrophysiological features.

## Installation

It is recommended to install Power Pixels in an Anaconda or Miniforge environment.
1. Install [Anaconda](https://www.anaconda.com/) or [Miniforge](https://github.com/conda-forge/miniforge) - Miniforge is the recommended option
2. Open the Anaconda or Miniforge prompt
3. Install git by typing ```conda install git``` or ```mamba install git``` depending on whether you use Anaconda or Miniforge, respectively
4. Navigate to the location on your computer you want the repository to be and clone the repository by typing ```git clone https://github.com/NeuroNetMem/PowerPixelsPipeline```
5. Additonally clone the iblapps repository ```git clone https://github.com/int-brain-lab/iblapps/```
6. Create the environment from the provided YML file by typing ```conda env create -f environment.yml``` or ```mamba env create -f environment.yml```
7. You can now activate the environment by typing ```conda activate powerpixels```
8. Install the iblapps repository by typing `pip install -e iblapps`

### Docker
SpikeInterface uses Docker to launch spike sorters in a docker container, this is great because it means that you don't need to tinker with grapic card drivers or have MATLAB installed. Instructions to set up Docker on Windows:
1. Install [Docker Desktop](https://www.docker.com/products/docker-desktop/)
2. Create an account on [Docker Hub](https://hub.docker.com/)
3. Install WSL2
4. Open a PowerShell terminal and type ```wsl --install```

### MATLAB 
At the moment a MATLAB license is necessary to convert the output of Universal Probe Finder to the input of the alignment GUI. The spike sorting is run in a Docker by SpikeInterface, which means that you can even run Kilosort without having MATLAB installed. The Bombcell package is optional and is turned off by default. To run Universal Probe Finder without MATLAB, it has the option to download a compiled version of the package (see [here](https://github.com/JorritMontijn/UniversalProbeFinder?tab=readme-ov-file#using-the-universal-probe-finder-without-a-matlab-license)).
#### Bombcell
If you want to use Bombcell as part of the pipeline you need to set up the MATLAB python engine so that python can run the toolbox. Follow [these](https://nl.mathworks.com/help/matlab/matlab_external/install-the-matlab-engine-for-python.html) instructions to set up the engine (tested with MATLAB 20223b).
#### Universal Probe Finder
Install the Universal Probe Finder by following the instructions [here](https://github.com/JorritMontijn/UniversalProbeFinder). In short: first do ```git clone --recursive https://github.com/JorritMontijn/UniversalProbeFinder``` and add with subfolders to the MATLAB path.

## First time use

After installing all the necessary components you can set up your pipeline for use.
1. Open your Anaconda or Miniforge prompt.
2. Activate your environment by typing ```conda activate powerpixels```
3. Navigate to where you cloned the repository.
4. Generate setting JSON files ```python PowerPixelsPipeline\generate_json_files.py```.
5. Open settings.json and fill in your settings (explanations of each item can be found in generate_setting_files.py).
6. Open nidq.wiring.json and fill in the synchronization channels you have in use.
7. (Optional) If you are planning on using a spike sorter other than Kilosort 2.5 or 3, you can generate the parameter file for this sorter by typing ```python PowerPixelsPipeline\get_default_sorter_params.py -s SPIKE_SORTER``` in which SPIKE_SORTER should be the sorter you want to use (see [here](https://spikeinterface.readthedocs.io/en/latest/modules/sorters.html#supported-spike-sorters) for all options). The default parameter file for your sorter will appear in the spikesorter_param_file folder and you can adjust any parameters you wish in there.

*Recommended parameters for Kilosort 2.5 and 3 are provided in the spikesorter_param_files folder, you can change these if you want but bear in mind that, because they come with the repository, they will be overwritten when you pull any new changes from the repo.*

## Folder structure
The pipeline is in principle agnostic to how your data is organized at a high level. The session folder can have any name and can be located anywhere in your top level data directory. However, each session folder does need to abide by some formatting requirements. Inside each session folder there needs to be a raw_ephys_data folder in which there are several folders for the different probes (named probe00, probe01, etc), these folders should be the output folders of SpikeGLX. For the pipeline to find which session folder to process you need to create a process_me.flag file and put it in the session folder.
```
├── session_folder
|   ├── raw_ephys_data
|   |   ├── probe00
|   |   ├── probe01
└── process_me.flag
```
To facilitate the process you can run the helper function `python PowerPixelsPipeline\prepare_sessions.py` which creates the folders and flags for you (optional).

## Usage workflow

1. Before starting a recording prepare the folder structure, either manually or by running `python PowerPixelsPipeline\prepare_sessions.py`. 
2. Perform your Neuropixel recording and make sure the output folder of SpikeGLX is one of the probe folders in raw_ephys_data.
3. Start the pipeline by running the command `python PowerPixelsPipeline\run_pipeline.py`, this will search your top-level data folder for any sessions that have the process_me.flag. The pipeline will take a long time to run so best to do it overnight. After the pipeline has run there will be new probe folders for each of the probes in the top-level of the session folder which contain the spike sorted data and other quality metrics.
4. After you've done your histology, launch Universal Probe Finder in MATLAB and do the Slice Prepper and Slice Finder steps to trace your probe insertions (you can skip Probe Finder).
5. To transform the tracks to something the alignment GUI can read run the `convert_histology_to_alignment_GUI.m` script in MATLAB. This will save .json files for all the tracks you traced in Universal Probe Finder.
6. Match these tracks to the recorded probes and move the .json files to the corresponsing probe folders that were created by the pipeline. Once it's in the correct probe folder, rename the .json file to `xyz_picks.json`.
7. Launch the alignment gui by typing `python iblapps\atlaselectrophysiology\ephys_atlas_gui.py -o True`, see instructions on how to use the GUI [here](https://github.com/int-brain-lab/iblapps/wiki/2.-Usage-instructions).
8. After the alignment is done click Upload and the final channel locations and brain regions will be saved in the `channel_locations.json` file.
9. You can do manual curation of the spike sorting output using [Phy](https://github.com/cortex-lab/phy) and launching it in the probe folder. If you've run BombCell the quality metrics should appear in Phy automatically.
10. That's it, enjoy your beautiful data!

*If you like this pipeline, you can star this repository and/or give me a shoutout on Twitter ([@guido_meijer](https://twitter.com/guido_meijer)).*



