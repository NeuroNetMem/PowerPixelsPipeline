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




