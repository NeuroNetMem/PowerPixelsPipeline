# Pixelzord Pipeline for pre/post-processing and spike sorting of Neuropixel recordings

Pixelzord pipeline combines the best parts of several pipelines (and even MATLAB) into one Python-based pipeline. It supports Neuropixel 1.0 probes (2.0 support coming soon) recorded on a National Instruments system (tested on NI PIXe-1071 with a BNC-2110 breakout board for synchronization channels) using SpikeGLX. OpenEphys GUI support can be implemented if there is interest. 

This pipeline is nothing new! It's all about combining existing modules and pipelines into one, which is especially usefull for people who are just starting out doing Neuropixel recordings and maybe have heard of all these tools but need some help getting them all integrated. Megazord relies on these amazing open-source projects:
- SpikeInterface (https://spikeinterface.readthedocs.io)
- ibllib (https://github.com/int-brain-lab/ibllib)
- Bombcell (https://github.com/Julie-Fabre/bombcell)
- Phy (https://github.com/cortex-lab/phy) 



