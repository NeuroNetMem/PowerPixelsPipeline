# -*- coding: utf-8 -*-
"""
Created on Wed Feb 25 10:40:54 2026

By Guido Meijer
"""

from pathlib import Path
import matplotlib.pyplot as plt
import spikeinterface.full as si
from spikeinterface.sortingcomponents.motion import estimate_motion, interpolate_motion


# Settings
PROBE_PATH = Path(r'path\to\raw_ephys_data\probe00')
N_CPUS = -1
si.set_global_job_kwargs(n_jobs=N_CPUS, progress_bar=True)

# Load in raw broadband data
rec_raw = si.read_spikeglx(PROBE_PATH, stream_id=si.get_neo_streams('spikeglx', PROBE_PATH)[0][0])

# %% Do CAR and save to disk

# Filter out LFP band
rec_lfp = si.bandpass_filter(rec_raw, freq_min=1, freq_max=400)

# Correct for inter-sample shift
rec_shifted = si.phase_shift(rec_lfp)    

# Do common average reference
rec_car = si.common_reference(rec_shifted)

# Downsample to 2500 Hz
rec_down = si.resample(rec_car, 2500)

# Save to disk
rec_down.save(folder=PROBE_PATH / 'lfp_car', format='binary', chunk_duration='1s',
              dtype='int16', n_jobs=N_CPUS)

# %% Do motion correction

# Filter out LFP band
rec_lfp = si.bandpass_filter(rec_raw, freq_min=2, freq_max=100)

# Correct for inter-sample shift
rec_shifted = si.phase_shift(rec_lfp)    

# Do common average reference
rec_car = si.common_reference(rec_shifted)

# Downsample to 250 Hz
rec_way_down = si.resample(rec_car, 250)

# Do DREDge motion detection
print('Do motion correction')
lfprec_1d = si.average_across_direction(rec_way_down, direction='y')
motion = estimate_motion(lfprec_1d, method='dredge_lfp', rigid=True, progress_bar=True)
rec_down_float = si.astype(rec_down, dtype='float32')  # convert to float for interpolation
rec_motion_corr = interpolate_motion(recording=rec_down_float, motion=motion, border_mode='force_zeros')
rec_motion_corr = si.astype(rec_motion_corr, dtype='int16')

# Save to disk
rec_motion_corr.save(folder=PROBE_PATH / 'lfp_motion_corr', format='binary', chunk_duration='1s',
                     dtype='int16', n_jobs=N_CPUS)

# %% Plot
# Plot motion
fig, ax = plt.subplots()
si.plot_motion(motion, mode='line', ax=ax)
ax.set_xlim(0, 2)

# Plot traces
channel_ids = rec_raw.get_channel_ids()
f, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 7))
w = si.plot_traces(recording=rec_down, backend="matplotlib", mode='line',
                   channel_ids=channel_ids[:20], figtitle='Before motion correction',
                   time_range=[10, 11], show_channel_ids=False, ax=ax1)
w = si.plot_traces(recording=rec_motion_corr, backend="matplotlib", mode='line',
                   channel_ids=channel_ids[:20], figtitle='After motion correction',
                   time_range=[10, 11], show_channel_ids=False, ax=ax2)
