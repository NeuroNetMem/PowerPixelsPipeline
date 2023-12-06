# -*- coding: utf-8 -*-
"""
Created on Tue Dec  5 12:15:39 2023

By Guido Meijer
"""

from spikeinterface.sorters import get_default_sorter_params
import argparse
import json
from os.path import join, dirname, realpath

parser = argparse.ArgumentParser()
parser.add_argument("-s", "--sorter", help="Spike sorter (for example 'kilosort2_5)'")
args = parser.parse_args()

# Get default sorter params
sorter_params = get_default_sorter_params(args.sorter)

# Save to disk
with open(join(dirname(realpath(__file__)),
               'spikesorter_param_files', f'{args.sorter}_params.json'), 'w') as outfile:
    outfile.write(json.dumps(sorter_params, indent=4))
