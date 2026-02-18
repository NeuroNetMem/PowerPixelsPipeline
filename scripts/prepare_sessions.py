# -*- coding: utf-8 -*-
"""
Created on Mon Aug 21 11:01:06 2023

@author: Guido Meijer
"""

import json
from datetime import date
from pathlib import Path

# Get data folder path
project_root = Path(__file__).parent.parent
settings_file = project_root / 'config' / 'settings.json'
with open(settings_file, 'r') as openfile:
    settings_dict = json.load(openfile)
data_folder = Path(settings_dict['DATA_FOLDER'])

# Get date of today
this_date = date.today().strftime('%Y%m%d')

# Get mouse name
subject_name = input('Subject name (q to quit): ')

while subject_name != 'q':
        
    # Make directories
    while not (data_folder / subject_name).is_dir():
        if not (data_folder / subject_name).is_dir():
            create_folder = input('Subject does not exist, create subject folder? (y/n) ')
            if create_folder == 'y':        
                (data_folder / subject_name).mkdir()
            else:
                subject_name = input('Subject name (q to quit): ')
            
    if not (data_folder / subject_name/ this_date).is_dir():
        (data_folder / subject_name / this_date).mkdir()
        (data_folder / subject_name / this_date / 'raw_behavior_data').mkdir()
        (data_folder / subject_name / this_date / 'raw_video_data').mkdir()
        (data_folder / subject_name / this_date / 'raw_ephys_data').mkdir()
        print(f'Created session {this_date} for {subject_name}')
        
    # Create flags
    if not (data_folder / subject_name / this_date / 'process_me.flag').is_file():
        with open(data_folder / subject_name / this_date / 'process_me.flag', 'w') as fp:
            pass
   
    # Get mouse name
    subject_name = input('Subject name (q to quit): ')
            

    

