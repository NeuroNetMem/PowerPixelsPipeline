# -*- coding: utf-8 -*-
"""
Created on Mon Aug 21 11:01:06 2023

@author: Guido Meijer
"""

import json
from os import mkdir
from os.path import join, dirname, realpath, isdir, isfile
from datetime import date

# Get data folder path
with open(join(dirname(realpath(__file__)), 'settings.json'), 'r') as openfile:
    settings_dict = json.load(openfile)

# Get date of today
this_date = date.today().strftime('%Y%m%d')

# Get mouse name
subject_name = input('Subject name (q to quit): ')

while subject_name != 'q':
        
    # Make directories
    while not isdir(join(settings_dict['DATA_FOLDER'], subject_name)):
        if not isdir(join(settings_dict['DATA_FOLDER'], subject_name)):
            create_folder = input('Subject does not exist, create subject folder? (y/n) ')
            if create_folder == 'y':        
                mkdir(join(settings_dict['DATA_FOLDER'], subject_name))
            else:
                subject_name = input('Subject name (q to quit): ')
            
    if not isdir(join(settings_dict['DATA_FOLDER'], subject_name, this_date)):
        mkdir(join(settings_dict['DATA_FOLDER'], subject_name, this_date))
        mkdir(join(settings_dict['DATA_FOLDER'], subject_name, this_date, 'raw_behavior_data'))
        mkdir(join(settings_dict['DATA_FOLDER'], subject_name, this_date, 'raw_video_data'))
        mkdir(join(settings_dict['DATA_FOLDER'], subject_name, this_date, 'raw_ephys_data'))    
        print(f'Created session {this_date} for {subject_name}')
        
    # Create flags
    if not isfile(join(settings_dict['DATA_FOLDER'], subject_name, this_date, 'process_me.flag')):
        with open(join(settings_dict['DATA_FOLDER'], subject_name, this_date, 'process_me.flag'), 'w') as fp:
            pass
   
    # Get mouse name
    subject_name = input('Subject name (q to quit): ')
            

    

