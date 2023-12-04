# -*- coding: utf-8 -*-
"""
Created on Mon Aug 21 11:01:06 2023

@author: Guido Meijer
"""

from os import mkdir
from os.path import join, isdir, isfile
from datetime import date
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-e", "--ephys", help="Ephys session", action='store_true')
args = parser.parse_args()

# Set path to save data
if args.ephys:
    PATH = 'D:\\NeuropixelData'  # SSD
else:
    PATH = 'K:\\Subjects'

# Get date of today
this_date = date.today().strftime('%Y%m%d')

# Get mouse name
subject_name = input('Subject name (q to quit): ')

while subject_name != 'q':
        
    # Make directories
    while not isdir(join(PATH, subject_name)):
        if not isdir(join(PATH, subject_name)):
            create_folder = input('Subject does not exist, create subject folder? (y/n) ')
            if create_folder == 'y':        
                mkdir(join(PATH, subject_name))
            else:
                subject_name = input('Subject name (q to quit): ')
            
    if not isdir(join(PATH, subject_name, this_date)):
        mkdir(join(PATH, subject_name, this_date))
        mkdir(join(PATH, subject_name, this_date, 'raw_behavior_data'))
        mkdir(join(PATH, subject_name, this_date, 'raw_video_data'))
        if args.ephys:
            mkdir(join(PATH, subject_name, this_date, 'raw_ephys_data'))    
            print(f'Created ephys session {this_date} for {subject_name}')
        else:
            print(f'Created session {this_date} for {subject_name}')
        
    # Create flags
    if not isfile(join(PATH, subject_name, this_date, 'extract_me.flag')):
        with open(join(PATH, subject_name, this_date, 'extract_me.flag'), 'w') as fp:
            pass
    if not isfile(join(PATH, subject_name, this_date, 'transfer_me.flag')):
        with open(join(PATH, subject_name, this_date, 'transfer_me.flag'), 'w') as fp:
            pass
    if args.ephys:
        if not isfile(join(PATH, subject_name, this_date, 'videotrack_me.flag')):
            with open(join(PATH, subject_name, this_date, 'videotrack_me.flag'), 'w') as fp:
                pass
    
    # Get mouse name
    subject_name = input('Subject name (q to quit): ')
            

    

