#!/usr/local/bin/python

####    Extract numbers from ngsAdmix log files    ####

## Copyright Marta Binaghi marta.binaghi at ips.unibe.ch
## 17th May 2019
## Last modified: 17th May 2019


import argparse
from os import listdir
from os.path import isfile, join
import re
import zipfile
import os
import shutil

def initialise_variables():
    global header
    header = ['K',
              'rep',
              'seed',
              'likelihood',
              'iterations']

def get_logfile_list( path ):     # list of log files in folder
    onlyfiles = [f for f in listdir(path) if isfile(join(path, f))]
    # look for files matching .log
    pattern = '\w+\.log'
    regexp = re.compile(pattern)
    log_files = [f for f in onlyfiles if regexp.search(f) is not None]
    del(pattern, regexp)
    return log_files

def initialise_csv_file( path ):
    # make a new empty file to store the csv
    global csvfile
    csvfile = path + 'log_stats.csv'
    open(csvfile, 'w').close()
    # add header
    with open(csvfile, 'a') as file:
        file.write(','.join(header) + '\n')

def parse_logfile( input_file_path ):
    thisfile = open(input_file_path,'r')
    new_line_dict = dict.fromkeys(header, '')
    for line in thisfile:
        if line[0:5] == 'Input':
            pattern = 'nPop=(?P<info>\d{1,2}),.+outfiles=.*run(?P<myrep>\d{1,4})$'
            regexp = re.compile(pattern)
            result = regexp.search(line)
            if result is not None:
                info = result.group('info')
                new_line_dict['K'] = info
                myrep = result.group('myrep')
                new_line_dict['rep'] = myrep
        elif line[0:5] == 'Setup':
            pattern = 'seed=(?P<info>-?\w+)\snThr'
            regexp = re.compile(pattern)
            result = regexp.search(line)
            if result is not None:
                info = result.group('info')
                new_line_dict['seed'] = info
        elif line[0:4] == 'best':
            pattern = 'like=(?P<info>-?\d+\.\d+) after (?P<myiter>\d+) iterations'
            regexp = re.compile(pattern)
            result = regexp.search(line)
            if result is not None:
                info = result.group('info')
                new_line_dict['likelihood'] = info
                myiter = result.group('myiter')
                new_line_dict['iterations'] = myiter
        else:
            continue
    thisfile.close()
    return new_line_dict

def append_line_to_csv( new_line_dict ):
    newline = []
    for column in header:
        newline = newline + [new_line_dict[column]]
    # append new line to file
    with open(csvfile, 'a') as file:
        file.write(','.join(newline) + '\n')

####
parser = argparse.ArgumentParser(description='Parse a set of ngsAdmix log files.')
parser.add_argument('-p', '--path', dest='path', type=str, nargs=1,
                    action='store',
                    help='path to directory where .log files are')

args = parser.parse_args()
path = args.path[0]

if path[len(path)-1] != '/':
    path = path + '/'

log_files = get_logfile_list( path )
initialise_variables()
initialise_csv_file( path )
for log_file in log_files:
    append_line_to_csv(parse_logfile(log_file))
