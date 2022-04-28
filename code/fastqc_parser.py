#!/usr/local/bin/python

####    Extract some info from fastqc .zip output    ####

## Copyright Marta Binaghi marta.binaghi <at> ips.unibe.ch
## 20th November 2018
## Last modified on 20th November 2018


import argparse
from os import listdir
from os.path import isfile, join
import re
import zipfile
import os
import shutil

def initialise_variables():
    global modules
    global header
    # store names of modules and their abbreviations
    modules = {'Basic Statistics': 'bs',
               'Per base sequence quality': 'pbsq',
               'Per tile sequence quality': 'ptsq',
               'Per sequence quality scores': 'psqs',
               'Per base sequence content': 'pbsc',
               'Per sequence GC content': 'psgc',
               'Per base N content': 'pbnc',
               'Sequence Length Distribution': 'slen',
               'Sequence Duplication Levels': 'sdup',
               'Overrepresented sequences': 'ors',
               'Adapter Content': 'adc'}
    header = ['filename',
              'totseq',
              'poor',
              'len',
              'bs',
              'pbsq',
              'ptsq',
              'psqs',
              'pbsc',
              'psgc',
              'pbnc',
              'slen',
              'sdup',
              'ors',
              'adc']

def extract_fastqc_data( path ):
    # list of extracted files
    data_files = []
    onlyfiles = [f for f in listdir(path) if isfile(join(path, f))]
    # look for files matching _fastqc.zip
    pattern = '\w+_fastqc.zip'
    regexp = re.compile(pattern)
    zipfiles = [f for f in onlyfiles if regexp.search(f) is not None]
    del(pattern, regexp)
    # make a temporary folder to store the data
    fqdir = path + r'fq_data_tmp/'
    if not os.path.exists(fqdir):
        os.makedirs(fqdir)
    # extract fastqc_data.txt
    for zf in zipfiles:
        filename = zf.strip('.zip')
        zip_data = zipfile.ZipFile(path + zf, 'r')
        zip_data.extract(filename + '/fastqc_data.txt', path = fqdir)
        zip_data.close()
        oldpath = fqdir + filename + r'/fastqc_data.txt'
        shutil.move(oldpath, fqdir + '/' + filename + r'_data.txt')
        os.rmdir(fqdir + filename)
        data_files = data_files + [fqdir + filename + r'_data.txt']
    return data_files

def initialise_summary_file( path ):
    # make a new empty file to store the csv summary
    global summaryfile
    summaryfile = path + 'fastqc_summary.csv'
    open(summaryfile, 'w').close()
    # add header
    with open(summaryfile, 'a') as file:
        file.write(','.join(header) + '\n')

def parse_fastqc_data( input_file_path ):
    thisfile = open(input_file_path,'r')
    new_line_dict = dict.fromkeys(header, '')
    current_module = None
    for line in thisfile:
        if line[0:2] == '>>':
            if line[0:12] == '>>END_MODULE':
                current_module = None
                continue
            for key in modules.keys():
                pattern = '>>(?P<section>' + key + ')\s(?P<info>\w{4})'
                regexp = re.compile(pattern)
                result = regexp.search(line)
                if result is not None:
                    current_module = result.group('section')
                    info = result.group('info')
                    new_line_dict[modules[current_module]] = info
        elif current_module == 'Basic Statistics':
            if line[0:9] == 'Filename\t':
                new_line_dict['filename'] = line.strip('Filename\t').strip('\n')
            elif line[0:5] == 'Total':
                new_line_dict['totseq'] = line.strip('Total Sequences\t').strip('\n')
            elif line[0:13] == 'Sequences fla':
                new_line_dict['poor'] = line.strip('Sequences flagged as poor quality\t').strip('\n')
            elif line[0:12] == 'Sequence len':
                new_line_dict['len'] = line.strip('Sequence length\t').strip('\n')
        else:
            continue
    thisfile.close()
    return new_line_dict

def append_line_to_csv( new_line_dict ):
    newline = []
    for column in header:
        newline = newline + [new_line_dict[column]]
    # append new line to file
    with open(summaryfile, 'a') as file:
        file.write(','.join(newline) + '\n')

####
parser = argparse.ArgumentParser(description='Parse a set of fastqc output files.')
parser.add_argument('-p', '--path', dest='path', type=str, nargs=1,
                    action='store',
                    help='path to directory where .zip input is')

args = parser.parse_args()
path = args.path[0]

if path[len(path)-1] != '/':
    path = path + '/'

data_files = extract_fastqc_data( path )
initialise_variables()
initialise_summary_file( path )
for data_file in data_files:
    append_line_to_csv(parse_fastqc_data(data_file))
shutil.rmtree(path + 'fq_data_tmp')
