#!/bin/bash

# Renaming the raw read files for
# natural ax ex hybrids project.
# Marta Binaghi
# Created April 4th, 2019
# Last modified March 22nd, 2022

# !! the script is tailored on the csv where the
# IDs are recorded. If you re-use the script
# you have to define the column ID and adapt it
# to your needs.

id_file=$1

file_folder=$2

cd $file_folder

for file in *.fastq.gz;
do
  sraId="$(cut -d'.' -f 1 <<< $file)"
#  echo $sraId
  line=$(grep -P "${sraId}" ${id_file})
#  echo $line
  newId="$(cut -f 6 -d ',' <<< $line)"  
#  echo $newId
  newLane="$(cut -f 12 -d ',' <<< $line)"
  newName="$(echo ${newId}_L${newLane}_${file})"
#  echo ${newName}
  echo "file ${file} moved to ${newName}" 
#  echo ""
  mv ${file} ${newName}  
done


