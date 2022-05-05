#!/usr/bin/python3

#-------------------------------------------------------------------

# GATK PL to Beagle GL format
#
# Take a vcf file from GATK or similar, extract PL values and
#   transform them into genotype likelihoods. Like:
#   marker,allele1,allele2,Ind0,Ind0,Ind0,Ind1,Ind1,Ind1
#
#
# Marta Binaghi <marta.binaghi <at> ips.unibe.ch>
#   Kuhlemeier group, University of Bern
#   Natural axillaris X exserta hybrids project
# Created: August 12th, 2019
# Last modified: August 12th, 2019

#-------------------------------------------------------------------

import os.path
import csv
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument("-i", "--input", type=str,
                    help="Input file (vcf)")
parser.add_argument("-o", "--output", type=str,
                    help="Output file name (Beagle format)")

args = parser.parse_args()

if args.input == None:
    print("\nAn input file must be supplied via the -i argument.\n")
    parser.print_help()
    parser.exit()
elif not os.path.isfile(args.input):
    print("\nThe input file specified is not valid. Is the path correct?\n")
    parser.print_help()
    parser.exit()

if args.output == None:
    args.output = args.input.strip(".vcf") + "_beagle.gz"

print("The input vcf file is:")
print(args.input)
print("The output will be saved in:")
print(args.output)

def pltogl(pls) :
    unscaledGl = list()
    scaledGl = list()
    for p in pls:
        unscaledGl.append( 10**(-(int(p)/10)) )
    for g in unscaledGl:
        scaledGl.append( g / sum(unscaledGl) )
    return(scaledGl)



f = open(args.input)
for line in f:
    # test for header
    if line[0:6] == "#CHROM":
        header = line.split()
        samples = header[9:len(header)]
        sampleRef = list()
        for smplidx in range(len(samples)) :
            sampleRef.append([smplidx,
                              samples[smplidx]
                              ])
        print("\n" + str(len(samples)) + " samples were detected in the vcf. Their new index is shown as row number:")
        for i in sampleRef:
            print(*i)
        #print()
        # initialise dictionary and a list that will be the header of the output csv file
        df = {
            "marker" : [None],
            "allele1" : [None],
            "allele2" : [None]
        }
        myheaderKeys = ["marker",
                        "allele1",
                        "allele2"]
        for smpl in samples:
            smplidx = "Ind" + str(samples.index(smpl))
            df[smplidx + "_00"] = [None]
            df[smplidx + "_01"] = [None]
            df[smplidx + "_11"] = [None]
            myheaderKeys.extend([smplidx + "_00",
                                 smplidx + "_01",
                                 smplidx + "_11"])
        continue
    elif line[0] == "#":
        continue
    # read line
    pos = line.split()
    df["marker"].append(pos[0] + "_" + pos[1] )
    df["allele1"].append(pos[3])
    df["allele2"].append(pos[4])
    # loop through samples in the line
    for smpl in samples:
        #print(samples.index(smpl))
        smplidx = "Ind" + str(samples.index(smpl))
        # test if sample is genotyped
        smplGT = pos[samples.index(smpl) + 9].split(sep=":")[pos[8].split(sep=":").index("GT")]
        if smplGT == "./." :
            df[smplidx + "_00"].append(0.333333)
            df[smplidx + "_01"].append(0.333333)
            df[smplidx + "_11"].append(0.333333)
        else :
            # get PL from format field into a three-items list
            smplPL = pos[samples.index(smpl) + 9].split(sep=":")[pos[8].split(sep=":").index("PL")]
            smplPL = smplPL.split(sep=",")
            smplPL = [ int(i) for i in smplPL ]
            #print(smplPL)
            # calculate GL
            GLs = pltogl(smplPL)
            #append GLs to posDic under keys smplidx+"00", "01", "11"
            df[smplidx + "_00"].append(GLs[0])
            df[smplidx + "_01"].append(GLs[1])
            df[smplidx + "_11"].append(GLs[2])
    #print(df)
f.close()

# remove first entry of each dictionary list because it was the None used to initialise
for key in df.keys():
    del df[key][0]

realHeader = myheaderKeys[0:3]
for ind in samples:
    realHeader.extend(["Ind" + str(samples.index(ind))]*3)
#print(realHeader)
with open(args.output, "w") as out :
    w = csv.writer(out, )
    w.writerow(realHeader)
    w.writerows(zip(*[df[key] for key in myheaderKeys]))

