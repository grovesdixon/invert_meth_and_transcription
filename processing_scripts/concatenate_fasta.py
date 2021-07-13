#!/usr/bin/env python
##concatenate_fasta.py
##written 2/23/15 by Groves Dixon
ProgramName = 'concatenate_fasta.py'
LastUpdated = '2/23/15'
By = 'Groves Dixon'
VersionNumber = '1.0'
print "\nRunning Program {}...".format(ProgramName)

VersionString = '{} version {} Last Updated {} by {}'.format(ProgramName, VersionNumber, LastUpdated, By)

Description = '''
Description:
This script concatenates all sequences in a fasta file into one single long sequence.

'''

AdditionalProgramInfo = '''
Additional Program Information:
The purpose of this is to run an entire set of coding sequences as one long sequence in dnaSP.
If the argument -sub is supplied then only a subset of sequences given as a table will be concatenated.
Otherwise the entire set of sequences are concatenated.
'''

##Import Modules 

import time
import argparse
from sys import argv
from sys import exit
import numpy as np
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio import SeqIO
Start_time = time.time() ##keeps track of how long the script takes to run


##Set Up Argument Parsing
parser = argparse.ArgumentParser(description=Description, epilog=AdditionalProgramInfo) ##create argument parser that will automatically return help texts from global variables above
parser.add_argument('-i', required = False, dest = 'input', help = 'The the input file')
parser.add_argument('-output', '-o', required = True, dest = 'out', help = 'The desired name for the output file')
parser.add_argument('-sub', required = False, default = "NONE", dest = 'sub', help = 'A subset of sequence names to pull and concatenate')
parser.add_argument('-checkCoding', required = False, default = "NONE", dest = 'checkCoding', help = 'Optional argument to check that the sequence lengths are divisible by three. Change this to "yes" if you are concatenating coding sequences for codon analyses')
args = parser.parse_args()

#Assign Arguments
InfileName = args.input
OutfileName = args.out
SubsetFile = args.sub
CheckCoding = args.checkCoding
if CheckCoding != "NONE":
    print "\nChecking that lengths of all sequences to be concatenated are divisible by three."
nucleotideList = ['A', 'T', 'C', 'G', 'a', 't', 'c', 'g', 'N']

def read_sub(SubsetFile):
    """Reads in the subset file"""
    subList = []
    with open(SubsetFile, 'r') as infile:
        for line in infile:
            subList.append(line.strip("\n"))
    print "\nLooking for {} sequences to concatenate".format(len(subList))
    return subList

def read_file(InfileName, SubsetFile):
    '''Function to read in a file as a list of lists
    '''
    #iterate through the seqs
    numSeqs = 0
    if SubsetFile != "NONE":
        subList = read_sub(SubsetFile)
    with open(OutfileName, 'w') as out:
        out.write('>{}\n'.format(OutfileName))
        fasSeqs = SeqIO.parse(open(InfileName), 'fasta')
        if SubsetFile == "NONE":
            print "\nNo subset file detected. Concatenating the entire fasta"
            for seq in fasSeqs:
                numSeqs += 1
                newSeq = str(seq.seq)
                #replace any noncononical nucleotides with N
                editSeq = ''
                for i in newSeq:
                   if i in nucleotideList:
                       editSeq = editSeq + i
                   else:
                       editSeq = editSeq + 'N'
                out.write(editSeq)
        else:
            for seq in fasSeqs:
                if seq.id in subList:
                    numSeqs += 1
                    newSeq = str(seq.seq)
                    #double check that all the sequences are indeed divisible by three
                    if CheckCoding != "NONE":
                        x = len(newSeq)
                        if x % 3 != 0:
                            print 'Problem! The sequence is not divisible by three!'
                            print newSeq
                            exit()  
                    #replace any noncononical nucleotides with N
                    editSeq = ''
                    for i in newSeq:
                        if i in nucleotideList:
                            editSeq = editSeq + i
                        else:
                            editSeq = editSeq + 'N'
                    out.write(editSeq)
            if numSeqs == len(subList):
                print "\nNice, found all {} of the sequences from the subset file".format(numSeqs)
            else:
                print "\nFound {} of the {} sequences in the subset file".format(numSeqs, len(subList))
 
                



read_file(InfileName, SubsetFile)
if CheckCoding != "NONE":
    print "\nAll concatenated sequences had length divisible by three"

#return time to run
Time = time.time() - Start_time
print('\nTime took to run: {}'.format(Time))


