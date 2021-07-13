#!/usr/bin/env python


from sys import argv

infileName = argv[1]
pairingMethod = argv[2]


#read the file
with open(infileName, 'r') as infile:
	lineList = [l.strip('\n') for l in infile]


#print out all
if pairingMethod == 'permute':
	pairList = []
	for i in lineList:
		for j in lineList:
			if i!=j:
				print('{}\t{}'.format(i,j))
else:
	print('Error, currently no other method than "permute"')


