# python 1by1.py 

import sys
import getopt
from subprocess import call
import os.path

opts, extraparams = getopt.getopt(sys.argv[1:], 'i:')
for o,p in opts:
	if o == '-i':
		print p
		inPath = p

with open("../data/"+inPath, "r") as f:
	fastaLines = f.readlines()

	for i in range(1,len(fastaLines),2):
		header = fastaLines[i-1]
		seq = fastaLines[i]
		with open("../data/tmp_"+inPath, 'w') as tempOut:
			tempOut.write(header+seq)
		
		call("perl microPred.pl tmp_"+inPath, shell=True)
		
		with open("../data/selected.tmp_"+inPath+".-21.features", "r") as inFile:
			with open("../data/selected."+inPath+".-21.features", 'a') as finalOut:
				features = inFile.readline()
				if features != '':
					finalOut.write(features)

# from subprocess import call

# with open("../data/hpnapnr-nofail.labels") as l:
# 	labelLines = l.readlines()

# with open("../data/hpnapnr-nofail.features") as f:
# 	featureLines = f.readlines()

# with open("../data/mirbase_19_hairpins.fa") as f:
# 	fastaLines = f.readlines()

# features = {}

# for i in range(len(labelLines)):
# 	features[labelLines[i]] = featureLines[i]

# doStuff = False
# for i in range(1,len(fastaLines),2):
# 	header = fastaLines[i-1].split()[0]+'\n'
# 	seq = fastaLines[i]
# 	if header == '>mtr-MIR2680d\n':
# 		doStuff = True
# 	if doStuff:
# 		if header in labelLines:
# 			featvec = features[header]
# 		else:	
# 			with open("../data/temp1.fasta", 'w') as tempOut:
# 				tempOut.write(header+seq)
# 			call("perl microPred.pl temp1.fasta", shell=True)
# 			with open("../data/selected.temp1.fasta.-21.features", "r") as inFile:
# 				featvec = inFile.readline()
# 		with open("../data/mirbase_19.mirdb", 'a') as finalOut:
# 			finalOut.write(header+seq+featvec)