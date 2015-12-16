from subprocess import call
# Local helper classes for fasta manipulation, feature set manipulation
import classes.FastaOperations as FastaOps
from classes.FeatureSet import *
# SKLearn classes
import numpy as np
from sklearn.ensemble import *
from sklearn.externals import joblib
# System classes for argument-passing
import sys
import getopt

opts, extraparams = getopt.getopt(sys.argv[1:], 'd:')
for o,p in opts:
	if o == '-d':
		mirdeepDir = p

# Chop the trailing / off of the directory, if necessary
if mirdeepDir[-1] == '/':
	mirdeepDir = mirdeepDir[:-1]

runtimeDir = mirdeepDir.split('/')[-1]
# Get paths for all input files from miRDeep2 directory structure
mdOutputPath = mirdeepDir+"/output.mrd"
mdStrPath = mirdeepDir+'/tmp/precursors.str'
mdFaPath = mirdeepDir+'/tmp/precursors.fa'
mdRandPath = mirdeepDir+'/tmp/precursors_for_randfold.rand'
expFeatsPath = mirdeepDir+'/tmp/expressionFeatures.csv'

# miRDeep spits out DNA sequences with odd mixes of upper and lower case. Fix those.
def convert_DNA_to_RNA(seq):
	seq = seq.upper()
	seq = seq.replace('T','U')
	return seq

# Extract ids of pre-processed precursors from the output.mrd file
ids = []
for line in open(mdOutputPath):
	if line[0] == ">":
		ids.append(line.strip()[1:])

# Extract sequences of pre-processed precursors from the precursors.fa file
seqs = {}
getSeq = False
for line in open(mdFaPath):
	if getSeq == True:
		seq = line.strip()
		seq = convert_DNA_to_RNA(seq)
		seqs[seqId] = seq
		getSeq = False
	for i in ids:
		if i in line:
			seqId = i
			getSeq = True

expFeats = {}
headerLine = True
for i, line in enumerate(open(expFeatsPath).readlines()[1:]):
	features = line.strip().split(',')
	expFeats[ids[i]] = features



with open("data/tmp.fa",'w') as outFile:
	for key in seqs.keys():
		outFile.write(">"+key+'\n')
		outFile.write(seqs[key]+'\n')

# call("python build_hmp_features.py -i tmp.fa -n 8", shell=True)

hmpHeaders = [l.strip()[1:] for l in open("data/tmp.fa.headers").readlines()]
hmpFeats = {}
headerNum = 0
for line in open("data/tmp.fa.hmp"):
	hmpFeats[hmpHeaders[headerNum]] = line.strip().split(',')[:-1]
	headerNum += 1

finalIds = []

for i in ids:
	if i in hmpHeaders:
		finalIds.append(i)


with open("data/tmp.csv",'w') as csvOut:
	featNames = []
	for i in range(len(hmpFeats[hmpHeaders[0]])+len(expFeats[ids[0]])):
		featNames.append("'feat"+str(i)+"'")
	csvOut.write(",".join(featNames) + ',"class"\n')
	for i in finalIds:
		csvOut.write(','.join(hmpFeats[i])+','+','.join(expFeats[i])+',miRNA\n')

fs = FeatureSet()
fs.load("data/tmp.csv", patternClass = "miRNA")
fs.select_features([12, 18, 24, 26, 38, 43, 64, 139, 163, 174, 176, 190, 192, 201, 202, 206, 208, 215, 218, 219])
fs.libsvm_scale()
fs.export("data/tmp_select.csv")

################################################################
# Use the RF model to predict on the new data
################################################################
rf = joblib.load("model/smirpdeep_model.pkl")

testLines = open("data/tmp_select.csv", 'r').readlines()[1:]
testData = np.array([[float(y) for y in x.split(',')[:-1]] for x in testLines])

allPredictions = rf.predict_proba(testData)
allPredictions = sorted(allPredictions, key = lambda x : x[0])
        
with open(mirdeepDir+"/OURMETHOD.results", 'w') as out:
    for i, p in enumerate(allPredictions):
        out.write(finalIds[i]+"\t"+str(p[0]) + "\t" + str(p[1])+ '\n')