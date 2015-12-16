import sys
import getopt

opts, extraparams = getopt.getopt(sys.argv[1:], 'i:')
for o,p in opts:
	if o == '-i':
		inPath = p

outPath = '20.' + inPath

featureIndices = [1,3,9,10,20,26,27,32,33,34,35,39,37,43,44,210,153,156,169,203]

with open(inPath, 'r') as inFile:
	with open(outPath, 'w') as outFile:
		for line in inFile:
			features = line.split(',')
			for index in featureIndices:
				outFile.write(features[index]+',')
			outFile.write(features[-1])