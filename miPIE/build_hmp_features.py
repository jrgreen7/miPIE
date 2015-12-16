import threading
import sys
import getopt
from subprocess import call
import os.path
import classes.FastaOperations as FastaOps
import os

patternClass = ''

# ##########################################################
# Parameters:
# -i input fasta file. File must be present in the data directory.
# -n Number of processor threads to be used by the program
# -c Class of data (either 'miRNA' or 'pseudo')
#
############################################################
opts, extraparams = getopt.getopt(sys.argv[1:], 'i:n:c:')
for o,p in opts:
	if o == '-i':
		inPath = p
	if o == '-n':
		numThreads = int(p)
	if o == '-c':
		patternClass = p

class myThread(threading.Thread):
	def __init__(self, inPath):
		threading.Thread.__init__(self)
		self.inPath = inPath

	def run(self):
		with open(self.inPath, "r") as f:
			fastaLines = f.readlines()

			for i in range(1,len(fastaLines),2):
				header = fastaLines[i-1]
				seq = fastaLines[i].strip()
				with open("tmp_"+self.inPath, 'w') as tempOut:
					tempOut.write(header+seq)
				
				call("perl HeteroMirPred.pl tmp_"+self.inPath, shell=True)
				
				with open("tmp_"+self.inPath+".csv", "r") as inFile:
					with open(self.inPath+".csv", 'a') as finalOut:
						inFile.readline()
						features = inFile.readline()
						if "Other" in features:
							exit()
						if features != '':
							finalOut.write(header)
							finalOut.write(features)
		call('rm *.ps', shell=True)

###########################################
# Main script starts here
###########################################

# Move data file into microPred directory
call('cp data/'+inPath+' progs/HeteroMirPred/'+inPath, shell=True)
exitFlag = 0

# Split input fasta into <numThreads> smaller fasta files
FastaOps.split_fasta('progs/HeteroMirPred/'+inPath, numThreads)

#Get a path to the input file that doesn't have '.fasta' on the end
threadPath = inPath.rsplit('.',1)[0]
threadExt = inPath.rsplit('.',1)[1]
# Move to microPred directory and run microPred on all the smaller files (multi-threaded)
os.chdir('progs/HeteroMirPred/')
threads = []
for i in range(numThreads):
	threads.append(myThread(threadPath+'.'+str(i)+'.'+threadExt))

for thread in threads:
	thread.start()

for thread in threads:
	thread.join()

os.chdir('../../')

# 
with open('progs/HeteroMirPred/'+inPath+'.csv', 'w') as finalOut:
	for i in range(numThreads):
		with open('progs/HeteroMirPred/'+threadPath+'.'+str(i)+'.'+threadExt+'.csv', 'r') as finalIn:
			for line in finalIn:
				finalOut.write(line)
		# call('rm progs/HeteroMirPred/'+threadPath+'.'+str(i)+'.'+threadExt+'.csv', shell=True)

with open('data/'+inPath+'.hmp', 'w') as newOut:
	with open('data/'+inPath+'.headers', 'w') as headerOut:
		with open('progs/HeteroMirPred/'+inPath+'.csv', 'r') as oldOut:
			for line in oldOut:
				if line[0] in '1234567890-':
					newOut.write(line)
				else:
					headerOut.write(line)

outPath20 = 'data/' + inPath + '.hmp20'

featureIndices = [1,3,9,10,20,26,27,32,33,34,35,39,37,43,44,210,153,156,169,203]

with open('data/'+inPath+'.hmp', 'r') as inFile:
	with open(outPath20, 'w') as outFile:
		for line in inFile:
			features = line.split(',')
			for index in featureIndices:
				outFile.write(features[index]+',')
			if patternClass == '':
				outFile.write(features[-1])
			else:
				outFile.write(patternClass+"\n")
call('rm progs/HeteroMirPred/*.csv', shell=True)
call('rm progs/HeteroMirPred/*.fasta', shell=True)
call('rm progs/HeteroMirPred/*.fa', shell=True)
call('rm progs/HeteroMirPred/tmp_*', shell=True)
