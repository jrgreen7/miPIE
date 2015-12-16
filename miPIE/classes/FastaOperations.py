# Removes newline characters from within sequences in a fasta format file.
# Newlines need to be removed before some programs (mostly RNALfold) can operate on a fasta file.
#
# inPath: input filename
# outPath: output filename
from subprocess import call

def remove_AU(inPath, outPath, numAUs):
	AUstring = 'AU'*numAUs
	with open(inPath, 'r') as inFile:
		with open(outPath, 'w') as outFile:
			with open(inPath+'.removed.AUs', 'w') as auFile:
				for line in inFile:
					if line[0] == '>':
						header = line
					else:
						seq = line
						if AUstring in seq:
							auFile.write(header)
							auFile.write(seq)
						else:
							outFile.write(header)
							outFile.write(seq)

def RNA_reverse_complement(inPath, outPath):
	with open(inPath, 'r') as inFile:
		with open(outPath, 'w') as outFile:
			inFile.seek(0)
			for line in inFile:
				if line[0] == '>':
					outFile.write(line)
				else:
					line = line.upper()
					line = line.replace('A','u')
					line = line.replace('U','a')
					line = line.replace('G','c')
					line = line.replace('C','g')
					line = line.upper()
					line = line[::-1].strip()
					outFile.write(line+'\n')

def convert_DNA_to_RNA(inPath, outPath):
	with open(inPath, 'r') as inFile:
		with open(outPath, 'w') as outFile:
			inFile.seek(0)
			for line in inFile:
				if line[0] == '>':
					outFile.write(line)
				else:
					line=line.upper()
					if 'U' in line:
						outFile.write(line)
					else:
						line = line.replace('A','u')
						line = line.replace('T','a')
						line = line.replace('G','c')
						line = line.replace('C','g')
						line = line.upper()
						line = line[::-1].strip()
						outFile.write(line+'\n')

def convert_RNA_to_DNA(inPath, outPath):
	with open(inPath, 'r') as inFile:
		with open(outPath, 'w') as outFile:
			inFile.seek(0)
			for line in inFile:
				if line[0] == '>':
					outFile.write(line)
				else:
					line=line.upper()
					if 'T' in line:
						outFile.write(line)
					else:
						line = line.replace('A','t')
						line = line.replace('U','a')
						line = line.replace('G','c')
						line = line.replace('C','g')
						line = line.upper()
						line = line[::-1].strip()
						outFile.write(line+'\n')

def remove_newlines(inPath, outPath):
	b = ''
	with open(inPath, 'r') as inFile:
		with open(outPath, 'w') as outFile:
			for line in inFile:
				if line[0] == '>':
					if b != '':
						outFile.write(b.strip() + '\n')
					outFile.write(line)
					b = ''
				else:
					b += line.strip()
			if b != '\n':
				outFile.write(b.strip())

# Splits a fasta format file into smaller fasta files.
# Output files will be named '<inPath>.1', '<inPath>.2', etc.
#
# inPath: input filename
# numOut: number of output files
def split_fasta(inPath, numOut):
	# Set lineCount to the number of lines in the input file
	remove_newlines(inPath, inPath+'.tmp_fixed')

	with open(inPath+'.tmp_fixed', 'r') as inFile:
		lineCount = 0
		for line in inFile:
			lineCount += 1

	with open(inPath+'.tmp_fixed', 'r') as inFile:
		for i in range(numOut-1):
			# Build path of outgoing file
			outPath = ""
			for text in inPath.split('.')[:-1]:
				outPath += text
				outPath += '.'
			outPath += str(i)
			outPath += '.'
			outPath += inPath.split('.')[-1]
			with open(outPath, 'w') as outFile:
				title = 'X'
				seq = 'X'
				for lineNum in range(lineCount/numOut/2):
					while title[0] != '>':
						title = inFile.readline()
					while seq[0] not in 'agtcunAGTCUN':
						seq = inFile.readline()
					outFile.write(title)
					outFile.write(seq)
					title = 'X'
					seq = 'X'
		# Build path of final outgoing file
		outPath = ""
		for text in inPath.split('.')[:-1]:
			outPath += text
			outPath += '.'
		outPath += str(numOut-1)
		outPath += '.'
		outPath += inPath.split('.')[-1]
		# Write end of incoming file to final outgoing file
		with open(outPath, 'w') as outFile:
			for line in inFile:
				outFile.write(line)
	call('rm '+inPath+'.tmp_fixed', shell=True)

def merge_fasta(outPath, numIn):
	with open(outPath, 'w') as outFile:
		inPath = outPath.rsplit('.',1)[0]
		inExt = outPath.rsplit('.',1)[1]
		for i in range(numIn):
			with open(inPath+'.'+str(i)+'.'+inExt, 'r') as inFile:
				for line in inFile:
					outFile.write(line)
			# call('rm '+outPath+'.'+str(i), shell=True)