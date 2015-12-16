###########################################################################
#
# Module: selfcontain.py
# Author: Miler T. Lee
# Date: 30 April 2008
# Version: 1.0a (unix/mac)
#
#--------------------------------------------------------------------------
#
# Copyright (c) 2008, Miler T. Lee and Junhyong Kim, University of
# Pennsylvania.  All Rights Reserved.
#
# You may not use this file except in compliance with the License
# located in the top directory of this distribution.  You may obtain
# a copy of the License at http://kim.bio.upenn.edu/software/
#
# Unless required by applicable law or agreed to in writing, this
# software is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR
# CONDITIONS OF ANY KIND, either express or implied.  See the License
# for the specific language governing permissions and limitations
# under the License.
#
###########################################################################

import sys
import random
import commands

RNAFOLD = 'RNAfold -noLP'
RNAFOLD_CHUNK_SIZE = 250 #number of seqs to pass to RNAfold at once
DEFAULT_N_CONTEXTS = 100
DEFAULT_CONTEXT_PROP = 1.0


def selfcontain(seq, n_contexts = DEFAULT_N_CONTEXTS, context_length_proportion = DEFAULT_CONTEXT_PROP, context_seq_list = None):
    """
    Given a sequence, calculates the self-containment index as described in
    Lee & Kim (2008).  Default parameters are 1000 uniform random contexts of
    length 1.0 * length(seq) on the left and the right.

    Optionally a list of user-specified context sequences can be input.  Every
    two sequences in this list are used as a pair of left/right contexts for
    one iteration of the algorithm, over n_contexts or
    length(context_seq_list)/2 iterations, whichever is smaller.  Context
    sequences are trimmed to the appropriate length as specified by
    context_length_proportion (left-hand contexts are trimmed from the 5' end,
    right-hand contexts are trimmed from the 3' end); contexts shorter than
    the specified context length are used as is without padding.

    Return is the self-containment index, which falls between 0 and 1.  If
    the input sequence is predicted to be unstructured, 0 is returned.
    """

    seq_length = len(seq)
    context_length = int(context_length_proportion * seq_length)
    struct = rna_struct(seq)
    if struct == None:
        print '0'
        return None
    trunc_struct, left, right = truncate_struct_coords(struct)
    if len(trunc_struct) == 0:
        return 0

    embedded_seqs = []

    if context_seq_list == None:
        for i in xrange(n_contexts):
            left_s = random_rna_sequence(context_length)
            right_s = random_rna_sequence(context_length)
            embedded = ''.join([left_s, seq, right_s])
            embedded_seqs.append(embedded)
    else:
        n_contexts = min(len(context_seq_list)/2, n_contexts)
        for i in xrange(n_contexts):
            left_s = context_seq_list[i*2][-1*context_length:]
            right_s = context_seq_list[i*2+1][0:context_length]
            embedded = ''.join([left_s, seq, right_s])
            embedded_seqs.append(embedded)

    embedded_structs = []
    for i in xrange(len(embedded_seqs) / RNAFOLD_CHUNK_SIZE + 1):
        seq_string = '\n'.join(embedded_seqs[i*RNAFOLD_CHUNK_SIZE:(i+1)*RNAFOLD_CHUNK_SIZE])
        output = commands.getoutput('echo "'+ seq_string + '" | ' + RNAFOLD)
        output = output.rstrip()
        embedded_structs += output.split('\n')
    match_prop_list = []

    if len(embedded_structs) < len(embedded_seqs) * 2:
        print '0.78'
        return 0

    for i in xrange(len(embedded_seqs)):
        embedded_struct = embedded_structs[i*2+1].split()[0]
        portion = embedded_struct[context_length+left:context_length+right+1]
        portion = repair_vienna(portion)
        match_prop_list.append(struct_match_proportion(portion, trunc_struct))

    return mean(match_prop_list)


def selfcontain_fasta(fasta_file, n_contexts = DEFAULT_N_CONTEXTS, context_length_proportion = DEFAULT_CONTEXT_PROP, context_seq_list_file = None):
    """
    Prints to stdout the self containment indices of each sequence in the
    fasta file, according to the input parameters.  Context sequence list is
    specified as a filename containing one context per line (None by default).
    """

    context_seq_list = None
    if context_seq_list_file != None:
        try:
            f = open(context_seq_list_file)
            context_seq_list = f.read().split('\n')
            f.close()
        except:
            print 'File ' + context_seq_list_file + ' not found'
            return

    fr = fasta_reader(fasta_file)
    if fr == None:
        print 'File ' + fasta_file + ' not found'
        return

    while fr.has_next():
        id, seq = fr.next()
       
        sc = selfcontain(seq, n_contexts, context_length_proportion, context_seq_list)
        if sc == None:
            return
        print sc
    fr.close()


def selfcontain_contextfile(seq, n_contexts = DEFAULT_N_CONTEXTS, context_length_proportion = DEFAULT_CONTEXT_PROP, context_seq_list_file = None):
    """
    Wrapper to call selfcontain using a context sequence file.
    """

    context_seq_list = None
    if context_seq_list_file != None:
        f = open(context_seq_list_file)
        context_seq_list = f.read().split('\n')
        f.close()
    return selfcontain(seq, n_contexts, context_length_proportion, context_seq_list)


def rna_struct(seq):
    """
    Given the input sequence, returns the parentheses
    structure as returned by RNAfold
    """

    output = commands.getoutput('echo '+ seq + ' | ' + RNAFOLD)

    if output.find('No such file') > -1:
        return None

    struct = output[output.find('\n')+1:output.find(' ')]
    return struct


def truncate_struct_coords(struct):
    """
    Takes a Vienna structure and trims unpaired bases
    from either end.  Returns the resulting trimmed struct, the left
    index, and the right index (inclusive) into the original structure
    corresponding to the substructure.
    """

    left = struct.find('(')
    right = struct[::-1].find(')')
    if right == 0:  #no unpaired bases on right
        truncated_struct = struct[left:]
    else:
        truncated_struct = struct[left:-1*right]

    return truncated_struct, left, len(struct) - right - 1


def random_rna_sequence(seq_length):
    """
    Returns a random string of [A, C, G, U] of length seq_length,
    uniform probability on the bases.
    """

    seq = []
    for i in xrange(seq_length):
        seq.append(random.choice(['A', 'C', 'G', 'U']))
    return ''.join(seq)


def repair_vienna(struct, substitute_char = '.'):
    """
    Takes a possibly inconsistent vienna rna structure
    and removes unmatched parentheses.  By default, these
    are replaced by '.' but optionally replaces them with '-'
    """

    ##create a stack of left parentheses' indices,
    ##pop whenever a right parenth is encountered
    paren_stack = []
    new_struct = list(struct)
    for i in xrange(len(struct)):
        if new_struct[i] == '(':
            paren_stack.append(i)
        elif new_struct[i] == ')':
            if len(paren_stack) == 0:
                new_struct[i] = substitute_char
            else:
                paren_stack.pop()
    for remaining in paren_stack:
        new_struct[remaining] = substitute_char
    return ''.join(new_struct)


def struct_match_proportion(orig_struct, new_struct):
    """
    Finds the longest stretch of matching aligned struct, divides by length
    does NOT check that the stretch is a legit structure
    """
    if len(orig_struct) == 0 or len(new_struct) == 0:
        return 0

    max_length = 0
    temp_length = 0
    for pair in zip(orig_struct, new_struct):
        if pair[0] == pair[1]:
            temp_length += 1
        else:
            if temp_length > max_length:
                max_length = temp_length
                temp_length = 0
    if temp_length > max_length:
        max_length = temp_length
    return 1.0 * max_length / len(orig_struct)


def mean(numbers):
    """
    Arithmetic mean.  Returns None if the input list is empty
    """

    if len(numbers) < 1:
        return None

    return 1.0 * sum(numbers) / len(numbers)


class fasta_reader:
    """
    Incrementally reads a fasta file.
    """

    file = None
    nextheader=''
    
    def __init__(self, filename):
        try:
            self.file = open(filename, 'r')
            # fast forward to the first entry
            while 1:
                line = self.file.readline()
                if line == '':
                    self.close()
                    return
                elif line[0] == '>':
                    self.nextheader = line[1:].rstrip()
                    return
        except IOError:
            print 'No such file: ' + filename
            sys.exit()
        
    def has_next(self):
        return len(self.nextheader) > 0

    def next(self):
        """Returns an (id, sequence) tuple, or () if file is finished"""
        #if global nextheader is empty, return empty
        #otherwise, the header is the nextheader
        try:
            identifier = self.nextheader
            total = []
            while 1:
                line = self.file.readline()
                if line == '' or line[0] == '>':  #EOF, end of entry
                    break
                total.append(line.rstrip())

            sequence = ''.join(total)

            if len(line) > 0:
                self.nextheader = line[1:].rstrip()
            else:
                self.nextheader = ''
                self.close()

            return (identifier, sequence)

        except:
            self.nextheader=''
            self.close()
            return ()

    def close(self):
        self.file.close()


def example():
    print selfcontain('UGAGGCUGAAACAUAGCAGGGCUCGCUCUGGAGGUAACGUUCCAGCUCUAGGCAGGUCCUGUUGCCAGUACUCCUCCA', n_contexts = 100, context_length_proportion = 1, context_seq_list = None)

    selfcontain_fasta('sample_files/sample.fa', n_contexts = 100, context_length_proportion = 1.0, context_seq_list_file = 'sample_files/random_seqs_200')



def main(argv = None):
    """
    Run at command line
    """

    if argv == None:
        argv = sys.argv
    
    switches = {'s':'', 'i':'', 'n': DEFAULT_N_CONTEXTS, 'p': DEFAULT_CONTEXT_PROP, 'c': None}

    for i in range(len(argv)/2):
        k = 2*i + 1   #key
        v = 2*i + 2   #value

        if argv[k] == '-s' or argv[k] == '--sequence (RNA)':
            switches['s'] = argv[v]
        elif argv[k] == '-i' or argv[k] == '-input_file':
            switches['i'] = argv[v]
        elif argv[k] == '-n' or argv[k] == '--num_contexts':
            switches['n'] = int(argv[v])
        elif argv[k] == '-p' or argv[k] == '--context_length_prop':
            switches['p'] = float(argv[v])
        elif argv[k] == '-c' or argv[k] == '--context_seq_list_file':
            switches['c'] = argv[v]
        else:
            print_help()
            sys.exit()

    if switches['i'] == '' and switches['s'] == '':
        print_help()
        sys.exit()
    elif switches['s'] != '':
        print selfcontain_contextfile(switches['s'], switches['n'], switches['p'], switches['c'])
    elif switches['i'] != '':
        selfcontain_fasta(switches['i'], switches['n'], switches['p'], switches['c'])




def print_help():
    print "Options:\n\
    -s, --sequence                  sequence\n\
    -i, --input_file                fasta sequence input file\n\
    -n, --num_contexts              number of contexts (default 100)\n\
    -p, --context_length_prop       context length as a fraction of seq length (default 1.0)\n\
    -c, --context_seq_list_file     file containing list of contexts to use (default\n\
    -h, --help                      help\n\
\n\
    if -s is specified, -i is ignored\n\
    "
    

if __name__ ==  "__main__":
    main()
