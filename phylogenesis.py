#   Author: Ryan Petit
#   CECS 660 - 01
#   Phylogenetic Tree Constructor
#   22 April 2020

# Import BioPython Modules
from Bio import Phylo

# Import from Custom Modules
from PairwiseSequencing import FASTA
from PairwiseSequencing import DP
from TreeConstruction import FitchMargoliash

# Import Helper Modules
from datetime import datetime
import matplotlib
import pylab
import os

# GLOBAL VARIABLES
ALIGNMENT_TYPE = 'l'
SEQUENCE_TYPE = 'dna'
SCORING_MATRIX = 'EPAM250'
GAP = -4
MISMATCH = -3
MATCH = 5
SEQUENCE_DIRECTORY = 'Sequences'
TREE_NAME = 'Unrooted_Phylogram'

# Read initialization file 
def process_initialization_file():
    global ALIGNMENT_TYPE
    global SEQUENCE_TYPE
    global GAP
    global MISMATCH
    global MATCH
    global SEQUENCE_DIRECTORY
    global TREE_NAME

    with open ("{}".format('initialization.txt'), "r") as file:
            data=file.readlines()

    for line in data:
        if line.split(':')[0].strip() == 'ALIGNMENT_TYPE' and line.split(':')[1].strip() is not None:
            ALIGNMENT_TYPE = str(line.split(':')[1].strip())
        if line.split(':')[0].strip() == 'SEQUENCE_TYPE' and line.split(':')[1].strip() is not None:
            SEQUENCE_TYPE = str(line.split(':')[1].strip())
        if line.split(':')[0].strip() == 'GAP' and line.split(':')[1].strip() is not None:
            GAP = int(line.split(':')[1].strip())
        if line.split(':')[0].strip() == 'MISMATCH' and line.split(':')[1].strip() is not None:
            MISMATCH = int(line.split(':')[1].strip())
        if line.split(':')[0].strip() == 'MATCH' and line.split(':')[1].strip() is not None:
            MATCH = int(line.split(':')[1].strip())
        if line.split(':')[0].strip() == 'SEQUENCE_DIRECTORY' and line.split(':')[1].strip() is not None:
            SEQUENCE_DIRECTORY = str(line.split(':')[1].strip())
        if line.split(':')[0].strip() == 'TREE_NAME' and line.split(':')[1].strip() is not None:
            TREE_NAME = str(line.split(':')[1].strip())
    
    return

process_initialization_file()

# Initialize Process Directories
now = str(datetime.now()).replace(' ','_').split('.')[0]
process_id = "P_{}".format(now)
os.mkdir('ProcessSummaries/{}'.format(process_id))
os.mkdir('ProcessSummaries/{}/PairwiseAlignments'.format(process_id))
os.mkdir('ProcessSummaries/{}/PhyloXML'.format(process_id))

# Input Sequence Processing
#   -   read FASTA files from directory
#   -   for each file, parse the sequence for header information and store
sequences = []
files = os.listdir(SEQUENCE_DIRECTORY)
for f in files:
    seq = FASTA(fileName="{}".format(f),directory=SEQUENCE_DIRECTORY)
    sequences.append(seq)

# Initialize Hamming Table
#   -   Hamming distance definition: number of mismatched positions in alignment
hamming_table = [ [ 'X' for i in range(len(files)) ] for j in range(len(files)) ] 
index = 0
for i in range(0,len(hamming_table)):
    for j in range(0,len(hamming_table)):
        if j > index:
            hamming_table[i][j] = ' '
    index += 1 

# Populate distance table using pairwise alignment
#
# if DNA Sequences:
if SEQUENCE_TYPE == 'dna':
    for i in range(0,len(hamming_table)):
        for j in range(0,len(hamming_table)):
            if hamming_table[i][j] is ' ':
                dp = DP(process_id=process_id,alignment_type=ALIGNMENT_TYPE,directory=SEQUENCE_DIRECTORY,file1=files[i],file2=files[j],sequence_type='dna',match=MATCH,mismatch=MISMATCH,gap=GAP,aa_score_name=SCORING_MATRIX)
                matrix_s, matrix_d = dp.initialize_matrices()
                matrix_s, matrix_d = dp.score(matrix_s=matrix_s,matrix_d=matrix_d)
                distance = dp.stacktrace(matrix_s=matrix_s,matrix_d=matrix_d)
                hamming_table[i][j] = str(distance)
else: 
    for i in range(0,len(hamming_table)):
        for j in range(0,len(hamming_table)):
            if hamming_table[i][j] is ' ':
                dp = DP(process_id=process_id,alignment_type=ALIGNMENT_TYPE,directory=SEQUENCE_DIRECTORY,file1=files[i],file2=files[j],sequence_type=SEQUENCE_TYPE,match=MATCH,mismatch=MISMATCH,gap=GAP,aa_score_name=SCORING_MATRIX)
                matrix_s, matrix_d = dp.initialize_matrices_aa()
                matrix_s, matrix_d = dp.score_aa(matrix_s=matrix_s,matrix_d=matrix_d)
                distance = dp.stacktrace(matrix_s=matrix_s,matrix_d=matrix_d)
                hamming_table[i][j] = str(distance)

point_dictionary = []
file_codes = []
alphabet = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
for i in range(0,len(hamming_table)):
    print("{} | {}".format(files[i],hamming_table[i]))
    point_dictionary.append(alphabet[i])

# Fitch and Margoliash Method to build tree structure in newick format
#
point_dict_copy = point_dictionary.copy()
fm = FitchMargoliash(hamming_table,point_dictionary)
handle = fm.run()
print(handle)
for i in range(len(files)):
    print(files[i],'  ',point_dict_copy[i])

# Phylogenetic Tree construction using phyloXML file (unrooted phylogram)

from io import StringIO

handle = StringIO(handle)
tree = Phylo.read(handle,'newick')
tree.name = TREE_NAME
tree.id = process_id
tree.ladderize()

def tabulate_names(tree):
    names = {}
    for idx, clade in enumerate(tree.find_clades()):
        if clade.name:
            new_name = sequences[point_dict_copy.index(clade.name)].name
            clade.name = '%d_%s' % (idx, new_name)
        else:
            clade.name = "{}_inner".format(idx)
        names[clade.name] = clade
    return names

tabulate_names(tree)
Phylo.draw(tree)

# Save the annotated phyloXML file
export = tree.as_phyloxml()
Phylo.write(export, 'ProcessSummaries/{}/PhyloXML/tree_{}.xml'.format(process_id,process_id), 'phyloxml')