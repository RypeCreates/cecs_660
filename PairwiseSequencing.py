#   Author: Ryan Petit
#   CECS 660 - 01
#   FASTA and Dynamic Programming (DP) Classes
#   23 April 2020


import numpy as np
import os
import re

class FASTA:

    def __init__(self, fileName=None, directory=None, header=None,sequence=None):
        super().__init__()

        if fileName is not None and directory is not None:
            self.fileName = str(fileName).split('.')[0]
            self.header,self.sequence = self.__read_sequence(fileName,directory)
        elif header is not None and sequence is not None:
            self.header = header
            self.sequence = sequence
        pattern = ">(.*?) "
        self.name = re.search(pattern,self.header).group(1)

    def __read_sequence(self,fileName, directory):
        with open ("{}/{}".format(directory,fileName), "r") as file:
            data=file.readlines()
        header = data[0].replace("\n","")
        seq = ""
        for i in range(1,len(data)):
            seq += data[i].replace("\n","").replace(" ","")
        return header,seq

class DP:

    # constructor
    def __init__(self,process_id,alignment_type,file1,file2,directory,sequence_type=None,match=None,mismatch=None,gap=None,aa_score_name=None):

        self.process_id = process_id
        self.alignment_type = alignment_type # either g or l
        self.seq1 = FASTA(fileName=file1,directory=directory)
        self.seq2 = FASTA(fileName=file2,directory=directory)

        if sequence_type is not None:
            self.sequence_type = sequence_type
        else: 

            self.sequence_type = "dna"
        if match is not None:
            self.match = float(match)
        else: 
            self.match = float(5)
        if mismatch is not None:
            self.mismatch = float(mismatch)
        else:
            self.mismatch = float(-3)
        if gap is not None:
            self.gap = float(gap)
        else:
            self.gap = float(-4)    
        if aa_score_name is not None:
            self.aa_score_name = aa_score_name
        else:
            self.aa_score_name = "EPAM250"
        if self.sequence_type == "aa":
            self.aa_score_sheet = self.__initialize_scoring_matrix(self.aa_score_name)

    # DNA and global methods
    def __initialize_scoring_matrix(self,matrix_name):
        with open ("Scoring Matrices/{}".format(matrix_name), "r") as file:
            data=file.readlines()

        values = []
        for row in data:
            if row[0] is '#':
                continue
            else:
                values.append(row.split())

        values[0].insert(0,"")
        return values

    def initialize_matrices(self):  

        # add space for gap logic
        self.seq1.sequence = " "+self.seq1.sequence
        self.seq2.sequence = " "+self.seq2.sequence

        matrix_s = np.array([[0 for i in range(len(self.seq1.sequence))] for j in range(len(self.seq2.sequence))]) # scoring matrix
        matrix_d = np.array([[0 for i in range(len(self.seq1.sequence))] for j in range(len(self.seq2.sequence))]) # direction matrix      

        if self.alignment_type is "g":
            # initialize gap penalties going down column 0
            current_mismatch = 0
            for cell in matrix_s:
                cell[0] = current_mismatch
                current_mismatch += self.gap

            # initialize gap penalties going across row 0
            current_mismatch = 0
            for i in range(0, len(matrix_s[0])):
                matrix_s[0][i] = current_mismatch
                current_mismatch += self.gap

        return matrix_s, matrix_d
    
    def score(self, matrix_s, matrix_d):
        
        sequence_1 = self.seq1.sequence
        sequence_2 = self.seq2.sequence

        for i in range(1,len(sequence_2)):            # for each row
            for j in range(1,len(sequence_1)):        # for each column
                
                if self.alignment_type is "l":
                    diag_score =    max(0,matrix_s.item(i-1,j-1) + self.__getMatch(sequence_2[i],sequence_1[j])) 
                    back_score =    max(0,matrix_s.item(i,j-1) + self.gap)
                    up_score =      max(0,matrix_s.item(i-1,j) + self.gap)
                elif self.alignment_type is "g":
                    diag_score = matrix_s.item(i-1,j-1) + self.__getMatch(sequence_2[i],sequence_1[j])
                    back_score = matrix_s.item(i,j-1) + self.gap
                    up_score = matrix_s.item(i-1,j) + self.gap         
                else: 
                    print("Parameter error in alignment")
                    return -1,-1           

                diag =  [diag_score,3] # match/mismatch
                back =  [back_score,2] # gap in sequence 2
                up   =  [up_score,1] # gap in sequence 1 
                
                moveset = [diag,up,back]                                # rearrange to prioritize gaps over mismatches
                optimal = moveset[0]

                for x in range(1,len(moveset)):
                    if optimal[0] < moveset[x][0]: optimal = moveset[x]

                matrix_s[i][j] = optimal[0]
                matrix_d[i][j] = optimal[1]
        return matrix_s, matrix_d

    def __getMatch(self, cell1, cell2):
        if cell1 is cell2: return self.match
        else: return self.mismatch

    def stacktrace(self, matrix_s,matrix_d):
        mismatched_position_count = 0.0

        sequence_1 = self.seq1.sequence
        sequence_2 = self.seq2.sequence

        s1_stack = []
        s2_stack = []
        align_stack = []

        # Procedure for Local
        if self.alignment_type == "l":
            max_value = 0
            max_coords = [0,0]
            for row in range(0,len(matrix_s)):
                for col in range(0,len(matrix_s[0])):
                    if matrix_s[row][col] > max_value:
                        max_value = matrix_s[row][col]
                        max_coords = [row,col]
            i = max_coords[0]
            j = max_coords[1]

            while(i is not 0 or j is not 0):
                if matrix_s.item(i,j) is 0:
                    break
                if matrix_d.item(i,j) is 3:
                    s1_stack.append(sequence_1[j])
                    s2_stack.append(sequence_2[i])
                    if sequence_1[j] is sequence_2[i]:
                        align_stack.append("|")
                    elif self.sequence_type == "aa": 

                        score = self.__getMatch_aa(sequence_1[j],sequence_2[i])
                        if score > 0: 
                            mismatched_position_count += 0.2
                            align_stack.append(":")
                            # print(':')
                        elif score == 0: 
                            mismatched_position_count += 0.5
                            align_stack.append(".")
                            # print('.')
                        else: 
                            mismatched_position_count += 1
                            align_stack.append(" ")
                    else: 
                        align_stack.append(" ")
                        mismatched_position_count += 1.0
                    i = i-1
                    j = j-1
                elif matrix_d.item(i,j) is 2:
                    s1_stack.append(sequence_1[j])
                    s2_stack.append("_")
                    align_stack.append(" ")
                    mismatched_position_count += 1
                    j = j-1
                elif matrix_d.item(i,j) is 1:
                    s1_stack.append("_")
                    s2_stack.append(sequence_2[i])
                    align_stack.append(" ")
                    mismatched_position_count += 1
                    i = i-1
                else:
                    print("End of Traceback.")
                    break

        elif self.alignment_type == "g":
            i = len(sequence_2)-1
            j = len(sequence_1)-1            
            
            while(i is not 0 or j is not 0):
                
                if matrix_d.item(i,j) is 3:
                    s1_stack.append(sequence_1[j])
                    s2_stack.append(sequence_2[i])
                    if sequence_1[j] == sequence_2[i]:
                        align_stack.append("|")
                    elif self.sequence_type == "aa": 

                        score = self.__getMatch_aa(sequence_1[j],sequence_2[i])
                        if score > 0: 
                            mismatched_position_count += 0.2
                            align_stack.append(":")
                        elif score == 0: 
                            mismatched_position_count += 0.5
                            align_stack.append(".")
                        else: 
                            mismatched_position_count += 1
                            align_stack.append(" ")
                    else: 
                        align_stack.append(" ")
                    i = i-1
                    j = j-1
                elif matrix_d.item(i,j) is 2 or (i == 0):
                    mismatched_position_count += 1
                    s1_stack.append(sequence_1[j])
                    s2_stack.append("_")
                    align_stack.append(" ")
                    j = j-1
                elif matrix_d.item(i,j) is 1 or (j == 0):
                    mismatched_position_count += 1
                    s1_stack.append("_")
                    s2_stack.append(sequence_2[i])
                    align_stack.append(" ")
                    i = i-1
                elif matrix_d.item(i,j) is 0:
                    print("Something went wrong in traceback.")
                    break
        stack_size = len(align_stack)       
        s1,al,s2,cushion = "","","",""
        print_index = 0

        # Create new file here
        fileName = '{}_align_{}'.format(self.seq1.fileName,self.seq2.fileName)
        with open('ProcessSummaries/{}/PairwiseAlignments/{}'.format(self.process_id,fileName),'w+') as file:
            # Create file header
            file.write('Process ID: {}'.format(self.process_id))
            file.write('\nSequence 1: {}'.format(self.seq1.header.replace('>','')))
            file.write('\nSequence 2: {}'.format(self.seq2.header.replace('>','')))
            file.write('\nAlignment Strategy: {}'.format(self.alignment_type))
            file.write('\n\tMatch: {}'.format(self.match))
            file.write('\n\tMismatch: {}'.format(self.mismatch))
            file.write('\n\tGap: {}\n'.format(self.gap))
            file.write('\nHamming Distance: {}\n'.format(mismatched_position_count))
            
            while len(s1_stack) > 0:
                s1 += s1_stack.pop()+" "
                al += align_stack.pop()+" "
                s2 += s2_stack.pop()+" "
                cushion += "##"
                
                print_index += 1
                if print_index is 40:
                    file.write('\n{}'.format(s1))
                    file.write('\n{}'.format(al))
                    file.write('\n{}'.format(s2))
                    file.write('\n{}'.format(cushion))
                    s1,al,s2,cushion = "","","",""
                    print_index = 0
            file.write('\n{}'.format(s1))
            file.write('\n{}'.format(al))
            file.write('\n{}'.format(s2))
            file.write('\n{}'.format(cushion))

        d = self.jukes_cantor_distance(stack_size,mismatched_position_count)
        return d

    def initialize_matrices_aa(self):  
        
        score_sheet = self.aa_score_sheet
        self.gap = self.__getMatch_aa('*','A')

        # add space for gap logic
        self.seq1.sequence = " "+self.seq1.sequence
        self.seq2.sequence = " "+self.seq2.sequence

        matrix_s = np.array([[0 for i in range(len(self.seq1.sequence))] for j in range(len(self.seq2.sequence))]) # scoring matrix
        matrix_d = np.array([[0 for i in range(len(self.seq1.sequence))] for j in range(len(self.seq2.sequence))]) # direction matrix      

        if self.alignment_type is "g":
            # initialize gap penalties going down column 0
            current_mismatch = 0
            for cell in matrix_s:
                cell[0] = current_mismatch
                current_mismatch += self.gap

            # initialize gap penalties going across row 0
            current_mismatch = 0
            for i in range(0, len(matrix_s[0])):
                matrix_s[0][i] = current_mismatch
                current_mismatch += self.gap

        return matrix_s, matrix_d

    def score_aa(self, matrix_s, matrix_d):
        score_sheet = self.aa_score_sheet

        sequence_1 = self.seq1.sequence
        sequence_2 = self.seq2.sequence

        for i in range(1,len(sequence_2)):            # for each row
            for j in range(1,len(sequence_1)):        # for each column
                
                if self.alignment_type is "l":
                    diag_score =    max(0,matrix_s.item(i-1,j-1) + self.__getMatch_aa(sequence_2[i],sequence_1[j])) 
                    back_score =    max(0,matrix_s.item(i,j-1) + self.__getMatch_aa('*',sequence_1[j]))
                    up_score =      max(0,matrix_s.item(i-1,j) + self.__getMatch_aa(sequence_2[i],'*'))
                elif self.alignment_type is "g":
                    diag_score = matrix_s.item(i-1,j-1) + self.__getMatch_aa(sequence_2[i],sequence_1[j])
                    back_score = matrix_s.item(i,j-1) + self.__getMatch_aa('*',sequence_1[j])
                    up_score = matrix_s.item(i-1,j) + self.__getMatch_aa(sequence_2[i],'*')
                else: 
                    print("Parameter error in alignment")
                    return -1,-1           

                diag =  [diag_score,3] # match/mismatch
                back =  [back_score,2] # gap in sequence 2
                up   =  [up_score,1] # gap in sequence 1 
                
                moveset = [diag,up,back]                                # rearrange to prioritize gaps over mismatches
                optimal = moveset[0]

                for x in range(1,len(moveset)):
                    if optimal[0] < moveset[x][0]: optimal = moveset[x]

                matrix_s[i][j] = optimal[0]
                matrix_d[i][j] = optimal[1]

        return matrix_s, matrix_d

    def __getMatch_aa(self, cell1, cell2):
        score_sheet = self.aa_score_sheet        
        i1 = score_sheet[0].index(cell1)

        for row in score_sheet:
            if cell2 is row[0]:
                return int(row[i1])    

        return int(score_sheet[1][len(score_sheet[1]-1)]) # else return gap value


    def jukes_cantor_distance(self, stack_size, mismatch_count):
        p = mismatch_count/stack_size
        d = (-3/4)*np.log(1-(4/3)*p)
        return d