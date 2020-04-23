#   Author: Ryan Petit
#   CECS 660 - 01
#   Fitch and Margoliash Algorithm Implementation
#   23 April 2020

import numpy as np

class FitchMargoliash:

    def __init__(self, distance_table, point_dictionary):
        super().__init__()
        self.distance_table = distance_table
        self.point_dictionary = point_dictionary

    def run(self):
        distance_table = self.distance_table
        point_dictionary = self.point_dictionary

        handle = ""

        z = 0
        r = 0.0
        case_1 = False
        while(True):
            min_distance, min_i, min_j = self.__minimum_distance(distance_table)
            # print('1: {} -> {} = {:.2f}'.format(point_dictionary[min_i],point_dictionary[min_j],min_distance))

            # find average distance between min_i and rest of points
            denominator = 0
            sum_distance_i = 0
            sum_distance_j = 0
            cluster = ''

            for i in range(0,len(distance_table)):
                if i is min_i:
                    continue
                elif i is min_j:
                    continue
                else:
                    denominator += 1

                    if distance_table[min_i][i] is 'X':
                        sum_distance_i += float(distance_table[i][min_i])
                    else: 
                        sum_distance_i += float(distance_table[min_i][i])
                    if distance_table[min_j][i] is 'X':
                        sum_distance_j += float(distance_table[i][min_j])
                    else: 
                        sum_distance_j += float(distance_table[min_j][i])
                    
                    cluster += point_dictionary[i]
            
            # HANDLE CODE
        
            if denominator is 0: 
                is_cluster = True
                if case_1: z = 0.0
                # add current z value here
                _handle = handle.rsplit(')',1)
                if _handle[1] is "":
                    for p in point_dictionary:
                        if p[0] in _handle[0]:
                            continue 
                        else: 
                            _handle[1] = p
                            is_cluster = False
                            break
                if is_cluster:
                    # print('End 1',_handle)
                    if _handle[0].endswith(')'): handle = "{}:{});".format(_handle[0],z)
                    else: handle = "{});".format(_handle[0],z)

                else: 
                    # print('End 2',_handle)
                    handle = "({}),{}:{});".format(_handle[0],_handle[1],z)
                break
            combinations = []
            combinations.append([[1,1,0],[1,0,1],[0,1,1]])         

            # print('a1',min_distance)
            # print('a2',sum_distance_i/denominator)
            # print('a3',sum_distance_j/denominator)

            a = np.array(combinations[0])
            b = np.array([min_distance,sum_distance_i/denominator,sum_distance_j/denominator])
            x,y,z = np.linalg.solve(a,b)
            # print(x,y,z)

            
            if x < 0 or y < 0: 
                print("ERROR - Distance table has not been normalized")

            
            # if we are not at the very beginning
            if handle != "":
                if point_dictionary[min_i][0] not in handle and point_dictionary[min_j][0] not in handle:
                    # print('Case 1')
                    handle = "({}:{},({}:{},{}:{}))".format(handle,z, point_dictionary[min_i],x,point_dictionary[min_j],y)
                    case_1 = True
                else:
                    # print('Case 2')
                    handle = "({}:{},{}:{})".format(handle,y,point_dictionary[min_i],x)
                    case_1 = False
            else: 
                # if nodes being added are not already in the handle
                handle = "({}:{},{}:{})".format(point_dictionary[min_i],x,point_dictionary[min_j],y)
            # print('>     ',handle)
            # print('2: {} -> {} = {:.2f}'.format(point_dictionary[min_i],cluster,sum_distance_i/denominator))
            # print('3: {} -> {} = {:.2f}'.format(point_dictionary[min_j],cluster,sum_distance_j/denominator))

            cluster = '{}{}'.format(point_dictionary[min_i],point_dictionary[min_j])

            new_table = distance_table
            new_point_dictionary = []
            x_lst = []
            for i in range(0,len(distance_table)):
                if i is min_i or i is min_j:
                    continue
                
                if distance_table[i][min_i] is 'X':
                    a = distance_table[min_i][i]
                else: a = distance_table[i][min_i]
                if distance_table[i][min_j] is 'X':
                    b = distance_table[min_j][i]
                else: b = distance_table[i][min_j]            
                x = (float(a) + float(b))/2 - min_distance/2
                x_lst.append(str(x))
                print('4: {} -> {} = {:.2f}'.format(point_dictionary[i],cluster,x))

            # Update point dictionary for next pass
            if min_i > min_j:
                point_dictionary.pop(point_dictionary.index(point_dictionary[min_i]))
                point_dictionary.pop(point_dictionary.index(point_dictionary[min_j]))
                new_table.pop(min_i)
                new_table.pop(min_j)
                for col in new_table:
                    col.pop(min_i)
                    col.pop(min_j)
            else:
                point_dictionary.pop(point_dictionary.index(point_dictionary[min_j]))
                point_dictionary.pop(point_dictionary.index(point_dictionary[min_i]))
                new_table.pop(min_j)
                new_table.pop(min_i)
                for col in new_table:
                    col.pop(min_j)
                    col.pop(min_i)

            point_dictionary.append(cluster)

            # Update data table for next pass
            for row in new_table:
                row.append(x_lst[0])
                x_lst.pop(0)

            new_table.append(['X' for i in range(len(new_table[0]))])
            distance_table = new_table
            # print('####################################################')
        return handle

    def __minimum_distance(self, hamming_table):
        min_i = 0
        min_j = len(hamming_table)-1
        min_distance = float(hamming_table[min_i][min_j])

        for i in range(0,len(hamming_table)):
            for j in range(0,len(hamming_table)):
                if hamming_table[i][j] is 'X':
                    continue
                distance = float(hamming_table[i][j])
                if distance < min_distance:
                    min_i = i
                    min_j = j
                    min_distance = distance
        return min_distance, min_i, min_j

