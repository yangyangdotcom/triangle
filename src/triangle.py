from hashlib import sha3_224
from operator import ge
from sys import getsizeof
import numpy as np
import pandas as pd
import math
from biopandas.pdb import PandasPdb
from itertools import combinations
from scipy.spatial import distance
import csv
ppdb = PandasPdb().fetch_pdb('3eiy')
pd.options.mode.chained_assignment = None
def calculate_dist(atom):

    atom = atom.reset_index(drop=True)
    atom_dis_dict = {}
    for i in range(len(atom.index)):
        for j in range(i+1,len(atom.index)):
            if(i != j):
                x1 = atom['x_coord'][i]
                y1 = atom['y_coord'][i]
                z1 = atom['z_coord'][i]
                a = (x1,y1,z1)

                #bundle the coord of the 2nd atom
                x2 = atom['x_coord'][j]
                y2 = atom['y_coord'][j]
                z2 = atom['z_coord'][j]
                b = (x2,y2,z2)

                atom_number_list = (atom.at[atom.index[i],'atom_number'],atom.at[atom.index[j],'atom_number'])
                atom_number_tuple = tuple(atom_number_list)
                #calculate the euclidean distance
                dist = distance.euclidean(a,b)
                atom_dis_dict[atom_number_tuple] = dist

    return atom_dis_dict

def cosine_rule(A, B, C, decider):
    if(decider == 0):
        numerator = A*A + B*B - C*C
        denominator = 2*A*B
        #return angle
        return math.degrees(math.acos(numerator/denominator))
    elif(decider == 1):
        cosC = math.cos(C)
        tmp = A*A + B*B - 2*A*B*cosC
        #return side
        return math.sqrt(tmp)

def length_angle_df(atom):
    atom_group = atom.groupby(['residue_number'])
    complete_final_df = pd.DataFrame()
    entire_atom_dis_dict = {}
    for atom, atom_df in atom_group:
        atom_number = atom_df['atom_number'].tolist()
        atom_dis_dict = calculate_dist(atom_df)
        tmp = combinations(atom_number,3)
        list_of_atom_combi = list(tmp)

        final_df = pd.DataFrame(list_of_atom_combi, columns=["atom1", "atom2", "atom3"])
        final_df['atom_name'] = atom
        final_df['D(1,2)'] = ""
        final_df['D(2,3)'] = ""
        final_df['D(1,3)'] = ""
        final_df['A(1)'] = ""
        final_df['A(2)'] = ""
        final_df['A(3)'] = ""
        
        for i in range(len(final_df.index)):
            a1 = (final_df["atom1"][i], final_df["atom2"][i])
            a2 = (final_df["atom2"][i], final_df["atom3"][i])
            a3 = (final_df["atom1"][i], final_df["atom3"][i])

            final_df.at[i, 'D(1,2)'] = atom_dis_dict[a1]
            final_df.at[i, 'D(2,3)'] = atom_dis_dict[a2]
            final_df.at[i, 'D(1,3)'] = atom_dis_dict[a3]

            d12 = float(final_df.loc[i, 'D(1,2)'])
            d13 = float(final_df.loc[i, 'D(1,3)'])
            d23 = float(final_df.loc[i, 'D(2,3)'])

            final_df.at[i, 'A(1)'] = cosine_rule(d12,d13,d23,0)
            final_df.at[i, 'A(2)'] = cosine_rule(d12,d23,d13,0)
            final_df.at[i, 'A(3)'] = cosine_rule(d13,d23,d12,0)
            # print(final_df)
            
        # print(final_df)
        complete_final_df = pd.concat([complete_final_df,final_df])

        entire_atom_dis_dict.update(atom_dis_dict)
        #remove the following break if you want it to run on the entire input file
        # break

    return complete_final_df, list_of_atom_combi, entire_atom_dis_dict

def read_lexical_to_dict(input_file):
    dict_from_csv = {}

    with open(input_file, mode='r') as inp:
        reader = csv.reader(inp)
        dict_from_csv = {rows[0]:rows[1] for rows in reader}
    
    return dict_from_csv

def get_key_from_val(my_dict,val):
    for key, value in my_dict.items():
        if val == value:
            return key
 
    return "key doesn't exist"

def determine_label(row, atom, entire_atom_dis_dict):
    label1 = 0
    label2 = 0
    label3 = 0
    i = 0

    a1 = row['atom1']
    a2 = row['atom2']
    a3 = row['atom3']
    atom_in_a_triangle = (a1, a2, a3)

    s1 = atom.loc[atom['atom_number'] == a1, 'seq'].item()
    s2 = atom.loc[atom['atom_number'] == a2, 'seq'].item()
    s3 = atom.loc[atom['atom_number'] == a3, 'seq'].item()
    seq_in_a_triangle = (s1, s2, s3)

    atom_seq_dict = dict(zip( atom_in_a_triangle, seq_in_a_triangle))

    #initialized the sequence
    sorted_keys = sorted(atom_seq_dict, key=atom_seq_dict.get)
    sorted_atom_seq_dict = {}
    for w in sorted_keys:
        sorted_atom_seq_dict[w] = atom_seq_dict[w]

    #l1,l2,l3 are the seq value
    l1_value = list(sorted_atom_seq_dict.values())[0]
    l2_value = list(sorted_atom_seq_dict.values())[1]
    l3_value = list(sorted_atom_seq_dict.values())[2]
    
    #these are the key from the dictionary
    l1_key = list(sorted_atom_seq_dict.keys())[0]
    l2_key = list(sorted_atom_seq_dict.keys())[1]
    l3_key = list(sorted_atom_seq_dict.keys())[2]

    #the lengths of the triangle
    dist1_2 = entire_atom_dis_dict.get((l1_key,l2_key))
    if dist1_2 is None:
        dist1_2 = entire_atom_dis_dict.get((l2_key,l1_key))

    dist1_3 = entire_atom_dis_dict.get((l1_key,l3_key))
    if dist1_3 is None:
        dist1_3 = entire_atom_dis_dict.get((l3_key,l1_key))

    dist2_3 = entire_atom_dis_dict.get((l2_key,l3_key))
    if dist2_3 is None:
        dist2_3 = entire_atom_dis_dict.get((l3_key,l2_key))

    if(l1_value != l2_value) and (l2_value != l3_value) and (l1_value != l3_value): 
        #i did not use the max function, since I know which one in the keys are largest and smallest
        label1 = l3_key
        label3 = l1_key
        label2 = l2_key

    elif(l1_value == l2_value) and (l2_value == l3_value):
        #the rules here follow the code that Dr. Xu gave me
        if(dist1_2) >= max((dist1_2,dist1_3,dist2_3)):
            label1 = l1_key
            label2 = l2_key
            label3 = l3_key
        elif(dist1_3) >= max((dist1_2,dist1_3,dist2_3)):
            label1 = l1_key
            label2 = l3_key
            label3 = l2_key
        else:
            label1 = l2_key
            label2 = l3_key
            label3 = l1_key

    elif(l1_value == l2_value) and (l2_value != l3_value):
        label3 = l3_key
        
        if dist1_3 >= dist2_3:
            label1 = l1_key
            label2 = l2_key
        else:
            label1 = l2_key
            label2 = l1_key

    elif(l1_value != l2_value) and (l2_value == l3_value):
        label1 = l1_key
        if dist1_2 > dist2_3:
            label2 = l2_key
            label3 = l3_key
        else:
            label2 = l3_key
            label3 = l2_key

    maxDist = max(dist1_2,dist1_3,dist2_3)
    return label1,label2,label3, maxDist

def calculate_d3_theta(label1, label2, label3, entire_atom_dis_dict):
    dist1_2 = entire_atom_dis_dict.get((label1,label2))
    if dist1_2 is None:
        dist1_2 = entire_atom_dis_dict.get((label2,label1))

    dist1_3 = entire_atom_dis_dict.get((label1,label3))
    if dist1_3 is None:
        dist1_3 = entire_atom_dis_dict.get((label3,label1))

    dist2_3 = entire_atom_dis_dict.get((label2,label3))
    if dist2_3 is None:
        dist2_3 = entire_atom_dis_dict.get((label3,label2))

    d12Mid = float(dist1_2)/2

    median = math.sqrt((2*dist1_3**2 + 2*dist2_3**2 - dist1_2**2)/4)

    s1 = dist1_3
    s2 = dist1_2/2
    s3 = median
    if(abs(median)<= 1e-9):
        Theta1 = 180
    elif abs((s1**2 - s2**2 - s3**2)/(2*s2*s3)) > 1.0:
        Theta1 = 0
    else:
        Theta1 = 180*(math.acos((s1**2-s2**2-s3**2)/(2*s2*s3)))/3.14

    if Theta1 <= 90:
        Theta = Theta1
    else:
        Theta = abs(180-Theta1)
    

    return median, Theta

def getBinNumber(theta, decider):
    binForThetaT = [0,12.11,17.32,21.53,25.21,28.54,31.64,34.55,37.34,40.03,42.64,45.17,47.64,50.05,52.43,54.77,57.08,59.38,61.64,63.87,66.09,68.30,70.5,72.69,
        79.2,81.36,83.51,85.67,87.8,90]
    binForDt = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52]

    if decider == 0:
        for i in range(len(binForThetaT)):
            if binForThetaT[i] > theta:
                thetaT = i
                print(index)
                break
        return thetaT
    elif decider == 1:
        for i in range(len(binForThetaT)):
            if binForDt[i] > theta:
                dT = i
                print(index)
                break
        return dT

    

def calculate_key(theta, maxDist, li1, li2, li3):
    thetaT = getBinNumber(theta, 0)
    dT = getBinNumber(maxDist, 0)
    m=109

    firstTerm = thetaT*dT*(int(li1)-1)*m*m
    secondTerm = thetaT*dT*(int(li2)-1)*m
    thirdTerm = thetaT*dT*(int(li3)-1)
    fourthTerm = thetaT*(maxDist-1)
    fifthTerm = thetaT-1

    return firstTerm + secondTerm + thirdTerm + fourthTerm + fifthTerm, thetaT, dT
 

ppdb.read_pdb("5zk1.pdb")
atom = ppdb.df['ATOM']
atom = atom[['atom_number', 'atom_name', 'residue_number', 'x_coord', 'y_coord', 'z_coord']]
# atom.to_csv("refer.csv")

length_and_angle, list_of_atom_combi, entire_atom_dis_dict = length_angle_df(atom)

#convert drug_atom_lexical to dictionary
lexical = read_lexical_to_dict("drug_atom_lexical.csv")

#adding new column seq
atom['seq'] = atom['atom_name'].map(lexical)

# atom.to_csv("atom.csv")

length_and_angle = length_and_angle.reset_index()

length_and_angle["label1"] = ""
length_and_angle["label2"] = ""
length_and_angle["label3"] = ""
length_and_angle["median"] = ""
length_and_angle["theta"] = ""
length_and_angle["key"] = ""

output = open('output.txt', 'w')
for index, row in length_and_angle.iterrows():
    label1,label2,label3, maxDist = determine_label(row, atom, entire_atom_dis_dict)
    length_and_angle.at[index, 'label1'] = label1
    length_and_angle.at[index, 'label2'] = label2
    length_and_angle.at[index, 'label3'] = label3

    median, theta = calculate_d3_theta(label1, label2, label3, entire_atom_dis_dict)
    length_and_angle.at[index, 'median'] = median
    length_and_angle.at[index, 'theta'] = theta

    li1 = atom.loc[atom.atom_number == label1, 'seq'].item()
    li2 = atom.loc[atom.atom_number == label2, 'seq'].item()
    li3 = atom.loc[atom.atom_number == label3, 'seq'].item()

    #I did not follow Dr. Xu way to calculate the theta, becau
    key,thetaT,dT = calculate_key(theta, maxDist, li1, li2, li3)
    length_and_angle.at[index, 'key'] = key

    atom_name1 = atom.loc[atom.atom_number == label1, 'atom_name'].item() 
    atom_name2 = atom.loc[atom.atom_number == label2, 'atom_name'].item() 
    atom_name3 = atom.loc[atom.atom_number == label3, 'atom_name'].item() 

    x1 = atom.loc[atom.atom_number == label1, 'x_coord'].item()
    y1 = atom.loc[atom.atom_number == label1, 'y_coord'].item()
    z1 = atom.loc[atom.atom_number == label1, 'z_coord'].item()

    x2 = atom.loc[atom.atom_number == label2, 'x_coord'].item()
    y2 = atom.loc[atom.atom_number == label2, 'y_coord'].item()
    z2 = atom.loc[atom.atom_number == label2, 'z_coord'].item()

    x3 = atom.loc[atom.atom_number == label3, 'x_coord'].item()
    y3 = atom.loc[atom.atom_number == label3, 'y_coord'].item()
    z3 = atom.loc[atom.atom_number == label3, 'z_coord'].item()

    line = (str(key)+"\t"+str(atom_name1)+"\t"+str(li1)+"\t"+str(atom_name2)+"\t"+str(li2)+"\t"+str(atom_name3)
                        +"\t"+str(li3)+"\t"+str(thetaT)+"\t"+str(theta)+"\t"+str(dT)+"\t"+str(maxDist)
                            +"\t"+str(x1)+"\t"+str(y1)+"\t"+str(z1)+"\t"+str(x2)+"\t"+str(y2)+"\t"+str(z2)
                                +"\t"+str(x3)+"\t"+str(y3)+"\t"+str(z3)+"\n")
    output.write(line)

# length_and_angle.to_csv("angles_and_length_test_median_and_theta_k.csv")