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
                # print('===========')
                # print(i)
                #bundle the coord of the 1st atom
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
                # atom_number_list = (i+1,j+1)
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
        # print(atom)
        atom_number = atom_df['atom_number'].tolist()
        atom_dis_dict = calculate_dist(atom_df)
        # print("---")
        # print(len(atom_dis_dict))

        # counter = 0
        tmp = combinations(atom_number,3)
        list_of_atom_combi = list(tmp)
        # for i in list_of_atom_combi:
        #     # counter = counter+1
        #     print(i)

        final_df = pd.DataFrame(list_of_atom_combi, columns=["atom1", "atom2", "atom3"])
        final_df['atom_name'] = atom
        final_df['D(1,2)'] = ""
        final_df['D(2,3)'] = ""
        final_df['D(1,3)'] = ""
        final_df['A(1)'] = ""
        final_df['A(2)'] = ""
        final_df['A(3)'] = ""
        #this process might be time consuming
        # final_df["D(1,2)"] = list(zip(final_df["atom1"],final_df["atom2"]))
        # final_df["D(2,3)"] = list(zip(final_df["atom2"],final_df["atom3"]))
        # final_df["D(1,3)"] = list(zip(final_df["atom1"],final_df["atom3"]))
        # print(final_df)
        
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
        break

    return complete_final_df, list_of_atom_combi, entire_atom_dis_dict

def read_lexical_to_dict(input_file):
    dict_from_csv = {}

    with open(input_file, mode='r') as inp:
        reader = csv.reader(inp)
        dict_from_csv = {rows[0]:rows[1] for rows in reader}
    
    # print(dict_from_csv)
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
    # print(atom_in_a_triangle)

    s1 = atom.loc[atom['atom_number'] == a1, 'seq'].item()
    s2 = atom.loc[atom['atom_number'] == a2, 'seq'].item()
    s3 = atom.loc[atom['atom_number'] == a3, 'seq'].item()
    seq_in_a_triangle = (s1, s2, s3)
    # print(seq_in_a_triangle)

    atom_seq_dict = dict(zip( atom_in_a_triangle, seq_in_a_triangle))
    # print(atom_seq_dict)

    #initialized the sequence
    sorted_keys = sorted(atom_seq_dict, key=atom_seq_dict.get)
    sorted_atom_seq_dict = {}
    for w in sorted_keys:
        sorted_atom_seq_dict[w] = atom_seq_dict[w]
    # print(sorted_atom_seq_dict)
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
        # print("case1")
        #redundant, just for reference
        # label1 = get_key_from_val(sorted_atom_seq_dict,l3_value)
        # label3 = get_key_from_val(sorted_atom_seq_dict,l1_value)
        # label2 = get_key_from_val(sorted_atom_seq_dict,l2_value)

        #i did not use the max function, since I know which one in the keys are largest and smallest
        label1 = l3_key
        label3 = l1_key
        label2 = l2_key
        # print(label1, label2, label3)

    elif(l1_value == l2_value) and (l2_value == l3_value):
        # print("case4")
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
            
        # print('-=-=-=-=')
        # print(label1, label2, label3)

    elif(l1_value == l2_value) and (l2_value != l3_value):
        # print("case3")
        label3 = l3_key
        
        if dist1_3 >= dist2_3:
            label1 = l1_key
            label2 = l2_key
        else:
            label1 = l2_key
            label2 = l1_key
        # print('-=-=-=-=')
        # print(label1, label2, label3)

    elif(l1_value != l2_value) and (l2_value == l3_value):
        # print("case2")
        label1 = l1_key

        if dist1_2 > dist2_3:
            label2 = l2_key
            label3 = l3_key
        else:
            label2 = l3_key
            label3 = l2_key
        # print('-=-=-=-=')
        # print(label1, label2, label3)

    return label1,label2,label3

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
    # print(entire_atom_dis_dict)
    # print(label1, label2, label3)

    # print(dist1_3)
    # print(dist1_2)
    # print(dist2_3)

    d12Mid = float(dist1_2)/2
    a1 = cosine_rule(dist1_3, dist1_2, dist2_3, 0)

    midD3 = cosine_rule(dist1_3, d12Mid, a1, 1)
    theta = cosine_rule(d12Mid, midD3, dist1_3, 0)

    return midD3, theta


ppdb.read_pdb("5zk1.pdb")
atom = ppdb.df['ATOM']
atom = atom[['atom_number', 'atom_name', 'residue_number', 'x_coord', 'y_coord', 'z_coord']]
# atom.to_csv("refer.csv")

length_and_angle, list_of_atom_combi, entire_atom_dis_dict = length_angle_df(atom)

#convert drug_atom_lexical to dictionary
lexical = read_lexical_to_dict("drug_atom_lexical.csv")
# print(lexical)

#adding new column seq
atom['seq'] = atom['atom_name'].map(lexical)

print(atom)
# atom.to_csv("atom.csv")

length_and_angle = length_and_angle.reset_index()

length_and_angle["label1"] = ""
length_and_angle["label2"] = ""
length_and_angle["label3"] = ""
length_and_angle["midD3"] = ""
length_and_angle["theta"] = ""

for index, row in length_and_angle.iterrows():
    # print(entire_atom_dis_dict)
    label1,label2,label3 = determine_label(row, atom, entire_atom_dis_dict)
    length_and_angle.at[index, 'label1'] = label1
    length_and_angle.at[index, 'label2'] = label2
    length_and_angle.at[index, 'label3'] = label3

    midD3, theta = calculate_d3_theta(label1, label2, label3, entire_atom_dis_dict)
    length_and_angle.at[index, 'midD3'] = midD3
    length_and_angle.at[index, 'theta'] = theta
    # break

print(length_and_angle)
length_and_angle.to_csv("angles_and_length.csv")