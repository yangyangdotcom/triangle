import numpy as np
import pandas as pd
from biopandas.pdb import PandasPdb
from itertools import combinations
from scipy.spatial import distance
import csv
ppdb = PandasPdb().fetch_pdb('3eiy')

def readSample(filename):
    df = pd.read_csv(filename)
    # print(df)
    protein_file = df["protein"]
    chain = df["chain"]
    lig_id = df["ligand_id"]
    return protein_file, chain, lig_id

# protein_file,chain,lig_id = readSample("sample_details.csv")
# protein_file, chain, lig_id = readSample("sample_details.csv")

def readPDBfile(filename):
    ppdb.read_pdb(filename)
    #split the file into ATOM and HETATM
    
    #ATOM part
    recordName = ppdb.df['ATOM']['record_name']
    atomNo = ppdb.df['ATOM']['atom_number']
    residueName = ppdb.df['ATOM']['residue_name']
    chain_id = ppdb.df['ATOM']['chain_id']
    residueNumber = ppdb.df['ATOM']['residue_number']
    xCor = ppdb.df['ATOM']['x_coord']
    yCor = ppdb.df['ATOM']['y_coord']
    zCor = ppdb.df['ATOM']['z_coord']

    d = {'record_name':recordName, 'atom_number':atomNo, 'residue_name':residueName,'chain_id':chain_id, 'residue_number':residueNumber
        ,'x_coord':xCor, 'y_coord':yCor, 'z_coord':zCor}
    atom = pd.DataFrame(data=d)

    #HETATM part
    recordName = ppdb.df['HETATM']['record_name']
    atomNo = ppdb.df['HETATM']['atom_number']
    residueName = ppdb.df['HETATM']['residue_name']
    chain = ppdb.df['HETATM']['chain_id']
    xCor = ppdb.df['HETATM']['x_coord']
    yCor = ppdb.df['HETATM']['y_coord']
    zCor = ppdb.df['HETATM']['z_coord']
    
    d = {'record_name':recordName, 'atom_number':atomNo, 'residue_name':residueName, 'chain_id':chain,'x_coord':xCor, 'y_coord':yCor, 'z_coord':zCor}
    hetatm = pd.DataFrame(data=d)
    
    return atom, hetatm

ppdb.read_pdb("5zk1.pdb")
atom = ppdb.df['ATOM']
atom = atom[['atom_number', 'x_coord', 'y_coord', 'z_coord']]
atom_number = atom['atom_number'].tolist()
# atom.to_csv('test.csv')

# atom_dis_dict = {}
# for i in range(len(atom.index)):
#     for j in range(i+1,len(atom.index)):
#         if(i != j):
#             #bundle the coord of the atom
#             x1 = atom['x_coord'][i]
#             y1 = atom['y_coord'][i]
#             z1 = atom['z_coord'][i]
#             a = (x1,y1,z1)

#             #bundle the coord of the hetatm
#             x2 = atom['x_coord'][j]
#             y2 = atom['y_coord'][j]
#             z2 = atom['z_coord'][j]
#             b = (x2,y2,z2)

#             atom_number_list = (i,j)
#             atom_number_tuple = tuple(atom_number_list)
#             #calculate the euclidean distance
#             dist = str(distance.euclidean(a,b))
#             atom_dis_dict[atom_number_tuple] = dist

# print(atom_dis_dict)

counter = 0

tmp = combinations(atom_number,3)
for i in list(tmp):
    counter = counter+1
    print(i)

print(counter)

