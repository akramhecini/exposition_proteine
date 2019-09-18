

import p7 as pp
import argparse 

parser = argparse.ArgumentParser()

parser.add_argument("-pdb", type = str, help = "this will generate a dataframe with the proteine coordinates", dest="pdb_file")

parser.add_argument("-s", type = float, help = "This is the VDW radius value of the solvant", dest="rayon_solvant")



dic_rayon = {"H": 1.2, "C":1.9 , "N": 1.5, "O":1.4, "F":1.47 , "P": 1.8, "S":1.85}

argmnt = parser.parse_args()

df = pp.coordonnees(argmnt.pdb_file)

data_new = pp.exposition_calculate(df,dic_rayon, argmnt.rayon_solvant)

data_new2 = pp.exposition_residues(df)

print(data_new2)

data_new.to_csv('expo_atome.csv', sep='\t')

data_new2.to_csv('expo_resid.csv', sep='\t')
