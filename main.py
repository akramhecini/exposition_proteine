#import p5 as pp
import p7 as pp
import argparse 

parser = argparse.ArgumentParser()

parser.add_argument("-pdb", type = str, help = "this will generate a dataframe with the proteine coordinates", dest="pdb_file")

parser.add_argument("-s", type = float, help = "This is the VDW radius value of the solvant", dest="rayon_solvant")



dic_rayon = {"H": 1.2, "C":1.7 , "N": 1.55, "O":1.52, "F":1.47 , "P": 1.8, "S":1.8}

argmnt = parser.parse_args()

df = pp.coordonnees(argmnt.pdb_file)

expo = pp.exposition_calculate(df,dic_rayon, argmnt.rayon_solvant)

df["Exposition"] = expo

print(df.head())

