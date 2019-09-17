import pandas as pd
import numpy as np
import math


####################


# fonction pour générer des points sur une spehere : centre (0,0,0) et r = 1

dic_rayon = {"H": 1.2, "C":1.7 , "N": 1.55, "O":1.52, "F":1.47 , "P": 1.8, "S":1.8}

def sample_spherical(npoints, ndim=3):

    vec = np.random.randn(ndim, npoints)

    vec /= np.linalg.norm(vec, axis=0)

    return vec

phi = np.linspace(0, np.pi, 20)

theta = np.linspace(0, 2 * np.pi, 40)

x = np.outer(np.sin(theta), np.cos(phi))

y = np.outer(np.sin(theta), np.sin(phi))

z = np.outer(np.cos(theta), np.ones_like(phi))


def coordonnees(fichier):

    """This function will return a DataFrame from a PDB file. The DataFrame contains 6 different columns : 
    Atomes, Premiere lettre de chaque Atome,  coordonnees ( x , y , z), Rayon de VDW """

    dic =  {} 

    with open(fichier, "r") as filin:

	    for ligne in filin:

		    if ligne[0:5].strip() == 'ATOM':

			    crd = [float(ligne[30:38].strip()), float(ligne[38:46].strip()), float(ligne[46:54].strip())]

			    signature = ligne[6:11].strip()+'_'+ligne[12:16].strip()+"_"+ligne[17:20].strip()

			    dic[signature] = crd
			

    df=pd.DataFrame.from_dict(dic,orient='index')

    df.columns = ['x', 'y','z']

    df.reset_index(level=0, inplace=True)


# new data frame with split value columns 

    new = df["index"].str.split("_",expand = True) 

  
# making separate index column from new data frame 

    df["ind"] = new[0] 
  
# making separate Atome column from new data frame 

    df["Atome"] = new[1] 

# making separate Residue column from new data frame 

    df["Res"] = new[2] 
  
# Dropping old Name columns 

    df.drop(columns =["index"], inplace = True) 

    #df.drop(columns =["ind"], inplace = True)

    df['At'] = df['Atome'].astype(str).str[0]

#adding rayon de van der waals 

    col_rayon = [0]*len(df.index)

    return(df)



def exposition_calculate(df, dic_rayon, rayon_solvant):

    """
This Function calculates the exposition of an atom to the solvant. 

First of all, it parses the DataFtame previosuly generated, line by line, and for each atom, it calls the function 

sample_spherical that create a sphere of 100 points around each atom. The coordinates of the 100 points are then translocated 

using the atom's coordinates as the new sphere center and the VDW radius as the sphere's new radius. 

The surface of the sphere is then calculated. The function 'Distance' is then called. 

exposition_calculate Function finally return the expostion of all atoms to the solvant. 
    """

    col_rayon = [0]*len(df.index)
	
    exposition = []

    for index, row in df.iterrows():

        if row["At"] in dic_rayon:

            col_rayon[index] = dic_rayon[row["At"]]

            xi, yi, zi = sample_spherical(100)  #create sphere 

            xi = (xi+row['x'])*dic_rayon[row["At"]] #translocation 

            yi = (yi+row['y'])*dic_rayon[row["At"]]

            zi = (zi+row['z'])*dic_rayon[row["At"]]

            surface = math.pi*4*dic_rayon[row["At"]] #je calcule ma surface

            dis_euc = distance(xi,yi,zi,index,surface,df)

            expp = (sum(d > rayon_solvant for d in dis_euc)/len(dis_euc))*surface #calculate the exposition

            exposition.append(expp)

    return exposition



def distance(xi, yi, zi, index, surface,df):

    """
This Function is used to calculate the distance between each point of the sphere and all the atoms of the protein. 

The distances calculated are then put in a vector. The values of distance are then compared to 3.92 which is the

VDW radius of a molecule of water (solvant). The number of values of distance greater than 3.92 are computed and 

their percentage is calculated. This percentage is then multiplied by the surface of the sphere which calculates the 

value of exposition to the solvant 
  
    """
    df = df.drop([index]) #I delete the sphere's center from the dataframe

    dis_euc = [] #a list containing the distance values

    for index2, row in df.iterrows():#dataframe parsing

        p2 = list(df.loc[index2,["x","y","z"]]) #coordinates of an atom 

        for ind in range(len(xi)): # for each point of the 100 points 

            p1 = [xi[ind], yi[ind], zi[ind]] #coordinates of the 100 points 

            dist_p1_p2 = np.linalg.norm(np.array(p1)-np.array(p2)) #calculating the distance between p1 & p2

            dis_euc.append(dist_p1_p2)#put the distance in a list

    return (dis_euc)




