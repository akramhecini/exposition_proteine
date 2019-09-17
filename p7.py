import pandas as pd
import numpy as np
import math


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

    dic =  {} #a dictionary that will contain atoms and their coordinates

    with open(fichier, "r") as filin: #reading the pdb file 

	    for ligne in filin:

		    if ligne[0:5].strip() == 'ATOM':

			    crd = [float(ligne[30:38].strip()), float(ligne[38:46].strip()), float(ligne[46:54].strip())] # extracting the atoms coordinates

			    signature = ligne[6:11].strip()+'_'+ligne[12:16].strip()+"_"+ligne[17:20].strip()+"_"+ligne[22:26].strip() #a signature that is specific to each atom 

			    dic[signature] = crd # the signature is the key and the coordinates are the values 
			

    df=pd.DataFrame.from_dict(dic,orient='index') # transoforming the dictionary into  a dataframe

    df.columns = ['x', 'y','z'] #changing the columns names 

    df.reset_index(level=0, inplace=True) #reset the indexing 


# new data frame with split value columns 

    new = df["index"].str.split("_",expand = True) 

  
# making separate index column from new data frame 

    df["ind"] = new[0] 
  
# making separate Atome column from new data frame 

    df["Atome"] = new[1] 

# making separate Residue column from new data frame 

    df["Res"] = new[2] 
    df["nbr"] = new[3]
  
# Dropping old Name columns 

    df.drop(columns =["index"], inplace = True) 

    df.drop(columns =["ind"], inplace = True)

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

exposition_calculate Function finally return the expostion of all atoms to the solvant and the percentage of  exposition
    """

    col_rayon = [0]*len(df.index) # a new list that will be transformed later into a dataframe column, it will contain the VDW radius 
	
    exposition = [] # a list of the exposition values 
    prc = []

    for index, row in df.iterrows(): #parsing the dataframe

        if row["At"] in dic_rayon: # if the atom exists in the dictionary 

            col_rayon[index] = dic_rayon[row["At"]] # I extract the VDW radius value 

            xi, yi, zi = sample_spherical(100)  #create sphere 

            rayon_atome = dic_rayon[row["At"]]

            xi = (xi+row['x'])*rayon_atome #translocation  of x values 

            yi = (yi+row['y'])*rayon_atome #translocation  of y values 

            zi = (zi+row['z'])*rayon_atome #translocation  of z values 

            surface = math.pi*4*rayon_atome #je calcule ma surface

            dis_euc = distance(xi,yi,zi,index,surface,df) # calling the distance function 

            pourcentage = (sum(d > (rayon_solvant*2+rayon_atome) for d in dis_euc)/len(dis_euc)) #calculate the percentage of exposition to the solvant

            expp = pourcentage*surface #calculate the exposition

            prc.append(pourcentage)

            exposition.append(expp) #adding the exposition value to the list 

    df["pourcentage"] = prc

    df["Exposition"] = exposition


    return df #returning the percentages and the exposition values 



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


def exposition_residues(df):

    df_new = df.groupby(['Res','nbr'])["Exposition"].sum()

    return(df_new)
        
    




