import scipy.sparse as spa
import numpy as np
import os
from ase.data import covalent_radii, atomic_numbers

"""
Function that reads an .xyz file and returns a list or the 
coordination numbers (CNs) corresponding to each atoms, and their
atop generalized coordination numbers (aGCNs). 
​
The length of the two lists is equal to the number of atoms in the 
cluster. 
​
The CN is the number of nearest neighbours (neighbours within 
a distance of r_cut_cn, calculated through the PDF) of each atom.
​
The aGCN is the sum of all the CNs of the neighboring atoms, 
divided by the CN of the atoms in the bulk, which is equal to 12 
for the atop sites.

"""

def agcn_generator(adj=None, NN=False):
    
    """
    
    Robert:
        
        Arguments:
            
            adj - The sparse matrix from the adjacency module. It contains
            only binary truth elements regarding two neighbours being adjacent
            or not.
            
            
        Returns:
            
            agcn - List of agcn values for each atom in a single trajectory frame
            
            Matrix - np.array of the number of nearest nieghbours each atom has
            at the given snapshot.
            
        Note that a frame is not specified as it is understood to be called in conjunction
        with a function which reads frame by frame meaning that it is never ambiguous as to
        which frame is being evaluated.
        
    """
    
    if adj is None:
        raise TypeError('You have not specified your adjacecny matrix.')
        
    Matrix = adj.sum(axis=1).getA1() #This is an ordered list of number of NN each atom has
    I_Row,_,_  = spa.find(adj)       #Indices of rows and columns with none-zero adjacency
    
    agcn=[]
    Tick=0                          #Allows us to run along the length of the bonds found in I_Row
    
    #In principle, the following routine is equivalent to that written by Elena

    for i in range(len(Matrix)):
        Temp_List=[];cc=Matrix[i];
        for j in range(cc):
            Temp = I_Row[Tick:(Tick+Matrix[i])]
            Temp_List.append(Matrix[Temp[j]])
        agcn.append("%.3f" % (sum(Temp_List)/12.0))
        Tick+=Matrix[i]
    if NN is True:    
        return(np.array(agcn, dtype = float),Matrix)
    elif NN is False:
        return np.array(agcn, dtype=float)
    
def Surface_Area(aGCN, Elements, Species = None, Homo = False):
    
    """
    Computes the approximate surface area of the cluster in accordance with 
    ACS Catal. 2020, 10, 6, 3911–3920
    
    A = sum_{atoms} (1/3) pi r_{atom}^{2} (12 - aGCN_{atom})
    
    This function will be passed directly into the Process Module and requires
    the aGCN value computed by the above function, the list of elements in the 
    system and the present atomic species.
    Will return a float value 
    """
    
    if Homo:
        if not Species:
            return None
        else:
            Radius =  covalent_radii[atomic_numbers[Homo]]
            Homo_aGCN = [ float(aGCN[i]) for i,x in enumerate(Elements) if x == Homo ]
            Temp = [ 12 - x for x in Homo_aGCN ]
            return (1/3) * np.pi * Radius**2 * sum(Temp)
        
    else:
        Radii = [ (x, covalent_radii[atomic_numbers[x]]) for x in Species ]
        Mod_aGCN = [ 12 - float(x) for x in aGCN ]
        T1 = []; T2 = []
        for i,x in enumerate(Mod_aGCN):
            if Radii[0][0] == Elements[i]:
                T1.append( (Radii[0][1]**2)*x )
            elif Radii[1][0] == Elements[i]:
                T1.append( (Radii[1][1]**2)*x )  
        return (1/3) * np.pi * (sum(T1) + sum(T2))

            
def Surface_Atoms(aGCN, Elements = None, Species = None, Homo = False):
    if Homo:
        Homo_aGCN = np.array(
            [ float(aGCN[i]) for i,x in enumerate(Elements) if x == Homo ], 
                             dtype = float
                             )
            
        Mask = Homo_aGCN < 9.1
        return sum(Mask)
    else:
        Temp = np.array(
            [ float(x) for x in aGCN ], dtype = float)
        Mask = Temp < 9.1
        return sum(Mask)
    
"""

Robert:
    Below is the code written by Elena which has asap3 as a dependency.
    I cannot get this module to run on Gravity and so this module is 
    unable to produce results at the present time. It still works, but 
    one should keep in mind its limited functionality on archaeic hardware.

def cn_generator (positions, r_cut_cn):
    #Creating two empty lists
    cn=[]
    agcn=[]
    for i, atoms in enumerate(positions):
        atoms.set_cell([[100, 0, 0], [0, 100, 0], [0, 0, 100]])
        #Creating an empty list for each atom, to which the 
        #indices of all its nearest neighbours will be appended
        ind=[]
        for j in np.arange(len(atoms)):
            nl = FullNeighborList(r_cut_cn, atoms=atoms)
            indices, positions, distances = nl.get_neighbors(j)
            ind.append([int(k) for k in indices])
            distancej=[]    
            distancej.extend(distances**0.5)
            #The CN of each atom corresponds to the length of the 
            distance array calculated with asap3
            cnj=len(distancej)
            cn.append(cnj)
        for l in np.arange(len(atoms)):
            #List of all the indices of the neighbours of each atom
            cc=ind[l][:]
            list=[]
            #Appending the CN of each of these neighbors to a list and 
            #calculating their sum divided by 12, giving the atop GCN
            for m in range(len(cc)):
                list.append(cn[ind[l][m]])
                sm=sum(list)/12
            agcn.append(sm)
    return (cn, agcn)
    
"""