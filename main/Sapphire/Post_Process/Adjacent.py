import numpy as np
import scipy.sparse as spa
import os
import functools
import operator
from ase.data import covalent_radii, atomic_numbers
from Sapphire.Post_Process import DistFuncs

class Adjacency_Matrix():

    def __init__(self, System = None, Positions = None, Distances = None, 
                 Adj = None, agcn = None, Surf_Area = None, Surf_Atoms = None, CN = None,
                 Elements = None, R_Cut = None, Type = None, Frame = None, Metals = None):

        """ Robert
            Args:
                Not yet fully implemented but in theory will take arguments of a single frame
                of an xyz trajectory to then be iterated over and the cut-off imposed to be the
                nearest nieghbour distances
                
            Returns:
                Distances:
                    N x N array (N being number of atoms present) containing all pairwise distances 
                    between atoms
                    I.e., Entry[i][j] is the absolute distance between atoms i and j
                
                Adjacent: N x N array of 0s and 1s identifying the "Truth" of a given atom
                pairing being neighbours based on the criterion of R_Cut
                I.e., Dist<R_Cut returns a 1 for that entry as the condition has been met
                All diagonal elements are set to 0 as it is meaningless to be your own neighbour.
                
                
            Writable objects:
                Adjacency matrix - Sparse scipy matrix written in NPZ format
                
                CN - Coordination number of a given atom written as an N length vector
                
                AGCN -  atop generalised CN, similar to above
                
                Surface area/atoms - Uses above value to determine surface-type properties
        """ 
        self.System = System
        self.Positions = Positions
        self.Distances = Distances
        self.R_Cut = R_Cut
        self.Type = Type
        self.Frame = Frame
        self.Metals = Metals #Species present
        self.Elements = Elements #List of atomic elements in 1:1 correspondance with coordinates
        self.NumAdj = np.zeros((len(self.Positions), len(self.Positions)), dtype=np.float64) #Instantiate the 'dense' matrix
        
        #Consider the calculable objects below
        self.Adjacent = [] #Primary point of this object
        self.CN = CN #Coordination number
        self.agcn = agcn #Atop generalised coordination number
        self.Surf_Area = Surf_Area #Self explanatory
        self.Surf_Atoms = Surf_Atoms #Number of surface atoms
        self.calculate_adj() #Calculates adjacency matrix, nothing can happen without this
        
        #Call the write objects to handle everything else
        #Note that subsequent quantities are computed via cascading if statements
        #It's inellegent but keeps everything related under one roof
        self.write()

    def ensure_dir(self, base_dir='', file_path=''):
        """

        Robert:

            A simple script to verify the existence of a directory
            given the path to it. If it does not exist, will create it.

        """

        directory = base_dir + file_path
        if not os.path.exists(directory):

            os.makedirs(directory)

    def MakeFile(self, Attributes):
        self.out = self.System['base_dir'] + Attributes['Dir'] + Attributes['File']

        if not os.path.isfile(self.out):
            with open(self.System['base_dir'] + Attributes['Dir'] + Attributes['File'], 'w') as out:
                out.close()
        else:
            pass

    def calculate_adj(self):
        try:
            if self.Type == 'Homo':
                self.Distances = DistFuncs.Euc_Dist(self.Positions, homo=True,
                                                    specie=self.Metals, elements=self.Elements)
                self.Positions = DistFuncs.get_subspecieslist(self.Metals, self.Elements, self.Positions)
                
            elif self.Type == 'Hetero':
                
                self.Positions = DistFuncs.Hetero(self.Positions, self.Metals, self.Elements)
                self.Distances = functools.reduce(operator.iconcat, self.Positions, [])
                
            Tick = 0
            for i in range(1,len(self.Positions)):
                self.Adjacent.append(self.Distances[Tick:Tick+len(self.Positions)-i])
                Tick += (len(self.Positions)-i)
            for i in range(len(self.Positions)):
                for j in range(len(self.Positions)-(i+1)):
                    self.NumAdj[i][j+i+1] = self.Adjacent[i][j]
                    self.NumAdj[j+i+1][i] = self.Adjacent[i][j]
            
            self.Adjacent=(self.NumAdj<self.R_Cut).astype(int) #Evaluate if a pair are within R_Cut of eachother
            np.fill_diagonal(self.Adjacent,0)
                
            self.Adjacent = spa.csc_matrix(self.Adjacent)
        except Exception as e:
            with open(self.System['base_dir'] + 'Sapphire_Errors.log', 'a') as f:
                f.write("Exception raised whilst computing the Adjacency Matrix:\n%s"%e)

    def get_coordination(self):
        try:
            Temp = spa.csr_matrix.todense(self.Adjacent)
            self.NN = np.array([ Temp[i].sum() for i in range(len(Temp)) ])
        except Exception as e:
            with open(self.System['base_dir'] + 'Sapphire_Errors.log', 'a') as f:
                f.write("Exception raised whilst computing the Coordination:\n%s"%e)

    def get_coordination_hetero(self):
        try:
            self.BoolAdj = (self.Adjacent<self.R_Cut).astype(int)
            self.CoordA = [ self.Adjacent[i].sum() for i in range(len(self.Adjacent)) ]
            self.CoordB = [ self.Adjacent[:,j].sum() for j in range(len(self.Adjacent[0])) ]
        except Exception as e:
            print("Exception raised whilst computing the Hetero Coordination:\n%s"%e)

    def ReturnAdj(self):
        """
        This function is mostly used to store the sparse matrix short term into
        result cache of the process object. Things like the mix
        if NN is True:    
            return(np.array(agcn, dtype = float),Matrix)
        elif NN is False:ing parameter
        have this as a depedendency requiring differently computed objects.
        I.e., Homo_Adj[1,2] & HeAdj 
        """
        return self.Adjacent

    def agcn_generator(self):
        
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
        try:
            Matrix = self.Adjacent.sum(axis=1).getA1() #This is an ordered list of number of NN each atom has
            I_Row,_,_  = spa.find(self.Adjacent)       #Indices of rows and columns with none-zero adjacency
            
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
                self.AGCN = np.array(agcn, dtype=float)
        except Exception as e:
            print("Exception raised whilst computing the agcn:\n%s"%e)

    def Surface_Area(self):
        
        """
        Computes the approximate surface area of the cluster in accordance with 
        ACS Catal. 2020, 10, 6, 3911–3920
        
        A = sum_{atoms} (1/3) pi r_{atom}^{2} (12 - aGCN_{atom})
        
        This function will be passed directly into the Process Module and requires
        the aGCN value computed by the above function, the list of elements in the 
        system and the present atomic species.
        Will return a float value 
        """
        try:
            if self.Type == 'Homo':
                Radius =  covalent_radii[atomic_numbers[self.Metals]]
                Homo_aGCN = [ float(self.AGCN[i]) for i,x in enumerate(self.Elements) if x == self.Metals ]
                Temp = [ 12 - x for x in Homo_aGCN ]
                self.Area = (1/3) * np.pi * Radius**2 * sum(Temp)
                
            elif self.Type == 'Full':
                Radii = [ (x, covalent_radii[atomic_numbers[x]]) for x in self.Metals ]
                print(Radii)
                Mod_aGCN = [ 12 - float(x) for x in self.AGCN ]
                print("Size of modified agcn is %s.\n" %len(Mod_aGCN))
                T1 = []; T2 = []
                for i,x in enumerate(Mod_aGCN):
                    if Radii[0][0] == self.Metals[i]:
                        T1.append( (Radii[0][1]**2)*x )
                    elif Radii[1][0] == self.Metals[i]:
                        T2.append( (Radii[1][1]**2)*x )  
                self.Area = (1/3) * np.pi * (sum(T1) + sum(T2))
        except Exception as e:
            print("Exception raised whilst computing the Surfce Area:\n%s"%e)

    def Surface_Atoms(self):
        try:
            if self.Type == 'Homo':
                Homo_aGCN = np.array(
                    [ float(self.AGCN[i]) for i,x in enumerate(self.Elements) if x == self.Metals ], 
                                     dtype = float
                                     )
                    
                Mask = Homo_aGCN < 9.1
                self.Surf_At = sum(Mask)
                
            elif self.Type == 'Full':
                Temp = np.array(
                    [ float(x) for x in self.AGCN ], dtype = float)
                Mask = Temp < 9.1
                self.Surf_At = sum(Mask)
        except Exception as e:
            print("Exception raised whilst computing the Surfce Atoms:\n%s"%e)

    def write(self):
    
        if self.Type == 'Full':
            from Sapphire.Utilities import OutputInfoFull as Out  # Case 1
            
            #Write object for the CoM
            Attributes = getattr(Out, str('adj')) #Loads in the write information for the object 
            OutFile = self.System['base_dir'] + Attributes['Dir'] + Attributes['File']
            self.ensure_dir(base_dir=self.System['base_dir'], file_path=Attributes['Dir'])   
            self.MakeFile(Attributes)
            self.filename = self.System['base_dir'] + Attributes['Dir'] + 'File%s' % str(self.Frame)
            Mat = spa.csr_matrix.todense(self.Adjacent)
            with open(self.filename, 'w') as f:
                for line in Mat:
                    np.savetxt(f, line, fmt='%d')
                
            if self.CN:
                self.get_coordination()
                #Write object for the CoMDistances
                Attributes = getattr(Out, str('nn')) #Loads in the write information for the object                
                OutFile = self.System['base_dir'] + Attributes['Dir'] + Attributes['File']
                self.ensure_dir(base_dir=self.System['base_dir'], file_path=Attributes['Dir'])   
                self.MakeFile(Attributes)
                with open(OutFile, 'a') as outfile:
                    outfile.write(str(self.Frame) + ' ' +  ' '.join(str(item) for item in self.NN) +'\n')

            if self.agcn:
                self.agcn_generator()
                #Write object for the homo CoM distances
                Attributes = getattr(Out, str('agcn')) #Loads in the write information for the object                  
                OutFile = self.System['base_dir'] + Attributes['Dir'] + Attributes['File']
                self.ensure_dir(base_dir=self.System['base_dir'], file_path=Attributes['Dir'])   
                self.MakeFile(Attributes)
                with open(OutFile, 'a') as outfile:
                    outfile.write(str(self.Frame) + ' ' +  ' '.join(str(item) for item in self.AGCN) +'\n')
                    
            if self.Surf_Area:
                self.Surface_Area()
                #Write object for the homo CoM distances
                Attributes = getattr(Out, str('surf_area')) #Loads in the write information for the object                  
                OutFile = self.System['base_dir'] + Attributes['Dir'] + Attributes['File']
                self.ensure_dir(base_dir=self.System['base_dir'], file_path=Attributes['Dir'])   
                self.MakeFile(Attributes)
                with open(OutFile, 'a') as outfile:
                    outfile.write(str(self.Frame) + ' ' +  ' '.join(str(item) for item in self.Area) +'\n')
                    
            if self.Surf_Atoms:
                self.Surface_Atoms()
                #Write object for the homo CoM distances
                Attributes = getattr(Out, str('surf_atoms')) #Loads in the write information for the object                  
                OutFile = self.System['base_dir'] + Attributes['Dir'] + Attributes['File']
                self.ensure_dir(base_dir=self.System['base_dir'], file_path=Attributes['Dir'])   
                self.MakeFile(Attributes)
                with open(OutFile, 'a') as outfile:
                    outfile.write(str(self.Frame) + ' ' +  ' '.join(str(item) for item in self.Surf_At) +'\n')                    

        elif self.Type == 'Homo':
            from Sapphire.Utilities import OutputInfoHomo as Out  # Case 2
            
            #Write object for the homo CoM 
            Attributes = getattr(Out, str('hoadj')) #Loads in the write information for the object 
            OutFile = self.System['base_dir'] + Attributes['Dir'] + Attributes['File']+self.Metals
            self.ensure_dir(base_dir=self.System['base_dir'], file_path=Attributes['Dir'])   
            self.MakeFile(Attributes)
            self.filename = self.System['base_dir'] + Attributes['Dir'] + Attributes['File']+self.Metals +'File%s' % str(self.Frame)
            self.Mat = spa.csr_matrix.todense(self.Adjacent)
            with open(self.filename, 'w') as f:
                for line in self.Mat:
                    np.savetxt(f, line, fmt='%d')
                    
            if self.CN:
                self.get_coordination()
                #Write object for the homo CoM distances
                Attributes = getattr(Out, str('honn')) #Loads in the write information for the object                  
                OutFile = self.System['base_dir'] + Attributes['Dir'] + Attributes['File']+self.Metals
                self.ensure_dir(base_dir=self.System['base_dir'], file_path=Attributes['Dir'])   
                self.MakeFile(Attributes)
                with open(OutFile, 'a') as outfile:
                    outfile.write(str(self.Frame) + ' ' +  ' '.join(str(item) for item in self.NN) +'\n')
                    
            if self.Surf_Area:
                self.Surface_Area()
                #Write object for the homo CoM distances
                Attributes = getattr(Out, str('surf_area')) #Loads in the write information for the object                  
                OutFile = self.System['base_dir'] + Attributes['Dir'] + Attributes['File']
                self.ensure_dir(base_dir=self.System['base_dir'], file_path=Attributes['Dir']) + self.Metals
                self.MakeFile(Attributes)
                with open(OutFile, 'a') as outfile:
                    outfile.write(str(self.Frame) + ' ' +  ' '.join(str(item) for item in self.Area) +'\n')
                    
            if self.Surf_Atoms:
                self.Surface_Atoms()
                #Write object for the homo CoM distances
                Attributes = getattr(Out, str('surf_atoms')) #Loads in the write information for the object                  
                OutFile = self.System['base_dir'] + Attributes['Dir'] + Attributes['File']  + self.Metals
                self.ensure_dir(base_dir=self.System['base_dir'], file_path=Attributes['Dir'])   
                self.MakeFile(Attributes)
                with open(OutFile, 'a') as outfile:
                    outfile.write(str(self.Frame) + ' ' +  ' '.join(str(item) for item in self.Surf_At) +'\n')   
                    
        elif self.Type == 'Hetero':
            from Sapphire.Utilities import OutputInfoHetero as Out  # Case 3
            
            Attributes = getattr(Out, str('headj')) #Loads in the write information for the object 
            OutFile = self.System['base_dir'] + Attributes['Dir'] + Attributes['File'] + str(self.Frame)
            self.ensure_dir(base_dir=self.System['base_dir'], file_path=Attributes['Dir'])   
            self.MakeFile(Attributes)
            self.filename = self.System['base_dir'] + Attributes['Dir'] + Attributes['File']+self.Metals +'File%s' % str(self.Frame)
            self.Mat = spa.csr_matrix.todense(self.Adjacent)
            with open(self.filename, 'w') as f:
                for line in self.Mat:
                    np.savetxt(f, line, fmt='%d')
   
            if self.CN:
                self.get_coordination_hetero()
                Attributes = getattr(Out, str('henn')) #Loads in the write information for the object
                
                OutFile = self.System['base_dir'] + Attributes['Dir'] + Attributes['File'] + self.Metals[0]
                self.ensure_dir(base_dir=self.System['base_dir'], file_path=Attributes['Dir'])   
                self.MakeFile(Attributes)
                with open(OutFile, 'a') as outfile:
                    outfile.write(str(self.Frame) + ' ' +  ' '.join(str(item) for item in self.CoordA) +'\n')
                    
                OutFile = self.System['base_dir'] + Attributes['Dir'] + Attributes['File'] + self.Metals[1]
                self.ensure_dir(base_dir=self.System['base_dir'], file_path=Attributes['Dir'])   
                self.MakeFile(Attributes)
                with open(OutFile, 'a') as outfile:
                    outfile.write(str(self.Frame) + ' ' +  ' '.join(str(item) for item in self.CoordB) +'\n')