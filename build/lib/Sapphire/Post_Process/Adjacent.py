import numpy as np
import scipy.sparse as spa
import os
import functools
import operator
from ase.data import covalent_radii, atomic_numbers
from Sapphire.Post_Process import DistFuncs

class Adjacency_Matrix():
    """ 
    This class function is desined for the fast evaluation of adjacency matrices
    for nanoparticles. Data required will typically be fed in from the DistFuncs.EucDist
    function as distances are required to evaluate the truth in 2 atoms being 
    considered to be adjacent.
    
    The object Adjacency_Matrix.Adjacent IS the matrix stored in sparse diagonal form.
    
    It is not explicitly returned as there are qrite flags which may be fed into the __init__ 
    task to determine if writing is required. Two forms of writing are available:
        1. Save as numpy object which must be appropriately decompressed by python
        handlers. This is the most preferable form of long-term storage given the 
        scaling of an NxN matrix consisting only boolean value.
        2. Writeable as a full matrix in .txt form for each desired frame. 
        This is strongly counter-indicated due to the heavy nature of storing/writing
        this objec. Though this feature has been left for those who wish to seamlessly maniplate
        the adjacency matrix in ways otherwise not anticipated by the Sapphire development team.
    
        Args:
            System : Type - Dict
                Description - Base system information regarding directories.
                Not necessary for separate use outside of Sapphire core, 
                but writing output is not possible without reference directories.
                
            Positions : Type - numpy array
                Description - N X 3 array of atomic positions to be passed.
                Generally will be handled by the ase.Atoms.positions scheme,
                though this can be handled manually by an experienced user.
                
            Distances : Type numpy array
                Description - N(N-1)/2 length numpy array of distances computed
                from the Sapphire.Post_Process.DistFuncs module.
                
            Adj : Type - Boolean
                Description - Whether or not to compute and write this quantity.
                
            agcn : Type - Boolean
                Description - Whether or not to compute and write this quantity.
                
            Surf_Area : Type - Boolean
                Description - Whether or not to compute and write this quantity.
                
            Surf_Atoms : Type - Boolean
                Description - Whether or not to compute and write this quantity.
                 
            CN : Type - Boolean
                Description - Whether or not to compute and write this quantity.
                
            Elements : Type - numpy array
                Description -  N length vector of strings containing the names
                of the atomic species considered.
                
            R_Cut : Type - Float
                Description - Cut-off distance for atoms to be considered neighbours.
                This value can either be passed directly, or more commonly computed
                from the Sapphire.Post_Process.Kernels module.
                
            Type : Type - String
                Description - Whether or not to do the full set of calculations
                including the writing of outputs to specific files. Alternate 
                computations are 'Homo' or 'Hetero'.
                
            Frame : Type - Integer
                Description - The frame being considered. All this really does
                is label the frame in written output files. Can be ignored if doing
                a single frame calculation.
                
            Metals : Type - List
                Description - List of metal species considered. Essentially just the
                set of objects passed by Elements.
            
        Returns : None
            
            
        Writable objects:
            Adjacency matrix - Sparse scipy matrix written in NPZ format
            
            CN - Coordination number of a given atom written as an N length vector
            
            AGCN -  atop generalised CN, similar to above
            
            Surface area/atoms - Uses above value to determine surface-type properties
    """ 
    
    def __init__(self, System = None, Positions = None, Distances = None, 
                 Adj = None, agcn = None, Surf_Area = None, Surf_Atoms = None, CN = None,
                 Elements = None, R_Cut = None, Type = None, Frame = 0, Metals = None):

        self.System = System
        self.Positions = Positions
        self.Distances = Distances
        self.R_Cut = R_Cut
        self.Type = Type
        self.Frame = Frame
        self.Metals = Metals #Species present
        self.Elements = Elements #List of atomic elements in 1:1 correspondance with coordinates
        
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
                                                    specie=self.Metals[0], elements=self.Elements)
                self.Positions = DistFuncs.get_subspecieslist(self.Metals[0], self.Elements, self.Positions)
                
                self.NumAdj = np.zeros((len(self.Positions), len(self.Positions)), dtype=np.float64) #Instantiate the 'dense' matrix
                
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
            elif self.Type == 'Hetero':
                
                self.Positions = DistFuncs.Hetero(self.Positions, self.Metals, self.Elements)
                self.Distances = np.asarray(functools.reduce(operator.iconcat, self.Positions, []))
                self.Adjacent = self.Distances.reshape(np.shape(self.Positions))
                self.Adjacent=(self.Adjacent<self.R_Cut).astype(int)
                self.Adjacent = spa.csc_matrix(self.Adjacent)
            
            elif self.Type == 'Full':
                self.NumAdj = np.zeros((len(self.Positions), len(self.Positions)), dtype=np.float64)
                
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
                Mod_aGCN = [ 12 - float(x) for x in self.AGCN ]
                T1 = []; T2 = []
                for i,x in enumerate(Mod_aGCN):
                    if Radii[0][0] == self.Elements[i]:
                        T1.append( (Radii[0][1]**2)*x )
                    elif Radii[1][0] == self.Elements[i]:
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
                self.Surf_At = Mask.astype(int)
        except Exception as e:
            print("Exception raised whilst computing the Surfce Atoms:\n%s"%e)

    def write(self):
    
        if self.Type == 'Full':
            from Sapphire.IO import OutputInfoFull as Out  # Case 1
            
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
                    outfile.write(str(self.Frame) + ' ' +  str(self.Area) +'\n')
                    
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
            from Sapphire.IO import OutputInfoHomo as Out  # Case 2
            
            #Write object for the homo CoM 
            Attributes = getattr(Out, str('hoadj')) #Loads in the write information for the object 
            OutFile = self.System['base_dir'] + Attributes['Dir'] + Attributes['File']+self.Metals[0]
            self.ensure_dir(base_dir=self.System['base_dir'], file_path=Attributes['Dir'])   
            self.MakeFile(Attributes)
            self.filename = self.System['base_dir'] + Attributes['Dir'] + Attributes['File']+self.Metals[0] +'File%s' % str(self.Frame)
            self.Mat = spa.csr_matrix.todense(self.Adjacent)
            with open(self.filename, 'w') as f:
                for line in self.Mat:
                    np.savetxt(f, line, fmt='%d')
                    
            if self.CN:
                self.get_coordination()
                #Write object for the homo CoM distances
                Attributes = getattr(Out, str('honn')) #Loads in the write information for the object                  
                OutFile = self.System['base_dir'] + Attributes['Dir'] + Attributes['File']+self.Metals[0]
                self.ensure_dir(base_dir=self.System['base_dir'], file_path=Attributes['Dir'])   
                self.MakeFile(Attributes)
                with open(OutFile, 'a') as outfile:
                    outfile.write(str(self.Frame) + ' ' +  ' '.join(str(item) for item in self.NN) +'\n')
                    
            if self.Surf_Area:
                self.Surface_Area()
                #Write object for the homo CoM distances
                Attributes = getattr(Out, str('surf_area')) #Loads in the write information for the object                  
                OutFile = self.System['base_dir'] + Attributes['Dir'] + Attributes['File']
                self.ensure_dir(base_dir=self.System['base_dir'], file_path=Attributes['Dir']) + self.Metals[0]
                self.MakeFile(Attributes)
                with open(OutFile, 'a') as outfile:
                    outfile.write(str(self.Frame) + ' ' +  ' '.join(str(item) for item in self.Area) +'\n')
                    
            if self.Surf_Atoms:
                self.Surface_Atoms()
                #Write object for the homo CoM distances
                Attributes = getattr(Out, str('surf_atoms')) #Loads in the write information for the object                  
                OutFile = self.System['base_dir'] + Attributes['Dir'] + Attributes['File']  + self.Metals[0]
                self.ensure_dir(base_dir=self.System['base_dir'], file_path=Attributes['Dir'])   
                self.MakeFile(Attributes)
                with open(OutFile, 'a') as outfile:
                    outfile.write(str(self.Frame) + ' ' +  ' '.join(str(item) for item in self.Surf_At) +'\n')   
                    
        elif self.Type == 'Hetero':
            from Sapphire.IO import OutputInfoHetero as Out  # Case 3
            
            Attributes = getattr(Out, str('headj')) #Loads in the write information for the object 
            OutFile = self.System['base_dir'] + Attributes['Dir'] + Attributes['File'] + str(self.Frame)
            self.ensure_dir(base_dir=self.System['base_dir'], file_path=Attributes['Dir'])   
            self.MakeFile(Attributes)
            self.filename = self.System['base_dir'] + Attributes['Dir'] + Attributes['File'] +'File%s' % str(self.Frame)
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