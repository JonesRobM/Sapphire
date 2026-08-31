import numpy as np
from scipy.spatial.distance import cdist, pdist
from Sapphire.Utilities import errors
import os
from Sapphire.Utilities.log import get_logger

log = get_logger('Sapphire.Post_Process.DistFuncs')


def distance(a, b):
    
    dx = abs(a[0] - b[0])
     
    dy = abs(a[1] - b[1])
     
    dz = abs(a[2] - b[2])
 
    return np.sqrt(dx**2 + dy**2 + dz**2)

def CoMDist(positions, CoM = None, homo = False, specie = None, elements = None):
    
    if homo == False:
        return [distance(x, CoM) for x in positions]
    elif homo:
        Temp = get_subspecieslist(specie, elements, positions)
        CoM = get_CoM(Temp)
        return [distance(x, CoM) for x in Temp]
        
def get_CoM(positions):
    """Geometric centre (unweighted) of an (N, 3) array."""
    return np.average(np.asarray(positions, dtype=float), axis=0)

def get_subspecieslist(specie, elements, positions):
    Temp = np.column_stack((elements,positions))
    Temp = [x for x in Temp if x[0] == specie]
    return np.array(np.delete(Temp,0,1), dtype = np.float64)

def Euc_Dist(positions, homo = False, specie = None, elements = None):
    """All pair distances |r_i - r_j| for i < j (row-major order), as a 1-D array.

    With ``homo=True`` only atoms of ``specie`` are considered; returns None if fewer than two.
    Equivalent to the original double loop but via ``scipy.spatial.distance.pdist``.
    """
    if homo:
        Temp = get_subspecieslist(specie, elements, positions)
        if len(Temp) < 2:
            return None
        return pdist(np.asarray(Temp, dtype=float))
    return pdist(np.asarray(positions, dtype=float))

def Hetero(positions, species, elements):
    """Distances between every atom of species[0] and every atom of species[1]: (N_A, N_B) array."""
    TempA = get_subspecieslist(species[0], elements, positions)
    TempB = get_subspecieslist(species[1], elements, positions)
    return cdist(np.atleast_2d(TempA), np.atleast_2d(TempB))


class CoM_Dist():
    
    def __init__(self, System, Positions, CoM = None, Type = False, Specie = None, Elements = None, Frame = None):
        self.System = System
        self.Positions = Positions
        # Reference point: the caller's CoM if given, else the geometric centre.
        # (Previously this was set to Positions itself, so every 'CoMDist' was a 3-vector, not a distance.)
        self.CoM = get_CoM(Positions) if CoM is None else np.asarray(CoM, dtype=float)
        self.Type = Type
        self.Specie= Specie
        self.Elements = Elements
        self.Frame = Frame
        self.calculate()
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
        
    def calculate(self):
        
        if self.Type == 'Full':
            self.Dist = np.array([distance(x, self.CoM) for x in self.Positions])
        elif self.Type == 'Homo':
            Temp = get_subspecieslist(self.Specie, self.Elements, self.Positions)
            self.CoM = get_CoM(Temp)  # sub-species centre first, then distances to it
            self.Dist = np.array([distance(x, self.CoM) for x in Temp])
            self.MidDist = np.array([distance(x, self.CoM) for x in Temp])
        
    def write(self):
    
        if self.Type == 'Full':
            from Sapphire.IO import OutputInfoFull as Out  # Case 1
            
            #Write object for the CoM
            Attributes = getattr(Out, str('com')) #Loads in the write information for the object 
            OutFile = self.System['base_dir'] + Attributes['Dir'] + Attributes['File']
            self.ensure_dir(base_dir=self.System['base_dir'], file_path=Attributes['Dir'])   
            self.MakeFile(Attributes)
            with open(OutFile, 'a') as outfile:
                outfile.write(str(self.Frame) + ' ' +  ' '.join(str(item) for item in self.CoM) +'\n')

            #Write object for the CoMDistances
            Attributes = getattr(Out, str('comdist')) #Loads in the write information for the object                
            OutFile = self.System['base_dir'] + Attributes['Dir'] + Attributes['File']
            self.ensure_dir(base_dir=self.System['base_dir'], file_path=Attributes['Dir'])   
            self.MakeFile(Attributes)
            with open(OutFile, 'a') as outfile:
                outfile.write(str(self.Frame) + ' ' +  ' '.join(str(item) for item in self.Dist) +'\n')

        elif self.Type == 'Homo':
            from Sapphire.IO import OutputInfoHomo as Out  # Case 2
            
            #Write object for the homo CoM 
            Attributes = getattr(Out, str('hocom')) #Loads in the write information for the object 
            OutFile = self.System['base_dir'] + Attributes['Dir'] + Attributes['File']+self.Specie
            self.ensure_dir(base_dir=self.System['base_dir'], file_path=Attributes['Dir'])   
            self.MakeFile(Attributes)
            with open(OutFile, 'a') as outfile:
                outfile.write(str(self.Frame) + ' ' +  ' '.join(str(item) for item in self.CoM) +'\n')
   
            #Write object for the homo CoM distances
            Attributes = getattr(Out, str('hocomdist')) #Loads in the write information for the object                  
            OutFile = self.System['base_dir'] + Attributes['Dir'] + Attributes['File']+self.Specie
            self.ensure_dir(base_dir=self.System['base_dir'], file_path=Attributes['Dir'])   
            self.MakeFile(Attributes)
            with open(OutFile, 'a') as outfile:
                outfile.write(str(self.Frame) + ' ' +  ' '.join(str(item) for item in self.Dist) +'\n')
    
            #Write object for the sub-cluster specific homo CoM distances
            Attributes = getattr(Out, str('homidcomdist')) #Loads in the write information for the object                  
            OutFile = self.System['base_dir'] + Attributes['Dir'] + Attributes['File']+self.Specie
            self.ensure_dir(base_dir=self.System['base_dir'], file_path=Attributes['Dir'])   
            self.MakeFile(Attributes)
            with open(OutFile, 'a') as outfile:
                outfile.write(str(self.Frame) + ' ' +  ' '.join(str(item) for item in self.MidDist) +'\n')

    
class RDF():
    

    def __init__(self, System, Positions, Res=100, R_Cut=10.0, Type = None, Species = None, Elements = None, Frame = None):
            
            
        """ Robert
        
        Args:
            Res: 
                int data type representing how finely you wish to make 
                the grid. Usually set in the order of 100
            
            positions: 
                Single frame of xyz coordinates for a set of atoms
                Is expected to be iterated over and so will only take a single frame of xyz
            
            R_Cut: 
                Float type variable which indicates how far you wish to create
                the distribution for.
                Good practice is to set it to ~0.5 Diameter of the cluster
                Tested with 10 Angstroms
        Returns:
            Radii:
                A numpy array of all the radii the distribution has been computed over
                Will have length of "Resolution" and is to be used as the x axis on
                an RDF plot.
            
            G:
                A numpy array of the (unnormalised) calculated RDF values corresponding 
                to the respective radius in Radii. To be set on the y axis in a given
                RDF plot.
                
        """
        
        self.R_Cut = R_Cut
        self.System = System
        self.Res = Res
        self.Positions = Positions
        self.Type = Type
        self.Species = Species
        self.Elements = Elements
        self.Frame = Frame
        self.dr = self.R_Cut / self.Res #Increments to grow the spheres by
        self.Radii = np.linspace(0, self.R_Cut, self.Res) #List of Sphere radii to use
        self.Volumes=np.zeros(self.Res)
        self.G=np.zeros(self.Res)
        self.calculate()
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
        
    def calculate(self):
        """Shell-count RDF: G[j] = 2 x (pairs with j*dr <= r < (j+1)*dr) / (N_ref x shell volume)."""
        if self.Type == 'Hetero':
            A = get_subspecieslist(self.Species[0], self.Elements, self.Positions)
            B = get_subspecieslist(self.Species[1], self.Elements, self.Positions)
            r = cdist(np.atleast_2d(A), np.atleast_2d(B)).ravel()
            n_ref = len(np.atleast_2d(A))
        else:
            P = np.asarray(self.Positions, dtype=float)
            r = pdist(P)
            n_ref = len(P)
        idx = (r / self.dr).astype(int)
        idx = idx[(idx > 0) & (idx < self.Res)]
        self.G = 2.0 * np.bincount(idx, minlength=self.Res).astype(float)
        j = np.arange(self.Res)
        self.Volumes = n_ref * (4.0 / 3.0) * np.pi * (((j + 1) * self.dr) ** 3 - (j * self.dr) ** 3)
        self.G = self.G / self.Volumes

    def write(self):
        
        if self.Type == 'Full':
            from Sapphire.IO import OutputInfoFull as Out  # Case 1
          
            Attributes = getattr(Out, str('rdf')) #Loads in the write information for the object 
            OutFile = self.System['base_dir'] + Attributes['Dir'] + Attributes['File']
            self.ensure_dir(base_dir=self.System['base_dir'], file_path=Attributes['Dir'])   
            self.MakeFile(Attributes)
            with open(OutFile, 'a') as outfile:
                outfile.write(str(self.Frame) + ' ' +  ' '.join(str(item) for item in self.G) +'\n')

        
            Attributes = getattr(Out, str('rdfspace')) #Loads in the write information for the object                
            OutFile = self.System['base_dir'] + Attributes['Dir'] + Attributes['File']
            self.ensure_dir(base_dir=self.System['base_dir'], file_path=Attributes['Dir'])   
            self.MakeFile(Attributes)
            with open(OutFile, 'a') as outfile:
                outfile.write(str(self.Frame) + ' ' +  ' '.join(str(item) for item in self.Radii) +'\n')

        elif self.Type == 'Homo':
            from Sapphire.IO import OutputInfoHomo as Out  # Case 2
          
            Attributes = getattr(Out, str('hordf')) #Loads in the write information for the object 
            OutFile = self.System['base_dir'] + Attributes['Dir'] + Attributes['File']+self.Species
            self.ensure_dir(base_dir=self.System['base_dir'], file_path=Attributes['Dir'])   
            self.MakeFile(Attributes)
            with open(OutFile, 'a') as outfile:
                outfile.write(str(self.Frame) + ' ' +  ' '.join(str(item) for item in self.G) +'\n')
        
            Attributes = getattr(Out, str('hordfspace')) #Loads in the write information for the object         
            
            OutFile = self.System['base_dir'] + Attributes['Dir'] + Attributes['File']+self.Species
            self.ensure_dir(base_dir=self.System['base_dir'], file_path=Attributes['Dir'])   
            self.MakeFile(Attributes)
            with open(OutFile, 'a') as outfile:
                outfile.write(str(self.Frame) + ' ' +  ' '.join(str(item) for item in self.Radii) +'\n')

        elif self.Type == 'Hetero':
            from Sapphire.IO import OutputInfoHetero as Out  # Case 3
          
            Attributes = getattr(Out, str('herdf')) #Loads in the write information for the object 
            OutFile = self.System['base_dir'] + Attributes['Dir'] + Attributes['File']
            self.ensure_dir(base_dir=self.System['base_dir'], file_path=Attributes['Dir'])   
            self.MakeFile(Attributes)
            with open(OutFile, 'a') as outfile:
                outfile.write(str(self.Frame) + ' ' +  ' '.join(str(item) for item in self.G) +'\n')
      
        
            Attributes = getattr(Out, str('herdfspace')) #Loads in the write information for the object                   
            OutFile = self.System['base_dir'] + Attributes['Dir'] + Attributes['File']
            self.ensure_dir(base_dir=self.System['base_dir'], file_path=Attributes['Dir'])   
            self.MakeFile(Attributes)
            with open(OutFile, 'a') as outfile:
                outfile.write(str(self.Frame) + ' ' +  ' '.join(str(item) for item in self.Radii) +'\n')

class Pair_Dist():
    
    def __init__(self, System, Positions, Type = None, Specie = None, Elements = None, Frame = None):
        self.System = System
        self.Positions = Positions
        self.Type = Type
        self.Specie = Specie
        self.Elements = Elements
        self.Frame = Frame
        self.calculate()
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
        
    def calculate(self):
        
        if self.Type == 'Homo':
            try:
                self.distances = Euc_Dist(self.Positions, True, self.Specie, self.Elements)
                #(positions, homo = False, specie = None, elements = None)
            except Exception as e:
                errors.report(e, 'Pair-distance histogram failed: %s', base_dir=self.System['base_dir'] if self.System else None)
        elif self.Type == 'Hetero':
            try:
                self.distances = Hetero(self.Positions, self.Specie, self.Elements)
            except Exception as e:
                errors.report(e, 'Pair-distance histogram failed: %s', base_dir=self.System['base_dir'] if self.System else None)
        else:
            self.distances = Euc_Dist(self.Positions)
        
        self.bins = int(round(200/(1+20*np.exp(-len(self.distances)/1000)))) #Wait, what the fuck???
        self.a, b = np.histogram(self.distances, self.bins)
        bin_width = b[1]-b[0]
        self.bin_cents = [ b[i]+ bin_width for i in range(len(b)-1) ]
        #bin_cents, a
    
    def write(self):
        
        if self.Type == 'Full':
            from Sapphire.IO import OutputInfoFull as Out  # Case 1
          
            Attributes = getattr(Out, str('pair_distance')) #Loads in the write information for the object 
            OutFile = self.System['base_dir'] + Attributes['Dir'] + Attributes['File']
            self.ensure_dir(base_dir=self.System['base_dir'], file_path=Attributes['Dir'])   
            self.MakeFile(Attributes)
            with open(OutFile, 'a') as outfile:
                outfile.write(str(self.Frame) + ' ' +  ' '.join(str(item) for item in self.a) +'\n')

        
            Attributes = getattr(Out, str('pair_distancespace')) #Loads in the write information for the object                
            OutFile = self.System['base_dir'] + Attributes['Dir'] + Attributes['File']
            self.ensure_dir(base_dir=self.System['base_dir'], file_path=Attributes['Dir'])   
            self.MakeFile(Attributes)
            with open(OutFile, 'a') as outfile:
                outfile.write(str(self.Frame) + ' ' +  ' '.join(str(item) for item in self.bin_cents) +'\n')

        elif self.Type == 'Homo':
            from Sapphire.IO import OutputInfoHomo as Out  # Case 2
          
            Attributes = getattr(Out, str('hopair_distance')) #Loads in the write information for the object 
            OutFile = self.System['base_dir'] + Attributes['Dir'] + Attributes['File']+self.Specie
            self.ensure_dir(base_dir=self.System['base_dir'], file_path=Attributes['Dir'])   
            self.MakeFile(Attributes)
            with open(OutFile, 'a') as outfile:
                outfile.write(str(self.Frame) + ' ' +  ' '.join(str(item) for item in self.a) +'\n')
        
            Attributes = getattr(Out, str('hopair_distancespace')) #Loads in the write information for the object         
            
            OutFile = self.System['base_dir'] + Attributes['Dir'] + Attributes['File']+self.Specie
            self.ensure_dir(base_dir=self.System['base_dir'], file_path=Attributes['Dir'])   
            self.MakeFile(Attributes)
            with open(OutFile, 'a') as outfile:
                outfile.write(str(self.Frame) + ' ' +  ' '.join(str(item) for item in self.bin_cents) +'\n')

        elif self.Type == 'Hetero':
            from Sapphire.IO import OutputInfoHetero as Out  # Case 3
          
            Attributes = getattr(Out, str('hepair_distance')) #Loads in the write information for the object 
            OutFile = self.System['base_dir'] + Attributes['Dir'] + Attributes['File']
            self.ensure_dir(base_dir=self.System['base_dir'], file_path=Attributes['Dir'])   
            self.MakeFile(Attributes)
            with open(OutFile, 'a') as outfile:
                outfile.write(str(self.Frame) + ' ' +  ' '.join(str(item) for item in self.a) +'\n')
      
        
            Attributes = getattr(Out, str('hepair_distancespace')) #Loads in the write information for the object                   
            OutFile = self.System['base_dir'] + Attributes['Dir'] + Attributes['File']
            self.ensure_dir(base_dir=self.System['base_dir'], file_path=Attributes['Dir'])   
            self.MakeFile(Attributes)
            with open(OutFile, 'a') as outfile:
                outfile.write(str(self.Frame) + ' ' +  ' '.join(str(item) for item in self.bin_cents) +'\n')