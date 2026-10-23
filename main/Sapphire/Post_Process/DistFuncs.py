import numpy as np
from scipy.spatial import ConvexHull, QhullError, cKDTree
from scipy.spatial.distance import cdist, pdist
from Sapphire.Utilities import errors
import os
from Sapphire.Utilities.log import get_logger

log = get_logger('Sapphire.Post_Process.DistFuncs')

# Above this many atoms nothing builds the all-pairs list or an N x N array: N(N-1)/2
# float64s is 1.6 GB at 20k atoms and 176 GiB for the 217k-atom Big_Rod. Larger systems
# take the KD-tree routes below, which find exactly the same pairs within each cutoff.
# At or below it the original all-pairs code runs unchanged.
LARGE_N = 20000


def is_large(n_atoms, n_other=0):
    return n_atoms + n_other > LARGE_N


def pairs_within(positions, r_max, other=None):
    """(i, j, d) for every pair closer than ``r_max``, from a KD-tree; memory ~ pairs found.

    Without ``other``: pairs within ``positions``, i < j. With it: i indexes ``positions``,
    j indexes ``other``. Strictly ``d < r_max``, as the all-pairs comparisons are.
    """
    a = np.asarray(positions, dtype=float)
    if other is None:
        ij = cKDTree(a).query_pairs(r_max, output_type='ndarray')
        i, j = ij[:, 0], ij[:, 1]
        d = np.linalg.norm(a[i] - a[j], axis=1)
    else:
        b = np.asarray(other, dtype=float)
        m = cKDTree(a).sparse_distance_matrix(cKDTree(b), r_max, output_type='coo_matrix')
        i, j, d = m.row, m.col, m.data
    keep = d < r_max
    return i[keep], j[keep], d[keep]


def max_pair_distance(positions):
    """Largest pair distance (the cluster's diameter). For large systems only the convex
    hull's vertices are compared: the farthest pair always lies on the hull."""
    p = np.asarray(positions, dtype=float)
    if is_large(len(p)):
        try:
            p = p[ConvexHull(p).vertices]
        except QhullError:      # flat or linear: fall through to all pairs of what's left
            pass
    return float(pdist(p).max())


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

def Euc_Dist(positions, homo = False, specie = None, elements = None, r_max = None):
    """All pair distances |r_i - r_j| for i < j (row-major order), as a 1-D array.

    With ``homo=True`` only atoms of ``specie`` are considered; returns None if fewer than two.
    Equivalent to the original double loop but via ``scipy.spatial.distance.pdist``.
    With ``r_max``, a large system returns only the distances below it, in no set order
    (enough for callers that discard the rest anyway, like the PDF's kernel estimate).
    """
    if homo:
        positions = get_subspecieslist(specie, elements, positions)
        if len(positions) < 2:
            return None
    p = np.asarray(positions, dtype=float)
    if r_max is not None and is_large(len(p)):
        return pairs_within(p, r_max)[2]
    return pdist(p)

def Hetero(positions, species, elements, r_max = None):
    """Distances between every atom of species[0] and every atom of species[1]: (N_A, N_B) array.

    With ``r_max``, a large system returns only the distances below it, as a 1-D array.
    """
    TempA = np.atleast_2d(get_subspecieslist(species[0], elements, positions))
    TempB = np.atleast_2d(get_subspecieslist(species[1], elements, positions))
    if r_max is not None and is_large(len(TempA), len(TempB)):
        return pairs_within(TempA, r_max, TempB)[2]
    return cdist(TempA, TempB)


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
            A = np.atleast_2d(get_subspecieslist(self.Species[0], self.Elements, self.Positions))
            B = np.atleast_2d(get_subspecieslist(self.Species[1], self.Elements, self.Positions))
            r = pairs_within(A, self.R_Cut, B)[2] if is_large(len(A), len(B)) else cdist(A, B).ravel()
            n_ref = len(A)
        else:
            P = np.asarray(self.Positions, dtype=float)
            r = pairs_within(P, self.R_Cut)[2] if is_large(len(P)) else pdist(P)
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