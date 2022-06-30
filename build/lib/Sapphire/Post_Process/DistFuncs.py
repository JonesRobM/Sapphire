import numpy as np
import os

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
    return (np.average(positions))

def get_subspecieslist(specie, elements, positions):
    Temp = np.column_stack((elements,positions))
    Temp = [x for x in Temp if x[0] == specie]
    return np.array(np.delete(Temp,0,1), dtype = np.float64)

def Euc_Dist(positions, homo = False, specie = None, elements = None):
    
    if homo == False:
        Distances=[]
        for i in range(len(positions)-1):
            for j in range(i+1,len(positions)):
                Euc = distance(positions[i],positions[j])

                Distances.append(Euc)
        return Distances
    
    elif homo:
        Distances = []
        Temp = get_subspecieslist(specie, elements, positions)
        if (len(Temp)>1) is False:
            return None
        else:
            for i in range(len(Temp)-1):
                for j in range(i+1,len(Temp)):
                    Euc = distance(Temp[i],Temp[j])

                    Distances.append(Euc)
            return Distances
    else:
        print("Variables used were:\n%s\n%s\n%s\n" %(homo, specie, (elements[0], elements[1])))
        raise TypeError("Euc_Dist function has encountered an error.\n")
        
        
def Hetero(positions, species, elements):
        
    """ Robert
    
    Note that no species need to be defined for this function as it is understood that LoDiS
    only has provision for mono/bimetallic systems (for the time being) although this
    function could be further generalised (albeit it a potential cost to computation time).
    """
    
    TempA = get_subspecieslist(species[0], elements, positions)
    TempB = get_subspecieslist(species[1], elements, positions)
    try:
        np.shape(TempA)[1]
        try:
            np.shape(TempB)[1]
            Dist=[]
            for a in TempA:
                Temp = [ distance(a,b) for b in TempB]
                Dist.append(Temp)
            return Dist
        except IndexError:
            Dist=[]
            for x in TempA:
                Dist.append( [distance(x, TempB) ])
            return Dist
            print("You have only one of a specific atom type in your simulation. I hope that this is correct.", "\n")
    except IndexError:
        try:
            np.shape(TempB)[1]           
            return [ distance(TempA, b) for b in TempB ]
            print("You have only one of a specific atom type in your simulation. I hope that this is correct.", "\n")
        except IndexError:
            print("You only have two atoms.\nIs this correct?", "\n")
            return None
        

class CoM_Dist():
    
    def __init__(self, System, Positions, CoM = None, Type = False, Specie = None, Elements = None, Frame = None):
        self.System = System
        self.Positions = Positions
        self.CoM =Positions
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
            self.Dist = np.array([distance(x, self.CoM) for x in Temp])
            self.CoM = get_CoM(Temp)
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
        
        if not self.Type == 'Hetero':
            for i, atom1 in enumerate(self.Positions):
                for j in range(self.Res):
                    r1 = j * self.dr #Inner radius for the spherical shell
                    r2 = r1 + self.dr #Outer radius increased by increment dr
                    v1 = 4.0 / 3.0 * np.pi * r1**3
                    v2 = 4.0 / 3.0 * np.pi * r2**3
                    self.Volumes[j] += v2 - v1 #Volume to consider when evaluating distribution
            
                for atom2 in self.Positions[i:]:
                    self.Distance = distance(atom1, atom2)
                    index = int(self.Distance / self.dr)
                    if 0 < index < self.Res:
                        self.G[index] += 2 #Identifies when there is an atom at this distance
                        
            for i, value in enumerate(self.G):
                self.G[i] = value / self.Volumes[i] #Rescaling the distribution with respect to enclosing volume
        elif self.Type == 'Hetero':
            TempA = get_subspecieslist(self.Species[0], self.Elements, self.Positions)
            TempB = get_subspecieslist(self.Species[1], self.Elements, self.Positions)
            for i, atom1 in enumerate(TempA):
                for j in range(self.Res):
                    r1 = j * self.dr #Inner radius for the spherical shell
                    r2 = r1 + self.dr #Outer radius increased by increment dr
                    v1 = 4.0 / 3.0 * np.pi * r1**3
                    v2 = 4.0 / 3.0 * np.pi * r2**3
                    self.Volumes[j] += v2 - v1 #Volume to consider when evaluating distribution
    
                for atom2 in TempB:
                    self.Distance = distance(atom1, atom2)
                    index = int(self.Distance / self.dr)
                    if 0 < index < self.Res:
                        self.G[index] += 2 #Identifies when there is an atom at this distance
            
        
            for i, value in enumerate(self.G):
                self.G[i] = value / self.Volumes[i] #Rescaling the distribution with respect to enclosing volume
                
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
                pass
        elif self.Type == 'Hetero':
            try:
                self.distances = Hetero(self.Positions, self.Specie, self.Elements)
            except Exception as e:
                pass
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