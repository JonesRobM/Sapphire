import os
import numpy as np

def distance(a, b):
    
    dx = abs(a[0] - b[0])
     
    dy = abs(a[1] - b[1])
     
    dz = abs(a[2] - b[2])
 
    return np.sqrt(dx**2 + dy**2 + dz**2)

def get_subspecieslist(specie, elements, positions):
    Temp = np.column_stack((elements,positions))
    Temp = [x for x in Temp if x[0] == specie]
    return np.array(np.delete(Temp,0,1), dtype = np.float64)

def CoMDist(positions, CoM = None, homo = False, specie = None, elements = None):
    
    if homo == False:
        return [distance(x, CoM) for x in positions]
    elif homo:
        Temp = get_subspecieslist(specie, elements, positions)
        CoM = get_CoM(Temp)
        return [distance(x, CoM) for x in Temp]

def get_CoM(positions):
    return (np.average(positions, axis = 0))

class Gyration():
    
    def __init__(self, System = None, Positions = None, Type = False, 
                 Metal = None, Elements = None, Masses = None, Frame = None):
        self.System = System
        self.Positions = Positions
        self.Type = Type
        self.Metal = Metal
        self.Elements = Elements
        self.Masses =  Masses
        self.Frame = Frame
        self.calculate()
        self.write()
        
    def ensure_dir(self, base_dir='', file_path=''):

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
            Distances = CoMDist(CoM = get_CoM(
                                    get_subspecieslist(self.Metal, self.Elements, self.Positions)), 
                                         homo=True, positions = self.Positions, 
                                         specie = self.Metal, elements=self.Elements)
            
            self.Gyrate = np.sqrt(sum(Distances)/len(Distances))
        
        elif self.Type == 'Full':
            Distances = np.array(CoMDist(self.Positions, CoM = get_CoM(self.Positions)))
            S = (Distances**2) *  np.array(self.Masses)
            self.Gyrate = np.sqrt( sum(S) / sum(self.Masses) )
    
    def write(self):
        
        if self.Type == 'Full':
            from Sapphire.IO import OutputInfoFull as Out  # Case 1
          
            Attributes = getattr(Out, str('gyration')) #Loads in the write information for the object 
            OutFile = self.System['base_dir'] + Attributes['Dir'] + Attributes['File']
            self.ensure_dir(base_dir=self.System['base_dir'], file_path=Attributes['Dir'])   
            self.MakeFile(Attributes)
            with open(OutFile, 'a') as outfile:
                outfile.write(str(self.Frame) + ' ' +  str(self.Gyrate) + '\n')

        if self.Type == 'Homo':
            from Sapphire.IO import OutputInfoHomo as Out  # Case 2
          
            Attributes = getattr(Out, str('hogyration')) #Loads in the write information for the object 
            OutFile = self.System['base_dir'] + Attributes['Dir'] + Attributes['File'] + self.Metal
            self.ensure_dir(base_dir=self.System['base_dir'], file_path=Attributes['Dir'])   
            self.MakeFile(Attributes)
            with open(OutFile, 'a') as outfile:
                outfile.write(str(self.Frame) + ' ' +  str(self.Gyrate) +'\n') 

class Stat_Radius():
    
    def __init__(self, System, Positions, Frame):
        self.System = System
        self.Positions = Positions
        self.Frame = Frame
        self.calculate()
        self.write()
    def ensure_dir(self, base_dir='', file_path=''):

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
        
        CoM = np.array( CoMDist(self.Positions, get_CoM(self.Positions)) )**2
        self.Rad = np.sqrt( (5/3)*(1/len(self.Positions))*sum(CoM) )

    def write(self):
        
        from Sapphire.IO import OutputInfoFull as Out  # Case 1
      
        Attributes = getattr(Out, str('stat_radius')) #Loads in the write information for the object 
        OutFile = self.System['base_dir'] + Attributes['Dir'] + Attributes['File']
        self.ensure_dir(base_dir=self.System['base_dir'], file_path=Attributes['Dir'])   
        self.MakeFile(Attributes)
        with open(OutFile, 'a') as outfile:
            outfile.write(str(self.Frame) + ' ' + str(self.Rad) +'\n')
