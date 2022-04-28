import numpy as np
from scipy.stats import norm
import os

class Gauss():
    
    """ Robert
    
    This class contains three flavours of kernel desnity estimators.
    Functionally, these attempt to estimate a distribution function 
    for a data set (at this stage, only able to work in 1D but could
    be extended to N dimensions {in theory}).
    
    Essentially, these take a given data point and assign it a weight
    at each point on a pre-defined grid. This weight is determined by the
    type of function used and the proximity of the data point to the position 
    on the pre-defined grid.
    
    As this code is developed, it is likely that each of these functions will
    take additional arguments to determine the cut-off and the desnity of
    grid points.
    
    Note that given the way in which these functions are created and designed, 
    are already normalised.
    
    Note that additional kernels are likely to be implemented into this module
    in the fullnes of time as and when I get around to developing and testing them.
    However; at the moment, these three below appear to be sufficient in efficiency
    and accuracy for the task at hand. But variety IS the spice of life.
    """
    
    def __init__(self, Data = None, Band = None, Ele = None,
                 Type=False, Space = None, System = None, Frame = None):
        self.Space = Space
        self.Band = Band
        self.Data = Data
        self.Type = Type
        self.System = System
        self.Frame = Frame
        self.Ele = Ele

        if self.Space is None:
            self.Space = np.linspace(0, 6.0, 200)

        self.Data=np.array(self.Data)
        self.Data = np.array([a for a in Data if a < 6.5])
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
        
        """ Robert
        
        Data: At this stage, this expects the input of a 1D vector containing
        the data to be iterated over. This is expected to be the output of one
        of the functions from the "Distances.py" module. Although, in theory,
        this could be an arbitrary vector of data.
        
        Band: The bandwidth assigned to this kernel. Increasing the bandwidth results 
        in a smoother looking distribution although runs the risk of missing important
        features of the true distribution. Decreasing it tends towards the behaviour
        of each data point being a dirac-delta peak at its precise location in the limit
        of Band tends to zero.
        This would have the effect of re-paramterising the raw data.
        Good results have come from using band ~ 0.05 for estimating a PDDF. Although,
        one may wish to tune this parameter depending on how they wish to present their
        data.
        
        This particular function assigns the weight according to an underlying
        Gaussian distribution. I.e., the weight that a given data point has is 
        Gaussian distributed about the position on the grid under consideration.
        """

        A=[] 
 
        for i in range(len(self.Data)):
            A.append(norm.pdf(self.Space, self.Data[i], self.Band))
        self.Density = np.asarray(np.sum(A, axis=0))
        self.Density = self.Density / np.trapz( self.Density, self.Space) #For normalisation purposes
        self.Density[np.where(self.Density < 0.01)] = 0 #This smooths near zero elements into zeroes so that minima may be found - change at your own peril

        Min = (np.diff(np.sign(np.diff(self.Density))) > 0).nonzero()[0] + 1 # local min
        try:
            self.R_Cut = self.Space[Min][np.where(self.Space[Min]>3)][0] #We expect a minimum in this region
            
        except Exception as e:
            return None
        
    def ReturnRCut(self):
        return float(self.R_Cut)
            
    def write(self):
        
        if self.Type == 'Full':
            from Sapphire.IO import OutputInfoFull as Out  # Case 1
          
            Attributes = getattr(Out, str('pdf')) #Loads in the write information for the object 
            OutFile = self.System['base_dir'] + Attributes['Dir'] + Attributes['File']
            self.ensure_dir(base_dir=self.System['base_dir'], file_path=Attributes['Dir'])   
            self.MakeFile(Attributes)
            with open(OutFile, 'a') as outfile:
                outfile.write(str(self.Frame) + ' ' +  ' '.join(str(item) for item in self.Density) +'\n')
     
            Attributes = getattr(Out, str('pdfspace')) #Loads in the write information for the object         
            
            OutFile = self.System['base_dir'] + Attributes['Dir'] + Attributes['File']
            self.ensure_dir(base_dir=self.System['base_dir'], file_path=Attributes['Dir'])   
            self.MakeFile(Attributes)
            with open(OutFile, 'a') as outfile:
                outfile.write(str(self.Frame) + ' ' +  ' '.join(str(item) for item in self.Space) +'\n')
          
            Attributes = getattr(Out, str('rcut')) #Loads in the write information for the object         
            
            OutFile = self.System['base_dir'] + Attributes['Dir'] + Attributes['File']
            self.ensure_dir(base_dir=self.System['base_dir'], file_path=Attributes['Dir'])   
            self.MakeFile(Attributes)
            with open(OutFile, 'a') as outfile:
                outfile.write(str(self.Frame) + ' ' + str(self.R_Cut) +'\n')                

        elif self.Type == 'Homo':
            from Sapphire.IO import OutputInfoHomo as Out  # Case 2
          
            Attributes = getattr(Out, str('hopdf')) #Loads in the write information for the object 
            OutFile = self.System['base_dir'] + Attributes['Dir'] + Attributes['File']+self.Ele
            self.ensure_dir(base_dir=self.System['base_dir'], file_path=Attributes['Dir'])   
            self.MakeFile(Attributes)
            with open(OutFile, 'a') as outfile:
                outfile.write(str(self.Frame) + ' ' +  ' '.join(str(item) for item in self.Density) +'\n')
        
            Attributes = getattr(Out, str('hopdfspace')) #Loads in the write information for the object         
            
            OutFile = self.System['base_dir'] + Attributes['Dir'] + Attributes['File']+self.Ele
            self.ensure_dir(base_dir=self.System['base_dir'], file_path=Attributes['Dir'])   
            self.MakeFile(Attributes)
            with open(OutFile, 'a') as outfile:
                outfile.write(str(self.Frame) + ' ' +  ' '.join(str(item) for item in self.Space) +'\n')

            Attributes = getattr(Out, str('hocut')) #Loads in the write information for the object         
            
            OutFile = self.System['base_dir'] + Attributes['Dir'] + Attributes['File']+self.Ele
            self.ensure_dir(base_dir=self.System['base_dir'], file_path=Attributes['Dir'])   
            self.MakeFile(Attributes)
            with open(OutFile, 'a') as outfile:
                outfile.write(str(self.Frame) + ' ' + str(self.R_Cut) +'\n')
       
        elif self.Type == 'Hetero':
            from Sapphire.IO import OutputInfoHetero as Out  # Case 3

            Attributes = getattr(Out, str('hepdf')) #Loads in the write information for the object 
            OutFile = self.System['base_dir'] + Attributes['Dir'] + Attributes['File']
            self.ensure_dir(base_dir=self.System['base_dir'], file_path=Attributes['Dir'])   
            self.MakeFile(Attributes)
            with open(OutFile, 'a') as outfile:
                outfile.write(str(self.Frame) + ' ' +  ' '.join(str(item) for item in self.Density) +'\n')
                
            
            Attributes = getattr(Out, str('hepdfspace')) #Loads in the write information for the object         
            
            OutFile = self.System['base_dir'] + Attributes['Dir'] + Attributes['File']
            self.ensure_dir(base_dir=self.System['base_dir'], file_path=Attributes['Dir'])   
            self.MakeFile(Attributes)
            with open(OutFile, 'a') as outfile:
                outfile.write(str(self.Frame) + ' ' +  ' '.join(str(item) for item in self.Space) +'\n')
 
            
            Attributes = getattr(Out, str('hecut')) #Loads in the write information for the object         
            
            OutFile = self.System['base_dir'] + Attributes['Dir'] + Attributes['File']
            self.ensure_dir(base_dir=self.System['base_dir'], file_path=Attributes['Dir'])   
            self.MakeFile(Attributes)
            with open(OutFile, 'a') as outfile:
                outfile.write(str(self.Frame) + ' ' + str(self.R_Cut) +'\n')
            

"""
class Uniform():
        
         Robert
        
        Data: At this stage, this expects the input of a 1D vector containing
        the data to be iterated over. This is expected to be the output of one
        of the functions from the "Distances.py" module. Although, in theory,
        this could be an arbitrary vector of data.
        
        Band: The bandwidth assigned to this kernel. Increasing the bandwidth results 
        in a smoother looking distribution although runs the risk of missing important
        features of the true distribution. Decreasing it tends towards the behaviour
        of each data point being a dirac-delta peak at its precise location in the limit
        of Band tends to zero.
        This would have the effect of re-paramterising the raw data.
        
        This assigns a uniform weight to a data point if it exists within an interval
        centred on the grid point under consideration. The weight and interval width are
        intrinsically linked for the purposes of distribution normalisation.
        
        Fine details of PDDFs (including peak splitting) has been best observed with 
        a bandwidth of ~ 0.25.
    
    
    def __init__(self, Data, Band, mon=False):
    
    
        Space=np.linspace(2,int(max(Data)/2),300); Density = []
        for i in Space:
            Density.append(minifunc(Data, Band, i))
        if mon == False:
            Min = (np.diff(np.sign(np.diff(Density))) > 0).nonzero()[0] + 1 # local min
            R_Cut = Space[Min[0]]
            return Space, Density, R_Cut
            
        elif mon == True:
            Min = (np.diff(np.sign(np.diff(Density))) > 0).nonzero()[0] + 1 # local min
            R_Cut = Space[Min[0]]
            return Space, Density, R_Cut
        
    def minifunc(Data,Band,i):
        X = Data - i
        A1 = -0.5*Band <= X
        A2 = X <= 0.5*Band
        Temp = np.multiply(A1, A2)
        Temp = Temp/Band
        return np.sum(Temp)/(len(Data)*Band)
            
    
class Epan():
    
    def __init__(self, Space, Density, Data, Band):
        self.Space = Space
        self.Density = Density
        
        
         Robert
        
        Data: At this stage, this expects the input of a 1D vector containing
        the data to be iterated over. This is expected to be the output of one
        of the functions from the "Distances.py" module. Although, in theory,
        this could be an arbitrary vector of data.
        
        Band: The bandwidth assigned to this kernel. Increasing the bandwidth results 
        in a smoother looking distribution although runs the risk of missing important
        features of the true distribution. Decreasing it tends towards the behaviour
        of each data point being a dirac-delta peak at its precise location in the limit
        of Band tends to zero.
        This would have the effect of re-paramterising the raw data.
        
        This particular function utilises the Epanechnikov convention for assigning weights
        to each data point. In essence, this creates a small semi-circle of weight around 
        each grid point to weight the surroudning data points by.
        
        Testing has had good results for a bandwidth of 0.25 when analysing PDDFs.
        
        
    def calculate(self):
   
        Space=np.linspace(0,8.0,400);Density=[]
        for i in Space:
            P=0
            for j in range(len(Data)):
                X = (Data[j]-i)/Band
                P+=0.75*max(1-X**2,0)
            
            Density.append(P/(len(Data)*Band))
        return Space, Density
"""