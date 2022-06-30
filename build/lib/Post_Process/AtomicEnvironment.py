import numpy as np
import scipy.sparse as spa
import os

class LAE():

    def __init__(self, System = None, Frame = None, HeAdj = None, Metal = None, Species = None):
        
        """
        Robert:
            This class function faacilitates the computation of heterogeneous 
            atomic quantities which requires adjacenency information to be fed in
            from the Sapphire.Post_Process.Adjacent module.
            
            Args:
                System : Type - Dict
                    Description - Base system information regarding directories.
                    Not necessary for separate use outside of Sapphire core, 
                    but writing output is not possible without reference directories.
                    
        """
        
        self.System = System
        self.Frame = Frame
        self.HeAdj = HeAdj.todense()
        self.Species = Species
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
        
    def LAE(self):
        
        a = np.shape(self.HeAdj)
        self.MatA = np.zeros(13) #Initialise w.r.t total time
        for i in range(a[0]):
            T = sum(self.HeAdj[i])
            try:
                self.MatA[T] += 1
            except IndexError:
                pass
            
        self.MatB = np.zeros(13) #Initialise w.r.t total time
        for i in range(a[1]):
            T = sum(self.HeAdj[:,i])
            try:
                self.MatB[T] += 1
            except IndexError:
                pass
        
    def write(self):

        from Sapphire.IO import OutputInfoHetero as Out  # Case 3
        #Write object for the homo CoM distances
        Attributes = getattr(Out, str('lae')) #Loads in the write information for the object                  
        OutFile = self.System['base_dir'] + Attributes['Dir'] + Attributes['File'] + self.Species[0]
        self.ensure_dir(base_dir=self.System['base_dir'], file_path=Attributes['Dir']) + self.Species[0]
        self.MakeFile(Attributes)
        with open(OutFile, 'a') as outfile:
            outfile.write(str(self.Frame) + ' ' +  ' '.join(str(item) for item in self.MatA) +'\n') 

        OutFile = self.System['base_dir'] + Attributes['Dir'] + Attributes['File'] + self.Species[1]
        self.ensure_dir(base_dir=self.System['base_dir'], file_path=Attributes['Dir']) + self.Species[1]
        self.MakeFile(Attributes)
        with open(OutFile, 'a') as outfile:
            outfile.write(str(self.Frame) + ' ' +  ' '.join(str(item) for item in self.MatB) +'\n') 
                
class Mix():

    def __init__(self, System = None, Frame = None, Adj1 = None, Adj2 = None, HeAdj = None, 
                 EleNN = None, lae = None, HomoBonds = None, HeteroBonds = None, Mix = None,
                 Metal = None, Species = None):
        
        """
        Robert:
            This class function faacilitates the computation of heterogeneous 
            atomic quantities which requires adjacenency information to be fed in
            from the Sapphire.Post_Process.Adjacent module.
            
            Args:
                System : Type - Dict
                    Description - Base system information regarding directories.
                    Not necessary for separate use outside of Sapphire core, 
                    but writing output is not possible without reference directories.

        
        """
        
        self.System = System
        self.Frame = Frame
        self.Adj1 = Adj1
        self.Adj2 = Adj2
        self.HeAdj = HeAdj
        self.HeteroBonds = HeteroBonds
        self.Species = Species
        self.Mix = Mix
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
        self.HoBonds = sum(self.Adj1)/2 + sum(self.Adj2)/2
        self.HeBonds = sum(self.HeAdj[0])
        self.Mix_Param =  (self.HoBonds - self.HeBonds) / (self.HoBonds + self.HeBonds)

    def write(self):

        from Sapphire.IO import OutputInfoHetero as Out  # Case 3
        #Write object for the CoM
        Attributes = getattr(Out, str('mix')) #Loads in the write information for the object 
        OutFile = self.System['base_dir'] + Attributes['Dir'] + Attributes['File']
        self.ensure_dir(base_dir=self.System['base_dir'], file_path=Attributes['Dir'])   
        self.MakeFile(Attributes)
        with open(OutFile, 'a') as outfile:
            outfile.write(str(self.Frame) + ' ' +  ' '.join(str(item) for item in self.Mix_Param) +'\n')

        if self.HeteroBonds:
            #Write object for the homo CoM distances
            Attributes = getattr(Out, str('hetero_bonds')) #Loads in the write information for the object                  
            OutFile = self.System['base_dir'] + Attributes['Dir'] + Attributes['File']
            self.ensure_dir(base_dir=self.System['base_dir'], file_path=Attributes['Dir'])   
            self.MakeFile(Attributes)
            with open(OutFile, 'a') as outfile:
                outfile.write(str(self.Frame) + ' ' +  ' '.join(str(item) for item in self.HeBonds) +'\n')
                
        if self.HomoBonds:
            from Sapphire.IO import OutputInfoHomo as Out
            #Write object for the CoMDistances
            Attributes = getattr(Out, str('homo_bonds')) #Loads in the write information for the object                
            OutFile = self.System['base_dir'] + Attributes['Dir'] + Attributes['File']
            self.ensure_dir(base_dir=self.System['base_dir'], file_path=Attributes['Dir'])   
            self.MakeFile(Attributes)
            with open(OutFile, 'a') as outfile:
                outfile.write(str(self.Frame) + ' ' +  ' '.join(str(item) for item in self.HoBonds) +'\n')
        
class Ele_NN():

    def __init__(self, System = None, Frame = None, Adj1 = None, Adj2 = None, HeAdj = None, 
                 EleNN = None, lae = None, HomoBonds = None, HeteroBonds = None, Mix = None,
                 Metal = None, Species = None):
        
        """
        Robert:
            This class function faacilitates the computation of heterogeneous 
            atomic quantities which requires adjacenency information to be fed in
            from the Sapphire.Post_Process.Adjacent module.
            
            System : Type - Dict
            Description - Full Sapphire calculation information regarding base directories and file composition.
        
        
        """
        
        self.System = System
        self.Frame = Frame
        self.Adj1 = Adj1
        self.Adj2 = Adj2
        self.HeAdj = HeAdj
        self.EleNN = EleNN
        self.lae = lae
        self.HomoBonds = HomoBonds
        self.HeteroBonds = HeteroBonds
        self.Species = Species        
        self.Metal = Metal
        self.Mix = Mix
        self.Metal_Index = self.Species.index(self.Metal)
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
        
    def ele_nn(self):
        if self.Metal_Index == 0:
            self.EleNN = self.Adj1 + self.HeAdj[self.Metal_Index]
        elif self.Metal_Index == 1:
            self.EleNN = self.Adj2 + self.HeAdj[self.Metal_Index]
        
    def write(self):

        from Sapphire.IO import OutputInfoHetero as Out  # Case 3