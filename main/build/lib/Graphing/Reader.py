import pickle
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import norm
import multiprocessing as mp
from functools import partial
import os
import linecache
import sys
import traceback 
from inspect import getmembers, isfunction

from Graphing import Plot_Funcs
from CNA.Utilities import Logo

"""
Supported=[
        'rdf', 'cna', 'adj', 'pdf', 'pdfhomo', 'agcn', 'nn',
        'SimTime', 'EPot', 'ETot', 'EKin', 'EDelta', 'MeanETot', 'Temp'
           ]
"""

def get_heights_asap3(CNA, Masterkey, frame):
    CNA[frame] = list(CNA[frame])
    CNA[frame][1] = list(CNA[frame][1])
    def getkey(item):
        return item[0]
    
    FullSample=[]; Heights=[]
    Temp1=CNA[frame][0]
    for x in Masterkey:
        if x not in Temp1:    
            CNA[frame][0][x] = 0
    Sample=[]
    for j in Masterkey:
        Sample.append((j, CNA[frame][0][j]))
    FullSample.append(Sample)
    FullSample.sort()
    A,B=zip(*FullSample[0])
    Heights.append(B/np.sum(B))
    Temp = FullSample[0]
    FullCNA = [ item[0] for i, item in enumerate(Temp) ]
    
    Sample = (FullCNA, Heights[0])
    Altern = []
    for x in range(len(FullCNA)):
        Altern.append( ( FullCNA[x],Heights[0][x] ) )
    Sample = sorted(Altern, key = getkey)
    
    return ( [ Sample[x][0] for x in range(len(Sample)) ], [ Sample[x][1] for x in range(len(Sample)) ])

def Get_Heights_Ovito(CNAs, Masterkey, Norm = False):
    
    """ Jones
    
    Arguments:
        filename: The string name of your input xyz file
            Normally something like 'movie.xyz'
        
        Metadata: The dictionary containng the time ordered CNA signatures 
        and the number of observed occurances.
        
        
        MasterKey: The output from calling the Master function.
        This is to do pairwise comparrison for creating full 
        distributions without having to know what the craic is.
            
        Norm: Default - False
        Whether or not the user wishes to normalise the distribution of 
        CNA signatures for each frame in order to perform meaningful
        statistical analysis.
        
        
    Returns:
        
        Heights: np.array(Frames/Skip, len(MasterKey)) The array containing 
        the (if desired) normalised
        distribution of CNA signature occurances. 
        
    """
    

    Heights=np.zeros((len(Metadata),len(Masterkey)))
    

    for frame in range(len(Metadata)):
        Temp = Metadata[frame].keys()
        for x in Masterkey:
            if x not in Temp:
        
                Metadata[frame][x] = 0

            Heights[frame][Masterkey.index(x)] = Metadata[frame][x]
            
            if Norm == True:
                Heights[frame] = Heights[frame]/sum(Heights[frame])
            
    return Heights

def ensure_dir(file_path=''):
    directory = os.path.dirname(file_path+'/')
    if not os.path.exists(directory):
        os.makedirs(directory)
        print("Made a new directory.")


class Read_Meta():
    def __init__(self, System = None):
        
        """        
        Robert
            
            Reading user defined inputs for where to find the simulation data,
            where it can be found and the names of files used.
            
            Alse provides the output directroy for plots to be sent to.
        """
        L = Logo().Logo()


        self.BigMeta = {}
        if System == None:
             self.System = None
             self.Base = ''
             with open(self.Base+'Plotting_Info.txt', "w") as f:
                 f.write(L)
                 f.write(" #"*50)
                 f.write("\n")
             with open(self.Base+'Plotting_Info.txt', "a") as f:
                 f.write("\nNo specified system info to read from. Attempting to load in defaults.\n")
             with open(self.Base+'Plotting_Info.txt', "a") as f:
                 f.write("Base directory set to current working directory.\n")
                 f.write(os.getcwd())
             self.Images = ''
             with open(self.Base+'Plotting_Info.txt', "a") as f:
                 f.write("\nAny images generated will be saved to this directory.\n")
             self.Meta = 'Metadata.csv'
             with open(self.Base+'Plotting_Info.txt', "a") as f:
                 f.write("Will attempt to read input data from 'Metadata.csv'.\n")
             self.single_file = True
             self.Iter = False
             with open(self.Base+'Plotting_Info.txt', "a") as f:
                 f.write("Will only attempt to read a single trajectory file.\n")
             
        else:
            self.System = System
            try:
                self.Base = System['base_dir']
                with open(self.Base+'Plotting_Info.txt', "w") as f:
                    f.write(L)
                    f.write(" #"*50)
                    f.write("\n")
                    f.write("\nBase directory set to %s.\n" %System['base_dir'])
            except KeyError:
                self.Base = ''
                with open(self.Base+'Plotting_Info.txt', "w") as f:
                    f.write(L)
                    f.write(" #"*50)
                    f.write("\n")                    
                    f.write("\nNo base directory given. Will attempt to read from current working directory.\n")
                    f.write(os.getcwd())
                    f.write("\n")
            try:
                self.Iter = System['iter_dir']
                if System['iter_dir'] is False:
                    self.single_file = True
                    self.Iter = False
                    with open(self.Base+'Plotting_Info.txt', "a") as f:
                        f.write("Will only attempt to read a single trajectory file.\n")
                else:
                    with open(self.Base+'Plotting_Info.txt', "a") as f:
                        f.write("Reading multiple simulations over the iterable directories:\n")
                        for obj in self.Iter:
                            f.write("\n%s"%obj)
            except KeyError:
                self.Iter = False
                with open(self.Base+'Plotting_Info.txt', "a") as f:
                    f.write("Will only attempt to read a single trajectory file.\n")                
            try:
                self.Images = System['plot_dir']
                ensure_dir(self.Base + self.Images)
                with open(self.Base+'Plotting_Info.txt', "a") as f:
                    f.write("Will be saving all figures into %s.\n"%(self.Base+self.Images))
            except KeyError:
                self.Images = ''
                with open(self.Base+'Plotting_Info.txt', "a") as f:
                    f.write("Saving all figures to the base directory.\n")
            
            try:
                self.Meta = System['meta_name']
                with open(self.Base+'Plotting_Info.txt', "a") as f:
                    f.write("Will be reading files named '%s'.\n"%self.Meta)
            except KeyError:
                self.Meta = 'Metadata.csv'    
                with open(self.Base+'Plotting_Info.txt', "a") as f:
                    f.write("Reading, by default, from files named '%s'.\n"%self.Meta)
                    
            """
            The next section will be to read through the iterable directories and 
            pull out the metadata files to be used in analysis.
            """

        
        """
        This section reads in the metadata files from each iterable
        directory and then creates a dictionary to average over.
        """
        ensure_dir(self.Base + self.Images)
        with open(self.Base+'Plotting_Info.txt', "a") as f:
            f.write("Verified the existence of a plotting directory.\n")
        self.Iter_Probs = []
        if self.Iter:
            with open(self.Base+'Plotting_Info.txt', "a") as f:
                f.write("Will be executing the software across all iterable directories to average over realisations.\n")
                
            for Object in self.Iter:
                try:
                    with open(self.Base + Object + '/' + self.Meta, "rb") as file:
                        self.BigMeta[str(Object)] = pickle.load(file)
                except FileNotFoundError:
                    self.Iter_Probs.append(Object)
                    with open(self.Base+'Plotting_Info.txt', "a") as f:
                        f.write("\n%s was not found and so this metadata will not be read.\n" %(self.Base + Object + '/' + self.Meta))
                except EOFError:
                    self.Iter_Probs.append(Object)
                    with open(self.Base+'Plotting_Info.txt', "a") as f:
                        f.write("\n%s was not found and so this metadata will not be read.\n" %(self.Base + Object + '/' + self.Meta))
                self.Keys = list(self.BigMeta[str(Object)].keys())
                
                    
            self.AverageMeta = {} #Averaged data across simulations
            self.Errors = {} #The standard errors which will be plotted
            for Probs in self.Iter_Probs:
                self.Iter.remove(Probs)
            
            #This front-loads the reading functions so that the user may either
            #settle with default (full simulation) or user-defined time ranges.
            
            for Item in ['Start', 'End','Step', 'Skip', 'Band',
                         'NSpecies', 'NFrames', 'NAtoms']:
                self.AverageMeta[Item] =  self.BigMeta[self.Iter[0]][Item] 
                self.Keys.remove(Item)
            
            self.Species = self.BigMeta[self.Iter[0]]['Species']
            self.Elements = self.BigMeta[self.Iter[0]]['Elements']
            
            for Item in ['pdftype', 'Elements', 'Species', 'euc', 'pos']:
                self.Keys.remove(Item)
                
        elif self.single_file:
            try:
                with open(self.Base + '/' + self.Meta, "rb") as file:
                    self.Metadata = pickle.load(file)                
            except FileNotFoundError:
                with open(self.Base + 'Plotting_Info.txt', 'a') as file:
                    file.write("\nNo metadata could be read. Terminating programme.\n")
                    sys.exit(0)
            self.Keys = list(self.Metadata.keys())
                    
            for Item in ['Start', 'End', 'Step', 'Skip', 'Band', 'NSpecies', 'NFrames', 'NAtoms']:
                self.Keys.remove(Item)
            
            self.Species = self.Metadata['Species']
            self.Elements = self.Metadata['Elements']
            
            for Item in ['pdftype', 'Elements', 'Species', 'euc', 'pos']:
                self.Keys.remove(Item)
            
    
    def Average(self):
        
        """
        This function takes in the dictionary of the metadata and averages over
        the quantities found to exist as the relevant keys.
        """
        if self.Iter:
            with open(self.Base+'Plotting_Info.txt', 'a') as file:
                file.write("\nAveraging over all of the iterable directories.\n")
            for obj in self.Keys:
            
                try:
                    
                    Truth = type(self.BigMeta[self.Iter[0]][obj][1]) is list
                    
                except TypeError:
                    with open(self.Base+'Plotting_Info.txt', "a") as f:
                        f.write("\nTypeError for %s as it is %s.\n"%(obj, type(self.BigMeta[self.Iter[0]][obj])))
                        
                except IndexError:
                    with open(self.Base+'Plotting_Info.txt', "a") as f:
                        f.write("\nIndexError for %s.\n"%(obj))
                    Truth = False
                    
                if Truth:
                    
                    if self.Range_Comp(obj):                        
                        TempDat, TempErr = self.Add_Quant_List(obj)
                        self.AverageMeta[obj] = TempDat
                        self.Errors[obj] = TempErr
                        Truth = False
                        
                    else:
                        try:
                            TempDat, TempErr = self.Add_Quant_List(obj, range(len(self.BigMeta[self.Iter[0]][obj])))
                        except IndexError:
                            with open(self.Base+'Plotting_Info.txt', "a") as f:
                                f.write("\IndexError raised for %s. Consider adding it yourself or review the metadata.\n"%(obj))                           
                    
                    self.AverageMeta[obj] = TempDat
                    self.Errors[obj] = TempErr
                    Truth = False
                    
                try:  
                    
                    Truth = type(self.BigMeta[self.Iter[0]][obj][2]) is tuple
                    
                except TypeError:
                    with open(self.Base+'Plotting_Info.txt', "a") as f:
                        f.write("\nTypeError for %s as it is %s.\n"%(obj, type(self.BigMeta[self.Iter[0]][obj])))
                    Truth = False
                    
                except IndexError:
                    with open(self.Base+'Plotting_Info.txt', "a") as f:
                        f.write("\nIndexError for %s.\n"%(obj))
                    Truth = False
    
                if Truth:
                    
                    if self.Range_Comp(obj):                        
                        TempDat, TempErr = self.Add_Quant_Tuple(obj)
                        
                    else:
                        TempDat, TempErr = self.Add_Quant_Tuple(obj, range(len(self.BigMeta[self.Iter[0]][obj])))
                        
                    self.AverageMeta[obj] = TempDat
                    self.Errors[obj] = TempErr
                    Truth = False
                    
                try:    
                    
                    Truth = 'float' in str(type(self.BigMeta[self.Iter[0]][obj][1]))
                    
                except TypeError:
                    with open(self.Base+'Plotting_Info.txt', "a") as f:
                        f.write("\nTypeError for %s as it is %s.\n"%(obj, type(self.BigMeta[self.Iter[0]][obj])))
                    Truth = False
                    
                except IndexError:
                    with open(self.Base+'Plotting_Info.txt', "a") as f:
                        f.write("\nIndexError for %s.\n"%(obj))
                    Truth = False
                    
                if Truth:
                    with open(self.Base+'Plotting_Info.txt', "a") as f:
                        f.write("\nCurrently averaging over %s.\n"%obj)
                    for It in self.Iter:
                        self.BigMeta[It][obj] = [ float(x) for x in self.BigMeta[It][obj] if x!= None ]
                    self.AverageMeta[obj] = np.average( [ self.BigMeta[x][obj] for x in self.Iter ], axis = 0)
                    self.Errors[obj] = np.std( [ self.BigMeta[x][obj] for x in self.Iter ], axis = 0)
                    Truth = False
                if 'CoM' in obj:
                    Truth = type(self.BigMeta[self.Iter[0]][obj][1]) is np.ndarray
                    if Truth:
                        with open(self.Base+'Plotting_Info.txt', "a") as f:
                            f.write("\nType found to be array for %s as it is %s.\n"%(obj, type(self.BigMeta[self.Iter[0]][obj])))
                        
                        if self.Range_Comp(obj):
                            try:
                                TempDat, TempErr = self.Add_Quant_List(obj)
                                self.AverageMeta[obj] = TempDat
                                self.Errors[obj] = TempErr
                                Truth = False
                            except IndexError:
                                with open(self.Base+'Plotting_Info.txt', "a") as f:
                                    f.write("\IndexError raised for %s. Consider adding it yourself or review the metadata.\n"%(obj))                                                  
                            
                        else:
                            try:
                                TempDat, TempErr = self.Add_Quant_List(obj, range(len(self.BigMeta[self.Iter[0]][obj])))                
                                self.AverageMeta[obj] = TempDat
                                self.Errors[obj] = TempErr
                                Truth = False                           
                            except IndexError:
                                with open(self.Base+'Plotting_Info.txt', "a") as f:
                                    f.write("\IndexError raised for %s. Consider adding it yourself or review the metadata.\n"%(obj))                                                  

            """
            
            The following code was historically used but is now obsolete:
 
            self.AverageMeta['cna'] = [] 
            self.Errors['cna'] = []
            try:
                for i in range(len(self.BigMeta[self.Iter[0]]['cna'])):
                    self.AverageMeta['cna'].append( np.average( [ self.BigMeta[x]['cna'][i][1] for x in self.Iter ], axis = 0) )
        
                    self.Errors['cna'].append( np.std( [ self.BigMeta[x]['cna'][i][1] for x in self.Iter ], axis = 0) )
            except KeyError:
                with open(self.Base+'Plotting_Info.txt', "a") as f:
                    f.write("\nKeyError raised for CNA averaging. Consider adding it yourself or review the metadata.\n")
            except TypeError:
                with open(self.Base+'Plotting_Info.txt', "a") as f:
                    f.write("\nTypeError raised for CNA averaging. Consider adding it yourself or review the metadata.\n")
    
            """
            
            """
            
            Robert:
                
                The following block gives the option for the user to save the averaged
                metadata and the associated errors for later analysis
                
                This means that they will not have to run the plotting script if they
                wish to revisit old data.
            """
            
            try:
                if self.System['save_meta'] is True:
                    with open(self.Base + 'Metadata.csv', "wb") as file:
                        pickle.dump( self.AverageMeta, file, protocol=pickle.HIGHEST_PROTOCOL )
            except KeyError:
                pass
            
            try:
                if self.System['save_errors'] is True:
                    with open(self.Base + 'Errors.csv', "wb") as file:
                        pickle.dump( self.Errors, file, protocol=pickle.HIGHEST_PROTOCOL )
            except KeyError:
                pass   
            
            return self.AverageMeta, self.Errors
        
        elif self.single_file:
            with open(self.Base+'Plotting_Info.txt', 'a') as file:
                file.write("\nAveraging over all of the iterable directories.\n")
            for obj in self.Keys:
            
                if type(self.Metadata[obj][0]) is list:
                    with open(self.Base+'Plotting_Info.txt', 'a') as file:
                        file.write("\nCurrently sanitising the data for %s.\n"%obj)
                    for t in range(len(self.Metadata[obj])):
                        self.Metadata[obj][t] = np.array(self.Metadata[obj][t], dtype = float)
                elif type(self.Metadata[obj][0]) is tuple:
                    with open(self.Base+'Plotting_Info.txt', 'a') as file:
                        file.write("\nCurrently sanitising the data for %s.\n"%obj)
                else:
                    pass
                    
            return self.Metadata, False
        else:
            with open(self.Base + 'Plotting_Info.txt', 'a') as file:
                file.write("\nNo metadata could be read. Terminating programme.\n")
        
    def Add_Quant_List(self, Quant, Range = None):
        with open(self.Base+'Plotting_Info.txt', "a") as f:
            f.write("\nCurrently adding %s to the metadata.\n"%Quant)
        
        Val, Err = [], []
        
        if Range is None:
            Range = range(int(self.AverageMeta['Start']), int(self.AverageMeta['End']/self.AverageMeta['Step']))
        
        for i in Range:
            try:
            
                TempDat = np.average( ( [self.BigMeta[It][Quant][i] for It in self.Iter ] ), axis = 0)
                
                TempErr = np.std( ( [self.BigMeta[It][Quant][i] for It in self.Iter ] ), axis = 0)
                
                Val.append(TempDat)
                Err.append(TempErr)
            except IndexError:
                pass
            except TypeError:
                TempDat = np.average( ( [np.array(self.BigMeta[It][Quant][i], dtype = float) for It in self.Iter ] ), axis = 0)
                
                TempErr = np.std( ( [np.array(self.BigMeta[It][Quant][i], dtype = float) for It in self.Iter ] ), axis = 0)
                
                Val.append(TempDat)
                Err.append(TempErr)
                    
        
        return Val, Err

    def Add_Quant_Tuple(self, Quant, Range = None):
        if Quant.lower is 'masterkey':
            return None
        with open(self.Base+'Plotting_Info.txt', "a") as f:
            f.write("\nCurrently adding %s to the metadata.\n"%Quant)
        
        if Range is None:
            Range = range(int(self.AverageMeta['Start']), int(self.AverageMeta['End']/self.AverageMeta['Step']))
            
        TempDat = np.empty( (len(Range),len(self.BigMeta[self.Iter[0]][Quant][0])), dtype = object )
        TempErr = np.empty( (len(Range),len(self.BigMeta[self.Iter[0]][Quant][0])), dtype = object )
        
        for i in Range:
            try:
                for j in range(len(self.BigMeta[self.Iter[0]][Quant][i])):
                    
                
                    TempDat[i][j] = np.average( ( [self.BigMeta[It][Quant][i][j] for It in self.Iter ] ), axis = 0)
                    
                    TempErr[i][j] = np.std( ( [self.BigMeta[It][Quant][i][j] for It in self.Iter ] ), axis = 0)

            except TypeError:
                with open(self.Base+'Plotting_Info.txt', "a") as f:
                    f.write("\IndexError raised for %s. Consider adding it yourself or review the metadata.\n"%(Quant))
            except IndexError:
                with open(self.Base+'Plotting_Info.txt', "a") as f:
                    f.write("\IndexError raised for %s. Consider adding it yourself or review the metadata.\n"%(Quant))
            
            
        return TempDat, TempErr
        
    def Range_Comp(self, obj):
        R = range(int(self.AverageMeta['Start']), int(self.AverageMeta['End']/self.AverageMeta['Step']))
        if len(self.BigMeta[self.Iter[0]][obj]) == len(R):
            return True
        else:
            return False           