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

"""
Supported=[
        'rdf', 'cna', 'adj', 'pdf', 'pdfhomo', 'agcn', 'nn',
        'SimTime', 'EPot', 'ETot', 'EKin', 'EDelta', 'MeanETot', 'Temp'
           ]
"""

def get_heights(CNA, Masterkey, frame):
    CNA[frame] = list(CNA[frame])
    CNA[frame][1] = list(CNA[frame][1])
    def getkey(item):
        return item[0]
    
    FullSample=[]; Heights=[]
    Temp1=CNA[frame][0]
    for x in Masterkey:
        if x not in Temp1:    
            CNA[frame][0].append(x)
            CNA[frame][1].append(0)
    Sample=[]
    for j in range(len(Masterkey)):
        Sample.append((CNA[frame][0][j], CNA[frame][1][j]))
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

def ensure_dir(file_path=''):
    directory = os.path.dirname(file_path)
    if not os.path.exists(directory):
        os.makedirs(directory)


class Plot_Meta():
    def __init__(self, System = None):
        
        """        
        Robert
            
            Reading user defined inputs for where to find the simulation data,
            where it can be found and the names of files used.
            
            Alse provides the output directroy for plots to be sent to.
        """
        
        self.BigMeta = {}
        
        if System == None:
             self.System = None
             self.Base = ''
             self.Iter = ''
             self.Images = ''
             self.Meta = 'MetaTrial.csv'
        else:
            self.System = System
            try:
                self.Base = System['base_dir']
            except KeyError:
                self.Base = ''
            try:
                self.Iter = System['iter_dir']
            except KeyError:
                self.Iter = ''
                
            try:
                self.Images = System['plot_dir']
                ensure_dir(self.Base + self.Images)
            except KeyError:
                self.Images = ''
            
            try:
                self.Meta = System['meta_name']
            except KeyError:
                self.Meta = 'MetaTrial.csv'
                

                
        with open(self.Base+'Plotting_Info.txt', "w") as f:
            f.write("""
                            
                      _____         _____  _____  _    _ _____ _____  ______ 
                     / ____|  /\   |  __ \|  __ \| |  | |_   _|  __ \|  ____|
                    | (___   /  \  | |__) | |__) | |__| | | | | |__) | |__   
                     \___ \ / /\ \ |  ___/|  ___/|  __  | | | |  _  /|  __|  
                     ____) / ____ \| |    | |    | |  | |_| |_| | \ \| |____ 
                    |_____/_/    \_\_|    |_|    |_|  |_|_____|_|  \_\______|
                                                                      
                     
                                              ____ 
                                             /\__/\ 
                                            /_/  \_\ 
                                            \ \__/ / 
                                             \/__\/ 
                                                                                                                                   
                                """)
        
    
                    
            """
            The next section will be to read through the iterable directories and 
            pull out the metadata files to be used in analysis.
            """

        
        """
        This section reads in the metadata files from each iterable
        directory and then creates a dictionary to average over.
        """
        
        for Object in self.Iter:
            with open(self.Base + Object + '/' + self.Meta, "rb") as file:
                self.BigMeta[str(Object)] = pickle.load(file)
                try:
                    self.BigMeta[str(Object)]['masterkey'].sort()
                except KeyError:
                    continue
                self.Keys = list(self.BigMeta[str(Object)].keys())
                
        self.AverageMeta = {} #Averaged data across simulations
        self.Errors = {} #The standard errors which will be plotted
        
        #This front-loads the reading functions so that the user may either
        #settle with default (full simulation) or user-defined time ranges.
        
        for Item in ['Start', 'End', 'Step', 'Skip', 'Band',
                     'NSpecies', 'NFrames', 'NAtoms']:
            self.AverageMeta[Item] = np.average( [self.BigMeta[x][Item] for x in self.Iter], axis = 0)
            self.Keys.remove(Item)
        
        self.Species = self.BigMeta[self.Iter[0]]['Species']
        self.Elements = self.BigMeta[self.Iter[0]]['Elements']
        
        for Item in ['pdftype', 'Elements', 'Species', 'euc', 'pos']:
            self.Keys.remove(Item)
        
    
    def Average(self):
        
        """
        This function takes in the dictionary of the metadata and averages over
        the quantities found to exist as the relevant keys.
        """
        
        for obj in self.Keys:
            
            if obj.lower() in ['cna', 'masterkey']:
                pass
            else:
                try:
                    
                    Truth = type(self.BigMeta[self.Iter[0]][obj][1]) is list
                    
                except TypeError:
                    print("\nTypeError for %s as it is %s.\n"%(obj, type(self.BigMeta[self.Iter[0]][obj])))
                    Truth = False
                    
                except IndexError:
                    print("\nIndexError for %s.\n"%(obj))
                    Truth = False
                    
                if Truth:
                    
                    if self.Range_Comp(obj):                        
                        TempDat, TempErr = self.Add_Quant_List(obj)
                        self.AverageMeta[obj] = TempDat
                        self.Errors[obj] = TempErr
                        Truth = False
                        
                    else:
                        TempDat, TempErr = self.Add_Quant_List(obj, range(len(self.BigMeta[self.Iter[1]][obj])))

                    
                    self.AverageMeta[obj] = TempDat
                    self.Errors[obj] = TempErr
                    Truth = False
                    
                try:  
                    
                    Truth = type(self.BigMeta[self.Iter[0]][obj][2]) is tuple
                    
                except TypeError:
                    print("\nTypeError for %s as it is %s.\n"%(obj, type(self.BigMeta[self.Iter[0]][obj])))
                    Truth = False
                    
                except IndexError:
                    print("\nIndexError for %s.\n"%(obj))
                    Truth = False

                if Truth:
                    
                    if self.Range_Comp(obj):                        
                        TempDat, TempErr = self.Add_Quant_Tuple(obj)
                        
                    else:
                        TempDat, TempErr = self.Add_Quant_Tuple(obj, range(len(self.BigMeta[self.Iter[1]][obj])))
                        
                    self.AverageMeta[obj] = TempDat
                    self.Errors[obj] = TempErr
                    Truth = False
                    
                try:    
                    
                    Truth = 'float' in str(type(self.BigMeta[self.Iter[0]][obj][1]))
                    
                except TypeError:
                    print("\nTypeError for %s as it is %s.\n"%(obj, type(self.BigMeta[self.Iter[0]][obj])))
                    Truth = False
                    
                except IndexError:
                    print("\nIndexError for %s.\n"%(obj))
                    Truth = False
                    
                if Truth:
                    print("\n%s\n"%obj)
                    for It in self.Iter:
                        self.BigMeta[It][obj] = [ float(x) for x in self.BigMeta[It][obj] if x!= None ]
                    self.AverageMeta[obj] = np.average( [ self.BigMeta[x][obj] for x in self.Iter ], axis = 0)
                    self.Errors[obj] = np.std( [ self.BigMeta[x][obj] for x in self.Iter ], axis = 0)
                    Truth = False
                
        self.AverageMeta['masterkey'] = []
        
        for x in self.BigMeta:
            try:
                for signature in self.BigMeta[x]['masterkey']:
                    if signature not in self.AverageMeta['masterkey']:
                        self.AverageMeta['masterkey'].append(signature)
            except KeyError:
                with open(self.Base+'Plotting_Info.txt', "a") as f:
                    f.write("\n%s\n"%x)
                    f.write(traceback.format_exc())
                        
        for x in self.BigMeta:
            for i in range(len(self.BigMeta[x]['cna'])):
                self.BigMeta[x]['cna'][i] = get_heights(  self.BigMeta[x]['cna'], self.AverageMeta['masterkey'], i)
                    
        self.AverageMeta['masterkey'].sort()


        self.AverageMeta['cna'] = []    
        for i in range(len(self.BigMeta[self.Iter[0]]['cna'])):
            self.AverageMeta['cna'].append( np.average( [ self.BigMeta[x]['cna'][i][1] for x in self.Iter ], axis = 0) )
            self.Errors['cna'] = []

            self.Errors['cna'].append( np.std( [ self.BigMeta[x]['cna'][i][1] for x in self.Iter ], axis = 0) )

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
        
        
        
    def Add_Quant_List(self, Quant, Range = None):
        print("\n%s\n"%Quant)
        
        Val, Err = [], []
        
        if Range is None:
            Range = range(int(self.AverageMeta['Start']), int(self.AverageMeta['End']), int(self.AverageMeta['Step']))
        
        for i in Range:
    
            for It in self.Iter:
                self.BigMeta[It][Quant][i] = [ float(x) for x in self.BigMeta[It][Quant][i] if x!= None ]
            
            TempDat = np.average( ( [self.BigMeta[It][Quant][i] for It in self.Iter ] ), axis = 0)
            
            TempErr = np.std( ( [self.BigMeta[It][Quant][i] for It in self.Iter ] ), axis = 0)
            
            Val.append(TempDat)
            Err.append(TempErr)
                    
        
        return Val, Err

    def Add_Quant_Tuple(self, Quant, Range = None):
        print("\n%s\n"%Quant)
        
        if Range is None:
            Range = range(self.AverageMeta['Start'], self.AverageMeta['End'], self.AverageMeta['Step'])
            
        TempDat = np.empty( (len(Range),len(self.BigMeta[self.Iter[0]][Quant][0])), dtype = object )
        TempErr = np.empty( (len(Range),len(self.BigMeta[self.Iter[0]][Quant][0])), dtype = object )
        
        for i in Range:

            for j in range(len(self.BigMeta[self.Iter[0]][Quant][i])):
            
                TempDat[i][j] = np.average( ( [self.BigMeta[It][Quant][i][j] for It in self.Iter ] ), axis = 0)
                
                TempErr[i][j] = np.std( ( [self.BigMeta[It][Quant][i][j] for It in self.Iter ] ), axis = 0)
            
            
        return TempDat, TempErr
        
    def Range_Comp(self, obj):
        R = range(int(self.AverageMeta['Start']), int(self.AverageMeta['End']), int(self.AverageMeta['Step']))
        if len(self.BigMeta[self.Iter[0]][obj]) == len(R):
            return True
        else:
            return False       

#########################################################################################################################
        
        
    def Lin_Func(Dist_List):
        return np.sqrt( np.average( [ a**2 for a in Dist_List ] ) - np.average(Dist_List)**2 )/ np.average(Dist_List)
    
    
    def Lin_List(Movie, index1, index2):
        List = []
        for T in range(0,len(Movie),10):
            atoms = Movie[T].get_positions()
            List.append(distance(atoms[index1],atoms[index2]))
        return Lin_Func(List)
    
    def Func(Movie, index1):
        A = [ Lin_List(Movie, index1, x) for x in range(0,891) if x!= index1 ]
        return sum(A)/len(A)
    
    def Main_Lind(Movie, N_Atoms, CPUs):
        iterable = range(0,N_Atoms)
        pool = mp.Pool(CPUs)
        function = partial(Func, Movie)
        New_Sample = pool.map(function, iterable)
        pool.close()
        pool.join()
        return New_Sample
    
    def CoMDist(Sim, Frame):
        Pos = read(path+'Sim-'+str(Sim)+'/movie.xyz', index = Frame*10).get_positions()
        CDist = [ distance(BigMeta[Sim]['CoM'][Frame], Pos[i]) for i in range(len(Pos)) ]
        return CDist
    
    def Main_CoM(Sim, CPUs):
        iterable = range(1000)
        pool = mp.Pool(CPUs)
        function = partial(CoMDist, Sim)
        New_Sample = pool.map(function, iterable)
        pool.close()
        pool.join()
        return New_Sample
    
        