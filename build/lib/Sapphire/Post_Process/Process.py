#Pre-designed pythonic libraries. Indepedent of Sapphire
import numpy as np
import time
import multiprocessing as mp
from inspect import getmembers, isfunction
from ase.io import read
import os

#Objects written specificaly for processing data
from Sapphire.Post_Process import Adjacent, Kernels, DistFuncs, Stats

#Common Neighbour Analysis specific functions
from Sapphire.CNA import FrameSignature, Utilities

#General purpose utility functions for parsing and tidying
from Sapphire.Utilities import Initial, Supported, System_Clean

class Process(object):

    def __init__(self, System=None, Quantities=None,
                 Pattern_Input=False):
        
        self.tick = time.time()

        self.System = System
        self.Quantities = Quantities
        if Pattern_Input:
            self.Pattern_Input = Pattern_Input

        self.filename = System['base_dir']+System['movie_file_name']
        Initial.Logo(self.System['base_dir'])._write_()
        Initial.Info(self.System['base_dir'])._write_()

        self.Tbar = False


        self.result_cache = {}
        """
        This list contains unique strings which are liable
        to require smaller storage objects in the metadata
        """
        self.T = time.time()

        """
        Initialise the CNA signature masterkey
        """

        try:
            if self.Pattern_Input['BULK_MASTERKEY']:
                self.Masterkey = Utilities.CNA_Masterkey().Key()

            else:
                self.Masterkey = []
        except KeyError:
            self.Masterkey = Utilities.CNA_Masterkey().Key()
        self.Initialising()
        #self.run_pdf()
        self.run_core()

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
        self.out = self.System.System['base_dir'] + Attributes['Dir'] + Attributes['File']

        if not os.path.isfile(self.out):
            with open(self.System.System['base_dir'] + Attributes['Dir'] + Attributes['File'], 'w') as out:
                out.close()
        else:
            pass

##############################################################################

    def Initialising(self):


        self.Calc_Quants = {}
        # This next sub-block instantiates the system object.

        self.System = System_Clean._Clean_System(self.System)


        self.Time = int(
            (self.System.System['End'] - self.System.System['Start']) / (self.System.System['Step']))
        self.Step = self.System.System['Step']
        self.Skip = self.System.System['Skip']
        self.End = self.System.System['End']
        self.Start = self.System.System['Start']
        self.Base = self.System.System['base_dir']

##############################################################################

        with open(self.System.System['base_dir']+'Sapphire_Info.txt', "a") as f:
            f.write("Initialising system environment took %.3f seconds.\n" %
                    (time.time()-self.tick))

##############################################################################

        tick = time.time()

        """
        Robert:

            This next block loads up the first frame of the trajectory and sets some initial file parameters and decides how to treat
            poly-metallic or mono-metallic systems depending on the user input.

            Note that this is NOT robust for growth simulations.

            It is a future goal to have the calculation react dynamically to a grand ensemble style of simulation
            but that is not present on this version.

        """

##############################################################################

        with open(self.System.System['base_dir']+'Sapphire_Info.txt', "a") as f:
            f.write("Loading in the dataset to be analysed.\n")
            f.write("Be aware that this may take a while for a large file.\n")
        Read_Time = time.time()
        self.Dataset = read(self.filename, index=':')
        with open(self.System.System['base_dir']+'Sapphire_Info.txt', "a") as f:
            f.write("Opened the dataset in %.3f seconds.\n" %
                    (time.time()-Read_Time))
        self.all_positions = self.Dataset[0].get_positions()
        self.max_dist = max(DistFuncs.Euc_Dist(self.all_positions))
        del(self.all_positions)
        
        #self.metadata['comspace'] = np.linspace(0, self.max_dist, 200)   #WRITE_REQUIRED
        
        #self.all_atoms contains the chemical symbols for all frames
        self.all_atoms = [ 
            self.Dataset[t].get_chemical_symbols() for t in range(
                self.System.System['Start'], self.System.System['End'], self.System.System['Step']
            )
        ] 

        used = set()
        self.Species = [x for x in self.all_atoms[0]
                        if x not in used and (used.add(x) or True)]

        self.NAtoms = [len(self.all_atoms[t]) for t in range(len(self.all_atoms))]

        tick = time.time()
        self.Ele = [] #This gives a convolved form of the elements - can be easier to store
        for atoms in self.all_atoms:
            Temp = []
            i = 0
            for atom in atoms:
                try:
                    if atom == Temp[i][0]:
                        Temp[i][1] += 1
                    else:
                        Temp.append([atom, 1])
                        i += 1
                except IndexError:
                    Temp.append([atom, 1])
            self.Ele.append(Temp)
        self.NSpecies = len(self.Species)
        
##############################################################################

        with open(self.System.System['base_dir']+'Sapphire_Info.txt', "a") as f:
            f.write("Checking user input for calculating homo properties in this run.\n")

        if self.System.System['Homo'] is None:
            with open(self.System.System['base_dir']+'Sapphire_Info.txt', "a") as f:
                f.write("No bimetallic properties will be calculated in this run.\n")

        else:
            with open(self.System.System['base_dir']+'Sapphire_Info.txt', "a") as f:
                f.write("Homo atom properties will be caluclated for %s in this run.\n" % (
                    self.System.System['Homo']))
      
            with open(self.System.System['base_dir']+'Sapphire_Info.txt', "a") as f:
                f.write("Checking user input for hetero atomic species.\n")
            

##############################################################################

            with open(self.System.System['base_dir']+'Sapphire_Info.txt', "a") as f:
                f.write("Initialising Metadata took %.3f seconds.\n" %
                        (time.time() - tick))

##############################################################################

            self.All_Times = list(range(
                self.Start, self.End,self.Step))
            self.Band = self.System.System['Band']


##############################################################################



    def calculate(self, i):
        """

        Robert:

            And we finally get to the meat and potatoes of the calculator.

            This is simply broken down into a series of small blocks.

            Each frame is instantiated by loading in the particular frame of the trajectory.
            While this is time intensive, it saves a lot on memory by not storing an entire trajectory
            in the memory for the entire duration.

            The positions of the atoms are stored in an array and analysis may begin from there.

            Each block is laid out in the following fashion:

                Are each of the required quantities calculated and is this wanted?
                    Y: Calculate and save to metadata by calling a calculator from an external module

                    N: Pass and continue.

            The premise being that the programme will be able to notice that you have not calculated a dependency for a given quantity
            E.g., no Homo quantities in a bimetallic situation
            And will not perform any future calculations which depend on this.

            These quantities are organised by their names as keys which are stored in frame wise metadata dictionaries.

            At the end of the calculation, these frame wise dictionaries are unloaded into a global dictionary and emptied for the next frame.

        """

        self.timer = time.time()

        with open(self.System.System['base_dir']+'Sapphire_Info.txt', "a") as f:
            f.write("\nLoading in atoms for frame %s.\n" % i)
        self.All_Atoms = self.Dataset[i]
        with open(self.System.System['base_dir']+'Sapphire_Info.txt', "a") as f:
            f.write("Loaded the atoms in %.3f seconds.\n" % (time.time()-self.timer))

        self.result_cache['pos'] = self.All_Atoms.get_positions()
        self.result_cache['euc'] = DistFuncs.Euc_Dist(self.result_cache['pos'])
        self.result_cache['syms'] = self.All_Atoms.get_chemical_symbols()
        if 'pdf' in self.Quantities['Full']:
            
            try:
                self.result_cache['FullCut'] = Kernels.Gauss(Data = self.result_cache['euc'], Band = self.Band, 
                                                                Ele = None, Type = 'Full', Space = None, 
                                                                System = self.System.System, Frame = i).ReturnRCut()
            except Exception as e:
                with open(self.Base + 'Sapphire_Errors.log', 'a') as f:
                    f.write('\nException raised while computing Full PDF was: \n%s' % e)

            if 'hepdf' in self.Quantities['Hetero']:
                
                self.result_cache['heteropos'] = DistFuncs.Hetero(self.result_cache['pos'], self.Species,
                                                                  self.result_cache['syms'])
                
                self.result_cache['HeCut'] = Kernels.Gauss(Data = self.result_cache['heteropos'][0], Band = self.Band, 
                                                          Ele = None, Type='Hetero', Space = None, 
                                                          System = self.System.System, Frame = i).ReturnRCut()

                

        if 'hopdf' in self.Quantities['Homo']:
            for x in self.System.System['Homo']:
                self.result_cache['homoed'+x] = DistFuncs.Euc_Dist(positions=self.result_cache['pos'], homo=True, specie=x, elements=self.result_cache['syms'])
                if self.result_cache['homoed'+x] is not None:
                    try:
                        self.result_cache[x+'Cut'] = Kernels.Gauss(self.result_cache['homoed'+x], 
                                              Band = self.Band, Ele = x, 
                                              Type='Homo', Space = None, 
                                              System = self.System.System, Frame = i).ReturnRCut()
                    except Exception as e:
                        with open(self.Base + 'Sapphire_Errors.log', 'a') as f:
                            f.write(
                                '\nException raised while computing Homo rcut was: \n%s' % e)


        if 'moi' in self.Quantities['Full']:
            self.moi = self.All_Atoms.get_moments_of_inertia()
            from Sapphire.Utilities import OutputInfoFull as Out
            Attributes = getattr(Out, str('moi')) #Loads in the write information for the object 
            OutFile = self.System.System['base_dir'] + Attributes['Dir'] + Attributes['File']
            self.ensure_dir(base_dir=self.System.System['base_dir'], file_path=Attributes['Dir'])   
            self.MakeFile(Attributes)
            with open(OutFile, 'a') as outfile:
                outfile.write(str(i) + ' ' +  ' '.join(str(item) for item in self.moi) +'\n')        


##############################################################################

        # All RDF calculations performed in the following block

##############################################################################


        if 'rdf' in self.Quantities['Full']:
            try:
                RDF = DistFuncs.RDF(System = self.System.System,
                    Positions = self.result_cache['pos'], Type = 'Full', Frame = i)
                
            except Exception as e:
                with open(self.Base + 'Sapphire_Errors.log', 'a') as f:
                    f.write('\nException raised while computing Full RDF was: \n%s' % e)
                    
                if self.System.System['Homo'] and 'hordf' in self.Quantities['Homo']:
                    try:
                        for x in self.System.System['Homo']:
                            self.result_cache['homopos'+x] = DistFuncs.get_subspecieslist(x, self.result_cache['syms'], self.result_cache['pos'])
                            HoRDF = DistFuncs.RDF(self.result_cache['homopos'+x], 
                                                  Type = 'Homo', Species = x, Frame = i, System = self.System.System)                          
                    except Exception as e:
                        with open(self.Base + 'Sapphire_Errors.log', 'a') as f:
                            f.write('\nException raised while computing Homo RDF was: \n%s' % e)

                if self.System.System['Hetero'] and 'herdf' in self.Quantities['Hetero']:
                    try:
                        HeRDF = DistFuncs.RDF(self.result_cache['pos'], Type = 'Hetero', System = self.System.System,
                                              Species=self.Species, Elements=self.result_cache['syms'], Frame = i)
                    except Exception as e:
                        with open(self.Base + 'Sapphire_Errors.log', 'a') as f:
                            f.write('\nException raised while computing Hetero RDF was: \n%s' % e)


##############################################################################

        # This block evaluates all of the CoM calculations.

##############################################################################


        if 'com' in self.Quantities['Full']:
            self.result_cache['com'] = self.All_Atoms.get_center_of_mass() #CoM of WHOLE cluster
            try:
                if 'comdist' in self.Quantities['Full']:
                    CoMDist = DistFuncs.CoM_Dist(Positions=self.result_cache['pos'],
                                                       CoM=self.result_cache['com'],
                                                       Type = 'Full', Frame = i, System = self.System.System)

            except Exception as e:
                with open(self.Base + 'Sapphire_Errors.log', 'a') as f:
                    f.write('\nException raised while computing Full CoM Distances was: \n%s' % e)

        try:
            if self.System.System['Homo'] and 'hocom' in self.Quantities['Homo']:
                for x in self.System.System['Homo']:
                    self.New_Temp = DistFuncs.CoM_Dist(Positions=self.result_cache['pos'], System = self.System.System,
                                                       CoM = self.result_cache['com'], Type = 'Homo',
                                                       Specie = x, Elements = self.result_cache['syms'], Frame = i) #W.r.t sub-system
        except Exception as e:
            with open(self.Base + 'Sapphire_Errors.log', 'a') as f:
                f.write(
                    '\nException raised while computing in the Homo CoM Distances block was: \n%s' % e)


##############################################################################

        # This block creates histograms for the global pair-distances.
        # Note that this is a coarser calculation than the pdf.

##############################################################################

            if 'pair_distance' in self.Quantities['Full']:
                try:
                    PD = DistFuncs.Pair_Dist(Positions = self.result_cache['pos'], 
                                             Type = 'Full', Frame = i, System = self.System.System)
                except Exception as e:
                    with open(self.Base + 'Sapphire_Errors.log', 'a') as f:
                        f.write('\nException raised while computing pair distances: \n%s' % e)


        if 'hopair_distance' in self.Quantities['Homo']:
            for x in self.System.System['Homo']:
                try:
                    HoPD = DistFuncs.Pair_Dist(System = self.System.System,
                        Positions = self.result_cache['pos'], Type = 'Homo', 
                        Specie=x, Elements=self.result_cache['syms'], Frame = i)
                except Exception as e:
                    with open(self.Base + 'Sapphire_Errors.log', 'a') as f:
                        f.write(
                            '\nException raised while computing homo pair distances: \n%s' % e)

        try:
            if('he_pair_distance' in self.Quantities['Hetero']):
                PDHe = DistFuncs.Pair_Dist(System = self.System.System,
                    Positions = self.result_cache['pos'], Specie=self.Species, 
                    Type = 'Hetero', Elements=self.result_cache['syms'], Frame = i
                )
        except Exception as e:
            
            with open(self.Base + 'Sapphire_Errors.log', 'a') as f:
                f.write('\nException raised while computing hetero pair distances: \n%s' % e)

                    
##############################################################################

        # This block evaluates the adjacency matrices
        # for the whole system, homo pair(s), & hetero atoms

##############################################################################

        #The following block computes full-system adjacency properties
        if 'adj' in self.Quantities['Full']:
            try:
                self.result_cache['Adj'] = Adjacent.Adjacency_Matrix(
                    System = self.System.System,
                    Adj = 'adj' in self.Quantities['Full'], 
                    agcn = 'agcn' in self.Quantities['Full'], 
                    Surf_Area = 'surf_area' in self.Quantities['Full'], 
                    Surf_Atoms = 'surf_atoms' in self.Quantities['Full'], 
                    CN = 'nn' in self.Quantities['Full'],
                    Positions = self.result_cache['pos'],
                    Distances = self.result_cache['euc'],
                    R_Cut = self.result_cache['FullCut'],
                    Type = 'Full', Frame = i, Metals = self.Species, Elements = self.result_cache['syms']
                    ).ReturnAdj()

            except Exception as e:
                with open(self.Base + 'Sapphire_Errors.log', 'a') as f:
                    f.write('\nException raised while computing adjacency properties: \n%s' % e)


        #The following block computes homo-type adjacency properties
        if 'hoadj' in self.Quantities['Homo']:
            for x in self.System.System['Homo']: #Considering metals in the system
                try:
                    self.result_cache['HoAdj'+x] = Adjacent.Adjacency_Matrix(
                        System = self.System.System,
                        Adj = 'adj' in self.Quantities['Homo'], 
                        agcn = 'agcn' in self.Quantities['Homo'], 
                        Surf_Area = 'surf_area' in self.Quantities['Homo'], 
                        Surf_Atoms = 'surf_atoms' in self.Quantities['Homo'], 
                        CN = 'honn' in self.Quantities['Homo'],
                        Positions = self.result_cache['pos'],
                        Distances = self.result_cache['euc'],
                        R_Cut = self.result_cache['FullCut'],
                        Type = 'Homo', Frame = i, Metals = x, Elements = self.result_cache['syms']
                        ).ReturnAdj()
                except Exception as e:
                    with open(self.Base + 'Sapphire_Errors.log', 'a') as f:
                        f.write('\nException raised while computing HoAdj%s properties: \n%s' %(x,e))
                    
        #This next block computes adjacency properties for hetero-type interactions
        if 'headj' in self.Quantities['Hetero']:
            try:
                self.result_cache['HeAdj'] = Adjacent.Adjacency_Matrix(
                    System = self.System.System,
                    Adj = 'adj' in self.Quantities['Hetero'],
                    CN = 'henn' in self.Quantities['Hetero'],
                    Positions = self.result_cache['pos'],
                    Distances = self.result_cache['euc'],
                    R_Cut = self.result_cache['FullCut'],
                    Type = 'Hetero', Frame = i, Metals = self.Species, Elements = self.result_cache['syms']
                    ).ReturnAdj() 
            except Exception as e:
                with open(self.Base + 'Sapphire_Errors.log', 'a') as f:
                    f.write('\nException raised while computing HeAdj properties: \n%s' % e)
        # This next section computes the mixing parameter

##############################################################################

        # This block calculates the CNA signatures
        # and cna patterns for the whole system, only

##############################################################################

        if 'cna_sigs' in self.Quantities['Full']:
            try:
                cna = FrameSignature.CNA(
                    System=self.System.System, Frame=i,
                    adj=self.result_cache['Adj'], Type = 'Full',
                    Masterkey=self.Masterkey,
                    Fingerprint='cna_patterns' in self.Quantities['Full']).calculate()
                
            except Exception as e:
                with open(self.Base + 'Sapphire_Errors.log', 'a') as f:
                    f.write('\nException raised while computing CNA properties: \n%s' % e)

##############################################################################

        # This block calculates the radius of
        # gyration for the whole system and sub-systems

##############################################################################

        if 'gyration' in self.Quantities['Full']:
            try:
                Gyr = DistFuncs.Gyration(
                    System = self.System.System, Positions = self.result_cache['pos'], 
                    Type = 'Full', Metal = None, Elements = None, 
                    Masses=self.All_Atoms.get_masses(), Frame = i)
            except Exception as e:
                with open(self.Base + 'Sapphire_Errors.log', 'a') as f:
                    f.write('\nException raised while computing Gyration properties: \n%s' % e)


        if 'hogyration' in self.Quantities['Homo']:
            for Metal in self.System.System['Homo']:
                try:
                    HoGyr = DistFuncs.Gyration(
                        System = self.System.System, Positions = self.result_cache['pos'], 
                        Type = 'Homo', Metal = Metal, Elements = self.result_cache['syms'], Frame = i)
                except Exception as e:
                    with open(self.Base + 'Sapphire_Errors.log', 'a') as f:
                        f.write('\nException raised while computing HomoGyration%s properties: \n%s' %(Metal,e))


        if 'stat_radius' in self.Quantities['Full']:
            try:
                Stat_Rad = DistFuncs.Stat_Radius(self.System.System, self.result_cache['pos'], Frame = i)
            except Exception as e:
                with open(self.Base + 'Sapphire_Errors.log', 'a') as f:
                    f.write('\nException raised while computing Stat Radius properties: \n%s' % e)
        """

        #############################################################################################

        Robert:
            And now we check to see if the users wishes to evaluate any of the quantities
            from the energy file and add them to the metadata.

        """
        if self.System.System['energy_file_name']:
            self.energy = np.loadtxt(
                self.System.System['base_dir'] + self.System.System['energy_file_name'])
            try:
                if 'simtime' in self.Quantities['Energy']:
                    self.metadata['simtime'][int(
                        i/self.Step)] = self.energy[:, 0][int(i)]
            except Exception as e:
                with open(self.Base + 'Sapphire_Errors.log', 'a') as f:
                    f.write('\nException raised while computing simulation time: \n%s' % e)

            try:
                if 'epot' in self.Quantities['Energy']:
                    self.metadata['epot'][int(i/self.Step)
                                          ] = self.energy[:, 1][int(i)]
            except Exception as e:
                with open(self.Base + 'Sapphire_Errors.log', 'a') as f:
        
                    f.write('\nException raised while computing potential energy: \n%s' % e)

            try:
                if 'etot' in self.Quantities['Energy']:
                    self.metadata['etot'][int(i/self.Step)
                                          ] = self.energy[:, 2][int(i)]
            except Exception as e:
                with open(self.Base + 'Sapphire_Errors.log', 'a') as f:
                    f.write('\nException raised while computing total energy: \n%s' % e)

            try:
                if 'ekin' in self.Quantities['Energy']:
                    self.metadata['ekin'][int(i/self.Step)
                                          ] = self.energy[:, 3][int(i)]
            except Exception as e:
                with open(self.Base + 'Sapphire_Errors.log', 'a') as f:
                    f.write('\nException raised while computing kinetic energy: \n%s' % e)

            try:
                if 'edelta' in self.Quantities['Energy']:
                    self.metadata['edelta'][int(
                        i/self.Step)] = self.energy[:, 4][int(i)]
            except Exception as e:
                with open(self.Base + 'Sapphire_Errors.log', 'a') as f:
                    f.write('\nException raised while computing excess energy: \n%s' % e)

            try:
                if 'meanetot' in self.Quantities['Energy']:
                    self.metadata['meanetot'][int(
                        i/self.Step)] = self.energy[:, 5][int(i)]
            except Exception as e:
                with open(self.Base + 'Sapphire_Errors.log', 'a') as f:
                    f.write('\nException raised while computing mean total energy: \n%s' % e)

            try:
                if 'temp' in self.Quantities['Energy']:
                    self.metadata['temp'][int(i/self.Step)
                                          ] = self.energy[:, 6][int(i)]
            except Exception as e:
                with open(self.Base + 'Sapphire_Errors.log', 'a') as f:
                    f.write('\nException raised while computing temperature: \n%s' % e)

        try:
    
            with open(self.System.System['base_dir']+'Sapphire_Info.txt', "a") as f:
                f.write("\nAnalysis of frame %s required %.3f seconds.\n" %
                        (i, time.time() - self.timer))
                f.write("\nThis is approximately %.3fms for each atom.\n" % (
                    1000*(time.time() - self.timer)/self.NAtoms[int(i/self.Step)]))
        except Exception as e:
            print(e)
            
            
    def run_core(self):
        """
        Robert:

        """
        
        for i in self.All_Times:
            self.calculate(i)
        with open(self.System.System['base_dir']+'Sapphire_Info.txt', "a") as f:
            self.T3 = time.time()
            f.write('Time for completion is %s.\n' %
                    (time.strftime("%H:%M:%S", time.gmtime((self.T3-self.T)))))
            

    def analyse(self, Stat_Tools=None):
        """
        Robert:

            In general, the user will define which functions they wish to call on each distribution in the input file.

            This will then create a function object for each tool to be used and enters it into a tuple with the metadata keys
            for quantities to be analysed by that function.

            E.g., Jenson Shannon Divergence on the Radial Distribution Function.

            New metadata keys are then created for these analysed distributions and their frame-wise analysis values are stored
            under these keys.

            Moreover, for H / C statistics, these have T-1 and T-2 entries respectively and so new storage arrays for them
            must be instantiated separately while the latter is dependent on the former.


            """

        for i in range(1, int((self.End - self.Start)/self.Step)):

            # This  block calculates the concertedness and collectivity of atom rearrangements
            try:
                if self.Quantities['Full']['collect']:
                    self.result_cache['r'] = Adjacent.R(
                        self.metadata['adj'][i], self.metadata['adj'][i-1])
                    self.metadata['collect'][i -
                                             1] = Adjacent.Collectivity(self.result_cache['r'])
                    if not(i < 3):
                        if self.Quantities['Full']['concert']:
                            self.metadata['concert'][i-2] = Adjacent.Concertedness(self.metadata['collect'][i-1],
                                                                                   self.metadata['collect'][i-3])
            except Exception as e:
                with open(self.Base + 'Sapphire_Errors.log', 'a') as f:
                    f.write('\nException raised while computing collecivity and concertednes:\n%s' % e)

        try:
            from CNA.Utilities import Pattern_Key as PK
            self.pattern_key = list(PK().Key().keys())

            if 'cna_patterns' in self.Quantities['Full']:
                self.metadata['pattern_indices'] = np.empty(
                    (len(self.metadata['cna_patterns']),), dtype=object)
                for j, t in enumerate(self.metadata['cna_patterns']):
                    Temp_List = np.zeros(len(t))
                    for i, atom in enumerate(t):
                        if not str(atom) in self.pattern_key:
                            self.pattern_key.append(str(atom))
                        Temp_List[i] += self.pattern_key.index(str(atom))
                    self.metadata['pattern_indices'][j] = Temp_List
                with open(self.System.System['base_dir'] + 'AllPatterns.txt', 'w') as outfile:
                    for i, thing in enumerate(self.pattern_key):
                        outfile.write(str(i) + ')\t' + str(thing)+'\n')

        except Exception as e:
            with open(self.Base + 'Sapphire_Errors.log', 'a') as f:
                f.write('\nException raised while cleaning cna patterns: \n%s' % e)
        """
        This next block creates a dictionary whose keys are the analysis tools to be implemented.
        The first entry is the function to be called.
        All of the subsequent entries are the keys, as they appear in the metadata, to be processed.
        """
        if Stat_Tools is not None:
            self.Stat_Tools = Stat_Tools
            self.functions_list = [o for o in getmembers(
                Stats.Dist_Stats) if isfunction(o[1])]
            self.Stat_Keys = self.Stat_Tools.keys()
            self.Meta_Keys = self.metadata.keys()
            self.Calc_Dict = {}
            for obj in self.Stat_Keys:
                for item in self.functions_list:
                    if obj.lower() in item[0].lower():
                        self.Calc_Dict[item[0]] = [item[1]]

            for A_Key in self.Stat_Keys:
                for M_Key in self.Meta_Keys:
                    for obj in self.Stat_Tools[A_Key]:
                        if obj.lower() in M_Key.lower():
                            if M_Key.lower() == 'pdftype':
                                pass
                            else:
                                self.Calc_Dict[A_Key].append(M_Key)
                try:
                    self.Calc_Dict[A_Key].remove('pdftype')
                except ValueError:
                    pass
            """
            This next block reads over the previously created dictionary and then doctors the relevant
            metadata entry to be ready for processing.

            That is to say, that the heights of the distributions are to be analysed as the x-axis are
            largely uniform across the sample.
            """

            for A_Key in self.Stat_Keys:
                for obj in self.Calc_Dict[A_Key][1:]:
                    try:
                        self.metadata[A_Key +
                                      obj] = np.empty((len(self.metadata[obj]),), dtype=object)
                        # This is the initial distribution to which we shall make comparrisons
                        Init = self.metadata[obj][0][1]
                        for frame in range(len(self.metadata[obj])):
                            try:
                                # This is the y-axis of the distribution under consideration
                                Temp = self.metadata[obj][frame][1]
                                self.metadata[A_Key +
                                              obj][frame] = self.Calc_Dict[A_Key][0](Init, Temp)
                            except TypeError:
                                continue
                    except TypeError:
                        with open(self.Base + 'Sapphire_Errors.log', 'a') as f:
                            f.write("Type error raised by %s when performing statistical analysis of %s."
                                    % (obj, A_Key))
            del(self.result_cache)
        from Utilities import Output
        Out = Output.Writer(self.System.System, self.metadata)
        Out.Run('Full')
        if not self.System.System['Homo'] is None:
            Out.Run('Homo')
        if not self.System.System['Hetero'] is None:
            Out.Run('Hetero')

        if self.System.System['extend_xyz'] is not None:
            from Utilities import ExtendXYZ
            Write = ExtendXYZ.Extend(
                Traj=self.Dataset,
                System=self.System.System,
                Metadata=self.metadata,
                Names=self.System.System['extend_xyz']
            )
