#Pre-designed pythonic libraries. Indepedent of Sapphire
import numpy as np
import time
from inspect import getmembers, isfunction
from ase.io import read
import os

#Objects written specificaly for processing data
from Sapphire.Post_Process import Adjacent, Kernels, DistFuncs, Stats, Radii, AtomicEnvironment

#Common Neighbour Analysis specific functions
from Sapphire.CNA import FrameSignature, Utilities

#General purpose utility functions for parsing and tidying
from Sapphire.Utilities import Initial, System_Clean, Pattern_Clean
from Sapphire.Utilities import errors
from Sapphire.Utilities.log import get_logger

log = get_logger('Sapphire.Process')

class Process(object):

    """
    Parameters.

    ----------
   
    System : python dictionary
        Dictionary of basic analysis parameters provided by both the
        user and the System_Clean module.
       
    Quantities : python dictionary
        Dictionary of basic analysis parameters provided by both the
        user and the System_Clean module.
       
    Pattern_Input : python dictionary
        This is the dictionary style outpute of Sapphire which contains
        all of the relevant analysed information for facile reading/writing.

    Returns
    -------
    Performs the full analysis of a given structure / trajectory given the input
    information fed into the Quantities and System arguments.
    Each of the sub-modules may be indivudually interrogated.
   
    An example of a given input scheme may be found in examples/run_analysis.py
   
    System = {
        'base_dir': '/path/to/directory/',
        'movie_file_name': 'path/from/directory/to/movie_file_name.xyz',
        'extend_xyz': ['', '', ''],

        'Homo': ['Element1', 'Element2'],

        'Hetero': True,

        'Start': 0, 'End': None, 'Step': 1, 'Skip': 50, 'UniformPDF': False, 'Band': 0.05
    }

    # Define the quantities you want calculating given the names
    # in the supporting documentation.

    Quantities = {
        'Full':
        {
            'euc': None, 'rdf': None, 'pos': None,  'comdist': None,
            'moi': None, 'adj': None, 'pdf': None, 'pair_distance': None,
            'agcn': {'Write_Movie': False},
            'nn': None, 'com': None, 'cna_sigs': None,
            'cna_patterns': {'Write_Movie': True},
            'gyration': None, 'stat_radius': None,
            'surf_area': None, 'surf_atoms': None
        },

            'Homo':
        {
            'hopdf': None, 'hordf': None,
            'hocom': None, 'hoadj': None,
            'hocomdist': None, 'homidcomdist': None,
            'euc': None, 'hocna_sigs': None,
            'hocna_patterns': None, 'hogyration': None,
            'hosurf_area': None, 'hosurf_atoms': None,
            'hopair_distance': None
        },

            'Hetero':
        {
            'hepdf': None, 'herdf': None,
            'headj': None, 'mix': None,
            'he_pair_distance': None
        }
    }

    CNA_Pattern_Settings = {
        'npz_dir': 'CNA_npz/',  # folder to save the npz files in
        'new_xyz_dir': 'CNA_XYZs/',
        'APPEND_DICTIONARY': False,
        'FROM_MEMORY': False,
        'BULK_MASTERKEY': True,
        'SAVING_XYZ': True,
        'PRINTING_PATTERNS': True
    }

    Data = Process.Process(System=System, Quantities=Quantities,
                           Pattern_Input=CNA_Pattern_Settings)

    """
    
    def __init__(self, System=None, Quantities=None,
                 Pattern_Input=None, strict=False, overwrite=True):
        
        self.tick = time.time()
        # strict=True re-raises any exception instead of only logging it to Sapphire_Errors.log.
        self.strict = strict
        errors.set_strict(strict)
        # Output files are appended to frame by frame; a fresh run must not inherit an old one's.
        self.overwrite = overwrite

        self.System = System
        self.Quantities = Quantities
        
        #Check for user CNAP info - If not present - allow the cleanup 
        #module to sanitise and pass apropriate input in the Initialising step.
        if Pattern_Input:
            self.Pattern_Input = Pattern_Input
        else:
            self.Pattern_Input = None
            
        self.filename = System['base_dir']+System['movie_file_name']
        Initial.Logo(self.System['base_dir'])._write_()
        Initial.Info(self.System['base_dir'])._write_()

        self.Tbar = False


        self.result_cache = {}
        # Legacy in-memory store. Since the 1.0 refactor results are written to
        # <base_dir>/Time_Dependent/<Quantity>; this stays so analyse()/write_meta() degrade gracefully.
        self.metadata = {}
        """
        This list contains unique strings which are liable
        to require smaller storage objects in the metadata
        """
        self.T = time.time()

        """
        Initialise the CNA signature masterkey
        """

        self.Initialising()
        #self.run_pdf()
        self.run_core()

    def _report(self, exc, message):
        """Record a caught exception in Sapphire_Errors.log; re-raise when strict."""
        errors.report(exc, message, base_dir=self.Base, strict=self.strict)

    OUTPUT_DIRS = ('Time_Dependent/', 'Time_Dependent/Stats/', 'Adjacency/', 'CNA/', 'Exec/')

    def _clear_outputs(self):
        """Remove previous result files (file by file) from the Sapphire-owned output folders."""
        removed = 0
        for d in self.OUTPUT_DIRS:
            folder = self.Base + d
            if not os.path.isdir(folder):
                continue
            for name in os.listdir(folder):
                path = os.path.join(folder, name)
                if os.path.isfile(path):
                    os.remove(path); removed += 1
        for name in ('Sapphire_Errors.log', 'Extended.xyz', 'Metadata.pkl'):
            if os.path.isfile(self.Base + name):
                os.remove(self.Base + name); removed += 1
        if removed:
            with open(self.Base + 'Sapphire_Info.txt', 'a') as f:
                f.write('Cleared %d result files from a previous run (overwrite=True).\n' % removed)

    def _write_series(self, name, values, subdir='Time_Dependent/Stats/'):
        """Persist a per-frame series so Reader / Graphing can find it (frame value per line)."""
        self.ensure_dir(self.Base, subdir)
        with open(self.Base + subdir + name, 'w') as f:
            for frame, val in zip(self.All_Times, np.atleast_1d(values)):
                f.write('%d %s\n' % (frame, val))

    def reader(self):
        """A :class:`Sapphire.IO.Reader.Reader` over this run's output directory."""
        from Sapphire.IO.Reader import Reader
        return Reader(self.Base)

    def load_metadata(self, keys=None):
        """Populate ``self.metadata`` from the files written during ``run_core``.

        Results are written to disk frame by frame; this reads them back so the
        statistical tools, ``write_meta`` and the extended-xyz writer can use them.
        """
        r = self.reader()
        available = r.available()
        for k in (available if keys is None else [k for k in keys if k in available]):
            self.metadata[k] = r.load(k)
        self.metadata['masterkey'] = r.masterkey()
        return self.metadata

    def ensure_dir(self, base_dir='', file_path=''):
        
        """Create ``base_dir + file_path`` if it does not exist."""

        directory = base_dir + file_path
        if not os.path.exists(directory):

            os.makedirs(directory)

    def MakeFile(self, Attributes):
        
        """Create the (empty) output file described by an ``OutputInfo`` entry if absent."""
        
        self.out = self.System['base_dir'] + Attributes['Dir'] + Attributes['File']

        if not os.path.isfile(self.out):
            with open(self.System['base_dir'] + Attributes['Dir'] + Attributes['File'], 'w') as out:
                out.close()
        else:
            pass

##############################################################################

    def Initialising(self):


        """Sanitise the user input, load the trajectory and derive per-run constants.

        Runs the ``System_Clean`` / ``Pattern_Clean`` validators, reads every frame with ASE,
        records species, atom counts and the CNA masterkey, and clears previous output files
        when ``overwrite`` is set.
        """

        self.Calc_Quants = {}

        # Absent quantity groups mean "nothing requested", not a crash.

        self.Quantities = dict(self.Quantities or {})

        for group in ('Full', 'Homo', 'Hetero'):

            self.Quantities.setdefault(group, {})
        # This next sub-block instantiates the system object.
        # It parses the user input into a clean-up module to sanitise arguments.
        self.System = System_Clean._Clean_System(self.System).System
        
        self.Pattern_Input = Pattern_Clean._Clean_Pattern(
            self.System,self.Pattern_Input).Pattern_Input
        
        try:
            if self.Pattern_Input['BULK_MASTERKEY']:
                self.Masterkey = Utilities.CNA_Masterkey().Key()

            else:
                self.Masterkey = []
        except KeyError:
            self.Masterkey = Utilities.CNA_Masterkey().Key()

        self.Time = int(
            (self.System['End'] - self.System['Start']) / (self.System['Step']))
        self.Step = self.System['Step']
        self.Skip = self.System['Skip']
        self.End = self.System['End']
        self.Start = self.System['Start']
        self.Base = self.System['base_dir']
        if self.overwrite:
            self._clear_outputs()

##############################################################################

        with open(self.System['base_dir']+'Sapphire_Info.txt', "a") as f:
            f.write("Initialising system environment took %.3f seconds.\n" %
                    (time.time()-self.tick))

##############################################################################

        tick = time.time()

        """
        Robert:

            This next block loads up the first frame of the trajectory and sets some initial file parameters
            deciding how to treat poly-metallic or mono-metallic systems depending on the user input.
        """

##############################################################################

        with open(self.System['base_dir']+'Sapphire_Info.txt', "a") as f:
            f.write("Loading in the dataset to be analysed.\n")
            f.write("Be aware that this may take a while for a large file.\n")
        Read_Time = time.time()
        self.Dataset = read(self.filename, index=':')
        with open(self.System['base_dir']+'Sapphire_Info.txt', "a") as f:
            f.write("Opened the dataset in %.3f seconds.\n" %
                    (time.time()-Read_Time))
        self.all_positions = self.Dataset[0].get_positions()
        self.max_dist = max(DistFuncs.Euc_Dist(self.all_positions))
        del(self.all_positions)
        
        #self.all_atoms contains the chemical symbols for all frames
        self.all_atoms = [ 
            self.Dataset[t].get_chemical_symbols() for t in range(
                self.System['Start'], self.System['End'], self.System['Step']
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

        with open(self.System['base_dir']+'Sapphire_Info.txt', "a") as f:
            f.write("Checking user input for calculating homo properties in this run.\n")

        if self.System['Homo'] is None:
            with open(self.System['base_dir']+'Sapphire_Info.txt', "a") as f:
                f.write("No bimetallic properties will be calculated in this run.\n")

        else:
            with open(self.System['base_dir']+'Sapphire_Info.txt', "a") as f:
                f.write("Homo atom properties will be caluclated for %s in this run.\n" % (
                    self.System['Homo']))
      
            with open(self.System['base_dir']+'Sapphire_Info.txt', "a") as f:
                f.write("Checking user input for hetero atomic species.\n")

##############################################################################

            with open(self.System['base_dir']+'Sapphire_Info.txt', "a") as f:
                f.write("Initialising Metadata took %.3f seconds.\n" %
                        (time.time() - tick))

##############################################################################

        # Needed for mono- and bi-metallic runs alike (was only set in the bimetallic branch).
        self.All_Times = list(range(self.Start, self.End, self.Step))
        self.Band = self.System['Band']

##############################################################################

    def calculate(self, i):
        """Analyse frame ``i``: compute every requested quantity and write it to disk.

        Quantities are evaluated in dependency order (pair distances → PDDF and cutoff →
        adjacency → coordination / aGCN / CNA → chemical ordering → radii and shape). Each
        block is guarded: a failure is recorded via :meth:`_report` and the remaining blocks
        still run, unless ``strict=True``.
        """

        self.timer = time.time()

        with open(self.System['base_dir']+'Sapphire_Info.txt', "a") as f:
            f.write("\nLoading in atoms for frame %s.\n" % i)
        self.All_Atoms = self.Dataset[i]
        with open(self.System['base_dir']+'Sapphire_Info.txt', "a") as f:
            f.write("Loaded the atoms in %.3f seconds.\n" % (time.time()-self.timer))

        self.result_cache['pos'] = self.All_Atoms.get_positions()
        self.result_cache['euc'] = DistFuncs.Euc_Dist(self.result_cache['pos'])
        self.result_cache['syms'] = self.All_Atoms.get_chemical_symbols()
        if 'pdf' in self.Quantities['Full']:
            
            try:
                self.result_cache['FullCut'] = Kernels.Gauss(Data = self.result_cache['euc'], Band = self.Band, 
                                                                Ele = None, Type = 'Full', Space = None, 
                                                                System = self.System, Frame = i).ReturnRCut()
            except Exception as e:
                self._report(e, '\nException raised while computing Full PDF: \n%s')
        if self.System['Hetero']:
            if 'hepdf' in self.Quantities['Hetero']:
                
                TempPos = DistFuncs.Hetero(self.result_cache['pos'], self.Species,
                                                                  self.result_cache['syms'])
                
                self.result_cache['heteropos'] = np.asarray(TempPos, dtype=float).ravel()     
                
                self.result_cache['HeCut'] = Kernels.Gauss(Data = self.result_cache['heteropos'], Band = self.Band, 
                                                          Ele = None, Type='Hetero', Space = None, 
                                                          System = self.System, Frame = i).ReturnRCut()

                

        if 'hopdf' in self.Quantities['Homo']:
            for x in self.System['Homo']:
                self.result_cache['homoed'+x] = DistFuncs.Euc_Dist(positions=self.result_cache['pos'], homo=True, specie=x, elements=self.result_cache['syms'])
                if self.result_cache['homoed'+x] is not None:
                    try:
                        self.result_cache[x+'Cut'] = Kernels.Gauss(self.result_cache['homoed'+x], 
                                              Band = self.Band, Ele = x, 
                                              Type='Homo', Space = None, 
                                              System = self.System, Frame = i).ReturnRCut()
                    except Exception as e:
                        self._report(e, '\nException raised while computing Homo rcut: \n%s')


        if 'moi' in self.Quantities['Full']:
            self.moi = self.All_Atoms.get_moments_of_inertia()
            from Sapphire.IO import OutputInfoFull as Out
            Attributes = getattr(Out, str('moi')) #Loads in the write information for the object 
            OutFile = self.System['base_dir'] + Attributes['Dir'] + Attributes['File']
            self.ensure_dir(base_dir=self.System['base_dir'], file_path=Attributes['Dir'])   
            self.MakeFile(Attributes)
            with open(OutFile, 'a') as outfile:
                outfile.write(str(i) + ' ' +  ' '.join(str(item) for item in self.moi) +'\n')        


##############################################################################

        # All RDF calculations performed in the following block

##############################################################################


        if 'rdf' in self.Quantities['Full']:
            try:
                DistFuncs.RDF(System = self.System,
                    Positions = self.result_cache['pos'], Type = 'Full', Frame = i)
                
            except Exception as e:
                self._report(e, '\nException raised while computing Full RDF was: \n%s')
                    
            if self.System['Homo'] and 'hordf' in self.Quantities['Homo']:
                try:
                    for x in self.System['Homo']:
                        self.result_cache['homopos'+x] = DistFuncs.get_subspecieslist(x, self.result_cache['syms'], self.result_cache['pos'])
                        DistFuncs.RDF(System = self.System, Positions = self.result_cache['homopos'+x],
                                      Type = 'Homo', Species = x, Frame = i)
                        
                except Exception as e:
                    self._report(e, '\nException raised while computing Homo RDF: \n%s')
                        
            if self.System['Hetero']:
                if self.System['Hetero'] and 'herdf' in self.Quantities['Hetero']:
                    try:
                        DistFuncs.RDF(System = self.System, Positions = self.result_cache['pos'], Type = 'Hetero',
                                              Species=self.Species, Elements=self.result_cache['syms'], Frame = i)
                    except Exception as e:
                        self._report(e, '\nException raised while computing Hetero RDF: \n%s')


##############################################################################

        # This block evaluates all of the CoM calculations.

##############################################################################


        if 'com' in self.Quantities['Full']:
            self.result_cache['com'] = self.All_Atoms.get_center_of_mass() #CoM of WHOLE cluster
            try:
                if 'comdist' in self.Quantities['Full']:
                    DistFuncs.CoM_Dist(Positions=self.result_cache['pos'],
                                                       CoM=self.result_cache['com'],
                                                       Type = 'Full', Frame = i, System = self.System)

            except Exception as e:
                self._report(e, '\nException raised while computing Full CoM Distances was: \n%s')

        try:
            if self.System['Homo'] and 'hocom' in self.Quantities['Homo']:
                for x in self.System['Homo']:
                    self.New_Temp = DistFuncs.CoM_Dist(Positions=self.result_cache['pos'], System = self.System,
                                                       CoM = self.result_cache['com'], Type = 'Homo',
                                                       Specie = x, Elements = self.result_cache['syms'], Frame = i) #W.r.t sub-system
        except Exception as e:
            self._report(e, '\nException raised while computing in the Homo CoM Distances block was: \n%s')


##############################################################################

        # This block creates histograms for the global pair-distances.
        # Note that this is a coarser calculation than the pdf.

##############################################################################

            if 'pair_distance' in self.Quantities['Full']:
                try:
                    DistFuncs.Pair_Dist(Positions = self.result_cache['pos'], 
                                             Type = 'Full', Frame = i, System = self.System)
                except Exception as e:
                    self._report(e, '\nException raised while computing pair distances: \n%s')


        if 'hopair_distance' in self.Quantities['Homo']:
            for x in self.System['Homo']:
                try:
                    DistFuncs.Pair_Dist(System = self.System,
                        Positions = self.result_cache['pos'], Type = 'Homo', 
                        Specie=x, Elements=self.result_cache['syms'], Frame = i)
                except Exception as e:
                    self._report(e, '\nException raised while computing homo pair distances: \n%s')

        try:
            if self.System['Hetero']:
                if('he_pair_distance' in self.Quantities['Hetero']):
                    DistFuncs.Pair_Dist(System = self.System,
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
                    System = self.System,
                    Adj = 'adj' in self.Quantities['Full'], 
                    agcn = 'agcn' in self.Quantities['Full'], 
                    Surf_Area = 'surf_area' in self.Quantities['Full'], 
                    Surf_Atoms = 'surf_atoms' in self.Quantities['Full'], 
                    CN = 'nn' in self.Quantities['Full'],
                    Positions = self.result_cache['pos'],
                    Distances = self.result_cache['euc'],
                    R_Cut = self.result_cache['FullCut'],
                    Type = 'Full', Frame = i, 
                    Metals = self.Species, 
                    Elements = self.result_cache['syms']
                    ).ReturnAdj()

            except Exception as e:
                self._report(e, '\nException raised while computing adjacency properties: \n%s')


        #The following block computes homo-type adjacency properties
        if 'hoadj' in self.Quantities['Homo']:
            for x in self.System['Homo']: #Considering metals in the system
                try:
                    self.result_cache['HoAdj'+x] = Adjacent.Adjacency_Matrix(
                        System = self.System,
                        Adj = 'adj' in self.Quantities['Homo'], 
                        agcn = 'agcn' in self.Quantities['Homo'], 
                        Surf_Area = 'surf_area' in self.Quantities['Homo'], 
                        Surf_Atoms = 'surf_atoms' in self.Quantities['Homo'], 
                        CN = 'honn' in self.Quantities['Homo'],
                        Positions = self.result_cache['pos'],
                        Distances = self.result_cache['euc'],
                        R_Cut = self.result_cache['FullCut'],
                        Type = 'Homo', Frame = i, 
                        Metals = [x], 
                        Elements = self.result_cache['syms']
                        ).ReturnAdj()
                except Exception as e:
                    self._report(e, '\nException raised while computing HoAdj%s properties: \n%s' %(x,e))

        #This next block computes adjacency properties for hetero-type interactions
        if self.System['Hetero']:
            if 'headj' in self.Quantities['Hetero']:
                try:
                    self.result_cache['HeAdj'] = Adjacent.Adjacency_Matrix(
                        System = self.System,
                        Adj = 'adj' in self.Quantities['Hetero'],
                        CN = 'henn' in self.Quantities['Hetero'],
                        Positions = self.result_cache['pos'],
                        Distances = self.result_cache['euc'],
                        R_Cut = self.result_cache['FullCut'],
                        Type = 'Hetero', Frame = i,
                        Metals = self.Species,
                        Elements = self.result_cache['syms']
                        ).ReturnAdj()
                except Exception as e:
                    self._report(e, '\nException raised while computing HeAdj properties: \n%s')
        # This next section computes the mixing parameter

##############################################################################

        # This block calculates the CNA signatures
        # and cna patterns for the whole system, only

##############################################################################

        if 'cna_sigs' in self.Quantities['Full']:
            try:
                cna = FrameSignature.CNA(self.System, self.result_cache['Adj'], self.Masterkey,
                                         'cna_patterns' in self.Quantities['Full'], Type = 'Full', Frame = i)
                cna.calculate()
                # carry every signature seen so far forward: later frames' columns then extend, never reorder
                self.Masterkey = tuple(cna.Sigs.keys())
                
            except Exception as e:
                self._report(e, '\nException raised while computing CNA properties: \n%s')

##############################################################################

        # This block calculates the radius of
        # gyration for the whole system and sub-systems

##############################################################################

        if 'gyration' in self.Quantities['Full']:
            try:
                Radii.Gyration(
                    System = self.System, Positions = self.result_cache['pos'], 
                    Type = 'Full', Metal = None, Elements = None, 
                    Masses=self.All_Atoms.get_masses(), Frame = i)
            except Exception as e:
                self._report(e, '\nException raised while computing Gyration properties: \n%s')


        if 'hogyration' in self.Quantities['Homo']:
            for Metal in self.System['Homo']:
                try:
                    Radii.Gyration(
                        System = self.System, Positions = self.result_cache['pos'], 
                        Type = 'Homo', Metal = Metal, Elements = self.result_cache['syms'], Frame = i)
                except Exception as e:
                    self._report(e, '\nException raised while computing HomoGyration%s properties: \n%s' %(Metal,e))

        if 'stat_radius' in self.Quantities['Full']:
            try:
                Radii.Stat_Radius(self.System, self.result_cache['pos'], Frame = i)
            except Exception as e:
                self._report(e, '\nException raised while computing Stat Radius properties: \n%s')

##############################################################################

        # This block evaluates local atomic environment variables such as the
        # mixing parameter, LAE, and homo / hetero "bond" analyses.

##############################################################################
        if self.System['Hetero'] and self.NSpecies == 2:
            HeQ, HoQ = self.Quantities['Hetero'], self.Quantities['Homo']
            hoadj = {x: self.result_cache[k] for x in self.System['Homo'] if (k := 'HoAdj' + x) in self.result_cache}
            if 'mix' in HeQ or 'heterobonds' in HeQ or 'homobonds' in HoQ:
                try:
                    AtomicEnvironment.Mix(System=self.System, Frame=i, HoAdj=hoadj,
                                          HeAdj=self.result_cache['HeAdj'],
                                          HomoBonds='homobonds' in HoQ, HeteroBonds='heterobonds' in HeQ)
                except Exception as e:
                    self._report(e, '\nException raised while computing the mixing parameter / bond counts: \n%s')
            if 'lae' in HeQ:
                try:
                    AtomicEnvironment.LAE(System=self.System, Frame=i, HeAdj=self.result_cache['HeAdj'],
                                          Species=self.Species)
                except Exception as e:
                    self._report(e, '\nException raised while computing the LAE: \n%s')
            if 'ele_nn' in HeQ:
                try:
                    AtomicEnvironment.Ele_NN(System=self.System, Frame=i, Adj=self.result_cache['Adj'],
                                             Elements=self.result_cache['syms'], Species=self.Species)
                except Exception as e:
                    self._report(e, '\nException raised while computing per-species neighbour counts: \n%s')

##############################################################################

        # This block simply writes a summary of the ith step of the
        # post-processing.

##############################################################################

        try:
    
            with open(self.System['base_dir']+'Sapphire_Info.txt', "a") as f:
                f.write("\nAnalysis of frame %s required %.3f seconds.\n" %
                        (i, time.time() - self.timer))
                f.write("\nThis is approximately %.3fms for each atom.\n" % (
                    1000*(time.time() - self.timer)/self.NAtoms[int(i/self.Step)]))
        except Exception as e:
            log.info(e)
            
            
    def run_core(self):
        """Loop :meth:`calculate` over the requested frames, logging timing to ``Sapphire_Info.txt``."""
        
        for i in self.All_Times:
            self.calculate(i)
        with open(self.System['base_dir']+'Sapphire_Info.txt', "a") as f:
            self.T3 = time.time()
            f.write('Time for completion is %s.\n' %
                    (time.strftime("%H:%M:%S", time.gmtime((self.T3-self.T)))))
            

    def analyse(self, Stat_Tools=None):
        """ In general, the user will define which functions they wish to call on each distribution in the input file.

            This will then create a function object for each tool to be used and enters it into a tuple with the metadata keys
            for quantities to be analysed by that function.

            E.g., Jenson Shannon Divergence on the Radial Distribution Function.

            New metadata keys are then created for these analysed distributions and their frame-wise analysis values are stored
            under these keys.

            Moreover, for H / C statistics, these have T-1 and T-2 entries respectively and so new storage arrays for them
            must be instantiated separately while the latter is dependent on the former.
        
        :param startHnd: Start index, defaults to 1
        :type startHnd: int, optional
        """

        self.load_metadata()
        if 'adj' in self.metadata and ('collect' in self.Quantities['Full'] or 'concert' in self.Quantities['Full']):
            n = len(self.metadata['adj'])
            self.metadata.setdefault('collect', np.zeros(max(n - 1, 0)))   # one value per consecutive pair
            self.metadata.setdefault('concert', np.zeros(max(n - 2, 0)))   # needs two pairs
        # one comparison per consecutive pair of analysed frames
        for i in range(1, len(self.metadata.get('adj', []))):

            # This  block calculates the concertedness and collectivity of atom rearrangements
            try:
                if 'collect' in self.Quantities['Full']:
                    self.result_cache['r'] = Stats.Mobility.R(
                        self.metadata['adj'][i], self.metadata['adj'][i-1])
                    self.metadata['collect'][i -
                                             1] = Stats.Mobility.Collectivity(self.result_cache['r'])
                    if not(i < 3):
                        if 'concert' in self.Quantities['Full']:
                            self.metadata['concert'][i-2] = Stats.Mobility.Concertedness(self.metadata['collect'][i-1],
                                                                                   self.metadata['collect'][i-3])
            except Exception as e:
                self._report(e, '\nException raised while computing collecivity and concertednes:\n%s')

        for key, attr in (('collect', 'collect'), ('concert', 'concert')):
            if key in self.metadata and key in self.Quantities['Full']:
                from Sapphire.IO import OutputInfoFull as Out
                A = getattr(Out, attr)
                self.ensure_dir(self.Base, A['Dir'])
                offset = len(self.All_Times) - len(self.metadata[key])   # series are shorter than the frame list
                with open(self.Base + A['Dir'] + A['File'], 'w') as f:
                    for frame, val in zip(self.All_Times[offset:], self.metadata[key]):
                        f.write('%d %s\n' % (frame, val))
        try:
            from Sapphire.CNA.Utilities import Pattern_Key as PK
            self.pattern_key = list(PK().Key().keys())

            if 'cna_patterns' in self.Quantities['Full'] and 'cna_patterns' in self.metadata:
                self.metadata['pattern_indices'] = np.empty(
                    (len(self.metadata['cna_patterns']),), dtype=object)
                for j, t in enumerate(self.metadata['cna_patterns']):
                    Temp_List = np.zeros(len(t))
                    for i, atom in enumerate(t):
                        if not str(atom) in self.pattern_key:
                            self.pattern_key.append(str(atom))
                        Temp_List[i] += self.pattern_key.index(str(atom))
                    self.metadata['pattern_indices'][j] = Temp_List
                with open(self.System['base_dir'] + 'AllPatterns.txt', 'w') as outfile:
                    for i, thing in enumerate(self.pattern_key):
                        outfile.write(str(i) + ')\t' + str(thing)+'\n')

        except Exception as e:
            self._report(e, '\nException raised while cleaning cna patterns: \n%s')
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
            self.Meta_Keys = list(self.metadata.keys())
            self.Calc_Dict = {}
            for obj in self.Stat_Keys:
                for item in self.functions_list:
                    if obj.lower() in item[0].lower():
                        self.Calc_Dict[item[0]] = [item[1]]

            for A_Key in self.Stat_Keys:
                for M_Key in self.Meta_Keys:
                    for obj in self.Stat_Tools[A_Key]:
                        # 'pdf' -> pdf, hopdfAu, hepdf ... but never the *space grids or pdftype
                        k = M_Key.lower()
                        derived = k.startswith(tuple(s.lower() for s in self.Stat_Keys))  # e.g. an earlier run's JSDpdf
                        if obj.lower() in k and 'space' not in k and k != 'pdftype' and not derived:
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
                func = self.Calc_Dict[A_Key][0]
                for obj in self.Calc_Dict[A_Key][1:]:
                    try:
                        dist = np.asarray(self.metadata[obj], dtype=float)
                        if dist.ndim != 2:
                            raise TypeError('%s is not a per-frame distribution' % obj)
                        Init = dist[0]  # reference distribution: the first frame
                        self.metadata[A_Key + obj] = np.array([func(Init, frame) for frame in dist])
                        self._write_series(A_Key + obj, self.metadata[A_Key + obj])
                    except Exception as e:
                        self._report(e, "\nError raised when performing %s analysis of %s: %%s" % (A_Key, obj))
            del(self.result_cache)


        if self.System['extend_xyz'] is not None:
            from Sapphire.Utilities import ExtendXYZ
            if set(self.metadata) <= {'masterkey'}:
                self.load_metadata()
            ExtendXYZ.Extend(
                Traj=self.Dataset,
                System=self.System,
                Metadata=self.metadata,
                Names=self.System['extend_xyz']
            )

    def write_meta(self, filename='Metadata.pkl'):
        """Pickle everything the run wrote (read back via :class:`Sapphire.IO.Reader.Reader`)."""
        import pickle
        if set(self.metadata) <= {'masterkey'}:
            self.load_metadata()
        path = self.Base + filename
        with open(path, 'wb') as f:
            pickle.dump(self.metadata, f, protocol=pickle.HIGHEST_PROTOCOL)
        return path
