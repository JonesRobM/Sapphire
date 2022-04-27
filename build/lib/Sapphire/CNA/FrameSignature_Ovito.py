# Import OVITO modules.
from ovito.io import *
from ovito.modifiers import *
from ovito.data import *

from ovito.io import import_file
from ovito.modifiers import CreateBondsModifier, CommonNeighborAnalysisModifier
from ovito.data import particles, BondsEnumerator

# Import standard Python and NumPy modules.
import numpy
import datetime
import platform
import getpass
import os

from Sapphire.CNA import Utilities

def cna_init(System = None, Pattern_Input = False):
    System = System
    if Pattern_Input:
        Pattern_Input = Pattern_Input
    filename = System['base_dir'] + System['movie_file_name']
    if (os.path.exists(System['base_dir']+'CNA_npz/')):
        pass
    else:
        os.mkdir(System['base_dir']+'CNA_npz/')
    npz_dir =  System['base_dir']+'CNA_npz/'
    __version__ = '1.0.0'
    Units = 'Angstrom & ev'
    if Pattern_Input:
        if not os.path.isfile(System['base_dir'] + 'CNA_Pattern_Info.txt'):
            with open(System['base_dir'] + 'CNA_Pattern_Info.txt', 'w') as f:
                f.write("""

              _____         _____  _____  _    _ _____ _____  ______
             / ____|  /\   |  __ \|  __ \| |  | |_   _|  __ \|  ____|     ____
            | (___   /  \  | |__) | |__) | |__| | | | | |__) | |__       /\__/\
             \___ \ / /\ \ |  ___/|  ___/|  __  | | | |  _  /|  __|     /_/  \_\
             ____) / ____ \| |    | |    | |  | |_| |_| | \ \| |____    \ \__/ /
            |_____/_/    \_\_|    |_|    |_|  |_|_____|_|  \_\______|    \/__\/
                        """
                                    "\nRunning version  -- %s --\n"
                                    "\nCurrent user is [ %s ]\n"
                                    "\nCalculation beginning %s\n"
                                    "\nArchitecture : [ %s ]\n"
                                    "\nUnits : [ %s ]\n"
                                    "\nThis file contains all of the user information regarding\nthe "
                                    "CNA Pattern Recognition and the support vector model.\n"
                                    %(__version__, getpass.getuser(), datetime.datetime.now().strftime("%a %d %b %Y %H:%M:%S"),
                                      platform.machine(), Units)
                    )
def row_histogram(a):
    ca = numpy.ascontiguousarray(a).view([('', a.dtype)] * a.shape[1])
    unique, indices, inverse = numpy.unique(ca, return_index=True, return_inverse=True)
    counts = numpy.bincount(inverse)
    return (a[indices], counts)
class Frame_CNA_Sigs():
    
    """
    Robert:
        
        THIS FUNCTION HAS SINCE BEEN DEPRACATED!
        PLEASE INSTEAD USE fRAMEsIGNATURE
        
        This function takes the following input arguments:
            
            frame: int - frame index of a movie.xyz trajectory.
            
            R_Cut: float - Interatomic separation. Can be set manually or read in 
            by a higher function. It is advised that this is set for each metal individually
            to ensure the most accurate and physically meaningful results.
            Will default to a value for the most abundant metal in a system otherwise - To be implemented -
            
            Masterkey: Tuple of tuples - Tells the programme which signatures are to be expected.
            There will be a mode in which this may be appended to and saved in an external file for reference.
            
            filename: str - self explanatory. where to find the file to be analysed.
            
            Homo: Boolean - Whether or not to search for a single atomic specie. Default of false means that the atoms will
            not be doctored prior to analysis.
            
            Metal: str - The name of the specie to be removed.
            
        Will output the following:
            
            signature_cna: a*N array where a is the length of the masterkey & N is the number of atoms considered.
            
            cna_pattern_indices: N dimensional list which tells you which pattern was found for each atom
            
            signature_cna_count: gives the frequency of a given pattern.
            
        May include a mode to write the patterns to a new file.
        
    """
    def __init__(self, System, frame, R_Cut=None, Masterkey=None, 
                   Type = None, Metal = None, Species = None,
                   Patterns = False, Fingerprint = None):
        self.System = System
        self.Frame = frame
        self.R_Cut = R_Cut
        self.Masterkey = Masterkey
        self.Type = Type
        self.Metal = Metal
        self.Patterns = Patterns
        self.Fingerprint = Fingerprint
        self.Masterkey = ((0,0,0),
                    (1,0,0),
                    (2,0,0),(2,1,1),
                    (3,0,0),(3,1,1),(3,2,2),
                    (4,0,0),(4,1,1),(4,2,1),(4,2,2),(4,3,3),(4,4,4),
                    (5,2,1),(5,2,2),(5,3,2),(5,3,3),(5,4,4),(5,5,5),
                    (6,6,6))
        self.Default = True
        self.npz_dir =  self.System['base_dir']+'CNA_npz/'
        self.filename = self.System['base_dir'] + self.System['movie_file_name']
        self.Pat_Key = Utilities.Pattern_Key().Key() #Calling dictionary of recognised patterns
        self.Keys = list(self.Pat_Key)
        self.Max_Label = max(len(str(label)) for label in self.Keys)
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
        
    def Ascii_Bars(self, Finger):
        with open(self.System['base_dir'] + 'CNA_Pattern_Info.txt', 'a', encoding='utf-8') as f:
            f.write('\nCNA Pattern distribution for full system at frame %s.\n' %self.Frame)
        Temp = numpy.zeros(len(self.Keys), int)
        for atom in Finger:
            if str(atom) in self.Keys:
                Temp[self.Keys.index(str(atom))] += 1
        with open(self.System['base_dir'] + 'CNA_Pattern_Info.txt', 'a', encoding='utf-8') as f:
            for i, count in enumerate(Temp):
                bar_chunks, remainder = divmod(int(count * 8 / (len(Finger) / 50)), 8)
                # First draw the full width chunks
                bar = '█' * bar_chunks
                # Then add the fractional part.  The Unicode code points for
                # block elements are (8/8), (7/8), (6/8), ... , so we need to
                # work backwards.
                if remainder > 0:
                    bar += chr(ord('█')+ (8 - remainder))
                # If the bar is empty, add a left one-eighth block
                bar = bar or  '|'
                f.write(f'{self.Keys[i].rjust(self.Max_Label)} | {count:#4d} {bar}\n')
    def calculate(self):
        """
        if self.Homo:
                
        
            self.pipeline.modifiers.append( SelectTypeModifier(property='Particle Type', 
                                                           types={self.Metal}) )
            
            self.pipeline.modifiers.append( DeleteSelectedModifier() )
        """
        pipeline = import_file(self.filename)

        pipeline.modifiers.append(CreateBondsModifier(cutoff = self.R_Cut))

        pipeline.modifiers.append( CommonNeighborAnalysisModifier(
            mode = CommonNeighborAnalysisModifier.Mode.BondBased))

        data = pipeline.compute(self.Frame)
        # The 'CNA Indices' bond property is a a two-dimensional array
        # containing the three CNA indices computed for each bond in the system.

        cna_indices = data.particles.bonds['CNA Indices']   

        # Used below for enumerating the bonds of each particle:
        bond_enumerator = BondsEnumerator(data.particles.bonds)

        particle_cnas = numpy.zeros((data.particles.count, len(self.Masterkey)), dtype=int)
        for particle_index in range(data.particles.count):
    
            # Create local list with CNA indices of the bonds of the current particle.
            bond_index_list = list(bond_enumerator.bonds_of_particle(particle_index))
            local_cna_indices = cna_indices[bond_index_list]
    
            # Count how often each type of CNA triplet occurred.
            unique_triplets, triplet_counts = row_histogram(local_cna_indices)

            for triplet, count in zip(unique_triplets, triplet_counts):
                try:
                    for index, signature in enumerate(self.Masterkey):
                        if(signature == tuple(triplet)):
                            particle_cnas[particle_index,index] = count
                                
                except KeyError:
                    pass
             
        Finger = numpy.zeros(len(particle_cnas), dtype=numpy.ndarray)
        for atom in range(len(particle_cnas)):
            Temp = []
            for i, x in enumerate(particle_cnas[atom]):
                if x > 0:
                    Temp.append((x, self.Masterkey[i]))
            
            Finger[atom] = Temp 
        
        Sigs = [ sum(particle_cnas[:,x]) for x in range(len(self.Masterkey)) ]
                    
        
        if (self.Fingerprint and self.Patterns):
            if self.Type == 'Homo':
                self.Ascii_Bars(Finger, True)
            else:
                self.Ascii_Bars(Finger)    
                
        self.Sigs = Sigs
        if self.Patterns:
            self.Finger = Finger
            
    def write(self):
        
        if self.Type == 'Full':
            from Sapphire.Utilities import OutputInfoFull as Out  # Case 1
            
            #Write object for the CoM
            Attributes = getattr(Out, str('cna_sigs')) #Loads in the write information for the object 
            OutFile = self.System['base_dir'] + Attributes['Dir'] + Attributes['File']
            self.ensure_dir(base_dir=self.System['base_dir'], file_path=Attributes['Dir'])   
            self.MakeFile(Attributes)
            with open(OutFile, 'a') as outfile:
                outfile.write(str(self.Frame) + ' ' +  ' '.join(str(item) for item in self.Sigs) +'\n')
                
            if self.Patterns:
   
                #Write object for the homo CoM distances
                Attributes = getattr(Out, str('pattern_indices')) #Loads in the write information for the object                  
                OutFile = self.System['base_dir'] + Attributes['Dir'] + Attributes['File']
                self.ensure_dir(base_dir=self.System['base_dir'], file_path=Attributes['Dir'])   
                self.MakeFile(Attributes)
                with open(OutFile, 'a') as outfile:
                    outfile.write(str(self.Frame) + ' ' +  ' '.join(str(item) for item in self.Finger) +'\n') 