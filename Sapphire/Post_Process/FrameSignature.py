# Import OVITO modules.
from ovito.io import *
from ovito.modifiers import *
from ovito.data import *

from ovito.io import import_file
from ovito.modifiers import CreateBondsModifier, CommonNeighborAnalysisModifier
from ovito.data import particles, BondsEnumerator

# Import standard Python and NumPy modules.
import sys
import numpy
import datetime
import platform
import getpass
import time

import glob
import os
import sys

from CNA import Utilities

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
    
    __version__ = '0.9.0'
    Units = 'Angstrom & ev'
    if Pattern_Input:
        if not os.path.isfile(System['base_dir'] + 'CNA_Pattern_Info.txt'):
            with open(System['base_dir'] + 'CNA_Pattern_Info.txt', 'w') as f:
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


def Frame_CNA_Sigs(System, frame, R_Cut=None, MasterKey=None, 
                   Homo = False, Metal = None, Patterns = False):
    
    def row_histogram(a):
        ca = numpy.ascontiguousarray(a).view([('', a.dtype)] * a.shape[1])
        unique, indices, inverse = numpy.unique(ca, return_index=True, return_inverse=True)
        counts = numpy.bincount(inverse)
        return (a[indices], counts)
    
    Masterkey = ((0,0,0),
            (1,0,0),
            (2,0,0),(2,1,1),
            (3,0,0),(3,1,1),(3,2,2),
            (4,0,0),(4,1,1),(4,2,1),(4,2,2),(4,3,3),(4,4,4),
            (5,2,1),(5,2,2),(5,3,2),(5,3,3),(5,4,4),(5,5,5),
            (6,6,6))
    
    npz_dir =  System['base_dir']+'CNA_npz/'
    filename = System['base_dir'] + System['movie_file_name']
    
    pipeline = import_file(filename)
    
    if MasterKey is None:
        MasterKey = Masterkey
    else:
        MasterKey = MasterKey
    if Patterns:
        Patterns = True
    else:
        Patterns = False

    
    #with open(System['base_dir'] + 'CNA_Pattern_Info.txt', 'a') as f:
     #   f.write("\nCurrently reading CNA Signatures and Patterns for frame [ %s ].\nR_Cut has been set to %s.\n" %(frame,R_Cut))
      #  if Metal is not None:
       #     f.write('Metal considered is [ %s ].\n'%Metal)
        #    f.write("Initialising system environment took %.3f seconds.\n" %(time.time()-tick))
    
    """
    Robert:
        
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
    tick = time.time()
    
    if bool((Homo is True)*(type(Metal) is str)):
            
    
        pipeline.modifiers.append( SelectTypeModifier(property='Particle Type', 
                                                       types={Metal}) )
        
        pipeline.modifiers.append( DeleteSelectedModifier() )

    pipeline.modifiers.append(CreateBondsModifier(cutoff = R_Cut))

    pipeline.modifiers.append( CommonNeighborAnalysisModifier(
        mode = CommonNeighborAnalysisModifier.Mode.BondBased))
    
    data = pipeline.compute(frame)
            
    
    # The 'CNA Indices' bond property is a a two-dimensional array
    # containing the three CNA indices computed for each bond in the system.
    
    cna_indices = data.particles.bonds['CNA Indices']   

    # Used below for enumerating the bonds of each particle:
    bond_enumerator = BondsEnumerator(data.particles.bonds)
    particle_cnas = numpy.zeros((data.particles.count, len(MasterKey)), dtype=int)
    # Loop over particles and print their CNA indices.
    for particle_index in range(data.particles.count):

        # Create local list with CNA indices of the bonds of the current particle.
        bond_index_list = list(bond_enumerator.bonds_of_particle(particle_index))
        local_cna_indices = cna_indices[bond_index_list]

        # Count how often each type of CNA triplet occurred.
        unique_triplets, triplet_counts = row_histogram(local_cna_indices)

        # Print list of triplets with their respective counts.
        for triplet, count in zip(unique_triplets, triplet_counts):
            try:
                for index, signature in enumerate(MasterKey):
                    if(signature == tuple(triplet)):
                        particle_cnas[particle_index,index] = count
                            
            except KeyError:
                pass

    if Patterns:
        
        signature_cna, CNA_pattern_indices, signature_cna_count = numpy.unique(
            particle_cnas, axis=0, return_inverse=True, return_counts=True)
        if bool((Homo is True)*(type(Metal) is str)):
            cna_path = npz_dir+'/CNA_'+System['movie_file_name'][:-4]+'-'+ Metal + '-' + str(frame)+'.npz'
            
        else:
            cna_path = npz_dir+'/CNA_'+System['movie_file_name'][:-4]+'-'+str(frame)+'.npz'
        
        
        numpy.savez(cna_path, 
                    particle_cnas = particle_cnas, 
                    signature_cna = signature_cna, 
                    signature_cna_count = signature_cna_count, 
                    signature_cna_indices = CNA_pattern_indices)
        

    print([ sum(particle_cnas[:,x]) for x in range(len(MasterKey)) ])
    return [ sum(particle_cnas[:,x]) for x in range(len(MasterKey)) ]