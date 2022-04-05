import pickle
import numpy
from ase.io import read

from FramePattern import patterns 
from FrameSignature import signature
import Utilities
import xyz_io


"""
Robert:
    
    This is a test script for running Armand's CNA Pattern code with the new 
    class structure.
    
    Generally speaking, the version of this that will be run in post-
    processing will borrow mmost of its input parameters from the PP input.
    
    However; everything will be explicitly written here.
"""

#for T in range(300,650,50):
#for sys in ['Sim-1345/', 'Sim-2783/', 'Sim-3987/', 'Sim-4009/']:
System = {
        'base_dir' : '/home/k1899676/Downloads/',
        'movie_file_name' : 'short_movie.xyz',
        'energy_file_name' : 'energy.out',
        #'Path_to_Pot' : '/path/to/potential/file/file.pot',
        'New_agcn_movie' : True,
        
        #'Homo' : ['Au', 'Pd'], 'HomoQuants' : [ 'HoPDF', 'HoRDF', 'CoM', 'HoAdj', 'CoMDist', 'MidCoMDist'], 
        #'Hetero' : True, 'HeteroQuants' : [ 'HePDF', 'HeRDF', 'HeAdj' ],
        
        'Start' : 0, 'End' : 300, 'Step' : 1, 'Skip' : 10, 'UniformPDF' : False, 'Band' : 0.05
        
        }

CNA_Pattern_Settings = {
                        'npz_dir' : 'CNA_npz/', #folder to save the npz files in
                        'new_xyz_dir' : 'CNA_XYZs/',
                        'APPEND_DICTIONARY' : True,
                        'FROM_MEMORY' : False,
                        'BULK_MASTERKEY' : True,
                        'SAVING_XYZ' : True,
                        'PRINTING_PATTERNS' : True,
                        'length_to_modify': 0
    }
try:
    with open('/home/k1899676/Downloads/Metadata.csv', "rb") as temp:
        R_Cut = pickle.load(temp)['R_Cut']
        R_Cut = numpy.array(R_Cut, dtype = float)[~numpy.isnan(numpy.array(R_Cut, dtype = float))]
except EOFError:
    R_Cut = 3.5
Ele = read(System['base_dir']+System['movie_file_name'], index=0).get_chemical_symbols()
with open(System['base_dir'] + 'Temp.xyz', "w+") as moviefile:
    moviefile.write(str(len(Ele)) + '\n')
    moviefile.write("Francesca's CNA Patterns \n")

    for i, t in enumerate(range(System['Start'], System['End'], System['Step'])):

        sigs = signature(System, CNA_Pattern_Settings)
        try:
            sigs.Frame_CNA_Sigs(frame=t, R_Cut = R_Cut[int(i/50)])
        except TypeError:
            sigs.Frame_CNA_Sigs(frame=t, R_Cut = 3.5)
        
        Pattern_Dict = Utilities.Bulk_Masterkey()
        Pattern_Dict = Pattern_Dict.Key()
        MasterKey = Utilities.CNA_Masterkey()
        MasterKey = MasterKey.Key()
        #frame = 0, System = None, Pattern_Input = None, MasterKey = None
        
        Pattern_Dict=patterns(frame = t, System = System, 
                              MasterKey = MasterKey, 
                              Pattern_Input = CNA_Pattern_Settings)
        
        
        XYZ = read(System['base_dir']+System['movie_file_name'], index = t).positions
        XYZ = numpy.column_stack((Ele,XYZ))
    
        Patterns = numpy.load(System['base_dir']+CNA_Pattern_Settings['npz_dir']+'/pattern_dictionary.npz', 
                              allow_pickle=True)['%s-%s'%(System['movie_file_name'][:-4],t)]
        Pats=numpy.zeros(len(XYZ))
        for i, atom in enumerate(Patterns):
            for j, val in enumerate(atom):
                if val:
                    Pats[i] = j+1
    
    
        Temp = numpy.column_stack((XYZ, Pats))
        for items in Temp:
            moviefile.write(' \t'.join(str(item) for item in items) + '\n')
        moviefile.write(str(len(XYZ))+'\n')
        moviefile.write('\n')
