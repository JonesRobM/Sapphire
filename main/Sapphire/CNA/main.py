# Import OVITO modules.
from ovito.io import *
from ovito.modifiers import *
from ovito.data import *

# Import standard Python and NumPy modules.
import sys
import numpy as np

#import glob and os
import glob
import os

#importing python modules
import cna_signatures
import cna_patterns
import xyz_io

"""
(Armand)
The idea of this code is to be a standalone CNA pattern code. It draws huge
inspiration from Robert's PP code, but due to issues with the r_cut, I ended up
making it a standalone code to go faster. The r_cut used for the CNA is now just
based on a default value based on the elemental composition. This has issues of
course, but I dont have time to fix it, I had to start getting results. Adapting
this code to use Robert PP's code r_cut from the PDF shouldn't be too hard, but
I leave it to another soul. I don't have time to do that.

The code works the following way:

-Looping over all xyz files found within System['xyz_dir'], and only looking at
    the first frame
-CNA signatures are calculated, and saved within an npz file (numpy archives)
    within the System['npz_dir'].
-Unique CNA patterns are then identified for all xyz files, and a masterkey of
    them is made. This masterkey can be chosen to NOT be updated, allowing to
    construct a masterkey from perfect geometries, and use it on messier xyz
    files.
-A Pattern Dictionary is made, which contains the cna pattern masterkey found
    and which contains boolean arrays identifying if an atom of the nanoparticle
    had a CNA pattern found within the masterkey. The keys are the filenames of
    the xyz files.
-The Pattern Dictionary is saved, as I needed to work on it without modifying it
    from past usages.
-A new xyz file is created with values between 0 and len(pattern_masterkey) for
    a recognized cna pattern, and a value of 999 for unfound values. This can be
    used with ovito's color coding to identify the CNA pattern identified.

This code is a bowl of spaghetti and a massive bodge to make it work asap, as
I needed results QUICKLY, to start using SVC on them.
"""

System = {
        'base_dir' : '',
        'xyz_dir' : 'XYZs_to_compute/',
        'npz_dir' : 'CNA_npz/', #folder to save the npz files in
        'new_xyz_dir' : 'CNA_XYZs/', #folder for new xyz files for ovito colouring
        'movie_file_name' : '', #LEAVE BLANK
        'length_to_modify': 0   #DONT TOUCH, INITIALIZED AUTOMATICALLY
        }

MasterKey = ((0,0,0),
            (1,0,0),
            (2,0,0),(2,1,1),
            (3,0,0),(3,1,1),(3,2,2),
            (4,0,0),(4,1,1),(4,2,1),(4,2,2),(4,3,3),(4,4,4),
            (5,2,1),(5,2,2),(5,3,2),(5,3,3),(5,4,4),(5,5,5),
            (6,6,6))

System['length_to_modify'] = len(System['base_dir'])+len(System['xyz_dir'])

print('Input the r_cut to use for all files (float)')
r_cut=input('(recommended value is 3.5 for Au) ')
try:
    r_cut=float(r_cut)
except ValueError:
    print('Error in value. Using a cutoff of 3.5.')
    r_cut=3.5


print('Should the Pattern Dictionary be appended ? (Y/N)')
a=input()
if(a == 'Y' or a == 'y'):
    APPEND_DICTIONARY = True
    FROM_MEMORY = True
    BULK_MASTERKEY = False
    SAVING_XYZ = False
elif(a == 'N' or a == 'n'):
    APPEND_DICTIONARY = False
else:
    print('Error in input, appending all new values to the dictionary.')
    APPEND_DICTIONARY = True

if(APPEND_DICTIONARY is False):
    print('\nShould the Pattern Dictionary be entirely remade? (Y/N)')
    a=input()
    if(a == 'Y' or a == 'y'):
        NEW_DICTIONARY = True
    elif(a == 'N' or a == 'n'):
        NEW_DICTIONARY = False
    else:
        print('Error in input, making a new pattern dictionary')
        NEW_DICTIONARY = True


    if(NEW_DICTIONARY is True):
        print('Should the Bulk Masterkey be used? (Y/N)')
        a=input()
        if(a =='Y' or a == 'y'):
            BULK_MASTERKEY = True
            FROM_MEMORY = False
        elif(a =='N' or a =='n'):
            BULK_MASTERKEY = False
        else:
            print('Error in input, using bulk masterkey')
            BULK_MASTERKEY = True
            FROM_MEMORY = False

        if(BULK_MASTERKEY is False):
            print('Should the Pattern Masterkey be loaded from before or remade? (L/R)')
            a=input()
            if(a == 'L' or a == 'l'):
                FROM_MEMORY = True
            elif(a =='R' or a =='r'):
                FROM_MEMORY = False
            else:
                print('Error in input, loading old pattern masterkey')
                FROM_MEMORY = False

    print('Should new xyz files be saved? (Y/N)')
    a=input()
    if(a =='Y' or a == 'y'):
        SAVING_XYZ = True
    elif(a =='N' or a == 'n'):
        SAVING_XYZ = False
    else:
        print('Error in input, not saving xyz files.')
        SAVING_XYZ = False


print('Should all found Patterns be printed to the console? (Y/N)')
a=input()
if(a =='Y' or a == 'y'):
    PRINTING_PATTERNS = True
elif(a =='N' or a == 'n'):
    PRINTING_PATTERNS = False
else:
    print('Error in input, not printing the patterns.')
    PRINTING_PATTERNS = False




xyz_folder_path = System['base_dir']+System['xyz_dir']
print('Calculating CNA signatures: ')
for filename in glob.glob(os.path.join(xyz_folder_path, '*.xyz')):
    System['movie_file_name']=filename[System['length_to_modify']:]
    xyz_file=System['base_dir']+System['xyz_dir']+System['movie_file_name']
    cna_signatures.CNA_calculations(xyz_file, MasterKey, System, cutoff=r_cut)
    print('CNA signatures of',System['movie_file_name'],'done.')


if(APPEND_DICTIONARY is True):
    #Loading the old Dictionary
    Temp_Dict=np.load(System['base_dir']+System['npz_dir']+'pattern_dictionary.npz', allow_pickle=True)
    #Rewriting it on a dictionary to be able to append/modify values
    Pattern_Dict={}
    for file in Temp_Dict.files:
        Pattern_Dict[file]=Temp_Dict[file]

elif(NEW_DICTIONARY is True):
    Pattern_Dict={}

Pattern_Dict=cna_patterns.pattern_dictionary_maker(System, MasterKey,
                                                    Pattern_Dict,
                                                    FROM_MEMORY=FROM_MEMORY,
                                                    BULK_MASTERKEY=BULK_MASTERKEY)
xyz_io.dictionary_saver(System, Pattern_Dict)



#reading all filenames
folder_path = System['base_dir']+System['xyz_dir']
for filename in glob.glob(os.path.join(folder_path, '*.xyz')):

    System['movie_file_name']=filename[System['length_to_modify']:]
    arrays_filename=System['npz_dir']+'/CNA_'+System['movie_file_name'][:-4]+'.npz'

    if (PRINTING_PATTERNS==True):
        print('\nCurrent file: '+System['movie_file_name'])
        cna_patterns.pattern_CNA_Reader(arrays_filename,MasterKey)

    #Grabbing all the movie file names
    PATH_TO_MOVIE = System['base_dir']+System['xyz_dir']+System['movie_file_name']
    PATH_TO_NEW_MOVIE = System['base_dir']+System['new_xyz_dir']+System['movie_file_name'][:-4]+'_CNA.xyz'

    #Saving the new XYZ files
    if(SAVING_XYZ==True):
        print('Saving: '+System['movie_file_name'])
        xyz_io.XYZ_data_maker(Pattern_Dict, System,
                                    PATH_TO_MOVIE, PATH_TO_NEW_MOVIE)
