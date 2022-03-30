import sys
import numpy as np

import os
import time

from CNA import Utilities

sys.path.append("../../")

class patterns():
    
    def __init__(self, frame = 0, System = None, Pattern_Input = None, MasterKey = None):
        
        tick = time.time()
        self.frame = frame
        self.System = System
        self.Pattern_Input = Pattern_Input
        if self. System is not None:
            self.filename = System['base_dir']+System['movie_file_name']
            self.npz_dir = System['base_dir'] + 'CNA_npz'
            with open(self.System['base_dir'] + 'CNA_Pattern_Info.txt', 'a') as f:
                f.write('\n')
                f.write(' # '*50)
                f.write('\nComputing CNA patterns for frame %s.\n'%frame)
                f.close()
        else:
            try:
                self.filename = 'movie.xyz'
            except FileNotFoundError:
                with open(self.System['base_dir'] + 'CNA_Pattern_Info.txt', 'a') as f:
                    f.write('\nCould not find a suitable file to examine.\n')
                    f.close()
        self.script_path = os.path.dirname(os.path.realpath(__file__))+'/'
        
        self.cwd = os.getcwd()
        
        if MasterKey is None:
            self.MasterKey = Utilities.CNA_Masterkey().Key()
        else:
            self.MasterKey = MasterKey
        
        if self.Pattern_Input['APPEND_DICTIONARY'] is True:
            os.chdir(self.script_path)
            os.chdir('../')
            self.Temp_Dict = np.load(
                                    'CNA_npz/pattern_dictionary.npz', 
                                    allow_pickle=True)
            self.Pattern_Dict = {}
            for file in self.Temp_Dict.files:
                self.Pattern_Dict[file] = self.Temp_Dict[file]
            
        elif (self.Pattern_Input['NEW_DICTIONARY'] is True):
            self.Pattern_Dict = {}

            
        self.Pattern_Dict = self.pattern_dictionary_maker()

        self.dictionary_saver()
        
        with open(self.System['base_dir'] + 'CNA_Pattern_Info.txt', 'a') as f:
            f.write("\nGenerating CNA Patterns took %.3f seconds.\n" %(time.time()-tick))
        
        os.chdir(self.cwd)
        
 
    def run(self):
        Info = self.Pattern_Dict[self.System['movie_file_name'][:-4]+'-'+str(self.frame)]
        Pats = np.zeros(len(Info))
        for i, atom in enumerate(Info):
            for j, val in enumerate(atom):
                if val:
                    Pats[i] = j+1
        return Pats
            

    def pattern_CNA_Reader(self):
    
        """
        Armand
        Formatting from the npz files, gives the cna patterns found and prints them.
        This isnt meant to be run normally, add it within the filename loop when you
        want to inspect a FEW files. Doesn't use the masterkey, so prepare to have a
        LOT of data printed at you at once.
        """
    
        self.CNA_arrays=np.load(self.npz_dir+'/CNA_'+self.filename[:-4]+'-'+str(self.frame)+'.npz', allow_pickle=True)
        
        with open(self.System['base_dir'] + 'CNA_Pattern_Info.txt', 'a') as f:
            f.write('\nTypes of CNA bonds found with each atom:\n')
            f.close()
    
        for i in range(len(self.CNA_arrays['signature_cna_count'])):
            with open(self.System['base_dir'] + 'CNA_Pattern_Info.txt', 'a') as f:
                f.write('\n%d Atoms had CNA patterns  (no: %d)\n'%(self.CNA_arrays['\nSignature_cna_count'][i], i))
                f.close()
                
            non_zero_values=np.nonzero(self.CNA_arrays['signature_cna'][i])
            for j in range(len(non_zero_values[0])):
                
                with open(self.System['base_dir'] + 'CNA_Pattern_Info.txt', 'a') as f:
                    f.write('\n%s on %s of its bonds.\n' %(self.MasterKey[non_zero_values[0][j]], 
                                                         self.CNA_arrays['signature_cna'][i][non_zero_values[0][j]]))
                    f.close()
                    
            with open(self.System['base_dir'] + 'CNA_Pattern_Info.txt', 'a') as f:
                f.write('\nCoordination number: %s\n'%np.sum(self.CNA_arrays['signature_cna'][i]))
                f.close()
    
    def cna_pattern_master_key_maker(self):
        
        """
        Armand
        This function creates a new cna pattern masterkey by running through ALL
        files within xyz_dir. This is meant for studying all cna patterns with the
        variable SAVING_XYZ == True, not for Support Vector Clustering.
        """
        
        self.CNA_arrays=np.load(self.npz_dir+'/CNA_'+self.filename[:-4]+'-'+str(self.frame)+'.npz', allow_pickle=True)
    
        self.cna_patterns=[]
        with open(self.System['base_dir'] + 'CNA_Pattern_Info.txt', 'a') as f:
            f.write("\nCreating the CNA pattern master key...\n")
            f.close()

        #Looping over the atoms
        for i in range(len(self.CNA_arrays['signature_cna_count'])):

            #Creating the atomistic pattern list
            self.atom_pattern=[]

            #Finding the non zero CNA signatures, and looping over them
            non_zero_values = np.nonzero(self.CNA_arrays['signature_cna'][i])
            for j in range(len(non_zero_values[0])):
                #Retrieving the CNA signature from the Master Key
                cna_sign = self.MasterKey[non_zero_values[0][j]]
                #Counting them
                count = self.CNA_arrays['signature_cna'][i][non_zero_values[0][j]]
                #Appending the tuples found within the list
                self.atom_pattern.append((cna_sign,count))
            #Checking if the atomic pattern is in the cna_pattern_masterkey
            if self.atom_pattern not in self.cna_patterns:
                self.cna_patterns.append(self.atom_pattern)
    
        #Ordering the pattern masterkey by the Coordination Number
        self.cna_pattern_array = np.asarray(self.cna_patterns)
        self.cna_pattern_master_key=np.copy(self.cna_pattern_array)
        a=[]
        for i in range(len(self.cna_pattern_array)):
            a.append(np.sum(self.cna_pattern_array[i],axis=0)[1])
        l=np.asarray(a).argsort()[::-1]
        for i in range(len(l)):
            self.cna_pattern_master_key[i]=self.cna_pattern_array[l[i]]
    
        #returning the pattern masterkey
        return self.cna_pattern_master_key
    
    def pattern_dictionary_maker(self):
        """
        Armand
        This is where the magic happens. The function first asks for a new MasterKey
        or receives one from memory. The function goes over all files within xyz_dir,
        and uses the npz files in npz_dir to find all of the atoms whose patterns
        are in MasterKey.
        """

    
        #READING THE MASTERKEY FROM PATTERN DICTIONARY NPZ
        if (self.Pattern_Input['FROM_MEMORY'] is True):
            
                self.Pattern_Dict['masterkey'] = np.load(
                    self.System['base_dir']  
                    + self.System['npz_dir']
                    + 'pattern_dictionary.npz', 
                    allow_pickle = True)['masterkey']
                
                with open(self.System['base_dir'] + 'CNA_Pattern_Info.txt', 'a') as f:
                    f.write()('\nKey CNA Patterns found in memory:\n')
                    f.close()
                
        #CREATING A NEW MASTERKEY
        elif (self.Pattern_Input['FROM_MEMORY'] is False):
            
            #USING THE MASTERKEY FOR SUPPORT VECTOR CLUSTERING
            
            if (self.Pattern_Input['BULK_MASTERKEY'] is True):
                
                self.Pattern_Dict = Utilities.Bulk_Masterkey(self.Pattern_Dict).Key()
                
                with open(self.System['base_dir'] + 'CNA_Pattern_Info.txt', 'a') as f:
                    f.write('\nUsing bulk pattern dictionary from Utilities.\n')                
                
            #CREATING A NEW MASTERKEY FROM ALL PATTERNS FOUND WITHIN THE THE XYZ_DIR
            
            if(self.Pattern_Input['BULK_MASTERKEY'] is False):
                
                self.Pattern_Dict['masterkey'] = self.cna_pattern_master_key_maker(self.System, self.MasterKey)
                
                with open(self.System['base_dir'] + 'CNA_Pattern_Info.txt', 'a') as f:
                    f.write('\nFound key CNA Patterns:\n')
                
        #printing it
        for key in self.Pattern_Dict['masterkey']:
            
            with open(self.System['base_dir'] + 'CNA_Pattern_Info.txt', 'a+') as f:
                f.write('\n')
                f.write('\t'.join(str(item) for item in key))
                f.write('\n')
                f.close()

        #Looping over all files again

        with open(self.System['base_dir'] + 'CNA_Pattern_Info.txt', 'a') as f:
            f.write('\nCalculating CNA Patterns of: '+self.System['movie_file_name']+'\n')
            f.write('\n Reading CNA arrays from:\n' + self.npz_dir)
            f.close()
        
        self.CNA_arrays=np.load(self.npz_dir+'/CNA_'+self.System['movie_file_name'][:-4]+'-'+str(self.frame)+'.npz', allow_pickle=True)

        #pattern_CNA_Reader(arrays_filename, MasterKey)

        #Loading the CNA arrays
        self.Pattern_Dict[self.System['movie_file_name'][:-4]+'-'+str(self.frame)] = np.zeros(
                                                                (len(self.CNA_arrays['particle_cnas']),
                                                                 len(self.Pattern_Dict['masterkey'])),dtype=bool)
        #Looping over the atoms
        for i in range(len(self.CNA_arrays['particle_cnas'])):

            #Creating the atomistic pattern list
            self.atom_pattern=[]

            #Finding the non zero CNA signatures, and looping over them
            self.non_zero_values = np.nonzero(self.CNA_arrays['particle_cnas'][i])
            for j in range(len(self.non_zero_values[0])):
                #Retrieving the CNA signature from the Master Key
                self.cna_sign = self.MasterKey[self.non_zero_values[0][j]]
                #Counting them
                self.count = self.CNA_arrays['particle_cnas'][i][self.non_zero_values[0][j]]
                #Appending the tuples found within the list
                self.atom_pattern.append((self.cna_sign, self.count))
            #Checking if the atomic pattern is in the cna_pattern_masterkey
            if self.atom_pattern in list(self.Pattern_Dict['masterkey']):
                k = list(self.Pattern_Dict['masterkey']).index(self.atom_pattern)
                self.Pattern_Dict[self.System['movie_file_name'][:-4]+'-'+str(self.frame)][i][k] = True
        
        return self.Pattern_Dict

    def dictionary_saver(self):
        #Saving the created dictionary in an npz file
        self.values_to_save={}
        for key in self.Pattern_Dict:
            #creating an argument dictionary to input in np.savez
            self.values_to_save[key]=self.Pattern_Dict[key]
        #saving the npz file
        os.chdir(self.script_path)
        os.chdir('../')
        self.path_to_npz = self.System['base_dir'] + 'CNA_npz/pattern_dictionary.npz'
        np.savez(self.path_to_npz, **self.values_to_save)
        
        with open(self.System['base_dir'] + 'CNA_Pattern_Info.txt', 'a') as f:
            f.write("\nPatterns saved in %s.\n"%('CNA_npz/pattern_dictionary.npz'))
            f.close()
            
    def movie_writer(self, Positions = None, Elements = None, Outfile = 'Pattern_Movie.xyz'):
        if not os.path.isfile(self.System['base_dir']+Outfile):
            with open(self.System['base_dir'] + Outfile, "w+") as moviefile:
                moviefile.write(str(len(Elements)) + '\n')
                moviefile.write("CNA Patterns \n")
                XYZ = np.column_stack((Elements, Positions))
                Pats = np.zeros(len(Positions))
                for i, atom in enumerate(self.Pattern_Dict[self.System['movie_file_name'][:-4]+'-'+str(self.frame)]):
                    for j, val in enumerate(atom):
                        if val:
                            Pats[i] = j+1
            
                Temp = np.column_stack((XYZ, Pats))
                for items in Temp:
                    moviefile.write(' \t'.join(str(item) for item in items) + '\n')
                moviefile.write(str(len(XYZ))+'\n')
                moviefile.write('\n')
        else:
            with open(self.System['base_dir'] + Outfile, "w+") as moviefile:
                   
                XYZ = np.column_stack((Elements, Positions))
                Pats = np.zeros(len(Positions))
                for i, atom in enumerate(self.Pattern_Dict[self.System['movie_file_name'][:-4]+'-'+str(self.frame)]):
                    for j, val in enumerate(atom):
                        if val:
                            Pats[i] = j+1
                       
                Temp = np.column_stack((XYZ, Pats))
                for items in Temp:
                    moviefile.write(' \t'.join(str(item) for item in items) + '\n')
                moviefile.write(str(len(XYZ))+'\n')
                moviefile.write('\n')