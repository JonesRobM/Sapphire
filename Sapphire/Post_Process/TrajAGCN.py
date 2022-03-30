import pickle
import numpy as np
from ase.io import read

class Trajagcn():


    def __init__(self, metadata, system):
        self.system = system
        self.metadata = metadata
        self.movie_file = self.system['base_dir']+self.system['movie_file_name']
        self.start = self.system['Start']
        self.end = self.system['End']
        self.step = self.system['Step']
        try:
            self.agcn = self.metadata['agcn']
        except KeyError:
            print("It would appear that you have not evaluated the agcn for this simulation, yet. \n")
    
    def edit_movie(self):
        self.New_Obj = []
        self.Frames = range(self.start, self.end, self.step)
        for i in self.Frames:
            temp = read(self.movie_file, index = i)
            c = temp.get_chemical_symbols()
            xyz = temp.get_positions()
            ag = self.agcn[Frames.index(i)]
            self.New_Obj.append( np.column_stack( (c,xyz,ag) ) )
        return self.New_Obj
        

    
    def New_File(self, new_movie='agcn_movie.xyz'):
        with open(self.system['base_dir'] + new_movie, 'w+') as self.movie:
            self.movie.write(str(self.metadata['NAtoms']) +'\n')
            self.movie.write('\t' + "This was made by Jones' post-processing code." + '\n')
            for Frame in self.New_Obj:
                for items in Frame:
                    self.movie.write(' \t'.join(str(item) for item in items) +'\n')
                self.movie.write(str(self.metadata['NAtoms']) + '\n')
                self.movie.write('\n')
        print("This movie has been saved as %s in %s.\n" %(new_movie, self.system['base_dir']))