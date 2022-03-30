import numpy as np
from ase.io import read

XYZ = read('/media/k1899676/Seagate/PhD/20/March/CMD/Au/1103/NVT/600/Sim-1345/movie.xyz').positions
Ele = read('/media/k1899676/Seagate/PhD/20/March/CMD/Au/1103/NVT/600/Sim-1345/movie.xyz').get_chemical_symbols()
XYZ = np.column_stack((Ele,XYZ))
with open('Temp.xyz', "w+") as moviefile:
    moviefile.write('1103' + '\n')
    moviefile.write('Core-shell ico AuPt + Delta Bader charges \n')

    Patterns = np.load('pattern_dictionary.npz', allow_pickle=True)['movie-0']
    Pats=np.zeros(1103)
    for i, atom in enumerate(Patterns):
        for j, val in enumerate(atom):
            if val:
                Pats[i] = j+1

         
    Temp = np.column_stack((XYZ, Pats))
    for items in Temp:
        moviefile.write(' \t'.join(str(item) for item in items) + '\n')
