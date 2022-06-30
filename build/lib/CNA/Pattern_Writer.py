from ase.io import read
import pickle
import numpy as np
import matplotlib.pyplot as plt

def Patterns(file='Metadata.csv', xyz = 'Strut.xyz', frame=0, outfile = 'Patterns'):
    with open(file, "rb") as infile:
        Data = pickle.load(infile)['cna_patterns'][frame]
    pos = read(xyz, index = frame).positions
    ele = read(xyz, index = frame).get_chemical_symbols()
    Temp = {}
    Patterns_Listed = []
    for x in Data:
        Temp[str(x)] = 0
    for x in Data:
        Temp[str(x)] += 1
    Numbered_Pats = Temp.copy()
    for i,x in enumerate(Numbered_Pats):
        Numbered_Pats[str(x)] = i
    for x in Data:
        Patterns_Listed.append(Numbered_Pats[str(x)])
    New = np.column_stack((ele,pos,Patterns_Listed))
    with open(outfile+'.xyz', 'w') as newfile:
        newfile.write(str(len(New)))
        newfile.write("\n \n")
        for atom in New:
            newfile.write(' \t'.join(str(item) for item in atom) + '\n')
    X_Pats = [ str(a) for a in Temp ]
    Plotting = [ Temp[x] for x in Temp ]
    cmap = plt.get_cmap('rainbow')
    #rescale = lambda Plotting: (Plotting - np.min(Plotting)) / (np.max(Plotting) - np.min(Plotting))
    fig, ax  = plt.subplots()
    fig.set_size_inches(12,8)
    X = np.array(range(len(X_Pats)))
    Plot = ax.barh(X, Plotting, color = cmap( X / (len(X)-1) ) )
    ax.set_ylabel('CNA Pattern', fontsize = 14)
    ax.set_xlabel('Frequency', fontsize = 14)
    plt.yticks(X, X_Pats)
    ax.set_title(xyz[3:-4], fontsize=14)
    plt.savefig(outfile+'.jpeg', dpi = 100, bbox_inches='tight')