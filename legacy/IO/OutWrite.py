from ase.io import read
import numpy as np
import pickle

def ExtendXYZ(Traj, Metadata, Quants, Names):
    with open('Extend.xyz', 'w') as movie:
        movie.write(str(Metadata['NAtoms'][0]) +'\n')
        movie.write('Extra columns are | \t')
        for name in Names:
            movie.write(str(name) + '\t')
        movie.write('\n')
        for i, Frame in enumerate(range(len(Metadata['agcn']))):
            items = np.column_stack(( Traj[i].get_chemical_symbols(), 
                                     Traj[i].positions ))
            for obj in Quants:
                items = np.column_stack((items, obj[i]))
            for atom in items:
                movie.write(' \t'.join(str(item) for item in atom) +'\n')
            movie.write(str(Metadata['NAtoms'][i]) + '\n')
            movie.write('\n')
            
def count(Input, Value):
    return(len([x for x in Input if x == Value]))

def PtInfo(Metadata):
    AvgPt = []
    Mat = np.zeros((len(Metadata['nn']),13))
    for t in range(len(Metadata['nn'])):
        for a in range(13):
            Mat[t][a] += count(Metadata['headj'][t][1], a)
        AvgPt.append(np.average(Metadata['nn'][t][-55:],0))
    N = np.column_stack((Mat,AvgPt))
    np.savetxt('PtInfo.dat', N)

def AuNN(Metadata):
    Temp = []
    for t in range(len(Metadata['headj'])):
        Temp.append(Metadata['hoadjAu'][t] + Metadata['headj'][t][1])
    return Temp

def PtNN(Metadata):
    Temp = []
    for t in range(len(Metadata['headj'])):
        Temp.append(Metadata['headj'][t][0] + Metadata['hoadjPt'][t])
    return Temp

def Output():
    Traj = read('NewMovie.xyz', index = ':')
    with open('Metadata.csv', 'rb') as file:
        Metadata = pickle.load(file)
    AuNeigh = AuNN(Metadata); PtNeigh = PtNN(Metadata); NN = Metadata['nn']
    ExtendXYZ(Traj, Metadata, [NN, AuNeigh, PtNeigh], ['Coordination', 'Au Neighbours', 'Pt Neighbours'])
    PtInfo(Metadata)
    
    np.savetxt('Gyration_Mix.dat', np.column_stack((np.column_stack((
    Metadata['gyration'], Metadata['gyrationPt'])), Metadata['mix'])))
    
    
if __name__ == '__main__':
    Output()