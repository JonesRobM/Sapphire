from ase.io import read
import numpy as np
from asap3.analysis import FullCNA


def get_all_cnas(filename, r_cut):
    """ (Claudio)
    Given a trajectory and a cutoff radius, returns a dictionary, sorted by value,
    with all the cnas that appear in the trajectory as keys and the number
    of times they appear as value.
    """
    
    
    """ Jones
    
    Arguments:
        filename: The string name of your input xyz file
            Normally something like 'movie.xyz'
            
        r_cut: The nearest neighbour spacing as a float
            This could be passed from a higher function
            or simply stated at the start of a run.
            
    Returns:
        sorted_cnas_dict: A dictionary whose keys are all of the identified CNA
        signtures appearing in the trajectory.
        The value associated with each key is the number of occurances.
        
    
    If one wishes to observe the entire distribution of CNA signatures over the 
    full production run. Simply call this function and plot its dictionary output.
    
    
    I personally use this to initially identify all of the signatures so that I may
    then identify frame-wise distributions for the cNA signatures.
    """
    
    all_cnas = {}
    traj=read(filename, index=':')
    for j, atoms in enumerate(traj):
        CNA = FullCNA(atoms, r_cut)
        atoms.set_cell([[100,0,0],[0,100,0],[0,0,100]])
        snapshot_cna = CNA.get_normal_cna()
        for i, atomic_cna in enumerate(snapshot_cna):
            for key in atomic_cna:
                try:
                    all_cnas[key] += atomic_cna[key]
                except KeyError:
                    all_cnas[key] = atomic_cna[key]

    sorted_cnas = sorted(all_cnas.items(), key=lambda kv: -kv[1])
    sorted_cnas_dict = {}
    for t in sorted_cnas:
        sorted_cnas_dict[t[0]] = t[1]

    return sorted_cnas_dict
 


def Master(filename, R_Cut):
    
    """ Jones
    
    Arguments: 
        filename: The string name of your input xyz file
            Normally something like 'movie.xyz'
        
        R_Cut: The nearest neighbour spacing as a float
            This could be passed from a higher function
            or simply stated at the start of a run.
           
    Returns:
        MasterKey:
            The sorted list containing all of the CNA signatures which appear in the 
            xyz file under consideration.
            
    """
        
            
    CNAS=get_all_cnas(filename,R_Cut)
    MasterKey=[]
    for keys in CNAS:
        MasterKey.append(keys)
    CNAS=0
    MasterKey.sort()
    return MasterKey
    



def get_cnas(frame, R_Cut, Masterkey=None, filename=None):
    """(Claudio)
    Given a trajectory and a cutoff radius, returns a dictionary, sorted by value,
    with all the cnas that appear in the trajectory as keys and the number
    of times they appear as value.
    """
    
    """ Jones
    
    Arguments:
        filename: The string name of your input xyz file
            Normally something like 'movie.xyz'
        
        R_Cut: The nearest neighbour spacing as a float
            This could be passed from a higher function
            or simply stated at the start of a run.
        
        j: Integer whicha specifies the frame of the trajectory file to be called
        
    Returns:
        (Key, Num) The tuple of the CNA signature alongside the number of times it has
        been observed in this given frame.
        
    In general, this function is to be called by the CNA_Sig_Frame and will proceed to be 
    processed by the follow-up guy.
    
    """
        
        
    
    all_cnas = {}
    traj=read(filename, index=frame)
    traj.set_cell([[100,0,0],[0,100,0],[0,0,100]])
    CNA = FullCNA(traj, R_Cut)
    snapshot_cna = CNA.get_normal_cna()
    for i, atomic_cna in enumerate(snapshot_cna):
        for key in atomic_cna:
            try:
                all_cnas[key] += atomic_cna[key]
            except KeyError:
                all_cnas[key] = atomic_cna[key]

    sorted_cnas = sorted(all_cnas.items())
    sorted_cnas_dict = {}
    for t in sorted_cnas:
        sorted_cnas_dict[t[0]] = t[1]

    Key=[]; Num=[]
    for keys in sorted_cnas_dict:
        Key.append(keys); Num.append(sorted_cnas_dict[keys])
        
        if keys not in Masterkey:
            Masterkey.append(keys)
    
    return (Key, Num)


def CNA_Sig_Frame(filename, MasterKey, R_Cut, Frames, Skip):
    
    """ Jones
    
    Arguments:
        filename: The string name of your input xyz file
            Normally something like 'movie.xyz'
        
        MasterKey: The output from calling the Master function.
        This is to do pairwise comparrison for creating full 
        distributions without having to know what the craic is.
            
        
        R_Cut: The nearest neighbour spacing as a float
            This could be passed from a higher function
            or simply stated at the start of a run.
        Frames: (int) The number of frames f your movie  file you wish to 
        consider up to
        
        Skip: (int) How many frames you wish to pass over before recording new data.
        
    Returns:
        
        Heights: np.array(Frames/Skip, len(MasterKey)) The array containing the normalised
        distribution of CNA signature occurances. 
        
    """
    
    FullList=[]
    FullSample=[]
    Heights=[]

    for frame in range(int(Frames/Skip)):
        FullList.append(get_cnas(filename,R_Cut,Skip*frame))
        Temp1=FullList[frame][0]
        for x in MasterKey:
            if x not in Temp1:
        
                FullList[frame][0].append(x)
                FullList[frame][1].append(0)
        Sample=[]
        for j in range(len(MasterKey)):
            Sample.append((FullList[frame][0][j],FullList[frame][1][j]))
        FullSample.append(Sample)
        FullSample[frame].sort()
        A,B=zip(*FullSample[frame])
        Heights.append(B/np.sum(B))
        
    return Heights
        

def get_heights(CNA, Masterkey, frame):
    
    def getkey(item):
        return item[0]
    
    FullSample=[]; Heights=[]
    Temp1=CNA[frame][0]
    for x in Masterkey:
        if x not in Temp1:    
            CNA[frame][0].append(x)
            CNA[frame][1].append(0)
    Sample=[]
    for j in range(len(Masterkey)):
        Sample.append((CNA[frame][0][j], CNA[frame][1][j]))
    FullSample.append(Sample)
    FullSample.sort()
    A,B=zip(*FullSample[0])
    Heights.append(B/np.sum(B))
    Temp = FullSample[0]
    FullCNA = [ item[0] for i, item in enumerate(Temp) ]
    
    Sample = (FullCNA, Heights[0])
    #Sample = sorted(Sample, key = getkey)        
    
    return (FullCNA, Heights[0])