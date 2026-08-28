"""Lindemann-index and CoM-distance helpers found pasted (unindented, no `self`) inside
Graphing/Plotter.py. Depend on undefined globals `distance`, `read`, `path`, `BigMeta`, `mp`, `partial`.
Moved here 2026-08-28."""

def Lin_Func(Dist_List):
    return np.sqrt( np.average( [ a**2 for a in Dist_List ] ) - np.average(Dist_List)**2 )/ np.average(Dist_List)


def Lin_List(Movie, index1, index2):
    List = []
    for T in range(0,len(Movie),10):
        atoms = Movie[T].get_positions()
        List.append(distance(atoms[index1],atoms[index2]))
    return Lin_Func(List)

def Func(Movie, index1):
    A = [ Lin_List(Movie, index1, x) for x in range(0,891) if x!= index1 ]
    return sum(A)/len(A)

def Main_Lind(Movie, N_Atoms, CPUs):
    iterable = range(0,N_Atoms)
    pool = mp.Pool(CPUs)
    function = partial(Func, Movie)
    New_Sample = pool.map(function, iterable)
    pool.close()
    pool.join()
    return New_Sample

def CoMDist(Sim, Frame):
    Pos = read(path+'Sim-'+str(Sim)+'/movie.xyz', index = Frame*10).get_positions()
    CDist = [ distance(BigMeta[Sim]['CoM'][Frame], Pos[i]) for i in range(len(Pos)) ]
    return CDist

def Main_CoM(Sim, CPUs):
    iterable = range(1000)
    pool = mp.Pool(CPUs)
    function = partial(CoMDist, Sim)
    New_Sample = pool.map(function, iterable)
    pool.close()
    pool.join()
    return New_Sample
