import numpy as np
import os
import networkx as nx
from Sapphire.CNA import Utilities
from Sapphire.Utilities import errors
from Sapphire.Utilities.log import get_logger

log = get_logger('Sapphire.CNA.FrameSignature')

class CNA(object):
    """
    RMJ 10/04/22
    Class template structure on calculating common neighbour analysis (CNA)
    signatures and patterns with the former of the form (r,s,t), and the latter
    [n_i (r_i, s_i, t_i)] over all recognised local atomic signature indices i.
    
    
    r is the number of nearest neighbours common to both atoms in the pair; 
    s is the number of bonds between shared neighbours;
    t is the longest chain which can be made from bonding s atoms if they are nearest neighbours.;
    
    Parameters
    ----------
    System : Full Sapphire calculation information regarding base directories and file composition.
    
    Adj : scipy sparse matrix - returned from Post_Process.Adjacent.ReturnAdj()
        the 1st param name Adj

    Masterkey : tuple - The user may provide their own cna masterkey if they wish to compare
                        against a theoretical cna signature list
    
    Fingerprint : boolean - Whether or not to compute the cna patterns
    
    Type : boolean - Whether or not to write out the full cna profile to an external file

    Returns
    -------
    numpy array
        2 x m matrix where m is the number of unique recognised signatures 
        Will be of the form (n, (r,s,t)) where n is the number of unique counts
        of the signature (r,s,t)
    
    list
        N dimensional list where N is the number of atoms considered.
        Each list element will be a tuple of the CNA pattern for that given atom

    """
    
    def __init__(self, System = None, Adj = None, Masterkey = None,
                         Fingerprint = True, Type = False, Frame = 0):
        
        self.System = System
        self.Type = Type
        self.Frame = Frame
        if Adj is not None:
            try:
                # np.asarray: .todense() gives np.matrix, whose 2-D column slices break scalar
                # conversion under numpy>=2 (and were O(N) per lookup anyway).
                self.Adj = np.asarray(Adj.todense())
            except Exception as e:
                errors.report(e, 'CNA adjacency densification failed: %s', base_dir=self.System['base_dir'] if self.System else None)
        else:
            pass
        if Fingerprint:
            self.Fingerprint = np.zeros(self.Adj.shape[0], dtype = object)
            self.Keys = np.zeros(self.Adj.shape[0], dtype = object)    
        else:
            self.Fingerprint = False

        if Masterkey is None:
            self.Masterkey = ((0,0,0),
                        (1,0,0),
                        (2,0,0),(2,1,1),
                        (3,0,0),(3,1,1),(3,2,2),
                        (4,0,0),(4,1,1),(4,2,1),(4,2,2),(4,3,3),(4,4,4),
                        (5,2,1),(5,2,2),(5,3,2),(5,3,3),(5,4,4),(5,5,5),
                        (6,6,6))
        else:
            self.Masterkey = Masterkey
        self.Sigs = {}
        for item in self.Masterkey:
            self.Sigs[item] = 0

        self.Pat_Key = Utilities.Pattern_Key().Key() #Calling dictionary of recognised patterns
        self.Keys = list(self.Pat_Key)
        self.Max_Label = max(len(str(label)) for label in self.Keys)

    def ensure_dir(self, base_dir='', file_path=''):
        """
            A simple script to verify the existence of a directory
            given the path to it. If it does not exist, will create it.

        """

        directory = base_dir + file_path
        if not os.path.exists(directory):

            os.makedirs(directory)

    def MakeFile(self, Attributes):
        
        """
            Checks for the existence of a file to write output to. Creates one if it doe not exist.

        """
        
        self.out = self.System['base_dir'] + Attributes['Dir'] + Attributes['File']

        if not os.path.isfile(self.out):
            with open(self.System['base_dir'] + Attributes['Dir'] + Attributes['File'], 'w') as out:
                out.close()
        else:
            pass
        
    def Ascii_Bars(self, Finger):
        with open(self.System['base_dir'] + 'CNA_Pattern_Info.txt', 'a', encoding='utf-8') as f:
            f.write('\nCNA Pattern distribution for full system at frame %s.\n' %self.Frame)
        Temp = np.zeros(len(self.Keys), int)
        for atom in Finger:
            if str(atom) in self.Keys:
                Temp[self.Keys.index(str(atom))] += 1
        with open(self.System['base_dir'] + 'CNA_Pattern_Info.txt', 'a', encoding='utf-8') as f:
            for i, count in enumerate(Temp):
                bar_chunks, remainder = divmod(int(count * 8 / (len(Finger) / 50)), 8)
                # First draw the full width chunks
                bar = '█' * bar_chunks
                # Then add the fractional part.  The Unicode code points for
                # block elements are (8/8), (7/8), (6/8), ... , so we need to
                # work backwards.
                if remainder > 0:
                    bar += chr(ord('█')+ (8 - remainder))
                # If the bar is empty, add a left one-eighth block
                bar = bar or  '|'
                f.write(f'{self.Keys[i].rjust(self.Max_Label)} | {count:#4d} {bar}\n')

    def NN(self, atom):
        """

        Parameters
        ----------
        atom : integer
            the atomic index being considered relative to the ordering or atoms in the frame of the trajectory

        Returns
        -------
        self.neigh : list
            indices of all atoms considered to be neighbours of the reference atom.

        """
        
        self.neigh = np.flatnonzero(self.Adj[:, atom] == 1).tolist()
        return self.neigh
        
    
    def R(self, atom, friend):
        
        """

        Parameters
        ----------
        atom : integer
            the reference atomic index being considered relative to the ordering or atoms in the frame of the trajectory
        friend : integer
            the neighbour atomic index being considered relative to the ordering or atoms in the frame of the trajectory

        Returns
        -------
        self.bonds : list
            indices of all atoms which are mutually bonded to both the atom and its friend

        """
        
        self.bonds = np.flatnonzero((self.Adj[:, atom] == 1) & (self.Adj[:, friend] == 1)).tolist()
        self.r = len(self.bonds)
        return self.r
                    
    
    def S(self):
        self.s = 0
        self.perm = []
        for i, b in enumerate(self.bonds):
            for j, c in enumerate(self.bonds[i:]):
                a = int(self.Adj[c, b])
                if a == 1:
                    self.s += a
                    self.perm.append((b,c))
        return self.s

    def T(self):
        """Longest chain among the bonds between common neighbours.

        Same definition as the original networkx implementation: the maximum over all simple
        paths (in edges) and, where the bond graph has cycles, the longest cycle in a cycle basis.
        The path search is a small DFS; networkx is only used for the (rare) cyclic case.
        """
        nodes = list(self.bonds)
        nbrs = {n: set() for n in nodes}
        for b1, b2 in self.perm:
            if b1 != b2:
                nbrs[b1].add(b2); nbrs[b2].add(b1)
        n_edges = sum(len(v) for v in nbrs.values()) // 2
        best = 0
        if n_edges:
            def dfs(node, visited, length):
                nonlocal best
                best = max(best, length)
                for m in nbrs[node]:
                    if m not in visited:
                        visited.add(m); dfs(m, visited, length + 1); visited.discard(m)
            for start in nodes:
                dfs(start, {start}, 0)
            if n_edges >= len(nodes):  # cyclic
                G = nx.Graph(); G.add_nodes_from(nodes); G.add_edges_from((a, b) for a in nbrs for b in nbrs[a] if a < b)
                cycles = [len(c) for c in nx.cycle_basis(G)]
                if cycles:
                    best = max(best, max(cycles))
        self.t = best
        return self.t


    def calculate(self):
        for i, atom in enumerate(self.Adj):
            self.particle_cnas = []
            self.NN(i)
            for neigh in self.neigh:
                sig = tuple((self.R(i,neigh), self.S(), self.T()))
                try:
                    self.Sigs[sig]+=1
                except KeyError:
                    self.Sigs[sig] = 1
                    
                self.particle_cnas.append(sig)
            if self.Fingerprint is not False:
                self.Fingerprint[i] = self.fingers()  
        self.write()
        
        
    def fingers(self):
        Temp = set(self.particle_cnas)
        self.Keys.append(Temp)
        return tuple((self.particle_cnas.count(x), x) for x in Temp) 

    """
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
    """
    
    def write(self):
        if self.System is None:  # standalone use (tutorials): nothing to write to
            return None
        if self.Type == 'Full':
            from Sapphire.IO import OutputInfoFull as Out  # Case 1
            
            #Write object for the CoM
            Attributes = getattr(Out, str('cna_sigs')) #Loads in the write information for the object 
            OutFile = self.System['base_dir'] + Attributes['Dir'] + Attributes['File']
            self.ensure_dir(base_dir=self.System['base_dir'], file_path=Attributes['Dir'])   
            self.MakeFile(Attributes)
            with open(OutFile, 'a') as outfile:
                outfile.write(str(self.Frame) + ' ' +  ' '.join(str(item) for item in np.array(list(self.Sigs.values()), dtype = int)) +'\n')
            
            if self.Fingerprint is not False:
   
                #Write object for the homo CoM distances
                Attributes = getattr(Out, str('pattern_indices')) #Loads in the write information for the object                  
                OutFile = self.System['base_dir'] + Attributes['Dir'] + Attributes['File']
                self.ensure_dir(base_dir=self.System['base_dir'], file_path=Attributes['Dir'])   
                self.MakeFile(Attributes)
                with open(OutFile, 'a') as outfile:
                    outfile.write(str(self.Frame) + ' ' +  ' '.join(str(item) for item in self.Fingerprint) +'\n') 
                    
        from Sapphire.IO import OutputInfoExec as Out
        
        Attributes = getattr(Out, str('masterkey')) #Loads in the write information for the object 
        OutFile = self.System['base_dir'] + Attributes['Dir'] + Attributes['File']
        self.ensure_dir(base_dir=self.System['base_dir'], file_path=Attributes['Dir'])   
        self.MakeFile(Attributes)
        keys = list(self.Sigs.keys())
        # overwrite: the masterkey only grows, so the latest frame's key is the complete one
        with open(OutFile, 'w') as outfile:
            for item in keys:
                outfile.write(''.join(str(index) for index in item ) +' ')
                