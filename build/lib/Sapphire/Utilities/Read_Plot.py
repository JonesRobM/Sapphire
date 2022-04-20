"""

This script is generally for the purpose of extracting quantities of interest
from the Sapphire-generated Metadata object as it exists  in V0.10.1
This will simply be a placeholder until a more robust way of writing and storing
the output data can be 

"""

import pickle
from ase.io import read
import numpy as np
import matplotlib.pyplot as plt

def distance(a, b):
    
    dx = abs(a[0] - b[0])
     
    dy = abs(a[1] - b[1])
     
    dz = abs(a[2] - b[2])
 
    return np.sqrt(dx**2 + dy**2 + dz**2)

def Collect_CNA(Data, Sig):
    """

    Parameters
    ----------
    Data : TYPE - List
        Must be of the form
        Read_Data = Reader(...)... etc...
        Data = Read_Data[2][Key] for whichever simulation
        you wish to extract the cna signatures from.
        Note that the index '2' corresponds to the extracted cna sigs.
    Sig : TYPE - Tuple
        Will be of the form (r, s, t) EXACTLY
        where r,s,t are the triplet of the desired signature

    Returns
    -------
    list
        A list of the normalised frequency of occurance for 
        a given signature to be observed in the given simulation.

    """
    try:
        Index = Data[0][1].index( Sig )
        #This will pull the index of the desired signature from the
        #first frame of the data.

        return [ Data[x][0][Index] for x in range(len(Data)) ]
    except Exception as e:
        print(e)
        return None

Sims = ['Sim-1345/', 'Sim-2783/', 'Sim-3987/', 'Sim-4009/']
Struts = ['Co/', 'Ih/']


def New_File(path, new_movie='Quantity_movie.xyz', Quantities = []):
    
    Reference = read(path, index = ':')
    
    """
    Robert:
        
        This function, at the moment, is only supporting the introduction of the aGCN
        to the new xyz file. 
        
        But this is easily appended to as needs dictate.     
    """
    
    with open(new_movie, 'w+') as movie:
        movie.write(str(len(Reference[0])) +'\n')
        movie.write('\t' + "This was made by Jones' post-processing code." + '\n')
        for i in range(len(Reference)):
            Ele = Reference[i].get_chemical_symbols()
            Pos = Reference[i].positions
            
            items = np.column_stack(( Ele, Pos))
            for obj in Quantities:
                items = np.column_stack((items, obj[i]))
            for atom in items:
                movie.write(' \t'.join(str(item) for item in atom) +'\n')
            movie.write(str(len(Ele)) + '\n')
            movie.write('\n')

def hebond(data):
    a = np.array([ sum(data[t][1]) for t in range(len(data)) ])
    
    b = np.array([ np.concatenate((data[t][0], data[t][1]), axis = 0) for t in range(len(data)) ])
    c = [ [0] * len(data[0][0]) ]; d = [ [0] * len(data[0][1]) ]
    for t in range(len(data)-1):
        c.append( [ data[t+1][0][x] - data[t][0][x] for x in data[t][0] ] )
        d.append( [data[t+1][0][x] - data[t][0][x] for x in data[t][1] ] )
    e = []
    f = []
    for t in range(len(data)):      
        e.append(np.concatenate((c[t], d[t])))
    e = np.array(e)
    c = [0]
    for t in range(len(data)-1):
        c.append( sum(data[t+1][1]) - sum(data[t][1]))
    c = np.array(c)
    d = np.array([data[t][1] for t in range(len(data)) ])
    return a,b,e,c,d

def Relative(headj, nn):
    return [ headj[t][1] / nn[t][309:] for t in range(len(nn)) ]

def Init():
    
    """
   

    Returns
    -------
    edelta : TYPE
        DESCRIPTION.
    comspace : TYPE
        DESCRIPTION.
    cna_sigs : TYPE
        DESCRIPTION.
    adj : TYPE
        DESCRIPTION.
    agcn : TYPE
        DESCRIPTION.
    com : TYPE
        DESCRIPTION.
    comdist : TYPE
        DESCRIPTION.
    surf_atoms : TYPE
        DESCRIPTION.
    comAu : TYPE
        DESCRIPTION.
    comPt : TYPE
        DESCRIPTION.
    hoadjAu : TYPE
        DESCRIPTION.
    hoadjPt : TYPE
        DESCRIPTION.
    comdistAu : TYPE
        DESCRIPTION.
    comdistPt : TYPE
        DESCRIPTION.
    midcomdistAu : TYPE
        DESCRIPTION.
    midcomdistPt : TYPE
        DESCRIPTION.
    surf_atomsPt : TYPE
        DESCRIPTION.
    headj : TYPE
        DESCRIPTION.
    mix : TYPE
        DESCRIPTION.

    """
    
    edelta = {}; comspace = {}; cna_sigs = {}
    com = {}; comdist = {}; surf_atoms = {} 
    comAu = {}; comPt = {}; hoadjAu = {}; hoadjPt = {} 
    comdistAu = {}; comdistPt = {}; midcomdistPt = {} ; nn = {}
    midcomdistAu = {}; surf_atomsPt = {}; headj = {}; mix = {}
    PtAu = {}; PtOnly = {}; AvgCoPt = {}; GyrationPt = {}; Gyration = {}
    return (edelta, comspace, cna_sigs, com, comdist, 
            surf_atoms, comAu, comPt, hoadjAu, hoadjPt, comdistAu, 
            comdistPt, midcomdistAu, midcomdistPt, surf_atomsPt, 
            headj, mix, nn, PtAu, PtOnly, AvgCoPt, Gyration, GyrationPt)

def Reader(T, Seed, Struts, Sims):
    """
    

    Parameters
    ----------
    Struts : TYPE
        DESCRIPTION.
    Sims : TYPE
        DESCRIPTION.

    Returns
    -------
    init : TYPE
        DESCRIPTION.

    """
    init = Init()
    for Strut in Struts:
        for Sim in Sims:
            try:
                with open(T+Seed+Strut+Sim+'Metadata.csv', 'rb') as infile:
                    Temp = pickle.load(infile)
                    init[0][Strut+Sim] = Temp['edelta'] #t-type: number
                    init[1][Strut+Sim] = Temp['comspace'] #t-type: array
                    init[2][Strut+Sim] = Temp['cna_sigs'] #t-type: number
                    init[3][Strut+Sim] = Temp['com'] #t-type: array 
                    init[4][Strut+Sim] = Temp['comdist'] #t-type: array
                    init[5][Strut+Sim] = Temp['surf_atoms'] #t-type: number
                    init[6][Strut+Sim] = Temp['comAu'] #t-type: array
                    init[7][Strut+Sim] = Temp['comPt'] #t-type: array
                    hoadjAu = Temp['hoadjAu']
                    init[8][Strut+Sim] = np.array([ x for x in hoadjAu ] ) #t-type: list
                    hoadjPt = Temp['hoadjPt']
                    init[9][Strut+Sim] = np.array([ x for x in hoadjPt ] ) #t-type: list
                    init[10][Strut+Sim] = Temp['comdistAu'] #t-type: array
                    init[11][Strut+Sim] = Temp['comdistPt'] #t-type: array
                    init[12][Strut+Sim] = Temp['midcomdistPt'] #t-type: array
                    init[13][Strut+Sim] = Temp['midcomdistAu'] #t-type: array
                    init[14][Strut+Sim] = Temp['surf_atomsPt'] #t-type: number
                    headj = Temp['headj'] #t-type: tuple #######
                    PtAuTemp = []
                    PtOnlyTemp = []
                    for t in range(len(headj)):
                        Temp1 = [ x for x in headj[t][1] if x == 0 ]
                        Temp2 = [ x for x in headj[t][1] if x > 9 ]
                        PtAuTemp.append(len(Temp2))
                        PtOnlyTemp.append(len(Temp1)/55)
                    init[15][Strut+Sim] = PtAuTemp
                    init[16][Strut+Sim] = PtOnlyTemp
                    init[17][Strut+Sim] = Temp['mix'] #t-type: number
                    Spare = []
                    for t in range(len(headj)):
                        c = np.average([ Temp['hoadjPt'][t][i] + headj[t][0][i] for i in range(len(Temp['hoadjPt'][t])) ])
                        Spare.append(c)
                    init[18][Strut+Sim] = Spare #t-type: number
                    init[19][Strut+Sim] = Temp['gyrationPt'] #t-type: number
                    init[20][Strut+Sim] = Temp['gyration'] #t-type: number
                    init[21][Strut+Sim] = headj
                    del(Temp)
                    print(Strut+Sim)
            except Exception as e:
                print(e)
    return init

def clean(data, strut):
    System = {
        'edetla' : np.zeros(len(data[0]), dtype = float),
        'comspace' : np.zeros(len(data[1]), dtype = float),
        '421' : np.zeros(len(data[2]), dtype = float),
        '422' : np.zeros(len(data[2]), dtype = float),
        '555' : np.zeros(len(data[2]), dtype = float),
        'com' : np.zeros(len(data[3]), dtype = object),
        'comdist' : np.zeros(len(data[4]), dtype = object),
        'surf_atoms' : np.zeros(len(data[5]), dtype = float),
        'comAu' : np.zeros(len(data[6]), dtype = object),
        'comPt' : np.zeros(len(data[7]), dtype = object),
        'hoadjAu' : np.zeros(len(data[8]), dtype = object),
        'hoadjPt' : np.zeros(len(data[9]), dtype = object),
        'comdistAu' : np.zeros(len(data[10]), dtype = object),
        'comdistPt' : np.zeros(len(data[11]), dtype = object),
        'midcomdistAu' : np.zeros(len(data[13]), dtype = object),
        'midcomdistPt' : np.zeros(len(data[12]), dtype = object),
        'surf_atomsPt' : np.zeros(len(data[14]), dtype = float),
        'mix' : np.zeros(len(data[16]), dtype = float),
        
        'headj' : np.zeros(len(data[15]), dtype = object),        
        'atombonds' : np.zeros(len(data[15]), dtype = object),
        'deltaatoms' : np.zeros(len(data[15]), dtype = object),
        'deltabonds' : np.zeros(len(data[15]), dtype = object),
        'nnadj' : np.zeros(len(data[15]), dtype = object),
        'nn' : np.zeros(len(data[15]), dtype = object),
        'PtOnly' : np.zeros(len(data[15]), dtype = object),
        'PtAu' : np.zeros(len(data[15]), dtype = object),
        'GyrPt' : np.zeros(len(data[15]), dtype = object),
        'Gyr' : np.zeros(len(data[15]), dtype = object),
        'AvgCoPt' : np.zeros(len(data[15]), dtype = object)
        }
    Keys = data[0].keys()
    print(Keys)
    
    Tempedelta = []; Tempcomspace = []; Temp421 = []; Temp422 = []; Temp555 = []
    Tempcom = []; Tempcomdist = []; Tempsurf_atoms = [] 
    TempcomAu = []; TempcomPt = []; TemphoadjAu = []; TemphoadjPt = [] 
    TempcomdistAu = []; TempcomdistPt = []; TempmidcomdistPt = [] 
    TempmidcomdistAu = []; Tempsurf_atomsPt = []; Tempmix = []
    Tempheadj = []; Tempatombonds = []; Tempdeltaatoms = []; Tempdeltabonds = []
    Tempnnadj = []; Tempnn = []; TempPtOnly = []; TempPtAu = []
    TempGyrPt = []; TempGyr = []; TempAvgCoPt = []
    
    for Key in Keys:
        try:
            Tempedelta.append(data[0][Key])
            Tempcomspace.append(data[1][Key])
            Temp421.append( Collect_CNA( data[2][Key], (4, 2, 1) ) )
            Temp422.append( Collect_CNA( data[2][Key], (4, 2, 2) ) )
            Temp555.append( Collect_CNA( data[2][Key], (5, 5, 5) ) )
            Tempcom.append(data[3][Key])
            Tempcomdist.append(data[4][Key])
            Tempsurf_atoms.append(data[5][Key])
            TempcomAu.append(data[6][Key])
            TempcomPt.append(data[7][Key])
            TemphoadjAu.append(data[8][Key])
            TemphoadjPt.append(data[9][Key])
            TempcomdistAu.append(data[10][Key])
            TempcomdistPt.append(data[11][Key])
            TempmidcomdistAu.append(data[13][Key])
            TempmidcomdistPt.append(data[12][Key])
            Tempsurf_atomsPt.append(data[14][Key])
            TempPtAu.append(data[15][Key])
            TempPtOnly.append(data[16][Key])
            Tempmix.append(data[17][Key])
            TempAvgCoPt.append(data[18][Key])
            TempGyrPt.append(data[19][Key])
            TempGyr.append(data[20][Key])
            HeAdj = hebond(data[21][Key])
            
            Tempheadj.append(HeAdj[0])
            Tempatombonds.append(HeAdj[1])
            Tempdeltaatoms.append(HeAdj[2])
            Tempdeltabonds.append(HeAdj[3]) 
            Tempnnadj.append(HeAdj[4])
            #New_File(Key+'NewMovie.xyz', new_movie=Key+'Quantity_movie.xyz', Quantities = [HeAdj[1], HeAdj[2]])

        except Exception as e:
            print(e)

    System['edetla'] = np.average(Tempedelta, axis = 0)
    System['comspace'] = np.average(Tempcomspace, axis = 0)
    System['421'] = np.average(Temp421, axis = 0)
    System['422'] = np.average(Temp422, axis = 0)
    System['555'] = np.average(Temp555, axis = 0)
    System['com'] = np.average(Tempcom, axis = 0)
    System['comdist'] = np.average(Tempcomdist, axis = 0)
    System['surf_atoms'] = np.average(Tempsurf_atoms, axis = 0)
    System['comAu'] = np.average(TempcomAu, axis = 0)
    System['comPt'] = np.average(TempcomPt, axis = 0)
    System['hoadjAu'] = np.average(TemphoadjAu, axis = 0)
    System['hoadjPt'] = np.average(TemphoadjPt, axis = 0)
    System['comdistAu'] = np.average(TempcomdistAu, axis = 0)
    System['comdistPt'] = np.average(TempcomdistPt, axis = 0)
    System['midcomdistAu'] = np.average(TempmidcomdistAu, axis = 0)
    System['midcomdistPt'] = np.average(TempmidcomdistPt, axis = 0)
    System['surf_atomsPt'] = np.average(Tempsurf_atomsPt, axis = 0)
    System['mix'] = np.average(Tempmix, axis = 0)
    
    System['headj'] = np.average(Tempheadj, axis = 0)
    System['atombonds'] = np.average(Tempatombonds, axis = 0)
    System['deltaatoms'] = np.average(Tempdeltaatoms, axis = 0)
    System['deltabonds'] = np.average(Tempdeltabonds, axis = 0)
    System['nnadj'] = np.average(Tempnnadj, axis = 0)
    System['PtAu'] = np.average(TempPtAu, 0)
    System['PtOnly'] = np.average(TempPtOnly, 0)
    System['nn'] = np.average(Tempnn, 0)
    System['GyrPt'] = np.average(TempGyrPt, 0)
    System['Gyr'] = np.average(TempGyr, 0)
    System['AvgCoPt'] = np.average(TempAvgCoPt, 0)

    return System

class Plot():
    def __init__(self):
        return None
    
    def plotting(self, System, Strut, T):
        self.System = System 
        self.Strut = Strut
        self.T = T
        
        self.Time = np.linspace(0, 500, len(self.System['edetla'])) #time in ns
        
    def cna(self, name = 'CNA_Sigs'):
        
        fig, ax = plt.subplots()
        fig.set_size_inches(8,3)
        ax.plot(self.Time, self.System['421'])
        ax.plot(self.Time, self.System['422'])
        ax2 = ax.twinx()
        x_filt = self.Time[self.System['555'] > 0]
        y_filt = self.System['555'][self.System['555'] > 0]
        ax2.scatter(x_filt, y_filt, color= 'g')
        ax2.set_ylim(0,max(y_filt))
        ax2.set_yticklabels([])
        ax.set_xlabel('Time (ns)', fontsize = 14)
        plt.subplots_adjust(right=0.8)
        scale = max(ax.get_yticks())/max(ax2.get_yticks())
        line_labels = ["(4 2 1)", "(4 2 2)", r"(5 5 5)$\times$%s" %int(scale)]
        fig.legend(labels=line_labels,loc = 'center right', title = 'CNA Sigs')
        plt.savefig(name+self.T+self.Strut+'.png', dpi=400, bbox_inches='tight')
    
    def hebonds(self, name = 'Bonds'):
        fig, ax = plt.subplots()
        fig.set_size_inches(8,3)
        ax.plot(self.Time, 
                [ self.System['headj'][t]/ self.System['headj'][0] for t in range(len(self.Time))], 
                label = '#Hetero Bonds')
        
        ax.set_xlabel('Time (ns)', fontsize = 14)
        ax2 = ax.twinx()
        ax2.plot(self.Time, 
                 [ self.System['surf_atomsPt'][t]/ self.System['surf_atomsPt'][0] for t in range(len(self.Time))], 
                 color='k')
        ax2.set_ylabel('#Surf(t) / #Surf(0)')
        ax.set_ylabel('#He(t) / #He(0)')
        labels = ['Hetero bonds', 'Pt Surface']
        fig.legend(labels = labels, loc = 'center', fontsize = 14)
        plt.savefig(name+self.T+self.Strut+'.png', dpi=400, bbox_inches='tight')
        
    def comtraj(self, name = 'Distance'):
        fig, ax = plt.subplots()
        fig.set_size_inches(8,3)
        A = A = [ distance(self.System['comPt'][t], self.System['comAu'][t]) for t in range(len(self.Time)) ]
        ax.plot(self.Time, max(A) - A)
        ax2 = ax.twinx()
        ax2.plot(self.Time, [ self.System['surf_atomsPt'][t]/self.System['surf_atomsPt'][0] for t in range(len(self.Time))], color='k')
        ax.set_ylabel(r'$\Delta$|CoM(Pt), CoM(Au)|')
        ax2.set_ylabel('#Surf(t) / #Surf(0)')
        ax.set_xlabel('Time (ns)', fontsize = 14)
        labels = ['Distance', 'Pt Surface']
        fig.legend(labels = labels, loc = 'center', fontsize = 14)
        plt.savefig(name+self.T+self.Strut+'.png', dpi=400, bbox_inches='tight')
        
    def midPtDist(self, name = 'MidPtDist'):
        imin = []; imax = []
        for t in range(len(self.Time)):
            s = np.flatnonzero(self.System['midcomdistPt'][t] > 0)
            imin.append(s[0])
            imax.append(s[-1])
        fig, ax = plt.subplots()
        fig.set_size_inches(8,3)
        ax.scatter(self.Time, self.System['comspace'][imin]/10, color= 'r', label = 'Min')
        ax.scatter(self.Time, self.System['comspace'][imax]/10, color= 'g', label = 'Max')
        fig.legend(loc = 'center')
        ax.set_xlabel('Time (ns)', fontsize = 14)
        ax.set_ylabel('Distance (nm)')
        plt.savefig(name+self.T+self.Strut+'.png', dpi=400, bbox_inches='tight')
        
    def PtDist(self, name = 'PtDist'):
        imin = []; imax = []
        for t in range(len(self.Time)):
            s = np.flatnonzero(self.System['comdistPt'][t] > 0)
            imin.append(s[0])
            imax.append(s[-1])
        fig, ax = plt.subplots()
        fig.set_size_inches(8,3)
        ax.scatter(self.Time, self.System['comspace'][imin]/10, color= 'r', label = 'Min')
        ax.scatter(self.Time, self.System['comspace'][imax]/10, color= 'g', label = 'Max')
        fig.legend(loc = 'center')
        ax.set_xlabel('Time (ns)', fontsize = 14)
        ax.set_ylabel('Distance (nm)')
        plt.savefig(name+self.T+self.Strut+'.png', dpi=400, bbox_inches='tight')
               
    def AuExpand(self, name = 'AuVolume'):
        Vol = []; imax = []
        for t in range(len(self.Time)):
            s = np.flatnonzero(self.System['comdistAu'][t] > 0)
            imax.append(s[-1])
            v1 = np.trapz(self.System['comdistAu'][t][:s[-1]], 
                         self.System['comspace'][:s[-1]], dx = 0.05) 
            v2 = np.trapz(self.System['comdistAu'][t][:imax[0]], 
                     self.System['comspace'][:imax[0]], dx = 0.05)
            Vol.append(abs(v1-v2))
        fig,ax = plt.subplots()
        fig.set_size_inches(8,3)
        ax.plot(self.Time, Vol, color = 'g')
        ax.set_ylabel('% Volume')
        ax.set_xlabel('Time (ns)', fontsize = 14)
        plt.savefig(name+self.T+self.Strut+'.png', dpi = 400, bbox_inches='tight')
        
    
    def AuNeigh(self, name = 'AuNN'):
        fig,ax = plt.subplots()
        fig.set_size_inches(10,6)
        
        for a in range(2,9):
            B = []
            for t in range(len(self.Time)):
                temp = [ x for x in self.System['nnadj'][t] if x == a]
                B.append(len(temp))
                
            ax.plot(self.Time,B, label = 'NN $\geq$ %s' %(a))
            
        ax.set_xlabel('Time (ns)')
        ax.set_ylabel('Au neighbours')
        fig.legend(loc = 'center right', fontsize = 13)
        plt.subplots_adjust(right=0.825)
        plt.savefig(name+self.T+self.Strut+'.png', dpi = 400, bbox_inches='tight')
        
    def CompBars(self, name = 'CompBars'):
        fig,ax = plt.subplots()
        fig.set_size_inches(10,6)
        
        for t in range(0,len(self.Time),10):
            bins = np.arange(0, 1.05, 0.05)
            a,b = np.histogram(self.System['PtOnly'][t], bins)
            bin_width = b[1]-b[0]
            bin_cents = [ b[i]+ bin_width for i in range(len(b)-1) ]
            fig,ax = plt.subplots()
            fig.set_size_inches(8,5)
            ax.bar(bin_cents, a, width= 0.05, color = 'k')
            ax.set_xlabel('NN$_{Au}$ / NN$_{Tpt}$', fontsize = 16)
            ax.set_ylabel('Frequency', fontsize = 16)
            ax.text(0.4, 0.5*max(a), 'Time | {:d} ns'.format(int(self.Time[t])), fontsize = 14)
            ax.set_xlim(0,1)
            plt.show()
            plt.savefig(name+self.T+self.Strut+'%s.png'%t, dpi = 400, bbox_inches='tight')
            
    def Panels(self, name = 'MultiPanels'):
        fig,axs = plt.subplots(2,1)
        fig.set_size_inches(14,6)
        
        ax1,ax2 = axs
        ax1.plot(self.Time, self.System['AvgCoPt'], label =r'$\langle$ NN$_{Pt}$ $\rangle$', color = 'k')
        ax1.set_ylabel(r'$\langle$ NN$_{Pt}$ $\rangle$', fontsize = 16)
        ax1.set_xticklabels([])
        ax12 = ax1.twinx()
        
        ax12.plot(self.Time, [self.System['PtOnly'][t]/55 for t in range(len(self.Time))], label = 'Pt | No Au NN', color='r')
        
        ax12.plot(self.Time, self.System['PtAu'], label = r'NN(Pt)|$_{Au>10}$', color = 'g')
        ax12.set_ylabel('#Pt', fontsize = 16)
        
        ax2.plot(self.Time, self.System['mix'], label = r'$\mu$', color = 'y')
        ax2.set_ylabel('Mixing', fontsize = 16)
        
        ax22 = ax2.twinx()
        ax22.plot(self.Time, [ self.System['GyrPt'][t] / self.System['GyrPt'][0] for t in range(len(self.Time))], label = 'Pt Radius of gyration', color = 'c')

        ax22.plot(self.Time, [ self.System['Gyr'][t] / self.System['Gyr'][0] for t in range(len(self.Time))], label = 'Radius of gyration', color = 'm')
        ax22.set_ylabel('RoG(t) / RoG(0)', fontsize = 16)
        
        ax2.set_xlabel('Time (ns)', fontsize = 16)
        lines_labels = [ax.get_legend_handles_labels() for ax in fig.axes]
        lines, labels = [sum(lol, []) for lol in zip(*lines_labels)]
        fig.legend(lines,labels,loc = 'upper center', fontsize = 14, ncol = 6)
        plt.savefig(name+self.T+self.Strut+'.png', dpi = 400, bbox_inches='tight')
    