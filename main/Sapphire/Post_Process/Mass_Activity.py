from ase import Atoms
from ase.io import read
import numpy as np
from matplotlib import pyplot as plt
import time
from scipy.interpolate import interp1d
import seaborn as sns
import pickle
import numpy as np
from itertools import groupby
from collections import namedtuple
from scipy.ndimage import gaussian_filter1d
import math

"""CONSTANTS TO BE UPDATED AND INSERTED MANUALLY BY THE USER !!!"""

def beta(T):
    return 1/((8.6173303E-5)*(T))

sigma=2 #how much you want to smoothen out your functions?
applied_V=1.1 

U_bins = 1299 #binning of the voltages (v)
 
#Change the temperature at whihc you want the catalytic activities to be performed at.
r=1.28*10**(-9) #atomic radius of your atoms (in cm)?
mass_cu = 1.0552e-19 #mass of your atoms (in mg)?
C = 12.56 # constant in equation, to be calculated from initial conditions

def heatmap(agcn):
    occ_long=[]
    for i in range(len(agcn)):
        (n, bins, patches)=plt.hist(agcn[i], bins=90, density=True)
        occ_long.append(n)
    occ_long=np.asarray(occ_long)
    corrected_occ=np.zeros((90, len(agcn)))
    for i in range(90):
        corrected_occ[-i-1]=np.asarray(occ_long[:, i])
    sns.heatmap(corrected_occ, vmax=1, xticklabels=False, yticklabels=False, cbar=False)
    #plt.savefig(str(folder)+’occurrency.png’, bbox_inches=’tight’, dpi=400)
    #plt.show()
    #plt.close()
    return corrected_occ



def catalytic_analysis (filename, agcn):
    corrected_occ = heatmap(agcn)
    traj = read(filename, index = ':')
    current=[]; mass_activity=[]; y=[]
    for j in range(0, len(agcn)):
        surf_area=0
        for i in range(len(agcn[j])):
            once = 4*np.pi*r**2*(1-agcn[j][i]/12)
            surf_area = surf_area + once
        y.append(surf_area)
        site=np.zeros((U_bins))
        for h in range(U_bins):
            spec=0
            for m in range(len(corrected_occ)):
                if m<31:
                    sitecurrent = C*np.exp(((0.162 * m/10 - 1.11)-(h*0.001))*beta)*m/10*corrected_occ[m][j]/len(agcn[j])
                elif 31<=m<81:
                    sitecurrent = C*np.exp(((-0.067 * m/10 - 0.416)-(h*0.001))*beta)*m/10*corrected_occ[m][j]/len(agcn[j])
                else:
                    sitecurrent = C*np.exp(((-0.222 * m/10 + 0.849)-(h*0.001))*beta)*m/10*corrected_occ[m][j]/len(agcn[j])
                spec=spec+sitecurrent
            site[h]=spec
        current.append(site)
        mass_activity.append(-site[int(1299-applied_V*1000)]*surf_area/mass_NP)
            
    """Plotting the current density at different applied potentials for different
    time steps."""
    potentials = np.linspace(-1.299, 0, U_bins)
    plt.plot(potentials, current[int(len(traj)-1)], color='orange',lw=3, label='final')
    plt.plot(potentials, current[int(len(traj)/2)-1], color='purple', label='middle')
    plt.plot(potentials, current[0], color='r', label='initial')
    
    plt.xlabel('V vs RHE')
    plt.ylabel('j (mA/cm^2)')
    plt.savefig('current densities.png', bbox_inches='tight', dpi=200)
    plt.show()
    plt.close()
    """Mass activity plots at desired applied_V."""
    plt.xlabel('Time (ns)')
    plt.ylabel('MA (mA/mg)')
    plt.plot(mass_activity, linestyle='dashed', color='k')
    plt.plot(gaussian_filter1d(mass_activity, sigma), linestyle='solid', color='k')
    plt.savefig('mass_activity.png', bbox_inches='tight', dpi=200)
    plt.show()
    plt.close()