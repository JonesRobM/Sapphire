import numpy as np
import matplotlib.pyplot as plt
import Epsilon_ExpClass as Epsilon
import os

from pyGDM2 import structures
from pyGDM2 import fields
from pyGDM2 import propagators

from pyGDM2 import core
from pyGDM2 import visu
from pyGDM2 import tools
from pyGDM2 import linear

from ase.io import read

## --------------- Setup structure

class Spectrum(object):
    
    def __init__(self, Structure, Step = 2.88, wl = np.linspace(300, 800, 101),
                 wavetype = 'planewave', angles = [0], n1 = 1, n2 = 1, scale = 1):
        
        self.wl = wl
        self.wave = wavetype
        self.anglekwargs = dict(theta = angles)
        self.scale = scale
        self.Step = Step*self.scale
        self.name = Structure

        Strut = read(Structure)
        ele = Strut.get_chemical_symbols()
        Ele = list(set(ele))
        Strut.positions -= Strut.get_center_of_mass()
        xyz = Strut.positions*self.scale
        
        if len(Ele)>1:
            self.g = [ xyz[i] for i, atom in enumerate(ele) if atom == Ele[0] ]
            self.m = len(self.g)*[getattr(Epsilon, Ele[0])()]
            for i in range(1, len(Ele)):
                gi = [ xyz[j] for j, atom in enumerate(ele) if atom == Ele[i] ]
                self.g = np.concatenate([self.g, gi])
                self.m += len(gi)*[getattr(Epsilon, Ele[i])()]
                
        else:   
            self.g = xyz
            self.m = len(self.g)*[getattr(Epsilon, Ele[0])()]
            
            
            
        ## structure instance
        self.Struct = structures.struct(Step, self.g, self.m)

        while len(self.Struct.geometry) < len(self.g):
            self.Optimise()
            
        
        ## incident field
        field_generator = getattr(fields, self.wave)
        
        efield = fields.efield(field_generator, wavelengths=self.wl, 
                               kwargs=self.anglekwargs)
        
        ## environment
        dyads = propagators.DyadsQuasistatic123(n1, n2)

        ## simulation initialization
        self.sim = core.simulation(self.Struct, efield, dyads)
        
        
        ## --------------- run scatter simulation
        self.efield = core.scatter(self.sim, verbose=False)
        visu.structure(self.sim, show=0)
        visu.vectorfield_by_fieldindex(self.sim, 0)
        
    def Optimise(self):
        self.Step -= 0.01*self.scale
        self.Struct = structures.struct(self.Step, self.g, self.m, verbose = False)


    def Plot_Spectrum_En(self):

        field_kwargs = tools.get_possible_field_params_spectra(self.sim)
        for i, conf in enumerate(field_kwargs):
            print("config", i, ":", conf)
        
      
            wl, spectrum = tools.calculate_spectrum(self.sim,
                                field_kwargs[i], linear.extinct)
            np.savetxt('Spectrum_'+self.name[:-4]+'.txt', np.column_stack((wl, spectrum.T[0], spectrum.T[1], spectrum.T[2])))
            
            ## --- linear.extinct returns 3-tuples, "spectrum" therefore consists
            ## --- of an array of 3-tuples, corresponding to extinction,
            ## --- scattering and absorption.
            
            plt.plot(1.2415285500000001e3/wl, spectrum.T[0], 'g-', label='ext.')
            plt.plot(1.2415285500000001e3/wl, spectrum.T[1], 'b-', label='scat.')
            plt.plot(1.2415285500000001e3/wl, spectrum.T[2], 'r-', label='abs.')
            
            plt.xlabel("Energy (ev)")
            plt.ylabel("cross section (nm$^2$)")
            plt.legend(loc='best', fontsize=8)
            plt.savefig('Spectrum_'+self.name[:-4]+'.png', dpi = 400, bbox_inches = 'tight')
            plt.show()
            
    def Plot_Spectrum_Wl(self):

        field_kwargs = tools.get_possible_field_params_spectra(self.sim)
        for i, conf in enumerate(field_kwargs):
            print("config", i, ":", conf)
        
      
            wl, spectrum = tools.calculate_spectrum(self.sim,
                                field_kwargs[i], linear.extinct)
            
            ## --- linear.extinct returns 3-tuples, "spectrum" therefore consists
            ## --- of an array of 3-tuples, corresponding to extinction,
            ## --- scattering and absorption.
            
            plt.plot(wl, spectrum.T[0], 'g-', label='ext.')
            plt.plot(wl, spectrum.T[1], 'b-', label='scat.')
            plt.plot(wl, spectrum.T[2], 'r-', label='abs.')
            
            plt.xlabel("Wavelength nm")
            plt.ylabel("cross section (nm$^2$)")
            plt.legend(loc='best', fontsize=8)
            plt.savefig('Spectrum_'+self.name[:-4]+'.png', dpi = 400, bbox_inches = 'tight')
            plt.show()

    def Plot_EField(self):
        
        r_probe = tools.generate_NF_map(2*min(self.g[:,0]),2*max(self.g[:,0]),101, 
                                        2*min(self.g[:,1]),2*max(self.g[:,1]),101, 
                                        Z0=self.g.T[2].max()+2*self.Step)
        Es, Et, Bs, Bt = linear.nearfield(self.sim, 0, r_probe)
        np.savetxt('Probe_'+self.name[:-4]+'.txt', np.column_stack((Es, Et, Bs, Bt)))
        visu.structure(self.g, show = 0, alpha = 0.5)
        im = visu.vectorfield_color(Es, tit='intenstiy of scattered field outside', show=0)
        plt.colorbar(im, label=r'$|E_s|^2 / |E_0|^2$')
        plt.xlabel("x (nm)")
        plt.ylabel("y (nm)")
        plt.savefig('E-Field'+self.name[:-4]+'.png', dpi = 400, bbox_inches = 'tight')
        plt.show()
  
if __name__ =='__main__':      
      
    for Strut in os.listdir():
        if Strut.endswith('.xyz'):
            Spec = Spectrum(Structure = Strut, scale = 0.1)
            Spec.Plot_EField()
            Spec.Plot_Spectrum_En()
            Spec.Plot_Spectrum_Wl()
    
    import matplotlib.colors as mcol
    import matplotlib.cm as cm
    
    # Make a user-defined colormap.
    cm1 = mcol.LinearSegmentedColormap.from_list("MyCmapName",["r","b"])
    cnorm = mcol.Normalize(vmin=91,vmax=112)
    cpick = cm.ScalarMappable(norm=cnorm,cmap=cm1)
    cpick.set_array([])
    
    os.chdir('../')
    fig, axs = plt.subplots(1,3, True)
    fig.set_size_inches(20/2.5, 4/2.5)
    fig.subplots_adjust(wspace=0, hspace=0)
    ax1,ax2,ax3= axs.flatten()
    os.chdir('Au32/')
    Temp = {}
    for Strut in os.listdir():
        if (Strut.endswith('.txt') and Strut.startswith('Spectrum')):
            Temp[Strut[-7:-4]] = np.loadtxt(Strut)
        New = dict(sorted(Temp.items()))
    for i, item in enumerate(New):
        gamma = str("{:.2f}".format(float(item)/100 - 1))
        if ('1' in gamma or '3' in gamma or '6' in gamma or '9' in gamma):
            ax1.plot(1.2415285500000001e3/New[item][:,0], New[item][:,1],
            color=cpick.to_rgba(float(item)))
    ax1.plot(1.2415285500000001e3/New['100'][:,0], New['100'][:,1],
    label = '0.00', color = 'g')
    os.chdir('../Cu32/')
    Temp = {}
    for Strut in os.listdir():
        if (Strut.endswith('.txt') and Strut.startswith('Spectrum')):
            Temp[Strut[-7:-4]] = np.loadtxt(Strut)
        New = dict(sorted(Temp.items()))
    for i, item in enumerate(New):
        gamma = str("{:.2f}".format(float(item)/100 - 1))
        if ('1' in gamma or '3' in gamma or '6' in gamma or '9' in gamma):
            ax2.plot(1.2415285500000001e3/New[item][:,0], New[item][:,1],
            color=cpick.to_rgba(float(item)))
    ax2.plot(1.2415285500000001e3/New['100'][:,0], New['100'][:,1], color = 'g')
    os.chdir('../Ag32/')
    Temp = {}
    for Strut in os.listdir():
        if (Strut.endswith('.txt') and Strut.startswith('Spectrum')):
            Temp[Strut[-7:-4]] = np.loadtxt(Strut)
        New = dict(sorted(Temp.items()))
    for i, item in enumerate(New):
        gamma = str("{:.2f}".format(float(item)/100 - 1))
        if ('1' in gamma or '3' in gamma or '6' in gamma or '9' in gamma):
            ax3.plot(1.2415285500000001e3/New[item][:,0], New[item][:,1],
            label = (str("{:.2f}".format(float(item)/100 - 1))),
            color=cpick.to_rgba(float(item)))
    ax3.plot(1.2415285500000001e3/New['100'][:,0], New['100'][:,1], color = 'g')
    ax1.set_yticklabels([])
    ax2.set_yticklabels([])
    ax3.set_yticklabels([])
    
    #ax1.set_ylabel('Ext. (a.u.)', fontsize = 10)
    fig.text(0.1, 0.33, 'Ext. (a.u.)', fontsize = 10, rotation = 90)
    fig.text(0.45, -0.1, 'Energy (eV)', fontsize = 10)
    fig.text(0.225, 0.75, 'Au', fontsize = 10)
    fig.text(0.475, 0.75, 'Cu', fontsize = 10)
    fig.text(0.725, 0.75,'Ag', fontsize = 10)
    fig.legend(bbox_to_anchor=(0.125, 0.9, 0.775, .12), loc = 3,
    ncol = 10, title = r'$\gamma$', fontsize = 8, mode="expand", borderaxespad=0.)
    plt.savefig('../ClassicalSpectra_Smol.png', dpi = 600, bbox_inches = 'tight')
