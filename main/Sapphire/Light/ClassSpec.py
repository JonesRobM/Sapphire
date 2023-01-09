import numpy as np
import matplotlib.pyplot as plt
import Epsilon_ExpClass as Epsilon

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

    
    """This is a conceptual class representation of a simple BLE device
    (GATT Server). It is essentially an extended combination of the
    :class:`bluepy.btle.Peripheral` and :class:`bluepy.btle.ScanEntry` classes

    :param client: A handle to the :class:`simpleble.SimpleBleClient` client
        object that detected the device
    :type client: class:`simpleb le.SimpleBleClient`
    :param addr: Device MAC address, defaults to None
    :type addr: str, optional
    :param addrType: Device address type - one of ADDR_TYPE_PUBLIC or
        ADDR_TYPE_RANDOM, defaults to ADDR_TYPE_PUBLIC
    :type addrType: str, optional
    :param iface: Bluetooth interface number (0 = /dev/hci0) used for the
        connection, defaults to 0
    :type iface: int, optional
    :param data: A list of tuples (adtype, description, value) containing the
        AD type code, human-readable description and value for all available
        advertising data items, defaults to None
    :type data: list, optional
    :param rssi: Received Signal Strength Indication for the last received
        broadcast from the device. This is an integer value measured in dB,
        where 0 dB is the maximum (theoretical) signal strength, and more
        negative numbers indicate a weaker signal, defaults to 0
    :type rssi: int, optional
    :param connectable: `True` if the device supports connections, and `False`
        otherwise (typically used for advertising ‘beacons’).,
        defaults to `False`
    :type connectable: bool, optional
    :param updateCount: Integer count of the number of advertising packets
        received from the device so far, defaults to 0
    :type updateCount: int, optional
    """

    def __init__(self, Structure, Frame = -1, Step = 2.88, wl = np.linspace(300, 800, 101),
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
