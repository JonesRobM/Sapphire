import os
import sys
import warnings
from inspect import getmembers, isfunction
from ase.io import read

none_template = '\nSystem property "%s" is bad. Typically, this is because the ' \
                'required information has not been provied by the user or is given incorrectly.\n' \
                'Reverting to System default "%s".\n'
"""
CNA_Pattern_Settings = {
    'npz_dir': 'CNA_npz/',  # folder to save the npz files in
    'new_xyz_dir': 'CNA_XYZs/',
    'APPEND_DICTIONARY': False,
    'FROM_MEMORY': False,
    'BULK_MASTERKEY': True,
    'SAVING_XYZ': True,
    'PRINTING_PATTERNS': True
}
"""

def ensure_dir(base_dir='', file_path=''):
    """

    Robert:

        A simple script to verify the existence of a directory
        given the path to it. If it does not exist, will create it.
        For this module we shall be looking to create, if requested, directories
        for saved pattern data to be passed to a classification module.

    """

    directory = base_dir + file_path
    if not os.path.exists(directory):

        os.makedirs(directory)

class _Clean_Pattern(object):

    """This is a cleanup class which serves 

    :param client: A handle to the :class:`simpleble.SimpleBleClient` client
        object that detected the device
    :type client: class:`simpleble.SimpleBleClient`
    :param addr: Device MAC address, defaults to None
    :type addr: str, optional
    """

    def __init__(self, System, Pattern_Input):
        
        """Constructor method
        """

        self.System = System

        self.Pattern_Input = Pattern_Input

        self.file = 'Sapphire_Info.txt' #Note that this file already exists assuming prescribed use of Sapphire

        #These settings are suggested arguments to use under most circumstances.
        #Tutorials will exist on how and when to use each argument, and what it does.
        self.Default = {
            'npz_dir': 'CNA_npz/',  # folder to save the npz files in
            'APPEND_DICTIONARY': False,
            'FROM_MEMORY': False,
            'BULK_MASTERKEY': True,
            'SAVING_XYZ': True,
            'PRINTING_PATTERNS': True
        }
        self.Keys = list(self.Default.keys())

        #Calls the alphebetised list of functions to pass through.
        #This is why each REQUIRED argument begins with an ordered captial letter.
        #We then proceed to execute each of these functions with the getattr builtin

        self.FunkList = [o for o in getmembers(_Clean_Pattern) if isfunction(o[1])]

        self.Funks = [x[0] for x in self.FunkList if not x[0].startswith('_')]
        for x in self.Funks:
            getattr(self, x)()
            

    def Anpz_dir(self):

        """Returns a list containing :class:`bluepy.btle.Characteristic`
        objects for the peripheral. If no arguments are given, will return all
        characteristics. If startHnd and/or endHnd are given, the list is
        restricted to characteristics whose handles are within the given range.
        """
        
        def _no_base():
            """Print info about device
            """  
            with open(self.System['base_dir']+self.file, 'a') as warn:
                warn.write(none_template % (self.Pattern_Input['npz_dir'], self.Default['npz_dir']))
            self.Pattern_Input['npz_dir'] = self.Defaults['npz_dir']
            ensure_dir(self.System['base_dir'], self.Pattern_Input['npz_dir'])

        try:
            self.Pattern_Input['npz_dir']
            if type(self.Pattern_Input['npz_dir']) is not str:
                self._no_base()

            else:
                if self.Pattern_Input['npz_dir'] == self.Defaults['npz_dir']:
                    ensure_dir(self.System['base_dir'], self.Pattern_Input['npz_dir'])
                else:
                    if not os.path.isdir(self.Pattern_Input['base_dir']):
                        _no_base()
        except KeyError:
            _no_base()
        with open(self.System['base_dir']+self.file, "a") as f:
            f.write("\nInitialising...\n")

    def CAPPEND_DICTIONARY(self):
        
        """Returns a list containing :class:`bluepy.btle.Characteristic`
        objects for the peripheral. If no arguments are given, will return all
        characteristics. If startHnd and/or endHnd are given, the list is
        restricted to characteristics whose handles are within the given range.
        """
        
        def _no_file():
            
            """Returns a list containing :class:`bluepy.btle.Characteristic`
            objects for the peripheral. If no arguments are given, will return all
            characteristics. If startHnd and/or endHnd are given, the list is
            restricted to characteristics whose handles are within the given range.
            """
            
            self.System['energy_file_name'] = self.Default['energy_file_name']
            with open(self.System['base_dir']+self.file, 'a') as warn:
                warn.write("\nNo energy file can be found at the specified location.\n'%s'\n"
                           "Please check your local directories and re-write your input file if you want energetic analysis.\n"
                           % (self.System['base_dir']))

        try:
            if type(self.System['energy_file_name']) is not str:
                self.System['energy_file_name'] = self.Default['energy_file_name']
                warnings.warn(none_template % ('energy_file_name', self.Default['energy_file_name']))
                with open(self.System['base_dir']+self.file, 'a') as warn:
                    warn.write(none_template % ('energy_file_name', self.Default['energy_file_name']))
                _no_file()
            else:
                if not os.path.isfile(self.System['base_dir']+self.System['energy_file_name']):
                    _no_file()

        except Exception as e:
            _no_file()


    def DFROM_MEMORY(self):
        
        """Returns a list containing :class:`bluepy.btle.Characteristic`
        objects for the peripheral. If no arguments are given, will return all
        characteristics. If startHnd and/or endHnd are given, the list is
        restricted to characteristics whose handles are within the given range.
        """

        try:

            if type(self.System['extend_xyz']) is not list:

                self.System['extend_xyz'] = self.Default['extend_xyz']
                warnings.warn(none_template % ('extend_xyz', self.Default['extend_xyz']))
                with open(self.System['base_dir']+self.file, 'a') as warn:
                    warn.write(none_template % ('extend_xyz', self.Default['extend_xyz']))
            else:

                with open(self.System['base_dir']+self.file, 'a') as warn:
                    warn.write("Will attempt to write the following quantities into an extended xyz file:\n")
                    for x in self.System['extend_xyz']:
                        warn.write("%s\n" % x)

        except KeyError:

            self.System['extend_xyz'] = self.Default['extend_xyz']
            warnings.warn(none_template % ('extend_xyz', self.Default['extend_xyz']))
            with open(self.System['base_dir']+self.file, 'a') as warn:
                warn.write(none_template % ('extend_xyz', self.Default['extend_xyz']))

    def _no_homo(self):
        
        """Returns a list containing :class:`bluepy.btle.Characteristic`
        objects for the peripheral. If no arguments are given, will return all
        characteristics. If startHnd and/or endHnd are given, the list is
        restricted to characteristics whose handles are within the given range.
        """

        with open(self.System['base_dir']+self.file, 'a') as warn:
            warn.write("\nNo specie-specific properties for homo species will be calculated in this run.\n")
        self.System['Homo'] = self.Default['Homo']

    def EBULK_MASTERKEY(self):
        
        """Returns a list containing :class:`bluepy.btle.Characteristic`
        objects for the peripheral. If no arguments are given, will return all
        characteristics. If startHnd and/or endHnd are given, the list is
        restricted to characteristics whose handles are within the given range.
        """

        with open(self.System['base_dir']+self.file, "a") as f:
            f.write("\nChecking user input for calculating homo properties in this run.\n")

        try:
            self.System['Homo']
            if self.System['Homo'] is None:
                self._no_homo()

            elif type(self.System['Homo']) is list:
                Temp = read(
                    self.System['base_dir']+self.System['movie_file_name'],
                    index=0).get_chemical_symbols()
            used = set()
            Species = [x for x in Temp
                       if x not in used and (used.add(x) or True)]
            Temp = []
            for x in self.System['Homo']:
                if x not in Species:
                    with open(self.System['base_dir']+self.file, "a") as f:
                        f.write("\nChemical specie %s not present in the trajectory."
                                "Consequently, this shall be discarded from Homo.\n" % x)
                else:
                    Temp.append(x)
            with open(self.System['base_dir']+self.file, "a") as f:
                f.write("\nSpecies being considered are:\n"+'\t'.join(str(x) for x in Temp))
            self.System['Homo'] = Temp
        except Exception as e:
            self._no_homo()

    def _no_hetero(self):
        
        """Returns a list containing :class:`bluepy.btle.Characteristic`
        objects for the peripheral. If no arguments are given, will return all
        characteristics. If startHnd and/or endHnd are given, the list is
        restricted to characteristics whose handles are within the given range.
        """
        
        with open(self.System['base_dir']+self.file, 'a') as warn:
            warn.write("\nNo specie-specific properties for homo species will be calculated in this run.\n")
        self.System['Hetero'] = self.Default['Hetero']

    def GSAVING_XYZ(self):
        
        """Returns a list containing :class:`bluepy.btle.Characteristic`
        objects for the peripheral. If no arguments are given, will return all
        characteristics. If startHnd and/or endHnd are given, the list is
        restricted to characteristics whose handles are within the given range.
        """

        with open(self.System['base_dir']+self.file, "a") as f:
            f.write("\nChecking user input for calculating homo properties in this run.\n")

        try:
            self.System['Hetero']

            if self.System['Hetero'] is None:
                self._no_hetero()
        except KeyError:
            self._no_hetero()

    def IPRINTING_PATTERNS(self):
        
        """Returns a list containing :class:`bluepy.btle.Characteristic`
        objects for the peripheral. If no arguments are given, will return all
        characteristics. If startHnd and/or endHnd are given, the list is
        restricted to characteristics whose handles are within the given range.
        """
        
        try:
            self.System['Start']
            if type(self.System['Start']) is not int or self.System['Start'] < 0:
                self.System['Start'] = 0
                warnings.warn(none_template % ('Start', self.Default['Start']))
                with open(self.System['base_dir']+self.file, 'a') as warn:
                    warn.write(none_template % ('Start', self.Default['Start']))

            else:
                with open(self.System['base_dir']+self.file, 'a') as file:
                    file.write("\nInitial frame has been set to %s.\n" % self.System['Start'])

        except KeyError:
            self.System['Start'] = 0
            warnings.warn(none_template % ('Start', self.Default['Start']))
            with open(self.System['base_dir']+self.file, 'a') as warn:
                warn.write(none_template % ('Start', self.Default['Start']))
