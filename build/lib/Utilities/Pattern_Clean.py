import os
import sys
import warnings
from inspect import getmembers, isfunction
from ase.io import read

none_template = '\nSystem property "%s" is bad. Typically, this is because the ' \
                'required information has not been provied by the user or is given incorrectly.\n' \
                'Reverting to System default "%s".\n'

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

    """This is a cleanup class which serves to preprocess user requests regarding CNAPs
    Defaults assume that if a user is requesting CNAPs, they wish for these only.
    Setting flags otherwise will signal to Sapphire to call the SVC method found
    in the CNA/Classification.py module.

    :param System: The system information for the data being analysed.
        Contains details on base directories to write log files to.
    :type System: dict, always passed by a previous cleanup module in
        Utilities/System_Clean.py
    
    :param Pattern_Input: Device MAC address, defaults to default arguments
        {
            'npz_dir': 'CNA_npz/',  # folder to save the npz files in
            'APPEND_DICTIONARY': False,
            'FROM_MEMORY': False,
            'BULK_MASTERKEY': True,
            'SAVING_XYZ': True,
            'PRINTING_PATTERNS': True
        }
    :type Pattern_Input: dict
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
            'PRINTING_PATTERNS': False
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

        """Parses and cleans the user input relating the the CNAP directory.
        """
        
        def _no_base():
            """Provides a warning that an input location is inappropriate and
                creates a suitable directory for patterns to be saved into.
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
                ensure_dir(self.System['base_dir'], self.Pattern_Input['npz_dir'])
                with open(self.System['base_dir']+self.file, 'a') as warn:
                    warn.write("Directory for storing patterns for this run has been set to ' %s '.\n" (self.Pattern_Input['npz_dir']))
                
        except Exception:
            self._no_base()


    def CAPPEND_DICTIONARY(self):
        
        """Parses and cleans user input on if an existing Pattern Dictionary should
            be appended to. In principle, this will exist within the user's home Sapphire
            directory. However, it is strongly advised that users do not enable this
            unless they know what they are doing!!!
        """

        try:
            if type(self.Pattern_Input['APPEND_DICTIONARY']) is not bool:
                warnings.warn(none_template % ('APPEND_DICTIONARY', self.Default['APPEND_DICTIONARY']))
                with open(self.System['base_dir']+self.file, 'a') as warn:
                    warn.write(none_template % ('APPEND_DICTIONARY', self.Default['APPEND_DICTIONARY']))
                    
                self.Pattern_Input['APPEND_DICTIONARY'] = self.Default['APPEND_DICTIONARY']
        except Exception:
            self.Pattern_Input['APPEND_DICTIONARY'] = self.Default['APPEND_DICTIONARY']
            with open(self.System['base_dir']+self.file, 'a') as warn:
                warn.write("Argument for appending to a known pattern dictionary has been set to ' %s '.\n"  % (self.Pattern_Input['APPEND_DICTIONARY']))


    def DFROM_MEMORY(self):
        
        """
        """

        try:

            if type(self.System['FROM_MEMORY']) is not bool:

                self.System['FROM_MEMORY'] = self.Default['FROM_MEMORY']
                warnings.warn(none_template % ('FROM_MEMORY', self.Default['FROM_MEMORY']))
                with open(self.System['base_dir']+self.file, 'a') as warn:
                    warn.write(none_template % ('FROM_MEMORY', self.Default['FROM_MEMORY']))
            else:

                with open(self.System['base_dir']+self.file, 'a') as warn:
                    warn.write("Argument for reading a known pattern dictionary has been set to ' %s '.\n" (self.Pattern_Input['FROM_MEMORY']))


        except Exception:

            self.System['FROM_MEMORY'] = self.Default['FROM_MEMORY']
            warnings.warn(none_template % ('FROM_MEMORY', self.Default['FROM_MEMORY']))
            with open(self.System['base_dir']+self.file, 'a') as warn:
                warn.write("Argument for reading a known pattern dictionary has been set to ' %s '.\n" (self.Pattern_Input['FROM_MEMORY']))

    def EBULK_MASTERKEY(self):
        
        """Returns a list containing :class:`bluepy.btle.Characteristic`
        objects for the peripheral. If no arguments are given, will return all
        characteristics. If startHnd and/or endHnd are given, the list is
        restricted to characteristics whose handles are within the given range.
        """
        
        try:

            if type(self.System['BULK_MASTERKEY']) is not bool:

                self.System['BULK_MASTERKEY'] = self.Default['BULK_MASTERKEY']
                warnings.warn(none_template % ('BULK_MASTERKEY', self.Default['BULK_MASTERKEY']))
                with open(self.System['base_dir']+self.file, 'a') as warn:
                    warn.write(none_template % ('BULK_MASTERKEY', self.Default['BULK_MASTERKEY']))
            else:

                with open(self.System['base_dir']+self.file, 'a') as warn:
                    warn.write("Argument for reading a known pattern dictionary has been set to ' %s '.\n" (self.Pattern_Input['BULK_MASTERKEY']))


        except Exception:

            self.System['BULK_MASTERKEY'] = self.Default['BULK_MASTERKEY']
            warnings.warn(none_template % ('BULK_MASTERKEY', self.Default['BULK_MASTERKEY']))
            with open(self.System['base_dir']+self.file, 'a') as warn:
                warn.write("Argument for reading a known pattern dictionary has been set to ' %s '.\n" (self.Pattern_Input['BULK_MASTERKEY']))


    def IPRINTING_PATTERNS(self):
        
        """Returns a list containing :class:`bluepy.btle.Characteristic`
        objects for the peripheral. If no arguments are given, will return all
        characteristics. If startHnd and/or endHnd are given, the list is
        restricted to characteristics whose handles are within the given range.
        """
        
        try:

            if type(self.System['PRINTING_PATTERNS']) is not bool:

                self.System['PRINTING_PATTERNS'] = self.Default['PRINTING_PATTERNS']
                warnings.warn(none_template % ('PRINTING_PATTERNS', self.Default['PRINTING_PATTERNS']))
                with open(self.System['base_dir']+self.file, 'a') as warn:
                    warn.write(none_template % ('PRINTING_PATTERNS', self.Default['PRINTING_PATTERNS']))
            else:

                with open(self.System['base_dir']+self.file, 'a') as warn:
                    warn.write("Argument for reading a known pattern dictionary has been set to ' %s '.\n" (self.Pattern_Input['PRINTING_PATTERNS']))


        except Exception:

            self.System['PRINTING_PATTERNS'] = self.Default['PRINTING_PATTERNS']
            warnings.warn(none_template % ('PRINTING_PATTERNS', self.Default['PRINTING_PATTERNS']))
            with open(self.System['base_dir']+self.file, 'a') as warn:
                warn.write("Argument for reading a known pattern dictionary has been set to ' %s '.\n" (self.Pattern_Input['PRINTING_PATTERNS']))
