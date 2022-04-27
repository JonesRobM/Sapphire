import os
import sys
import warnings
from inspect import getmembers, isfunction
from ase.io import read

unsupported_template = '\nProperty "%s" not available. Please verify which features' \
    'in Sapphire are supported first by calling\n' \
    '\nfrom Utilities.Supported import Supported\n' \
    'print(Supported().Full(), Supported().Homo(), Supported().Hetero())\n'

none_template = '\nSystem property "%s" is bad. Typically, this is because the ' \
                'required information has not been provied by the user or is given incorrectly.\n' \
                'Reverting to System default "%s".\n'


class _Clean_System(object):

    def __init__(self, System={}):
        self.System = System
        self.file = 'Sapphire_Info.txt'
        self.Default = {
            'base_dir': '',
            'movie_file_name': 'movie.xyz',
            'energy_file_name': None,
            'extend_xyz': None,
            'Homo': None,
            'Hetero': None,
            'Start': 0, 'End': None, 'Step': 1, 'Skip': 50,
            'UniformPDF': False, 'Band': 0.05
        }
        self.Keys = list(self.Default.keys())

        self.FunkList = [o for o in getmembers(_Clean_System) if isfunction(o[1])]

        self.Funks = [x[0] for x in self.FunkList if not x[0].startswith('_')]
        for x in self.Funks:
            getattr(self, x)()

    def Abase_dir(self):

        def _no_base():
            self.System['base_dir'] = ''
            with open(self.System['base_dir']+self.file, 'a') as warn:
                warn.write(none_template % (self.System['base_dir'], self.Default['base_dir']))

        try:
            self.System['base_dir']
            if type(self.System['base_dir']) is not str:
                self._no_base()

            else:
                if self.System['base_dir'] == '':
                    pass
                else:
                    if not os.path.isdir(self.System['base_dir']):
                        _no_base()
        except KeyError:
            _no_base()
        with open(self.System['base_dir']+self.file, "a") as f:
            f.write("\nInitialising...\n")

    def Bmovie_file_name(self):

        def _exit():
            try:
                if not os.path.isfile(self.System['base_dir']+self.System['movie_file_name']):
                    with open(self.System['base_dir']+self.file, 'a') as warn:
                        warn.write("\nNo trajectory file can be found at the specified location.\n"
                                   "Please check your local directories and re-write your input file.\n"
                                   "Sapphire will now terminate.\n")
                    raise SystemExit("No trajectory found at '%s'.\n" % (
                        self.System['base_dir']+self.System['movie_file_name']))
                    _exit()
            except Exception as e:
                sys.exit('\nCannot find this file.\nExiting now due to error rasied as:\n.%s' % e)

        try:
            _exit()
            if type(self.System['movie_file_name']) is not str:
                self.System['movie_file_name'] = self.Default['movie_file_name']
                _exit()
                warnings.warn(none_template % ('movie_file_name', self.Default['movie_file_name']))
                with open(self.System['movie_file_name']+self.file, 'a') as warn:
                    warn.write(none_template % ('movie_file_name', self.Default['movie_file_name']))
                _exit()

            else:
                if not os.path.isfile(self.System['base_dir']+self.System['movie_file_name']):
                    self.System['movie_file_name'] = self.Default['movie_file_name']
                    warnings.warn(none_template % ('movie_file_name', self.Default['movie_file_name']))
                    with open(self.System['base_dir']+self.file, 'a') as warn:
                        warn.write(none_template % ('movie_file_name', self.Default['movie_file_name']))
                    _exit()

        except Exception as e:
            self.System['movie_file_name'] = self.Default['movie_file_name']
            warnings.warn(none_template % ('movie_file_name', self.Default['movie_file_name']))
            with open(self.System['base_dir']+self.file, 'a') as warn:
                warn.write(
                    none_template % (
                        self.System['movie_file_name'], self.Default['movie_file_name']
                    )
                )
            _exit()

        with open(self.System['base_dir']+self.file, "a") as f:
            f.write('\nReading from the %s file.\n' % (self.System['base_dir']+self.System['movie_file_name']))
            
    """

    def Cenergy_file_name(self):
        
        """""""
        
        Please note that this command has since been removed due to being obsolete.

        """""""
        def _no_file():
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
    """

    def Dextend_xyz(self):

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

        with open(self.System['base_dir']+self.file, 'a') as warn:
            warn.write("\nNo specie-specific properties for homo species will be calculated in this run.\n")
        self.System['Homo'] = self.Default['Homo']

    def EHomo(self):

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
        with open(self.System['base_dir']+self.file, 'a') as warn:
            warn.write("\nNo specie-specific properties for homo species will be calculated in this run.\n")
        self.System['Hetero'] = self.Default['Hetero']

    def GHetero(self):

        with open(self.System['base_dir']+self.file, "a") as f:
            f.write("\nChecking user input for calculating homo properties in this run.\n")

        try:
            self.System['Hetero']

            if self.System['Hetero'] is None:
                self._no_hetero()
        except KeyError:
            self._no_hetero()

    def IStart(self):
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

    def JEnd(self):
        try:
            if not type(self.System['End']) is int or self.System['End'] < self.System['Start']:
                Temp = read(self.System['base_dir']+self.System['movie_file_name'], index=':')
                self.Default['End'] = len(Temp)
                self.System['End'] = len(Temp)
                del(Temp)
                warnings.warn(none_template % ('End', self.Default['End']))
                with open(self.System['base_dir']+self.file, 'a') as warn:
                    warn.write(none_template % (self.System['End'], self.Default['End']))

            elif self.System['End'] < self.System['Start']:
                Temp = read(self.System['base_dir']+self.System['movie_file_name'], index=':')
                self.Default['End'] = len(Temp)
                self.System['End'] = len(Temp)
                del(Temp)
                warnings.warn(none_template % ('End', self.Default['End']))
                with open(self.System['base_dir']+self.file, 'a') as warn:
                    warn.write(none_template % ('End', self.Default['End']))

            else:
                with open(self.System['base_dir']+self.file, 'a') as file:
                    file.write("\nFinal frame has been set to %s.\n" % self.System['End'])

        except KeyError:
            Temp = read(self.System['base_dir']+self.System['movie_file_name'], index=':')
            self.Default['End'] = len(Temp)
            self.System['End'] = len(Temp)
            del(Temp)
            warnings.warn(none_template % ('End', self.Default['End']))
            with open(self.System['base_dir']+self.file, 'a') as warn:
                warn.write(none_template % ('End', self.Default['End']))

    def KStep(self):

        try:
            if not type(self.System['Step']) is int or self.System['Step'] < 1:
                self.System['Step'] = self.Default['Step']
                warnings.warn(none_template % ('Step', self.Default['Step']))
                with open(self.System['base_dir']+self.file, 'a') as warn:
                    warn.write(none_template % ('Step', self.Default['Step']))

        except KeyError:
            self.System['Step'] = self.Default['Step']
            warnings.warn(none_template % ('Step', self.Default['Step']))
            with open(self.System['base_dir']+self.file, 'a') as warn:
                warn.write(none_template % ('Step', self.Default['Step']))

    def LSkip(self):

        try:
            if not type(self.System['Skip']) is int or self.System['Skip'] < 1:
                self.Default['Skip'] = int(self.System['End']-self.System['Start']/25.0)
                if self.Default['Skip'] < 1:
                    self.Default['Skip'] = 1
                warnings.warn(none_template % ('Skip', self.Default['Skip']))
                with open(self.System['base_dir']+self.file, 'a') as warn:
                    warn.write(none_template % ('Skip', self.Default['Skip']))
                self.System['Skip'] = self.Default['Skip']

        except KeyError:
            self.Default['Skip'] = int(self.System['End']-self.System['Start']/25.0)
            if self.Default['Skip'] < 1:
                self.Default['Skip'] = 1
            warnings.warn(none_template % ('Step', self.Default['Step']))
            with open(self.System['base_dir']+self.file, 'a') as warn:
                warn.write(none_template % ('Step', self.Default['Step']))
            self.System['Skip'] = self.Default['Skip']

    def MUniformPDF(self):

        try:
            if type(self.System['UniformPDF']) is not bool:
                warnings.warn(none_template % ('UniformPDF', self.Default['UniformPDF']))
                with open(self.System['base_dir']+self.file, 'a') as warn:
                    warn.write(none_template % ('UniformPDF', self.Default['UniformPDF']))
                self.System['UniformPDF'] = self.Default['UniformPDF']

        except KeyError:
            warnings.warn(none_template % ('UniformPDF', self.Default['UniformPDF']))
            with open(self.System['base_dir']+self.file, 'a') as warn:
                warn.write(none_template % ('UniformPDF', self.Default['UniformPDF']))
            self.System['UniformPDF'] = self.Default['UniformPDF']

    def NBand(self):

        try:
            if type(self.System['Band']) is not float:
                self.Default['Band'] = self.Default['Band']
                warnings.warn(none_template % ('Band', self.Default['Band']))
                with open(self.System['base_dir']+self.file, 'a') as warn:
                    warn.write(none_template % ('Band', self.Default['Band']))
                self.System['Band'] = self.Default['Band']

        except KeyError:
            warnings.warn(none_template % ('Band', self.Default['Band']))
            with open(self.System['base_dir']+self.file, 'a') as warn:
                warn.write(none_template % ('Band', self.Default['Band']))
            self.System['Band'] = self.Default['Band']
