import getpass
import datetime
import platform


class Logo():
    def __init__(self, base=''):

        self.Info_file = base+'Sapphire_Info.txt'
        self.Error_file = base+'Sapphire_Errors.log'

        self.ENDC = ''

    def Logo(self):

        return """\n
      _____         _____  _____  _    _ _____ _____  ______ 
     / ____|  /\   |  __ \|  __ \| |  | |_   _|  __ \|  ____|     ____ 
    | (___   /  \  | |__) | |__) | |__| | | | | |__) | |__       /\__/\ 
     \___ \ / /\ \ |  ___/|  ___/|  __  | | | |  _  /|  __|     /_/  \_\ 
     ____) / ____ \| |    | |    | |  | |_| |_| | \ \| |____    \ \__/ / 
    |_____/_/    \_\_|    |_|    |_|  |_|_____|_|  \_\______|    \/__\/ \n"""

    def _write_(self):
        with open(self.Info_file, 'w') as NewSim:
            NewSim.write(self.Logo())
        with open(self.Error_file, 'w') as NewSim:
            NewSim.write(self.Logo())


class Info():
    def __init__(self, base=''):

        self.file = base+'Sapphire_Info.txt'
        self._Version_ = "1.0.0"  # Increments by 0.0.1 with every minor patch, 0.1 with every large addition
        self._Units_ = "ev angstrom"  # Currently the only supported units
        self.quants = [
            '_version_', '_arch_', '_node_', '_user_',
            '_init_time_', '_units_', '_quote_'
        ]

    def _version_(self):
        return "\nRunning version  -- %s --\n" % (self._Version_)

    def _arch_(self):
        return "\nArchitecture : [ %s ]\n" % platform.machine()

    def _node_(self):
        return "\nSapphire is shining on [ %s ]\n" % (platform.node())

    def _user_(self):
        return "\nCurrent user is [ %s ]\n" % (getpass.getuser())

    def _init_time_(self):
        return "\nCalculation beginning %s\n" % (datetime.datetime.now().strftime("%a %d %b %Y %H:%M:%S"))

    def _units_(self):
        return "\nUnits : [ %s ]\n" % self._Units_

    def _quote_(self):
        try:
            import wikiquote
            Wiki = True
        except ModuleNotFoundError:
            Wiki = False
            return "\nNo random quote today.\n"
        if Wiki:
            return str(wikiquote.quotes(wikiquote.random_titles(max_titles=1)[0]))+"\n"

    def _write_(self):
        with open(self.file, 'a') as Sim:
            for x in self.quants:
                Sim.write(getattr(self, x)())
