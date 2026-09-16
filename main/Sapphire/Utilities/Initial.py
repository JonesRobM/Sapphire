import getpass
import datetime
import os
import platform
import socket


class Logo():
    def __init__(self, base=''):

        self.Info_file = base+'Sapphire_Info.txt'
        self.Error_file = base+'Sapphire_Errors.log'

        self.ENDC = ''

    def Logo(self):

        return r"""
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
        from Sapphire import __version__ as _v
        self._Version_ = _v
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

    _NO_QUOTE = "\nNo random quote today.\n"

    def _quote_(self):
        """A random quote for the log header.

        This reaches out to the network, which a batch node may not be able to do.
        It is therefore opt-out (``SAPPHIRE_NO_QUOTE=1``, set by the CLI's
        ``--no-quote``), bounded by a short socket timeout, and never allowed to
        raise: a decoration must not be able to fail an analysis run.
        """
        if os.environ.get('SAPPHIRE_NO_QUOTE'):
            return self._NO_QUOTE
        try:
            import wikiquote
        except ImportError:
            return self._NO_QUOTE
        try:
            timeout = float(os.environ.get('SAPPHIRE_QUOTE_TIMEOUT', '3'))
        except ValueError:
            timeout = 3.0
        previous = socket.getdefaulttimeout()
        socket.setdefaulttimeout(timeout)
        try:
            return str(wikiquote.quotes(wikiquote.random_titles(max_titles=1)[0]))+"\n"
        except Exception:
            return self._NO_QUOTE
        finally:
            socket.setdefaulttimeout(previous)

    def _write_(self):
        with open(self.file, 'a') as Sim:
            for x in self.quants:
                Sim.write(getattr(self, x)())
