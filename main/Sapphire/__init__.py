#Your face here
"""

This is the main Sapphire library hub.
From here, you may access all available Sapphire modules.

CNA performs common neighbour analyses and offers classification tools.

Graphing provides a suite of pre-made matplotlib codes for presenting your data.

IO handles all input / output calls for Sapphire.

Post_Process contains the primary set of analysis tools

Potentials contains a range of classical inter-atomic potentials

Tutorials contains example on how to run various types of Sapphire codes.

Utilities contains useful features and main module docstrings.

"""

from os.path import dirname, basename, isfile, join
import glob
modules = glob.glob(join(dirname(__file__), "*.py"))
__all__ = [ basename(f)[:-3] for f in modules if isfile(f) and not f.endswith('__init__.py')]