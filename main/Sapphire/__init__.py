"""

This is the main Sapphire library hub.
From here, you may access all available Sapphire modules.

CNA performs common neighbour analyses and offers classification tools.

Graphing provides a suite of pre-made matplotlib codes for presenting your data.

IO handles all input / output calls for Sapphire.

Post_Process contains the primary set of analysis tools

Potentials contains a range of classical inter-atomic potentials

Light computes optical spectra (optional pyGDM2 dependency).

Tutorials contains example on how to run various types of Sapphire codes.

Utilities contains useful features and main module docstrings.

"""

__version__ = "1.2.0.dev0"

__all__ = ['api', 'CNA', 'Graphing', 'IO', 'Light', 'Post_Process', 'Potentials',
           'Process', 'Tutorials', 'Utilities']
