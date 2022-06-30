import os
import sys
import warnings
from inspect import getmembers, isfunction
import inspect
import numpy as np
from ase.io import read
import scipy.sparse as sp

from Sapphire.Utilities import Initial

no_dir_template = "\nThere does not exist a suitable directory in which to place these" \
    "quantities.\n\nInstead, we shall generate one at '%s'.\n"

no_file_template = "\nThere does not exist a file in which to write the quantity %s.\n" \
    "\nInstead, we shall create the file '%s' at location '%s'."

AttrErr = "Unable to find a write object for {0}:\n"\
    "\nException traceback:\n{1}.\n"


class Writer():

    """

    Robert:

        This class object has been written with the purpose of handling the
        creation and distribution of Sapphire Output.
        In version 0.10.1, the pickle function is inadequate to facilitate
        the entirity of the metadata.

        In principle, all of the handling of output should be handled out of
        sight of the user.

    """

    def __init__(self, System, Metadata):

        self.output_info_file = System['base_dir']+'Output_Info.txt'
        self.output_error_file = System['base_dir']+'Output_Errors.txt'
        self.Quants = {
            'Dir': 'Time_Dependent/', 'File': 'R_Cut', 'Iterate': False, 'Bool': False,
            'Skip': True, 'Energy': False, 'Homo': False, 'Hetero': False, 'xyz': False
        }

        self.Metadata = Metadata  # This is the data provided to the user by Sapphire after post processing
        self.System = System  # Significant system information regarding I/O streams
        self.Logo = Initial.Logo().Logo()

        with open(self.output_info_file, 'w') as outfile:
            outfile.write(self.Logo)
            outfile.write('\n')

        with open(self.output_error_file, 'w') as outfile:
            outfile.write(self.Logo)
            outfile.write('\n')

        """
        
        This provides a dictionary with the function names as keys and the 
        function itself.
        This allows us to have 1-1-1 mapping between the output p
        
        """

        self.functions_list = [o for o in getmembers(Writer) if isfunction(o[1])]
        self.Functions = {}

        for x in self.functions_list:
            if x in self.Quants.keys():
                self.Functions[x[0]] = inspect.getfullargspec(x[1])[0][1:]

    def ensure_dir(self, base_dir='', file_path=''):
        """

        Robert:

            A simple script to verify the existence of a directory
            given the path to it. If it does not exist, will create it.

        """

        directory = base_dir + file_path
        if not os.path.exists(directory):

            os.makedirs(directory)
            with open(self.output_info_file, 'w') as outfile:
                outfile.write(no_dir_template % (base_dir+file_path))

    def MakeFile(self, Attributes):
        self.out = self.System['base_dir'] + Attributes['Dir'] + Attributes['File']

        if not os.path.isfile(self.out):
            with open(self.System['base_dir'] + Attributes['Dir'] + Attributes['File'], 'w') as out:
                out.close()
        else:
            pass

    def Masterkey(self, Quantity):
        try:
            with open(self.out, 'w') as f:
                for item in self.Metadata[self.x]:
                    f.write(str(item)+'\n')
        except Exception as e:
            with open(self.output_error_file, 'a') as outfile:
                outfile.write(AttrErr % (self.x, e))

    def Adj(self, Quantity):
        self.out = self.System['base_dir'] + Quantity['Dir'] + Quantity['File']
        self.ensure_dir(base_dir=self.System['base_dir'], file_path=Quantity['Dir'])
        for i, t in enumerate(self.Metadata[self.x]):
            try:
                self.filename = self.System['base_dir'] + Quantity['Dir'] + 'File%s' % i
                self.Mat = sp.csr_matrix.todense(t)
                with open(self.filename, 'w') as f:
                    for line in self.Mat:
                        np.savetxt(f, line, fmt='%d')
            except Exception as e:
                with open(self.output_error_file, 'a') as outfile:
                    outfile.write(AttrErr % (self.x, e))

    def Ele(self, Quantity):
        self.out = self.System['base_dir'] + Quantity['Dir'] + Quantity['File']
        self.ensure_dir(base_dir=self.System['base_dir'], file_path=Quantity['Dir'])
        with open(self.out, 'w') as file:
            for i, t in enumerate(self.Metadata[self.x]):
                try:
                    self.filename = self.System['base_dir'] + Quantity['Dir'] + 'File%s' % i
                    file.write('\t|\t'.join(str(item) for item in t[0])+'\n')
                except Exception as e:
                    with open(self.output_error_file, 'a') as outfile:
                        outfile.write(AttrErr % (self.x, e))

    def HeAdj(self, Quantity):
        self.Homo = self.System['Homo']
        for Ele in self.Homo:
            if len(self.Metadata[self.x]) > 1:
                Temp = np.column_stack((
                    self.Metadata[self.x][0][self.Homo.index(Ele)],
                    self.Metadata[self.x][1][self.Homo.index(Ele)]
                ))
                for t in range(2, len(self.Metadata[self.x])):
                    Temp = np.column_stack((
                        Temp, np.array(self.Metadata[self.x][t][self.Homo.index(Ele)], int)
                    ))
                np.savetxt(
                    self.System['base_dir'] + Quantity['Dir'] + Quantity['File']+Ele,
                    Temp.transpose(), fmt='%d')
            else:
                np.savetxt(
                    self.System['base_dir'] + Quantity['Dir'] + Quantity['File']+Ele,
                    np.array(self.Metadata[self.x][0][self.Homo.index(Ele)]).transpose(),
                    fmt='%d')

    def Write_Homo(self, Quantity):
        # self.MakeFile(Quantity) #See if the file already exists
        for Ele in self.System['Homo']:
            File = str(self.x)[:-2]+Ele
            self.out = self.System['base_dir'] + Quantity['Dir'] + Quantity['File']+Ele
            self.ensure_dir(base_dir=self.System['base_dir'], file_path=Quantity['Dir'])
            try:

                if not Quantity['Iterate'] and not Quantity['Bool'] and not Quantity['array']:
                    try:
                        np.savetxt(self.out, self.Metadata[File], fmt='%s')
                    except Exception as e:
                        with open(self.output_error_file, 'a') as error:
                            error.write(AttrErr.format(str(File), str(e)))
                        try:
                            with open(self.out, 'a') as CurrentOut:
                                CurrentOut.write(str(File)+str(self.Metadata[File]))
                                CurrentOut.write('\n')
                        except Exception as e:
                            with open(self.output_error_file, 'a') as outfile:
                                outfile.write(AttrErr % (File, e))

                elif Quantity['Iterate'] and Quantity['array']:
                    try:
                        if len(self.Metadata[File]) > 1:
                            Temp = np.column_stack((self.Metadata[File][0], self.Metadata[File][1]))
                            for t in range(2, len(self.Metadata[File])):
                                Temp = np.column_stack((Temp, self.Metadata[File][t]))
                            np.savetxt(self.out, Temp.transpose(), fmt='%f')
                        else:
                            np.savetxt(
                                self.out,
                                np.array(self.Metadata[File][0]).transpose(),
                                fmt='%f')
                    except Exception as e:
                        with open(self.output_error_file, 'a') as outfile:
                            outfile.write(AttrErr % (File, e))

                elif Quantity['Iterate'] and not Quantity['array']:
                    try:
                        np.savetxt(self.out, np.array(self.Metadata[File], dtype=float).transpose(), fmt='%f')
                    except Exception as e:
                        with open(self.output_error_file, 'a') as outfile:
                            outfile.write(AttrErr % (File, e))

            except Exception as e:
                with open(self.output_error_file, 'a') as error:
                    error.write(AttrErr.format(str(File), str(e)))

    def Write(self, Quantity):
        self.out = self.System['base_dir'] + Quantity['Dir'] + Quantity['File']
        self.ensure_dir(base_dir=self.System['base_dir'], file_path=Quantity['Dir'])  # See if the directory already exists
        # self.MakeFile(Quantity) #See if the file already exists

        if Quantity['Exec']:
            try:
                with open(self.out, 'a') as CurrentOut:
                    CurrentOut.write(str(self.x)+'\t|\t'+str(self.Metadata[self.x]))
                    CurrentOut.write('\n')
            except Exception as e:
                with open(self.output_error_file, 'a') as outfile:
                    outfile.write(AttrErr % (self.x, e))
        else:

            try:
                if Quantity['Bool']:
                    try:
                        with open(self.out, 'a') as CurrentOut:
                            CurrentOut.write(str(self.x) + '\t|\t' + str(self.Metadata[self.x]))
                            CurrentOut.write('\n')
                    except Exception as e:
                        with open(self.output_error_file, 'a') as outfile:
                            outfile.write(AttrErr % (self.x, e))

                elif not Quantity['Iterate'] and not Quantity['Bool'] and not Quantity['array']:
                    try:
                        np.savetxt(self.out, self.Metadata[self.x], fmt='%s')
                    except Exception as e:
                        with open(self.output_error_file, 'a') as error:
                            error.write(AttrErr.format(str(self.x), str(e)))
                        try:
                            with open(self.out, 'a') as CurrentOut:
                                CurrentOut.write(str(self.x)+str(self.Metadata[self.x]))
                                CurrentOut.write('\n')
                        except Exception as e:
                            with open(self.output_error_file, 'a') as outfile:
                                outfile.write(AttrErr % (self.x, e))

                elif Quantity['Iterate'] and Quantity['array']:
                    try:
                        if len(self.Metadata[self.x]) > 1:
                            Temp = np.column_stack((self.Metadata[self.x][0], self.Metadata[self.x][1]))
                            for t in range(2, len(self.Metadata[self.x])):
                                Temp = np.column_stack((Temp, self.Metadata[self.x][t]))
                            np.savetxt(self.out, Temp.transpose(), fmt='%f')
                        else:
                            np.savetxt(
                                self.out,
                                np.array(self.Metadata[self.x][0]).transpose(),
                                fmt='%f')
                    except Exception as e:
                        with open(self.output_error_file, 'a') as outfile:
                            outfile.write(AttrErr % (self.x, e))

                elif Quantity['Iterate'] and not Quantity['array']:
                    try:
                        np.savetxt(self.out, np.array(self.Metadata[self.x], dtype=float).transpose(), fmt='%f')
                    except Exception as e:
                        with open(self.output_error_file, 'a') as outfile:
                            outfile.write(AttrErr % (self.x, e))

            except Exception as e:
                with open(self.output_error_file, 'a') as error:
                    error.write(AttrErr.format(str(self.x), str(e)))

    def Run(self, Output_Type):
        """
        Robert.

            This will need to be handled internally delicately so as to not confuse
            the user.

            I would like to be able to determine whether or not to call a given
            output file type based on it being part of the Full, Homo, or Hetero
            sub-systems.

            In principle, the User is at liberty (not now, but soon) to pre-select their
            own output parameters. Though deviating from the defaults could be dangerous.


            At present, one of three string-types can be assigned to the 'Output_Type'
            free variable:

                Full - Loads in the OutputInfoFull.py file for its attributes to be read.
                Homo - Loads in the OutputInfoHomo.py file for its attributes to be read.
                Hetero - Loads in the OutputInfoHetero.py file for its attributes to be read.

        """
        if Output_Type == 'Full':
            from Utilities import OutputInfoFull as Out  # Case 1

        elif Output_Type == 'Homo':
            from Utilities import OutputInfoHomo as Out  # Case 2

        elif Output_Type == 'Hetero':
            from Utilities import OutputInfoHetero as Out  # Case 3

        self.Write_List = []

        for self.x in self.Metadata.keys():  # Things such as: 'pdf', 'R_Cut', ...

            try:
                if Output_Type == 'Homo' and self.x.startswith('ho'):
                    Attributes = getattr(Out, str(self.x[:-2]))  # Pulls dictionaries with names corresponding to x as above
                    with open(self.output_info_file, 'a') as outfile:
                        outfile.write('Working now with %s and placing it in %s with file name %s.\n' % (self.x, Attributes['Dir'], Attributes['File']))
                    try:
                        self.Write_Homo(Attributes)
                    except Exception as e:
                        with open(self.output_error_file, 'a') as error:
                            error.write(AttrErr.format(str(self.x), str(e)))

                else:
                    Attributes = getattr(Out, str(self.x))  # Pulls dictionaries with names corresponding to x as above
                    if self.x == 'adj':
                        try:
                            self.Adj(Attributes)
                        except Exception as e:
                            with open(self.output_error_file, 'a') as error:
                                error.write(AttrErr.format(str(self.x), str(e)))

                    elif self.x == 'Elements':
                        try:
                            self.Ele(Attributes)
                        except Exception as e:
                            with open(self.output_error_file, 'a') as error:
                                error.write(AttrErr.format(str(self.x), str(e)))

                    elif self.x == 'headj':
                        try:
                            self.HeAdj(Attributes)
                        except Exception as e:
                            with open(self.output_error_file, 'a') as error:
                                error.write(AttrErr.format(str(self.x), str(e)))

                    elif self.x == 'master':
                        try:
                            self.Masterkey(Attributes)
                        except Exception as e:
                            with open(self.output_error_file, 'a') as error:
                                error.write(AttrErr.format(str(self.x), str(e)))

                    else:
                        self.Write(Attributes)
                    with open(self.output_info_file, 'a') as outfile:
                        outfile.write('Working now with %s and placing it in %s with file name %s.\n' % (self.x, Attributes['Dir'], Attributes['File']))
            except Exception as e:
                with open(self.output_error_file, 'a') as error:
                    error.write(AttrErr.format(str(self.x), str(e)))
        try:
            from CNA.Utilities import Pattern_Key as PK
            self.pattern_key = PK().Key()
            with open(self.System['base_dir'] + 'RecognisedPatterns.txt', 'w') as outfile:
                for i, thing in enumerate(self.pattern_key.keys()):
                    outfile.write(str(i) + ')\t' + str(thing)+':\t')
                    for item in self.pattern_key[thing]:
                        outfile.write(str(item) + ':\t' + str(self.pattern_key[thing][item])+'\t|\t')
                    outfile.write('\n\n')
        except Exception as e:
            with open(self.output_error_file, 'a') as error:
                error.write(AttrErr.format('CNA_Patterns', e))
