"""

Robert:
    
    Below is a series of dictionaries which contain the output writing
    information for each quantity in the metadate. In principle, each entry in
    the respective dictionary will determine how that quantity is to be handled.
    
    
    It will be used in the following format during production runs...
    
    import OutputInfo as OI
    
    for x in Metadata.keys():
    try:
        Attrs = getattr(OI, str(x)) 
        #This will get each item's output parameters and use them to determine
        #the correct way to handle the output.
        for attr in Attrs:
            use Attrs[attr] to do something interesting.
        
        
    except Exception as e:
        #Handle the exception as necessary
        
        
    Each dictionary will have elements of the form - 
    
            'Dir' : '', #Where to place the output
            'File' : '', #What to call the file
            'Iterate' : , #Whether to expect to write multiple times
            'Bool' : , #Will quantity be boolean? T - Make empty file, F - No
            'xyz' : , If the quantiy should be available to written into an extended xyz
            'Array' : , #Prepare the write object for each element to be an array
            'Exec' : , #If the quantity is written info for the sysinfo file

"""

#########################################################################################
#########################################################################################
#########################################################################################

masterkey = {
            'Dir' : 'Exec/', 'File' : 'Masterkey', 'Iterate' : False, 'Bool' : False, 
            'xyz' : False, 'array' : True, 'Exec' : False
            }

base_dir = {
            'Dir' : 'Exec/', 'File' : 'SysInfo', 'Iterate' : False, 'Bool' : False,
            'xyz' : False, 'array' : False, 'Exec' : True
            }

movie_file_name = {
            'Dir' : 'Exec/', 'File' : 'SysInfo', 'Iterate' : False, 'Bool' : False,
            'xyz' : False, 'array' : False, 'Exec' : True
            }

energy_file_name = {
            'Dir' : 'Exec/', 'File' : 'SysInfo', 'Iterate' : False, 'Bool' : False, 
            'xyz' : False, 'array' : False, 'Exec' : True
            }

extend_xyz = {
            'Dir' : 'Exec/', 'File' : 'SysInfo', 'Iterate' : False, 'Bool' : True,
            'xyz' : False, 'array' : True, 'Exec' : True
            }

Homo = {
            'Dir' : 'Exec/', 'File' : 'SysInfo', 'Iterate' : False, 'Bool' : False,
            'xyz' : False, 'array' : True, 'Exec' : True
            }

Hetero = {
            'Dir' : 'Exec/', 'File' : 'SysInfo', 'Iterate' : False, 'Bool' : True,
            'xyz' : False, 'array' : False, 'Exec' : True
            }

Start = {
            'Dir' : 'Exec/', 'File' : 'SysInfo', 'Iterate' : False, 'Bool' : False,
            'xyz' : False, 'array' : False, 'Exec' : True
            }

End = {
            'Dir' : 'Exec/', 'File' : 'SysInfo', 'Iterate' : False, 'Bool' : False,
            'xyz' : False, 'array' : False, 'Exec' : True
            }

Step = {
            'Dir' : 'Exec/', 'File' : 'SysInfo', 'Iterate' : False, 'Bool' : False,
            'xyz' : False, 'array' : False, 'Exec' : True
            }

NFrames = {
            'Dir' : 'Exec/', 'File' : 'SysInfo', 'Iterate' : False, 'Bool' : False,
            'xyz' : False, 'array' : False, 'Exec' : True
            }

Skip = {
            'Dir' : 'Exec/', 'File' : 'SysInfo', 'Iterate' : False, 'Bool' : False,
            'xyz' : False, 'array' : False, 'Exec' : True
            }
UniformPDF = {
            'Dir' : 'Exec/', 'File' : 'SysInfo', 'Iterate' : False, 'Bool' : False,
            'xyz' : False, 'array' : False, 'Exec' : True
            }

Band = {
            'Dir' : 'Exec/', 'File' : 'SysInfo', 'Iterate' : False, 'Bool' : False,
            'xyz' : False, 'array' : False, 'Exec' : True
            }

Species = {
            'Dir' : 'Exec/', 'File' : 'SysInfo', 'Iterate' : False, 'Bool' : False,
            'xyz' : False, 'array' : False, 'Exec' : True
            }

NSpecies = {
            'Dir' : 'Exec/', 'File' : 'SysInfo', 'Iterate' : False, 'Bool' : False,
            'xyz' : False, 'array' : False, 'Exec' : True
            }

comspace = {
            'Dir' : 'Exec/', 'File' : 'comspace', 'Iterate' : False, 'Bool' : False,
            'xyz' : False, 'array' : True, 'Exec' : False
            }

#########################################################################################
#########################################################################################
#########################################################################################

rcut = {
            'Dir' : 'Time_Dependent/', 'File' : 'RCut', 'Iterate' : True, 'Bool' : False,
            'xyz' : False, 'array' : False, 'Exec' : False
            }


moi = {
            'Dir' : 'Time_Dependent/', 'File' : 'MoI', 'Iterate' : True, 'Bool' : False,
            'xyz' : False, 'array' : True, 'Exec' : False
            }

gyration = {
            'Dir' : 'Time_Dependent/', 'File' : 'Gyration', 'Iterate' : True, 'Bool' : False,
            'xyz' : False, 'array' : False, 'Exec' : False
            }

stat_radius = {
            'Dir' : 'Time_Dependent/', 'File' : 'Stat_Radius', 'Iterate' : True, 'Bool' : False,
            'xyz' : False, 'array' : False, 'Exec' : False
            }

surf_area = {
            'Dir' : 'Time_Dependent/', 'File' : 'Surf_Area', 'Iterate' : True, 'Bool' : False,
            'xyz' : False, 'array' : False, 'Exec' : False
            }

surf_atoms = {
            'Dir' : 'Time_Dependent/', 'File' : 'Surf_Atoms', 'Iterate' : True, 'Bool' : False,
            'xyz' : False, 'array' : False, 'Exec' : False
            }

com = {
            'Dir' : 'Time_Dependent/', 'File' : 'CoM', 'Iterate' : True, 'Bool' : False,
            'xyz' : False, 'array' : True, 'Exec' : False
            }

Elements = {
            'Dir' : 'Time_Dependent/', 'File' : 'Elements', 'Iterate' : True, 'Bool' : False,
            'xyz' : False, 'array' : True, 'Exec' : False
            }

NAtoms = {
            'Dir' : 'Time_Dependent/', 'File' : 'NAtoms', 'Iterate' : True, 'Bool' : False,
            'xyz' : False, 'array' : False, 'Exec' : False
            }

#########################################################################################
#########################################################################################
#########################################################################################

"""

Robert:
    
    Below are the quantities which may be written into an extended xyz file
    if this is requested by the user. 
    
    It is likely that I will set an additional flag in the dictionaries to enable
    this mode to be compatable with the "extend_xyz" logical variable.
    
"""


pattern_indices = {
            'Dir' : 'CNA/', 'File' : 'Patterns', 'Iterate' : True, 'Bool' : False,
            'xyz' : True, 'array' : True, 'Exec' : False
            }

agcn = {
            'Dir' : 'Time_Dependent/', 'File' : 'AGCN', 'Iterate' : True, 'Bool' : False,
            'xyz' : True, 'array' : True, 'Exec' : False
            }

nn = {
            'Dir' : 'Time_Dependent/', 'File' : 'NN', 'Iterate' : True, 'Bool' : False,
            'xyz' : True, 'array' : True, 'Exec' : False
            }

#########################################################################################
#########################################################################################
#########################################################################################

"""

Robert:
    
    Below are the quantities which require undisclosed array sizes and so are
    likely to be awkward to write as meaningful and useful outputs.
    
"""

rdf = {
            'Dir' : 'Time_Dependent/', 'File' : 'RDF', 'Iterate' : True, 'Bool' : False,
            'xyz' : False, 'array' : True, 'Exec' : False
            }

rdfspace = {
            'Dir' : 'Time_Dependent/', 'File' : 'RDFSpace', 'Iterate' : True, 'Bool' : False,
            'xyz' : False, 'array' : True, 'Exec' : False
            }

cna_sigs = {
            'Dir' : 'CNA/', 'File' : 'Signatures', 'Iterate' : True, 'Bool' : False,
            'xyz' : False, 'array' : True, 'Exec' : False
            }

adj = {
            'Dir' : 'Adjacency/', 'File' : 'Full', 'Iterate' : True, 'Bool' : False,
            'xyz' : False, 'array' : True, 'Exec' : False
            }

pdf = {
            'Dir' : 'Time_Dependent/', 'File' : 'PDF', 'Iterate' : True, 'Bool' : False,
            'xyz' : False, 'array' : True, 'Exec' : False
            }

pdfspace = {
            'Dir' : 'Time_Dependent/', 'File' : 'PDFSpace', 'Iterate' : True, 'Bool' : False,
            'xyz' : False, 'array' : True, 'Exec' : False
            }

comdist = {
            'Dir' : 'Time_Dependent/', 'File' : 'CoMDist', 'Iterate' : True, 'Bool' : False,
            'xyz' : False, 'array' : True, 'Exec' : False
            }

pair_distance = {
            'Dir' : 'Time_Dependent/', 'File' : 'PairDistances', 'Iterate' : True, 'Bool' : False,
            'xyz' : False, 'array' : True, 'Exec' : False
            }

pair_distancespace = {
            'Dir' : 'Time_Dependent/', 'File' : 'PairDistancesSpace', 'Iterate' : True, 'Bool' : False,
            'xyz' : False, 'array' : True, 'Exec' : False
            }