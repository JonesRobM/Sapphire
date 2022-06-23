"""

Robert:
    
    Below are quantities which are to be used only in the HOMO calculations.
    It is likely that these shall soon be seperated for ease of implementation.
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

hocut = {
            'Dir' : 'Time_Dependent/', 'File' : 'HomoCut', 'Iterate' : True, 'Bool' : False,
            'xyz' : False, 'array' : False, 'Exec' : False
            }

hocom = {
            'Dir' : 'Time_Dependent/', 'File' : 'HomoCoM', 'Iterate' : True, 'Bool' : False,
            'xyz' : False, 'array' : True, 'Exec' : False
            }

hogyration = {
            'Dir' : 'Time_Dependent/', 'File' : 'HomoGyration', 'Iterate' : True, 'Bool' : False,
            'xyz' : False, 'array' : False, 'Exec' : False
            }

hocomdist = {
            'Dir' : 'Time_Dependent/', 'File' : 'HomoCoMDist', 'Iterate' : True, 'Bool' : False,
            'xyz' : False, 'array' : True, 'Exec' : False
            }

hosurf_area = {
            'Dir' : 'Time_Dependent/', 'File' : 'HomoSurf_Area', 'Iterate' : True, 'Bool' : False,
            'xyz' : False, 'array' : True, 'Exec' : False
            }

hosurf_atoms = {
            'Dir' : 'Time_Dependent/', 'File' : 'HomoSurf_Atoms', 'Iterate' : True, 'Bool' : False,
            'xyz' : False, 'array' : False, 'Exec' : False
            }

hordf = {
            'Dir' : 'Time_Dependent/', 'File' : 'HomoRDF', 'Iterate' : True, 'Bool' : False,
            'xyz' : False, 'array' : True, 'Exec' : False
            }

hordfspace = {
            'Dir' : 'Time_Dependent/', 'File' : 'HomoRDFSpace', 'Iterate' : True, 'Bool' : False,
            'xyz' : False, 'array' : True, 'Exec' : False
            }

hocna_patterns = {
            'Dir' : 'CNA/', 'File' : 'HomoCNA_Patterns', 'Iterate' : True, 'Bool' : False,
            'xyz' : False, 'array' : True, 'Exec' : False
            }

hocna_sigs = {
            'Dir' : 'CNA/', 'File' : 'HomoCNA_Signatures', 'Iterate' : True, 'Bool' : False,
            'xyz' : False, 'array' : True, 'Exec' : False
            }

hoadj = {
            'Dir' : 'Adjacency/', 'File' : 'HomoAdj', 'Iterate' : True, 'Bool' : False,
            'xyz' : False, 'array' : True, 'Exec' : False
            }

hopdf = {
            'Dir' : 'Time_Dependent/', 'File' : 'HomoPDF', 'Iterate' : True, 'Bool' : False,
            'xyz' : False, 'array' : True, 'Exec' : False
            }

hopdfspace = {
            'Dir' : 'Time_Dependent/', 'File' : 'HomoPDFSpace', 'Iterate' : True, 'Bool' : False,
            'xyz' : False, 'array' : True, 'Exec' : False
            }

hocomdist = {
            'Dir' : 'Time_Dependent/', 'File' : 'HomoCoMDist', 'Iterate' : True, 'Bool' : False,
            'xyz' : False, 'array' : True, 'Exec' : False
            }

homidcomdist = {
            'Dir' : 'Time_Dependent/', 'File' : 'HomoMidCoMDist', 'Iterate' : True, 'Bool' : False,
            'xyz' : False, 'array' : True, 'Exec' : False
            }

hopair_distance = {
            'Dir' : 'Time_Dependent/', 'File' : 'HomoPair_Distance', 'Iterate' : True, 'Bool' : False,
            'xyz' : False, 'array' : True, 'Exec' : False
            }

honn = {
            'Dir' : 'Time_Dependent/', 'File' : 'HomoCoordination', 'Iterate' : True, 'Bool' : False,
            'xyz' : False, 'array' : True, 'Exec' : False
            }

hopair_distancespace = {
            'Dir' : 'Time_Dependent/', 'File' : 'HomoPair_DistanceSpace', 'Iterate' : True, 'Bool' : False,
            'xyz' : False, 'array' : True, 'Exec' : False
            }

homo_bonds = {
            'Dir' : 'Time_Dependent/', 'File' : 'Homo_Bonds', 'Iterate' : True, 'Bool' : False,
            'xyz' : False, 'array' : False, 'Exec' : False
            }