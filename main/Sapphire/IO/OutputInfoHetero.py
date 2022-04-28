"""

Robert:
    
    Below are quantities which are to be used only in the HETERO calculations.
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

hecut = {
            'Dir' : 'Time_Dependent/', 'File' : 'HeCut', 'Iterate' : True, 'Bool' : False,
            'xyz' : False, 'array' : False, 'Exec' : False
            }

headj = {
            'Dir' : 'Time_Dependent/', 'File' : 'HeAdj', 'Iterate' : True, 'Bool' : False,
            'xyz' : False, 'array' : True, 'Exec' : False
            }

henn = {
            'Dir' : 'Time_Dependent/', 'File' : 'HeCoordination', 'Iterate' : True, 'Bool' : False,
            'xyz' : False, 'array' : True, 'Exec' : False
            }

hepdf = {
            'Dir' : 'Time_Dependent/', 'File' : 'HePDF', 'Iterate' : True, 'Bool' : False,
            'xyz' : False, 'array' : True, 'Exec' : False
            }

hepdfspace = {
            'Dir' : 'Time_Dependent/', 'File' : 'HePDF_Space', 'Iterate' : True, 'Bool' : False,
            'xyz' : False, 'array' : True, 'Exec' : False
            }

herdf = {
            'Dir' : 'Time_Dependent/', 'File' : 'HeRDF', 'Iterate' : True, 'Bool' : False,
            'xyz' : False, 'array' : True, 'Exec' : False
            }

herdfspace = {
            'Dir' : 'Time_Dependent/', 'File' : 'HeRDF_Space', 'Iterate' : True, 'Bool' : False,
            'xyz' : False, 'array' : True, 'Exec' : False
            }

mix = {
            'Dir' : 'Time_Dependent/', 'File' : 'Mixing', 'Iterate' : True, 'Bool' : False,
            'xyz' : False, 'array' : False, 'Exec' : False
            }
ele_nn = {
            'Dir' : 'Time_Dependent/', 'File' : 'Ele_NN', 'Iterate' : True, 'Bool' : False,
            'xyz' : False, 'array' : False, 'Exec' : False
            }

lae = {
            'Dir' : 'Time_Dependent/', 'File' : 'LAE', 'Iterate' : True, 'Bool' : False,
            'xyz' : False, 'array' : False, 'Exec' : False
            }

hetero_bonds = {
            'Dir' : 'Time_Dependent/', 'File' : 'Hetero_Bonds', 'Iterate' : True, 'Bool' : False,
            'xyz' : False, 'array' : False, 'Exec' : False
            }

hepair_distance = {
            'Dir' : 'Time_Dependent/', 'File' : 'He_Pair_Distances', 'Iterate' : True, 'Bool' : False,
            'xyz' : False, 'array' : True, 'Exec' : False
            }

hepair_distancespace = {
            'Dir' : 'Time_Dependent/', 'File' : 'He_Pair_Distances_Space', 'Iterate' : True, 'Bool' : False,
            'xyz' : False, 'array' : True, 'Exec' : False
            }