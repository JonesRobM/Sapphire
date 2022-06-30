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