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
            
            'Skip' : , #Expect to write 1/Skip steps
            'Energy' : , #Requires a given energy file?
            'Homo' : , #Homo only quantity
            'Hetero' : #Hetero only quantity

"""

#########################################################################################
#########################################################################################
#########################################################################################

masterkey = {
            'Dir' : 'Exec/', 'File' : 'Masterkey', 'Iterate' : False, 'Bool' : False, 
            'Skip' : False, 'Energy' : False, 'Homo' : False, 'Hetero' : False
            }

movie_file_name = {
            'Dir' : 'Exec/', 'File' : 'SysInfo', 'Iterate' : False, 'Bool' : False,
            'Skip' : False, 'Energy' : False, 'Homo' : False, 'Hetero' : False
            }

energy_file_name = {
            'Dir' : 'Exec/', 'File' : 'SysInfo', 'Iterate' : False, 'Bool' : False, 
            'Skip' : False, 'Energy' : True, 'Homo' : False, 'Hetero' : False
            }

extend_xyz = {
            'Dir' : 'Exec/', 'File' : 'SysInfo', 'Iterate' : False, 'Bool' : True,
            'Skip' : False, 'Energy' : False, 'Homo' : False, 'Hetero' : False
            }

Homo = {
            'Dir' : 'Exec/', 'File' : 'SysInfo', 'Iterate' : False, 'Bool' : False,
            'Skip' : False, 'Energy' : False, 'Homo' : True, 'Hetero' : False
            }

Hetero = {
            'Dir' : 'Exec/', 'File' : 'SysInfo', 'Iterate' : False, 'Bool' : True,
            'Skip' : False, 'Energy' : False, 'Homo' : False, 'Hetero' : True
            }

Start = {
            'Dir' : 'Exec/', 'File' : 'SysInfo', 'Iterate' : False, 'Bool' : False,
            'Skip' : False, 'Energy' : False, 'Homo' : False, 'Hetero' : False
            }

End = {
            'Dir' : 'Exec/', 'File' : 'SysInfo', 'Iterate' : False, 'Bool' : False,
            'Skip' : False, 'Energy' : False, 'Homo' : False, 'Hetero' : False
            }

Step = {
            'Dir' : 'Exec/', 'File' : 'SysInfo', 'Iterate' : False, 'Bool' : False,
            'Skip' : False, 'Energy' : False, 'Homo' : False, 'Hetero' : False
            }

NFrames = {
            'Dir' : 'Exec/', 'File' : 'SysInfo', 'Iterate' : False, 'Bool' : False,
            'Skip' : False, 'Energy' : False, 'Homo' : False, 'Hetero' : False
            }

Skip = {
            'Dir' : 'Exec/', 'File' : 'SysInfo', 'Iterate' : False, 'Bool' : False,
            'Skip' : False, 'Energy' : False, 'Homo' : False, 'Hetero' : False
            }
KDE = {
            'Dir' : 'Exec/', 'File' : 'SysInfo', 'Iterate' : False, 'Bool' : False,
            'Skip' : False, 'Energy' : False, 'Homo' : False, 'Hetero' : False
            }

Band = {
            'Dir' : 'Exec/', 'File' : 'SysInfo', 'Iterate' : False, 'Bool' : False,
            'Skip' : False, 'Energy' : False, 'Homo' : False, 'Hetero' : False
            }

Species = {
            'Dir' : 'Exec/', 'File' : 'SysInfo', 'Iterate' : False, 'Bool' : False,
            'Skip' : False, 'Energy' : False, 'Homo' : False, 'Hetero' : False
            }

NAtoms = {
            'Dir' : 'Exec/', 'File' : 'SysInfo', 'Iterate' : False, 'Bool' : False,
            'Skip' : False, 'Energy' : False, 'Homo' : False, 'Hetero' : False
            }

comspace = {
            'Dir' : 'Exec/', 'File' : 'comspace', 'Iterate' : False, 'Bool' : False,
            'Skip' : False, 'Energy' : False, 'Homo' : False, 'Hetero' : False
            }

#########################################################################################
#########################################################################################
#########################################################################################

R_Cut = {
            'Dir' : 'Time_Dependent/', 'File' : 'R_Cut', 'Iterate' : False, 'Bool' : False,
            'Skip' : True, 'Energy' : False, 'Homo' : False, 'Hetero' : False
            }


moi = {
            'Dir' : 'Time_Dependent/', 'File' : 'R_Cut', 'Iterate' : False, 'Bool' : False,
            'Skip' : True, 'Energy' : False, 'Homo' : False, 'Hetero' : False
            }

Gyration = {
            'Dir' : 'Time_Dependent/', 'File' : 'R_Cut', 'Iterate' : False, 'Bool' : False,
            'Skip' : True, 'Energy' : False, 'Homo' : False, 'Hetero' : False
            }

stat_radius = {
            'Dir' : 'Time_Dependent/', 'File' : 'R_Cut', 'Iterate' : False, 'Bool' : False,
            'Skip' : True, 'Energy' : False, 'Homo' : False, 'Hetero' : False
            }

surf_area = {
            'Dir' : 'Time_Dependent/', 'File' : 'R_Cut', 'Iterate' : False, 'Bool' : False,
            'Skip' : True, 'Energy' : False, 'Homo' : False, 'Hetero' : False
            }

surf_atoms = {
            'Dir' : 'Time_Dependent/', 'File' : 'R_Cut', 'Iterate' : False, 'Bool' : False,
            'Skip' : True, 'Energy' : False, 'Homo' : False, 'Hetero' : False
            }

simtime = {
            'Dir' : 'Time_Dependent/', 'File' : 'R_Cut', 'Iterate' : False, 'Bool' : False,
            'Skip' : True, 'Energy' : False, 'Homo' : False, 'Hetero' : False
            }

epot = {
            'Dir' : 'Time_Dependent/', 'File' : 'R_Cut', 'Iterate' : False, 'Bool' : False,
            'Skip' : True, 'Energy' : False, 'Homo' : False, 'Hetero' : False
            }

etot = {
            'Dir' : 'Time_Dependent/', 'File' : 'R_Cut', 'Iterate' : False, 'Bool' : False,
            'Skip' : True, 'Energy' : False, 'Homo' : False, 'Hetero' : False
            }

ekin = {
            'Dir' : 'Time_Dependent/', 'File' : 'R_Cut', 'Iterate' : False, 'Bool' : False,
            'Skip' : True, 'Energy' : False, 'Homo' : False, 'Hetero' : False
            }

edelta = {
            'Dir' : 'Time_Dependent/', 'File' : 'R_Cut', 'Iterate' : False, 'Bool' : False,
            'Skip' : True, 'Energy' : False, 'Homo' : False, 'Hetero' : False
            }

meanetot = {
            'Dir' : 'Time_Dependent/', 'File' : 'R_Cut', 'Iterate' : False, 'Bool' : False,
            'Skip' : True, 'Energy' : False, 'Homo' : False, 'Hetero' : False
            }

temp = {
            'Dir' : 'Time_Dependent/', 'File' : 'R_Cut', 'Iterate' : False, 'Bool' : False,
            'Skip' : True, 'Energy' : False, 'Homo' : False, 'Hetero' : False
            }

concert = {
            'Dir' : 'Time_Dependent/', 'File' : 'R_Cut', 'Iterate' : False, 'Bool' : False,
            'Skip' : True, 'Energy' : False, 'Homo' : False, 'Hetero' : False
            }

collect = {
            'Dir' : 'Time_Dependent/', 'File' : 'R_Cut', 'Iterate' : False, 'Bool' : False,
            'Skip' : True, 'Energy' : False, 'Homo' : False, 'Hetero' : False
            }

com = {
            'Dir' : 'Time_Dependent/', 'File' : 'R_Cut', 'Iterate' : False, 'Bool' : False,
            'Skip' : True, 'Energy' : False, 'Homo' : False, 'Hetero' : False
            }


#########################################################################################
#########################################################################################
#########################################################################################