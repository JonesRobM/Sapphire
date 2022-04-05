from Graphing import Reader, Plot_Funcs
import numpy as np

#The system dictionary directs the Plotter module to find the completed simulations
#and then creates (if it does not exist) an image directory for plotted figures.

System = {
        'base_dir' : '/mnt/d/PhD/2020/November/Sapphire_Test/',
        'single_file' : True,
        'iter_dir' : False,
        'plot_dir' : 'Images/',
        'meta_name' : 'Metadata.csv',
        'save_meta' : False, 'save_errors' : False
        }

#These are the quantities which the user would like to be plotted. There will be
#a set name to call a given plotter and the arguments for that function can be 
#passed to the dictionary entry.

#One many find a list of supported plotters and their arguments in the documentation.


Quantities = { 
               'agcn_heat' : 
                      {'Name' : 'agcn_Heat.png'},
                      
               'prdf_plot' : 
                       {'Names' : ['pdf', 'rdf'], 'Frames' : range(0,15)},
                       
               'plot_stats' : 
                   {'Stats' : ['Kullback', 'JSD'], 'Quants' : ['pdf', 'rdf'], 'Temp' : False, 'Frames': range(0,15)},
                       
               'com_plot_bi' : 
                       {'Dists' : ['CoMDist', 'MidCoMDist'], 'Species' : ['Ag'], 'Frames' : range(0,300,1)},
                       
               'cna_plot' : 
                   { 'Frames' : range(0,300,20)},
                       
               'agcn_histo' : 
                       {'Frames' : range(0,300,20)},
                       
               'com_full_plot' : 
                       {'Frames' : range(0,300,20)},
                       
               'cum_com' : 
                       {'Frames' : range(0,300,50)},
                       
               'cna_traj' :    
                       {'Sigs' : [(4,2,2), (2,0,0), (5,5,5), (3,1,1), (4,2,1), (1,0,0)]},
                'h_c' : {}
        }

Pipeline = Reader.Read_Meta(System)
Metadata, Errors = Pipeline.Average()
Figures = Plot_Funcs.Plot_Funcs(Metadata, Quantities=Quantities, System=System, Errors = None)
Figures.Make_Plots()