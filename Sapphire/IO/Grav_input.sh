#!/bin/bash

#SBATCH --job-name=Sapphire
#SBATCH --output=output.out
#SBATCH --mail-type=END
#SBATCH --partition=henryd
#SBATCH --cpus-per-task=24
#SBATCH --mem=10G

srun hostname
srun sleep 1

        cat > Process.py << EOF
import sys
sys.path.append('/path/to/Sapphire/')

from Post_Process import ProcessV2
import pickle
#Define the environment variables for your analysis.

System = {
        'base_dir' : '/path/to/directory/',
        'movie_file_name' : 'path/from/directory/to/movie_file_name.xyz',
        'energy_file_name' : 'path/from/directory/to/energy_file_name.out',
        'extend_xyz' : True,

        'Homo' : ['Element1', 'Element2'],  
        
        'Hetero' : True,
        
        'Start' : 0, 'End' : None, 'Step' : 1, 'Skip' : 50, 'UniformPDF' : False, 'Band' : 0.05
        }

#Define the quantities you want calculating given the names
#in the supportinf documentation.

Quantities = { 
        'Full' :
              {
                'euc' : None, 'rdf' : None, 'pos' : None,  'comdist' : None, 
                'moi' : None, 'adj' : None, 'pdf' : None, 
                'agcn' : {'Write_Movie' : False}, 
                'nn' : None, 'com' : None, 'cna_sigs' : None,
                'cna_patterns' : {'Write_Movie' : True}, 
                'gyration' : None, 'stat_radius' : None,
                'surf_area' : None, 'surf_atoms' : None,     
                'SimTime': None, 'epot': None, 'etot' : None, 
                'ekin' : None, 'edelta' : None, 'eanetot' : None, 'temp' : None
                },
        'Homo' :
                {
                    'hopdf' : None, 'hordf' : None,
                    'com' : None, 'hoadj' : None, 
                    'comdist' : None, 'midcomdist' : None, 
                    'euc' : None, 'cna_sigs' : None, 
                    'cna_patterns' : None, 'gyration' : None,
                    'surf_area' : None, 'surf_atoms' : None
                    },
        'Hetero' : 
                {
                    'hepdf' : None, 'herdf' : None, 
                    'headj' : None, 'mix' : None
                    }
                }

CNA_Pattern_Settings = {
                        'npz_dir' : 'CNA_npz/', #folder to save the npz files in
                        'new_xyz_dir' : 'CNA_XYZs/',
                        'APPEND_DICTIONARY' : True,
                        'FROM_MEMORY' : False,
                        'BULK_MASTERKEY' : True,
                        'SAVING_XYZ' : True,
                        'PRINTING_PATTERNS' : True
    }

Analysis = {
    "JSD" : ["rdf", "pdf"],
    "Kullback" : ["rdf", "pdf"],
    "PStat" : ["rdf", "pdf"]
    }

Data = ProcessV2.Process(
    System = System, Quantities = Quantities, 
    Pattern_Input = CNA_Pattern_Settings
    )

Data.Initialising()
Data.run_pdf(10)
Data.clean_pdf()
Data.run_core(10)
Data.clean_core()

#Data.New_File(Quantities = ['nn','agcn'])

Meta = Data.analyse(Analysis)

with open(System['base_dir']+"Metadata.csv", "wb") as file:
    pickle.dump(Data.metadata, file, protocol=4)

EOF

python Process.py