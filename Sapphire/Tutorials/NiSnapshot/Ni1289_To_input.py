from Post_Process import Process
import pickle

#Define the environment variables for your analysis.

System = {
    'base_dir': '',
    'movie_file_name': 'Ni1289_To.xyz',
    'Homo': ['Ni'],
    'Hetero': False,
    'Start': 0, 'End': None, 'Step': 1, 'Skip': 1
}

# Define the quantities you want calculating given the names
# in the supporting documentation.

Quantities = {
    'Full':
    {
        'euc': None, 'rdf': None, 'pos': None,  'comdist': None,
        'moi': None, 'adj': None, 'pdf': None, 'pair_distance': None,
        'agcn': {'Write_Movie': False},
        'nn': None, 'com': None, 'cna_sigs': None,
        'cna_patterns': {'Write_Movie': True},
        'gyration': None, 'stat_radius': None,
        'surf_area': None, 'surf_atoms': None
    }
}

CNA_Pattern_Settings = {
    'npz_dir': 'CNA_npz/',  # folder to save the npz files in
    'new_xyz_dir': 'CNA_XYZs/',
    'APPEND_DICTIONARY': False,
    'FROM_MEMORY': False,
    'BULK_MASTERKEY': True,
    'SAVING_XYZ': True,
    'PRINTING_PATTERNS': True
}


Data = Process.Process(
    System=System, Quantities=Quantities,
    Pattern_Input=CNA_Pattern_Settings, Cores=1
)

Meta = Data.analyse()

with open(System['base_dir']+"Metadata.csv", "wb") as file:
    pickle.dump(Data.metadata, file, protocol=pickle.HIGHEST_PROTOCOL)
