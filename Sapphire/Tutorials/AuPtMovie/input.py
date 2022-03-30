from Post_Process import Process
import pickle
# Define the environment variables for your analysis.

System = {
    'base_dir': '',
    'movie_file_name': 'movie.xyz',
    'energy_file_name': 'energy.out',

    'Homo': ['Au', 'Pt'],

    'Hetero': False,

    'Start': 0, 'End': None, 'Step': 1, 'Skip': 50
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
    },

        'Energy':
    {
        'simtime': None, 'epot': None, 'etot': None,
        'ekin': None, 'edelta': None, 'meanetot': None, 'temp': None
    },

        'Homo':
    {
        'hopdf': None, 'hordf': None,
        'hocom': None, 'hoadj': None,
        'hocomdist': None, 'homidcomdist': None,
        'euc': None, 'hocna_sigs': None,
        'hocna_patterns': None, 'hogyration': None,
        'hosurf_area': None, 'hosurf_atoms': None,
        'hopair_distance': None
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
    Pattern_Input=CNA_Pattern_Settings, Cores=6
)

Meta = Data.analyse()

with open(System['base_dir']+"Metadata.csv", "wb") as file:
    pickle.dump(Data.metadata, file, protocol=pickle.HIGHEST_PROTOCOL)
