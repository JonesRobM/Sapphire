from Post_Process import Process
import pickle
# Define the environment variables for your analysis.

System = {
    'base_dir': '/path/to/directory/',
    'movie_file_name': 'path/from/directory/to/movie_file_name.xyz',
    'extend_xyz': ['', '', ''],

    'Homo': ['Element1', 'Element2'],

    'Hetero': True,

    'Start': 0, 'End': None, 'Step': 1, 'Skip': 50, 'UniformPDF': False, 'Band': 0.05
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

        'Homo':
    {
        'hopdf': None, 'hordf': None,
        'hocom': None, 'hoadj': None,
        'hocomdist': None, 'homidcomdist': None,
        'euc': None, 'hocna_sigs': None,
        'hocna_patterns': None, 'hogyration': None,
        'hosurf_area': None, 'hosurf_atoms': None,
        'hopair_distance': None
    },

        'Hetero':
    {
        'hepdf': None, 'herdf': None,
        'headj': None, 'mix': None,
        'he_pair_distance': None
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

Analysis = {
    # "JSD" : ["rdf", "pdf"],
    # "Kullback" : ["rdf", "pdf"],
    # "PStat" : ["rdf", "pdf"]
}

Data = Process.Process(System=System, Quantities=Quantities,
                       Pattern_Input=CNA_Pattern_Settings)

Meta = Data.analyse(Analysis)

with open(System['base_dir']+"Metadata.csv", "wb") as file:
    pickle.dump(Data.metadata, file, protocol=pickle.HIGHEST_PROTOCOL)
