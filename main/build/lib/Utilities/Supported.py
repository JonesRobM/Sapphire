class Supported(object):
    
    def __init__(self):
        
        return None
    
    def Full(self):
        
        self.Supported_Full=[
                'rdf', 'cna_sigs', 'cna_patterns', 'adj', 'pdf',  'agcn', 'nn', 'com', 
                'comdist', 'moi', 'gyration', 'stat_radius', 'surf_area', 'surf_atoms',
                'concert', 'collect', 'pair_distance'
                ]
        return self.Supported_Full
        
    def Homo(self):
        self.Supported_Homo=[ 
            'hopdf', 'hordf', 'hocom', 'hoadj', 'hocomdist', 'homidcomdist', 'hopair_distance',
         'euc', 'hocna_sigs', 'hocna_patterns', 'hogyration', 'hosurf_area', 'hosurf_atoms', 'homobonds'
         ]
        return self.Supported_Homo
        
    def Hetero(self):
        self.Supported_Hetero = [ 
            'hepdf', 'herdf', 'headj', 'he_pair_distance', 'mix', 'lae', 'ele_nn',
            'heterobonds'
            ]
        return self.Supported_Hetero