import numpy as np

class Bulk_Masterkey():
    def __init__(self, Pattern_Dict = {}):
        
        self.Pattern_Dict = Pattern_Dict
        
    def Key(self):

        self.Pattern_Dict['masterkey']=np.empty((4,), dtype=object)
        self.Pattern_Dict['masterkey'][0]=[((5,5,5),12)]
        self.Pattern_Dict['masterkey'][1]=[((4,2,2),10), ((5,5,5),2)]
        self.Pattern_Dict['masterkey'][2]=[((4,2,1),6), ((4,2,2),6)]
        self.Pattern_Dict['masterkey'][3]=[((4,2,1),12)]

        return self.Pattern_Dict
    
class Pattern_Key():
    def __init__(self, Patterns = {}):
        
        self.Patterns = Patterns
        
    def Key(self):

        self.Patterns['[(12, (5, 5, 5))]']={
            'Bulk' : 'Yes', 'Description' : 'Ih Centre'
            }
        self.Patterns['[(10, (4, 2, 2)), (2, (5, 5, 5))]'] = {
            'Bulk': 'Yes', 'Description' : 'HCP Twinning Planes'
            }
        self.Patterns['[(6, (4, 2, 1)), (6, (4, 2, 2))]'] = {
            'Bulk': 'Yes', 'Description' : 'HCP Bulk'
            }
        self.Patterns['[(12, (4, 2, 1))]'] = {
            'Bulk': 'Yes', 'Description' : 'FCC Bulk'
            }

        self.Patterns['[(4, (2, 1, 1)), (4, (4, 2, 1))]']={
            'Bulk' : 'No', 'Description' : '(1 0 0) Facet'
            }
        self.Patterns['[(3, (2, 1, 1)), (2, (3, 1, 1)), (2, (4, 2, 1))]'] = {
            'Bulk': 'No', 'Description' : '(1 0 0) & (1 1 1) Bounding edges'
            }
        self.Patterns['[(1, (2, 0, 0)), (2, (2, 1, 1)), (2, (3, 1, 1)), (1, (4, 2, 1))]'] = {
            'Bulk': 'No', 'Description' : '(1 0 0) & (1 1 1) Bounding corners'
            }
        self.Patterns['[(2, (2, 0, 0)), (4, (3, 1, 1)), (1, (4, 2, 1))]'] = {
            'Bulk': 'No', 'Description' : 'MDh Re-entrance vertex (1 0 0) to re-entrance centre'
            }
        self.Patterns['[(2, (3, 0, 0)), (4, (3, 1, 1)), (2, (4, 2, 1)), (2, (4, 2, 2))]']={
            'Bulk' : 'No', 'Description' : 'Centre of MDh re-entrance'
            }
        self.Patterns['[(2, (2, 0, 0)), (1, (3, 0, 0)), (2, (3, 1, 1)), (1, (3, 2, 2)), (1, (4, 2, 1))]']={
            'Bulk' : 'No', 'Description' : 'MDh: vertices connecting re-entrance lines to 5-fold symmetry axes'
            }
        self.Patterns['[(4, (3, 1, 1)), (2, (3, 2, 2)), (2, (4, 2, 2))]']={
            'Bulk' : 'No', 'Description' : '5-Fold symmetry axis dividing (1 1 1) facets'
            }
        self.Patterns['[(5, (3, 2, 2)), (1, (5, 5, 5))]']={
            'Bulk' : 'No', 'Description' : 'Tip of 5-fold symmetry axis'
            }
        self.Patterns['[(6, (3, 1, 1)), (3, (4, 2, 1))]']={
            'Bulk' : 'No', 'Description' : '(1 1 1) Facets'
            }
        return self.Patterns   
    
class CNA_Masterkey():
    def __init__(self):
        return None
    
    def Key(self):

        self.MasterKey = ((0,0,0),
                    (1,0,0),
                    (2,0,0),(2,1,1),
                    (3,0,0),(3,1,1),(3,2,2),
                    (4,0,0),(4,1,1),(4,2,1),(4,2,2),(4,3,3),(4,4,4),
                    (5,2,1),(5,2,2),(5,3,2),(5,3,3),(5,4,4),(5,5,5),
                    (6,6,6))
        return self.MasterKey
    
class Logo():
    def __init__(self):
    
        return None
    
    def Logo(self):
        self.Sapphire = """

      _____         _____  _____  _    _ _____ _____  ______
     / ____|  /\   |  __ \|  __ \| |  | |_   _|  __ \|  ____|     ____
    | (___   /  \  | |__) | |__) | |__| | | | | |__) | |__       /\__/\
     \___ \ / /\ \ |  ___/|  ___/|  __  | | | |  _  /|  __|     /_/  \_\
     ____) / ____ \| |    | |    | |  | |_| |_| | \ \| |____    \ \__/ /
    |_____/_/    \_\_|    |_|    |_|  |_|_____|_|  \_\______|    \/__\/
                """
        return self.Sapphire
"""  
class parse_sigs():
    
    def __init__(self, system, sig_file_root, master_file_root):
        
        
        self.sf_root = sig_file_root
        self.mk_root = master_file_root
        self.system = system
        
    def parse_frame(self,i):
"""