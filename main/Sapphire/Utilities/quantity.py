class Quantity():
    name = ''
    displayname = ''
    dependencies = ''
    number = ''
    
    def __init__(self, system, metadata):
        pass
    
    def not_implemented_string(self):
        return "Sorry, pal! We don't support " + self.name
    
    def calculate(self, i_frame, result_cache, metadata):
        raise NotImplementedError(self.not_implemented_string())
        
    def get_dimensions(self, n_atoms, n_elements):
        raise NotImplementedError(self.not_implemented_string())
        
    def cleanup(self):
        pass
    
    
        