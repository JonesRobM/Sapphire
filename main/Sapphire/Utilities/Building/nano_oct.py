from ase.cluster.octahedron import Octahedron
from ase.io import write
import numpy as np

class Nanoct:

    @classmethod    
    def _regular_octahedron(cls, symbol, length, lattice_const):
        """
        Returns a cordinates of an icosahedron cluster 

        Parameters
        ----------
        symbol: :class: `str.list`
                A string list with values indicating the elements to construct the nanocluster 
        length: :class: 'np.int'
                Number of atoms on the square edges of the complete octahedron.
        lattice_const: :class: 'np.float'
                The value of the lattice constant of the corresponding bulk structure

        Returns
        ---------
        positions : :class:'np.array'
              An array containing the spatial positions of the selected layer 
        """
        length=int(length)
        lat_con=float(lattice_const)
        cutoff = 0.0
        nanoparticle = Octahedron(symbol,length, cutoff, lat_con)
        return nanoparticle

    @classmethod
    def _truncated_octahedron(cls, symbol, length, cutoff, lattice_const):
        """
        Returns a cordinates of an icosahedron cluster 

        Parameters
        ----------
        symbol: :class: `str`
                A sequence (list) of chemical elements to construct the nanocluster 
        length: :class: 'np.int'
                Number of atoms on the square edges of the complete octahedron.
        cutoff: :class: 'np.int'
                Number of layers cut at each vertex.
        lattice_const: :class: 'np.float'
                The value of the lattice constant of the corresponding bulk structure

        Returns
        ---------
        positions : :class:'np.array'
              An array containing the spatial positions of the selected layer 
        """
        assert cutoff > 0, 'cutoff must be greater than 0'
        length=int(length)
        lat_con=float(lattice_const)
        cutoff = float(cutoff)
        nanoparticle = Octahedron(symbol,length, cutoff, lat_con)
        return nanoparticle

    @classmethod
    def _reg_truncated_octahedron(cls, symbol, cutoff, lattice_const):
        """
        Returns a cordinates of an icosahedron cluster 

        Parameters
        ----------
        symbol: :class: `str`
                A sequence (list) of chemical elements to construct the nanocluster 
        cutoff: :class: 'np.int'
                Number of layers cut at each vertex.
        lattice_const: :class: 'np.float'
                The value of the lattice constant of the corresponding bulk structure

        Returns
        ---------
        positions : :class:'np.array'
              An array containing the spatial positions of the selected layer 
        """
        assert cutoff > 0, 'cutoff must be greater than 0'
        lat_con=float(lattice_const)
        cutoff = float(cutoff)
        length = float(3*cutoff + 1)
        nanoparticle = Octahedron(symbol,length, cutoff, lat_con)
        return nanoparticle

    @classmethod
    def _cuboctahedron(cls,symbol, cutoff, lattice_const):
       """
       Returns a cordinates of an cuboctahedron cluster 

       Parameters
       ----------
       symbol: :class: `str`
                A sequence (list) of chemical elements to construct the nanocluster 
       cutoff: :class: 'np.int'
                Number of layers cut at each vertex.
       lattice_const: :class: 'np.float'
                The value of the lattice constant of the corresponding bulk structure

       Returns
       ---------
       positions : :class:'np.array'
                    An array containing the spatial positions of the selected layer 
       """
       assert cutoff > 0, 'cutoff must be greater than 0'
       lat_con=float(lattice_const)
       cutoff = float(cutoff)
       length = float(2*cutoff + 1)
       nanoparticle = Octahedron(symbol,length, cutoff, lat_con)
       return nanoparticle

    @classmethod
    def alloy_octahedron(cls, elements, length, cutoff, lattice_const, cluster_type, name):
        """
        Returns a cordinates of an cuboctahedron cluster 

        Parameters
        ----------
        elements: :class: `str`
                A sequence (list) of chemical elements to construct the nanocluster 
        length: :class: 'np.int'
                Number of atoms on the square edges of the complete octahedron.
        cutoff: :class: 'np.int'
                Number of layers cut at each vertex.
        lattice_const: :class: 'np.float'
                The value of the lattice constant of the corresponding bulk structure
        cluster_type: :class: 'str'
                A list of keywords to select four different geometries of octahedron clusters:
                1.) regular
                2.) truncated  
                3.) reg_trunc
                4.) cuboct

        name: :class:: 'str'
               Name of the file created after the formation of the nanocluster

        Returns
        ---------
        positions : :class:'np.array'
              An array containing the spatial positions of the selected layer 
        """
        length=int(length)
        lat_con=float(lattice_const)
        cutoff = float(cutoff)
        symbols=[]
        for i in range(len(elements)):
            symbols.append(elements[i])     

        if cluster_type == 'regular' :
           nano = cls._regular_octahedron(symbols[0],length,lat_con)

        elif cluster_type == 'truncated':
           nano = cls._truncated_octahedron(symbols[0],length, cutoff, lat_con)    

        elif cluster_type == 'reg_trunc':
           nano = cls._reg_truncated_octahedron(symbols[0], cutoff, lat_con)

        elif cluster_type == 'cuboct':
           nano = cls._cuboctahedron(symbols[0], cutoff, lat_con)
       
        else:
           raise AssertionError("Geometry not implemented")
       
        #Creating the alloy of the regular_octahedron
        nanoparticle = nano.get_positions()
        Ap = np.around(np.mean(nanoparticle, axis=0), decimals=3)
        dist = np.linalg.norm(nanoparticle - Ap, ord=2, axis=1)
        new_dist=[]
        for i in range(len(dist)):
           new_dist.append(round(dist[i],3))

        values, counts = np.unique(new_dist, return_counts=True)
        freq = np.dstack((values,counts)).reshape(-1,2)
        sorted_nano = nanoparticle[np.argsort(new_dist)]

        # spliting the nanoparticle in the number of declared elements
        if len(symbols) > len(values):
            raise AssertionError("Number of layers do not match Number of symbols")
        else:
          layers = np.array_split(counts,len(symbols))
          elem_col = []
          for i in range(len(layers)):
               for j in range(len(layers[i])):
                  elem_col.append(np.full(layers[i][j],symbols[i]))
          elem_col=np.array(elem_col)
    
          # Unpacking the matrix into a X,Y,Z column 
          x,y,z = sorted_nano.T
          x = np.around(x,decimals=7)
          y = np.around(y,decimals=7)
          z = np.around(z,decimals=7)

          # Unpacking the symbols to be 
          elem_tmp = elem_col.flatten()
          elements_final = np.concatenate(elem_tmp).ravel()

          cls._print_xyz(elements_final,x,y,z,name)

    @classmethod
    def _print_xyz(cls, a, b, c, d, output_file):
        """
        Returns a file a xyz file of the computed icosahedron

        Parameters
        ----------
        a : string ::class:'np.array'
              An array with the chemical symbols   

        b : float ::class:'np.array'
              An array with the x positions 

        c : float ::class:'np.array'
            An array with the y positions
    
        d : float ::class:'np.array'
            An array with the z positions

        output: 'string'
            Name of the corresponding file where atoms are printed in .xyz format         

         Returns
         ---------
         None : :class:'NoneType'
        
        """    
        with open(output_file,'w') as f:
                f.write("{}".format(len(a)))
                f.write("\n")
                f.write("\n")
                for i in range(len(a)):
                    f.write("{} {} {} {}\n".format(a[i], b[i], c[i], d[i]))
        print ("The file:", output_file, "has been created")

    

