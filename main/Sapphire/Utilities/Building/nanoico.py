import numpy as np


class Nanoalloy:
    
    @classmethod
    def _icosahedron_alloy(cls, layer, lattice_constant):
        """
        Returns a cordinates of an icosahedron cluster 

        Parameters
        ----------

        layer: :class:`np.int`
                An integer indicating the creation of the corresponding coordinates of a selected layer 

        lattice_const: :class: 'np.float'
                The value of the lattice constant of the corresponding bulk structure

        Returns
        ---------
        positions : :class:'np.array'
              An array containing the spatial positions of the selected layer  
        
        """   
        layer=int(layer)
        lattice_constant=float(lattice_constant)

        # golden ratio
        t = (1 + np.sqrt(5)) / 2

        vertices = np.array([[t, 0., 1.],
                         [t, 0., -1.],
                         [-t, 0., 1.],
                         [-t, 0., -1.],
                         [1., t, 0.],
                         [-1., t, 0.],
                         [1., -t, 0.],
                         [-1., -t, 0.],
                         [0., 1., t],
                         [0., -1., t],
                         [0., 1., -t],
                         [0., -1., -t]])

        positions = []
        tags = []
        positions.append(np.zeros(3))
        tags.append(1)
      
        #Construct square edges (6)
        for k in range(0, 12, 2):
                v1 = vertices[k]
                v2 = vertices[k+1]
                for i in range(layer+1):
                    pos = i*v1 + (layer-i)*v2
                    positions.append(pos)
                    tags.append(layer + 1)

        #Construct triangle planes (12)
        if layer > 1:
                map = {0: (8, 9), 1: (10, 11),
                       2: (8, 9), 3: (10, 11),
                       4: (0, 1), 5: (2, 3),
                       6: (0, 1), 7: (2, 3),
                       8: (4, 5), 9: (6, 7),
                      10: (4, 5), 11: (6, 7)}

                for k in range(0, 12):
                    v0 = layer*vertices[k]
                    v1 = (vertices[map[k][0]] - vertices[k])
                    v2 = (vertices[map[k][1]] - vertices[k])
                    for i in range(layer):
                        for j in range(layer-i):
                            if i == 0 and j == 0:
                                continue
                            pos = v0 + i*v1 + j*v2
                            positions.append(pos)
                            tags.append(layer + 1)

        #Fill missing triangle planes (8)
        if layer > 2:
                map = {0: (9, 6, 8, 4,),
                       1: (11, 6, 10, 4),
                       2: (9, 7, 8, 5,),
                       3: (11, 7, 10, 5)}

                for k in range(0, 4):
                    v0 = layer*vertices[k]
                    v1 = (vertices[map[k][0]] - vertices[k])
                    v2 = (vertices[map[k][1]] - vertices[k])
                    v3 = (vertices[map[k][2]] - vertices[k])
                    v4 = (vertices[map[k][3]] - vertices[k])
                    for i in range(1, layer):
                        for j in range(1, layer-i):
                            pos = v0 + i*v1 + j*v2
                            positions.append(pos)
                            tags.append(layer + 1)
                            pos = v0 + i*v3 + j*v4
                            positions.append(pos)
                            tags.append(layer + 1)

        # Scale the positions
        scaling_factor = lattice_constant / np.sqrt(2*(1 + t**2))
        positions = np.array(positions) * scaling_factor

        return positions

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

                
    @classmethod
    def create_alloy_ico(cls, element, noshells, lattice_const, output_file):
        """
        Returns a file a .xyz file for an alloyed nanoparticle with icosahedron geometry

        Parameters
        ----------
        noshells: string :class:`np.array`
                  An array of strings assigning the chemical element to each computed layer 
  
  
        noshells: :class:`np.int`
                  An integer indicating the number of shells to be created  


        lattice_const: :class: 'np.float'
                 The value of an estimated lattice constant of the corresponding structure    


        output: 'string'
                 Name of the corresponding file where atoms are printed in .xyz format  


        Returns
        ---------
        None : :class:'NoneType'
        
        """    
        assert int(len(element)+1) == int(noshells), "Number of outer shells (excluding the center) must be equal to number of defined elements"

        noshells=int(noshells)
        lattice_const=float(lattice_const)

        alloy_temp = []
        for i in range(1,noshells):
            alloy_temp.append(cls._icosahedron_alloy(i,lattice_const))
        alloy_temp=np.array(alloy_temp)

        alloy = []
        for i in range(0,len(alloy_temp)):
            if i == 0:
                alloy.append(alloy_temp[i][0:])
            elif i >= 1:
                alloy.append(alloy_temp[i][1:])

        symbols = []
        for i in range(0,len(element)):
           symbols.append(np.full(len(alloy[i]), element[i]))

        coord_tot = []
        for i in range(len(symbols)):
            coord_tot.append(np.c_[symbols[i],alloy[i]])

        element = []; x_temp = []; y_temp=[]; z_temp=[]
        for i in range(len(symbols)):    
            for j in range(len(symbols[i])):
                element.append(coord_tot[i][j][0])
                x_temp.append(coord_tot[i][j][1])
                y_temp.append(coord_tot[i][j][2])
                z_temp.append(coord_tot[i][j][3])

        cls._print_xyz(element, x_temp, y_temp, z_temp, output_file)

