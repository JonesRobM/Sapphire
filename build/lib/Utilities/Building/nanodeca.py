import numpy as np
import itertools 
import collections
from scipy import ndimage

class Nanodeca:
    
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
    def _Decahedron(cls, p, q , r , lattice_const):
        lattice_constant = float(lattice_const)
        p = int(p)
        q = int(q)
        r = int(r)

        # Checking if p and q values are allowed
        assert (p >= 1.0 and q >= 1.0), "p and q must be greater than 0"

        # Checking if r is an allowed value
        assert (r >= 0.0), "r must be greater than or equal to 0"
    
        t = 2.0*np.pi/5.0
        b = lattice_constant/np.sqrt(2.0)
        a = b*np.sqrt(3.0)/2.0

        vertices = a * np.array([[np.cos(np.pi/2.), np.sin(np.pi/2.), 0.],
                                  [np.cos(t*1. + np.pi/2.), np.sin(t*1. + np.pi/2.), 0.],
                                  [np.cos(t*2. + np.pi/2.), np.sin(t*2. + np.pi/2.), 0.],
                                  [np.cos(t*3. + np.pi/2.), np.sin(t*3. + np.pi/2.), 0.],
                                  [np.cos(t*4. + np.pi/2.), np.sin(t*4. + np.pi/2.), 0.]])

        # Number of atoms on the five fold axis and a nice constant
        h = p + q + 2*r - 1
        g = h - q + 1 # p + 2*r

        positions = []
        # Make the five fold axis
        for j in range(h):
            pos = np.array([0.0, 0.0, j*b - (h-1)*b/2.0])
            positions.append(pos)

        # Make pentagon rings around the five fold axis
        for n in range(1, h):
            # Condition for (100)-planes
            if n < g:
                for m in range(5):
                    v1 = vertices[m-1]
                    v2 = vertices[m]
                    for i in range(n):
                        # Condition for marks re-entrence
                        if n - i < g - r and i < g - r:
                            for j in range(h-n):
                                pos = (n-i)*v1 + i*v2
                                pos += np.array([0.0, 0.0, j*b - (h-n-1)*b/2.0])
                                positions.append(pos)

        # Prepare the array for printing in xyz format
        positions = np.array(positions)
    
        x_temp = []; y_temp=[]; z_temp=[]
        for k in range(len(positions)):
           x_temp.append(positions[k][0])
           y_temp.append(positions[k][1])
           z_temp.append(positions[k][2])
       
        x = np.array(x_temp)
        y = np.array(y_temp)
        z = np.array(z_temp)

        coord = np.stack((x,y,z), axis=1)
    
        return coord

    @classmethod
    def alloy_decahedron(cls, element, p, q, r, lattice_constant, output_file):
        p = int(p)
        q = int(q)
        r = int(r)
        lat_con = float(lattice_constant)
        symbols = []
        for i in range(0,len(element)):
           symbols.append(element[i])       

        # Reorganizing the nanoparticle in terms of distances     
    
        if p == 1 and q == 1 and r == 0:
             nanoparticle = cls._Decahedron(p, q, r, lat_con)
             result = open(output_file, 'w')
             print (len(nanoparticle), file=result)
             print("", file=result)
             print (symbols[0],nanoparticle[0][0],nanoparticle[0][1],nanoparticle[0][2],file=result)
             result.close()

        if p == 1 and q == 2 and r == 0:
             nanoparticle=cls._Decahedron(p,q,r,lat_con)
             result = open(output_file, 'w')
             print (len(nanoparticle), file=result)
             print("", file=result)
             print (symbols[0],nanoparticle[0][0],nanoparticle[0][1],nanoparticle[0][2],file=result)
             print (symbols[1],nanoparticle[1][0],nanoparticle[1][1],nanoparticle[1][2],file=result)
             result.close()

        else :   
         # Computing distances from the geometrical center 
         nanoparticle = cls._Decahedron(p, q, r, lat_con)
         Ap = np.around(np.mean(nanoparticle, axis=0), decimals=3)
         dist = np.linalg.norm(nanoparticle - Ap, ord=2, axis=1)
         new_dist=[]
         for i in range(len(dist)):
             new_dist.append(round(dist[i],3))
         values, counts = np.unique(new_dist, return_counts=True)
         freq = np.dstack((values,counts)).reshape(-1,2)
         sorted_nano = nanoparticle[np.argsort(new_dist)]


        # Spliting the nanoparticle in the number of declared elements
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
       
        cls._print_xyz(elements_final, x, y, z, output_file)
    
