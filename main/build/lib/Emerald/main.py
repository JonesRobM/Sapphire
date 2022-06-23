import math
import time
import matplotlib.pyplot as plt
import numpy as np
from math import sqrt

#Input: file xyz
#Output: matrix of i atoms in row with columns of symbol, x, y, z,
#number of atoms and title

class Emerald(object):
    
    """
    Robert:
        In principle, this will serve as the class descriptor for the key actions
        that Emerald purports to support.
        
        Ideally this shouold be either a stand-alone module or to be called from
        the primary Sapphire routine.

    Parameters
    ----------
    filename : TYPE
        DESCRIPTION.

    Returns
    -------
    coordinates1 : TYPE
        DESCRIPTION.

    """    
    
    def __init__(self, *kwargs):
        
        return None

def read_file(filename):
    """
    Robert:
        This simply returns the atoms-like object which may be supported by ase.
        I do not believe that it is necessary.
    """
    
    coordinates1 = []
    xyz = open(filename)
    n_atoms1 = int(xyz.readline())
    title1 = xyz.readline()
    for line in xyz:
        atom, x, y, z = line.split()
        coordinates1.append([atom, float(x), float(y), float(z)])
    xyz.close()
    print("filename: %s" % filename)
    print("title: %s" % title1)
    print("number of atoms: %d" % n_atoms1)
    print("number of coordinates: %d" % len(coordinates1))
    return coordinates1, n_atoms1, title1

#Input: matrix of coordinates Symbol, x, y, z
#and the number of atoms (len(coordinates)).
#Output: coordinates of center of mass
def get_coordinatesCM(coordinates, n):
    """
    This function assumes that all species are of equivilent mass!
    Do not use - core Sapphire does this better.
    """
    
    xcm = 0.
    ycm = 0.
    zcm = 0.
    for i in range(len(coordinates)):
        xcm = coordinates[i][1]+xcm
    xcm1 = xcm/n
    for i in range(len(coordinates)):
        ycm = coordinates[i][2]+ycm
    ycm1 = ycm/n
    for i in range(len(coordinates)):
        zcm = coordinates[i][3]+zcm
    zcm1 = zcm/n
    return xcm1, ycm1, zcm1

#Input: coordinates x, y, z
#Output: matrix of Euclidean distance
def Euc_Dist(coordinates):
    """
    Robert:
        This function already exists more efficiently in DistFuncs.
    """
    Distance = []
    for i in range(len(coordinates)-1):
        for j in range(i+1, len(coordinates)):
            Euc = get_distance(coordinates, i, j)
            Distance.append(Euc)
    return Distance


#Input: matrix of coordinates: Symbol, x, y, z
#the coordinates of mass center get_coordinatesCM
#Output: matrix of coordinates rescaled: Symbol, x, y, z
def riscale_coordination(coordinates, x, y, z):
    """
    Robert:
        I do not understand why you would want this quantity.
    """
    coordinatesCM = []
    for i in range(len(coordinates)):
        coordinatesCM.append(
            [str(coordinates[i][0]), coordinates[i][1]-x,
             coordinates[i][2]-y, coordinates[i][3]-z])
    return coordinatesCM

#Input: matrix of coordinates: Symbol, x, y, z
#Output: distance i j
def get_distance(coordinates, i, j):
    return math.sqrt(
        pow(coordinates[i][1] - coordinates[j][1], 2)
        + pow(coordinates[i][2] - coordinates[j][2], 2)
        + pow(coordinates[i][3] - coordinates[j][3], 2))

#Input: matrix of coordinates (use the riscale_coordination: Symbol, x, y, z)
#Output: radius
def get_radius(coordinates):
    radius = []
    for i in range(len(coordinates)):
        radius.append(
            math.sqrt(pow(coordinates[i][1], 2)
                      + pow(coordinates[i][2], 2)
                      + pow(coordinates[i][3], 2)))
    return radius

#Input: coordinates and cutoff
#Output: coordination number of all atoms in coordinates
#and a matrix quali1 which contain the nearest neighbour
#index
def get_coordination_number(coordinates, cutoff1):
    count1 = 0
    count = []
    quale = 0
    quali = []
    quali1 = []
    for i in range(len(coordinates)):
        for j in range(len(coordinates)):
            if get_distance(coordinates, i, j) <= cutoff1 and \
                    i != j:
                count1 = count1 + 1
                quale = j
                quali.append(quale)
        quali1.append(quali)
        quali = []
        count.append(count1)
        count1 = 0
    return count, quali1

#Input: coordinates and cutoff
#Output: general coordination number of all atoms
def get_generalized_CN(coordinates, cutoff1):
    count2 = 0
    generalcount = []
    count, quali = get_coordination_number(coordinates, cutoff1)
    #For every i we use the line quali[i]
    for i in range(len(coordinates)):
        for j in range(len(quali[i])):
            count2 = count2 + count[quali[i][j]]
        #FCC 12 max CN
        generalcount.append(count2 / 12)
        count2 = 0
    return generalcount

#Input: coordinates and cutoff
#Output: solid angol of all atoms
def get_solid_angol(coordinates, cutoff1):
    temp = []
    solidangle = []
    rho = 0.
    temp1 = 0.
    count, quali = get_coordination_number(coordinates, cutoff1)
    for i in range(len(coordinates)):
        for j in range(len(quali[i])):
            temp1 = get_distance(coordinates, i, quali[i][j])
            temp.append(temp1)
        m = len(temp)
        for q in range(len(temp)):
            rho += temp[q]
        rho1 = rho/(m-2)
        rho = 0.
        solidangle.append(
            (pow(rho1, 2)*math.pi)/pow(math.sqrt(
                pow(coordinates[i][1], 2)
                + pow(coordinates[i][2], 2)
                + pow(coordinates[i][3], 2)), 2))
    return solidangle

#Input: coordinates and cutoff
#Output: the list of the surface atoms and the
#bulk list as: Symbol, x, y, z, count, general count,
#radius and solid angle.
#The output is also the individual
#columns of radius, count and general count
def get_which_surface(coordinates, cutoff1):
    bulk = []
    surface = []
    radiuss = []
    countt = []
    gcn = []
    count, quali = get_coordination_number(coordinates, cutoff1)
    radius = get_radius(coordinates)
    solidangle = get_solid_angol(coordinates, cutoff1)
    generalcount = get_generalized_CN(coordinates, cutoff1)
    for i in range(len(coordinates)):
        if count[i] <= 10:
            surface.append(
                [str(coordinates[i][0]), coordinates[i][1],
                 coordinates[i][2], coordinates[i][3],
                 count[i], generalcount[i],
                 radius[i], solidangle[i]])
            radiuss.append(radius[i])
            countt.append(count[i])
            gcn.append(generalcount[i])
        else:
            bulk.append(
                [str(coordinates[i][0]), coordinates[i][1],
                 coordinates[i][2], coordinates[i][3],
                 count[i], generalcount[i],
                 radius[i], solidangle[i]])

    return surface, radiuss, countt, bulk, gcn

def get_which_bulk(coordinates, cutoff1):
    bulk = []
    count, quali = get_coordination_number(coordinates, cutoff1)
    radius = get_radius(coordinates)
    solidangle = get_solid_angol(coordinates, cutoff1)
    generalcount = get_generalized_CN(coordinates, cutoff1)
    for i in range(len(coordinates)):
        if count[i] > 10:
            bulk.append(
                [str(coordinates[i][0]), coordinates[i][1],
                 coordinates[i][2], coordinates[i][3],
                 count[i], generalcount[i],
                 radius[i], solidangle[i]])
    return bulk

#Input: coordinates, cutoff and atomic radius
#Output: the surface in Å^2
def get_surface(coordinates, cutoff1, rAT):
    surface, radius, count, \
    bulk, generalcount = get_which_surface(coordinates, cutoff1)
    surface1 = 0.
    for i in range(len(surface)):
        #For FCC max 12
        surface1 += (1-(generalcount[i]/12))
    return surface1*4*math.pi*pow(rAT, 2)

def get_diff_surface(coordinates, cutoff1, rAT):
    surfacemin = 0.
    surfacemax = 0.
    radiuss = get_radius(coordinates)
    surface, radius, generalcount, \
    bulk, generalcount = get_which_surface(coordinates, cutoff1)
    rmin = min(radius)
    rmax = max(radiuss)
    surfacemin = 4 * math.pi*pow(rmin+rAT, 2)
    surfacemax = 4 * math.pi * pow(rmax+rAT, 2)
    return surfacemax, surfacemin

def get_medium_radius_surface(coordinates, cutoff1):
    rad1 = 0.
    r = []
    deltarad1 = 0.
    surface, radius, count,\
    bulk, generalcount = get_which_surface(coordinates, cutoff1)
    rad = 0.
    for i in range(len(radius)):
        rad += radius[i]
    rad1 = rad/len(radius)
    for i in range(len(radius)):
        r.append(abs(radius[i] - rad1))
    deltarad1 = max(r)
    return rad1, deltarad1

def get_medium_radius_peeledsurface(coordinates, cutoff1):
    rad1 = 0.
    r = []
    deltarad1 = 0.
    surface, radius, count,\
    bulk, generalcount = get_which_surface(coordinates, cutoff1)
    surfacepeel, radius3, \
    generalcount3, bulk4 = get_which_surface_peeling(bulk, cutoff)
    rad = 0.
    for i in range(len(surfacepeel)):
        rad += radius3[i]
    rad1 = rad/len(radius3)
    for i in range(len(radius3)):
        r.append(abs(radius3[i] - rad1))
    deltarad1 = max(r)
    return rad1, deltarad1

#Input: coordinates and cutoff
#Output: the list of the peeled surface atoms and the
#bulk list as: Symbol, x, y, z, count, general count,
#radius and solid angle.
#The output is also the individual
#columns of radius, count and general count
def get_which_surface_peeling(coordinates, cutoff1):
    surface = []
    radiuss = []
    countt = []
    bulk = []
    count, quali = get_coordination_number(coordinates, cutoff1)
    radius = get_radius(coordinates)
    solidangle = get_solid_angol(coordinates, cutoff1)
    generalcount = get_generalized_CN(coordinates, cutoff1)
    for i in range(len(coordinates)):
        if count[i] <= 10:
            surface.append(
                [str(coordinates[i][0]), coordinates[i][1],
                 coordinates[i][2], coordinates[i][3],
                 count[i], generalcount[i],
                 radius[i], solidangle[i]])
            radiuss.append(radius[i])
            countt.append(count[i])
        else:
            bulk.append(
                [str(coordinates[i][0]), coordinates[i][1],
                 coordinates[i][2], coordinates[i][3],
                 count[i], generalcount[i],
                 radius[i], solidangle[i]])

    return surface, radiuss, countt, bulk

#Input: coordinates, atomic radius of the element
# and the cutoff
#Output: the total volume of the cluster
def get_volume(coordinates, rAT, cutoff1):
    r2 = []
    r20 = 0.
    volumemin = 0.
    volumemancante = 0.
    surface, radius, count,\
    bulk, generalcount = get_which_surface(coordinates, cutoff1)
    rmin = min(radius) - rAT
    volumemin = (4/3)*math.pi*pow(rmin, 3)
    volume = 0.
    volumemancante = (2/3)*math.pi*pow(rAT, 3)*len(radius)
    for i in range(len(radius)):
        r20 = rAT-sqrt(pow(pow(pow(radius[i], 2) +
                                 pow(rAT, 2), 0.5) - rmin, 2) -
                         pow(radius[i] - rmin, 2))
        r2.append(r20)
    for j in range(len(r2)):
        volume += (1/3)*math.pi*(pow(rAT, 2) +
                                 rAT*r2[j] + pow(r2[j], 2))*\
                  (radius[j] - rmin)
    return volume + volumemin + volumemancante

#Input: the coordinates and the Wiener radius
#Output: the volume 4/3*pi*r^3
def get_volume_WS(coordinates, rWS1):
    N = 0.
    volume = 0.
    N = len(coordinates)
    volume = (4/3)*math.pi*N*pow(rWS1, 3)
    return volume

#Input: coordinates and the radius
#Output: the count of atoms i found in the
#shell r-0.2, r+0.2
def get_RDF(coordinates, r):
    count = 0.
    k = []
    radius = get_radius(coordinates)
    for i in range(len(coordinates)):
        if i not in k:
            if radius[i]-0.2 < r and radius[i]+0.2 > r:
                count += 1
                k.append(i)
    return count

#From Robert Michael Jones King's College London (UK)
#Input: coordinates and the bins, is defined out
#of the function so it will be the same for each
#cases of surface, peeled surface and core
#Output: the x and y for the plot
def get_PDDF(coordinates, bins):
    distances = Euc_Dist(coordinates)
    a, b = np.histogram(distances, bins)
    bin_width = b[1] - b[0]
    bin_cents = [b[i] + bin_width for i in range(len(b)-1)]
    return bin_cents, a


#Input: coordinates and bins
#Output: array of radius and the count for
#any radius for the histogram
def get_histo(coordinates, bins):
    lunghezza = []
    radius2primo = []
    count5 = 0.
    radius = get_radius(coordinates)
    rmin = min(radius)
    rmax = max(radius)
    Nint = int((rmax - rmin)/bins)
    rmin2 = round(rmin, 1)
    for j in range(0, Nint+1):
        for i in range(len(radius)):
            if rmin + bins * j  <= \
                    radius[i] and radius[i] <= \
                    rmin + bins * j + bins:
                count5 += 1
        lunghezza.append(count5)
        count5 = 0.
        radius2primo.append(rmin2 + (bins/2)*(2*j+1))
    return radius2primo, lunghezza

def number_count(coordinates, cutoff1):
    d = []
    count = 0.
    count2 = []
    surface, radius, count1, \
    bulk, generalcount = get_which_surface(coordinates, cutoff1)
    for element in count1:
        if element not in d:
            d.append(element)
    for i in range(len(d)):
        for j in range(len(surface)):
            if d[i] == count1[j]:
                count += 1
        count2.append(count)
        count = 0.

    return d, count2

def number_count_total(coordinates, cutoff1):
    d = []
    count = 0
    count2 = []
    count1, quali = get_coordination_number(coordinates, cutoff1)
    for element in count1:
        if element not in d:
            d.append(element)
    for i in range(len(d)):
        for j in range(len(coordinates)):
            if d[i] == count1[j]:
                count += 1
        count2.append(count)
        count = 0.

    return d, count2


def number_generalcount(coordinates, cutoff1):
    d = []
    count = 0.
    count2 = []
    surface, radius, count1, bulk, generalcount = get_which_surface(coordinates, cutoff1)
    for element in generalcount:
        if element not in d:
            d.append(element)
    for i in range(len(d)):
        for j in range(len(surface)):
            if d[i] == generalcount[j]:
                count += 1
        count2.append(count)
        count = 0.

    return d, count2

def colored_GCN(coordinates, cutoff1):
    surface, radius, count1, \
    bulk, generalcount = get_which_surface(coordinates, cutoff1)
    d = []
    c = []
    total = []
    for element in generalcount:
        if element not in d:
            d.append(element)
    d.sort()
    n = 0
    for i in range(len(d)):
        c.append(str(n))
        n += 1
    for i in range(len(surface)):
        for j in range(len(d)):
            if generalcount[i] == d[j]:
                total.append([c[j], surface[i][1], surface[i][2], surface[i][3]])
    return total, d

def colored_GCN2(coordinates, cutoff1):
    surface, radius, count1, \
    bulk, generalcount = get_which_surface(coordinates, cutoff1)
    generalcount2 = get_generalized_CN(surface, cutoff1)
    d = []
    c = []
    total = []
    for element in generalcount2:
        if element not in d:
            d.append(element)
    n = 0
    for i in range(len(d)):
        c.append(str(n))
        n += 1
    for i in range(len(surface)):
        for j in range(len(d)):
            if generalcount2[i] == d[j]:
                total.append([c[j], surface[i][1], surface[i][2], surface[i][3]])
    return total

def get_plane_colored_GCN(coordinates, cutoff1):
    surface, radius, count1, \
    bulk, generalcount = get_which_surface(coordinates, cutoff1)
    uno = []
    a = []
    b = []
    c = []
    d = []
    a1 = []
    b1 = []
    c1 = []
    d1 = []
    matrixplane = []
    matrixplane1 = []
    dd = []
    cc = []
    total = []
    j = 0
    m = 0
    k = 0
    n = 0
    for element in generalcount:
        if element not in dd:
            dd.append(element)
    dd.sort()
    for i in range(len(surface)):
        if generalcount[i] == dd[0]:
            uno.append([surface[i][0], surface[i][1], surface[i][2], surface[i][3]])
        if generalcount[i] == dd[1]:
            uno.append([surface[i][0], surface[i][1], surface[i][2], surface[i][3]])

    for i in range(len(uno)-2):
        for j in range(i+1, len(uno)-1):
            for m in range(j+1, len(uno)):
                if i != j and j != m:
                    a.append(round((uno[j][2] - uno[i][2]) *
                                   (uno[m][3] - uno[i][3]) -
                                   (uno[j][3] - uno[i][3]) *
                                   (uno[m][2] - uno[i][2]), 2))
                    b.append(round((uno[j][3] - uno[i][3]) *
                             (uno[m][1] - uno[i][1]) -
                             (uno[j][1] - uno[i][1]) *
                             (uno[m][3] - uno[i][3]), 2))
                    c.append(round((uno[j][1] - uno[i][1]) *
                             (uno[m][2] - uno[i][2]) -
                             (uno[j][2] - uno[i][2]) *
                             (uno[m][1] - uno[i][1]), 2))
                    d.append(round(-a[i] * uno[i][1] - b[i] *
                             uno[i][2] - c[i] * uno[i][3], 2))
    for i in range(len(a)):
        if a[i] != 0. and b[i] != 0. and c[i] != 0. and d[i] != 0.:
            matrixplane.append([a[i], b[i], c[i], d[i]])

    return matrixplane, uno, dd


def get_GCN_atom_plane1(coordinates, cutoff1, rAT):
    matrixplane, uno, dd = get_plane_colored_GCN(coordinates, cutoff1)
    surface, radius, count1,\
    bulk, generalcount = get_which_surface(coordinates, cutoff1)
    numero = []
    quali = []
    quali1 = []
    quale = 0
    d = 0.
    n = 0
    for i in range(len(matrixplane)):
        for j in range(len(surface)):
            d = abs(matrixplane[i][0] * round(surface[j][1], 8) +
                    matrixplane[i][1] * surface[j][2] +
                    matrixplane[i][2] * surface[j][3] +
                    matrixplane[i][3]) / (pow((pow(matrixplane[i][0], 2) +
                                        pow(matrixplane[i][1], 2) +
                                        pow(matrixplane[i][2], 2)), 0.5))
            if d < 0.3*rAT:
                n += 1
                quale = j
                quali.append(quale)
        quali1.append(quali)
        quali = []
        numero.append(n)
        n = 0

    return numero, quali1

def get_thickness(coordinates, cutoff1):
    r1 = 0.
    r2 = 0.
    surface, radius, count1, \
    bulk, generalcount = get_which_surface(coordinates, cutoff1)
    surfacepeel, radius3, \
    generalcount3, bulk4 = get_which_surface_peeling(bulk, cutoff)
    r1 = min(radius) - min(radius3)
    r2 = max(radius) - max(radius3)
    return r1, r2

#Input: coordinates and cutoff
#Output: the faceting ratio
def get_FEratio(coordinates, cutoff1):
    ratio1 = 0.
    surface, radius, count1, \
    bulk, generalcount = get_which_surface(coordinates, cutoff1)
    surfacepeel, radius3, \
    generalcount3, bulk4 = get_which_surface_peeling(bulk, cutoff)
    ratio1 = (len(surface)/len(coordinates))*100
    return ratio1

