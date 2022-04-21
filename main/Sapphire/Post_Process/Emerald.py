import math
import time
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import functools
import operator
import multiprocessing as mp
from ase.calculators.emt import EMT
from ase.io import read, write
from scipy.stats import norm
import sys
from math import sqrt, exp, log
from ase.data import chemical_symbols, atomic_numbers
from ase.units import Bohr
from ase.neighborlist import NeighborList
from ase.calculators.calculator import (Calculator, all_changes, PropertyNotImplementedError)
from ase import Atom, Atoms
from ase.build import bulk
from ase.calculators.lammpsrun import LAMMPS
import sys
from importlib import import_module
from ase.io import read

cutoff = 3.5 #first mininum of PDDF Å
rAtEl = 1.442 #raggio atomico oro Å
nearst = 2.882 #Å
rC1 = 4.09 #Å
rC2 = 10.50 #Å
epsilon = 3.81 #eV
rws = 1.593 #Å
#constant for energy Au
p = 10.229
q = 4.036
A = 0.2061 #eV
E = 1.790 #eV

#Input: file xyz
#Output: matrix of i atoms in row with columns of symbol, x, y, z,
#number of atoms and title
def read_file(filename):
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

#Input: coordinates, cutoff, nearest neighbor distance,
#p, q, A, E that are defined on the top
#Output: total energy
def get_energy(coordinates, nearst1, p, q, A, E):
    energy = 0.
    energy1 = 0.
    energyb = []
    energyr = []
    energytot = 0.
    for i in range(len(coordinates)):
        for j in range(len(coordinates)):
            if i != j:
                energy += pow(E, 2)*\
                          exp(-2*q*((get_distance(coordinates, i, j) / nearst1)-1))
                energy1 += A*\
                           exp(-p*((get_distance(coordinates, i, j) / nearst1)-1))
        energyb.append(-pow(energy, 0.5))
        energyr.append(energy1)
        energy = 0.
        energy1 = 0.
    for i in range(len(energyb)):
        energytot += energyb[i] + energyr[i]
    return energytot

def get_energy_surface(coordinates, cutoff1, nearst1, p, q, A, E):
    surface, radius, count1,\
    bulk, generalcount = get_which_surface(coordinates, cutoff1)
    count, quali = get_coordination_number(coordinates, cutoff1)
    surfacepeel, radius3, \
    generalcount3, bulk4 = get_which_surface_peeling(bulk, cutoff)
    energy = 0.
    energy1 = 0.
    energyb = []
    energyr = []
    energytot = 0.
    subcoordinates = []
    subcoordinates2 = []
    for i in range(len(surface)):
        subcoordinates.append(surface[i])
    for i in range(len(surfacepeel)):
        subcoordinates.append(surfacepeel[i])
    for i in range(len(surface)):
        for j in range(len(subcoordinates)):
            if i != j:
                energy += pow(E, 2)*\
                          exp(-2*q*((get_distance(subcoordinates, i, j) / nearst1)-1))
                energy1 += A*\
                           exp(-p*((get_distance(subcoordinates, i, j) / nearst1)-1))
        energyb.append(-pow(energy, 0.5))
        energyr.append(energy1)
        energy = 0.
        energy1 = 0.
    for i in range(len(energyb)):
        energytot += energyb[i] + energyr[i]
    return energytot

#Input: coordinates, cutoff, nearest neighbor distance,
#p, q, A, E that are defined on the top
#Output: total energy
def get_energy_correct(coordinates, nearst1, rC1, rC2, p, q, A, E):
    energy = 0.
    energy1 = 0.
    energy2 = 0.
    energyb = []
    energyr = []
    energytot = 0.
    ar = -A*exp(-p*(rC1/nearst1-1))/(pow(rC2-rC1, 3))
    br = -(p/nearst1)*A*exp(-p*(rC1/nearst1-1))/(pow(rC2-rC1, 2))
    cr = -pow((p/nearst1), 2)*A*exp(-p*(rC1/nearst1-1))/(rC2-rC1)
    a3 = (20.*ar-8.*br+cr)/2.
    a4 = (15.*ar-7.*br+cr)/(rC2-rC1)
    a5 = (12.*ar-6.*br+cr)/(2.*pow(rC2-rC1, 2))
    for i in range(len(coordinates)):
        for j in range(len(coordinates)):
            if i != j:
                if get_distance(coordinates, i, j) < rC1:
                    energy += pow(E, 2)*exp(-2*q*((get_distance(coordinates, i, j) / nearst1)-1))
                    energy1 += A*exp(-p*((get_distance(coordinates, i, j) / nearst1)-1))
                if get_distance(coordinates, i, j) < rC2 and \
                        get_distance(coordinates, i, j) > rC1:
                    energy2 += a3*pow(get_distance(coordinates, i, j) - rC2, 3) + \
                               a4*pow(get_distance(coordinates, i, j) - rC2, 4) + \
                               a5*pow(get_distance(coordinates, i, j) - rC2, 5)
        energyb.append(-pow(energy, 0.5))
        energyr.append(energy1)
        energy = 0.
        energy1 = 0.
    for i in range(len(energyb)):
        energytot += energyb[i] + energyr[i]
    return energytot + energy2

#Input: coordinates, cutoff, nearest neighbor distance,
#p, q, A, E that are defined on the top
#Output: total energy
def get_energy_surface_correct(coordinates, cutoff1, nearst1, rC1, rC2, p, q, A, E):
    surface, radius, count1, bulk, gcn = get_which_surface(coordinates, cutoff1)
    count, quali = get_coordination_number(coordinates, cutoff1)
    surfacepeel, radius3, generalcount3, bulk4 = get_which_surface_peeling(bulk, cutoff)
    energy = 0.
    energy1 = 0.
    energy2 = 0.
    energyb = []
    energyr = []
    energytot = 0.
    subcoordinates = []
    subcoordinates2 = []
    ar = -A * exp(-p * (rC1 / nearst1 - 1)) / (pow(rC2 - rC1, 3))
    br = -(p / nearst1) * A * exp(-p * (rC1 / nearst1 - 1)) / (pow(rC2 - rC1, 2))
    cr = -pow((p / nearst1), 2) * A * exp(-p * (rC1 / nearst1 - 1)) / (rC2 - rC1)
    a3 = (20. * ar - 8. * br + cr) / 2.
    a4 = (15. * ar - 7. * br + cr) / (rC2 - rC1)
    a5 = (12. * ar - 6. * br + cr) / (2. * pow(rC2 - rC1, 2))
    for i in range(len(surface)):
        subcoordinates.append(surface[i])
    for i in range(len(surfacepeel)):
        subcoordinates.append(surfacepeel[i])
    for i in range(len(surface)):
        for j in range(len(subcoordinates)):
            if i != j:
                if get_distance(coordinates, i, j) < rC1:
                    energy += pow(E, 2)*exp(-2*q*((get_distance(coordinates, i, j) / nearst1)-1))
                    energy1 += A*exp(-p*((get_distance(coordinates, i, j) / nearst1)-1))
                if get_distance(coordinates, i, j) < rC2 and \
                        get_distance(coordinates, i, j) > rC1:
                    energy2 += a3*pow(get_distance(coordinates, i, j) - rC2, 3) + \
                               a4*pow(get_distance(coordinates, i, j) - rC2, 4) + \
                               a5*pow(get_distance(coordinates, i, j) - rC2, 5)
        energyb.append(-pow(energy, 0.5))
        energyr.append(energy1)
        energy = 0.
        energy1 = 0.
    for i in range(len(energyb)):
        energytot += energyb[i] + energyr[i]
    return energytot + energy2

#Input: coordinates, cutoff, energy bulk const
#Output: total energy
def get_bulk_energy(coordinates, cutoff1, epsilon1):
    surface, radius, count1,\
    bulk, generalcount = get_which_surface(coordinates, cutoff1)
    N = len(bulk)
    N1 = len(surface)
    return -N*epsilon1, -N1*epsilon1

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



start0 = time.time()
start10 = time.time()
coordinates, n_atoms, title = read_file("Au887_MDh453_rlx.xyz")
end10 = time.time()

start11 = time.time()
xcm2, ycm2, zcm2 = get_coordinatesCM(coordinates, n_atoms)
end11 = time.time()

start12 = time.time()
coordinatesCM = riscale_coordination(coordinates, xcm2, ycm2, zcm2)
end12 = time.time()

start13 = time.time()
radius = get_radius(coordinatesCM)
end13 = time.time()

start = time.time()
count, quali = get_coordination_number(coordinatesCM, cutoff)
end = time.time()

start14 = time.time()
generalcount = get_generalized_CN(coordinatesCM, cutoff)
end14 = time.time()

start1 = time.time()
solidangle = get_solid_angol(coordinatesCM, cutoff)
end1 = time.time()

start2 = time.time()
surface = get_surface(coordinatesCM, cutoff, rAtEl)
end2 = time.time()

start15 = time.time()
surfacediffmax, surfacediffmin = get_diff_surface(coordinatesCM, cutoff, rAtEl)
end15 = time.time()

#volume = get_volume(coordinatesCM, cutoff, rAtEl)

#volumediffmax, volumediffmin = get_diff_volume(coordinatesCM, cutoff, rAtEl)
start3 = time.time()
surface2, radius2, generalcount2, bulk1, generalcount1 = get_which_surface(coordinatesCM, cutoff)
end3 = time.time()

count2, quali2 = get_coordination_number(bulk1, cutoff)

start4 = time.time()
volume = get_volume(coordinatesCM, rAtEl, cutoff)
end4 = time.time()

volumemin = (4/3)*math.pi*pow(min(radius2), 3)

volumemax = (4/3)*math.pi*pow(max(radius2), 3)

volumecirca = get_volume_WS(coordinatesCM, rws)

generalcountpeeling = get_generalized_CN(bulk1, cutoff)

angolsolidpeeling = get_solid_angol(bulk1, cutoff)

radmedium, deltaradmedium = get_medium_radius_surface(coordinatesCM, cutoff)

start16 = time.time()
surfacepeel, radius3, generalcount3, bulk4 = get_which_surface_peeling(bulk1, cutoff)
end16 = time.time()

radius4 = get_radius(bulk4)

radmediumpeel, deltaradmediumpeel = get_medium_radius_peeledsurface(coordinatesCM, cutoff)

#def bins for histo
distances = Euc_Dist(coordinates)

bins = int(round(200/(1+20*np.exp(-len(distances)/1000))))

start5 = time.time()
bin_cent, a = get_PDDF(coordinates, bins)
end5 = time.time()

start6 = time.time()
bin_cent_bulk, a_bulk = get_PDDF(bulk1, bins)
end6 = time.time()
start7 = time.time()
bin_cent_surface, a_surface = get_PDDF(surface2, bins)
end7 = time.time()
start8 = time.time()
bin_cent_surface_peel, a_surface_peel = get_PDDF(surfacepeel, bins)
end8 = time.time()

#energy = get_EMT(coordinatesCM)

#energysurface = get_EMT_surface(coordinatesCM, cutoff)
start9 = time.time()
energy = get_energy(coordinatesCM, nearst, p, q, A, E)
end9 = time.time()

energy2 = get_energy_correct(coordinatesCM, nearst, rC1, rC2, p, q, A, E)

energy3 = get_energy_surface_correct(coordinatesCM, cutoff, nearst, rC1, rC2, p, q, A, E)

start17 = time.time()
energysurface = get_energy_surface(coordinatesCM, cutoff, nearst, p, q, A, E)
end17 = time.time()

energybulk = energy - energysurface

energybulkcorrect = energy2 - energy3

energysommabulk, energysommasurface = get_bulk_energy(coordinatesCM, cutoff, epsilon)

#matrix = get_plane1(coordinatesCM, cutoff)

#numero, quali4 = get_count_atom_plane1(coordinatesCM, cutoff, rAtEl)

#rpdf2, countrdf2 = get_RDF2(coordinatesCM)

start18 = time.time()
cn, number2 = number_count(coordinatesCM, cutoff)
end18 = time.time()

start19 = time.time()
cn1, number3 = number_count_total(coordinatesCM, cutoff)
end19 = time.time()

start20 = time.time()
gcn, number1 = number_generalcount(coordinatesCM, cutoff)
end20 = time.time()
#gcnforplane = get_gcn_for_planes(coordinatesCM, cutoff, rAtEl)

start21 = time.time()
coloredGCN, d3 = colored_GCN(coordinatesCM, cutoff)
end21 = time.time()
start22 = time.time()
coloredGCN2 = colored_GCN2(coordinatesCM, cutoff)
end22 = time.time()
#gcnpiani, uno, dd = get_plane_colored_GCN(coordinatesCM, cutoff)

#quantinumero, qualinumero = get_GCN_atom_plane1(coordinatesCM, cutoff, rAtEl)

thicknessmin, thicknessmax = get_thickness(coordinatesCM, cutoff)

ratio1 = get_FEratio(coordinatesCM, cutoff)

thickness = (thicknessmax + thicknessmin)/2

end0 = time.time()
print("Xcm: %s" % xcm2)
print("Ycm: %s" % ycm2)
print("Zcm: %s" % zcm2)

print("Radius Max: %s" % max(radius2))
print("Radius Min: %s" % min(radius2))

print("Surface: %s" % surface)
print("Surface Max: %s" % surfacediffmax)
print("Surface Min: %s" % surfacediffmin)

print("Medium Radius: %s" % radmedium)
print("Dev Standard: %s" % deltaradmedium)

print("Medium Radius: %s" % radmediumpeel)
print("Dev Standard: %s" % deltaradmediumpeel)

print("Volume: %s" % volume)
print("Volume massimo: %s" % volumemax)
print("Volume minimo: %s" % volumemin)
print("Volume WS: %s" % volumecirca)

print("Energia: %s" % energy)

print("Energia superficie : %s" % energysurface)

print("Energia bulk : %s" % energybulk)

print("EnergiaCorrect: %s" % energy2)

print("Energia surface Correct: %s" % energy3)

print("Energia bulk Correct: %s" % energybulkcorrect)

print("Energia come somma core: %s" % energysommabulk)

print("Energia come somma surface: %s" % energysommasurface)

print("Time for CN: %s" % (end-start))
print("Time for Solid Angle: %s" % (end1-start1))
print("Time for Surface: %s" % (end2-start2))
print("Time for Which surface: %s" % (end3-start3))
print("Time for Volume: %s" % (end4-start4))
print("Time for PDDF total: %s" % (end5-start5))
print("Time for PDDF surface: %s" % (end6-start6))
print("Time for PDDF subsurface: %s" % (end7-start7))
print("Time for PDDF core: %s" % (end8-start8))

print("Time for PDDF core: %s" % (end10-start10))
print("Time for PDDF core: %s" % (end12-start12))
print("Time for PDDF core: %s" % (end13-start13))
print("Time for PDDF core: %s" % (end14-start14))
print("Time for PDDF core: %s" % (end15-start15))

print("Thickness min: %s" %thicknessmin)
print("Thickness max: %s" %thicknessmax)

print("Thickness: %s" %thickness)

print("FE ratio surface-bulk: %s" %ratio1)

print("Time for Code: %s" % (end0-start0))

print(d3)

time1 = [60.7, 179.8, 78.4, 70.5, 27.4, 1237.5, 62.0, 4259]
timex1 = []
timex = []
for i in range(1, 9):
    timex1.append(i)
n, bins, patches = plt.hist(time1, 8, density=True, facecolor='b', alpha=0.75)
plt.xlabel('Clusters')
plt.ylabel('Time Coding [s]')
plt.xlim(0, 8)
plt.ylim(0, 4270)
plt.show()

fig, ax12 = plt.subplots()
ax12.bar(timex1, time1, width=0.2)
ax12.set(xlim=(0, 9),
       xticks=np.arange(0, 9, 1),
       ylim=(0, 1400), yticks=np.arange(0, 1400, 200))
plt.xlabel('Clusters')
plt.ylabel('Time Coding [s]')
plt.show()

time = [end-start, end1-start1, end2-start2, end3-start3,
             end4-start4, end5-start5, end6-start6, end7-start7,
             end8-start8, end9-start9, end10-start10, end11-start11, end12-start12, end13-start13,
             end14-start14, end15-start15, end16-start16, end17-start17,
             end18-start18, end19-start19, end20-start20,
             end21-start21, end22-start22]
timex = []
for i in range(1, 24):
    timex.append(i)

fig23, ax23 = plt.subplots()
ax23.bar(timex, time, width=0.4, edgecolor='white', linewidth=0.7)
ax23.set(xlim=(0, 25),
       xticks=np.arange(0, 25, 1),
       ylim=(0, 15), yticks=np.arange(0, 15, 1))
plt.xlabel('Code number')
plt.ylabel('Time [s]')
plt.show()
plt.savefig('time.png')

#SURFACE HISTO
perhistox, perhistoy = get_histo(surface2, 0.2)
perhistorad = get_radius(surface2)
fig, ax = plt.subplots()
ax.bar(perhistox, perhistoy, width=0.2, edgecolor='white', linewidth=0.7)
ax.set(xlim=(round(min(perhistorad), 1)-1, round(max(perhistorad), 1)+1),
       xticks=np.arange(round(min(perhistorad), 1) - 1,
                        round(max(perhistorad), 1) + 1, 0.2),
       ylim=(0, max(perhistoy)+4), yticks=np.arange(0, max(perhistoy)+4, 4))
plt.xlabel('Radius [Å]')
plt.ylabel('Number')
plt.show()
plt.savefig('historadiussurface1212.png')

#PEELED HISTO
perhistox1, perhistoy1 = get_histo(surfacepeel, 0.2)
perhistorad1 = get_radius(surfacepeel)
fig1, ax = plt.subplots()
ax.bar(perhistox1, perhistoy1, width=0.2, edgecolor='white', linewidth=0.7)
ax.set(xlim=(round(min(perhistorad1), 1)-1, round(max(perhistorad1), 1)+1),
       xticks=np.arange(round(min(perhistorad1), 1) - 1,
                        round(max(perhistorad1), 1) + 1, 0.2),
       ylim=(0, max(perhistoy1)+4), yticks=np.arange(0, max(perhistoy1)+4, 4))
plt.xlabel('Radius [Å]')
plt.ylabel('Number')
plt.show()
plt.savefig('historadiussurfacepeeled1212.png')

#CORE HISTO
perhistox2, perhistoy2 = get_histo(bulk4, 0.2)
perhistorad2 = get_radius(bulk4)
fig2, ax = plt.subplots()
ax.bar(perhistox2, perhistoy2, width=0.2, edgecolor='white', linewidth=0.7)
ax.set(xlim=(round(min(perhistorad2), 1)-1, round(max(perhistorad2), 1)+1),
       xticks=np.arange(round(min(perhistorad2), 1) - 1,
                        round(max(perhistorad2), 1) + 1, 0.5),
       ylim=(0, max(perhistoy2)+4), yticks=np.arange(0, max(perhistoy2)+4, 4))
plt.xlabel('Radius [Å]')
plt.ylabel('Number')
plt.show()
plt.savefig('historadiuscore1212.png')

#TOTALE HISTO
perhistox3, perhistoy3 = get_histo(coordinatesCM, 0.2)
perhistorad2 = get_radius(coordinatesCM)
fig3, ax3 = plt.subplots()
ax3.plot(perhistox3, perhistoy3, color='grey')
plt.xlabel('Radius [Å]')
plt.ylabel('Number')
plt.show()
plt.savefig('historadiustotale1212.png')

#add in legend
plt.figure(1)
fig4, ax3 = plt.subplots()
ax3.bar(perhistox, perhistoy, color='blue', width=0.2, edgecolor='white', linewidth=0.7)
ax3.bar(perhistox1, perhistoy1, color='green', width=0.2, edgecolor='black', linewidth=0.5, align='edge')
ax3.bar(perhistox2, perhistoy2, color='red', width=0.2, edgecolor='white', linewidth=0.5)
ax3.plot(perhistox3, perhistoy3, color='grey')
plt.xlabel('Radius [Å]')
plt.ylabel('Number')
plt.show()
plt.savefig('Rdf.png')

#add in legend
plt.figure(1)
fig50, ax50 = plt.subplots()
ax3.bar(perhistox, perhistoy, color='blue', width=0.2, edgecolor='white', linewidth=0.7)
ax3.bar(perhistox1, perhistoy1, color='green', width=0.2, edgecolor='black', linewidth=0.7)
ax3.bar(perhistox2, perhistoy2, color='red', width=0.2, edgecolor='white', linewidth=0.7)
ax3.plot()
plt.xlabel('Radius [Å]')
plt.ylabel('Number')
plt.show()
plt.savefig('Rdf.png')

#GRAPH PDDF ALL
plt.figure(1)
plt.plot(bin_cent, a, color='black', linestyle="-", label="Totale")
plt.plot(bin_cent_bulk, a_bulk, color='red', linestyle="-", label="Bulk")
plt.plot(bin_cent_surface, a_surface, color='blue', linestyle="-", label="Surface only totale")
plt.plot(bin_cent_surface_peel, a_surface_peel, color='green', linestyle="-", label="Surface peeled")
plt.axvline(3.5, 0, 2000, ls='--', color='grey')
plt.xlabel('Radius [Å]')
plt.ylabel('PDDF')
plt.show()
plt.savefig("pdftutti.png")


#GRAPH PDDF ALL ZOOMED
plt.figure(1)
plt.plot(bin_cent, a, color='black', linestyle="-", label="Totale")
plt.plot(bin_cent_bulk, a_bulk, color='red', linestyle="-", label="Bulk")
plt.plot(bin_cent_surface, a_surface, color='blue', linestyle="-", label="Surface only totale")
plt.plot(bin_cent_surface_peel, a_surface_peel, color='green', linestyle="-", label="Surface peeled")
plt.axvline(3.5, 0, 2000, ls='--', color='grey')
plt.xlim(2.0, 10.0)
plt.ylim(0.0, 8000)
plt.xlabel('Radius [Å]')
plt.ylabel('PDDF')
plt.show()
plt.savefig('pdftuttizoom.png')

#GRAPH PDDF SURFACE
plt.figure(1)
plt.plot(bin_cent, a)
plt.xlabel('Radius [Å]')
plt.ylabel('PDDF')
plt.show()
plt.savefig('pdftotale.png')


#GRAPH PDDF SURFACE ZOOMED
plt.figure(1)
plt.plot(bin_cent, a)
plt.xlim(2, 10.0)
plt.ylim(0.0, 8000)
plt.axvline(3.5, 0, 2000, ls='--', color='grey')
plt.xlabel('Radius [Å]')
plt.ylabel('PDDF')
plt.show()
plt.savefig('pdftotalezoom.png')


plt.figure(1)
plt.plot(bin_cent_bulk, a_bulk)
plt.xlabel('Radius [Å]')
plt.ylabel('PDDF')
plt.show()
plt.savefig('pdfbulk.png')


plt.figure(1)
plt.plot(bin_cent_surface, a_surface)
plt.xlabel('Radius [Å]')
plt.ylabel('PDDF')
plt.show()
plt.savefig('pdfsurface.png')


plt.figure(1)
plt.plot(bin_cent_surface_peel, a_surface_peel)
plt.xlabel('Radius [Å]')
plt.ylabel('PDDF')
plt.ylim(0.0, 900)
plt.show()
plt.savefig('pdfpeeled.png')

plt.figure(1)
fig5, ax5 = plt.subplots()
ax5.bar(cn, number2, width=0.2, edgecolor='white', linewidth=0.7)
ax5.set(xlim=(round(min(cn), 1)-1, round(max(cn), 1)+1),
       xticks=np.arange(round(min(cn), 1)-1, round(max(cn), 1)+1, 0.5),
       ylim=(0, max(number2)+4), yticks=np.arange(0, max(number2)+4, 10))
plt.xlabel('CN')
plt.ylabel('Number')
plt.show()
plt.savefig('histoCN.png')

plt.figure(1)
fig6, ax6 = plt.subplots()
ax6.bar(gcn, number1, width=0.1, edgecolor='white', linewidth=0.7)
ax6.set(xlim=(round(min(gcn), 1)-1, round(max(gcn), 1)+1),
       xticks=np.arange(round(min(gcn), 1)-1, round(max(gcn), 1)+1, 0.5),
       ylim=(0, max(number1)+4), yticks=np.arange(0, max(number1)+4, 10))
plt.xlabel('GCN')
plt.ylabel('Number')
plt.show()
plt.savefig('histoGCN.png')

plt.figure(1)
fig20, ax20 = plt.subplots()
ax20.bar(cn1, number3, width=0.2, edgecolor='white', linewidth=0.7)
ax20.set(xlim=(3, 13),
       xticks=np.arange(3, 13, 1),
       ylim=(0, max(number3)+4), yticks=np.arange(0, max(number3)+20, 30))
plt.xlabel('CN')
plt.ylabel('Number')
plt.show()
plt.savefig('histoCNtotale.png')

temp1 = []
temp2 = []

coordinatesCM2 = []
coordinatesSf = []
coordinatesSf2 = []
coordinatesSf3 = []
coordinatesSf4 = []
coordinatesSf5 = []
coordinatesSf6 = []
coordinatesSf90 = []
angolsurface = []
peeling = []
radiuspeel = []

'''
#add in legend
plt.figure(1)
plt.plot(rpdf2, countrdf2)
plt.xlabel('Radius [A]')
plt.ylabel('RDF')
plt.show()
plt.savefig('Rdf.png')
'''
for i in range(len(coordinatesCM)):
    if count[i] <= 9:
        coordinatesSf.append(
            [str(coordinatesCM[i][0]), coordinatesCM[i][1],
             coordinatesCM[i][2], coordinatesCM[i][3],
             count[i], generalcount[i],
             radius[i], solidangle[i]])

for i in range(len(coordinatesCM)):
    if count[i] <= 10:
        coordinatesSf2.append(
            [str(coordinatesCM[i][0]), coordinatesCM[i][1],
             coordinatesCM[i][2], coordinatesCM[i][3],
             count[i], generalcount[i],
             radius[i], solidangle[i]])
    else:
        coordinatesSf3.append(
            [str(coordinatesCM[i][0]), coordinatesCM[i][1],
             coordinatesCM[i][2], coordinatesCM[i][3],
             count[i], generalcount[i],
             radius[i], solidangle[i]])

for i in range(len(coordinatesCM)):
    if count[i] <= 11:
        coordinatesSf90.append(
            [str(coordinatesCM[i][0]), coordinatesCM[i][1],
             coordinatesCM[i][2], coordinatesCM[i][3],
             count[i], generalcount[i],
             radius[i], solidangle[i]])

for i in range(len(coordinatesCM)):
    if count[i] <= 10 or solidangle[i] < 0.155 and generalcount[i] < 10.3:
        coordinatesSf4.append(
            [str(coordinatesCM[i][0]), coordinatesCM[i][1],
             coordinatesCM[i][2], coordinatesCM[i][3],
             count[i], generalcount[i],
             radius[i], solidangle[i]])
    else:
        coordinatesSf5.append(
            [str(coordinatesCM[i][0]), coordinatesCM[i][1],
             coordinatesCM[i][2], coordinatesCM[i][3],
             count[i], generalcount[i],
             radius[i], solidangle[i]])

for i in range(len(coordinatesSf2)):
    if coordinatesSf2[i][5] < 6.5:
        temp1.append(
            [str(coordinatesSf2[i][0]), coordinatesSf2[i][1],
             coordinatesSf2[i][2], coordinatesSf2[i][3],
             coordinatesSf2[i][4], coordinatesSf2[i][5],
             coordinatesSf2[i][6], coordinatesSf2[i][7]])

for i in range(len(coordinatesCM2)):
    if coordinatesCM2[i][5] >= 6.5:
        temp2.append(
            [str(coordinatesCM2[i][0]), coordinatesCM2[i][1],
             coordinatesCM2[i][2], coordinatesCM2[i][3],
             coordinatesCM2[i][4], coordinatesCM2[i][5],
             coordinatesCM2[i][6], coordinatesCM2[i][7]])

for i in range(len(coordinatesCM)):
    if count[i] <= 10:
        coordinatesSf6.append(
            [str(coordinatesCM[i][0]), coordinatesCM[i][1],
             coordinatesCM[i][2], coordinatesCM[i][3],
             count[i], generalcount[i],
             radius[i], solidangle[i]])


for i in range(len(bulk1)):
    if count2[i] <= 10:
        peeling.append(
            [str(bulk1[i][0]), bulk1[i][1],
             bulk1[i][2], bulk1[i][3],
             count2[i], generalcountpeeling[i],
             radius[i], angolsolidpeeling[i]])
        radiuspeel.append(radius[i])

with open("fileverticiedge.xyz", 'w') as xyz_file:
    xyz_file.write("%d\n%s" % (len(temp1), title))
    for i in range(len(temp1)):
        xyz_file.write(
            "{:4} {:11.6f} {:11.6f} {:11.6f} {:11.6f} {:11.6f} {:11.6f} {:11.6f}\n".format(
                temp1[i][0], temp1[i][1], temp1[i][2],
                temp1[i][3], temp1[i][4], temp1[i][5],
                temp1[i][6], temp1[i][7]))

with open("surfaceCN10.xyz", 'w') as xyz_file:
    xyz_file.write("%d\n%s" % (len(coordinatesSf2), title))
    for i in range(len(coordinatesSf2)):
        xyz_file.write(
            "{:4} {:11.6f} {:11.6f} {:11.6f} {:11.6f} {:11.6f} {:11.6f} {:11.6f}\n".format(
                coordinatesSf2[i][0], coordinatesSf2[i][1],
                coordinatesSf2[i][2], coordinatesSf2[i][3],
                coordinatesSf2[i][4], coordinatesSf2[i][5],
                coordinatesSf2[i][6], coordinatesSf2[i][7]))

with open("surfaceCN10Else.xyz", 'w') as xyz_file:
    xyz_file.write("%d\n%s" % (len(coordinatesSf3), title))
    for i in range(len(coordinatesSf3)):
        xyz_file.write(
            "{:4} {:11.6f} {:11.6f} {:11.6f} {:11.6f} {:11.6f} {:11.6f} {:11.6f}\n".format(
                coordinatesSf3[i][0], coordinatesSf3[i][1],
                coordinatesSf3[i][2], coordinatesSf3[i][3],
                coordinatesSf3[i][4], coordinatesSf3[i][5],
                coordinatesSf3[i][6], coordinatesSf3[i][7]))

with open("surfaceCN10AndElse.xyz", 'w') as xyz_file:
    xyz_file.write("%d\n%s" % (len(coordinatesSf4), title))
    for i in range(len(coordinatesSf4)):
        xyz_file.write(
            "{:4} {:11.6f} {:11.6f} {:11.6f} {:11.6f} {:11.6f} {:11.6f} {:11.6f}\n".format(
                coordinatesSf4[i][0], coordinatesSf4[i][1],
                coordinatesSf4[i][2], coordinatesSf4[i][3],
                coordinatesSf4[i][4], coordinatesSf4[i][5],
                coordinatesSf4[i][6], coordinatesSf4[i][7]))

with open("surfaceCN10AndElseElse.xyz", 'w') as xyz_file:
    xyz_file.write("%d\n%s" % (len(coordinatesSf5), title))
    for i in range(len(coordinatesSf5)):
        xyz_file.write(
            "{:4} {:11.6f} {:11.6f} {:11.6f} {:11.6f} {:11.6f} {:11.6f} {:11.6f}\n".format(
                coordinatesSf5[i][0], coordinatesSf5[i][1],
                coordinatesSf5[i][2], coordinatesSf5[i][3],
                coordinatesSf5[i][4], coordinatesSf5[i][5],
                coordinatesSf5[i][6], coordinatesSf5[i][7]))

with open("surfaceCN9.xyz", 'w') as xyz_file:
    xyz_file.write("%d\n%s" % (len(coordinatesSf), title))
    for i in range(len(coordinatesSf)):
        xyz_file.write(
            "{:4} {:11.6f} {:11.6f} {:11.6f} {:11.6f} {:11.6f} {:11.6f} {:11.6f}\n".format(
                coordinatesSf[i][0], coordinatesSf[i][1],
                coordinatesSf[i][2], coordinatesSf[i][3],
                coordinatesSf[i][4], coordinatesSf[i][5],
                coordinatesSf[i][6], coordinatesSf[i][7]))

with open("Au887_MDh453_rlxAdd.xyz", 'w') as xyz_file:
    xyz_file.write("%d\n%s" % (len(coordinatesCM), title))
    for i in range(len(coordinatesCM)):
        xyz_file.write(
            "{:4} {:11.6f} {:11.6f} {:11.6f} {:11.6f} {:11.6f} {:11.6f} {:11.6f}\n".format(
                coordinatesCM[i][0], coordinatesCM[i][1],
                coordinatesCM[i][2], coordinatesCM[i][3],
                int(count[i]), generalcount[i],
                radius[i], solidangle[i]))

with open("surfaceSA.xyz", 'w') as xyz_file:
    xyz_file.write("%d\n%s" % (len(angolsurface), title))
    for i in range(len(angolsurface)):
        xyz_file.write(
            "{:4} {:11.6f} {:11.6f} {:11.6f} {:11.6f} {:11.6f} {:11.6f} {:11.6f} {:11.6f}\n".format(
                angolsurface[i][0], angolsurface[i][1],
                angolsurface[i][2], angolsurface[i][3],
                angolsurface[i][4], angolsurface[i][5],
                angolsurface[i][6], angolsurface[i][7],
                angolsurface[i][8]))

with open("surfaceCN10AndElse2.xyz", 'w') as xyz_file:
    xyz_file.write("%d\n%s" % (len(coordinatesSf6), title))
    for i in range(len(coordinatesSf6)):
        xyz_file.write(
            "{:4} {:11.6f} {:11.6f} {:11.6f} {:11.6f} {:11.6f} {:11.6f} {:11.6f}\n".format(
                coordinatesSf6[i][0], coordinatesSf6[i][1],
                coordinatesSf6[i][2], coordinatesSf6[i][3],
                coordinatesSf6[i][4], coordinatesSf6[i][5],
                coordinatesSf6[i][6], coordinatesSf6[i][7]))

with open("surfaceCN10AndElseElse2.xyz", 'w') as xyz_file:
    xyz_file.write("%d\n%s" % (len(bulk1), title))
    for i in range(len(bulk1)):
        xyz_file.write(
            "{:4} {:11.6f} {:11.6f} {:11.6f} {:11.6f} {:11.6f} {:11.6f} {:11.6f}\n".format(
                bulk1[i][0], bulk1[i][1],
                bulk1[i][2], bulk1[i][3],
                bulk1[i][4], bulk1[i][5],
                bulk1[i][6], bulk1[i][7]))

with open("peelingCN10.xyz", 'w') as xyz_file:
    xyz_file.write("%d\n%s" % (len(peeling), title))
    for i in range(len(peeling)):
        xyz_file.write(
            "{:4} {:11.6f} {:11.6f} {:11.6f} {:11.6f} {:11.6f} {:11.6f} {:11.6f}\n".format(
                peeling[i][0], peeling[i][1],
                peeling[i][2], peeling[i][3],
                peeling[i][4], peeling[i][5],
                peeling[i][6], peeling[i][7]))

with open("Whichsurface.xyz", 'w') as xyz_file:
    xyz_file.write("%d\n%s" % (len(surface2), title))
    for i in range(len(surface2)):
        xyz_file.write(
            "{:4} {:11.6f} {:11.6f} {:11.6f} {:11.6f} {:11.6f} {:11.6f} {:11.6f}\n".format(
                surface2[i][0], surface2[i][1],
                surface2[i][2], surface2[i][3],
                surface2[i][4], surface2[i][5],
                surface2[i][6], surface2[i][7]))

with open("Whichbulk.xyz", 'w') as xyz_file:
    xyz_file.write("%d\n%s" % (len(bulk1), title))
    for i in range(len(bulk1)):
        xyz_file.write(
            "{:4} {:11.6f} {:11.6f} {:11.6f} {:11.6f} {:11.6f} {:11.6f} {:11.6f}\n".format(
                bulk1[i][0], bulk1[i][1],
                bulk1[i][2], bulk1[i][3],
                bulk1[i][4], bulk1[i][5],
                bulk1[i][6], bulk1[i][7]))

with open("Whichpeel.xyz", 'w') as xyz_file:
    xyz_file.write("%d\n%s" % (len(surfacepeel), title))
    for i in range(len(surfacepeel)):
        xyz_file.write(
            "{:4} {:11.6f} {:11.6f} {:11.6f} {:11.6f} {:11.6f} {:11.6f} {:11.6f}\n".format(
                surfacepeel[i][0], surfacepeel[i][1],
                surfacepeel[i][2], surfacepeel[i][3],
                surfacepeel[i][4], surfacepeel[i][5],
                surfacepeel[i][6], surfacepeel[i][7]))

with open("surfaceCN11.xyz", 'w') as xyz_file:
    xyz_file.write("%d\n%s" % (len(coordinatesSf90), title))
    for i in range(len(coordinatesSf90)):
        xyz_file.write(
            "{:4} {:11.6f} {:11.6f} {:11.6f} {:11.6f} {:11.6f} {:11.6f} {:11.6f}\n".format(
                coordinatesSf90[i][0], coordinatesSf90[i][1],
                coordinatesSf90[i][2], coordinatesSf90[i][3],
                coordinatesSf90[i][4], coordinatesSf90[i][5],
                coordinatesSf90[i][6], coordinatesSf90[i][7]))

with open("Whichbulkpeeled.xyz", 'w') as xyz_file:
    xyz_file.write("%d\n%s" % (len(bulk4), title))
    for i in range(len(bulk4)):
        xyz_file.write(
            "{:4} {:11.6f} {:11.6f} {:11.6f} {:11.6f} {:11.6f} {:11.6f} {:11.6f}\n".format(
                bulk4[i][0], bulk4[i][1],
                bulk4[i][2], bulk4[i][3],
                bulk4[i][4], bulk4[i][5],
                bulk4[i][6], bulk4[i][7]))
'''
with open("Plane.xyz", 'w') as xyz_file:
    xyz_file.write("%d\n%s" % (len(matrix), title))
    for i in range(len(matrix)):
        xyz_file.write(
            "{:4} {:11.6f} {:11.6f} {:11.6f} {:11.6f} {:11.6f} {:11.6f} {:11.6f}\n".format(
                matrix[i][0], matrix[i][1],
                matrix[i][2], matrix[i][3],
                matrix[i][4], matrix[i][5],
                matrix[i][6], matrix[i][7], numero[i]))

'''
#with open("Plane.xyz", 'w') as xyz_file:
 #   xyz_file.write("%d\n%s" % (len(matrix), title))
  #  for i in range(len(matrix)):
   #     xyz_file.write(
    #        "{:4} {:11.6f} {:11.6f} {:11.6f} {:11.6f}\n".format(
     #           matrix[i][0], matrix[i][1],
      #          matrix[i][2], matrix[i][3], numero[i]))

'''
with open("Plane1.xyz", 'w') as xyz_file:
    xyz_file.write("%d\n%s" % (len(matrix), title))
    for i in range(len(matrix)):
        if numero[i] == 20:
            xyz_file.write(
                "{:4} {:11.6f} {:11.6f} {:11.6f} {:11.6f} {:11.6f}\n".format(
                    matrix[i][0], matrix[i][1],
                    matrix[i][2], matrix[i][3], numero[i], gcnforplane[i]))
'''
ag = []
N1 = len(surface2)
N2 = len(surfacepeel)
for i in range(N1):
    ag.append(str("Ag"))

pt = []
for i in range(N2):
    pt.append(str("Cu"))


with open("Ovito.xyz", 'w') as xyz_file:
    xyz_file.write("%d\n%s" % (len(coordinates), title))
    for i in range(len(surface2)):
        xyz_file.write(
            "{:4} {:11.6f} {:11.6f} {:11.6f}\n".format(
                ag[i], surface2[i][1],
                surface2[i][2], surface2[i][3]))
    for i in range(len(surfacepeel)):
        xyz_file.write(
            "{:4} {:11.6f} {:11.6f} {:11.6f}\n".format(
                pt[i], surfacepeel[i][1],
                surfacepeel[i][2], surfacepeel[i][3]))
    for i in range(len(bulk4)):
        xyz_file.write(
            "{:4} {:11.6f} {:11.6f} {:11.6f}\n".format(
                bulk4[i][0], bulk4[i][1],
                bulk4[i][2], bulk4[i][3]))

with open("OvitoGCN.xyz", 'w') as xyz_file:
    xyz_file.write("%d\n%s" % (len(coloredGCN), title))
    for i in range(len(coloredGCN)):
        xyz_file.write(
            "{:4} {:11.6f} {:11.6f} {:11.6f}\n".format(
                coloredGCN[i][0], coloredGCN[i][1],
                coloredGCN[i][2], coloredGCN[i][3]))

with open("OvitoGCN1.xyz", 'w') as xyz_file:
    xyz_file.write("%d\n%s" % (len(coloredGCN2), title))
    for i in range(len(coloredGCN2)):
        xyz_file.write(
            "{:4} {:11.6f} {:11.6f} {:11.6f}\n".format(
                coloredGCN2[i][0], coloredGCN2[i][1],
                coloredGCN2[i][2], coloredGCN2[i][3]))
'''
with open("Plane1.xyz", 'w') as xyz_file:
    xyz_file.write("%d\n%s" % (len(gcnpiani), title))
    for i in range(len(gcnpiani)):
        if quantinumero[i] >= 20:
            xyz_file.write(
                "{:4} {:11.6f} {:11.6f} {:11.6f} {:11.6f}\n".format(
                    gcnpiani[i][0], gcnpiani[i][1],
                    gcnpiani[i][2], gcnpiani[i][3], quantinumero[i]))

with open("Plane3.xyz", 'w') as xyz_file:
    xyz_file.write("%d\n%s" % (len(uno), title))
    for i in range(len(uno)):
            xyz_file.write(
                "{:4} {:11.6f} {:11.6f} {:11.6f}\n".format(
                    uno[i][0], uno[i][1],
                    uno[i][2], uno[i][3]))

with open("Plane4.xyz", 'w') as xyz_file:
    xyz_file.write("%d\n%s" % (len(dd), title))
    for i in range(len(dd)):
            xyz_file.write(
                "{:4} \n".format(
                    dd[i]))
'''