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