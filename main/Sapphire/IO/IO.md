The purpose of this set of scripts will be to allocate an output environment
for each of the possible output quantites. 

At this point, I'm thinking of breaking everything up into one of four data
types and taking it from there:
    
1 - System info:
    
base_dir, masterkey, movie_file_name, energy_file_name,
extend_xyz*, Homo, Hetero*, Start, End, Step, Skip, 
uniform_pdf*, Band, Species, NFrames, NAtoms, CoMSpace
    
2 - Time series + ( single datum / frame ):
    
R_Cut^, comdist, moi, gyration, stat_radius, surf_area, surf_atoms, 
simtime`, edelta`, meanetot`, temp`, epot`, etot`, edkin`, edelta`, 
concert, colect, CutE^¬, comE¬, gyrationE¬, surf_areaE¬, surf_atomsE¬, 
HeCut@, mix@, elements

3 - Time series - 1 data point per atom per frame (Fuel for extended XYZ):
    
cna_patterns, agcn, nn, cna_patternsE¬, (Hydrostatic Pressure??)
    
4 - Variable sized arrays (SUPER AWKWARD!):
    
rdf^, cna_sigs(Needs ReWorking), adj, pdf^, comdist, comdistE¬, midcomdist¬,
cna_sigsE¬, herdf^@, hepdf^@, hordfE^¬, hopdf^¬, headj@, hoadjE¬, pair_distance,
he_pair_distance@, pair_distanceE¬
    
Symbol Key:

* = Boolean

^ = Computed every 1/Skip frames

` = With Energy File only

¬ = With Homo enabled only

@ = With Hetero enabled only

