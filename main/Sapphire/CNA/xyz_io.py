import pickle
import numpy as np
from itertools import groupby
from collections import namedtuple

PATH_TO_MOVIE = '/home/annte/Master_Project/PP/Structures_to_compute/Au309_Ih.xyz'
PATH_TO_NEW_MOVIE = '/home/annte/Master_Project/PP/test/Au_309_Ih_CNA.xyz'

#DONT TOUCH THIS readMovieFileXYZ !

def readMovieFileXYZ(path_to_movie):
    """
    Reads a LoDiS movie.xyz file and fetches the coordinates
    for each atom for each frame.
    Input:
        path_to_movie: path to movie.xyz
    Returns:
        Named tuple read_movie:
            - read_movie.Frames: list of frames; each is an array of atoms each described by [Atom, x, y, z, Col]
            - read_movie.Headers: list of the movie frames headers

    """
    read_file_chars = []
    with open(path_to_movie, 'r') as file:
        for line in file:
            read_file_chars.append(line)
    # 1. Delete line jump
    read_file_chars = [line[:-1] for line in read_file_chars]
    read_file_chars
    # 2. Separate line by line
    grouped_lines = [([list(group) for k, group in groupby(line,lambda x: x == " ") if not k]) for line in read_file_chars]
    # 3. Concatenate charaters
    joined_string = [[''.join(info_elem) for info_elem in grouped_line] for grouped_line in grouped_lines]
    # 4. Regroup into list of lists. Elements of outerlist are movie frames
    merged_frames = []
    current_frame = []
    for line in joined_string:
        if(line==joined_string[0]):
            if len(current_frame)!=0:
                merged_frames.append(current_frame)
            current_frame=[]
        else:
            current_frame.append(line)
    merged_frames.append(current_frame)
    # 5. Removing second line of header
    movie_headers_all_frames = [frame[0] for frame in merged_frames]
    merged_frames = [frame[1:] for frame in merged_frames]
    # 6. Converting coordinates and pressure to floats
    for frame in merged_frames:
        for line in frame:
            line[1] = float(line[1]) # x coord
            line[2] = float(line[2]) # y coord
            line[3] = float(line[3]) # z coord
    Movie = namedtuple('Movie', 'Frames Headers')
    read_movie = Movie(merged_frames, movie_headers_all_frames)
    return(read_movie)

def XYZ_data_maker(Pattern_Dict, System, PATH_TO_MOVIE, PATH_TO_NEW_MOVIE):

    """ (Armand)

    This function takes a npz file location, finds the xyz with the same name,
    within a different folder, and rewrites it to be used with ovito to
    visualize the CNA patterns found.

    input:

    arrays_filename (str): the location of the npz filenames

    output:

    xyz file within the PATH_TO_NEW_MOVIE location, rewritten to include
    where each unique CNA pattern was found.

    """

    with open(PATH_TO_NEW_MOVIE, 'a+') as newmovie: # Mode chosen: append
        movie = readMovieFileXYZ(PATH_TO_MOVIE)
        NATOM = int(len(movie[0][0]))
        open(PATH_TO_NEW_MOVIE, 'w').close() #Clear old movie pressure file
        
        for frame_num, current_frame in enumerate(movie.Frames):
            print('Analyzing frame: {}/{}'.format(frame_num+1, len(movie.Frames)))
            
            with open(PATH_TO_NEW_MOVIE, 'a+') as newmovie: # Mode chosen: append
                num_lines = sum(1 for line in open(PATH_TO_NEW_MOVIE))
                
                if (num_lines==0): # No newline for first line -- bugs Ovito if there is newline at beginning
                    newmovie.write(str(NATOM)+'\n')
                    
                else:
                    newmovie.write('\n' + str(NATOM)+'\n')
                newmovie.write(System['movie_file_name'][:-4])
                
                for atom_index, atom_info in enumerate(current_frame):
                    if atom_index in np.nonzero(Pattern_Dict[System['movie_file_name']])[0]:
                        k=list(np.nonzero(Pattern_Dict[System['movie_file_name']])[0]).index(atom_index)
                        atom_info[-1] = np.nonzero(Pattern_Dict[System['movie_file_name']])[1][k] + 1
                    
                    else:
                        atom_info[-1] = 0
                    newmovie.write('\n')
                    newmovie.write('  \t'.join(str(item) for item in atom_info))


