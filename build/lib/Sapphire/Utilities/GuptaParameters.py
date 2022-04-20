"""Potential parameters for the Gupta potential."""

# Parameter sequence: [p, q, A(eV), xi(eV), r0(Angs)]

Ag_parameters = {
    'Ag': [10.85, 3.18, 0.1031, 1.1895, 2.89],
}

Pd_parameters = {
    'Pd': [11.00, 3.794, 0.1715, 1.7019, 2.75],
}

Cu_parameters = {
    'Cu': [10.55, 2.43, 0.0894, 1.2799, 2.56],
}

Ni_parameters = {
    'Ni': [11.34, 2.27, 0.0958, 1.5624, 2.49],
}

Pt_parameters = {
    'Pt': [10.612, 4.004, 0.2975, 2.695, 2.77],
}

Au_parameters = {
    'Au': [10.53, 4.30, 0.2197, 1.855, 2.89],
}

AuCu_parameters = {
    'Au': [10.229, 4.0360, 0.2061, 1.7900, 2.884],
    'Cu': [10.96, 2.2780, 0.0855, 1.224, 2.556],
    ('Au','Cu'): [11.05, 3.047, 0.1539, 2.2, 2.83],
}

AgCu_parameters = {
    'Ag': [10.85, 3.18, 0.1031, 1.1895, 2.89],
    'Cu': [10.55, 2.43, 0.0894, 1.279, 2.725],
    ('Ag','Cu'): [11.05, 3.047, 0.1539, 1.5605, 2.8075],
}


AuPt_parameters = {
    'Au': [10.53, 4.30, 0.2197, 1.855, 2.89],
    'Pt': [10.612, 4.004, 0.2975, 2.695, 2.77],
    ('Au','Pt'): [10.42, 4.02, 0.1539, 1.5605, 2.8075],
}
AgAu_parameters = {
    'Ag': [10.85, 3.18, 0.1031, 1.1895, 2.89],
    'Au': [10.139, 4.03, 0.2097, 1.815, 2.89],
    ('Ag','Au'): [10.494, 3.607, 0.149, 1.4874, 2.87],
}

AgPd_parameters = {
    'Ag': [10.85, 3.18, 0.1031, 1.1895, 2.89],
    'Pd': [11.00, 3.794, 0.1715, 1.7019, 2.75],
    ('Ag','Pd'): [10.895, 3.492, 0.161, 1.5597, 2.82],
}

AgNi_parameters = {
    'Ag': [10.85, 3.18, 0.1031, 1.1895, 2.89],
    'Ni': [11.34, 2.27, 0.0958, 1.5624, 2.49],   
    ('Ag','Ni'): [11.095, 2.725, 0.096, 1.34, 2.69],
}

PdPt_parameters = {
    'Pd': [10.867, 3.742, 0.1746, 1.718, 2.75],
    'Pt': [10.612, 4.004, 0.2975, 2.695, 2.77],
    ('Pd','Pt'): [10.74, 3.87, 0.23, 2.2, 2.76],

}