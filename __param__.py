import numpy as np

coordinates_path="/esnas/releases/models/nemo/v3.6_ecearth/inidata/nemo/"

# Constant definitions
grav = 9.81                  # acceleration due to gravity (m.s-2)
omega = 7.292115083046061e-5 # earth rotation rate (s-1)
earthrad = 6371229           # mean earth radius (m)
deg2rad = np.pi / 180.
rau0=1025.
rcp=3990.

input_chunks_ORCA1L75={
    "time_counter": 12, "t": 12,
    "deptht": 25, "depthu": 25, "depthv": 25, "depthw": 25,
    "y": 145, "x": 60
    }
target_chunks_ORCA1L75={
    "t": 12,
    "z_c": 25, "z_l": 25,
    "y_c": 145, "y_r": 145,
    "x_c": 60, "x_r": 60
    }

update_orca_variables={
    'uo': {'dims': ['t', 'z_c', 'y_c', 'x_r']},
    'u2o': {'dims': ['t', 'z_c', 'y_c', 'x_r']},
    'umo': {'dims': ['t', 'z_c', 'y_c', 'x_r']},
    'vo': {'dims': ['t', 'z_c', 'y_r', 'x_c']},
    'v2o': {'dims': ['t', 'z_c', 'y_r', 'x_c']},
    'vmo': {'dims': ['t', 'z_c', 'y_r', 'x_c']},
    'tauuo': {'dims': ['t', 'y_c', 'x_r']},
    'tauvo': {'dims': ['t', 'y_r', 'x_c']},
    'thkcello' : {'dims': ['t', 'z_c', 'y_c', 'x_c']},
    'pbo' : {'dims': ['t', 'y_c', 'x_c']},
    'zos' : {'dims': ['t', 'y_c', 'x_c']},
    'zossq' : {'dims': ['t', 'y_c', 'x_c']},
    'thetao' : {'dims': ['t', 'z_c', 'y_c', 'x_c']},
    'tos' : {'dims': ['t', 'y_c', 'x_c']},
    'tossq' : {'dims': ['t', 'y_c', 'x_c']},
    'so' : {'dims': ['t', 'z_c', 'y_c', 'x_c']},
    'sos' : {'dims': ['t', 'y_c', 'x_c']},
    'rhopoto' : {'dims': ['t', 'z_c', 'y_c', 'x_c']},
    'mlotst' : {'dims': ['t', 'y_c', 'x_c']},
    'tmaskgin' : {"dims": ["y_c", "x_c"]},
    'umaskgin' : {"dims": ["y_c", "x_r"]},
    'vmaskgin' : {"dims": ["y_r", "x_c"]},
    'fmaskgin' : {"dims": ["y_r", "x_r"]},
    'tmasknna' : {"dims": ["y_c", "x_c"]},
    'umasknna' : {"dims": ["y_c", "x_r"]},
    'vmasknna' : {"dims": ["y_r", "x_c"]},
    'fmasknna' : {"dims": ["y_r", "x_r"]},
    'tmaskmed' : {"dims": ["y_c", "x_c"]},
    'umaskmed' : {"dims": ["y_c", "x_r"]},
    'vmaskmed' : {"dims": ["y_r", "x_c"]},
    'fmaskmed' : {"dims": ["y_r", "x_r"]},
    'tmasklab' : {"dims": ["y_c", "x_c"]},
    'umasklab' : {"dims": ["y_c", "x_r"]},
    'vmasklab' : {"dims": ["y_r", "x_c"]},
    'fmasklab' : {"dims": ["y_r", "x_r"]},
    'tmaskwed' : {"dims": ["y_c", "x_c"]},
    'umaskwed' : {"dims": ["y_c", "x_r"]},
    'vmaskwed' : {"dims": ["y_r", "x_c"]},
    'fmaskwed' : {"dims": ["y_r", "x_r"]},
    'tmaskarc' : {"dims": ["y_c", "x_c"]},
    'umaskarc' : {"dims": ["y_c", "x_r"]},
    'vmaskarc' : {"dims": ["y_r", "x_c"]},
    'fmaskarc' : {"dims": ["y_r", "x_r"]},
    'tmasksoc' : {"dims": ["y_c", "x_c"]},
    'umasksoc' : {"dims": ["y_c", "x_c"]},
    'vmasksoc' : {"dims": ["y_c", "x_c"]},
    'fmasksoc' : {"dims": ["y_c", "x_c"]},
    'tmaskindpac' : {"dims": ["y_c", "x_c"], "old_names": ["indpacmsk", ]},
    'umaskindpac' : {"dims": ["y_c", "x_r"], "old_names": ["indpacmsk", ]},
    'vmaskindpac' : {"dims": ["y_r", "x_c"], "old_names": ["indpacmsk", ]},
    'fmaskindpac' : {"dims": ["y_r", "x_r"], "old_names": ["indpacmsk", ]},
    "e1f": {"dims": ["y_r", "x_r"]},
    "e2f": {"dims": ["y_r", "x_r"]},
    }
