# Initializes FREE object. Includes path and label generation.

from datetime import datetime
import sys
import os
import json


if sys.version[0] == '3':  # if using Python 3, we must import the reload module (it is built in to Python 2)
    from importlib import reload

import free_fea.user_input
reload(free_fea.user_input)
from free_fea.user_input import *

import free_fea.geometry_setup
reload(free_fea.geometry_setup)
from free_fea.geometry_setup import *

# initialize FREE object using parameters defined in free_fea.user_input()
def initialize_FREE(label='user_input'):
    if label != 'user_input':  # if label specified, then load input file
        myFREE = load_FREE(label)
    else:  # otherwise, create a new FREE using the parameters in free_fea.user_input
        myFREE = FREE()
        myFREE = apply_user_input(myFREE)

        myFREE.path = generate_FREE_path(myFREE)

        TNB_T, TNB_N, TNB_B = find_TNB_frame(myFREE.path)  # calculate TNB frame for FREE path
        myFREE.path = reorient_path(myFREE.path, TNB_T[0], TNB_N[0], TNB_B[0])  # orient curve to initially point along z axis and curve about the y-axis

        myFREE.length = find_path_length(myFREE.path)  # initial FREE length (mm)

        # check FREE radius - will trigger an error message if curvature is too tight
        check_FREE_curvature(myFREE.path, myFREE.outer_radius)

        myFREE.label = generate_label(myFREE)

        save_FREE(myFREE)  # saves FREE class object in the input folder

    return myFREE


# Generate descriptive label for output files
def generate_label(FREE):
    if FREE.path_type == 'straight':
        label = 'str'
    elif FREE.path_type == 'arc':  # 'arc': [arc_radius(mm), arc_angle(deg)]
        label = 'arc'
    elif FREE.path_type == 'quadratic_Bezier':
        label = 'qbez'
    elif FREE.path_type == 'cubic_Bezier':
        label = 'cbez'
    elif FREE.path_type == 'spiral':
        label = 'spi'
    elif FREE.path_type == 'from_2D_experiment':
        label = '2Dexp'
        label += datetime.now().strftime("%Y%m%d") + "_" + datetime.now().strftime("%H%M")  # experimental results are labeled with the current date and time
    else:
        raise ValueError('FREE path type not recognized.')

    if FREE.path_type != 'from_2D_experiment':
        if FREE.path_parameters.ndim == 1:
            label += "_"
            for i in range(np.shape(FREE.path_parameters)[0]):
                label += str(int(FREE.path_parameters[i]))
                if i != np.shape(FREE.path_parameters)[0] - 1:
                    label += "-"
        elif FREE.path_parameters.ndim == 2:
            for i in range(np.shape(FREE.path_parameters)[0]):
                label += "_"
                for j in range(np.shape(FREE.path_parameters)[1]):
                    label += str(int(FREE.path_parameters[i][j]))
                    if j != np.shape(FREE.path_parameters)[1] - 1:
                        label += "-"
        else:
            raise ValueError('FREE parameters not a recognized shape. Check input.')

    label += ('_r' + str(int(FREE.outer_radius)) +
              '_a' + str(int(FREE.fiber_a_angle_deg)) +
              '_b' + str(int(FREE.fiber_b_angle_deg)))

    if FREE.has_def_fib_angle:
        label += '_DFA'

    if hasattr(FREE, 'label_addon'):
        if FREE.label_addon != '':
            label += '_' + FREE.label_addon

    return label


# Save FREE object to json file
def save_FREE(saveFREE):
    # convert numpy arrays in FREE to lists for JSON compatibility
    if isinstance(saveFREE.path, np.ndarray):
        np_path = saveFREE.path  # store as np array for re-conversion
        saveFREE.path = saveFREE.path.tolist()
    if isinstance(saveFREE.path_parameters, np.ndarray):
        np_path_parameters = saveFREE.path  # store as np array for re-conversion
        saveFREE.path_parameters = saveFREE.path_parameters.tolist()

    # navigate to output folder. Abaqus Work Directory and curved_FREEs must be in the same folder
    current_directory = os.getcwd()  # get current work directory
    parent_directory = os.path.split(current_directory)[0]  # go up one folder
    save_directory = os.path.join(parent_directory, 'curved_FREEs', 'input')  # navigate to curved_FREEs\output\(myFREE.label)
    save_path = os.path.join(save_directory, saveFREE.label + '.json')

    save_file = open(save_path, 'w')
    json.dump(vars(saveFREE), save_file, indent=4)  # converts FREE class instance to dict for storage
    save_file.close()

    # convert back to numpy arrays
    if type(saveFREE.path) == list:
        saveFREE.path = np_path
    if type(saveFREE.path_parameters) == list:
        saveFREE.path_parameters = np_path_parameters


# load npz file containing FREE object
def load_FREE(label):
    # trim filetype from label if included
    if label[-4] == '.':
        label = label[:-4]
    elif label[-5] == '.':
        label = label[:-5]

    # navigate to input folder. Abaqus Work Directory and curved_FREEs must be in the same folder
    current_directory = os.getcwd()  # get current work directory
    parent_directory = os.path.split(current_directory)[0]  # go up one folder
    load_directory = os.path.join(parent_directory, 'curved_FREEs', 'input')  # navigate to curved_FREEs\input
    load_path = os.path.join(load_directory, label + '.json')

    load_file = open(load_path)
    load_dict = json.load(load_file)  # converts FREE class instance to dict for storage
    load_file.close()

    if sys.version[0] == '2':  # convert unicode to strings after loading json - necessary in Python2 (i.e. Abaqus)
        load_dict = byteify(load_dict)

    loadFREE = FREE()
    loadFREE.Dict2Class(load_dict)

    # convert lists back to numpy arrays
    if type(loadFREE.path) == list:
        loadFREE.path = np.array(loadFREE.path)
    if type(loadFREE.path_parameters) == list:
        loadFREE.path_parameters = np.array(loadFREE.path_parameters)

    return loadFREE

# convert unicode to strings after loading json - necessary in Python2 (i.e. Abaqus)
# https://stackoverflow.com/questions/956867/how-to-get-string-objects-instead-of-unicode-from-json
def byteify(input):
    if isinstance(input, dict):
        return {byteify(key): byteify(value)
                for key, value in input.items()}
    elif isinstance(input, list):
        return [byteify(element) for element in input]
    elif isinstance(input, unicode):
        return input.encode('utf-8')
    else:
        return input


# FREE class for storing all FREE parameters
class FREE:

    def __init__(self):
        # empty fields to be populated in free_fea.geometry_setup
        self.has_var_fib_angle = False  # FREE has variable fiber angle
        self.path = np.array([])  # Path will be built using path_parameters.
        self.label = ''  # Label will be generated using path_parameters.
        self.length = 0.0  # length is calculated using path.
        self.label_addon = ''  # customizable label addon

        # FREE parameters
        self.path_type = ''  # FREE path type. See geometry_setup.define_FREE_path.
        self.path_parameters = np.array([])  # FREE path parameters. Dependent on path_type.
        self.material = ''  # see abaqus_setup.define_Abaqus_materials() for options
        self.outer_radius = 0.0 # FREE outside radius (mm)
        self.thickness = 0.0  # FREE thickness (mm)
        self.fiber_a_angle_deg = 0.0  # initial fiber angle for fiber direction A (deg), i.e. alpha
        self.fiber_b_angle_deg = 0.0  # initial fiber angle for fiber direction B (deg), i.e. beta
        self.num_fibers = 0  # number of fibers in each direction
        self.endcap_thickness = 0.0  # endcap thickness (mm)
        self.final_pressure = 0.0  # FREE pressure (MPa)

        # Path parameters
        self.num_path_points = 0  # Number of points along FREE path. Affects path fidelity and output resolution.

        # Abaqus mesh parameters
        self.fiber_seed_size = 0.0  # fiber seed size
        self.bladder_seed_size = 0.0  # bladder seed size
        self.endcap_seed_size = 0.0  # endcap seed size
        self.num_output_nodes_circ = 0  # number of output nodes in bladder mesh circumference
        self.fiber_element = ''  # fiber element type. Use 'beam' when alpha = -beta. Use 'truss' for alpha != -beta.

    # Turns a dictionary into a class. Used for loading FREEs from json files.
    # https://www.geeksforgeeks.org/how-to-change-a-dictionary-into-a-class/
    def Dict2Class(self, my_dict):
        for key in my_dict:
            setattr(self, key, my_dict[key])

