# Save output data after running the FEA job created by 1_Abaqus_FREE_setup.py
# This file can only be run in the Abaqus Python environment (Python 2.7).
# 1_Abaqus_FREE_setup.py MUST be run in Abaqus before this script
# saves output as .npz files in a folder named [FREE.label] in curved_FREEs\output

import os
import sys

from part import *
from material import *
from section import *
from assembly import *
from step import *
from interaction import *
from load import *
from mesh import *
from optimization import *
from job import *
from sketch import *
from visualization import *
from connectorBehavior import *
from odbAccess import *

# import custom functions
# Navigate to curved_FREEs to import custom packages. Abaqus Work Directory and curved_FREEs must be in the same folder
current_directory = os.getcwd()  # get current work directory
parent_directory = os.path.split(current_directory)[0]  # go up one folder
new_path = os.path.join(parent_directory, 'curved_FREEs')  # navigate to curved_FREEs to import custom packages
if all(sys_path != new_path for sys_path in sys.path):
    sys.path.insert(1, new_path)

import free_fea.initialization
reload(free_fea.initialization)
from free_fea.initialization import *

########## USER INPUT ##########
# initialize FREE object
# label='user_input' for a new FREE from free_fea.user_input.
# label='[input file name]' to load a FREE from an input file
myFREE = initialize_FREE(label='user_input')
# myFREE = initialize_FREE(label='arc_30-220_r8_a30_b-30')
################################

odb = openOdb(myFREE.label + '.odb')

# main code block
def main():
    num_frames = len(odb.steps['Step-1'].frames)

    # extract endcap output coordinates
    initial_endcap_nodes = mdb.models['Model-1'].rootAssembly.instances['ENDCAP_OUTPUT-1'].sets['all'].nodes
    endcap_output_set = odb.rootAssembly.instances['ENDCAP_OUTPUT-1'].nodeSets['ALL']
    endcap_coordinates = get_all_endcap_coordinates(initial_endcap_nodes, endcap_output_set, num_frames)

    # extract inner bladder output coordinates
    initial_inner_bladder_nodes = mdb.models['Model-1'].rootAssembly.instances['INNER_BLADDER_OUTPUT-1'].sets['all'].nodes
    inner_bladder_output_set = odb.rootAssembly.instances['INNER_BLADDER_OUTPUT-1'].nodeSets['ALL']
    inner_bladder_output_nodes = odb.rootAssembly.instances['INNER_BLADDER_OUTPUT-1'].nodeSets['ALL'].nodes

    # bladder_coordinates[f][i][j] = point at bladder circumferential index j at s = s[i] in frame f
    inner_bladder_coordinates = get_all_bladder_coordinates(inner_bladder_output_nodes, inner_bladder_output_set, num_frames)

    # extract fiber output coordinates
    fiber_a_initial_nodes = mdb.models['Model-1'].rootAssembly.instances['FIBER_GROUP_A_OUTPUT-1'].sets['all'].nodes
    fiber_a_output_set = odb.rootAssembly.instances['FIBER_GROUP_A_OUTPUT-1'].nodeSets['ALL']
    fiber_a_output_nodes = odb.rootAssembly.instances['FIBER_GROUP_A_OUTPUT-1'].nodeSets['ALL'].nodes
    fiber_b_initial_nodes = mdb.models['Model-1'].rootAssembly.instances['FIBER_GROUP_B_OUTPUT-1'].sets['all'].nodes
    fiber_b_output_set = odb.rootAssembly.instances['FIBER_GROUP_B_OUTPUT-1'].nodeSets['ALL']
    fiber_b_output_nodes = odb.rootAssembly.instances['FIBER_GROUP_B_OUTPUT-1'].nodeSets['ALL'].nodes

    # fiber_x_coordinates[f][i][j] = np.array([x, y, z]) for fiber j at s = s[i] in frame f (group X)
    fiber_a_coordinates = get_all_fiber_coordinates(fiber_a_initial_nodes, fiber_a_output_set, num_frames)
    fiber_b_coordinates = get_all_fiber_coordinates(fiber_b_initial_nodes, fiber_b_output_set, num_frames)

    pressure = get_pressure_array(num_frames)  # pressure in each frame (MPa)

    # navigate to output folder. Abaqus Work Directory and curved_FREEs must be in the same folder
    current_directory = os.getcwd()  # get current work directory
    parent_directory = os.path.split(current_directory)[0]  # go up one folder
    output_directory = os.path.join(parent_directory, 'curved_FREEs', 'output',
                                    myFREE.label)  # navigate to curved_FREEs\output\(myFREE.label)
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)

    # write output data to npz file
    output_label = 'output.npz'
    np.savez(os.path.join(output_directory, output_label),
             FREE_label=myFREE.label,
             endcap_coordinates=endcap_coordinates,
             inner_bladder_coordinates=inner_bladder_coordinates,
             fiber_a_coordinates=fiber_a_coordinates,
             fiber_b_coordinates=fiber_b_coordinates,
             pressure=pressure)


# get pressure array
def get_pressure_array(num_frames):
    pressure = np.zeros(num_frames)

    for f in range(num_frames):
        frame = odb.steps['Step-1'].frames[f]
        pressure[f] = myFREE.final_pressure * frame.frameValue  # pressure applied linearly, so incremental pressure can be found with frameValue

    return pressure


# extract coordinates from all fibers' nodes in initial frame
def get_initial_fiber_coordinates(initial_fiber_nodes):
    initial_fiber_coordinates = np.zeros([myFREE.num_path_points, myFREE.num_fibers, 3])

    for i in range(myFREE.num_path_points):

        for j in range(myFREE.num_fibers):
            n = i + j * myFREE.num_path_points

            # initial_fiber_coordinates[i][j] = point on fiber j at s = s[i]
            initial_fiber_coordinates[i][j] = np.array([initial_fiber_nodes[n].coordinates])

    return initial_fiber_coordinates


# get fiber node locations for all frames
def get_all_fiber_coordinates(initial_fiber_nodes, fiber_output_set, num_frames):
    # initial_fiber_coordinates[i][j] = point on fiber j at s = s[i]
    initial_fiber_coordinates = get_initial_fiber_coordinates(initial_fiber_nodes)

    # fiber_coordinates[i][j] = point on fiber j at s = s[i] in frame f
    fiber_coordinates = np.zeros([num_frames, myFREE.num_path_points, myFREE.num_fibers, 3])

    for f in range(num_frames):
        frame = odb.steps['Step-1'].frames[f]
        displacement = frame.fieldOutputs['U'].getSubset(region=fiber_output_set)

        for i in range(myFREE.num_path_points):

            for j in range(myFREE.num_fibers):
                n = i + myFREE.num_path_points * j

                node_x_disp = displacement.values[n].data[0]
                node_y_disp = displacement.values[n].data[1]
                node_z_disp = displacement.values[n].data[2]

                x = initial_fiber_coordinates[i][j][0] + node_x_disp
                y = initial_fiber_coordinates[i][j][1] + node_y_disp
                z = initial_fiber_coordinates[i][j][2] + node_z_disp

                fiber_coordinates[f][i][j] = np.array([x, y, z])

    return fiber_coordinates


# extract coordinates from the bladder output nodes in initial frame
def get_initial_bladder_coordinates(initial_bladder_nodes):
    initial_bladder_coordinates = np.zeros([myFREE.num_path_points, myFREE.num_output_nodes_circ, 3])

    for i in range(myFREE.num_path_points):

        for j in range(myFREE.num_output_nodes_circ):
            n = i + j * myFREE.num_path_points

            # initial_bladder_coordinates[i][j] = point at circumferential index j at s = s[i]
            initial_bladder_coordinates[i][j] = np.array([initial_bladder_nodes[n].coordinates])

    return initial_bladder_coordinates


# get bladder node locations for all frames
def get_all_bladder_coordinates(initial_bladder_nodes, bladder_output_set, num_frames):
    # initial_bladder_coordinates[i][j] = point at circumferential index j at s = s[i]
    initial_bladder_coordinates = get_initial_bladder_coordinates(initial_bladder_nodes)

    # bladder_coordinates[f][i][j] = point at circumferential index j at s = s[i] in frame f
    bladder_coordinates = np.zeros([num_frames, myFREE.num_path_points, myFREE.num_output_nodes_circ, 3])

    for f in range(num_frames):
        frame = odb.steps['Step-1'].frames[f]
        displacement = frame.fieldOutputs['U'].getSubset(region=bladder_output_set)

        for i in range(myFREE.num_path_points):

            for j in range(myFREE.num_output_nodes_circ):
                n = i + myFREE.num_path_points * j

                node_x_disp = displacement.values[n].data[0]
                node_y_disp = displacement.values[n].data[1]
                node_z_disp = displacement.values[n].data[2]

                x = initial_bladder_coordinates[i][j][0] + node_x_disp
                y = initial_bladder_coordinates[i][j][1] + node_y_disp
                z = initial_bladder_coordinates[i][j][2] + node_z_disp

                bladder_coordinates[f][i][j] = np.array([x, y, z])

    return bladder_coordinates


# get endcap node locations for all frames (offset from center by e1 in the first frame)
def get_all_endcap_coordinates(initial_endcap_nodes, endcap_output_set, num_frames):
    initial_endcap_coordinates = np.zeros([2,3])
    for n in range(2):
        # initial_endcap_coordinates[n] = point n on endcap in first frame
        initial_endcap_coordinates[n] = np.array([initial_endcap_nodes[n].coordinates])

    # endcap_coordinates[f][n] = position in frame f of point n on free endcap at local origin and along e1 in first frame
    endcap_coordinates = np.zeros([num_frames, 2, 3])

    for f in range(num_frames):

        for n in range(np.shape(initial_endcap_coordinates)[0]):
            frame = odb.steps['Step-1'].frames[f]
            displacement = frame.fieldOutputs['U'].getSubset(region=endcap_output_set)

            node_x_disp = displacement.values[n].data[0]
            node_y_disp = displacement.values[n].data[1]
            node_z_disp = displacement.values[n].data[2]

            x = initial_endcap_coordinates[n][0] + node_x_disp
            y = initial_endcap_coordinates[n][1] + node_y_disp
            z = initial_endcap_coordinates[n][2] + node_z_disp

            endcap_coordinates[f][n] = np.array([x, y, z])

    return endcap_coordinates

# run main()
if __name__ == "__main__":
    main()
