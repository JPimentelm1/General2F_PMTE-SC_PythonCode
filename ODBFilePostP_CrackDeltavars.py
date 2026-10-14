# -*- coding: utf-8 -*-
"""
Created on Mon Feb 17 10:08:20 2025

@author: Jose M Pimentel
"""

#-------------------------------------------------------
#||||	 module imports   ||||   	   
#-------------------------------------------------------
import os
from odbAccess import*
from sys import argv
import sys
import math
import numpy as np
from odbAccess import openOdb
from abaqusConstants import READ_WRITE
# Redirect stdout to the command prompt
sys.stdout = sys.__stdout__

working_directory = os.getcwd()
os.chdir(working_directory)
# type the input file name:
name_INP = 'twofiberInfMtrx0deg'
inp_file_path = name_INP+'.inp'
with open(inp_file_path, 'r') as inp_file:
    inp_lines = inp_file.readlines()
    
def fiber2XYcoord(inp_lines):
    for line in inp_lines:
        if line.startswith('*Instance, name=Interface-2, '):
            # Assuming the first parameter is on the next line
            next_line = inp_lines[inp_lines.index(line) + 1]
            # Split the line by commas (assuming comma-separated values)
            parameters = next_line.split(', ')
            if len(parameters) > 0:
                try:
                    xcoordf2 = float(parameters[0])
                    ycoordf2 = float(parameters[1])
                    return xcoordf2, ycoordf2
                except ValueError:
                    print("Error: First parameter is not numeric.")
                    return None
    print("Error: Second instance not found.")
    return None
# call of function to obtain the x, y, coordinates of the secondary fibre centre:
xcf2, ycf2 = fiber2XYcoord(inp_lines)

working_directory_FFM=working_directory+'/FFM_optimizada01_2F0degmu9.0_PMTESCcriTen2Fibers-t2tcRandp0.0_adaptivefprop'          
working_directory_odbs=working_directory_FFM+'/odbS'
odbDir = working_directory_odbs

# Define the radius R
R_o = 7.51
R_i = 7.50
toler = 1E-6
# Loop through the files in the directory
for filename in os.listdir(odbDir):
    # Check if the file has a .lck extension 
    if filename.endswith('.odb'):
        # Construct the full file path 
        file_path = os.path.join(odbDir, filename)
        name_files = filename
        odb = openOdb(file_path, mode=READ_WRITE)
        myAssembly = odb.rootAssembly
        key_step = odb.steps.keys()
        lastFrame = odb.steps[key_step[0]].frames[-1]
        InterElementSet = odb.rootAssembly.elementSets['NINTERFACE']
             
        # Access the assembly instances
        instance1 = odb.rootAssembly.instances['INTERFACE1-1']
        instance2 = odb.rootAssembly.instances['INTERFACE2-1']
        
        interfU1 = lastFrame.fieldOutputs['U'].getScalarField(componentLabel='U1',)
        interfU2 = lastFrame.fieldOutputs['U'].getScalarField(componentLabel='U2',)
        interfc1 = lastFrame.fieldOutputs['COORD'].getScalarField(componentLabel='COOR1',)
        interfc2 = lastFrame.fieldOutputs['COORD'].getScalarField(componentLabel='COOR2',)

        # Initialize a list to store node labels that satisfy the condition
        node_labelsi1 = []
        node_labelsi2 = []
        # Loop through all nodes in the instance
        for nodei1, nodei2 in zip(instance1.nodes, instance2.nodes):
            x, y, z = nodei1.coordinates  # Get the coordinates of the node1
            x2, y2, z2 = nodei2.coordinates
            
            # Check if x^2 + y^2 equals R^2 (with a small tolerance)
            if instance1.name=='INTERFACE1-1' and \
                (abs(math.sqrt(x**2 + y**2) - R_o) < toler or abs(math.sqrt(x**2 + y**2) - R_i) < toler):  
                node_labelsi1.append(nodei1.label)  # Add the node label to the list
                print(nodei1.label)
            if instance2.name=='INTERFACE2-1' and \
                (abs(math.sqrt((x2-xcf2)**2 + (y2-ycf2)**2) - R_o) < toler or abs(math.sqrt((x2-xcf2)**2 + (y2-ycf2)**2) - R_i) < toler):  
                node_labelsi2.append(nodei2.label)  # Add the node label to the list
                print(nodei2.label)
        # Create the node sets and assemble array of coordinates and CODisp:
        # if node_labelsi1 :
        #     node_set_name = 'NODES_AT_RADIUS_R'
        #     assembly.NodeSetFromNodeLabels(name=node_set_name, nodeLabels=node_labels, instance=instance)
        #     print(f"Node set '{node_set_name}' created with {len(node_labels)} nodes.")
        # else:
        #     print(f"No nodes found at radius R = {R}.")

        
        odb.close()


