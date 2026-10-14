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
# os.chdir(working_directory_FFM)

# Loop through the files in the directory
for filename in os.listdir(odbDir):
    # Check if the file has a .lck extension 
    if filename.endswith('.odb'):
        # Construct the full file path 
        file_path = os.path.join(odbDir, filename)
        name_files = filename
        odb = openOdb(file_path)
        myAssembly = odb.rootAssembly
        key_step = odb.steps.keys()
        lastFrame = odb.steps[key_step[0]].frames[-1]
        InterElementSet = odb.rootAssembly.elementSets['NINTERFACE']
        path_SDV1=odb.steps[key_step[0]].frames[-1].fieldOutputs['SDV1'].\
         	getSubset(region=InterElementSet).values
        path_SDV2=odb.steps[key_step[0]].frames[-1].fieldOutputs['SDV2'].\
         	getSubset(region=InterElementSet).values
        path_SDV4=odb.steps[key_step[0]].frames[-1].fieldOutputs['SDV4'].\
         	getSubset(region=InterElementSet).values
        path_SDV5=odb.steps[key_step[0]].frames[-1].fieldOutputs['SDV5'].\
         	getSubset(region=InterElementSet).values
        path_SDV6=odb.steps[key_step[0]].frames[-1].fieldOutputs['SDV6'].\
         	getSubset(region=InterElementSet).values   
        path_SDV8=odb.steps[key_step[0]].frames[-1].fieldOutputs['SDV8'].\
         	getSubset(region=InterElementSet).values  
        path_SDV10=odb.steps[key_step[0]].frames[-1].fieldOutputs['SDV10'].\
         	getSubset(region=InterElementSet).values 
        path_SDV11=odb.steps[key_step[0]].frames[-1].fieldOutputs['SDV11'].\
         	getSubset(region=InterElementSet).values
        path_SDV12=odb.steps[key_step[0]].frames[-1].fieldOutputs['SDV12'].\
         	getSubset(region=InterElementSet).values
             
        # Access the assembly instances
        instance1 = odb.rootAssembly.instances['INTERFACE-1']
        instance2 = odb.rootAssembly.instances['INTERFACE-2']
        
        # 2 FIBERS: LOOP ITERATION TO SAVE THE INTERFACE DEBOND ANGLES andS (Sigma_rr, Tau)
        ERRf1_array=[]
        ERRf2_array=[]
        intpoint = len(path_SDV1)
        for value in range(intpoint):
            # Gtot=path_SDV4[value].data
            # GcT=path_SDV5[value].data
            x=path_SDV10[value].data
            y=path_SDV11[value].data
            GtE=path_SDV12[value].data
            G2E=path_SDV8[value].data
            GcE=path_SDV6[value].data
            GcE=min(GcE, 1E-1)
            # sigmarr=path_SDV13[value].data
            # sigmart=path_SDV14[value].data
            
            NeL2=str(path_SDV2[value].elementLabel)+'_'+str(path_SDV2[value].instance.name)
            if NeL2.endswith('_INTERFACE-1'):
                if y >= 0 and x >= 0:
                    thetad1=(math.atan(y/x)*(180/(math.pi)))
                    ERRf1_array.append([thetad1, G2E, GtE, GcE])
                    # Gt1v.append(Gtot)
                    # factorcoord.append([Gtot, GcT])
                elif y >= 0  and x < 0:
                    thetad1=float(180) - (math.atan(y/abs(x))*(180/(math.pi)))
                    ERRf1_array.append([thetad1, G2E, GtE, GcE])
                    # Gt2v.append(Gtot)
                    # factorcoord.append([Gtot, GcT])
                if y < 0 and x >= 0:
                    thetad1=float(-1)*math.atan(abs(y)/x)*(180/(math.pi))
                    ERRf1_array.append([thetad1, G2E, GtE, GcE])
                elif y < 0  and x < 0:
                    thetad1=float(-180) + (math.atan(abs(y)/abs(x))*(180/(math.pi)))
                    ERRf1_array.append([thetad1, G2E, GtE, GcE])
            if NeL2.endswith('_INTERFACE-2') and (x-xcf2)>0:
                if (y-ycf2) >= 0 :
                    thetad2=math.atan((y-ycf2)/(x-xcf2))*(180/(math.pi))
                    ERRf2_array.append([thetad2, G2E, GtE, GcE])
                    # Gt1v.append(Gtot)
                    # factorcoord.append([Gtot, GcT])
                elif (y-ycf2) < 0:
                    thetad2=float(-1)*math.atan(abs(y-ycf2)/(x-xcf2))*(180/(math.pi))
                    ERRf2_array.append([thetad2, G2E, GtE, GcE])
            if NeL2.endswith('_INTERFACE-2') and (x-xcf2)<0:
                if (y-ycf2) >= 0 :
                    thetad2=float(180) - math.atan((y-ycf2)/abs(x-xcf2))*(180/(math.pi))
                    ERRf2_array.append([thetad2, G2E, GtE, GcE])

                elif (y-ycf2) < 0:
                    thetad2=float(-180) + math.atan(abs(y-ycf2)/abs(x-xcf2))*(180/(math.pi))
                    ERRf2_array.append([thetad2, G2E, GtE, GcE])

        # Sort the ERR lists in ascending order of thetad:
        ERRf1_array.sort(key=lambda ERRf1_array:ERRf1_array[:][0], reverse=False)
        ERRf2_array.sort(key=lambda ERRf2_array:ERRf2_array[:][0], reverse=False)
        ERRf1np = np.array(ERRf1_array); ERRf2np = np.array(ERRf2_array);
        # Save the arrays as text files:
            # Specify the format for each column
        header_fmt = ['%.8e', '%.8e', '%.8e', '%.8e']
        np.savetxt(working_directory_FFM+'//'+name_files+'_ERRf1.txt', ERRf1np, fmt=header_fmt, delimiter='\t')
        np.savetxt(working_directory_FFM+'//'+name_files+'_ERRf2.txt', ERRf2np, fmt=header_fmt, delimiter='\t')
        odb.close()


