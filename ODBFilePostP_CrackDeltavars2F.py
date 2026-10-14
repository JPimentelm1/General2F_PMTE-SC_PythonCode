# -*- coding: utf-8 -*-
"""
Created on Mon Feb 17 10:08:20 2025

@author: Jose M Pimentel
"""

#-------------------------------------------------------
#||||	 module imports   ||||   	   
#-------------------------------------------------------
import os
from sys import argv
import sys
import math
import numpy as np
from odbAccess import openOdb
import abaqusConstants
import datetime
# Redirect stdout to the command prompt
sys.stdout = sys.__stdout__

start = datetime.datetime.now()
working_directory = os.getcwd()
os.chdir(working_directory)
# type the input file name:
name_INP = 'twofiberInfMtrx0deg'
working_directory_FFM=working_directory+'/FFM_optimizada01_2F0degmu9M2el_PMTESCcriTen2Fibers-t2tcRandp0.0_adaptivefprop'         
working_directory_odbs=working_directory_FFM+'/odbS'
odbDir = working_directory_odbs
inp_file_path = name_INP+'.inp'

with open(inp_file_path, 'r') as inp_file:
    inp_lines = inp_file.readlines()
    inp_file.close()
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

def sigmac_parameter(inp_lines):
    for line in inp_lines:
        if line.startswith('*User Material'):
            # Assuming the first parameter is on the next line
            next_line = inp_lines[inp_lines.index(line) + 1]
            # Split the line by commas (assuming comma-separated values)
            parameters = next_line.split(', ')
            if len(parameters) > 0:
                try:
                    first_param = float(parameters[0])
                    brittleness_num=float(parameters[-1])
                    GIc=float(parameters[1])
                    h=float(parameters[-2])
                    xi=float(parameters[3])
                    lambdhs=float(parameters[2])
                    
                    return first_param, GIc, brittleness_num, h, xi, lambdhs
                except ValueError:
                    print("Error: First parameter is not numeric.")
                    return None
    print("Error: User material section not found.")
    return None
# call of function to obtain the x, y, coordinates of the secondary fibre centre:
xcf2, ycf2 = fiber2XYcoord(inp_lines)
first_param, GIc, mu_const, interfh, xi, lambdhs = sigmac_parameter(inp_lines)

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
        print(file_path)
        odb = openOdb(file_path, readInternalSets=True)
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
        path_SDV13=odb.steps[key_step[0]].frames[-1].fieldOutputs['SDV13'].\
        	getSubset(region=InterElementSet).values
        path_SDV14=odb.steps[key_step[0]].frames[-1].fieldOutputs['SDV14'].\
        	getSubset(region=InterElementSet).values
        # Access the assembly instances
        instance1 = odb.rootAssembly.instances['INTERFACE-1']
        instance2 = odb.rootAssembly.instances['INTERFACE-2']
        
        interfU1 = lastFrame.fieldOutputs['U'].getScalarField(componentLabel='U1',)
        interfU2 = lastFrame.fieldOutputs['U'].getScalarField(componentLabel='U2',)
        # interfc1 = lastFrame.fieldOutputs['COORD'].getScalarField(componentLabel='COOR1',)
        # interfc2 = lastFrame.fieldOutputs['COORD'].getScalarField(componentLabel='COOR2',)

        # Initialize a list to store node labels that satisfy the condition
        node_labelsi1int = []; node_labelsi1ext = []
        node_labelsi2int = []; node_labelsi2ext = []
        # Loop through all nodes in the instance
        for nodei1, nodei2 in zip(instance1.nodes, instance2.nodes):
            x, y, z = nodei1.coordinates  # Get the coordinates of the node1
            x2, y2, z2 = nodei2.coordinates
            if x==0: x=x+1E-6
            elif abs(x2-xcf2)==0: x2=abs(xcf2)+1E-6
            
            if instance1.name=='INTERFACE-1' :
                # Extract the displacement values for the target node
                u1_values = interfU1.getSubset(region=nodei1).values
                u2_values = interfU2.getSubset(region=nodei1).values
                # Check if x^2 + y^2 equals R^2 (with a small tolerance)
                if abs(math.sqrt(x**2 + y**2) - R_o) < toler :  
                    
                    interf1C1 = x
                    interf1C2 = y
                    # umag = math.sqrt(u1_values[0].data**2 + u2_values[0].data**2)
                    interfext_defc1 = interf1C1 + u1_values[0].data
                    interfext_defc2 = abs(interf1C2 + u2_values[0].data)
                    deltacoor1_ext = abs(interfext_defc1)
                    URext = math.sqrt(math.pow(float(deltacoor1_ext),2) + math.pow(interfext_defc2,2))
                    if x>=0 and y>=0:
                        node_labelsi1ext.append([math.atan(interfext_defc2/interfext_defc1)*(180.0/math.pi), interf1C1, interf1C2, URext])  
                        # print([math.atan(y/x)*(180.0/math.pi), interf1C1[0].data, interf1C2[0].data, URext])
                    elif x>=0 and y<0: node_labelsi1ext.append([float(-1)*math.atan((interfext_defc2)/interfext_defc1)*(180.0/math.pi), interf1C1, interf1C2, URext])
                    if x<0 and y>=0: node_labelsi1ext.append([float(180) - math.atan(interfext_defc2/abs(interfext_defc1))*(180.0/math.pi), interf1C1, interf1C2, URext])
                    elif x<0 and y<0: node_labelsi1ext.append([float(-180) + math.atan((interfext_defc2)/abs(interfext_defc1))*(180.0/math.pi), interf1C1, interf1C2, URext])
                # Check if x^2 + y^2 equals R^2 (with a small tolerance)
                if abs(math.sqrt(x**2 + y**2) - R_i) < toler :
                    
                    interf1C1 = x
                    interf1C2 = y
                    # umag = math.sqrt(u1_values[0].data**2 + u2_values[0].data**2)
                    interfint_defc1 = interf1C1 + u1_values[0].data
                    interfint_defc2 = abs(interf1C2 + u2_values[0].data)
                    deltacoor1_int = abs(interfint_defc1)
                    URint = math.sqrt(float(deltacoor1_int**2) + interfint_defc2**2)
                    if x>=0 and y>=0:
                        node_labelsi1int.append([math.atan(interfint_defc2/interfint_defc1)*(180.0/math.pi), interf1C1, interf1C2, URint])  
                        # print([math.atan(y/x)*(180.0/math.pi), interf1C1[0].data, interf1C2[0].data, URext])
                    elif x>=0 and y<0: node_labelsi1int.append([float(-1)*math.atan((interfint_defc2)/interfint_defc1)*(180.0/math.pi), interf1C1, interf1C2, URint])
                    if x<0 and y>=0: node_labelsi1int.append([float(180) - math.atan(interfint_defc2/abs(interfint_defc1))*(180.0/math.pi), interf1C1, interf1C2, URint])
                    elif x<0 and y<0: node_labelsi1int.append([float(-180) + math.atan((interfint_defc2)/abs(interfint_defc1))*(180.0/math.pi), interf1C1, interf1C2, URint])
                    
                    
            if instance2.name=='INTERFACE-2' :
                # Extract the displacement values for the target node
                u1_values = interfU1.getSubset(region=nodei2).values
                u2_values = interfU2.getSubset(region=nodei2).values
                # umag = math.sqrt(u1_values[0].data**2 + u2_values[0].data**2)
                if abs(math.sqrt((x2-xcf2)**2 + (y2-ycf2)**2) - R_o) < toler : 
                    interf2C1 = (x2-xcf2)
                    interf2C2 = (y2-ycf2)
                    interfext_defc1 = interf2C1 + u1_values[0].data
                    interfext_defc2 = abs(interf2C2 + u2_values[0].data)
                    deltacoor1_ext = (interfext_defc1)
                    URext = math.sqrt((deltacoor1_ext**2) + (interfext_defc2**2))
                    if (x2-xcf2)>=0 and (y2-ycf2)>=0:
                        node_labelsi2ext.append([math.atan((interfext_defc2)/(interfext_defc1))*(180.0/math.pi), interf2C1, interf2C2, URext])  
                        # print([math.atan(y/x)*(180.0/math.pi), interf1C1[0].data, interf1C2[0].data, URext])
                    elif (x2-xcf2)>=0 and (y2-ycf2)<0: node_labelsi2ext.append([float(-1)*math.atan(abs(interfext_defc2)/(interfext_defc1))*(180.0/math.pi), interf2C1, interf2C2, URext])
                    if (x2-xcf2)<0 and (y2-ycf2)>=0: node_labelsi2ext.append([float(180) - math.atan((interfext_defc2)/abs(interfext_defc1))*(180.0/math.pi), interf2C1, interf2C2, URext])
                    elif (x2-xcf2)<0 and (y2-ycf2)<0: node_labelsi2ext.append([float(-180) + math.atan(abs(interfext_defc2)/abs(interfext_defc1))*(180.0/math.pi), interf2C1, interf2C2, URext])

                    # node_labelsi2ext.append([180.0 - math.atan(y2/abs(x2))*(180.0/math.pi), interf2C1[0].data, interf2C2[0].data, URext])  # Add the node label to the list
                    # print(nodei2.coordinates)
                if abs(math.sqrt((x2-xcf2)**2 + (y2-ycf2)**2) - R_i) < toler :
                    interf2C1 = (x2-xcf2)
                    interf2C2 = (y2-ycf2)
                    interfint_defc1 = (interf2C1) + u1_values[0].data
                    interfint_defc2 = abs((interf2C2) + u2_values[0].data)
                    deltacoor1_int = abs(interfint_defc1)
                    URint = math.sqrt((deltacoor1_int**2) + (interfint_defc2**2))
                    if (x2-xcf2)>=0 and (y2-ycf2)>=0:
                        node_labelsi2int.append([math.atan((y2-ycf2)/(interfint_defc1))*(180.0/math.pi), interf2C1, interf2C2, URint])  
                        # print([math.atan(y/x)*(180.0/math.pi), interf1C1[0].data, interf1C2[0].data, URext])
                    elif (x2-xcf2)>=0 and (y2-ycf2)<0: node_labelsi2int.append([float(-1)*math.atan(interfint_defc2/(interfint_defc1))*(180.0/math.pi), interf2C1, interf2C2, URint])
                    if (x2-xcf2)<0 and (y2-ycf2)>=0: node_labelsi2int.append([float(180) - math.atan((interfint_defc2)/abs(interfint_defc1))*(180.0/math.pi), interf2C1, interf2C2, URint])
                    elif (x2-xcf2)<0 and (y2-ycf2)<0: node_labelsi2int.append([float(-180) + math.atan((interfint_defc2)/abs(interfint_defc1))*(180.0/math.pi), interf2C1, interf2C2, URint])

                    # node_labelsi2int.append([180.0 - math.atan(y2/abs(x2))*(180.0/math.pi), interf2C1[0].data, interf2C2[0].data, URint]) # Add the node label to the list
                    # print(nodei2.coordinates)
        node_labelsi1int.sort(key=lambda node_labelsi1int:node_labelsi1int[0], reverse=False)        
        node_labelsi1ext.sort(key=lambda node_labelsi1ext:node_labelsi1ext[0], reverse=False)      
        node_labelsi2ext.sort(key=lambda node_labelsi2ext:node_labelsi2ext[0], reverse=False)
        node_labelsi2int.sort(key=lambda node_labelsi2int:node_labelsi2int[0], reverse=False)
        
        thetaUR_arr = []; thetaUR_f2arr = [];
        for nodes_int, nodes_ext, nodes2int, nodes2ext in zip(node_labelsi1int, node_labelsi1ext, node_labelsi2int, node_labelsi2ext):
            DeltaUR = nodes_ext[-1] - nodes_int[-1]
            DeltaURi2 = nodes2ext[-1] - nodes2int[-1]
            thetaUR_arr.append([nodes_ext[0], DeltaUR])
            thetaUR_f2arr.append([nodes2ext[0], DeltaURi2])
        thetaUR_arr.sort(key=lambda thetaUR_arr:thetaUR_arr[0], reverse=False)
        thetaUR_f2arr.sort(key=lambda thetaUR_f2arr:thetaUR_f2arr[0], reverse=False)
        thetaURnp = np.array(thetaUR_arr); thetaURf2np = np.array(thetaUR_f2arr)
        # Save the array into a text file:
            # Specify the format for each column
        fmt = ['%.8e', '%.8e']
        np.savetxt(working_directory_FFM+'\\'+name_files+'_thetaURdeltaf1.txt', thetaURnp, fmt=fmt, delimiter='\t')
        np.savetxt(working_directory_FFM+'\\'+name_files+'_thetaURdeltaf2.txt', thetaURf2np, fmt=fmt, delimiter='\t')
        
        # 2 FIBERS: LOOP ITERATION TO SAVE THE INTERFACE DEBOND ANGLES andS (Sigma_rr, Tau)
        ERRf1_array=[]
        ERRf2_array=[]
        intpoint=len(path_SDV1)
        for value in range(intpoint):
            # Gtot=path_SDV4[value].data
            # GcT=path_SDV5[value].data
            ip=path_SDV2[value].integrationPoint
            NeG=path_SDV2[value].data
            psi_g=math.atan2(path_SDV14[value].data*(math.sqrt(1/xi)), path_SDV13[value].data)
            GcE=float(GIc)*(1+(math.tan(psi_g*(1-lambdhs)))**2)*(mu_const/interfh); psi_gc=(math.pi)/(2*(1-lambdhs));
            x=path_SDV10[value].data
            y=path_SDV11[value].data
            GtE=path_SDV12[value].data
            G2E=path_SDV8[value].data
            if abs(psi_g)>=psi_gc: GcE=GIc*1E8
            GcE=min(GcE, 0.5E-1)
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
                    # Gt2v.append(Gtot)
                    # factorcoord.append([Gtot, GcT])
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
        np.savetxt(working_directory_FFM+'\\'+name_files+'_ERRf1.txt', ERRf1np, fmt=header_fmt, delimiter='\t')
        np.savetxt(working_directory_FFM+'\\'+name_files+'_ERRf2.txt', ERRf2np, fmt=header_fmt, delimiter='\t')

        odb.close()
end = datetime.datetime.now()
timediff=end-start
print('\n')
print('Script has finished in %.2f' %(timediff.total_seconds()/60)+' min')

