# -*- coding: utf-8 -*-
"""
Created on Tue Feb 18 18:31:50 2025

@author: wz10
"""

import numpy as np
from scipy.interpolate import interp1d

def read_airfoil_coordinates(airfoil_file, num_points):
    data = np.loadtxt(airfoil_file, delimiter=',')
    x, y = data[:, 0], data[:, 1]
    
    
    # Interpolation functions
    refined_x = np.flip((np.cos(np.linspace(-np.pi,-3*np.pi/2,num_points)))+1)
    f = interp1d(x, y, kind='quadratic')
    
    # Generate refined points
    refined_y = f(refined_x)
    
    points = np.array([refined_x, refined_y])
    return points
def compute_first_layer_height(reynolds, y_plus,mach,reference_length=1.0, nu=1.5e-5):
    # Estimate skin friction coefficient
    cf = 0.026 /(reynolds ** (1/7))
    
    # Estimate friction velocity
    u_tau = np.sqrt(cf*1.225*(mach*340)**2/2/1.225)
    
    # Compute first layer height
    h1 = (y_plus * nu) /( u_tau*1.225)
    return h1

def write_gmsh_geo(airfoil_fileU,airfoil_fileD, reynolds, y_plus, growth_rate, bl_layers, farfield_size, boundary_conditions,points,mach):
    geo_filename = airfoil_fileU.replace('U.csv', '.geo')
    first_layer_height = compute_first_layer_height(reynolds, y_plus,mach)
    with open(geo_filename, 'w') as geo:
        geo.write(f'// Gmsh script for 2D airfoil mesh\n')
        geo.write(f'// Airfoil: {airfoil_fileU}\n')
        #geo.write('SetFactory("OpenCASCADE");\n')
        # Read airfoil coordinates
        num_farfield_points = 100
        farfield_points = []
        for j in range(num_farfield_points):
            angle = 2 * np.pi * j / num_farfield_points
            x_far = farfield_size * np.cos(angle)
            y_far = farfield_size * np.sin(angle)
            pid = j + 1
            geo.write(f'Point({pid}) = {{{x_far}, {y_far}, 0, 1.0}};\n')
            farfield_points.append(pid)
            counter=0
        for i in range(num_farfield_points-1):
            geo.write(f'Line({counter+1})={{{i+1},{i+2}}};\n')
            counter=counter+1
        geo.write(f'Line({counter+1})={{{num_farfield_points},1}};\n')
        geo.write(f'Curve Loop(1) = {{{", ".join(str(p) for p in farfield_points)}}};\n')
        pointsU =read_airfoil_coordinates(airfoil_fileU,points)
        pointsD =  np.flip(read_airfoil_coordinates(airfoil_fileD,points),axis=1)
        #solve LE location for curvature continuity
        pUp = pointsU[:,-4:-1]
        pDp = pointsD[:,1:4]
        pp = np.hstack((pUp,pDp))
        cur = interp1d(pp[1,:], pp[0,:], kind='cubic')
        xp = cur(0)
        points = np.hstack((pointsU[:,:-1],[[xp],[0]],pointsD[:,1:],[[pointsU[0,0]],[(pointsU[1,0]+pointsD[1,-1])/2]]))
        for i in range(len(points[0])):
            X = points[0,i]
            Y = points[1,i]
            geo.write(f'Point({num_farfield_points+i+1}) = {{{X}, {Y},0, 1.0}};\n')
        # Define airfoil curve
        counter=num_farfield_points
        for i in range(len(points[0])-1):
            geo.write(f'Line({counter+1})={{{num_farfield_points+i+1},{num_farfield_points+i+2}}};\n')
            counter=counter+1
        geo.write(f'Line({counter+1})={{{num_farfield_points+len(points[0])},{num_farfield_points+1}}};\n')
        geo.write(f'Curve Loop(2) = {{{", ".join(str(p+1+num_farfield_points) for p in range(len(points[0])))}}};\n')
        geo.write('Plane Surface(1) = {1, 2};\n')
        # Boundary layer settings
        geo.write(f"hwall_n = {first_layer_height};\n")
        geo.write(f"thickness = {bl_layers};\n")
        geo.write(f"ratio = {growth_rate};\n")
        geo.write("use_quads = 1;\n\n")
        geo.write("// Apply boundary layer mesh near airfoil\n")
        geo.write("Field[1] = BoundaryLayer;\n")
        geo.write(f'Field[1].CurvesList = {{{", ".join(str(p+1+num_farfield_points) for p in range(len(points[0])))}}}; \n')
        geo.write("Field[1].hwall_n = hwall_n;\n")
        geo.write("Field[1].thickness = thickness;\n")
        geo.write("Field[1].ratio = ratio;\n")
        geo.write("Field[1].Quads = use_quads;\n")
        
        # geo.write("Field[2] = Wake;\n")
        # geo.write(f'Field[2].CurvesList = {{{len(points[0])},{len(points[0])+1}}}; \n')
        # geo.write("Field[2].hwall_n = 0.01;\n")
        # geo.write("Field[2].thickness = 3;\n")
        # geo.write("Field[2].ratio = 1;\n")
        # geo.write("Field[2].Quads = use_quads;\n")
        geo.write("BoundaryLayer Field = 1;\n\n")
        #geo.write("BoundaryLayer Field = 2;\n\n")
        #
        
        geo.write(f'Physical Curve("Airfoil") = {{{", ".join(str(p+1+num_farfield_points) for p in range(len(points[0])-1))}}};\n')
        geo.write(f'Physical Curve("Farfield") = {{{", ".join(str(p) for p in farfield_points)}}};\n')
        geo.write('Mesh.MeshSizeFactor = 100;\n')
        geo.write("Mesh 1;\n")
        geo.write("Mesh 2;\n")
        
        
    print(f'Gmsh .geo file written to {geo_filename}')
    return geo_filename

# Test case for NACA0012
airfoil_fileU = 'NACA64204U.csv'
airfoil_fileD = 'NACA64204D.csv'
reynolds = 2.8e6
y_plus = 1
growth_rate = 1.2
num_bl_layers = 0.08
farfield_size = 230
boundary_conditions = {"airfoil": 1, "farfield": 2}
mach=0.78
write_gmsh_geo(airfoil_fileU,airfoil_fileD, reynolds, y_plus, growth_rate, num_bl_layers, farfield_size, boundary_conditions,80,mach)
