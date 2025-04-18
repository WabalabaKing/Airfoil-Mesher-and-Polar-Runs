# -*- coding: utf-8 -*-
"""
Created on Tue Feb 18 22:18:33 2025

@author: wz10
"""

import os
import subprocess
import csv
import numpy as np
def modify_su2_cfg(cfg_file, aoa, restart):
    """Modify the SU2 configuration file for a given angle of attack and restart setting."""
    with open(cfg_file, 'r') as file:
        lines = file.readlines()
    
    with open(cfg_file, 'w') as file:
        for line in lines:
            if line.startswith('AOA='):
                file.write(f'AOA=  {aoa}\n')
            elif line.startswith('RESTART_SOL'):
                file.write(f'RESTART_SOL= {"YES" if restart else "NO"}\n')
            else:
                file.write(line)

def run_su2(cfg_file, angles, num_procs=8):
    """Run SU2 for a range of angles of attack and collect CL, CD, CMz."""
    results = []

    for i, aoa in enumerate(angles):
        restart = i > 0  # Restart from the second case onwards
        modify_su2_cfg(cfg_file, aoa, restart)
        
        # Run SU2 with MPI
        command = ["mpiexec", "-n", str(num_procs), "SU2_CFD", config_file]
        process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
        # Stream the output to terminal in real-time
        for line in iter(process.stdout.readline, ""):
            print(line, end="")  # Print line by line without extra newlines

        process.stdout.close()
        process.wait()  # Ensure SU2_CFD finishes before proceeding
        # Extract aerodynamic coefficients
        cl, cd, cmz = extract_forces("forces_breakdown.dat")
        results.append((aoa, cl, cd, cmz))
        print(f"AoA: {aoa}, CL: {cl}, CD: {cd}, CMz: {cmz}")

    return results
def extract_forces(filename):
    """Extract CL, CD, and CMz from forces_breakdown.dat using string splitting."""
    cl, cd, cmz = None, None, None
    
    if not os.path.exists(filename):
        print(f"Warning: {filename} not found!")
        return cl, cd, cmz

    with open(filename, 'r') as file:
        for line in file:
            if line.startswith("Total CL:"):
                parts = line.split("|")
                cl = float(parts[0].split(":")[-1].strip())
            elif line.startswith("Total CD"):
                parts = line.split("|")
                cd = float(parts[0].split(":")[-1].strip())
            elif line.startswith("Total CMz:"):
                parts = line.split("|")
                cmz = float(parts[0].split(":")[-1].strip())

    return cl, cd, cmz
# Example usage
angles_of_attack = np.linspace(-4,4,9)  # Define angles of attack to test
config_file = "turb_NACA0012.cfg"

results = run_su2(config_file, angles_of_attack)
with open("NACA0012_results.csv", "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["AoA", "CL", "CD", "CMz"])
    writer.writerows(results)

print("Results saved to NACA0012_results.csv")