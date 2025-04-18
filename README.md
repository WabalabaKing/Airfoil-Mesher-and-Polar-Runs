1. break airfoil into upper and lower surface and store coordinates in csv file

2. Run airfoilMesher.py

3. open Gmsh, open the output .geo file from previous python script, export mesh into .su2

4. Edit turb_NACA0012.cfg for your Mach and Reynolds number, then run RunPolars.py. DO NOT RUN FOR EXCESSIVE AOA!

5. CL,CD CM data will be stored in a .csv file
