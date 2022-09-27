To build the executable on Linux (I use WSL2 with Ubuntu 20.04): 

1. Install Eigen (http://eigen.tuxfamily.org/index.php?title=Main_Page#Download): 
download the source files and add a symlink to the directory, e.g.

``ln -s /mnt/c/users/lenny/documents/simulation/eigen/Eigen /usr/local/include/Eigen``.

2. Install MFEM (https://mfem.org/building/). I moved the /mfem/ directory to the project
folder; the g++ command can be modified accordingly (in particular, update the `-Lmfem` flag
to point to the correct directory where MFEM is locatied, e.g. `-L..` if you placed it
in the parent folder).

Then, under the project directory, run the command 

``g++ *.cpp -o name_of_your_executable -O3 -std=c++17 -I. -Lmfem -lmfem -lrt``

I used Gmsh (https://gmsh.info/) to generate 2D mesh in 3D space. However, exporting
into .vtk format is tricky:
  1. You may have to remove boundary lines on the mesh in the Gmsh GUI.
  2. In order to properly read the mesh data with the MFEM library, you need to 
  manually remove cells on the vertices (cell type == 1 in the .vtk file)
  and update the corresponding counts in the .vtk file. 
