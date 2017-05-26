# HoleDetection
HoleDetection, a program for detecting holes and computing their boundaries from the outer boundary detected Delaunay triangulation of a planar point set.
Copyright (C) 2017 Subhasree Methirumangalath, Shyam Sundar Kannan, Amal Dev Parakkat, and Ramanathan Muthuganapathy, Advanced Geometric Computing Lab, Department of Engineering Design, IIT Madras, Chennai, India.

For comments or questions, please contact Subhasree Methirumangalath at subhasree.rajiv@gmail.com or Ramanathan Muthuganapathy at mraman@iitm.ac.in.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

--------------------------------------------------------------------------------

DESCRIPTION:

HoleDetection is a program for computing the hole boundaries of an outer boundary reconstructed Delaunay Triangulation using the hole detection algorithm proposed in "Hole Detection of a Planar Point Set: An Empty Disk Approach" by Subhasree Methirumangalath, Shyam Sundar Kannan, Amal Dev Parakkat, and Ramanathan Muthuganapathy to appear at the Shape Modeling International (SMI-2017), with the outer boundary reconstructed using a Delaunay Triangulation method described in "A unified approach towards reconstruction of a planar point set" by Subhasree Methirumangalath, Amal Dev Parakkat and Ramanathan Muthuganapathy.

--------------------------------------------------------------------------------

REQUIREMENTS:

- Ubuntu 16.04 (earlier versions of Ubuntu and other Linux environments will likely also work).
- g++-5                  (required; compiler).
- Boost                  (required; C++ library, http://www.boost.org/).
- CGAL                   (required; C++ library, http://www.cgal.org/).
- OpenGL and GLUT        (required; for visualization).

All necessary packages can be obtained through Ubuntu's package system via apt-get:

sudo apt-get update

sudo apt-get install build-essential

sudo apt-get install cmake

sudo apt-get install libboost-dev

sudo apt-get install libglu1-mesa

sudo apt-get install libxmu-dev libxi-dev

sudo apt-get install freeglut3 freeglut3-dev

sudo apt-get install libcgal-dev

All the packages can be installed using the "install.sh" script provided. The script also compiles the program.

--------------------------------------------------------------------------------

USAGE:

Run format:

./HoleDetection --in filename [--p outer_boundary_parameter]

Description:

  --in - the input file name (.xyz extension is not needed) (required)
  
  --p  - the parameter for outer boundary detection (between 0 to 1) (optional)

Ex.
./HoleDetection --in sampleInput --p 0.7


Input:

The input file specifies n input points in the plane, (x_0, y_0, z_0), ..., (x_{n-1}, y_{n-1}, z_{n-1}) where x_i, y_i, z_i are real numbers and the z value is a constant since the points are planar.

The input point set can be either a dot pattern or a boundary sample.

Input file format:

x_0 y_0 z_0

x_1 y_1 z_1

...

x_{n-1} y_{n-1} z_{n-1}

Input file example :

300 300 1

900 300 1

600 900 1

.

.

.

1000 1000 1


See the input files in "data" folder for samples.

Output:

By default, the output will be shown in an OpenGL window.

--------------------------------------------------------------------------------

NOTES:

The implementation of the HoleDetection algorithm takes the reconstructed outer boundary of the planar point set as the input. The outer boundary reconstruction in based on "A unified approach towards reconstruction of a planar point set" by Subhasree Methirumangalath, Amal Dev Parakkat and Ramanathan Muthuganapathy. The parameter required by the program is for the outer boundary detection and not for the hole detection algorithm.
The parameter is used to define the extent of removal of exterior triangles for the outer boundary reconstruction. (For proper reconstruction of the hole boundaries, it is necessary that the outer boundary is reconstructed correctly.)

