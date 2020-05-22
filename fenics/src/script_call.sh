#!/bin/bash

#gmsh -3 -format msh40 -o sphere_msh40.msh sphere.geo
gmsh -3 -format vtk -o sphere.vtk sphere.geo
