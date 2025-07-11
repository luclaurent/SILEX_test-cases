#!/bin/bash

find . -type f -name \*.geo -exec gmsh -3 -format msh22 {} \;
