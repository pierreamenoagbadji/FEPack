#!/bin/bash

# Colors
RED='\033[0;31m'
NC='\033[0m' # No Color
YELLOW='\033[1;33m'

# Path to +FEPack
chemin='/home/pierre/Documents/Etudes/these/codes/Chantier\ matlab/FEPack/+FEPack';
cd $chemin

# Launch matlab 
numFloquetPoints=50;

# matlab -nosplash -nodesktop -r "addpath(genpath('../FEPack')); addpath(genpath('+FEPack')); init_script($numFloquetPoints); quit;"

# # Solve each guide problem
# for ((idFB=1; idFB<=$numFloquetPoints; idFB++))
# do

#   # Solve waveguide problem
#   matlab -nosplash -nodesktop -r "addpath(genpath('../FEPack')); addpath(genpath('+FEPack')); solve_TFB_waveguide($idFB); quit;"

# done

