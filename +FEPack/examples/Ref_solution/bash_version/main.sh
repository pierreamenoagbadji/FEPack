#!/bin/bash

# Colors
RED='\033[0;31m'
YELLOW='\033[1;33m'
GREEN='\033[0;92m'
NC='\033[0m' # No Color

# Path to +FEPack
chemin='/home/pierre/Documents/Etudes/these/codes/FEPack/+FEPack';

# Folder for outputs 
cheminDonnees='/home/pierre/Documents/Etudes/these/codes/FEPack/+FEPack/outputs';
cd $chemin

# Delete outputs/* if requested
while getopts ":r" opt; do
  case $opt in
    r)
      echo -e "${GREEN}-r was triggered: deleting $cheminDonnees/*${NC}" >&2
      rm $cheminDonnees/*
      ;;
    \?)
      echo -e "${RED}Invalid option: -$OPTARG. Exit code${NC}" >&2
      exit 1
      ;;
  esac
done

# Initialization
numCellsZ=40; # Total number of cells
sizeClust=10;  # Number of points per cluster
(( numClusters=numCellsZ/sizeClust )); # Number of clusters
numCellsZsave=20;
sizeClustSave=5;  # Number of points per cluster
(( numSaveClusters=numCellsZsave/sizeClustSave )); # Number of clusters

matlab -nosplash -nodesktop -nojvm -r "addpath(genpath('../../FEPack')); addpath(genpath('../+FEPack')); numCellsZ=$numCellsZ; cheminDonnees='$cheminDonnees'; init_script; quit;"

# ///////////////////////////////////////////////////////////#
# Positive side                                              #
# ///////////////////////////////////////////////////////////#
matlab -nosplash -nodesktop -nojvm -r "addpath(genpath('../../FEPack')); addpath(genpath('../+FEPack')); half_guide_script_avant('pos', '$cheminDonnees'); quit;"

for (( idCluster=0; idCluster<$numClusters; idCluster++ ))
do
  # Show progress
  echo -e "${RED}/////////////////  Cluster $idCluster  /////////////////${NC}"

  # Solve the cell problems associated to the current cluster
  matlab -nosplash -nodesktop -nojvm -r "addpath(genpath('../../FEPack')); addpath(genpath('../+FEPack')); local_cell_script($idCluster*$sizeClust+1, ($idCluster+1)*$sizeClust, '$cheminDonnees'); quit;" &
done
wait # for all the local cell problems to be solved

matlab -nosplash -nodesktop -nojvm -r "addpath(genpath('../../FEPack')); addpath(genpath('../+FEPack')); cheminDonnees='$cheminDonnees'; half_guide_script_apres; quit;"

rm $cheminDonnees/inputs_half_guide.mat

# ///////////////////////////////////////////////////////////#
# Negative side                                              #
# ///////////////////////////////////////////////////////////#
matlab -nosplash -nodesktop -nojvm -r "addpath(genpath('../../FEPack')); addpath(genpath('../+FEPack')); half_guide_script_avant('neg', '$cheminDonnees'); quit;"

for (( idCluster=0; idCluster<$numClusters; idCluster++ ))
do
  # Show progress
  echo -e "${RED}/////////////////  Cluster $idCluster  /////////////////${NC}"

  # Solve the cell problems associated to the current cluster
  matlab -nosplash -nodesktop -nojvm -r "addpath(genpath('../../FEPack')); addpath(genpath('../+FEPack')); local_cell_script($idCluster*$sizeClust+1, ($idCluster+1)*$sizeClust, '$cheminDonnees'); quit;" &
done
wait # for all the local cell problems to be solved

matlab -nosplash -nodesktop -nojvm -r "addpath(genpath('../../FEPack')); addpath(genpath('../+FEPack')); cheminDonnees='$cheminDonnees'; half_guide_script_apres; quit;"

rm $cheminDonnees/inputs_half_guide.mat

# ///////////////////////////////////////////////////////////#
# Transmission condition and guide solution                  #
# ///////////////////////////////////////////////////////////#
matlab -nosplash -nodesktop -nojvm -r "addpath(genpath('../../FEPack')); addpath(genpath('../+FEPack')); cheminDonnees='$cheminDonnees'; transmission_condition_script; quit;"

for (( idCluster=0; idCluster<$numSaveClusters; idCluster++ ))
do
  # Show progress
  echo -e "${RED}/////////////////  Cluster $idCluster  /////////////////${NC}"

  # Solve the cell problems associated to the current cluster
  matlab -nosplash -nodesktop -nojvm -r "addpath(genpath('../../FEPack')); addpath(genpath('../+FEPack')); compute_guide_solution_script($idCluster*$sizeClustSave+1+($numCellsZ - $numCellsZsave)/2, ($idCluster+1)*$sizeClustSave+($numCellsZ - $numCellsZsave)/2, '$cheminDonnees'); quit;" &
done
wait # for all the local cell problems to be solved

# ///////////////////////////////////////////////////////////#
# End script                                                 #
# ///////////////////////////////////////////////////////////#
matlab -nosplash -nodesktop -r "addpath(genpath('../../FEPack')); addpath(genpath('../+FEPack')); numCellsZsave = $numCellsZsave; cheminDonnees='$cheminDonnees'; end_script; quit;"

# Delete some now useless outputs
rm $cheminDonnees/guide_sol_*
rm $cheminDonnees/local_cell_sol_*
rm $cheminDonnees/half_guide_solution_*