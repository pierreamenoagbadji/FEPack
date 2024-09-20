#!/bin/bash

# Arg 1: image prefix (with path)
# Arg 2: number of tiles along X
# Arg 3: number of tiles along Y
# Arg 4: output image (with path)

# Extract directory name from path
outDIR="$(dirname "$4")";
outFILE="$(basename "$4")";
mkdir ${outDIR}/joinImages_temp

# for ((idY=1; idY<=$3; idY++))
# do
#   echo $idY

#   for ((idX=1; idX<=$2; idX++))
#   do
#     # echo $1_X_${idX}_Y_${idY}.png
#     mv $1_X_${idX}_Y_${idY}.png $1_X_$(printf '%03d\n' ${idX})_Y_${idY}.png
#   done
# done


for ((idY=1; idY<=$3; idY++))
do
  echo $idY
  idYout="$(($3-$idY+1))"
  convert $1_X_*_Y_$idY.png -resize 2% +append ${outDIR}/joinImages_temp/out_$(printf '%03d\n' $idYout).png
done

convert ${outDIR}/joinImages_temp/out_*.png -append $4
