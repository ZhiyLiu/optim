

#To use this file execute: sh execImage2Distance.sh <Filename binary image>
#example sh execImage2Distance.sh Right_Hippocampus/Right_Hippocampus.mhd
echo ${1}
#/home/north/src/MATLAB/Image2DistanceTransform/Image2SignedDistanceMap ${1}
./Image2SignedDistanceMap ${1}

fn=${1%.*}

cd MATLAB

matlab_exec=/usr/local/MATLAB/R2010b/bin/matlab
X="antiAliasWrapper('${fn}-ddm.mhd','${fn}',[0.5 0.5]);
   quit;
  "
echo ${X} > temp.m
cat temp.m
${matlab_exec} -nodisplay -r temp.m

rm temp.m 

