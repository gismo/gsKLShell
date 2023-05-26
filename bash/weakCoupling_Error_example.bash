#!/bin/bash
# Declare an array of string with type
declare -a Runname="../../../build/bin/weakCoupling_Error_example"
#declare -a Filenames=("1p_hyperboloid" "3p_hyperboloid" "4p_hyperboloid" "6p_hyperboloid" "4p_hyperboloid_hole_deg4deg4")
declare -a Filenames=("4p_hyperboloid_hole_deg4deg4")
declare -a p=4

mkdir -p Output
  
for file in ${Filenames[@]}; do
#    for (( p=2; p<5; p++)) do
        echo "$Runname" -p $p -s $(($p-1)) -G ../filedata/pde/"$file"_geom.xml -B ../filedata/pde/"$file"_bvp.xml -C 1e3 -r 5
            
        eval "$Runname" -p $p -s $(($p-1)) -G ../filedata/pde/"$file"_geom.xml -B ../filedata/pde/"$file"_bvp.xml -C 1e3 -r 5 > Output/"$file"_p"$p"_weak.log
#    done
done
