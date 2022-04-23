#!/bin/bash
# Declare an array of string with type
declare -a Runname="../../../build/bin/strongCoupling_Error_example"
declare -a Filenames=("1p_hyperboloid" "3p_hyperboloid" "4p_hyperboloid" "6p_hyperboloid" "4p_hyperboloid_hole_deg4deg4")
declare -a Methods=(1 2 3 4)

mkdir -p Output
  
for file in ${Filenames[@]}; do
    for (( p=2; p<4; p++)) do  
        for m in ${Methods[@]}; do
            echo "$Runname" -m $m -p $p -s $(($p-1)) -G ../extensions/gsKLShell/filedata/pde/"$file"_geom.xml -B ../filedata/pde/"$file"_bvp.xml -C 1e3 -r 5
            
            "$Runname" -m $m -p $p -s $(($p-1)) -G ../extensions/gsKLShell/filedata/pde/"$file"_geom.xml -B ../filedata/pde/"$file"_bvp.xml -C 1e3 -r 5 1> Output/"$file"_p"$p"_m"$m".log 2> Output/"$file"_p"$p"_m"$m".log
        done
    done
done
