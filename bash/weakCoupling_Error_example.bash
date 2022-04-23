#!/bin/bash
# Declare an array of string with type
declare -a Runname="../../../build/bin/weakCoupling_Error_example"
declare -a Filenames=("1p_hyperboloid" "3p_hyperboloid" "4p_hyperboloid" "6p_hyperboloid" )

mkdir -p Output
  
for file in ${Filenames[@]}; do
    for (( p=2; p<4; p++)) do  
        echo "$Runname" -p $p -s $(($p-1)) -G ../extensions/gsKLShell/filedata/pde/"$file"_geom.xml -B ../extensions/gsKLShell/filedata/pde/"$file"_bvp.xml -C 1e3 -r 3
            
        "$Runname" -p $p -s $(($p-1)) -G ../extensions/gsKLShell/filedata/pde/"$file"_geom.xml -B ../extensions/gsKLShell/filedata/pde/"$file"_bvp.xml -C 1e3 -r 3 1> Output/"$file"_p"$p"_weak.log 2> Output/"$file"_p"$p"_weak.log
    done
done
