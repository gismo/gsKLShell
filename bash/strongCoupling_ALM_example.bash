#!/bin/bash
# Declare an array of string with type
declare -a Runname="../../../build/bin/strongCoupling_ALM_example"
declare -a Outname="../../../build/ArcLengthResults"
declare -a Filenames=("17p_Lshape_4holes" )
declare -a Methods=(1 2 3 4)
  
mkdir -p Output
for file in ${Filenames[@]}; do
    for (( p=2; p<4; p++)) do  
        for m in ${Methods[@]}; do
            mkdir -p "$Outname"/"$file"_p"$p"_m"$m"            
        
            echo "$Runname" -m $m -p $p -s $(($p-1)) -r 2 -L 5e-3 -l 2e-1 -N 20 --plot --write -F 1e1 -C 1e3 -G ../filedata/pde/"$file"_geom.xml -B ../filedata/pde/"$file"_bvp.xml -o "$Outname"/"$file"_p"$p"_m"$m"
            
            "$Runname" -m $m -p $p -s $(($p-1)) -r 2 -L 5e-3 -l 2e-1 -N 20 --plot --write -F 1e1 -C 1e3 -G ../filedata/pde/"$file"_geom.xml -B ../filedata/pde/"$file"_bvp.xml -o "$Outname"/"$file"_p"$p"_m"$m"  1> Output/"$file"_p"$p"_m"$m".log 2> Output/"$file"_p"$p"_m"$m".log
        done
    done
done
