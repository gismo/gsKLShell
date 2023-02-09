#!/bin/bash
# Declare an array of string with type
declare -a Runname="../../../build_deb/bin/strongCoupling_Error_example"
declare -a Filenames=("3p_hyperboloid" "4p_hyperboloid2" "6p_hyperboloid" "6p_hyperboloid2")
# 1p_hyperboloid
declare -a Methods=(3)

for file in ${Filenames[@]}; do
    for (( p=4; p<6; p++)) do
        for m in ${Methods[@]}; do
            declare -a s=$((p-2))
            declare -a R0=2
            echo "$Runname" -m $m -p $p -s $s -G ../filedata/pde/"$file"_geom.xml -B ../filedata/pde/"$file"_bvp.xml -C 1e13 -r 5 -R $R0
            eval "$Runname" -m $m -p $p -s $s -G ../filedata/pde/"$file"_geom.xml -B ../filedata/pde/"$file"_bvp.xml -C 1e13 -r 5 -R $R0 -w Output/"$file"_p"$p"_s"$s"_m"$m".csv >  Output/"$file"_p"$p"_s"$s"_m"$m".log
        done
    done
done
