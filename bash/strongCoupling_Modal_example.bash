#!/bin/bash
# Declare an array of string with type
declare -a Runname="./bin/strongCoupling_Modal_example"
declare -a Meshes=("2x2" "3x3" "4x4")
declare -a BVP="../extensions/gsKLShell/filedata/pde/car_bvp.xml"
declare -a Methods=(1 4)
declare -a p=2
declare -a s=1

for (( r=0; r<4; r++)) do
    for mesh in ${Meshes[@]}; do
        for method in ${Methods[@]}; do
            echo "$Runname" -G ../filedata/surfaces/neon_side_fit_"$mesh".xml -B $BVP -m $method -r $r -p $p -s $s --plot -N 16 -S 1e-5 --write -o ModalResults/Car_"$mesh"_r"$r"_p"$p"_s"$s"_m"$method"
            eval "$Runname" -G ../filedata/surfaces/neon_side_fit_"$mesh".xml -B $BVP -m $method -r $r -p $p -s $s --plot -N 16 -S 1e-5 --write -o ModalResults/Car_"$mesh"_r"$r"_p"$p"_s"$s"_m"$method"
        done
    done
done
