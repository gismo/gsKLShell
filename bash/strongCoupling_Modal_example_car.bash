#!/bin/bash
# Declare an array of string with type
declare -a Runname="./bin/strongCoupling_Modal_example"
declare -a Meshes=("0x0" "1x1" "2x2" "3x3" "4x4")
declare -a BVP="../optional/gsKLShell/filedata/pde/car_bvp.xml"
declare -a Methods=(1 4)
declare -a p=2
declare -a s=1
declare -a S=-1e-6
declare -a out="/mnt/sda2/hverhelst/"

for (( r=4; r<6; r++)) do
    for mesh in ${Meshes[@]}; do
        for method in ${Methods[@]}; do
            echo "$Runname" -G ../optional/gsUnstructuredSplines/filedata/surfaces/neon/neon_side_fit_"$mesh"_p2.xml -B $BVP -m $method -r $r -p $p -s $s --plot -N 16 -S $S --write -o "$out"ModalResults/Car_"$mesh"_r"$r"_p"$p"_s"$s"_m"$method"
            eval "$Runname" -G ../optional/gsUnstructuredSplines/filedata/surfaces/neon/neon_side_fit_"$mesh"_p2.xml -B $BVP -m $method -r $r -p $p -s $s --plot -N 16 -S $S --write -o "$out"ModalResults/Car_"$mesh"_r"$r"_p"$p"_s"$s"_m"$method" 1> "$out"ModalResults/Car_"$mesh"_r"$r"_p"$p"_s"$s"_m"$method".log 2> "$out"ModalResults/Car_"$mesh"_r"$r"_p"$p"_s"$s"_m"$method".err
        done
    done
done

#declare -a Methods=(2)
#declare -a p=3
#declare -a s=2
#for (( r=1; r<2; r++)) do
#    for mesh in ${Meshes[@]}; do
#        for method in ${Methods[@]}; do
#            echo "$Runname" -G ../extensions/gsUnstructuredSplines/filedata/surfaces/neon/neon_side_split_fit_"$mesh"_smooth.xml -B $BVP -m $method -r $r -p $p -s $s --plot -N 16 -S $S --write -o "$out"ModalResults/Car_"$mesh"_r"$r"_p"$p"_s"$s"_m"$method"
#            eval "$Runname" -G ../extensions/gsUnstructuredSplines/filedata/surfaces/neon/neon_side_split_fit_"$mesh"_smooth.xml -B $BVP -m $method -r $r -p $p -s $s --plot -N 16 -S $S --write -o "$out"ModalResults/Car_"$mesh"_r"$r"_p"$p"_s"$s"_m"$method" 1> "$out"ModalResults/Car_"$mesh"_r"$r"_p"$p"_s"$s"_m"$method".log 2> "$out"ModalResults/Car_"$mesh"_r"$r"_p"$p"_s"$s"_m"$method".err
#        done
#    done
#done
