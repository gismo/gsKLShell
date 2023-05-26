#!/bin/bash
# Declare an array of string with type
declare -a Runname="./bin/weakCoupling_Modal_example"
declare -a Meshes=("0x0" "1x1" "2x2" "3x3" "4x4")
declare -a BVP="../optional/gsKLShell/filedata/pde/car_pillar_extended_bvp.xml"
declare -a p=2
declare -a s=1
declare -a S=-1e-6
declare -a out="/mnt/sda2/hverhelst/"
declare -a d=1e5

for (( r=0; r<3; r++)) do
    for mesh in ${Meshes[@]}; do
            echo "$Runname" -G ../optional/gsUnstructuredSplines/filedata/surfaces/neon/neon_side_split_pillar_extended_"$mesh"_p2.xml -B $BVP -d $d -r $r -p $p -s $s --plot -N 16 -S $S --write -o "$out"ModalResults/Weak_Car_pillar_extended_"$mesh"_r"$r"_p"$p"_s"$s"_d"$d"
            eval "$Runname" -G ../optional/gsUnstructuredSplines/filedata/surfaces/neon/neon_side_split_pillar_extended_"$mesh"_p2.xml -B $BVP -d $d -r $r -p $p -s $s --plot -N 16 -S $S --write -o "$out"ModalResults/Weak_Car_pillar_extended_"$mesh"_r"$r"_p"$p"_s"$s"_d"$d" 1> "$out"ModalResults/Weak_Car_pillar_extended_"$mesh"_r"$r"_p"$p"_s"$s"_d"$d".log 2> "$out"ModalResults/Weak_Car_pillar_extended_"$mesh"_r"$r"_p"$p"_s"$s"_d"$d".err
    done
done
