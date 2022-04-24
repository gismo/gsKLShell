#!/bin/bash
# Declare an array of string with type
declare -a Runname="../../../build/bin/strongCoupling_Error_example"
declare -a Filenames=("1p_hyperboloid" "3p_hyperboloid" "4p_hyperboloid" "6p_hyperboloid" "4p_hyperboloid_hole_deg4deg4")
declare -a Methods=(1 3 4)

mkdir -p Output
  
for file in ${Filenames[@]}; do
    for (( p=2; p<5; p++)) do
        for m in ${Methods[@]}; do
            # Define initial refinement level
            if (($m==3))
            then
                declare -a R0=2
            else
                declare -a R0=0
            fi

            # Options per method
            if (($m==1))
            then
                continue
            # elif (($m==2))
            # then

            elif (($m==3))
            then
                if (($p<3))
                then
                    continue
                fi
            elif (($m==4))
            then
                if (($p>2))
                then
                    continue
                fi
            fi

            echo "$Runname" -m $m -p $p -s $(($p-1)) -G ../filedata/pde/"$file"_geom.xml -B ../filedata/pde/"$file"_bvp.xml -C 1e3 -r 6 -R $R0
            
            eval "$Runname" -m $m -p $p -s $(($p-1)) -G ../filedata/pde/"$file"_geom.xml -B ../filedata/pde/"$file"_bvp.xml -C 1e3 -r 6 -R $R0 > Output/"$file"_p"$p"_m"$m".log
        done
    done
done
