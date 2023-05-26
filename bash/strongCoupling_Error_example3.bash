#!/bin/bash
# Declare an array of string with type
declare -a Runname="../../../build/bin/strongCoupling_Error_example"
declare -a Filenames=("4p_hyperboloid_hole_deg4deg4")
declare -a Methods=(2 3 4)

mkdir -p Output
  
for file in ${Filenames[@]}; do
    for (( p=4; p<5; p++)) do
        for (( s=$(($p-2)); s<$p; s++ )) do
            for m in ${Methods[@]}; do
                # Set smoothness
                if (($m==3 && $s>$(($p-2))))
                then
                    continue
                fi

                # Define initial refinement level
                if ((($m==2 || $m==3) && ($p==3))) ## Only pre-refine for p=4 ---- or make an initial knot insertion flag
                then
                    declare -a R0=2
                else
                    declare -a R0=0
                fi

                # Options per method
    #            if (($m==1))
    #            then
    #               if  (("$file"=="3p_hyperboloid" || "$file"=="4p_hyperboloid_hole_deg4deg4"))
    #               then
    #                   continue
    #               fi
                if (($m==2))
                then
                    if (($p<3))
                    then
                        continue
                    fi
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

                echo "$Runname" -m $m -p $p -s $s -G ../filedata/pde/"$file"_geom.xml -B ../filedata/pde/"$file"_bvp.xml -C 1e3 -r 6 -R $R0

                eval "$Runname" -m $m -p $p -s $s -G ../filedata/pde/"$file"_geom.xml -B ../filedata/pde/"$file"_bvp.xml -C 1e3 -r 6 -R $R0 > Output/"$file"_p"$p"_s"$s"_m"$m".log
            done
        done
    done
done