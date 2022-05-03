#!/bin/bash
# Declare an array of string with type
declare -a Runname="../../../build/bin/strongCoupling_ALM_example"
declare -a Outname="../../../build/ArcLengthResults"
declare -a Filenames=("17p_Lshape_4holes" )
declare -a Methods=(4)
  
    # use $GISMO_BUILD_DIR??

mkdir -p Output
for file in ${Filenames[@]}; do
    for (( p=2; p<5; p++)) do
        for m in ${Methods[@]}; do
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

            mkdir -p "$Outname"/"$file"_p"$p"_m"$m"            
        
            echo "$Runname" -m $m -p $p -s $(($p-1)) -r 3 -L 5e-3 -l 2e-1 -N 20 --plot --write -F 1e1 -C 1e3 -G ../filedata/pde/"$file"_geom.xml -B ../filedata/pde/"$file"_bvp.xml -o "$Outname"/"$file"_p"$p"_m"$m"
            
            eval "$Runname" -m $m -p $p -s $(($p-1)) -r 3 -L 5e-3 -l 2e-1 -N 20 --plot --write -F 1e1 -C 1e3 -G ../filedata/pde/"$file"_geom.xml -B ../filedata/pde/"$file"_bvp.xml -o "$Outname"/"$file"_p"$p"_m"$m"  > Output/"$file"_p"$p"_m"$m".log
        done
    done
done
