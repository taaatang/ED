#!/usr/bin/env bash

# env=env.sh
# set environment
cwd=$(pwd)
# cd ~
# source ./${env}
# cd ${cwd}

JobDir=Job/36/Raman
keyword1=J2idx
keyword2=Jkidx
app=Raman.out
script=run_haswell.sh
input=params_input.txt
for i in 0 15 30
    for j in 0 1 2 3 10 20 30 40 
        do
            dir=${JobDir}/J2_${i}_Jk_${j}
            appDir=${dir}/App
            inputDir=${dir}/Input
            mkdir -p ${appDir}
            mkdir -p ${inputDir}
            cp build/${app}  ${appDir}/${app}
            sed -e "s/${keyword1}/${i}/g;s/${keyword2}/${j}/g" Input/${input} > ${inputDir}/${input}
            sed -e "s/${keyword1}/${i}/g;s/${keyword2}/${j}/g;s/APP/${app}/g" ${script} > ${appDir}/${script}
            cd ${appDir}
            rm -rf job.out
            rm -rf job.err
            sbatch ${script}
            cd ${cwd}
        done
    done