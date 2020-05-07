#!/usr/bin/env bash

env=env.sh
# set environment
cwd=$(pwd)
cd ~
source ./${env}
cd ${cwd}

JobDir=Job/sq8_4u4d/Basis
keyword=kIndex
app=build/genBasis.out
script=run_haswell.sh
input=Input/symm_input.txt
for i in 0 1 2 3 4 5 6 7
    do
        dir=${JobDir}/k${i}
        appDir=dir/App
        inputDir=dir/Input
        mkdir -p ${appDir}
        mkdir -p ${inputDir}
        cp ${app}  ${appDir}/${app}
        sed -e "s/${keyword}/${i}/g" ${input} > ${inputDir}/${input}
        sed -e "s/${keyword}/${i}/g" ${script} > ${dir}/${script}
        cd ${dir}
        rm -rf job.out
        rm -rf job.err
        sbatch ${script}
        cd ${cwd}
    done