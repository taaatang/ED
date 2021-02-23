#!/usr/bin/env bash

# env=env.sh
# set environment
cwd=$(pwd)
# cd ~
# source ./${env}
# cd ${cwd}

JobDir=../Job/sq2x2_2u2d/pump/k0
keyword1=inputw
keyword2=inputflu
app=main.out
script=run_sherlock.sh
input=pulse.txt

for i in $(seq 0 0.1 5)
    do
        for j in $(seq 0 0.2 2)
            do
                appDir=${JobDir}/w_${i}/flu_${j}
                mkdir -p ${appDir}
                inputDir=${appDir}/Input
                mkdir -p ${inputDir}
                cp Input/path.txt ${inputDir}/path.txt
                cp Input/Hubbard.txt ${inputDir}/Hubbard.txt
                cp Input/lattice.txt ${inputDir}/lattice.txt
                cp Input/measure.txt ${inputDir}/measure.txt

                sed -e "s/${keyword1}/${i}/g;s/${keyword2}/${j}/g" Input/${input} > ${inputDir}/${input}
                sed -e "s/${keyword1}/${i}/g;s/${keyword2}/${j}/g" ${script} > ${appDir}/${script}
                cd ${appDir}
                rm -rf job.out
                rm -rf job.err
                sbatch ${script}
                cd ${cwd}
            done
    done