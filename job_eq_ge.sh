#!/usr/bin/env bash

# env=env.sh
# set environment
cwd=$(pwd)
# cd ~
# source ./${env}
# cd ${cwd}

JobDir=../Job/sq2x2_2u2d/eq/k0
app=main.out
script=run_sherlock.sh

appDir=${JobDir}
mkdir -p ${appDir}
inputDir=${appDir}/Input
mkdir -p ${inputDir}

cp Input/path.txt ${inputDir}/path.txt
cp Input/Hubbard.txt ${inputDir}/Hubbard.txt
cp Input/lattice.txt ${inputDir}/lattice.txt
cp Input/pulse.txt ${inputDir}/pulse.txt
cp Input/measure.txt $inputDir}/measure.txt

cp ${script} ${appDir}/${script}

cd ${appDir}
rm -rf job.out
rm -rf job.err
sbatch ${script}
cd ${cwd}
