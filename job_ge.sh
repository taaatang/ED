#!/usr/bin/env bash

N1=3
keyword=N2
app=Spectra.out
script=run.sh
input=spectra_input.txt
for i in 11
    do
     dir=${N1}by${i}/spectra
     mkdir -p ${dir}
     cp ${app}  ${dir}/${app}
     sed -e "s/${keyword}/${i}/g" ${script} > ${dir}/${script}
     sed -e "s/${keyword}/${i}/g" ${input} > ${dir}/${input}
     cd ${dir}
     rm -rf job.out
     rm -rf job.err
     sbatch run.sh
     cd ../../
    done