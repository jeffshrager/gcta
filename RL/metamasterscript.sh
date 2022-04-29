#!/bin/bash
params=$(ls 'exparams')
timestamp=$(date +%s)
echo $timestamp
mkdir 'experiments/'$timestamp
for p in $params
do
    echo 'exparams/'$p
    cp 'exparams/'$p './parameters.json'
    python3 scikit-master-script.py
done
