#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  1 18:32:23 2023

@author: jachym
"""

def write_new_file(type, tau, H, He, pot, name="gamma-TPC1223"):
    bash_file=f"""#!/bin/bash
#PBS -l nodes=1:soroban-node-02:ppn=2
#PBS -q long
#PBS -N {type}T{tau}H{H}He{He}P{pot}

module purge
module load intelmpi/18

export I_MPI_SHM_LMT=shm
ulimit -s unlimited

BIN_DIR=/compass/Shared/Common/IT/SW/SPICE/soroban/spice2/bin
EXE_NAME=spice-2.14-release.bin

INP_DIR=/compass/home/buben/SPICE2/simulations/gamma-TPC1223
INP_NAME={type}T{tau}H{H}He{He}P{pot}
SCRATCH_DIR=/scratch/buben/{name}

DATA_DIR=$SCRATCH_DIR/d
LOG_NAME=$SCRATCH_DIR/$INP_NAME.log

mkdir -p $SCRATCH_DIR
mkdir -p $DATA_DIR

cp -f $INP_DIR/$INP_NAME.inp $SCRATCH_DIR/$INP_NAME.inp
chmod a+rw $SCRATCH_DIR/$INP_NAME.inp

echo "job $INP_NAME is running">$SCRATCH_DIR/$INP_NAME.running

cd $BIN_DIR

mpirun -np 2 $BIN_DIR/$EXE_NAME -v -i $SCRATCH_DIR/$INP_NAME.inp -t $DATA_DIR/$INP_NAME-t -o $DATA_DIR/$INP_NAME-o > $LOG_NAME

rm -f $SCRATCH_DIR/$INP_NAME.running
"""
    with open(f"{name}/{type}T{tau}H{H}He{He}P{pot}.sh", "w+") as output:
        output.write(bash_file)

taus = ["05", "10", "20", "35", "50"]
Hs = ["02", "05", "10", "40", "60", "90", "95", "98"]
Hes = ["98", "95", "90", "60", "40", "10", "05", "02"]
Pots = ["01", "03", "06", "10"]
for tau in taus:
    for Pot in Pots:
        for i in range(len(Hs)):
            write_new_file("1", tau, Hs[i], Hes[i], Pot)
            write_new_file("2", tau, Hs[i], Hes[i], Pot)        