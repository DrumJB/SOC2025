#!/bin/bash
#PBS -l nodes=1:ppn=2
#PBS -q long
#PBS -N 1T10H100He0P05

module purge
module load intelmpi/18

export I_MPI_SHM_LMT=shm
ulimit -s unlimited

BIN_DIR=/compass/Shared/Common/IT/SW/SPICE/soroban/spice2/bin
EXE_NAME=spice-2.14-release.bin

INP_DIR=/compass/home/buben/SPICE2/simulations/gamma-TH0124
INP_NAME=1T10H100He0P05
SCRATCH_DIR=/scratch/buben/gamma-TH0124

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
