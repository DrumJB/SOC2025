#!/bin/bash
#PBS -l nodes=1:ppn=4
#PBS -q long
#PBS -N jb_test-0

module purge
module load intelmpi/18

export I_MPI_SHM_LMT=shm
ulimit -s unlimited

BIN_DIR=/compass/Shared/Common/IT/SW/SPICE/soroban/spice2/bin
EXE_NAME=spice-2.14-release.bin

INP_DIR=/compass/home/buben/SPICE2/simulations/test
INP_NAME=jb-test01
SCRATCH_DIR=/compass/home/buben/SPICE2/scratch

DATA_DIR=$SCRATCH_DIR/d
BACKUP_DIR=$DATA_DIR-$RANDOM
LOG_NAME=$SCRATCH_DIR/$INP_NAME.log

mkdir -p $SCRATCH_DIR
mkdir -p $DATA_DIR
mkdir -p $BACKUP_DIR

cp -vf $DATA_DIR/$INP_NAME*.mat $BACKUP_DIR
cp -vf $SCRATCH_DIR/$INP_NAME.log $BACKUP_DIR

cp -f $INP_DIR/$INP_NAME.inp $SCRATCH_DIR/$INP_NAME.inp
chmod a+rw $SCRATCH_DIR/$INP_NAME.inp

echo "job $INP_NAME is running">$SCRATCH_DIR/$INP_NAME.running

cd $BIN_DIR

mpirun -np 4 $BIN_DIR/$EXE_NAME -v -i $SCRATCH_DIR/$INP_NAME.inp -t $DATA_DIR/$INP_NAME-t -o $DATA_DIR/$INP_NAME-o > $LOG_NAME

rm -f $SCRATCH_DIR/$INP_NAME.running