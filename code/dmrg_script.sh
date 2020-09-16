#!/bin/bash
#SBATCH -N 1
#SBATCH --tasks-per-node 1
#SBATCH -A quant
#SBATCH -p public
#SBATCH -o dmrg.out
#SBATCH -e dmrg.err
#SBATCH -J dmrg
#SBATCH --mem=3750
#SBATCH -t 1:00:00
#####
#
#   Load the modules needed to run block-1.5.3
#
####
source /home/configfiles/bashrc.default
module purge
module load slurm/slurm
module load intel/intel16
module load openmpi/openmpi-2.0.i17
module load /usr/local/quant/modules/boost-1-55-0.i16
####
#
#   dir where the data to run block are located : $DATADIR
#
####

export DATADIR=$PWD
export EXE=/usr/local/quant/block-1.5.3/block.spin_adapted

####
#
#   define a $SCRATCH with the jobid in the $WORKDIR
#   copy the data (.conf) and FCIDUMP file in the $WORKDIR
#   run block on 1 core the output will be created in the $DATADIR
#
####

export WORKDIR=/b/work/quant/bsenjean
export SCRATCH=${WORKDIR}/block-${SLURM_JOBID}
mkdir $SCRATCH
cd $SCRATCH
cp ${DATADIR}/FCIDUMP .
cp ${DATADIR}/dmrg_onepdm.conf .
cp ${DATADIR}/dmrg_twopdm.conf .
mpirun -np 1 $EXE dmrg_onepdm.conf > ${DATADIR}/dmrg.out 2>${DATADIR}/dmrg.err
wait
mpirun -np 1 $EXE dmrg_twopdm.conf > ${DATADIR}/dmrg.out 2>${DATADIR}/dmrg.err
wait

####
#
#  erase the scratch dir
#
####
cp ${SCRATCH}/node0/spatial_onepdm.0.0.txt ${DATADIR}
cp ${SCRATCH}/node0/spatial_twopdm.0.0.txt ${DATADIR}
cd  ${SCRATCH}
rm -rf $SCRATCH
