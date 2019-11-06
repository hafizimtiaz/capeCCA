#!/bin/bash

## The submit script for MATLAB multicore run.
## It runs matlab threads on one node only.
## 1. Specify the job name in '--job-name='.
## 2. Specify the number of requested CPUs/Cores in '--ntasks-per-node=' 
## 3. Load matlab module: "module purge; module load matlab" 
## 4. Submit the script to the cluster through SLURM: "sbatch matlab_multicore_batch.sh"  


#SBATCH --job-name=mnist_del
#SBATCH --output=mnist_del.out
#SBATCH --error=mnist_del.err
#SBATCH --partition=SOE_main
#SBATCH --ntasks-per-node=4
#SBATCH --mail-type=ALL														
#SBATCH --mail-user=hafiz.imtiaz@rutgers.edu											

MYHDIR=$SLURM_SUBMIT_DIR													
MYTMP="/tmp/$USER/$SLURM_JOB_ID"    												
mkdir -p $MYTMP                     												
cp $MYHDIR/mnist_vs_delta.m  $MYTMP	
cp $MYHDIR/myCHindex.m  $MYTMP
cp $MYHDIR/myCCAscore.m  $MYTMP
cp $MYHDIR/mySTindex.m  $MYTMP
cp $MYHDIR/myClusterPerf.m  $MYTMP
cp $MYHDIR/kmeans.m  $MYTMP
cp $MYHDIR/myDistCCA.m  $MYTMP
cp $MYHDIR/MNIST_preprocessed_d100_N50k_new.mat  $MYTMP

cd $MYTMP                           										

matlab -nodisplay -nosplash -r "mnist_vs_delta, exit"

cp $MYTMP/* $MYHDIR                 										
rm -rf  $MYTMP                      										
