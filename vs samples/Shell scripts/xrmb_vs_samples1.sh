#!/bin/bash

## The submit script for MATLAB multicore run.
## It runs matlab threads on one node only.
## 1. Specify the job name in '--job-name='.
## 2. Specify the number of requested CPUs/Cores in '--ntasks-per-node=' 
## 3. Load matlab module: "module purge; module load matlab" 
## 4. Submit the script to the cluster through SLURM: "sbatch matlab_multicore_batch.sh"  


#SBATCH --job-name=xrmb_sample1
#SBATCH --output=xrmb_sample1.out
#SBATCH --error=xrmb_sample1.err
#SBATCH --partition=SOE_main
#SBATCH --ntasks-per-node=4
#SBATCH --mail-type=ALL														
#SBATCH --mail-user=hafiz.imtiaz@rutgers.edu											

MYHDIR=$SLURM_SUBMIT_DIR													
MYTMP="/tmp/$USER/$SLURM_JOB_ID"    												
mkdir -p $MYTMP                     												
cp $MYHDIR/xrmb_vs_samples1.m  $MYTMP
cp $MYHDIR/myCHindex.m  $MYTMP
cp $MYHDIR/myCCAscore.m  $MYTMP
cp $MYHDIR/mySTindex.m  $MYTMP
cp $MYHDIR/myClusterPerf.m  $MYTMP
cp $MYHDIR/kmeans.m  $MYTMP
cp $MYHDIR/myDistCCA.m  $MYTMP
cp $MYHDIR/XRMB_preprocessed_d50_p10_new.mat  $MYTMP
cp $MYHDIR/XRMB_preprocessed_d50_p20_new.mat  $MYTMP
cp $MYHDIR/XRMB_preprocessed_d50_p30_new.mat  $MYTMP
cp $MYHDIR/XRMB_preprocessed_d50_p40_new.mat  $MYTMP
cp $MYHDIR/XRMB_preprocessed_d50_p50_new.mat  $MYTMP

cd $MYTMP                           										

matlab -nodisplay -nosplash -r "xrmb_vs_samples1, exit"

cp $MYTMP/* $MYHDIR                 										
rm -rf  $MYTMP                      										