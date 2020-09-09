#!/bin/bash
#
# SLURM control file automatically created by launch
#
# Created on: 2020-09-08 09:10:28.252651
# Launching single command: -j test_job.txt
#
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -J test
#SBATCH -o test.o%j
#SBATCH -p normal
#SBATCH -t 3:00:00
#SBATCH -A LewPea_MRI_Analysis

umask 2

echo " Starting at $(date)"
start=$(date +%s)
echo " WORKING DIR: /opt/tacc_mounts/home1/05426/ach3377/gPPI/"
echo " JOB ID:      $SLURM_JOB_ID"
echo " JOB NAME:    $SLURM_JOB_NAME"
echo " NODES:       $SLURM_NODELIST"
echo " N NODES:     $SLURM_NNODES"
echo " N TASKS:     $SLURM_NTASKS"
set -x
-j test_job.txt
set +x
echo " "
echo " Job complete at $(date)"
echo " "
finish=$(date +%s)
printf "Job duration: %02d:%02d:%02d (%d s)
" $(((finish-start)/3600)) $(((finish-start)%3600/60)) $(((finish-start)%60)) $((finish-start))
