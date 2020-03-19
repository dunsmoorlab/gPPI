#!/bin/bash
#had to install dev of freesurfer and matlab runtime on my own

export FREESURFER_HOME=/scratch/05426/ach3377/freesurfer
source $FREESURFER_HOME/SetUpFreeSurfer.sh
export SUBJECTS_DIR=/scratch/05426/ach3377/fc-bids/derivatives/freesurfer
export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=48

#for SUBI in 001 002 003 004 005 006 007 008 010 012 013 014 015 016 017 018 019 020 021 023 024 025 026 101 102 103 104 105 106 107 108 109 110 111 112 113 114 115 116 117 118 120 121 122 123 124 125
for SUBI in 002
do
    segmentHA_T1.sh "sub-FC${SUBI}"
done
