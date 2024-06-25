#!/bin/bash
#SBATCH --job-name=01_PP_VAR_01
#SBATCH -N 1         # number of nodes
#SBATCH -n 1         # number of cores
#SBATCH --mem=1gb    # memory pool for all cores
#SBATCH -t 0-12:00   # time (D-HH:MM)
#SBATCH -o 0_log/zlog_01_pp_VAR_01_%j_out.txt # STDOUT
#SBATCH -e 0_log/zlog_01_pp_VAR_01_%j_err.txt # STDERR

#################################################################
printf "\n\n config"
#################################################################
export LC_ALL=en_US.utf-8
export   LANG=en_US.utf-8
printf "\n\n run time: $(date)"
set -x

#################################################################
printf "\n\n config run"
#################################################################

PROJECT_ID=${1}
SAMPLE__ID=${2}
FASTQ_PATH=${3}
NORMALBULK=${4}

WORK___DIR=/PATH/TO/WORK_DIR/${SAMPLE__ID}
mkdir -p ${WORK___DIR}/0_log ${WORK___DIR}/0_scripts
cd       ${WORK___DIR}

#################################################################
printf "\n\n source conda"
#################################################################
set +x && conda activate dna && set -x

#################################################################
printf "\n\n cut sequence adaptor"
#################################################################
cutadapt -j 2 -m 50 \
-a  CTGTCTCTTATACACATCT \
-A  CTGTCTCTTATACACATCT \
-o  03_tmm_R1.fastq.gz \
-p  03_tmm_R2.fastq.gz \
    ${FASTQ_PATH}/*_R1*.f*q.gz \
    ${FASTQ_PATH}/*_R2*.f*q.gz \
> trimming_report.txt

#rm 02_ext_R1.fastq.gz 02_ext_R2.fastq.gz

#################################################################
printf "\n\n run s1_CompDup_02_Mapping.sh"
#################################################################
SCRIPT_DIR=PATH/TO/SCRIPT_DIR

cd ${WORK___DIR}
sed "s/VAR_02/${SAMPLE__ID}/g" \
    ${SCRIPT_DIR}/s1_CompDup_02_Mapping.sh \
>       ./0_scripts/s1_CompDup_02_Mapping_${SAMPLE__ID}.sh
sbatch  ./0_scripts/s1_CompDup_02_Mapping_${SAMPLE__ID}.sh ${PROJECT_ID} ${SAMPLE__ID} ${NORMALBULK}
