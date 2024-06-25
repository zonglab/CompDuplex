#!/bin/bash
#SBATCH --job-name=01_PP_VAR_01
#SBATCH -N 1         # number of nodes
#SBATCH -n 1         # number of cores
#SBATCH --mem=1gb    # memory pool for all cores
#SBATCH -t 0-12:00   # time (D-HH:MM)
#SBATCH -o 0_logs/zlog_01_pp_VAR_01_%j_out.txt # STDOUT
#SBATCH -e 0_logs/zlog_01_pp_VAR_01_%j_err.txt # STDERR

#################################################################
printf "config"
#################################################################
export LC_ALL=en_US.utf-8
export   LANG=en_US.utf-8
printf "run time: $(date)"
set -x

#################################################################
printf "config run \n"
#################################################################

WORK__PATH=${1}     # absolute path to the final working directory
FASTQ_PREF=${2}     # path to the raw fastq files and prefix e.g. </PATH/TO/FASTQ> then R1_001.fq.gz
NORM__BULK=${3}     # path to the normal bulk, includes both WGS bam and variant calls
SAMPLE__ID=$(basename "${WORK__PATH}")

SCRIPT_DIR=PATH_TO_SCRIPT_DIR
mkdir -p ${WORK__PATH}/0_logs ${WORK__PATH}/0_scripts
cd       ${WORK__PATH}

#################################################################
printf "source conda \n"
#################################################################
source "PATH_TO_CONDA/etc/profile.d/conda.sh"
set +x && conda activate <ENVIRONMENT> && set -x

#################################################################
printf "cut sequence adaptor \n"
#################################################################
cutadapt -j 2 -m 50 \
-a  CTGTCTCTTATACACATCT \
-A  CTGTCTCTTATACACATCT \
-o  03_tmm_R1.fastq.gz \
-p  03_tmm_R2.fastq.gz \
    ${FASTQ_PREF}*_R1*.f*q.gz \
    ${FASTQ_PREF}*_R2*.f*q.gz \
> trimming_report.txt

#################################################################
printf "run s1_CompDup_02_Mapping.sh \n"
#################################################################
cd ${WORK__PATH}
sed "s/VAR_02/${SAMPLE__ID}/g" \
    ${SCRIPT_DIR}/s1_CompDup_02_Mapping.sh \
>       ./0_scripts/s1_CompDup_02_Mapping_${SAMPLE__ID}.sh
sbatch  ./0_scripts/s1_CompDup_02_Mapping_${SAMPLE__ID}.sh ${WORK__PATH} ${NORM__BULK}
