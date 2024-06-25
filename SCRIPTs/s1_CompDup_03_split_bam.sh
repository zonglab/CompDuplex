#!/bin/bash
#SBATCH --job-name=03_SB_VAR_03
#SBATCH -N 1         # number of nodes
#SBATCH -n 2         # number of cores
#SBATCH --mem=4gb    # memory pool for all cores
#SBATCH -t 0-12:00   # time (D-HH:MM)
#SBATCH -o 0_logs/zlog_03_sp_VAR_03_%j_out.txt # STDOUT
#SBATCH -e 0_logs/zlog_03_sp_VAR_03_%j_err.txt # STDERR

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

WORK__PATH=${1}     # absolute path to the final working directory
NORM__BULK=${2}     # path to the normal bulk, includes both WGS bam and variant calls
SAMPLE__ID=$(basename "${WORK__PATH}")

SCRIPT_DIR=PATH_TO_SCRIPT_DIR
mkdir -p ${WORK__PATH}/0_logs ${WORK__PATH}/0_scripts
cd       ${WORK__PATH}

#################################################################
printf "\n\n source conda"
#################################################################
source "PATH_TO_CONDA/etc/profile.d/conda.sh"
set +x && conda activate <ENVIRONMENT> && set -x

#################################################################
printf "\n\n split bam"
#################################################################
mkdir -p 1_split
cd 1_split

for chr in {1..22} X Y
do
    samtools  view ../07_sorted.bam  chr${chr} -b > 07_chr${chr}.bam
    samtools index 07_chr${chr}.bam
done

#################################################################
printf "\n\nrun 04_Variant_call.sh"
#################################################################
cd ${WORK__PATH}
### for chr in 1 3 5 7 9 11 13 15 17 19 21 X 2 4 6 8 10 12 14 16 18 20 22 Y
for chr in 21
do
    echo ${chr}
    sed  "s/VAR_04_1/${SAMPLE__ID}/g" \
         ${SCRIPT_DIR}/s1_CompDup_04_Variant_call.sh | \
    sed    "s/VAR_04_2/chr${chr}/g"  - \
    >      ./0_scripts/s1_CompDup_04_Variant_call_${SAMPLE__ID}_${chr}.sh
    sbatch ./0_scripts/s1_CompDup_04_Variant_call_${SAMPLE__ID}_${chr}.sh  \
           ${WORK__PATH}
    sleep  5s     # to avoid some clash of multi-threading setting
done

#################################################################
printf "\n\nrun 05_Call_summary.sh"
#################################################################
cd ${WORK__PATH}
sed   "s/VAR_05/${SAMPLE__ID}/g" \
      ${SCRIPT_DIR}/s1_CompDup_05_Bulk_intersect.sh \
>       ./0_scripts/s1_CompDup_05_Bulk_intersect_${SAMPLE__ID}.sh
sbatch  ./0_scripts/s1_CompDup_05_Bulk_intersect_${SAMPLE__ID}.sh  \
        ${WORK__PATH} ${NORM__BULK}