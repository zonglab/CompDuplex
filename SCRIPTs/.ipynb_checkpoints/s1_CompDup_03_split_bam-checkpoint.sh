#!/bin/bash
#SBATCH --job-name=03_SB_VAR_03
#SBATCH -N 1         # number of nodes
#SBATCH -n 2         # number of cores
#SBATCH --mem=4gb    # memory pool for all cores
#SBATCH -t 0-12:00   # time (D-HH:MM)
#SBATCH -o 0_log/zlog_03_sp_VAR_03_%j_out.txt # STDOUT
#SBATCH -e 0_log/zlog_03_sp_VAR_03_%j_err.txt # STDERR

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
NORMALBULK=${3}

SCRIPT_DIR=PATH/TO/SCRIPT_DIR
WORK___DIR=PATH/TO/PARENTAL_DIR/${PROJECT_ID}/${SAMPLE__ID}
mkdir -p ${WORK___DIR}/0_log ${WORK___DIR}/0_scripts
cd       ${WORK___DIR}

#################################################################
printf "\n\n source conda"
#################################################################
set +x && conda activate dna && set -x

#################################################################
printf "\n\n split bam"
#################################################################
mkdir -p split
cd split

for i in {1..22} X Y
do
    samtools  view ../06_sorted.bam  ${chr} -b > 07_${chr}.bam
    samtools index 07_${chr}.bam
done

#################################################################
printf "\n\nrun 04_Variant_call.sh"
#################################################################
cd ${WORK___DIR}
for chr in 1 3 5 7 9 11 13 15 17 19 21 X 2 4 6 8 10 12 14 16 18 20 22 Y
do
    echo ${chr}
    sed  "s/VAR_04_1/${SAMPLE__ID}/g" \
         ${SCRIPT_DIR}/s1_CompDup_04_Variant_call.sh | \
    sed    "s/VAR_04_2/chr${chr}/g"  - \
    >      ./0_scripts/s1_CompDup_04_Variant_call_${SAMPLE__ID}_${chr}.sh
    sbatch ./0_scripts/s1_CompDup_04_Variant_call_${SAMPLE__ID}_${chr}.sh  \
           ${PROJECT_ID} ${SAMPLE__ID} ${NORMALBULK}
    sleep  5s     # to avoid some clash of multi-threading setting
done

#################################################################
printf "\n\nrun 05_Call_summary.sh"
#################################################################
cd ${WORK___DIR}
sed   "s/VAR_05/${SAMPLE__ID}/g" \
      ${SCRIPT_DIR}/s1_CompDup_05_Bulk_intersect.sh \
>       ./0_scripts/s1_CompDup_05_Bulk_intersect_${SAMPLE__ID}.sh
sbatch  ./0_scripts/s1_CompDup_05_Bulk_intersect_${SAMPLE__ID}.sh  \
        ${PROJECT_ID} ${SAMPLE__ID} ${NORMALBULK}
