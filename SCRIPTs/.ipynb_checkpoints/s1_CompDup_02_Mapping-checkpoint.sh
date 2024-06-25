#!/bin/bash
#SBATCH --job-name=02_Mp_VAR_02
#SBATCH -N 1         # number of nodes
#SBATCH -n 12        # number of cores
#SBATCH --mem=10gb   # memory pool for all cores
#SBATCH -t 0-12:00   # time (D-HH:MM)
#SBATCH -o 0_log/zlog_02_Mp_VAR_02_%j_out.txt # STDOUT
#SBATCH -e 0_log/zlog_02_Mp_VAR_02_%j_err.txt # STDERR

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
printf "\n\n map reads"
# note that we're using samtools-1.12 and biobuilds-2016.04/bwa
#################################################################
seqtk mergepe 03_tmm_R1.fastq.gz 03_tmm_R2.fastq.gz | \
pigz        > 03_merged.fastq.gz

bwa mem -t 8 \
-p   PATH/TO/REFSEQ.fa \
     03_merged.fastq.gz | \
samtools view -bS - \
>    04_alignment.bam

samtools sort -@ 12 \
     04_alignment.bam \
>    05_sorted.bam

#################################################################
printf "\n\n remove optical duplicates"
#################################################################
### duplicates are marked by +1024 to sam flag
### dup  types are tagged as DT:SQ/LB
### optical dup pixel dist: 100-2500-12000
### only remove seq dup by REMOVE_SEQUENCING_DUPLICATES=true
### to remove all or none, REMOVE_DUPLICATE=true/false
#################################################################
picard MarkDuplicates \
    I=05_sorted.bam \
    O=06_rm_opt_dup.bam \
    M=06_rm_opt_dup.txt \
    TAGGING_POLICY=All \
    OPTICAL_DUPLICATE_PIXEL_DISTANCE=12000 \
    REMOVE_SEQUENCING_DUPLICATES=true

#################################################################
printf "\n\n quality filtering"
# note that we're using samtools-1.12
#################################################################
samtools view  -@ 12 -b -q 50 -f 3 06_rm_opt_dup.bam | \
samtools sort  -@ 12 - >           07_sorted.bam
samtools index -@ 12               07_sorted.bam

#################################################################
printf "\n\n run s1_CompDup_03_split_bam.sh"
#################################################################
SCRIPT_DIR=PATH/TO/SCRIPT_DIR

cd ${WORK___DIR}
sed "s/VAR_03/${SAMPLE__ID}/g" \
    ${SCRIPT_DIR}/s1_CompDup_03_split_bam.sh \
>       ./0_scripts/s1_CompDup_03_split_bam_${SAMPLE__ID}.sh
sbatch  ./0_scripts/s1_CompDup_03_split_bam_${SAMPLE__ID}.sh  ${PROJECT_ID} ${SAMPLE__ID} ${NORMALBULK}
