#!/bin/bash
#SBATCH --job-name=04_VC_VAR_04_1_VAR_04_2
#SBATCH -N 1         # number of nodes
#SBATCH -n 1         # number of cores
#SBATCH --mem=4gb    # memory pool for all cores
#SBATCH -t 0-12:00   # time (D-HH:MM)
#SBATCH -o 0_log/zlog_04_VC_VAR_04_1_VAR_04_2_%j_out.txt # STDOUT
#SBATCH -e 0_log/zlog_04_VC_VAR_04_1_VAR_04_2_%j_err.txt # STDERR

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
chr=VAR_04_2
PROJECT_ID=${1}
SAMPLE__ID=${2}

SCRIPT_DIR=PATH/TO/SCRIPT_DIR
WORK___DIR=PATH/TO/PARENTAL_DIR/${PROJECT_ID}/${SAMPLE__ID}/1_split
mkdir -p ${WORK___DIR}/0_log ${WORK___DIR}/0_scripts
cd       ${WORK___DIR}

#################################################################
printf "\n\n source conda"
#################################################################
set +x && conda activate dna && set -x

#################################################################
printf "\n\n variant calling"
# note that we're using samtools-1.12 and biobuilds-2016.04/bwa
#################################################################
#export OMP_NUM_THREADS=1
samtools index 07_${chr}.bam

time \
samtools mpileup --ff QCFAIL \
-uf /storage/zong/ref/hg19/hg19.fa \
        07_${chr}.bam | \
bcftools filter -i '(I16[2]+I16[3]) > 3' \
>       08_${chr}_raw.vcf

grep ^# 08_${chr}_raw.vcf \
>       08_${chr}_header.txt

time bcftools norm -m-any \
-f  /storage/zong/ref/hg19/hg19.fa \
        08_${chr}_raw.vcf | \
grep -v "<*>" \
>       09_${chr}_norm.vcf

#################################################################
printf "\n\n remove black-list"
#################################################################


bedtools intersect -v \
-a  09_${chr}_norm.vcf  \
-b  hg19_bad_merged.bed \
>   10_${chr}__clean.bdy

cat 10_${chr}__clean.bdy  |  \
awk '{print $1"\t"$2-10"\t"$2+10}' \
>   11_${chr}__local.bed

/storage/zong/opt/seqtk/seqtk subseq -t \
    /storage/zong/ref/hg19/hg19.fa \
    11_${chr}__local.bed \
>   12_${chr}__local.seq

get_bulk_search=${SCRIPT_DIR}/SCPT_3_local_seq_filter/localseqfilter_py3_24.05.01.py
python3 ${get_bulk_search} \
    10_${chr}__clean.bdy \
    12_${chr}__local.seq \
    13_${chr}_fltred.bdy

#################################################################
printf "\n\n a4s2 calling"
#################################################################
aMsN_script=${SCRIPT_DIR}/SCPT_4_aMsN_call/Dup_Mutation_call_aMsN_current.py

rm      13_${chr}_fltred.bdy.gz
pigz -k 13_${chr}_fltred.bdy
python3 ${aMsN_script} \
    ./ \
    13_${chr}_fltred.bdy.gz \
    07_${chr}.bam \
    4 2 \
    14_${chr}_a4s2_call.txt 150

mkdir -p complete
echo "Complete" > complete/${chr}.txt
