#!/bin/bash
#SBATCH --job-name=04_VC_VAR_04_1_VAR_04_2
#SBATCH -N 1         # number of nodes
#SBATCH -n 1         # number of cores
#SBATCH --mem=4gb    # memory pool for all cores
#SBATCH -t 0-12:00   # time (D-HH:MM)
#SBATCH -o 0_logs/zlog_04_VC_VAR_04_1_VAR_04_2_%j_out.txt # STDOUT
#SBATCH -e 0_logs/zlog_04_VC_VAR_04_1_VAR_04_2_%j_err.txt # STDERR

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
chr=VAR_04_2
SAMPLE__ID=$(basename "${WORK__PATH}")

SCRIPT_DIR=PATH_TO_SCRIPT_DIR
mkdir -p ${WORK__PATH}/0_logs ${WORK__PATH}/0_scripts ${WORK__PATH}/1_split
cd       ${WORK__PATH}/1_split

#################################################################
printf "\n\n source conda"
#################################################################
source "PATH_TO_CONDA/etc/profile.d/conda.sh"
set +x && conda activate <ENVIRONMENT> && set -x

#################################################################
printf "\n\n variant calling"
# note that we're using samtools-1.12 and biobuilds-2016.04/bwa
#################################################################
#export OMP_NUM_THREADS=1
ref_seq=PATH_TO_REFSEQ/hg19.fa
samtools index 07_${chr}.bam

time \
samtools mpileup --ff QCFAIL \
-uf     ${ref_seq} \
        07_${chr}.bam | \
bcftools filter -i '(I16[2]+I16[3]) > 3' \
>       08_${chr}_raw.vcf

grep ^# 08_${chr}_raw.vcf \
>       08_${chr}_header.txt

time /opt/biobuilds-2016.04/bin/bcftools norm -m-any \
-f      ${ref_seq} \
        08_${chr}_raw.vcf | \
grep -v "<*>" \
>       09_${chr}_norm.vcf

#################################################################
printf "\n\n remove black-list"
#################################################################
bedtools intersect -v \
-a  09_${chr}_norm.vcf  \
-b  ${SCRIPT_DIR}/hg19_bad_merged.bed \
>   10_${chr}__clean.bdy

cat 10_${chr}__clean.bdy  |  \
awk '{print $1"\t"$2-10"\t"$2+10}' \
>   11_${chr}__local.bed

seqtk subseq -t \
    ${ref_seq} \
    11_${chr}__local.bed \
>   12_${chr}__local.seq

get_bulk_search=${SCRIPT_DIR}/py_localseqfilter_current.py
python3 ${get_bulk_search} \
    10_${chr}__clean.bdy \
    12_${chr}__local.seq \
    13_${chr}_fltred.bdy

#################################################################
printf "\n\n a4s2 calling"
#################################################################
aMsN_script=${SCRIPT_DIR}/py_Dup_Mutation_call_aMsN_current.py

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
