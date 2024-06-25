#!/bin/bash
#SBATCH --job-name=05_SM_VAR_05
#SBATCH -N 1         # number of nodes
#SBATCH -n 2         # number of cores
#SBATCH --mem=10gb   # memory pool for all cores
#SBATCH -t 0-12:00   # time (D-HH:MM)
#SBATCH -o 0_logs/zlog_05_cs_VAR_05_%j_out.txt # STDOUT
#SBATCH -e 0_logs/zlog_05_cs_VAR_05_%j_err.txt # STDERR

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
printf "\n\n Filter for detected germline mutations"
#################################################################
python ${SCRIPT_DIR}/py_Detected_Mutation_Demultiplex.py \
    ${WORK__PATH} \
    ${NORM__BULK}

#################################################################
printf "\n\n Get Germline and Somatic"
#################################################################
BULK__BAM=${NORM__BULK}/05_add_RG.bam

python3 ${SCRIPT_DIR}/py_gtbulksearch_current.py \
    15_3_not_germlines.bdy ${BULK__BAM} 15_4  0

#################################################################
printf "\n\n Filter for dbSNP"
#################################################################
# annotating SNP
SnpSift=PATH/TO/snpEff_dir/SnpSift.jar
Snp__db=PATH/TO/dbSNP.vcf

# add a header
cat ${vcf_header} \
    15_4_somatc.bdy \
>   15_4_somatc.vcf

# annotate dbSNP mutations
java -jar ${SnpSift} annotate  \
    ${Snp__db}  \
    15_4_somatc.vcf \
>   15_5_somatc_sifted.vcf

# remove dbSNP mutations
grep -v rs \
    15_5_somatc_sifted.vcf \
>   15_6_somatc_fltred.vcf
