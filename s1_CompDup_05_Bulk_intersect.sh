#!/bin/bash
#SBATCH --job-name=05_SM_VAR_05
#SBATCH -N 1         # number of nodes
#SBATCH -n 2         # number of cores
#SBATCH --mem=10gb   # memory pool for all cores
#SBATCH -t 0-12:00   # time (D-HH:MM)
#SBATCH -o 0_log/zlog_05_cs_VAR_05_%j_out.txt # STDOUT
#SBATCH -e 0_log/zlog_05_cs_VAR_05_%j_err.txt # STDERR

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
NORML_BULK=${3}

SCRIPT_DIR=/storage/zongarchive/Yang/2_duplex_bulk/0_scripts/
WORK___DIR=/storage/zongarchive/Yang/2_duplex_bulk/${PROJECT_ID}/${SAMPLE__ID}
BULK__PATH=/storage/zongarchive/Yang/1_normal_bulk/${NORML_BULK}
mkdir -p ${WORK___DIR}/0_log ${WORK___DIR}/0_scripts
cd       ${WORK___DIR}

#################################################################
printf "\n\n source conda"
#################################################################
source "/storage/zong/home/u244338/miniconda3/etc/profile.d/conda.sh"
set +x && conda activate dna && set -x

#################################################################
printf "\n\n Filter for detected germline mutations"
#################################################################
python Detected_Mutation_Demultiplex.py \
${WORK___DIR} \
${NORML_BULK} \
${BULK__PATH}

#################################################################
printf "\n\n Wait for bam_add_META_tag to finish"
#################################################################
## Wait for bam_add_META_tag to finish
for chr in 1 3 5 7 9 11 13 15 17 19 21 X 2 4 6 8 10 12 14 16 18 20 22 Y
do
    FILE=split/complete/${chr}.txt
    while true; do
        if [[ -f "$FILE" ]]; then
           echo "$FILE exists." && break
          else
           printf '$(date +"%T") \n waiting for $FILE' && sleep 3m
         fi
    done
done

#################################################################
printf "\n\n Get Germline and Somatic"
#################################################################
BULK__BAM=${BULK__PATH}/05_${NORML_BULK}_add_RG.bam
BULK__TOT=${BULK__PATH}/15_${NORML_BULK}_fin_allvar.vcf
BULK__HET=${BULK__PATH}/15_${NORML_BULK}_fin_hetero.vcf
BULK__HOM=${BULK__PATH}/15_${NORML_BULK}_fin___homo.vcf

bulk_comp=${SCRIPT_DIR}/py_gtbulksearch_current.py
python3 ${bulk_comp} 15_3_not_germlines.bdy ${BULK__BAM} 15_4  0

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
