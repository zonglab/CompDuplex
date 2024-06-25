#################################################################
printf "config run \n"
#################################################################
WORK__PATH=${1}     # absolute path to the final working directory
FASTQ_PREF=${2}     # path to the raw fastq files and prefix e.g. </PATH/TO/FASTQ> then R1_001.fq.gz
NORM__BULK=${3}     # path to the normal bulk, includes both WGS bam and variant calls
SAMPLE__ID=$(basename "${WORK__PATH}")

SCRIPT_DIR=PATH_TO_SCRIPT_DIR
mkdir -p ${WORK__PATH}/0_logs ${WORK__PATH}/0_scripts ${WORK__PATH}/1_split
cd       ${WORK__PATH}

#################################################################
printf "run s1_CompDup_01_PreProcess.sh \n"
#################################################################
sed "s/VAR_01/${SAMPLE__ID}/g" \
    ${SCRIPT_DIR}/s1_CompDup_01_PreProcess.sh \
>       ./0_scripts/s1_CompDup_01_PreProcess_${SAMPLE__ID}.sh
sbatch  ./0_scripts/s1_CompDup_01_PreProcess_${SAMPLE__ID}.sh \
        ${WORK__PATH} ${FASTQ_PREF} ${NORM__BULK}

