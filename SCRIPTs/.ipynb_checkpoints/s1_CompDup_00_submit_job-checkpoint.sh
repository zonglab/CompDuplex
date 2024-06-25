PROJECT_ID=${1}     # grouping of samples
SAMPLE__ID=${2}     # sample ID
FASTQ_PATH=${3}     # path to the raw fastq files
NORMALBULK=${4}     # name of the normal bulk file

WORK___DIR=/storage/zongarchive/Yang/2_duplex_bulk/1_SMaHT/${SAMPLE__ID}
mkdir -p ${WORK___DIR}/0_log ${WORK___DIR}/0_scripts
cd       ${WORK___DIR}

#################################################################
printf "run s1_CompDup_01_PreProcess.sh"
#################################################################
SCRIPT_DIR=PATH/TO/SCRIPT_DIR

sed "s/VAR_01/${SAMPLE__ID}/g" \
    /storage/zongarchive/Yang/2_duplex_bulk/0_scripts/s1_CompDup_01_PreProcess.sh \
>       ./0_scripts/s1_CompDup_01_PreProcess_${SAMPLE__ID}.sh
sbatch  ./0_scripts/s1_CompDup_01_PreProcess_${SAMPLE__ID}.sh ${PROJECT_ID} ${SAMPLE__ID} ${DATA__PATH} ${NORMALBULK}

