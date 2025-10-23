#!/bin/bash
NCORES=4
MEM=2
WALLTIME="0:30"
MERPY="/sc/arion/projects/ji_lab/TOOLS/sif/merpy_latest.sif"

while read RUN_NAME; do

IN_DIR="/sc/arion/projects/ji_lab/DATA/MERFISH/merscope_outs/"${RUN_NAME}
WORK_DIR="/sc/arion/projects/ji_lab/DATA/MERFISH/segmentation/NS-Atlas/baysor_0.7.1/"${RUN_NAME}

mkdir -p $WORK_DIR
cd $IN_DIR

for REGION in region_*; do
echo "RUNNING " $REGION "FOR RUN NAME" $RUN_NAME "!"

OUTPUT_DIR=$WORK_DIR/$REGION
mkdir -p $OUTPUT_DIR 

echo "#!/bin/bash
#BSUB -J ${RUN_NAME}.${REGION}.QC
#BSUB -P acc_ji_lab
#BSUB -q express
#BSUB -n ${NCORES}
#BSUB -R \"rusage[mem=${MEM}GB]\"
#BSUB -R \"span[hosts=1]\"
#BSUB -W ${WALLTIME}
#BSUB -Ep 'mailx -s \"LSF Job \${LSB_JOBNAME} \${LSB_JOBID} exited with status \${LSB_JOB_STATUS}\" paula.restrepo@icahn.mssm.edu'
#BSUB -o ${OUTPUT_DIR}/log/qc_lsf.%J.log
#BSUB -eo ${OUTPUT_DIR}/log/qc_lsf.%J.err
#BSUB -L /bin/bash

set -e # exit whole thing if one command fails

cd ${OUTPUT_DIR} 
ml singularity
export MALLOC_TRIM_THRESHOLD_=0
export DASK_TEMPORARY_DIRECTORY=${OUTPUT_DIR}/tmp

# get qc metrics
singularity run ${MERPY} \\
python /sc/arion/projects/ji_lab/TOOLS/merfish-analysis/bin/py/generate_qc_metrics.py \\
--counts-matrix baysor_cell_by_gene.csv \\
--cell-metadata baysor_cell_metadata.csv \\
--detected-transcripts baysor_detected_transcripts.csv \\
--output-directory ${OUTPUT_DIR}/qc \\
--sample-name ${RUN_NAME}_${REGION}

">${OUTPUT_DIR}/make_qc.lsf
done

done < /sc/arion/projects/ji_lab/DATA/MERFISH/segmentation/NS-Atlas/baysor_0.7.1/run_names.txt