#!/bin/bash
NCORES=4
MEM=2

VPT="/sc/arion/projects/ji_lab/TOOLS/sif/vpt_1.3.sif"

while read RUN_NAME; do

CELLPOSE_DIR="/sc/arion/projects/ji_lab/DATA/MERFISH/merscope_outs/"${RUN_NAME}
IN_DIR="/sc/arion/projects/ji_lab/DATA/MERFISH/merscope_outs/"${RUN_NAME}
WORK_DIR="/sc/arion/projects/ji_lab/DATA/MERFISH/segmentation/NS-Atlas/baysor_0.7.1/"${RUN_NAME}

mkdir -p $WORK_DIR
cd $IN_DIR

for REGION in region_*; do
echo "RUNNING " $REGION "FOR RUN NAME" $RUN_NAME "!"

OUTPUT_DIR=$WORK_DIR/$REGION
mkdir -p $OUTPUT_DIR 
mkdir -p $OUTPUT_DIR/tmp
mkdir -p $OUTPUT_DIR/log
mkdir -p $OUTPUT_DIR/qc

echo "#!/bin/bash
#BSUB -J ${RUN_NAME}.${REGION}.REPORT
#BSUB -P acc_ji_lab
#BSUB -q express
#BSUB -n ${NCORES}
#BSUB -R \"rusage[mem=${MEM}GB]\"
#BSUB -R \"span[hosts=1]\"
#BSUB -W 0:30
#BSUB -Ep 'mailx -s \"LSF Job \${LSB_JOBNAME} \${LSB_JOBID} exited with status \${LSB_JOB_STATUS}\" paula.restrepo@icahn.mssm.edu'
#BSUB -o ${OUTPUT_DIR}/log/report_lsf.%J.log
#BSUB -eo ${OUTPUT_DIR}/log/report_lsf.%J.err
#BSUB -L /bin/bash

set -e # exit whole thing if one command fails

cd ${OUTPUT_DIR} 

ml singularity
export MALLOC_TRIM_THRESHOLD_=0
export DASK_TEMPORARY_DIRECTORY=${OUTPUT_DIR}/tmp

singularity run ${VPT} vpt \\
--verbose --processes ${NCORES} \\
--log-file log/sum_signals.log \\
sum-signals --overwrite \\
--input-images ${IN_DIR}/${REGION}/images \\
--input-micron-to-mosaic ${IN_DIR}/${REGION}/images/micron_to_mosaic_pixel_transform.csv \\
--input-boundaries baysor_micron_space.parquet \\
--output-csv baysor_cell_signals.csv

singularity run ${VPT} vpt \\
--verbose --processes ${NCORES} \\
--log-file log/update_vzg.log \\
generate-segmentation-metrics --overwrite \\
--input-entity-by-gene baysor_cell_by_gene.csv \\
--input-metadata baysor_cell_metadata.csv \\
--input-transcripts baysor_detected_transcripts.csv \\
--input-images ${IN_DIR}/${REGION}/images \\
--input-boundaries baysor_micron_space.parquet \\
--input-micron-to-mosaic ${IN_DIR}/${REGION}/images/micron_to_mosaic_pixel_transform.csv \\
--output-csv qc/segmentation_metrics.csv \\
--input-z-index 4 \\
--output-report qc/segmentation_report_baysor_${REGION}.html \\
--output-clustering qc/ \\
--transcript-count-filter-threshold 10 \\
--volume-filter-threshold 100

">${OUTPUT_DIR}/make_report.lsf
done

done < /sc/arion/projects/ji_lab/DATA/MERFISH/segmentation/NS-Atlas/baysor_0.7.1/run_names.txt