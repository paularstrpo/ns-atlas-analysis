#!/bin/bash

RUN_NAME="202404241001_D165--7-12--04242024_VMSC02301"
NCORES=32
MEM=8
IN_DIR="/sc/arion/projects/ji_lab/DATA/MERFISH/merscope_outs/"${RUN_NAME}
WORK_DIR="/sc/arion/projects/ji_lab/DATA/MERFISH/segmentation/D165-benchmark/baysor_0.7.1/"${RUN_NAME}
BAYSOR="/sc/arion/projects/ji_lab/TOOLS/sif/baysor_v0.7.1.sif"
CONFIG="/sc/arion/projects/ji_lab/DATA/MERFISH/segmentation/D165-benchmark/baysor_0.7.1/input_config_conf_0.95_autoscale.toml"
VPT="/sc/arion/projects/ji_lab/TOOLS/sif/vpt_1.3.sif"
MERPY="/sc/arion/projects/ji_lab/TOOLS/sif/merpy_latest.sif"
CELLPOSE_DIR="/sc/arion/projects/ji_lab/DATA/MERFISH/merscope_outs/"${RUN_NAME}

mkdir -p $WORK_DIR
cd $IN_DIR

for REGION in region_*; do
echo "RUNNING " $REGION "FOR RUN NAME" $RUN_NAME "!"

OUTPUT_DIR=$WORK_DIR/$REGION
mkdir -p $OUTPUT_DIR 

echo "#!/bin/bash
#BSUB -J ${RUN_NAME}.${REGION}.runBAYSOR 
#BSUB -P acc_ji_lab
#BSUB -q express
#BSUB -n ${NCORES}
#BSUB -R \"rusage[mem=${MEM}GB]\"
#BSUB -R \"span[hosts=1]\"
#BSUB -W 2:00
#BSUB -Ep 'mailx -s \"LSF Job \${LSB_JOBNAME} \${LSB_JOBID} exited with status \${LSB_JOB_STATUS}\" paula.restrepo@icahn.mssm.edu'
#BSUB -o ${OUTPUT_DIR}/baysor_lsf.%J.log
#BSUB -eo ${OUTPUT_DIR}/baysor_lsf.%J.err
#BSUB -L /bin/bash

set -e # exit whole thing if one command fails

cd ${OUTPUT_DIR} 
mkdir -p tmp
mkdir -p log
mkdir -p qc
ml singularity
export MALLOC_TRIM_THRESHOLD_=0
export DASK_TEMPORARY_DIRECTORY=${OUTPUT_DIR}/tmp

singularity run ${MERPY} \\
 python /bin/prepare_baysor_inputs.py \\
 --cellpose-boundaries ${CELLPOSE_DIR}/${REGION}/cellpose_micron_space.parquet \\
 --cellpose-transcripts ${CELLPOSE_DIR}/${REGION}/detected_transcripts.csv \\
 --output-directory ${OUTPUT_DIR} 

# Step 2: Run Baysor
export JULIA_NUM_THREADS=${NCORES}
singularity run ${BAYSOR} run \\
--plot \\
--output baysor_segmentation.csv \\
--config ${CONFIG} \\
input_transcripts.csv :cell_id

# Step 3: post process baysor polygons, will output temp.geojson 
singularity run ${MERPY} \\
python /sc/arion/projects/ji_lab/TOOLS/baysor/baysor2vizgen.py \\
--input baysor_segmentation_polygons_2d.json

#Step 4: convert baysor geometries
singularity run ${VPT} vpt \\
--verbose --processes ${NCORES} \\
--log-file log/convert_baysor_geometry.log  \\
convert-geometry --overwrite \\
--input-boundaries temp.geojson \\
--convert-to-3D \\
--number-z-planes 8 \\
--spacing-z-planes 1.5 \\
--output-boundaries baysor_micron_space.parquet 

#Step 4: Partition Transcripts
singularity run ${VPT} vpt \\
--verbose --processes ${NCORES} \\
--log-file log/partition_transcripts.log \\
partition-transcripts --overwrite \\
--input-boundaries baysor_micron_space.parquet \\
--input-transcripts input_transcripts.csv \\
--output-entity-by-gene baysor_cell_by_gene.csv \\
--output-transcripts baysor_detected_transcripts.csv

#Get metadata
singularity run ${VPT} vpt \\
--verbose --processes ${NCORES} \\
--log-file log/derive_entity_metadata.log \\
derive-entity-metadata --overwrite \\
--input-boundaries baysor_micron_space.parquet \\
--input-entity-by-gene baysor_cell_by_gene.csv \\
--output-metadata baysor_cell_metadata.csv

# make new VZG
singularity run ${VPT} vpt \\
--verbose --processes ${NCORES} \\
--log-file log/update_vzg.log \\
update-vzg --overwrite \\
--input-vzg ${IN_DIR}/${REGION}/${RUN_NAME}_${REGION}.vzg \\
--input-boundaries baysor_micron_space.parquet \\
--input-entity-by-gene baysor_cell_by_gene.csv \\
--input-metadata baysor_cell_metadata.csv \\
--output-vzg "${OUTPUT_DIR}/${RUN_NAME}_${REGION}_baysor.vzg" \\
--temp-path tmp/ 


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

">${OUTPUT_DIR}/run_baysor.lsf

done