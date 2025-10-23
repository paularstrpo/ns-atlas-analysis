#!/bin/bash
NCORES=4
MEM=2
WALLTIME="0:30"
VPT="/sc/arion/projects/ji_lab/TOOLS/sif/vpt_1.3.sif"


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
#BSUB -J ${RUN_NAME}.${REGION}.MAKEVZG
#BSUB -P acc_ji_lab
#BSUB -q express
#BSUB -n ${NCORES}
#BSUB -R \"rusage[mem=${MEM}GB]\"
#BSUB -R \"span[hosts=1]\"
#BSUB -W ${WALLTIME}
#BSUB -Ep 'mailx -s \"LSF Job \${LSB_JOBNAME} \${LSB_JOBID} exited with status \${LSB_JOB_STATUS}\" paula.restrepo@icahn.mssm.edu'
#BSUB -o ${OUTPUT_DIR}/log/vzg_lsf.%J.log
#BSUB -eo ${OUTPUT_DIR}/log/vzg_lsf.%J.err
#BSUB -L /bin/bash

set -e # exit whole thing if one command fails

cd ${OUTPUT_DIR} 
ml singularity
export MALLOC_TRIM_THRESHOLD_=0
export DASK_TEMPORARY_DIRECTORY=${OUTPUT_DIR}/tmp

# make new VZG
singularity run ${VPT} vpt \\
--verbose --processes ${NCORES} \\
--log-file log/update_vzg.log \\
update-vzg \\
--input-vzg ${IN_DIR}/${REGION}/${RUN_NAME}_${REGION}.vzg \\
--input-boundaries baysor_micron_space.parquet \\
--input-entity-by-gene baysor_cell_by_gene.csv \\
--input-metadata baysor_cell_metadata.csv \\
--output-vzg "${OUTPUT_DIR}/${RUN_NAME}_${REGION}_baysor.vzg" \\
--temp-path tmp/ 

">${OUTPUT_DIR}/make_vzg.lsf
done

done < /sc/arion/projects/ji_lab/DATA/MERFISH/segmentation/NS-Atlas/baysor_0.7.1/run_names.txt