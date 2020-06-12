#!/bin/bash
################################################################################
################################################################################
# File name: souporcell.sh
# Author:
# Created: mar. 28 avril 2020 09:03:29 CEST
# Last modified: ven. 12 juin 2020 08:44:12 CEST
################################################################################
################################################################################

################################################################################
# SLURM options
#SBATCH --job-name=spx1
#SBATCH --mail-user=yann.aquino@pasteur.fr
#SBATCH --mail-type=END
#SBATCH --error=/pasteur/homes/yaaquino/LOGS/scATACseqMonocytes/%j.err
#SBATCH --output=/pasteur/homes/yaaquino/LOGS/scATACseqMonocytes/%j.log
#SBATCH -c 8
#SBATCH --mem=100000
#SBATCH --partition=common
#SBATCH --qos=normal
################################################################################

################################################################################
# Load modules
source /local/gensoft2/adm/etc/profile.d/modules.sh # module loading scripts
module purge # remove modules that may have been loaded by mistake
module load souporcell/2.0
################################################################################

################################################################################
SAMPLE="ATAC2"
SAMPLEDIR="/pasteur/projets/policy01/evo_immuno_pop/Mary/scATACseq_pilot"
TASKDIR=/pasteur/projets/policy01/evo_immuno_pop/Yann/souporcell/internalValidity

cd ${TASKDIR}

# singularity exec \
#   ${TASKDIR}/souporcell/updatedVersion2/souporcell.sif \
#   souporcell_pipeline.py \

souporcell_pipeline.py \
  -i ${SAMPLEDIR}/cellRanger/${SAMPLE}/outs/possorted_bam.bam \
  -b ${SAMPLEDIR}/souporcell/passing/${SAMPLE}_custom_qc.tsv \
  -f ${SAMPLEDIR}/souporcell/genome.fa \
  -o ${TASKDIR}/${SAMPLE}_test1 \
  -k 8 \
  --skip_remap True \
  --no_umi True \
  --ignore True \
  -t ${SLURM_CPUS_PER_TASK}
################################################################################
