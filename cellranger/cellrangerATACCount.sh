#!/bin/bash
################################################################################
################################################################################
# File name: cellrangerATACCount.sh
# Author:
# Created: mar. 28 avril 2020 09:27:52 CEST
# Last modified: mar. 28 avril 2020 09:36:56 CEST
################################################################################
################################################################################

################################################################################
# SLURM options
#SBATCH --job-name=ATAC2Count
#SBATCH --mail-user=yann.aquino@pasteur.fr
#SBATCH --mail-type=END
#SBATCH --error=/pasteur/homes/yaaquino/LOGS/scATACseqMonocytes/%j.err
#SBATCH --output=/pasteur/homes/yaaquino/LOGS/scATACseqMonocytes/%j.log
#SBATCH --ntasks=16
#SBATCH --mem=40000
#SBATCH --qos=geh
#SBATCH --partition=geh
################################################################################

################################################################################
# Load modules

################################################################################

################################################################################
SAMPLE="ATAC2"
TOOLDIR='/pasteur/projets/policy01/evo_immuno_pop/sc_resources'
source ${TOOLDIR}/cellranger-atac-1.2.0/sourceme.bash
FASTQPATH='/pasteur/projets/policy01/evo_immuno_pop/Mary/scATACseq_pilot/HN00121869_10X_RawData_Outs'

cd /pasteur/projets/policy01/evo_immuno_pop/Yann/monocytes/scATACseq/

cellranger-atac count --id=${SAMPLE} \
                      --fastqs=/pasteur/projets/policy01/evo_immuno_pop/Mary/scATACseq_pilot/HN00121869_10X_RawData_Outs/${SAMPLE}_100120/H73GYCCX2/ \
      		      --reference=/pasteur/projets/policy01/evo_immuno_pop/sc_resources/atac-reference/hsapiens_atac \
		      #--sample=${SAMPLE} \
		      --localcores=$SLURM_NTASKS \
		      --localmem=40 \
		      --expect-cells=10000
################################################################################
