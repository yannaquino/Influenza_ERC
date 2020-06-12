#!/bin/bash
################################################################################
################################################################################
# File name: mkref.sh
# Author:
# Created: mar. 28 avril 2020 09:37:53 CEST
# Last modified: mar. 28 avril 2020 09:41:55 CEST
################################################################################
################################################################################

################################################################################
# SLURM options
#SBATCH --job-name=mkrefATAC
#SBATCH --mail-user=yann.aquino@pasteur.fr
#SBATCH --mail-type=END
#SBATCH --error=/pasteur/homes/yaaquino/LOGS/scATACseqMonocytes/%j.err
#SBATCH --output=/pasteur/homes/yaaquino/LOGS/scATACseqMonocytes/%j.log
#SBATCH --ntasks=8
#SBATCH --mem=40000
#SBATCH --qos=geh
#SBATCH --partition=geh
################################################################################

################################################################################
# Load modules

################################################################################

################################################################################
TOOLDIR='/pasteur/projets/policy01/evo_immuno_pop/sc_resources'
source ${TOOLDIR}/cellranger-atac-1.2.0/sourceme.bash

cd ${TOOLDIR}/atac-reference

cellranger-atac mkref hsapiens_atac --config hsapiens_atac.config
################################################################################
