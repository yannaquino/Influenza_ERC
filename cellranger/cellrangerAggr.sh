#!/bin/bash
################################################################################
################################################################################
# File name: cellrangerAggr.sh
# Author:
# Created: mar. 28 avril 2020 09:20:44 CEST
# Last modified: mar. 28 avril 2020 09:28:02 CEST
################################################################################
################################################################################

################################################################################
# SLURM options
#SBATCH --job-name=ATACAggr
#SBATCH --mail-user=yann.aquino@pasteur.fr
#SBATCH --mail-type=END
#SBATCH --error=/pasteur/homes/yaaquino/LOGS/scRNAseqPBMCs/%j.err
#SBATCH --output=/pasteur/homes/yaaquino/LOGS/scRNAseqPBMCs/%j.log
#SBATCH --ntasks=16
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

cd /pasteur/projets/policy01/evo_immuno_pop/Yann/monocytes/scATACseq/

cellranger-atac aggr --id=ATACAggrCustom \
                     --csv=aggr.csv \
                     --reference=/pasteur/projets/policy01/evo_immuno_pop/sc_resources/atac-reference/hsapiens_atac \
                     --normalize=depth \
                     --dim-reduce=lsa
################################################################################
