#!/bin/bash
################################################################################
################################################################################
# File name: demuxlet.sh
# Author:
# Created: mar. 28 avril 2020 09:48:38 CEST
# Last modified: jeu. 11 juin 2020 19:00:55 CEST
################################################################################
################################################################################

################################################################################
# SLURM options
#SBATCH --job-name=dmx1
#SBATCH --mail-user=yann.aquino@pasteur.fr
#SBATCH --mail-type=END
#SBATCH --error=/pasteur/homes/yaaquino/LOGS/scATACseqMonocytes/%j.err
#SBATCH --output=/pasteur/homes/yaaquino/LOGS/scATACseqMonocytes/%j.log
#SBATCH --ntasks=16
#SBATCH --mem=128000
#SBATCH --partition=common
#SBATCH --qos=normal
################################################################################

################################################################################
# Load modules
source /local/gensoft2/adm/etc/profile.d/modules.sh # module loading scripts
module purge # remove modules that may have been loaded by mistake
module load popscle/0.1-beta
# module load demuxet
################################################################################

################################################################################
SAMPLEDIR=/pasteur/projets/policy01/evo_immuno_pop/Mary/scATACseq_pilot
SAMPLE=ATAC2
TASKDIR=/pasteur/projets/policy01/evo_immuno_pop/Yann/demuxlet/internalValidity

popscle demuxlet --sam ${SAMPLEDIR}/cellRanger/${SAMPLE}/outs/possorted_bam.bam \
		 --tag-group CB --tag-UMI UB \
		 --vcf /pasteur/projets/policy01/evo_immuno_pop/Mary/scATACseq_pilot/demuxlet/popscle/ERC_continued_lifted38_8ind.vcf \
		 --field GT --cap-BQ 63 \
		 --out ${TASKDIR}/${SAMPLE}_test1

# demuxlet --sam ${SAMPLEDIR}/cellRanger/${SAMPLE}/outs/possorted_bam.bam \
# 		 --tag-group CB --tag-UMI UB \
# 		 --vcf /pasteur/projets/policy01/evo_immuno_pop/Mary/scATACseq_pilot/demuxlet/popscle/ERC_continued_lifted38_8ind.vcf \
# 		 --field GT --cap-BQ 63 \
# 		 --out ${TASKDIR}/${SAMPLE}_demuxletV2New
################################################################################
