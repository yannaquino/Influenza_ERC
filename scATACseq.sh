#!/bin/bash
################################################################################
################################################################################
# File name: scATACseq.sh
# Author:
# Created: jeu. 30 avril 2020 11:41:51 CEST
# Last modified: mer. 06 mai 2020 11:56:56 CEST
################################################################################
################################################################################

################################################################################
# SLURM options
#SBATCH --job-name=scATAC2
#SBATCH --mail-user=yann.aquino@pasteur.fr
#SBATCH --mail-type=END
#SBATCH --error=/pasteur/homes/yaaquino/LOGS/scATACseqMonocytes/%j.err
#SBATCH --output=/pasteur/homes/yaaquino/LOGS/scATACseqMonocytes/%j.log
#SBATCH -c 8
#SBATCH --mem=100000
#SBATCH --partition=geh
#SBATCH --qos=geh
################################################################################

################################################################################
# Load modules
source /local/gensoft2/adm/etc/profile.d/modules.sh # module loading scripts
module purge # remove modules that may have been loaded by mistake
module load R/3.6.0
################################################################################

################################################################################
Rscript ~/SCRIPTS/Influenza_ERC/scATACseq/scATACseq.R
################################################################################
