#!/bin/bash
################################################################################
################################################################################
# File name: bigWigTracks.sh
# Author:
# Created: ven. 05 juin 2020 10:56:06 CEST
# Last modified: mar. 09 juin 2020 19:26:58 CEST
################################################################################
################################################################################

################################################################################
# SLURM options
#SBATCH --job-name=bigwigs
#SBATCH --mail-user=yann.aquino@pasteur.fr
#SBATCH --mail-type=END
#SBATCH --error=/pasteur/homes/yaaquino/LOGS/scATACseqMonocytes/%j.err
#SBATCH --output=/pasteur/homes/yaaquino/LOGS/scATACseqMonocytes/%j.log
#SBATCH --ntasks=16
#SBATCH --mem=100000
#SBATCH --qos=geh
#SBATCH --partition=geh
################################################################################

################################################################################
# Load modules
module load samtools
module load bedtools
module load BigWig_tools
module load UCSC-tools/v373
module load bx-python
module load deepTools/2.4.2 
################################################################################

################################################################################
individuals=(AFB108 AFB022 AFB176 AFB131 EUB058 EUB078 EUB060 EUB061)

cd /pasteur/projets/policy01/evo_immuno_pop/Mary/scATACseq_pilot/cellRanger/ATAC2/outs

for i in ${individuals[@]}

do

mkdir ${i}

echo "Created directory for ${i}"

sed -e s/^/CB:Z:/g ${i}Barcodes.tsv > ${i}/${i}BarcodesMod.tsv

#samtools view -H possorted_bam.bam > header.sam

awk -F "\t" -v OFS="\t" '{if (NR==FNR) {a[$1]} else {if ($19 in a) {print $0}}}' ${i}/${i}BarcodesMod.tsv <(samtools view possorted_bam.bam) | cat header.sam - > ${i}/${i}.bam

echo "Created BAM for ${i}"

cd ${i}

samtools sort -o ${i}Sorted.bam ${i}.bam

echo "Sorted ${i}.bam"

samtools index -b ${i}Sorted.bam

echo "Indexed ${i}.bam"

bamCoverage -b ${i}Sorted.bam --binSize 1 -o ${i}.bw

echo "Created ${i}.bw"

cd ..

done
################################################################################
