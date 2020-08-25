#!/bin/bash
#SBATCH -t 03:00:00 
#SBATCH -p normal 
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=m.galland@uva.nl
#SBATCH --cpus-per-task=20
#SBATCH --time=10:00:00
#SBATCH --mem=128G
#SBATCH --job-name=snp_calling


# If specified, has to match the cpus-per-task required
# --nodes=1

mkdir -p "$TMPDIR"/marc_scratch


# build the singularity image first
# singularity build freebayes.simg docker://bleekerlab/freebayes:latest 
singularity run freebayes.simg --cores 30 \
  --config workdir="$TMPDIR"/marc_scratch/ \
           resultdir="$HOME"/workspace/freebayes_snp_calling/results/ 