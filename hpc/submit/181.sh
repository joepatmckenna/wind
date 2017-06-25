#!/bin/bash
#SBATCH -J "wind181"
#SBATCH -N 1
#SBATCH --ntasks-per-node=32
#SBATCH -p backfill
#SBATCH --mail-type="ALL"
#SBATCH -t 02:00:00
python ../write_flow.py 2005-09-05-15