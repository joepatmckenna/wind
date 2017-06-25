#!/bin/bash
#SBATCH -J "wind100"
#SBATCH -N 1
#SBATCH --ntasks-per-node=32
#SBATCH -p backfill
#SBATCH --mail-type="ALL"
#SBATCH -t 02:00:00
python ../write_flow.py 2005-08-26-12