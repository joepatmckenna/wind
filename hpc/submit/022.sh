#!/bin/bash
#SBATCH -J "wind022"
#SBATCH -N 1
#SBATCH --ntasks-per-node=32
#SBATCH -p backfill
#SBATCH --mail-type="ALL"
#SBATCH -t 02:00:00
python ../write_flow.py 2005-08-16-18