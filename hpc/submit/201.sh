#!/bin/bash
#SBATCH -J "wind201"
#SBATCH -N 1
#SBATCH --ntasks-per-node=32
#SBATCH -p backfill
#SBATCH --mail-type="ALL"
#SBATCH -t 02:00:00
python ../write_flow.py 2005-09-08-03