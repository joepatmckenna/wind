#!/bin/bash
files=submit/*.sh
for f in $files
do
	  (sbatch $f)
done
