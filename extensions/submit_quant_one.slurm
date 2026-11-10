#!/bin/bash
#SBATCH -J q1
#SBATCH -p AMD_9A14,AMD_9654
#SBATCH -c 12
#SBATCH -t 4:00:00
bash /public/home/xuyuchen/data/RNA-seq_public/01.public1/scripts/quant_one.sh "$1" "$2"
