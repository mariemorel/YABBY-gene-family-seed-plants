#!/bin/bash

#SBATCH --job-name=baliphy_AUplot_aa
#SBATCH --output=baliphy_AUplot_aa.out
#SBATCH --error=baliphy_AUplot_aa.err
#SBATCH --nodes=1                     # Number of nodes
#SBATCH --ntasks-per-node=1           # Number of tasks per node
#SBATCH --cpus-per-task=1             # Number of CPU cores per task
#SBATCH --time=3:00:00                # Maximum runtime (D-HH:MM:SS)

#gymno_codon
cut-range results_baliphy_angio_codon/pbil_angio_codon-1/C1.P1.fastas results_baliphy_angio_codon/pbil_angio_codon-2/C1.P1.fastas results_baliphy_angio_codon/pbil_angio_codon-3/C1.P1.fastas results_baliphy_angio_codon/pbil_angio_codon-4/C1.P1.fastas | alignment-chop-internal --tree results_baliphy_angio_codon/Results/greedy.PP.tree | alignment-max > results_baliphy_angio_codon/P1-max.fasta

cut-range results_baliphy_angio_codon/pbil_angio_codon-1/C1.P1.fastas results_baliphy_angio_codon/pbil_angio_codon-2/C1.P1.fastas results_baliphy_angio_codon/pbil_angio_codon-3/C1.P1.fastas results_baliphy_angio_codon/pbil_angio_codon-4/C1.P1.fastas | alignment-chop-internal --tree results_baliphy_angio_codon/Results/greedy.PP.tree | alignment-gild results_baliphy_angio_codon/P1-max.fasta results_baliphy_angio_codon/Results/greedy.PP.tree > results_baliphy_angio_codon/alignment-AU.prob

alignment-draw results_baliphy_angio_codon/P1-max.fasta --AU results_baliphy_angio_codon/alignment-AU.prob > results_baliphy_angio_codon/alignment-AU.html




