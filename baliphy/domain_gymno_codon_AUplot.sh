#!/bin/bash

#SBATCH --job-name=baliphy_AUplot_aa
#SBATCH --output=baliphy_AUplot_aa.out
#SBATCH --error=baliphy_AUplot_aa.err
#SBATCH --nodes=1                     # Number of nodes
#SBATCH --ntasks-per-node=1           # Number of tasks per node
#SBATCH --cpus-per-task=1             # Number of CPU cores per task
#SBATCH --time=3:00:00                # Maximum runtime (D-HH:MM:SS)

#gymno_codon
cut-range angio_domain_codon/pbil_angio_domain_codon-1/C1.P1.fastas angio_domain_codon/pbil_angio_domain_codon-2/C1.P1.fastas angio_domain_codon/pbil_angio_domain_codon-3/C1.P1.fastas angio_domain_codon/pbil_angio_domain_codon-4/C1.P1.fastas | alignment-chop-internal --tree angio_domain_codon/Results/greedy.PP.tree | alignment-max > angio_domain_codon/P1-max.fasta

cut-range angio_domain_codon/pbil_angio_domain_codon-1/C1.P1.fastas angio_domain_codon/pbil_angio_domain_codon-2/C1.P1.fastas angio_domain_codon/pbil_angio_domain_codon-3/C1.P1.fastas angio_domain_codon/pbil_angio_domain_codon-4/C1.P1.fastas | alignment-chop-internal --tree angio_domain_codon/Results/greedy.PP.tree | alignment-gild angio_domain_codon/P1-max.fasta angio_domain_codon/Results/greedy.PP.tree > angio_domain_codon/alignment-AU.prob

alignment-draw angio_domain_codon/P1-max.fasta --AU angio_domain_codon/alignment-AU.prob > angio_domain_codon/alignment-AU.html

