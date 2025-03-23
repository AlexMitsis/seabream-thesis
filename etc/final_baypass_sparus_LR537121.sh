#!/bin/bash
#SBATCH --job-name=final_LR537121
#SBATCH --time=36:00:00               # Run time (days-hh:mm:ss) - (max 7days)
#SBATCH --ntasks-per-node=12
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=amitsb@bio.auth.gr

# Run command directly with file paths
./g_baypass -gfile ../Sparus_WG_baypass_per_chr/Final_Test_Alex/maf_LR537121.genotype -efile ../Sparus_WG_baypass_per_chr/sparus_wg_maf.covariate -contrastfile ../Sparus_WG_baypass_per_chr/sparus_wg_maf.contrast -poolsizefile ../Sparus_WG_baypass_per_chr/sparus_wg_maf.poolsize -outprefix Final_Sparus_LR537121_baypass -nthreads 12 -nval 10000 -burnin 10000
