#!/bin/bash
#SBATCH -p cpu
#SBATCH --job-name=LCS_MPI
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G              # cabe confortavelmente nos 15GB do nó
#SBATCH --time=01:00:00
#SBATCH --output=log_lcs_%j.out


mpirun -np $SLURM_NTASKS ./lcs_mpi testes/inputs/teste_10000A.in testes/inputs/teste_10000B.in

