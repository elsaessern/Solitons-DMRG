#!/bin/bash

# HOW TO USE:
# Put the onePoint.sh in the same folder as 
# OnePointFuncConsole.py and run the following command
# 		sbatch onePoint.sh 

# You can also run this without sbatch to run it on your 
# computer (using onePoint.sh $N $L $delta ${init conditions})

# EXAMPLE COMMAND:
# 	sbatch 



# SBATCH COMMANDS 
# name of the job
#SBATCH -J solitonJulia

#job resource specifications
#SBATCH -p share
#SBATCH --mem=8G
#SBATCH -c 4
#SBATCH --time=24:00:00

# files to write to
#SBATCH -o "solitonJuliaStdOut (job %j).txt"
#SBATCH -e "solitonJuliaStdError (job %j).txt"


# load python yay 
module load julia


# run the thing 
julia soliton.jl $@
