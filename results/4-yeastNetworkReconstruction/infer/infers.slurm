#!/usr/bin/env bash
#SBATCH -N 1
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=4
#SBATCH --time=0-00:50:00

# loadr
module load foss/2020b 2> /dev/null
module load r/4.1.0

time srun -n 1 -c 4 ./infer.sh 'inner+-strict+1.0+0.0++'
time srun -n 1 -c 4 ./infer.sh 'inner+-strict+0.1+0.0++'
time srun -n 1 -c 4 ./infer.sh 'inner+-strict+0.0+1.0++'
time srun -n 1 -c 4 ./infer.sh 'inner+-strict+0.0+0.1++'
time srun -n 1 -c 4 ./infer.sh 'outer+-strict+1.0+0.0++'
time srun -n 1 -c 4 ./infer.sh 'outer+-strict+0.1+0.0++'
time srun -n 1 -c 4 ./infer.sh 'outer+-strict+0.0+1.0++'
time srun -n 1 -c 4 ./infer.sh 'outer+-strict+0.0+0.1++'
