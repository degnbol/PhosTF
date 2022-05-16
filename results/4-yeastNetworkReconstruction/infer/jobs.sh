#!/usr/bin/env zsh

for params in \
'inner+-strict+1.0+0.0++' \
'inner+-strict+0.1+0.0++' \
'inner+-strict+0.0+1.0++' \
'inner+-strict+0.0+0.1++' \
'outer+-strict+1.0+0.0++' \
'outer+-strict+0.1+0.0++' \
'outer+-strict+0.0+1.0++' \
'outer+-strict+0.0+0.1++'; do

echo "#!/usr/bin/env bash
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=0-20:00:00
#SBATCH --name=$params

# loadr
module load foss/2020b 2> /dev/null
module load r/4.1.0

time ./infer.sh '$params'
" > job-$params.slurm
sbatch job-$params.slurm

done
