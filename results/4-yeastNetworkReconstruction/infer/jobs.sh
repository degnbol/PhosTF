#!/usr/bin/env zsh

for params in \
'inner+-strict_norm+1.0+0.0++' \
'inner+-strict_norm+0.1+0.0++' \
'inner+-strict_norm+0.0+1.0++' \
'inner+-strict_norm+0.0+0.1++' \
'outer+-strict_norm+1.0+0.0++' \
'outer+-strict_norm+0.1+0.0++' \
'outer+-strict_norm+0.0+1.0++' \
'outer+-strict_norm+0.0+0.1++' \
'inner+-strict_norm+1.0+0.0+FDR20+' \
'inner+-strict_norm+0.1+0.0+FDR20+' \
'inner+-strict_norm+0.0+1.0+FDR20+' \
'inner+-strict_norm+0.0+0.1+FDR20+' \
'outer+-strict_norm+1.0+0.0+FDR20+' \
'outer+-strict_norm+0.1+0.0+FDR20+' \
'outer+-strict_norm+0.0+1.0+FDR20+' \
'outer+-strict_norm+0.0+0.1+FDR20+' \
'inner+-strict_norm+1.0+0.0++FDR20_abssign_mask' \
'inner+-strict_norm+0.1+0.0++FDR20_abssign_mask' \
'inner+-strict_norm+0.0+1.0++FDR20_abssign_mask' \
'inner+-strict_norm+0.0+0.1++FDR20_abssign_mask' \
'outer+-strict_norm+1.0+0.0++FDR20_abssign_mask' \
'outer+-strict_norm+0.1+0.0++FDR20_abssign_mask' \
'outer+-strict_norm+0.0+1.0++FDR20_abssign_mask' \
'outer+-strict_norm+0.0+0.1++FDR20_abssign_mask'
do
echo "#!/usr/bin/env bash
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=50gb
#SBATCH --time=0-20:00:00
#SBATCH --job-name=$params
#SBATCH --error=jobs/$params.err
#SBATCH --output=jobs/$params.out

# loadr
module load foss/2020b 2> /dev/null
module load r/4.1.0

time ./infer.sh '$params'
" > job-$params.slurm
sbatch job-$params.slurm
mv job-$params.slurm jobs/$params.slurm
done
