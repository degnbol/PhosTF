# run 5 samples for a given dataset (argument 1 to this bash script)
for i in {1..5}; do
../../../dream_test.jl run $1 $i > ${1}_${i}/training.log &
done
