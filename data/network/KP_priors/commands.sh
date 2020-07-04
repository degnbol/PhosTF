./KP2TF_wilcoxon_setup.R
./split_KPs.sh
# do in parallel:
./wilcoxon_part.sh KP_01.txt
./wilcoxon_part.sh ...

./after_parts.sh
./KP2TF.R
./WP.R
