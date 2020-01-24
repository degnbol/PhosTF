while read KP; do
	echo $KP
	./KP2TF_wilcoxon.R $KP
done < $1
