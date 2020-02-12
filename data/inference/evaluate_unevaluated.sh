for file in ??/WP*.mat; do dir=$(dirname $file); [ -f $dir/score.txt ] || echo $dir; done
