head -n20 KP.txt > KP_01.txt
head -n37 KP.txt | tail -n17 > KP_02.txt
head -n54 KP.txt | tail -n17 > KP_03.txt
head -n71 KP.txt | tail -n17 > KP_04.txt
head -n88 KP.txt | tail -n17 > KP_05.txt
head -n105 KP.txt | tail -n17 > KP_06.txt
head -n122 KP.txt | tail -n17 > KP_07.txt
head -n139 KP.txt | tail -n17 > KP_08.txt
head -n156 KP.txt | tail -n17 > KP_09.txt
tail -n17 KP.txt > KP_10.txt
# was split correctly
cat KP_??.txt | cmp - KP.txt
