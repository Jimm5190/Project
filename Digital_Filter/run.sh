#!/bin/bash
#!/usr/bin/env python3

gcc simple_filter.c -o simple_filter
M=(8 32 1024)
hL=("hL-8.txt" "hL-32.txt" "hL-1024.txt")
hR=("hR-8.txt" "hR-32.txt" "hR-1024.txt")
YL=("YL-8.txt" "YL-32.txt" "YL-1024.txt")
YR=("YR-8.txt" "YR-32.txt" "YR-1024.txt")
input="blue_giant_fragment.wav"
output=("output-8.wav" "output-32.wav" "output-1024.wav")

for i in {0..2};do
	./simple_filter ${M[i]} ${hL[i]} ${hR[i]} ${YL[i]} ${YR[i]} ${input} ${output[i]}
	done

python3 show_data.py 	    
