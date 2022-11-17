#!/bin/zsh

algos=("aoso2" "aoso4" "aoso6" "kmp" "bndm" "bndmq2" "bndmq4" "bsdm" "bsdm2" "bsdm3" "bsdm4" "sbndm" "sbndmq2" "sbndmq4" "bm" "ebom" "fbom" "hash3" "hash5" "hor");
for i in "${algos[@]}"
do
	echo $i
	bin/$i $1
done
