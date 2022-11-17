import os
import sys

#algos=["chct", "ct", "aoso2", "aoso4", "aoso6", "kmp", "bndm", "bndmq2", "bndmq4", "bndmq6", "bsdm", "bsdm2", "bsdm3", "bsdm4", "sbndm", "sbndm2", "sbndmq2", "sbndmq4", "bm", "ebom", "fbom", "hash3", "hash5", "hor", "simdv", "simdv2"];
#nums=["5", "9", "17", "33", "65"];
nums=["7", "13"]
#algos=["ct", "chct", "sbndm", "aoso4", "simdv", "simdv2"];
#algos=["simdv", "sksop"]
algos=["ct", "chct", "sbndmq2", "sbndmq4", "sbndmq6","shor4", "shor8", "shor12", "shor16", "ssk4", "ssk8", "ssk12", "ssk16"]
names=["KMPCT", "KMPCT2", "SBNDMCT2", "SBNDMCT4", "SBNDMCT6", "BMHCT4", "BMHCT8", "BMHCT12", "BMHCT16", "SKSCT4","SKSCT8", "SKSCT12", "SKSCT16"]
#algos=["shor4", "shor8", "shor12", "shor16", "ssk4", "ssk8", "ssk12", "ssk16"]
#names=["BMHCT4", "BMHCT8", "BMHCT12", "BMHCT16", "SKSCT4","SKSCT8", "SKSCT12", "SKSCT16"]
for k in range(1):
    for j in nums:
        print('&', end=' ')
        print(j, end=' ')
        sys.stdout.flush()
        for i in algos:
            print('&', end=" ")
            sys.stdout.flush()
            #cmd = "taskset 8 bin/"+i+" "+j+" > dummy"
            cmd = "taskset 8 bin/"+i+" "+j
            os.system(cmd)
            print(" ", end="")
        print("\\\\")
        sys.stdout.flush()
