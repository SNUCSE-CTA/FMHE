# Homomorphic Computations of Optimal and Approximate Local Alignments

A tool that computes optimal and approximate local alignments between two encrypted genomic sequences with affine gap penalty.

## Requirement

OpenMP<br/>
[TFHE](https://tfhe.github.io/)<br/>
[HEaan.STAT](https://www.cryptolab.co.kr/eng/product/heaan.php)<br/>
Python 3

## Run

```
python3 homla.py (app|opt) inputFileX inputFileY [reflen] [--scheme M,S,Go,Ge] [--lib (TFHE|HEAAN)]
```

`M`-Match score, `S`-Mismatch score, `Go`-Gap opening penalty, and `Ge`-Gap extending penalty.

## Default parameter

HE scheme: TFHE<br/>
Scoring scheme: 5/-3, -9, -1

## Sample input

For optimal local alignment: x.fasta, y.fasta<br/>
For approximate local alignment: x.vcf, y.vcf

## Sample run

Compute optimal local alignment:

```
python3 homla.py opt x.fasta y.fasta --scheme 5,-3,-9,-1
```

Output:

```
Scoring scheme: 5/-3, -9, -1
Input: x.fasta y.fasta
Lengths : 7 9
Encryption time: 1.01s
Homomorphically evaluating alignment...
Homomorphic computation time: 80.64s
Decryption time: 0.43s
RESULT
------------------
Score: 11
Starting pos: 2 3
Ending pos: 5 7
```

Compute approximate local alignment:

```
python3 homla.py app x.vcf y.vcf 26 --scheme 5,-3,-9,-1 --lib TFHE
```

Output:

```
Scoring scheme: 5/-3, -9, -1
Input: x.vcf y.vcf
Number of varaints: 5
Length of Reference Genome: 26
Encryption time: 1.04s
Homomorphically computing alignment...
Homomorphic computation time: 31.30s
Decryption time: 0.29s
RESULT
------------------
Score: 33
Starting pos: 5
Ending pos: 16
```

Or

```
python3 homla.py app x.vcf y.vcf 26 --scheme 5,-3,-9,-1 --lib HEAAN
```

Output:

```
Scoring scheme: 5/-3, -9, -1
Input: x.vcf y.vcf
Number of varaints: 5
Length of Reference Genome: 26
Encryption time: 3.54s
Homomorphically computing alignment...
Homomorphic computation time: 623.59s
Decryption time: 0.06s
RESULT
------------------
Score: 32.82
```

## References

[1] M. Bataa, S. Song, K. Park, M. Kim, J.H. Cheon, and S. Kim. <b>Homomorphic Computation of Local Alignment</b>, Proceedings of the IEEE BIBM 2020, pp. 2167-2174, 2020. <a href="https://doi.org/10.1109/BIBM49941.2020.9313199">DOI: 10.1109/BIBM49941.2020.9313199</a>
