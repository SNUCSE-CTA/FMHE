# Compressed/Searchable Index and Genome Analysis with Homomorphic Encryption

# (1) FM-index of Alignment

Stored in `FMA` directory.

## Requirements

[Succinct Data Structures Library 2.0](https://github.com/simongog/sdsl-lite) (SDSL)

## Installation

To complie, set `SDSL_DIR` in the makefile to point to your SDSL directory.

```
git clone https://github.com/SNUCSE-CTA/FMHE
cd FMHE/FMA
make all
```

## Running example in [4]

**Build FM-index of alignment**  
Set `saa_dens` and `isaa_dens` in buildFMA.cpp, locate.cpp, and extract.cpp same as `sampling_rate` of SAA. (default: 32)

```
mkdir example/output

#./buildSAA sampling_rate reference_file vcf_file output_directory [num_threads]
./buildSAA 32 example/example_reference example/example.vcf ./example/output/ 8

#./buildFMA saa_directory num_strings output_directory
./buildFMA ./example/output/ 4 ./example/output/
```

**Locate**

```
#./locate saa_directory fma_directory num_strings pattern_file
./locate ./example/output/ ./example/output/ 4 ./example/locate_example.txt > ./example/output/result_locate.txt
```

**Extract**

```
#./extract saa_directory fma_directory num_strings pattern_file pattern_size
./extract ./example/output/ ./example/output/ 4 ./example/extract_example.txt 4 > ./example/output/result_extract.txt
```

## Data

[chr1.tar.gz](https://drive.google.com/open?id=1fX4Re8hhidHLYtrYovYGyi9UXVVAobhD) (61MB)  
[3.chr1.vcf.tar.gz](https://drive.google.com/open?id=1iqNvOYUndey_PatDBIf8hA3kyMBSfk1S) (20MB)  
[100.chr1.vcf.tar.gz](https://drive.google.com/open?id=1bCgA8bQVf2sPPgapQe7ZzAlXxWxrZ2cG) (495MB)  
[250.chr1.vcf.tar.gz](https://drive.google.com/open?id=1Xk3ioGgcKp1bDBU1cAujFdee-nlNyG5j) (1.5GB)  
[500.chr1.vcf.tar.gz](https://drive.google.com/open?id=1T27ahsMWGytraDqpzhjNGWftl8TJaokX) (3.8GB)  
[locate pattern](https://drive.google.com/open?id=1eobVHrqMCAX6J5gPfawNg8au3EdiLgsm)  
[extract pattern](https://drive.google.com/open?id=1nJ5tODAB--S0IJE73BR5ZvVmSP-5pS4Q)

## References

[1] Joong Chae Na, Heejin Park, Maxime Crochemore, Jan Holub, Costas S. Iliopoulos, Laurent Mouchard, and Kunsoo Park. **Suffix tree of alignment: An efficient index for similar data**, In International Workshop on Combinatorial Algorithms, 8288:337-348, Springer, Berlin, Heidelberg, 2013. [DOI: 10.1007/978-3-642-45278-9_29](https://doi.org/10.1007/978-3-642-45278-9_29)

[2] Joong Chae Na, Heejin Park, Sunho Lee, Minsung Hong, Thierry Lecroq, Laurent Mouchard, and Kunsoo Park. **Suffix array of alignment: A practical index for similar data**, In International Symposium on String Processing and Information Retrieval, 8214:243-254, Springer, Cham, 2013. [DOI: 10.1007/978-3-319-02432-5_27](https://doi.org/10.1007/978-3-319-02432-5_27)

[3] Joong Chae Na, Hyunjoon Kim, Heejin Park, Thierry Lecroq, Martine Léonard, Laurent Mouchard, and Kunsoo Park. **FM-index of alignment: A compressed index for similar strings**, Theoretical Computer Science 638:159-170, 2016. [DOI: 10.1016/j.tcs.2015.08.008](https://doi.org/10.1016/j.tcs.2015.08.008)

[4] Joong Chae Na, Hyunjoon Kim, Seunghwan Min, Heejin Park, Thierry Lecroq, Martine Léonard, Laurent Mouchard, and Kunsoo Park. **FM-index of alignment with gaps**, Theoretical Computer Science 710:148-157, 2018. [DOI: 10.1016/j.tcs.2017.02.020](https://doi.org/10.1016/j.tcs.2017.02.020)

# (2) Homomorphic Computation of Local Alignment

A tool that computes optimal and approximate local alignments between two encrypted genomic sequences with affine gap penalty.</br>
Stored in `HLA` directory.

## Requirements

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
