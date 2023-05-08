## Kallisto to matrix

🙈: yukaiquan
<br/>
📧: 1962568272@qq.com

### 1. Introduction

This is a simple script to convert kallisto output to matrix file. The input file is a list of kallisto output directory. The output file is a matrix file. The matrix file is a tab-delimited file. The first column is gene id, the first row is sample name, and the other cells are TPM values.

### 2. Usage

```
kallisto2matrix -i samples.txt -o test
```

samples.txt is a list of kallisto output directory. The format is as follows:

```
./example/BB313-01T0001_sfs/abundance.tsv,01T0001_sfs
./example/BB313-01T0002_sfs/abundance.tsv,01T0002_sfs
```

The first column is kallisto output file, the second column is sample name. The output file is test_tpm.matrix and test_count.matrix.txt.

test_tpm.matrix.txt:

```
gene_id	01T0001_sfs	01T0002_sfs
A.satnudSFS6C01G000491.1	0	0
A.satnudSFS3D01G003755.1	6.44158	9.49274
A.satnudSFS6A01G001901.1	61.5822	90.448
A.satnudSFS1A01G000402.1	2.57207	2.24161
A.satnudSFS7C01G002031.1	0	0
A.satnudSFS6D01G000164.1	0	0
```

### 3. test

```bash
cd example
../kallisto2matrix -i samples.txt -o test


Welcome to use kallisto2matrix! Author: Yu kaiquan <1962568272@qq.com>
read kallisto output file: ./BB313-01T0001_sfs/abundance.tsv
read kallisto output file: ./BB313-01T0002_sfs/abundance.tsv
write count matrix file: test_count_matrix.txt
write count matrix file: test_count_matrix.txt done!
write tpm matrix file: test_tpm_matrix.txt
write tpm matrix file: test_tpm_matrix.txt done!
Done!Goodbye!
Total elapsed time: 582.4223ms

```
