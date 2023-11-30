## Kallisto to matrix

<img data-toggle="modal" data-target="[data-modal='10.5281/zenodo.8285146']" src="https://zenodo.org/badge/DOI/10.5281/zenodo.8285146.svg" alt="10.5281/zenodo.8285146">

Kaiquan Yu. (2023). Kallisto2matrix: a simple script to convert kallisto/salmon output to matrix files (0.0.1). Zenodo. https://doi.org/10.5281/zenodo.8285146

🙈: yukaiquan
<br/>
📧: 1962568272@qq.com

### 1. Introduction

This is a simple script to convert kallisto/salmon output to matrix file. The input file is a list of kallisto/salmon output directory. The output file is a matrix file. The matrix file is a tab-delimited file. The first column is gene id, the first row is sample name, and the other cells are TPM values.

### 2. Installation (no requirement)

```bash
# windowse or linux
cargo build --release
# Win10 cross-compile to linux application
cargo build --release --target=x86_64-unknown-linux-musl
```

you can find the executable file in target/release/kallisto2matrix or target/release/kallisto2matrix.exe

### 3. Usage

**Now, you can use the executable file to convert kallisto/salmon output to matrix file.**

````
## linux
kallisto2matrix -i samples.txt -o test
## windowse
kallisto2matrix.exe -i samples.txt -o test

```

samples.txt is a list of kallisto output directory. The format is as follows:

```

./example/BB313-01T0001_sfs/abundance.tsv,01T0001_sfs
./example/BB313-01T0002_sfs/abundance.tsv,01T0002_sfs

```

The first column is kallisto/salmon output file, the second column is sample name. The output file is test_tpm.matrix and test_count.matrix.txt.

test_tpm.matrix.txt:

```

gene_id 01T0001_sfs 01T0002_sfs
A.satnudSFS6C01G000491.1 0 0
A.satnudSFS3D01G003755.1 6.44158 9.49274
A.satnudSFS6A01G001901.1 61.5822 90.448
A.satnudSFS1A01G000402.1 2.57207 2.24161
A.satnudSFS7C01G002031.1 0 0
A.satnudSFS6D01G000164.1 0 0

````

### 4. test

```bash
cd example
## linux
../kallisto2matrix -i samples.txt -o test
## windowse
../kallisto2matrix.exe -i samples.txt -o test

Welcome to use kallisto2matrix! Author: Yu kaiquan <1962568272@qq.com>
read kallisto output file: ./BB313-01T0001_sfs/abundance.tsv
read kallisto output file: ./BB313-01T0002_sfs/abundance.tsv
write count matrix file: test_count_matrix.txt
write count matrix file: test_count_matrix.txt done!
write tpm matrix file: test_tpm_matrix.txt
write tpm matrix file: test_tpm_matrix.txt done!
Done!Goodbye!
Total elapsed time: 582.4223ms

# convert salmon output to matrix
## windows
..\kallisto2matrix.exe -i .\samples_salmon.txt -o salmon -t salmon
## linux
..\kallisto2matrix -i .\samples_salmon.txt -o salmon -t salmon

```

### 5. threads

I want to use 4 threads to run kallisto2matrix.

#### 5.1 threads

```bash
../kallisto2matrix.exe -i ./samples_salmon.txt -o salmon -t salmon
Welcome to use kallisto2matrix! Author: Yu kaiquan <1962568272@qq.com>
read salmon output file: ./SRR22937051_sfs/quant.sf
read salmon output file: ./SRR22937052_sfs/quant.sf
write count matrix file: salmon_count_matrix.txt
write count matrix file: salmon_count_matrix.txt done!
write tpm matrix file: salmon_tpm_matrix.txt
write tpm matrix file: salmon_tpm_matrix.txt done!
write fpkm matrix file: salmon_fpkm_matrix.txt
write fpkm matrix file: salmon_fpkm_matrix.txt done!
Done!Goodbye!
Total elapsed time: 1.0790253s
```

#### 5.2 4 threads

```bash
../kallisto2matrix.exe -i ./samples_salmon.txt -o salmon -t salmon
Welcome to use kallisto2matrix! Author: Yu kaiquan <1962568272@qq.com>
read salmon output file: ./SRR22937051_sfs/quant.sf
read salmon output file: ./SRR22937052_sfs/quant.sf
write count matrix file: salmon_count_matrix.txt
write tpm matrix file: salmon_tpm_matrix.txt
write fpkm matrix file: salmon_fpkm_matrix.txt
write count matrix file: salmon_count_matrix.txt done!
write tpm matrix file: salmon_tpm_matrix.txt done!
write fpkm matrix file: salmon_fpkm_matrix.txt done!
Done!Goodbye!
Total elapsed time: 805.5274ms
```

It is much faster than 1 thread.(800ms vs 1s)
