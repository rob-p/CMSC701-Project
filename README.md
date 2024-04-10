# CMSC701-Project
Class project for CMSC701

## Installing zlib
```
git clone git@github.com:madler/zlib.git
cd zlib

./configure
make test

# if the above succeeds
make install
```

## Generating test data
```
mkdir test-data
cd test-data

# 1.5 MB compressed
wget https://raw.githubusercontent.com/umd-cmsc701/project_0_test_data/main/test_input/salmonella.fa
gzip salmonella.fa

# 28 MB compressed
wget https://www.be-md.ncbi.nlm.nih.gov/Traces/sra-reads-be/fastq?acc=SRR28596297 -O ecoli.fq.gzip

# 595 MB compressed
wget https://www.be-md.ncbi.nlm.nih.gov/Traces/sra-reads-be/fastq?acc=SRR28592514 -O nematode.fq.gzip
```

## Running this project

1. Build everything
```
make
```

2. Generate index file
```
./zran.out test-data/salmonella.fa.gz
```

3. Use generated index file to extract `LEN` (16384) bytes from an offset:

```
./zran.out test-data/salmonella.fa.gz 400 ./test-data/salmonella.fa.gz.index
```