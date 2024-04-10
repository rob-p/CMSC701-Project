# Benchmarks to run
1. FQFeeder (https://github.com/rob-p/FQFeeder) - This is a lightweight class that
    provided single and paired-end Fasta/q parsing when there are multiple threads that
    parallelize the work done over the reads. 

**Task to be done is count number of bases**

# Datasets
## [seqkit-benchmark-data.tar.gz](http://app.shenwei.me/data/seqkit/seqkit-benchmark-data.tar.gz) - [details](https://bioinf.shenwei.me/seqkit/benchmark/)
### dataset_A.fa.gz
Size on disk = 867.2MB
Number of records = 67748

| Base | # of Base       |
|------|-----------------|
| A    |  761422333|
| C    | 634445262 | 
| G    | 634920403 |
| T    | 763678013 |

| Method                       | #Reader Threads | #Consumer Threads | Time to read (ms) | Total Time (ms) |
|------------------------------|-----------------|-------------------|-------------------|-----------------|
| Sequencial Read using kseq++ | NA              | NA                | 8607              | 25860              |
| FQFeeder                     | 1               | 1                 | NA                | 18627              |
| FQFeeder                     | 1               | 2                 | NA                | 9891              |
| FQFeeder                     | 1               | 3                 | NA                | 8852              |
| FQFeeder                     | 1               | 4                 | NA                | 8771              |
| FQFeeder                     | 1               | 5                 | NA                | 8828              |
| FQFeeder                     | 1               | 6                 | NA                | 8832              |
### dataset_B.fa.gz
Size on disk = 941.7MB
Number of records = 194

| Base | # of Base       |
|------|-----------------|
| A    |  409086287|
| C    | 277979602 | 
| G    | 278285676 |
| T    | 409609143 |

| Method                       | #Reader Threads | #Consumer Threads | Time to read (ms) | Total Time (ms) |
|------------------------------|-----------------|-------------------|-------------------|-----------------|
| Sequencial Read using kseq++ | NA              | NA                | 9351                | 19111           |
| FQFeeder                     | 1               | 1                 | NA                | 20715              |
| FQFeeder                     | 1               | 2                 | NA                | 20376              |
| FQFeeder                     | 1               | 3                 | NA                | 20559              |
| FQFeeder                     | 1               | 4                 | NA                | 20543              |
| FQFeeder                     | 1               | 5                 | NA                | 20453              |
| FQFeeder                     | 1               | 6                 | NA                | 20674              |

### dataset_C.fq.gz
Size on disk = 544.9MB
Number of records = 9186045

| Base | # of Base       |
|------|-----------------|
| A    |  235421932|
| C    | 223263521 | 
| G    | 223597960 |
| T    | 236193960 |

| Method                       | #Reader Threads | #Consumer Threads | Time to read (ms) | Total Time (ms) |
|------------------------------|-----------------|-------------------|-------------------|-----------------|
| Sequencial Read using kseq++ | NA              | NA                | 6499                | 11992              |
| FQFeeder                     | 1               | 1                 | NA                | 6439              |
| FQFeeder                     | 1               | 2                 | NA                | 6144              |
| FQFeeder                     | 1               | 3                 | NA                | 6138              |
| FQFeeder                     | 1               | 4                 | NA                | 6111              |
| FQFeeder                     | 1               | 5                 | NA                | 6120              |
| FQFeeder                     | 1               | 6                 | NA                | 6137              |


## Salmonella.fa.gz
Size on disk = 1.5MB
Number of records = 1

| Base | # of Base       |
|------|-----------------|
| A    |  1160904|
| C    | 1268422 | 
| G    | 1268221 |
| T    | 1159903 |

| Method                       | #Reader Threads | #Consumer Threads | Time to read (ms) | Total Time (ms) |
|------------------------------|-----------------|-------------------|-------------------|-----------------|
| Sequencial Read using kseq++ | NA              | NA                | 6                 | 50              |
| FQFeeder                     | 1               | 1                 | NA                | 61              |
| FQFeeder                     | 1               | 2                 | NA                | 62              |
| FQFeeder                     | 1               | 3                 | NA                | 63              |
| FQFeeder                     | 1               | 4                 | NA                | 63              |
| FQFeeder                     | 1               | 5                 | NA                | 60              |
| FQFeeder                     | 1               | 6                 | NA                | 60              |


## Ecoli.fq.gz
Size on disk = 41.3MB
Number of records = 1495900

| Base | # of Base       |
|------|-----------------|
| A    |  90216441|
| C    | 92370117 | 
| G    | 92760925 |
| T    | 90441732 |

| Method                       | #Reader Threads | #Consumer Threads | Time to read (ms) | Total Time (ms) |
|------------------------------|-----------------|-------------------|-------------------|-----------------|
| Sequencial Read using kseq++ | NA              | NA                | 1147                | 1991              |
| FQFeeder                     | 1               | 1                 | NA                | 1000              |
| FQFeeder                     | 1               | 2                 | NA                | 1004              |
| FQFeeder                     | 1               | 3                 | NA                | 1016              |
| FQFeeder                     | 1               | 4                 | NA                | 998              |
| FQFeeder                     | 1               | 5                 | NA                | 1002              |
| FQFeeder                     | 1               | 6                 | NA                | 1002              |
## Nematode
Size on disk = 686.9MB
Number of records = 11535642

| Base | # of Base       |
|------|-----------------|
| A    |  474118508|
| C    | 298365660 | 
| G    | 394364281 |
| T    | 563453165 |

| Method                       | #Reader Threads | #Consumer Threads | Time to read (ms) | Total Time (ms) |
|------------------------------|-----------------|-------------------|-------------------|-----------------|
| Sequencial Read using kseq++ | NA              | NA                | 9443                | 19311              |
| FQFeeder                     | 1               | 1                 | NA                | 10821              |
| FQFeeder                     | 1               | 2                 | NA                | 8533              |
| FQFeeder                     | 1               | 3                 | NA                | 8446              |
| FQFeeder                     | 1               | 4                 | NA                | 8319              |
| FQFeeder                     | 1               | 5                 | NA                | 8331              |
| FQFeeder                     | 1               | 6                 | NA                | 8371              |

# Commands to compile

1. FQFeeder  
```unix
cd $PRJECT_ROOT/benchmarks
g++ -std=c++17 -Wall -O3 -o fq_feeder BenchmarkFQFeeder.cpp ./FQFeeder/src/FastxParser.cpp -I ./FQFeeder/include -L ./ -lz
```

2. ../Scripts/CountBases.cpp
```unix
cd $PROJECT_ROOT
g++ -std=c++17 -Wall -O3 -o countbases scripts/CountBases.cpp -I ./ -I ./include/ -L ./ -lz
```