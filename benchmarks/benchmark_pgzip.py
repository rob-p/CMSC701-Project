import pgzip
import gzip
import time


def compress_fastq(filename):
    data = None
    size = None
    with gzip.open(filename, "rb") as file:
        data = file.read()
        data = data.decode('utf-8')
        size = len(data)
    
    with pgzip.open("test.txt.gz", "wt", thread=8, blocksize=2*10**8) as fw:
        fw.write(data)
    return size

def benchmark_pgzip(filename="test.txt.gz", size=1):
    t1 = time.perf_counter()
    A = 0
    C = 0
    G = 0
    T = 0
    with pgzip.open("test.txt.gz", "rt", thread=8) as fr:
        data = fr.read(size)
        for _x in data:
            if _x == "A":
                A = A + 1
            if _x == "C":
                C = C + 1
            if _x == "G":
                G = G + 1
            if _x == "T":
                T = T + 1
    print(f"A= {A}")
    print(f"C= {C}")
    print(f"G= {G}")
    print(f"T= {T}")

    t2 = time.perf_counter()
    print(f"took {t2-t1} ms")

if __name__ == "__main__":
    # size = compress_fastq("/nfshomes/sbharti/CMSC701-Project/data/nematode.fq.gzip")
    # print(size)//5174354312
    benchmark_pgzip("test.txt.gz", 5174354312)