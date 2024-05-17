import pgzip
import gzip
import time
import concurrent.futures


def count_bases(segment):
    count = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
    for char in segment:
        if char in count:
            count[char] += 1
    return count

def combine_counts(total_count, segment_count):
    for key in total_count:
        total_count[key] += segment_count[key]


def threaded_count_dna_bases(dna_string, num_threads=8):
    # Determine the length of each segment
    segment_length = len(dna_string) // num_threads
    segments = [dna_string[i * segment_length:(i + 1) * segment_length] for i in range(num_threads)]
    # Handle any remainder by adding it to the last segment
    if len(dna_string) % num_threads:
        segments[-1] += dna_string[num_threads * segment_length:]
    
    total_count = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
    with concurrent.futures.ThreadPoolExecutor(max_workers=num_threads) as executor:
        results = executor.map(count_bases, segments)
        for result in results:
            combine_counts(total_count, result)
    
    return total_count

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
    
    total_count = threaded_count_dna_bases(data, num_threads=8)
    # for _x in data:
    #     if _x == "A":
    #         A = A + 1
    #     if _x == "C":
    #         C = C + 1
    #     if _x == "G":
    #         G = G + 1
    #     if _x == "T":
    #         T = T + 1
    print(f"A= {total_count['A']}")
    print(f"C= {total_count['C']}")
    print(f"G= {total_count['G']}")
    print(f"T= {total_count['T']}")
    """	# of Base
A	474118508
C	298365660
G	394364281
T	563453165
"""

    t2 = time.perf_counter()
    print(f"took {t2-t1} secs")

if __name__ == "__main__":
    # size = compress_fastq("/nfshomes/sbharti/CMSC701-Project/data/nematode.fq.gzip")
    # print(size)//5174354312
    benchmark_pgzip("test.txt.gz", 5174354312)