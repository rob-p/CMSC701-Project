CC=c++
CFLAGS=-I.

ffindex: main.cpp include/zran.hpp
	$(CC) -std=c++17 -ggdb -Wall -O3 -o ffindex main.cpp -I ./ -I ./include/ -L ./ -lz

offsets: offsets.cpp
	$(CC) -std=c++17 -Wall -O3 -o offsets.out offsets.cpp -I ./ -I ./include/ -L ./ -lz

test_parser: test_parser.cpp
	$(CC) -std=c++17 -Wall -O3 -o test_parser.out test_parser.cpp -I ./ -I ./include/ -L ./ -lz -lpthread

test_serial_parser: test_serial_parser.cpp
	$(CC) -std=c++17 -Wall -O3 -o test_serial_parser.out test_serial_parser.cpp -I ./ -I ./include/ -L ./ -lz -lpthread

baseline:
	$(CC) -std=c++17 -Wall -O3 -o countbases.out scripts/CountBases.cpp -I ./ -I ./include/ -L ./ -lz

fqfeeder:
	cd benchmarks && $(CC) -std=c++17 -Wall -O3 -o fqfeeder.out BenchmarkFQFeeder.cpp ./FQFeeder/src/FastxParser.cpp -I ./FQFeeder/include -L ./ -lz -lpthread && mv fqfeeder.out ..

all: main offsets

clean:
	rm -f zran.out offsets.out main.out countbases.out fqfeeder.out test_parser.out
