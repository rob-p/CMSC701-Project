CC=c++
CFLAGS=-I.

main: main.cpp
	$(CC) -std=c++17 -Wall -O3 -o main.out main.cpp -I ./ -I ./include/ -L ./ -lz

offsets: offsets.cpp
	g++ -std=c++17 -Wall -O3 -o offsets.out offsets.cpp -I ./ -I ./include/ -L ./ -lz

test_parser: test_parser.cpp
	g++ -std=c++17 -Wall -O3 -o test_parser.out test_parser.cpp -I ./ -I ./include/ -L ./ -lz -lpthread

baseline:
	g++ -std=c++17 -Wall -O3 -o countbases.out scripts/CountBases.cpp -I ./ -I ./include/ -L ./ -lz

fqfeeder:
	cd benchmarks && g++ -std=c++17 -Wall -O3 -o fqfeeder.out BenchmarkFQFeeder.cpp ./FQFeeder/src/FastxParser.cpp -I ./FQFeeder/include -L ./ -lz -lpthread && mv fqfeeder.out ..

all: main offsets

clean:
	rm -f zran.out offsets.out main.out countbases.out fqfeeder.out test_parser.out
