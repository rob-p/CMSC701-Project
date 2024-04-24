CC=g++
CFLAGS=-I.

main: main.cpp
	$(CC) -std=c++17 -Wall -O3 -o main.out main.cpp -I ./ -I ./include/ -L ./ -lz

offsets: offsets.cpp
	g++ -std=c++17 -Wall -O3 -o offsets.out offsets.cpp -I ./ -I ./include/ -L ./ -lz

all: main offsets

clean:
	rm -f zran.out offsets.out main.out