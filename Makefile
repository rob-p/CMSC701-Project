CC=g++
CFLAGS=-I.

zran: zran.cpp
	$(CC) -o zran.out zran.cpp -lz

all: zran

clean:
	rm -f zran.out