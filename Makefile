CC=gcc
CFLAGS=-I.

zran: zran.c
	$(CC) -o zran zran.c -lz

all: zran

clean:
	rm -f zran