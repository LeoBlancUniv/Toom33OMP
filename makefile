CC = gcc
FLAG = -Wall -Wextra -pedantic -O3 -pthread -march=native -fopenmp
#FLAG = -Wall -Wextra -pedantic -O3 
LINKER = -lgmp

all:main.exe bench.exe

bench.exe:bench.o toom33_mul_mpn.o
	$(CC) $^ $(FLAG) -o $@ $(LINKER)

main.exe:main.o
	$(CC) $^ $(FLAG) -o $@ $(LINKER)

bench.o:bench.c
	$(CC) $< $(FLAG) -c -o $@ $(LINKER)

main.o:main.c
	$(CC) $< $(FLAG) -c -o $@ $(LINKER)

toom33_mul_mpn.o:toom33_mul_mpn.c
	$(CC) $< $(FLAG) -c -o $@ $(LINKER)