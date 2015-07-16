all:
	mpicc main.c -Wall -O2 -s -I/usr/local/include -L/usr/local/lib  -Wl,-Bstatic -lflint -lmpfr -lmpir -lgmp -lm -o output -Wl,-Bdynamic
