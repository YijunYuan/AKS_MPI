all:
	mpic++ main.cpp -Wall -O2 -I/usr/local/include -L/usr/local/lib  -Wl,-Bstatic -lflint -lmpfr -lmpir -lgmp -o output -Wl,-Bdynamic
