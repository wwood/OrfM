all:
	gcc -O2 -Wall orfm.c ext/ac.c -Iext -o orfm -lz
