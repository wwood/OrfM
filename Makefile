all: orfm

orfm: orfm.c
	gcc -O2 -Wall orfm.c ext/ac.c -Iext -o orfm -lz

test: orfm test/orfm_spec.rb
	rspec test/orfm_spec.rb

profile:
	gcc -O2 -Wall orfm.c ext/ac.c -Iext -o orfm_gprof -lz -pg -g -fprofile-arcs -ftest-coverage
