all: orfm

orfm: orfm.c
	gcc -O2 -Wall orfm.c ext/ac.c -Iext -o orfm -lz

test: orfm test/orfm_spec.rb
	rspec test/orfm_spec.rb
