CFLAGS := -O2 -Wall
LIBRARIES := -lz
VERSION := 0.2.0

all: orfm

orfm: orfm.c
	gcc $(CFLAGS) orfm.c ext/ac.c -Iext -o orfm $(LIBRARIES)

test: orfm test/orfm_spec.rb
	rspec test/orfm_spec.rb

profile:
	gcc $(CFLAGS) orfm.c ext/ac.c -Iext -o orfm_gprof $(LIBRARIES) -pg -g -fprofile-arcs -ftest-coverage

static_linux:
	mkdir -p orfm-$(VERSION)_Linux_x86_64
	gcc $(CFLAGS) --static orfm.c ext/ac.c -Iext -o orfm-$(VERSION)_Linux_x86_64/orfm-$(VERSION)_Linux_x86_64 $(LIBRARIES)
	cd orfm-$(VERSION)_Linux_x86_64 && ln -s orfm-$(VERSION)_Linux_x86_64 orfm && cd -
	tar czf orfm-$(VERSION)_Linux_x86_64.tar.gz orfm-$(VERSION)_Linux_x86_64
