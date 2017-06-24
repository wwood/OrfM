AUTOMAKE_OPTIONS=foreign

CFLAGS= -O2 -Wall -D ORFM_VERSION=\"$(VERSION)\" -Iext
LDFLAGS= -lz

bin_PROGRAMS=orfm
orfm_SOURCES= src/orfm.c ext/ac.c ext/ac.h ext/seqtk/kseq.h test/orfm_spec.rb test/20_random100bp.fna

DEF := -D ORFM_VERSION=$(ORFM_VERSION)

debug: clean
	gcc $(CFLAGS) -g src/orfm.c ext/ac.c -Iext -o orfm $(LDFLAGS)

check: orfm test/orfm_spec.rb
	rspec test/orfm_spec.rb

profile:
	gcc $(CFLAGS) src/orfm.c ext/ac.c -Iext -o orfm_gprof $(LDFLAGS) -pg -g -fprofile-arcs -ftest-coverage

static_linux:
	mkdir -p orfm-$(VERSION)_Linux_x86_64
	gcc $(CFLAGS) --static src/orfm.c ext/ac.c -Iext -o orfm-$(VERSION)_Linux_x86_64/orfm-$(VERSION)_Linux_x86_64 $(LDFLAGS)
	cd orfm-$(VERSION)_Linux_x86_64 && ln -s orfm-$(VERSION)_Linux_x86_64 orfm && cd -
	tar czf orfm-$(VERSION)_Linux_x86_64.tar.gz orfm-$(VERSION)_Linux_x86_64
