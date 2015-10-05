AUTOMAKE_OPTIONS=foreign

CFLAGS= -O2 -Wall -D ORFM_VERSION=\"$(VERSION)\" -Iext
LDFLAGS= -lz

bin_PROGRAMS=orfm
orfm_SOURCES= orfm.c ext/ac.c ext/ac.h ext/seqtk/kseq.h

DEF := -D ORFM_VERSION=$(ORFM_VERSION)

debug: clean
	gcc $(CFLAGS) -g orfm.c ext/ac.c -Iext -o orfm $(LDFLAGS)

test: orfm test/orfm_spec.rb
	rspec test/orfm_spec.rb

profile:
	gcc $(CFLAGS) orfm.c ext/ac.c -Iext -o orfm_gprof $(LDFLAGS) -pg -g -fprofile-arcs -ftest-coverage

static_linux:
	mkdir -p orfm-$(VERSION)_Linux_x86_64
	gcc $(CFLAGS) --static orfm.c ext/ac.c -Iext -o orfm-$(VERSION)_Linux_x86_64/orfm-$(VERSION)_Linux_x86_64 $(LDFLAGS)
	cd orfm-$(VERSION)_Linux_x86_64 && ln -s orfm-$(VERSION)_Linux_x86_64 orfm && cd -
	tar czf orfm-$(VERSION)_Linux_x86_64.tar.gz orfm-$(VERSION)_Linux_x86_64
