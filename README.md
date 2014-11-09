OrfM
====

NOTE: not really tested properly yet!

A simple and not slow ORF caller. No bells or whistles like frameshift detection, just a straightforward goal 
of returning a FASTA file of open reading frames over a certain length from a FASTA/Q file of nucleotide sequences. 

Install
----
```sh
make
```

Running
----
To find all reading frames greater than 99 nucleotides in length:
```sh
orfm -m 99 <seq_file> >orfs.fa
```
The `<seq_file>` can be a FASTA or FASTQ file, gzipped or uncompressed.

Not too slow
-----
It runs in reasonable time compared to e.g. `getorf` from the `emboss` toolkit. For a 300MB compressed fastq file:
```
orfm -m 33 the.fq.gz >orfs.fa
  #=> 42 seconds
  
pigz -cd the.fq.gz |fq2fa |getorf -sequence /dev/stdin -minsize 33 -outseq >orfs.fa
  #=> 3 min 17 sec
```
The vast majority of time is used by getorf rather than the first two commands in the pipe.

Credits
----
Compiled into the code is `kseq.h` from [seqtk](https://github.com/lh3/seqtk) and an 
implementation of the [Alo-Corasick algorithm](https://en.wikipedia.org/wiki/Aho%E2%80%93Corasick_string_matching_algorithm)
from [strmat](http://web.cs.ucdavis.edu/~gusfield/strmat.html) modified [slightly](https://github.com/aurelian/ruby-ahocorasick).
Both are MIT licenced. A few GNU `libc` libraries are used too.

Software by Ben J. Woodcroft, currently unpublished.

