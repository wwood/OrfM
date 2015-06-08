OrfM
====

A simple and not slow open reading frame (ORF) caller. No bells or whistles like frameshift detection, just a straightforward goal 
of returning a FASTA file of open reading frames over a certain length from a FASTA/Q file of nucleotide sequences. 

Install
----
Pre-compiled binaries are available at https://github.com/wwood/OrfM/releases

Or, you can install from source:
```sh
git clone --recursive https://github.com/wwood/OrfM
cd OrfM
make
```

Running
----
To find all reading frames greater than 96 nucleotides in length:
```sh
orfm <seq_file> >orfs.fa
```
The `<seq_file>` can be a FASTA or FASTQ file, gzipped or uncompressed. The default is 96
because this is the correct number for 100bp so that each of the 6 frames can be translated.
Using 99 would mean that the third frame forward (and the corresponding reverse frame) cannot 
possibly returned as an ORF because this would entail it encapsulating bases 2-101, and 101>100.

Output
---
The output ORFs fasta file contains any stretch of continuous codons which does not include a stop codon. 
There is no requirement for a start codon to be included in the ORF. One could say that OrfM is an ORF caller, not a gene caller (like say prodigal or genscan).

The output ORFs are named in a straitforward manner. The name of the sequence (i.e. anything before a space) is followed by `_startPosition_frameNumber_orfNumber` and then 
the comment of the sequence (i.e. anything after the space) is given after a space, if one exists. For example,
```
$ cat eg.fasta
>abc|123|name some comment
ATGTTA
$ orfm -m 3 eg.fasta
>abc|123|name_1_1_1 some comment
ML
```
The `startPosition` of reverse frames is the left-most position in the original sequence, not the codon where the ORF starts.

Not too slow
-----
It runs in reasonable time compared to e.g. `translate` from Sean Eddy's `squid` (available as part of the Ubuntu  [biosquid package](https://launchpad.net/ubuntu/+source/biosquid)), `getorf` from the `emboss` toolkit, and `prodigal`, a more nuanced gene caller. For a 463MB fasta file of 100bp sequences:
```
orfm -m 96 the.fa >orfm.fa
  #=> 7 seconds

translate -l 32 the.fa >biosquid.m33.txt
  #=> 29 seconds
  
getorf -sequence the.fa -minsize 96 -outseq getorf.fa
  #=> 38 sec

pigz -cd 110811_E_1_D_nesoni_single.fq.gz |fq2fa |prodigal -q -p meta -i /dev/stdin -a 110811_E_1_D_nesoni_single.prodigal.faa -o /dev/null
  #=> 16 min 6 sec
```
`translate` also does not appear to be able to handle fastq files (even piped in on `stdin` as fasta), and does not output a standard FASTA format file.

Contributing to OrfM
----
Patches most welcome. There is a few tests, which can be tested after installing `ruby`, as well as the `rspec` and `bio-commandeer` rubygems.
```sh
make test
```

Credits
----
Compiled into the code is `kseq.h` from [seqtk](https://github.com/lh3/seqtk) and an 
implementation of the [Alo-Corasick algorithm](https://en.wikipedia.org/wiki/Aho%E2%80%93Corasick_string_matching_algorithm)
from [strmat](http://web.cs.ucdavis.edu/~gusfield/strmat.html) modified [slightly](https://github.com/aurelian/ruby-ahocorasick).
Both are MIT licenced. A few GNU `libc` libraries are used too.

Software by Ben J. Woodcroft, currently unpublished. Released under LGPL - see the LICENSE.txt for licensing details.

