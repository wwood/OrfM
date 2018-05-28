OrfM
====

A simple and not slow open reading frame (ORF) caller. No bells or whistles like frameshift detection, just a straightforward goal 
of returning a FASTA file of open reading frames over a certain length from a FASTA/Q file of nucleotide sequences. 

Install
----
OrfM can be installed in 3 ways.
### 1) Install from pre-compiled binaries
OrfM can be installed by downloading pre-compiled binaries available at https://github.com/wwood/OrfM/releases. Once you have downloaded the package, extract and run it e.g. for GNU/Linux:
```sh
tar xzf orfm-x.x.x_Linux_x86_64.tar.gz
cd orfm-x.x.x_Linux_x86_64
./orfm -h
```
### 2) Install from source
If you desire, OrfM can also be installed from source. Download the `orfm-x.x.x.tar.gz` from the [releases](https://github.com/wwood/OrfM/releases) page (_not_ the 'Source code' or the 'Download zip') and then follow the usual protocol for compilation and installation:
```sh
tar xzf orfm-x.x.x.tar.gz
cd orfm-x.x.x
./configure
make
```
To run `make check` you need Ruby and as well as the `rspec` and `bio-commandeer` rubygems. This step is optional.
```
gem install rspec bio-commandeer # may require 'sudo'
make check
```
Then finally to install OrfM
```
sudo make install
orfm -h
```
### 3) Install with GNU Guix
Or, you can install through [guix](http://www.gnu.org/software/guix/):
```
guix package -i orfm
```

### 4) Install with brew
Thanks to Torsten Seemann (@tseemann), OrfM can be installed through homebrew:
```
brew install brewsci/bio/orfm
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
----
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

FAQ
----
### `bash: ./configure: No such file or directory`

This can happen when trying to build OrfM from source. It might mean that the original source code has been downloaded, rather than the 'dist' archive. Download `orfm-x.x.x.tar.gz` from the [releases page](https://github.com/wwood/OrfM/releases) which contains the `configure` script (not the 'Source code'), and then follow the instructions for building from source above.

Contributing to OrfM
----
Patches most welcome. To get started:
```sh
git clone --recursive https://github.com/wwood/OrfM
cd OrfM
./autogen.sh
./configure
make check
```

Credits
----
Compiled into the code is `kseq.h` from [seqtk](https://github.com/lh3/seqtk) and an 
implementation of the [Alo-Corasick algorithm](https://en.wikipedia.org/wiki/Aho%E2%80%93Corasick_string_matching_algorithm)
from [strmat](http://web.cs.ucdavis.edu/~gusfield/strmat.html) modified [slightly](https://github.com/aurelian/ruby-ahocorasick).
Both are MIT licenced. A few GNU `libc` libraries are used too.

Citing OrfM
----
Software (c) Ben J. Woodcroft, released under LGPL - see the LICENSE.txt for licensing details.

A peer-reviewed manuscript describing OrfM has been published. If you use OrfM in your work then please help us out by citing it - thank you.

Ben J. Woodcroft, Joel A. Boyd, and Gene W. Tyson. [_OrfM: A fast open reading frame predictor for metagenomic data_](http://bioinformatics.oxfordjournals.org/content/32/17/2702). (2016). Bioinformatics. doi:10.1093/bioinformatics/btw241.

