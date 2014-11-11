#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <getopt.h>    /* for getopt */

#include <zlib.h>
#include <seqtk/kseq.h>

#include <ac.h>

KSEQ_INIT(gzFile, gzread)


//codons listed in alphabetical order
char codonTable1[64] = {
  'K', 'N', 'K', 'N', 'T', 'T', 'T', 'T', 'R', 'S', 'R', 'S', 'I',
  'I', 'M', 'I', 'Q', 'H', 'Q', 'H', 'P', 'P', 'P', 'P', 'R', 'R',
  'R', 'R', 'L', 'L', 'L', 'L', 'E', 'D', 'E', 'D', 'A', 'A', 'A',
  'A', 'G', 'G', 'G', 'G', 'V', 'V', 'V', 'V', '*', 'Y', '*', 'Y',
  'S', 'S', 'S', 'S', '*', 'C', 'W', 'C', 'L', 'F', 'L', 'F'};


inline char translate_codon(char* codon, char* codonTable){
  int i;
  int index=0;
  bool is_n = false;
  for (i=0; i<=2; i++){
    switch (codon[i]){
      case 'A':
      case 'a':
        break;
      case 'C':
      case 'c':
        index += 1 << ((2-i)*2);
        break;
      case 'G':
      case 'g':
        index += 2 << ((2-i)*2);
        break;
      case 'T':
      case 't':
        index += 3 << ((2-i)*2);
        break;
      case 'N':
      case 'n':
        is_n = true;
        break;
      default:
        fprintf(stderr, "Detected unexcepted codon: %s\n", codon);
        exit(1);
    }
  }

  if (is_n){
    return 'X';
  } else {
    return codonTable[index];
  }
}


inline char translate_codon_revcom(char* codon, char* codonTable){
  int i;
  int index=0;
  bool is_n = false;
  for (i=0; i<=2; i++){
    switch (codon[i]){
      case 'A':
      case 'a':
        index += 3 << (i*2);
        break;
      case 'C':
      case 'c':
        index += 2 << (i*2);
        break;
      case 'G':
      case 'g':
        index += 1 << (i*2);
        break;
      case 'T':
      case 't':
        break;
      case 'N':
      case 'n':
        is_n = true;
        break;
      default:
        fprintf(stderr, "Detected unexcepted codon: %s\n", codon);
        exit(1);
    }
  }

  if (is_n){
    return 'X';
  } else {
    return codonTable[index];
  }
}


void translate(char* begin, int num, bool reverse, char* codonTable){
//   char* tmp = malloc(10000); printf("%s %i\n", strncpy(tmp, begin, num), reverse); free(tmp);
  if (reverse){
    num -= 3;
    for(;num >= 0; num-=3){
      putchar(translate_codon_revcom(begin+num, codonTable));
    }
  } else {
    int i;
    for(i=0; i<num; i+=3){
      putchar(translate_codon(begin+i, codonTable));
    }
  }
  puts("");
}

inline void print_sequence_header(kseq_t* seq_struct, int start_position, int frame, int *orf_counter){
  if (seq_struct->comment.l > 0){
    printf(">%s_%i_%i_%i %s\n", seq_struct->name.s, start_position+1, frame, (*orf_counter)++, seq_struct->comment.s);
  } else {
    printf(">%s_%i_%i_%i\n", seq_struct->name.s, start_position+1, frame, (*orf_counter)++);
  }
}


void process_sequence_file(char *path, int min_length, char* codonTable){
  //create string searching structure made of stop codons in fwd mode
  //create reverse searching structure
  AC_STRUCT * ac;
  ac = ac_alloc();
  ac_add_string(ac, "TAA", 3, 1);
  ac_add_string(ac, "TGA", 3, 2);
  ac_add_string(ac, "TAG", 3, 3);
  ac_add_string(ac, "TTA", 3, 4); //revcom
  ac_add_string(ac, "TCA", 3, 5); //revcom
  ac_add_string(ac, "CTA", 3, 6); //revcom
  ac_prep(ac);
  char *search_result;
  int length_out, id_out, ends_at;
  int mod3;
  int length_to_string_start;

  //setup kseq reading
  gzFile fp;
  if (path == NULL){
    fp = gzopen("/dev/stdin","r");
  } else {
    fp = gzopen(path, "r");
  }

  //read in one sequence at a time, looping
  kseq_t *seq;
  int l;
  seq = kseq_init(fp);

  while ((l = kseq_read(seq)) >= 0) {
    //printf("Processing sequence %s\n",seq->seq.s);

    //set current positions for each of the frames
    int last_found_positions[6] = {0,1,2,0,1,2};
    int orf_counter = 1;

    //search for forward facing stop codons
    ac_search_init(ac, seq->seq.s, seq->seq.l);

    //for each found position,
    while ((search_result = ac_search(ac, &length_out, &id_out, &ends_at)) != NULL){
      //find the frame of the position using mod 3 operation
      length_to_string_start = search_result - seq->seq.s;
      mod3 = length_to_string_start % 3;
//       printf("Found potential ORF at position %i given found string id %i. Current positions are %i,%i,%i %i,%i,%i and mod3 %i from total sequence length %i and ORF length %i\n",
//              length_to_string_start,
//              id_out,
//              last_found_positions[0],
//              last_found_positions[1],
//              last_found_positions[2],
//              last_found_positions[3],
//              last_found_positions[4],
//              last_found_positions[5],
//              mod3,
//              (int)seq->seq.l,
//              length_to_string_start - last_found_positions[mod3]
//             );

      //if current position - last position >= min_length, translate and spit out translated sequence
      if (id_out <= 3) {
        //in fwd direction is this ORF
        if (length_to_string_start - last_found_positions[mod3] >= min_length){
          //spit out sequence
          print_sequence_header(seq, last_found_positions[mod3], mod3+1, &orf_counter);
          translate(seq->seq.s+last_found_positions[mod3],
                         length_to_string_start - last_found_positions[mod3],
                         false, codonTable);
        }
        //last_position[frame] = current_position (for next time)
        last_found_positions[mod3] = length_to_string_start+3;

      } else { //if id of ac hit is 3 or more, spit out the reverse complement
        if (length_to_string_start - last_found_positions[mod3+3] >= min_length){
          //spit out revcom sequence
          print_sequence_header(seq, last_found_positions[mod3+3], mod3+4, &orf_counter);
          translate(seq->seq.s+last_found_positions[mod3+3],
                         length_to_string_start - last_found_positions[mod3+3],
                         true, codonTable);

        }
        last_found_positions[mod3+3] = length_to_string_start+3;
      }
    } //end while finding matches



    //then check for the last ORF in each of the frames
    mod3 = seq->seq.l % 3;
    //     printf("Final positions: %i,%i,%i %i,%i,%i and mod3 %i from length %i\n",
    //            last_found_positions[0],
    //            last_found_positions[1],
    //            last_found_positions[2],
    //            last_found_positions[3],
    //            last_found_positions[4],
    //            last_found_positions[5],
    //            mod3,
    //            (int)seq->seq.l
    //           );
    switch (mod3){
      case 0:
        //translate each of the frames in order
        //for the 1st and 4th frames translate (length-last)bp
        //for the 2,3,5,6th frames translate (length-last-3-{2,1})bp
        if ((l = seq->seq.l - last_found_positions[0]) >= min_length){
          print_sequence_header(seq, last_found_positions[0], 1, &orf_counter);
          translate(seq->seq.s+last_found_positions[0], l, false, codonTable);
        }
        if ((l = seq->seq.l - last_found_positions[1] - 2) >= min_length){
          print_sequence_header(seq, last_found_positions[1], 2, &orf_counter);
          translate(seq->seq.s+last_found_positions[1], l, false, codonTable);
        }
        if ((l = seq->seq.l - last_found_positions[2] - 1) >= min_length){
          print_sequence_header(seq, last_found_positions[2], 3, &orf_counter);
          translate(seq->seq.s+last_found_positions[2], l, false, codonTable);
        }
        if ((l = seq->seq.l - last_found_positions[3]) >= min_length){
          print_sequence_header(seq, last_found_positions[3], 4, &orf_counter);
          translate(seq->seq.s+last_found_positions[3], l, true, codonTable);
        }
        if ((l = seq->seq.l - last_found_positions[4] - 2) >= min_length){
          print_sequence_header(seq, last_found_positions[4], 5, &orf_counter);
          translate(seq->seq.s+last_found_positions[4], l, true, codonTable);
        }
        if ((l = seq->seq.l - last_found_positions[5] - 1) >= min_length){
          print_sequence_header(seq, last_found_positions[5], 6, &orf_counter);
          translate(seq->seq.s+last_found_positions[5], l, true, codonTable);
        }
        break;
      case 1:
        if ((l = seq->seq.l - last_found_positions[0] - 1) >= min_length){
          print_sequence_header(seq, last_found_positions[0], 1, &orf_counter);
          translate(seq->seq.s+last_found_positions[0], l, false, codonTable);
        }
        if ((l = seq->seq.l - last_found_positions[1]) >= min_length){
          print_sequence_header(seq, last_found_positions[1], 2, &orf_counter);
          translate(seq->seq.s+last_found_positions[1], l, false, codonTable);
        }
        if ((l = seq->seq.l - last_found_positions[2] - 2) >= min_length){
          print_sequence_header(seq, last_found_positions[2], 3, &orf_counter);
          translate(seq->seq.s+last_found_positions[2], l, false, codonTable);
        }
        if ((l = seq->seq.l - last_found_positions[3] - 1) >= min_length){
          print_sequence_header(seq, last_found_positions[3], 4, &orf_counter);
          translate(seq->seq.s+last_found_positions[3], l, true, codonTable);
        }
        if ((l = seq->seq.l - last_found_positions[4]) >= min_length){
          print_sequence_header(seq, last_found_positions[4], 5, &orf_counter);
          translate(seq->seq.s+last_found_positions[4], l, true, codonTable);
        }
        if ((l = seq->seq.l - last_found_positions[5] - 2) >= min_length){
          print_sequence_header(seq, last_found_positions[5], 6, &orf_counter);
          translate(seq->seq.s+last_found_positions[5], l, true, codonTable);
        }
        break;
      case 2:
        if ((l = seq->seq.l - last_found_positions[0] - 2) >= min_length){
          print_sequence_header(seq, last_found_positions[0], 1, &orf_counter);
          translate(seq->seq.s+last_found_positions[0], l, false, codonTable);
        }
        if ((l = seq->seq.l - last_found_positions[1] - 1) >= min_length){
          print_sequence_header(seq, last_found_positions[1], 2, &orf_counter);
          translate(seq->seq.s+last_found_positions[1], l, false, codonTable);
        }
        if ((l = seq->seq.l - last_found_positions[2]) >= min_length){
          print_sequence_header(seq, last_found_positions[2], 3, &orf_counter);
          translate(seq->seq.s+last_found_positions[2], l, false, codonTable);
        }
        if ((l = seq->seq.l - last_found_positions[3] - 2) >= min_length){
          print_sequence_header(seq, last_found_positions[3], 4, &orf_counter);
          translate(seq->seq.s+last_found_positions[3], l, true, codonTable);
        }
        if ((l = seq->seq.l - last_found_positions[4] - 1) >= min_length){
          print_sequence_header(seq, last_found_positions[4], 5, &orf_counter);
          translate(seq->seq.s+last_found_positions[4], l, true, codonTable);
        }
        if ((l = seq->seq.l - last_found_positions[5]) >= min_length){
          print_sequence_header(seq, last_found_positions[5], 6, &orf_counter);
          translate(seq->seq.s+last_found_positions[5], l, true, codonTable);
        }
        break;
    }
  } //onto next sequence

  //cleanup
  kseq_destroy(seq);
  ac_free(ac);
  gzclose(fp);
}


int main(int argc, char *argv[]){
  int min_length = 96;
  char c;
  char* codonTable = codonTable1;

  while ((c = getopt (argc, argv, "hvm:")) != -1){
    switch (c){
      case 'v':
        printf("OrfM version 0.1\n");
        exit(0);
      case 'h':
        printf("\n  Usage: orfm [options] <seq_file>\n\n");
        printf("  The <seq_file> can be a FASTA or FASTQ file, gzipped or uncompressed.\n\n");
        printf("  Options:\n");
        printf("   -m LENGTH   minimum number of nucleotides (not amino acids) to call an ORF on [default: %i]\n", min_length);
        printf("   -v          show version information\n");
        printf("   -h          show this help\n");
        printf("\n");
        exit(0);
      case 'm':
        min_length = atoi(optarg);
        if (min_length < 3){
          fprintf(stderr, "ERROR: -m minimum length must be 3 or more\n");
          exit(1);
        }
        if (min_length % 3 != 0){
          fprintf(stderr, "ERROR: -m minimum length argument must be a multiple of 3\n");
          exit(2);
        }
        break;
    }
  }
  //printf("Processing sequence file %s with min_length %i\n", argv[optind], min_length);
  if (argc-optind == 0){
    process_sequence_file(NULL, min_length, codonTable);
  } else if (argc-optind == 1){
    process_sequence_file(argv[optind], min_length, codonTable);
  } else {
    fprintf(stderr, "ERROR: one file at most can be given as an argument, found %i\n", argc-optind);
    exit(3);
  }

  return 0;
}
