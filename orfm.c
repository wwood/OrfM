#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <getopt.h>    /* for getopt */
#include <string.h>

#include <zlib.h>
#include <seqtk/kseq.h>

#include <ac.h>

KSEQ_INIT(gzFile, gzread)


char codonTable1[64] = {
  'K', 'N', 'K', 'N', 'T', 'T', 'T', 'T', 'R', 'S', 'R', 'S', 'I',
  'I', 'M', 'I', 'Q', 'H', 'Q', 'H', 'P', 'P', 'P', 'P', 'R', 'R',
  'R', 'R', 'L', 'L', 'L', 'L', 'E', 'D', 'E', 'D', 'A', 'A', 'A',
  'A', 'G', 'G', 'G', 'G', 'V', 'V', 'V', 'V', '*', 'Y', '*', 'Y',
  'S', 'S', 'S', 'S', '*', 'C', 'W', 'C', 'L', 'F', 'L', 'F'};


char* for_printing;


void translate(char* begin, int num, bool reverse, char* codonTable){
//   char* tmp = malloc(10000); printf("%s %i\n", strncpy(tmp, begin, num), reverse); free(tmp);
  for_printing = realloc(for_printing, (num/3) + 1);
  int j=0;
  if (reverse){
    num -= 3;
    for(;num >= 0; num-=3){
      for_printing[j] = codonTable[ 128 + (*(begin+num+2)<<16) + (*(begin+num+1)<<8) + *(begin+num)];
      j++;
    }
  } else {
    int i;
    for(i=0; i<num; i+=3){
      for_printing[j] = codonTable[ (*(begin+i)<<16) + (*(begin+i+1)<<8) + *(begin+i+2)];
      j++;
    }
  }
  for_printing[j] = 0;
  puts(for_printing);
}

inline void print_sequence_header(kseq_t* seq_struct, int start_position, int frame, int *orf_counter){
  if (seq_struct->comment.l > 0){
    printf(">%s_%i_%i_%i %s\n", seq_struct->name.s, start_position+1, frame, (*orf_counter)++, seq_struct->comment.s);
  } else {
    printf(">%s_%i_%i_%i\n", seq_struct->name.s, start_position+1, frame, (*orf_counter)++);
  }
}

/** Create an index that maps 3 base codons (the bits from 3 chars)
    to the amino acid that they encode, with X for non-DNA bases.
    The reverse complements are encoded as +128 since A, T, G, C, a, t, g, c are
    all less than 128 bits.
    */
char* create_codon_table(char* input_codon_table){
  int table_size = (1<<8)*(1<<8)*(1<<8) * sizeof(char); //3 bases per codon indexed to another char
  //printf("table size %i\n", table_size);
  char* faster_table = (char*) malloc(table_size);

  char* order = "ACGT";
  char* low_order = "acgt";
  char* complements = "TGCA";
  char* low_complements = "tgca";
  int i,j,k,count;
  for (i=0; i<table_size; i++)
    faster_table[i] = 'X';
  count = 0;
  for (i=0; i<4; i++){
    for (j=0; j<4; j++){
      for (k=0; k<4; k++){
        //fwd
        faster_table[(order[i]<<16) + (order[j]<<8) + order[k]] = input_codon_table[count];
        faster_table[(order[i]<<16) + (order[j]<<8) + low_order[k]] = input_codon_table[count];
        faster_table[(order[i]<<16) + (low_order[j]<<8) + order[k]] = input_codon_table[count];
        faster_table[(order[i]<<16) + (low_order[j]<<8) + low_order[k]] = input_codon_table[count];
        faster_table[(low_order[i]<<16) + (low_order[j]<<8) + order[k]] = input_codon_table[count];
        faster_table[(low_order[i]<<16) + (low_order[j]<<8) + low_order[k]] = input_codon_table[count];
        faster_table[(low_order[i]<<16) + (order[j]<<8) + order[k]] = input_codon_table[count];
        faster_table[(low_order[i]<<16) + (order[j]<<8) + low_order[k]] = input_codon_table[count];

        //encode the reverse complement of the codon matching to the same amino acid
        faster_table[128 + (complements[i]<<16) + (complements[j]<<8) + complements[k]] = input_codon_table[count];
        faster_table[128 + (complements[i]<<16) + (complements[j]<<8) + low_complements[k]] = input_codon_table[count];
        faster_table[128 + (complements[i]<<16) + (low_complements[j]<<8) + complements[k]] = input_codon_table[count];
        faster_table[128 + (complements[i]<<16) + (low_complements[j]<<8) + low_complements[k]] = input_codon_table[count];
        faster_table[128 + (low_complements[i]<<16) + (complements[j]<<8) + complements[k]] = input_codon_table[count];
        faster_table[128 + (low_complements[i]<<16) + (complements[j]<<8) + low_complements[k]] = input_codon_table[count];
        faster_table[128 + (low_complements[i]<<16) + (low_complements[j]<<8) + complements[k]] = input_codon_table[count];
        faster_table[128 + (low_complements[i]<<16) + (low_complements[j]<<8) + low_complements[k]] = input_codon_table[count];

        count += 1;
      }
    }
  }

  return faster_table;
}


void process_sequence_file(char *path, int min_length, char* codonTable, int position_limit){
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

  for_printing = (char *) malloc(1);

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
  int read_position_limit;
  codonTable = create_codon_table(codonTable);

  while ((l = kseq_read(seq)) >= 0) {
    //if there is no chance of getting an ORF here, don't even start
    if (seq->seq.l < min_length){
      continue;
    }
    //printf("Processing sequence %s\n",seq->seq.s);
    if (!position_limit || seq->seq.l < position_limit){
      read_position_limit = seq->seq.l;
    } else {
      read_position_limit = position_limit;
    }

    //set current positions for each of the frames
    int last_found_positions[6] = {0,1,2,0,1,2};
    int orf_counter = 1;

    //search for forward facing stop codons
    ac_search_init(ac, seq->seq.s, read_position_limit);

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
      //              (int)read_position_limit,
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
    mod3 = read_position_limit % 3;
    //     printf("Final positions: %i,%i,%i %i,%i,%i and mod3 %i from length %i\n",
    //            last_found_positions[0],
    //            last_found_positions[1],
    //            last_found_positions[2],
    //            last_found_positions[3],
    //            last_found_positions[4],
    //            last_found_positions[5],
    //            mod3,
    //            (int)read_position_limit
    //           );
    switch (mod3){
      case 0:
        //translate each of the frames in order
        //for the 1st and 4th frames translate (length-last)bp
        //for the 2,3,5,6th frames translate (length-last-3-{2,1})bp
        if ((l = read_position_limit - last_found_positions[0]) >= min_length){
          print_sequence_header(seq, last_found_positions[0], 1, &orf_counter);
          translate(seq->seq.s+last_found_positions[0], l, false, codonTable);
        }
        if ((l = read_position_limit - last_found_positions[1] - 2) >= min_length){
          print_sequence_header(seq, last_found_positions[1], 2, &orf_counter);
          translate(seq->seq.s+last_found_positions[1], l, false, codonTable);
        }
        if ((l = read_position_limit - last_found_positions[2] - 1) >= min_length){
          print_sequence_header(seq, last_found_positions[2], 3, &orf_counter);
          translate(seq->seq.s+last_found_positions[2], l, false, codonTable);
        }
        if ((l = read_position_limit - last_found_positions[3]) >= min_length){
          print_sequence_header(seq, last_found_positions[3], 4, &orf_counter);
          translate(seq->seq.s+last_found_positions[3], l, true, codonTable);
        }
        if ((l = read_position_limit - last_found_positions[4] - 2) >= min_length){
          print_sequence_header(seq, last_found_positions[4], 5, &orf_counter);
          translate(seq->seq.s+last_found_positions[4], l, true, codonTable);
        }
        if ((l = read_position_limit - last_found_positions[5] - 1) >= min_length){
          print_sequence_header(seq, last_found_positions[5], 6, &orf_counter);
          translate(seq->seq.s+last_found_positions[5], l, true, codonTable);
        }
        break;
      case 1:
        if ((l = read_position_limit - last_found_positions[0] - 1) >= min_length){
          print_sequence_header(seq, last_found_positions[0], 1, &orf_counter);
          translate(seq->seq.s+last_found_positions[0], l, false, codonTable);
        }
        if ((l = read_position_limit - last_found_positions[1]) >= min_length){
          print_sequence_header(seq, last_found_positions[1], 2, &orf_counter);
          translate(seq->seq.s+last_found_positions[1], l, false, codonTable);
        }
        if ((l = read_position_limit - last_found_positions[2] - 2) >= min_length){
          print_sequence_header(seq, last_found_positions[2], 3, &orf_counter);
          translate(seq->seq.s+last_found_positions[2], l, false, codonTable);
        }
        if ((l = read_position_limit - last_found_positions[3] - 1) >= min_length){
          print_sequence_header(seq, last_found_positions[3], 4, &orf_counter);
          translate(seq->seq.s+last_found_positions[3], l, true, codonTable);
        }
        if ((l = read_position_limit - last_found_positions[4]) >= min_length){
          print_sequence_header(seq, last_found_positions[4], 5, &orf_counter);
          translate(seq->seq.s+last_found_positions[4], l, true, codonTable);
        }
        if ((l = read_position_limit - last_found_positions[5] - 2) >= min_length){
          print_sequence_header(seq, last_found_positions[5], 6, &orf_counter);
          translate(seq->seq.s+last_found_positions[5], l, true, codonTable);
        }
        break;
      case 2:
        if ((l = read_position_limit - last_found_positions[0] - 2) >= min_length){
          print_sequence_header(seq, last_found_positions[0], 1, &orf_counter);
          translate(seq->seq.s+last_found_positions[0], l, false, codonTable);
        }
        if ((l = read_position_limit - last_found_positions[1] - 1) >= min_length){
          print_sequence_header(seq, last_found_positions[1], 2, &orf_counter);
          translate(seq->seq.s+last_found_positions[1], l, false, codonTable);
        }
        if ((l = read_position_limit - last_found_positions[2]) >= min_length){
          print_sequence_header(seq, last_found_positions[2], 3, &orf_counter);
          translate(seq->seq.s+last_found_positions[2], l, false, codonTable);
        }
        if ((l = read_position_limit - last_found_positions[3] - 2) >= min_length){
          print_sequence_header(seq, last_found_positions[3], 4, &orf_counter);
          translate(seq->seq.s+last_found_positions[3], l, true, codonTable);
        }
        if ((l = read_position_limit - last_found_positions[4] - 1) >= min_length){
          print_sequence_header(seq, last_found_positions[4], 5, &orf_counter);
          translate(seq->seq.s+last_found_positions[4], l, true, codonTable);
        }
        if ((l = read_position_limit - last_found_positions[5]) >= min_length){
          print_sequence_header(seq, last_found_positions[5], 6, &orf_counter);
          translate(seq->seq.s+last_found_positions[5], l, true, codonTable);
        }
        break;
    }
  } //onto next sequence

  //cleanup
  free(codonTable);
  free(for_printing);
  kseq_destroy(seq);
  ac_free(ac);
  gzclose(fp);
}

/* Return true if the current_version (i.e. ORFM_VERSION) is >= required_version.
The input strings are modified. Not thread-safe since strtok is not thread-safe. */
bool compare_version(char* required_version, char* current_version){
  char *s1;
  char *s2;
  char *s3;
  int r1, r2, r3, c1, c2, c3;
  s1 = (char *) strtok(required_version, ".");
  s2 = (char *) strtok(NULL, ".");
  s3 = (char *) strtok(NULL, ".");
  if (s1 == NULL || s2 == NULL || s3 == NULL || strtok(NULL, ".") != NULL){
    fprintf(stderr, "Unexpected format of version string, I required something like 0.1.4\n");
    exit(4);
  }
  r1 = atoi(s1);
  r2 = atoi(s2);
  r3 = atoi(s3);

  char *current_version2 = malloc((1+strlen(current_version))*sizeof(char));
  strcpy(current_version2, current_version);
  c1 = atoi(strtok(current_version2, "."));
  c2 = atoi(strtok(NULL, "."));
  c3 = atoi(strtok(NULL, "."));
  free(current_version2);

  if (r1 < c1) return true;
  else if (r1 > c1) return false;
  if (r2 < c2) return true;
  else if (r2 > c2) return false;
  if (r3 <= c3) return true;
  else return false;
}


int main(int argc, char *argv[]){
  int min_length = 96;
  char c;
  char* codonTable = codonTable1;
  int position_limit = 0;
  char* required_version;

  while ((c = getopt(argc, argv, "hvm:l:r:")) != -1){
    switch (c){
      case 'v':
        printf("OrfM version %s\n",ORFM_VERSION);
        exit(0);
      case 'h':
        printf("\n  Usage: orfm [options] <seq_file>\n\n");
        printf("  The <seq_file> can be a FASTA or FASTQ file, gzipped or uncompressed.\n\n");
        printf("  Options:\n");
        printf("   -m LENGTH   minimum number of nucleotides (not amino acids) to call an ORF on [default: %i]\n", min_length);
        printf("   -l LENGTH   ignore the sequence of the read beyond this, useful when comparing reads from with different read lengths [default: none]\n");
        printf("   -v          show version information\n");
        printf("   -h          show this help\n");
        printf("   -r VERSION  do not run unless this version of OrfM is at least this version number (e.g. %s)\n",ORFM_VERSION);
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
      case 'l':
        position_limit = atoi(optarg);
        if (position_limit < 3){
          fprintf(stderr, "ERROR: -l must be 3 or more\n");
          exit(1);
        }
        break;
      case 'r':
        required_version = malloc((1+strlen(optarg))*sizeof(char));
        strcpy(required_version, optarg);
        if (!compare_version(required_version, ORFM_VERSION)){
          fprintf(stderr, "ERROR: this version of OrfM is older than the required version; please upgrade OrfM\n");
          free(required_version);
          exit(3);
        }
        free(required_version);
        break;
    }
  }

  if (position_limit != 0 && min_length > position_limit){
    fprintf(stderr, "-l cannot be greater than -m, otherwise no ORFs are possible\n");
    exit(4);
  }
  //printf("Processing sequence file %s with min_length %i\n", argv[optind], min_length);
  if (argc-optind == 0){
    process_sequence_file(NULL, min_length, codonTable, position_limit);
  } else if (argc-optind == 1){
    process_sequence_file(argv[optind], min_length, codonTable, position_limit);
  } else {
    fprintf(stderr, "ERROR: one file at most can be given as an argument, found %i\n", argc-optind);
    exit(3);
  }

  return 0;
}
