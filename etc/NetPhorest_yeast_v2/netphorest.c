#include <assert.h>
#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Modules to include */
#define INSR
#define NN
#define PSSM

/* Define size of the alphabet */
#define na 21

/* Define max number of neurons */
#define maxwindow 99
#define center_of_maxwindow 49
#define maxhidden 20
#define maxname 1024


/* nn_data.h contains static_assert calls to verify the length of the arrays.
 * static_assert is available since C++11 only and must be replaced by assert
 * in C and C++98 and C++03.
 * Compiling as C++11 allows the static_asserts to be verified at compile time.
 * Compile with GCC >= 4.7 to activate it:
 * g++ -std=c++11 -O3 -o netphorest netphorest.c -lm
 * Xavier - 2017-05-01
 */
#if __cplusplus < 201103L /* __cplusplus will be replaced by 0 in C */
#define STATIC_ASSERT(x, y) assert(x)
#else
#define STATIC_ASSERT(x, y) static_assert(x, y)
#endif


char alphabet[] = "FIVWMLCHYAGNRTPDEQSK-UXZJBO"; /* To deal with unsupported amino acid characters - jinho 2014/01/30 */

float sigmoid_data[256] = {
  0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
  0.000000, 0.000000, 0.000000, 0.000000, 0.000001, 0.000001, 0.000001, 0.000001,
  0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000002, 0.000002, 0.000002,
  0.000002, 0.000003, 0.000003, 0.000003, 0.000004, 0.000004, 0.000005, 0.000005,
  0.000006, 0.000007, 0.000008, 0.000009, 0.000010, 0.000011, 0.000013, 0.000015,
  0.000017, 0.000019, 0.000021, 0.000024, 0.000028, 0.000031, 0.000035, 0.000040,
  0.000045, 0.000051, 0.000058, 0.000066, 0.000075, 0.000085, 0.000096, 0.000109,
  0.000123, 0.000140, 0.000158, 0.000180, 0.000203, 0.000231, 0.000261, 0.000296,
  0.000335, 0.000380, 0.000431, 0.000488, 0.000553, 0.000626, 0.000710, 0.000804,
  0.000911, 0.001032, 0.001170, 0.001325, 0.001501, 0.001701, 0.001927, 0.002183,
  0.002473, 0.002801, 0.003173, 0.003594, 0.004070, 0.004610, 0.005220, 0.005911,
  0.006693, 0.007577, 0.008577, 0.009708, 0.010987, 0.012432, 0.014064, 0.015906,
  0.017986, 0.020332, 0.022977, 0.025957, 0.029312, 0.033086, 0.037327, 0.042088,
  0.047426, 0.053403, 0.060087, 0.067547, 0.075858, 0.085099, 0.095349, 0.106691,
  0.119203, 0.132964, 0.148047, 0.164516, 0.182426, 0.201813, 0.222700, 0.245085,
  0.268941, 0.294215, 0.320821, 0.348645, 0.377541, 0.407333, 0.437823, 0.468791,
  0.500000, 0.531209, 0.562177, 0.592667, 0.622459, 0.651355, 0.679179, 0.705785,
  0.731059, 0.754915, 0.777300, 0.798187, 0.817574, 0.835484, 0.851953, 0.867036,
  0.880797, 0.893309, 0.904651, 0.914901, 0.924142, 0.932453, 0.939913, 0.946597,
  0.952574, 0.957912, 0.962673, 0.966914, 0.970688, 0.974043, 0.977023, 0.979668,
  0.982014, 0.984094, 0.985936, 0.987568, 0.989013, 0.990292, 0.991423, 0.992423,
  0.993307, 0.994089, 0.994780, 0.995390, 0.995930, 0.996406, 0.996827, 0.997199,
  0.997527, 0.997817, 0.998073, 0.998299, 0.998499, 0.998675, 0.998830, 0.998968,
  0.999089, 0.999196, 0.999290, 0.999374, 0.999447, 0.999512, 0.999569, 0.999620,
  0.999665, 0.999704, 0.999739, 0.999769, 0.999797, 0.999820, 0.999842, 0.999860,
  0.999877, 0.999891, 0.999904, 0.999915, 0.999925, 0.999934, 0.999942, 0.999949,
  0.999955, 0.999960, 0.999965, 0.999969, 0.999972, 0.999976, 0.999979, 0.999981,
  0.999983, 0.999985, 0.999987, 0.999989, 0.999990, 0.999991, 0.999992, 0.999993,
  0.999994, 0.999995, 0.999995, 0.999996, 0.999996, 0.999997, 0.999997, 0.999997,
  0.999998, 0.999998, 0.999998, 0.999998, 0.999999, 0.999999, 0.999999, 0.999999,
  0.999999, 0.999999, 0.999999, 0.999999, 0.999999, 1.000000, 1.000000, 1.000000,
  1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000
};


int gDistance_to_unsupported_amino_acid = center_of_maxwindow + 1;
int gbHas_unsupported_amino_acid = 0;

/*
 *  Fast sigmoid function that uses a lookup table.
 */

float sigmoid(float x) {

  if (x <= -16) {
    return(0);
  } else if (x >= 16) {
    return(1);
  } else {
    return(sigmoid_data[(int)(8*x+128)]);
  }

}




/*
 *  Perform feed forward evaluation of a neural network on a sequence window.
 */

float feed_forward(int const *s, float const w[], int nw, int nh) {

  float h[maxhidden], o[2], x;
  int i, j;

  /* To deal with unsupported amino acid characters - jinho 2014/01/30 */
  if (gDistance_to_unsupported_amino_acid <= (nw-1)/2){
    gbHas_unsupported_amino_acid = 1;
    return 0;
  }
  s += (maxwindow-nw)/2;

  if (nh > 0) {
    for (i = 0; i < nh; ++i) {
      x = w[(na*nw+1)*(i+1)-1];
      for (j = 0; j < nw; ++j) {
        x += w[(na*nw+1)*i+na*j+s[j]];
      }
      h[i] = sigmoid(x);
    }
    for (i = 0; i <= 1; ++i) {
      x = w[(na*nw+1)*nh+(nh+1)*(i+1)-1];
      for (j = 0; j < nh; ++j) {
        x += w[(na*nw+1)*nh+(nh+1)*i+j]*h[j];
      }
      o[i] = sigmoid(x);
    }
  } else {
    for (i = 0; i <= 1; ++i) {
      x = w[(na*nw+1)*(i+1)-1];
      for (j = 0; j < nw; ++j) {
        x += w[(na*nw+1)*i+na*j+s[j]];
      }
      o[i] = sigmoid(x);
    }
  }

  /* Combine the scores from the two output neurons */
  return((o[0]+1-o[1])/2);

}




/*
 *
 */

float pssm(int const *s, float const w[], int nw) {

  float o;
  int i;

  /* To deal with unsupported amino acid characters - jinho 2014/01/30 */
  if (gDistance_to_unsupported_amino_acid <= (nw-1)/2){
    return 0;
  }
  
  s += (maxwindow-nw)/2;

  o = 1;
  for (i = 0; i < nw; ++i) {
    o *= w[na*(nw-i-1)+s[i]];
  }

  return(o);

}



/*
 *  Print peptide sequence.
 */

void print_peptide (int const *s) {

  int i;

  s += (maxwindow-11)/2;

  for (i = 10; i >= 6; --i) {
    fputc(toupper(alphabet[s[i]]), stdout);
  }
  fputc(tolower(alphabet[s[5]]), stdout);
  for (i = 4; i >= 0; --i) {
    fputc(toupper(alphabet[s[i]]), stdout);
  }

}

/* To deal with unsupported amino acid characters - jinho 2014/01/30 */
int GetDistanceToUnsupportedAminoAcid(int const *s){
  int i;
  for (i = 1; i < center_of_maxwindow+1; i++){
    if (s[center_of_maxwindow+i] >= na || s[center_of_maxwindow-i] >= na){
      return i;
    }
  }
  return center_of_maxwindow + 1;
}

/*
 *  Calculate and print scores for a sequence window.
 */

void predict(char *name, int *n, int const *s) {

  #ifdef NETPHOS
  #include "netphos_data.h"
  #endif

  #ifdef NETPHOSK
  #include "netphosk_data.h"
  #endif

  #ifdef INSR
  #include "insr_data.h"
  #endif

  #ifdef NN
  #include "nn_data.h"
  #endif

  #ifdef PSSM
  #include "pssm_data.h"
  #endif

  float o;
  int c;

  if (s[center_of_maxwindow] == na-1) {
    return;
  }

  (*n)++;

  c = alphabet[s[center_of_maxwindow]];

  if (c != 'S' && c != 'T' && c != 'Y') {
    return;
  }
  
  gDistance_to_unsupported_amino_acid = GetDistanceToUnsupportedAminoAcid(s); /* To deal with unsupported amino acid characters - jinho 2014/01/30 */

  #ifdef NETPHOS
  #include "netphos_code.h"
  #endif

  #ifdef NETPHOSK
  #include "netphosk_code.h"
  #endif

  #ifdef INSR
  #include "insr_code.h"
  #endif

  #ifdef NN
  #include "nn_code.h"
  #endif

  #ifdef PSSM
  #include "pssm_code.h"
  #endif

}




int main(int ARGC, char *ARGV[]) {

  char name[maxname + 1], *p;
  int c, i, j, n, pos, s[maxwindow];

  fputs("# Name\tPosition\tResidue\tPeptide\tMethod\tTree\tClassifier\tPosterior\tPrior\n", stdout);

  c = fgetc(stdin);
  while (!feof(stdin)) {

    /* Skip to fasta header */
    while (!feof(stdin) && c != '>') {
      c = fgetc(stdin);
    }

    /* Read name */
    n = 0;
    c = fgetc(stdin);
    while (!feof(stdin) && !isspace(c) && n < maxname) {
      name[n] = c;
      n++;
      c = fgetc(stdin);
    }
    name[n] = '\0';

    /* Skip rest of header line */
    while (!feof(stdin) && c != '\n') {
      c = fgetc(stdin);
    }

    /* Initialize sequence window */
    for (i = 0; i < maxwindow; ++i) {
      s[i] = na-1;
    }

    /* Get first character after header eol */
    c = fgetc(stdin);

    /* Read and process sequence */
    n = 0;
    pos = 0; /* first character is previous \n */
    while (!feof(stdin) && c != '>') {
      ++pos;
      if (!isascii(c)) {
          fprintf(stderr, "Non-ascii character found in sequence %s, position %d\n", name, pos);
	  return(1);
      }
      if (isalpha(c)) {
        p = strchr(alphabet, c);
        if (p != NULL) {
          for (i = 1; i < maxwindow; i++) {
            s[maxwindow-i] = s[maxwindow-i-1];
          }
          s[0] = p-alphabet;
          predict(name, &n, s);
        }
      }
      c = fgetc(stdin);
    }

    /* Process sequence tail */
    for (i = 0; i < (maxwindow-1)/2; i++) {
      for (j = 1; j < maxwindow; j++) {
        s[maxwindow-j] = s[maxwindow-j-1];
      }
      s[0] = na-1;
      predict(name, &n, s);
    }

  }

  return(0);

}
