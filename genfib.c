#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

#include <math.h>

#include "convol.h"
#include <fftw3.h>



void print_vec(double *v) {
  size_t i;

  for(i = 0; i < 2 * N; i++)
    printf("  %3d  %12f\n", i, v[i]);
  printf("\n");
}

void print_2vec(double *v1, double sl1, double *v2, double sl2) {
  size_t i;

  for(i = 0; i < N; i++)
    printf("  %3d  %12f  %12f\n", i, v1[i] * exp(sl1 * i), v2[i] * exp(sl2 * i));
  printf("\n");
}


void init_weights(size_t N, double *w) {
  size_t i;

  for(i = 1; i < N; i++) 
    w[i] = 1./i;
  w[0] = 0;
}

double sla=0, slb=0;

void new_slope(double *sl, double *v, double new_sl) {
  size_t i;

  *sl += new_sl;
  for(i=0; i<N; i++)
    v[i] = v[i] * exp(-new_sl*i);
}

void fix_slope(double *sl, double *v, int nr, int nl) {
  double a = (log(v[nr]) - log(v[nl]));
  double b = nr-nl;
  printf("[%i %i %f %f %f]\n", nr-nl, nr, a , b, a / b);
  new_slope(sl, v, (log(v[nr]) - log(v[nl])) / (nr-nl));
}

int main() {
  size_t i,valid;
  double b0;

  init_conv(128);

  init_weights(N, inb);  // b = w
  memset(ina, 0, sizeof(ina[0]) * 2 * N); ina[0] = 1;  // a = delta
  valid = 1;

  while(valid < N) {
    // inv: a = sum_{k=0}^{valid-1} w^{*k}
    // inv: forall_{i=1}^{valid-1} a[i]=sum_{k=0}^i w[k] a[i-k]

    b0 = inb[0]; inb[0] = b0 + 1;
    conv();  // a=a*(b+delta)
    inb[0] = b0;

    autoconv();  // b=b*b
    valid *= 2;

    if(valid > 2)
      fix_slope(&sla, ina, 1, valid-1);
    new_slope(&slb, inb, sla-slb); slb=sla;

    printf("valid=%i\n", valid);
    printf("slope=%f,%f\n", sla, slb);
    print_2vec(ina, sla, inb, slb);
    print_2vec(ina, 0, inb, 0);

    // if(valid < 2 * N) memset(&ina[valid], 0, sizeof(ina[0]) * (2 * N - valid));
    // if(valid < 2 * N) memset(&inb[valid], 0, sizeof(inb[0]) * (2 * N - valid));
  }

  deinit_conv(128);
}
