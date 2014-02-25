#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

#include <math.h>

#include "convol.h"
#include <fftw3.h>



void print_vec(double *v, double sl) {
  size_t i;

  for(i = 0; i < 2 * N; i++)
    printf("  %3d  %12f  %12f\n", i, v[i] * exp(sl * i), v[i]);
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
    w[i] = (double)1. / i;
    //w[i] = pow((double)i/N, i/10.);
  w[0] = 0;
}

double f_df(double x, double *df, double *w) {
  size_t i;
  double v, dv;

  v = dv = 0;
  for(i = N-1; i >= 1; i--) {
    //printf("[%i:%f; v=%f; x=%f]\n", i,w[i], v,x);
    v = x*v + w[i];
    dv = x*dv + i*w[i];
  }
  //exit(1);

  if(df != NULL) *df = dv;
  return x*v - 1;
}

double asympt_slope(double *w) {
  double cx, cf, cdf;
  double x0;

  if(w[N-1] != 0) x0 = pow(w[1]/w[N-1], 1. / (N-1.));
  else x0 = 1;

  cx=x0; cf = f_df(cx, &cdf, w);
  printf("f(%f)=%f\n", cx, cf);

  while(cf < 0) {
    cx += 0.25 * x0;
    cf = f_df(cx, &cdf, w);
    printf("f(%f)=%f\n", cx, cf);
  }

  while(fabs(cf) >= 1e-2 / N * cdf) {
    cx -= cf/cdf;
    cf = f_df(cx, &cdf, w);
    printf("f(%f)=%f\n", cx, cf);
  }

  return -log(cx);
}

double sla=0, slb=0;

void new_slope(double *sl, double *v, double new_sl) {
  size_t i;

  *sl += new_sl;
  for(i=0; i<N; i++)
    v[i] = v[i] * exp(-new_sl*i);
}

void fix_slope(double *sl, double *v, size_t nl, size_t nr) {
  double a = (log(v[nr]) - log(v[nl]));
  double b = nr-nl;
  printf("[%i %i %f %f %f]\n", nr-nl, nr, a , b, a / b);
  new_slope(sl, v, (log(v[nr]) - log(v[nl])) / (nr-nl));
}

int main() {
  size_t i, valid;
  double b0;

  double *savea;

  init_conv(1024*1024);
  savea = fftw_malloc(sizeof(double) * 2 * N);

  init_weights(N, b);  // b = w
  memset(a, 0, sizeof(a[0]) * 2 * N); a[0] = 1;  // a = delta
  valid = 1;


  new_slope(&sla, a, asympt_slope(b));
  new_slope(&slb, b, sla-slb); slb=sla;

  while(valid < N) {
    // inv: a = sum_{k=0}^{valid-1} w^{*k}
    // inv: forall_{i=1}^{valid-1} a[i]=sum_{k=0}^i w[k] a[i-k]

    // a=a*b+a
    memcpy(savea, a, sizeof(savea[0]) * N); // We only need a[i] for i<N
    conv();  // a=a*b
    for(i = 0; i < N; i++)  // a=a+a_old
      a[i] += savea[i];

    //printf("[valid=%i]\n", valid);
    //print_2vec(a, 0, b, 0);

    //printf("[valid=%i] before\n", valid);
    //print_vec(b, slb);

    // b=b*b
    autoconv();

    //printf("[valid=%i] after\n", valid);
    //print_vec(b, slb);

    valid *= 2;

    /*if(valid > 2)
      fix_slope(&sla, a, 1, valid-1);
    else {
      new_slope(&sla, a, log(2));
    }
    new_slope(&slb, b, sla-slb); slb=sla;*/

    printf("valid=%i\n", valid);
    printf("slope=%f,%f\n", sla, slb);
    //print_2vec(a, 0, b, 0);
    printf("N=%i\nlog_2(a_{N-1}) = %.18f\n", valid, (log(a[valid-1])+sla*(valid-1))/log(2.));

    // if(valid < 2 * N) memset(&a[valid], 0, sizeof(a[0]) * (2 * N - valid));
    // if(valid < 2 * N) memset(&b[valid], 0, sizeof(b[0]) * (2 * N - valid));
  }

    //print_2vec(a, sla, b, slb);
    //print_2vec(a, 0, b, 0);
  printf("N=%i\nlog_2(a_{N-1}) = %.18f\n", N, (log(a[N-1])+sla*(N-1))/log(2.));

  deinit_conv();
  fftw_free(savea);


  return 0;
}
