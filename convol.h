#ifndef _CONVOL_H_
#define _CONVOL_H_

#include <string.h>
#include <fftw3.h>

size_t N;
double *ina;
double *inb;
fftw_complex *outa;
fftw_complex *outb;
fftw_plan plan_fw_a, plan_fw_b, plan_bw_a, plan_bw_b;

void init_conv(size_t _N);
void deinit_conv();
void conv();  // a[n:2n-1]=b[n:2n-1]=0; a=a*b; 
void autoconv();  // b[n:2n-1]=0; b=b*b;

#endif  // _CONVOL_H_