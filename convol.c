#include "convol.h"

void uninit_conv();

void init_conv(size_t _N) {
  if(N != 0) deinit_conv();
  N = _N;

  ina = fftw_malloc(sizeof(double) * 2 * N);
  inb = fftw_malloc(sizeof(double) * 2 * N);
  outa = fftw_malloc(sizeof(fftw_complex) * (N + 1));
  outb = fftw_malloc(sizeof(fftw_complex) * (N + 1));

  plan_fw_a = fftw_plan_dft_r2c_1d(2 * N, ina, outa, FFTW_MEASURE);
  plan_bw_a = fftw_plan_dft_c2r_1d(2 * N, outa, ina, FFTW_MEASURE);
  plan_fw_b = fftw_plan_dft_r2c_1d(2 * N, inb, outb, FFTW_ESTIMATE);
  plan_bw_b = fftw_plan_dft_c2r_1d(2 * N, outb, inb, FFTW_ESTIMATE);
}

void deinit_conv() {
  fftw_destroy_plan(plan_fw_a);
  fftw_destroy_plan(plan_bw_a);
  fftw_destroy_plan(plan_fw_b);
  fftw_destroy_plan(plan_bw_b);

  fftw_free(ina);
  fftw_free(inb);
  fftw_free(outa);
  fftw_free(outb);

  N = 0;
}

void conv() {  // a[n:2n-1]=b[n:2n-1]=0; a=a*b;
  size_t i;

  memset(&(ina[N]), 0, sizeof(ina[0]) * N);
  memset(&(inb[N]), 0, sizeof(inb[0]) * N);

  fftw_execute(plan_fw_a);
  fftw_execute(plan_fw_b);

  for(i=0; i<=N; i++) {
    double a0 = outa[i][0];
    double a1 = outa[i][1];

    outa[i][0] = (a0 * outb[i][0] - a1 * outb[i][1]) / (2.0 * N);
    outa[i][1] = (a0 * outb[i][1] + a1 * outb[i][0]) / (2.0 * N);
  }

  fftw_execute(plan_bw_a);
}

void autoconv() {  // b[n:2n-1]=0; b=b*b;
  size_t i;

  memset(&(inb[N]), 0, sizeof(inb[0]) * N);

  fftw_execute(plan_fw_b);

  for(i=0; i<=N; i++) {
    double b0 = outb[i][0];
    double b1 = outb[i][1];

    outb[i][0] = (b0 * b0 - b1 * b1) / (2.0 * N);
    outb[i][1] = (2 * b0 * b1      ) / (2.0 * N);
  }

  fftw_execute(plan_bw_b);
}
