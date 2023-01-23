//
// Created by 欧阳康的macbook on 2022/2/10.
//
#ifndef LPC_PITCH_H
#define LPC_PITCH_H

#include "Platform.h"


void lpc_analyze(float *_lpc, const float *ac, int p);

int lpc_autocorr(const float *x, float *ac,\
                   const float *window, int overlap, int lag, int n);


void pitch_residual(float *x[], float *x_lp,\
                      int len, int C);

void pitch_search(const float *x_lp, float *y,\
                  int len, int max_pitch, int *pitch);

float remove_doubling(float *x, int maxperiod, int minperiod,\
                           int N, int *T0, int prev_period, float prev_gain);

void pitch_xcorr(const float *_x, const float *_y,\
                      float *xcorr, int len, int max_pitch);


#endif /* LPC_PITCH_H */


