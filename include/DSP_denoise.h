//
// Created by 欧阳康的macbook on 2022/2/11.
//
#ifndef DSP_DENOISE_H
#define DSP_DENOISE_H
#include "define.h"
#include "denoise_function.h"
float amp_min[FREQ_SIZE], amp[FREQ_SIZE], amp_tmp[FREQ_SIZE];
float arr_temp[FREQ_SIZE];
float init_p[FREQ_SIZE];
float post_snr[FREQ_SIZE], current_snr[FREQ_SIZE], pr_snr[FREQ_SIZE];
float v_int[FREQ_SIZE];
int  integra[FREQ_SIZE];
float m_int[FREQ_SIZE];
float gh1[FREQ_SIZE], g[FREQ_SIZE];
float q_est[FREQ_SIZE], p_est[FREQ_SIZE];
float hann_win[WIN_ALL];
float snr[FREQ_SIZE], old_snr[FREQ_SIZE];
float snr_local[FREQ_SIZE], snr_global[FREQ_SIZE];
float p_local[FREQ_SIZE], p_global[FREQ_SIZE];

unsigned int amp_pr[FREQ_SIZE], amp_last[FREQ_SIZE];
float noise_est[FREQ_SIZE];

void mcra_init();
void NoiseEstimation(int blockInd);
void SpeechAbsenceEstm();
float G_calculate_process(DenoiseState *st, complex *X, complex *DSP_X, int Block);


#endif