//
// Created by 欧阳康的macbook on 2022/2/12.
//

#ifndef OM_DEBUG_DENOISE_FUNCTION_H
#define OM_DEBUG_DENOISE_FUNCTION_H
#include "fft.h"
#include "rnn_weight.h"
struct DenoiseState {
    float analysis_mem[FRAME_SIZE];
    float cepstral_mem[CEPS_MEM][NB_BANDS];
    int memid;
    float synthesis_mem[FRAME_SIZE];
    float pitch_buf[PITCH_BUF_SIZE];
    float pitch_enh_buf[PITCH_BUF_SIZE];
    float last_gain;
    int last_period;
    float mem_hp_x[2];
    float lastg[NB_BANDS];
    RNNState rnn;
};



typedef struct DenoiseState DenoiseState;
void compute_band_energy(float *bandE, const complex *X);
void compute_band_corr(float *bandE, const complex *X, const complex *P);
void interp_band_gain(float *g, const float *bandE);
void check_init();
void dct(float *out, const float *in);
void forward_transform(complex *out, const float *in);
void inverse_transform(float *out, const complex *in);
void apply_window(float *x);
int get_size();
int model_init(DenoiseState *st);
DenoiseState *model_create();
void mem_destroy(DenoiseState *st);
void frame_analysis(DenoiseState *st, complex *X, float *Ex, const float *in);
int compute_frame_features(DenoiseState *st, complex *X, complex *P,
                           float *Ex, float *Ep, float *Exp, float *features, const float *in);
void frame_synthesis(DenoiseState *st, float *out, const complex *y);
void biquad(float *y, float mem[2], const float *x, const float *b, const float *a, int N);
void pitch_filter(complex *X, const complex *P, const float *Ex, const float *Ep,
                  const float *Exp, const float *g);
#endif //OM_DEBUG_DENOISE_FUNCTION_H
