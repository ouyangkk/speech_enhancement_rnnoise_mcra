#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include "fft.h"
#include "lpc_pitch.h"
#include "define.h"
#include "Platform.h"

#include "rnn.h"
#include "rnn_weight.h"
#include "denoise_function.h"
typedef struct {
    int init;
    _fft_state *kfft;
    float half_window[FRAME_SIZE];
    float dct_table[NB_BANDS*NB_BANDS];
} CommonState;
CommonState common;



int eband5ms[] = {
        0,  1,  2,  3,  4,  5,  6,  7,  8, 10, 12, 14, 16, 20, 24, 28, 34, 40
};

void compute_band_energy(float *bandE, const complex *X) {
    int i;
    float sum[NB_BANDS] = {0};
    for (i=0;i<NB_BANDS-1;i++)
    {
        int j;
        int band_size;
        band_size = (eband5ms[i+1]-eband5ms[i])<<SHIFT;
        for (j=0;j<band_size;j++) {
            float tmp;
            float frac = (float)j/band_size;
            tmp = SQUARE(X[(eband5ms[i]<<SHIFT) + j].r);
            tmp += SQUARE(X[(eband5ms[i]<<SHIFT) + j].i);
            sum[i] += (1-frac)*tmp;
            sum[i+1] += frac*tmp;
        }
    }
    sum[0] *= 2;
    sum[NB_BANDS-1] *= 2;
    for (i=0;i<NB_BANDS;i++)
    {
        bandE[i] = sum[i];
    }
}

void compute_band_corr(float *bandE, const complex *X, const complex *P) {
    int i;
    float sum[NB_BANDS] = {0};
    for (i=0;i<NB_BANDS-1;i++)
    {
        int j;
        int band_size;
        band_size = (eband5ms[i+1]-eband5ms[i])<<SHIFT;
        for (j=0;j<band_size;j++) {
            float tmp;
            float frac = (float)j/band_size;
            tmp = X[(eband5ms[i]<<SHIFT) + j].r * P[(eband5ms[i]<<SHIFT) + j].r;
            tmp += X[(eband5ms[i]<<SHIFT) + j].i * P[(eband5ms[i]<<SHIFT) + j].i;
            sum[i] += (1-frac)*tmp;
            sum[i+1] += frac*tmp;
        }
    }
    sum[0] *= 2;
    sum[NB_BANDS-1] *= 2;
    for (i=0;i<NB_BANDS;i++)
    {
        bandE[i] = sum[i];
    }
}

void interp_band_gain(float *g, const float *bandE) {
    int i;
    memset(g, 0, FREQ_SIZE);
    for (i=0;i<NB_BANDS-1;i++)
    {
        int j;
        int band_size;
        band_size = (eband5ms[i+1]-eband5ms[i])<<SHIFT;
        for (j=0;j<band_size;j++) {
            float frac = (float)j/band_size;
            g[(eband5ms[i]<<SHIFT) + j] = (1-frac)*bandE[i] + frac*bandE[i+1];
        }
    }
}


void check_init() {
    int i;
    if (common.init) return;
    common.kfft = fft_alloc_twiddles(2*FRAME_SIZE, NULL, NULL, NULL, 0);
    for (i=0;i<FRAME_SIZE;i++)
        common.half_window[i] = sin(.5*PI*sin(.5*PI*(i+.5)/FRAME_SIZE) * sin(.5*PI*(i+.5)/FRAME_SIZE));
    for (i=0;i<NB_BANDS;i++) {
        int j;
        for (j=0;j<NB_BANDS;j++) {
            common.dct_table[i*NB_BANDS + j] = cos((i+.5)*j*PI/NB_BANDS);
            if (j==0) common.dct_table[i*NB_BANDS + j] *= sqrt(.5);
        }
    }
    common.init = 1;
}

void dct(float *out, const float *in) {
    int i;
    check_init();
    for (i=0;i<NB_BANDS;i++) {
        int j;
        float sum = 0;
        for (j=0;j<NB_BANDS;j++) {
            sum += in[j] * common.dct_table[j*NB_BANDS + i];
        }
        out[i] = sum*sqrt(2./22);
    }
}



void forward_transform(complex *out, const float *in) {
    int i;
    complex x[WINDOW_SIZE];
    complex y[WINDOW_SIZE];
    check_init();
    for (i=0;i<WINDOW_SIZE;i++) {
        x[i].r = in[i] * 0.8;
        x[i].i = 0;
    }
    fft_c(common.kfft, x, y);
    for (i=0;i<FREQ_SIZE;i++) {
        out[i] = y[i];
    }
}

void inverse_transform(float *out, const complex *in) {
    int i;
    complex x[WINDOW_SIZE];
    complex y[WINDOW_SIZE];
    check_init();
    for (i=0;i<FREQ_SIZE;i++) {
        x[i] = in[i];
    }
    for (;i<WINDOW_SIZE;i++) {
        x[i].r = x[WINDOW_SIZE - i].r;
        x[i].i = -x[WINDOW_SIZE - i].i;
    }
    fft_c(common.kfft, x, y);
    /* output in reverse order for IFFT. */
    out[0] = WINDOW_SIZE*y[0].r;
    for (i=1;i<WINDOW_SIZE;i++) {
        out[i] = WINDOW_SIZE*y[WINDOW_SIZE - i].r;
    }
}

void apply_window(float *x) {
    int i;
    check_init();
    for (i=0;i<FRAME_SIZE;i++) {
        x[i] *= common.half_window[i];
        x[WINDOW_SIZE - 1 - i] *= common.half_window[i];
    }
}

int get_size() {
    return sizeof(DenoiseState);
}

int model_init(DenoiseState *st) {
    memset(st, 0, sizeof(*st));
    return 0;
}

DenoiseState *model_create() {
    DenoiseState *st;
    st = malloc(get_size());
    model_init(st);
    return st;
}

void mem_destroy(DenoiseState *st) {
    free(st);
}

void frame_analysis(DenoiseState *st, complex *X, float *Ex, const float *in) {
    int i;
    float x[WINDOW_SIZE];
    RNN_COPY(x, st->analysis_mem, FRAME_SIZE);
    for (i=0;i<FRAME_SIZE;i++) x[FRAME_SIZE + i] = in[i];
    RNN_COPY(st->analysis_mem, in, FRAME_SIZE);
    apply_window(x);
    forward_transform(X, x);
    compute_band_energy(Ex, X);
}

int compute_frame_features(DenoiseState *st, complex *X, complex *P,
                                  float *Ex, float *Ep, float *Exp, float *features, const float *in) {
    int i;
    float E = 0;
    float *ceps_0, *ceps_1, *ceps_2;
    float spec_variability = 0;
    float Ly[NB_BANDS];
    float p[WINDOW_SIZE];
    float pitch_buf[PITCH_BUF_SIZE>>1];
    int pitch_index;
    float gain;
    float *(pre[1]);
    float tmp[NB_BANDS];
    float follow, logMax;
    frame_analysis(st, X, Ex, in);
    RNN_MOVE(st->pitch_buf, &st->pitch_buf[FRAME_SIZE], PITCH_BUF_SIZE-FRAME_SIZE);
    RNN_COPY(&st->pitch_buf[PITCH_BUF_SIZE-FRAME_SIZE], in, FRAME_SIZE);
    pre[0] = &st->pitch_buf[0];
    pitch_residual(pre, pitch_buf, PITCH_BUF_SIZE, 1);
    pitch_search(pitch_buf+(PITCH_MAX_PERIOD>>1), pitch_buf, PITCH_FRAME_SIZE,
                 PITCH_MAX_PERIOD-3*PITCH_MIN_PERIOD, &pitch_index);
    pitch_index = PITCH_MAX_PERIOD-pitch_index;

    gain = remove_doubling(pitch_buf, PITCH_MAX_PERIOD, PITCH_MIN_PERIOD,
                           PITCH_FRAME_SIZE, &pitch_index, st->last_period, st->last_gain);
    st->last_period = pitch_index;
    st->last_gain = gain;
    for (i=0;i<WINDOW_SIZE;i++)
        p[i] = st->pitch_buf[PITCH_BUF_SIZE-WINDOW_SIZE-pitch_index+i];
    apply_window(p);
    forward_transform(P, p);
    compute_band_energy(Ep, P);
    compute_band_corr(Exp, X, P);
    for (i=0;i<NB_BANDS;i++) Exp[i] = Exp[i]/sqrt(.001+Ex[i]*Ep[i]);
    dct(tmp, Exp);
    for (i=0;i<NB_DELTA_CEPS;i++) features[NB_BANDS+2*NB_DELTA_CEPS+i] = tmp[i];
    features[NB_BANDS+2*NB_DELTA_CEPS] -= 1.3;
    features[NB_BANDS+2*NB_DELTA_CEPS+1] -= 0.9;
    features[NB_BANDS+3*NB_DELTA_CEPS] = .01*(pitch_index-100);
    logMax = -2;
    follow = -2;
    for (i=0;i<NB_BANDS;i++) {
        Ly[i] = log10(1e-2+Ex[i]);
        Ly[i] = MAX(logMax-7, MAX(follow-1.5, Ly[i]));
        logMax = MAX(logMax, Ly[i]);
        follow = MAX(follow-1.5, Ly[i]);
        E += Ex[i];
    }
    if (!TRAINING && E < 0.03) {
        /* If there's no audio, avoid messing up the state. */
        RNN_CLEAR(features, NB_FEATURES);
        return 1;
    }
    dct(features, Ly);
    features[0] -= 12;
    features[1] -= 4;
    ceps_0 = st->cepstral_mem[st->memid];
    ceps_1 = (st->memid < 1) ? st->cepstral_mem[CEPS_MEM+st->memid-1] : st->cepstral_mem[st->memid-1];
    ceps_2 = (st->memid < 2) ? st->cepstral_mem[CEPS_MEM+st->memid-2] : st->cepstral_mem[st->memid-2];
    for (i=0;i<NB_BANDS;i++) ceps_0[i] = features[i];
    st->memid++;
    for (i=0;i<NB_DELTA_CEPS;i++) {
        features[i] = ceps_0[i] + ceps_1[i] + ceps_2[i];
        features[NB_BANDS+i] = ceps_0[i] - ceps_2[i];
        features[NB_BANDS+NB_DELTA_CEPS+i] =  ceps_0[i] - 2*ceps_1[i] + ceps_2[i];
    }
    /* Spectral variability features. */
    if (st->memid == CEPS_MEM) st->memid = 0;
    for (i=0;i<CEPS_MEM;i++)
    {
        int j;
        float mindist = 1e15f;
        for (j=0;j<CEPS_MEM;j++)
        {
            int k;
            float dist=0;
            for (k=0;k<NB_BANDS;k++)
            {
                float tmp;
                tmp = st->cepstral_mem[i][k] - st->cepstral_mem[j][k];
                dist += tmp*tmp;
            }
            if (j!=i)
                mindist = MIN(mindist, dist);
        }
        spec_variability += mindist;
    }
    features[NB_BANDS+3*NB_DELTA_CEPS+1] = spec_variability/CEPS_MEM-2.1;
    return TRAINING && E < 0.1;
}

void frame_synthesis(DenoiseState *st, float *out, const complex *y) {
    float x[WINDOW_SIZE];
    int i;
    inverse_transform(x, y);
    apply_window(x);
    for (i=0;i<FRAME_SIZE;i++) out[i] = x[i] + st->synthesis_mem[i];
    RNN_COPY(st->synthesis_mem, &x[FRAME_SIZE], FRAME_SIZE);
}

void biquad(float *y, float mem[2], const float *x, const float *b, const float *a, int N) {
    int i;
    for (i=0;i<N;i++) {
        float xi, yi;
        xi = x[i];
        yi = x[i] + mem[0];
        mem[0] = mem[1] + (b[0]*(double)xi - a[0]*(double)yi);
        mem[1] = (b[1]*(double)xi - a[1]*(double)yi);
        y[i] = yi;
    }
}

void pitch_filter(complex *X, const complex *P, const float *Ex, const float *Ep,
                  const float *Exp, const float *g) {
    int i;
    float r[NB_BANDS];
    float rf[FREQ_SIZE] = {0};
    for (i=0;i<NB_BANDS;i++) {
#if 0 
    if (Exp[i]>g[i]) r[i] = 1;
    else r[i] = Exp[i]*(1-g[i])/(.001 + g[i]*(1-Exp[i]));
    r[i] = MIN(1, MAX(0, r[i]));
#else
        if (Exp[i]>g[i]) r[i] = 1;
        else r[i] = SQUARE(Exp[i])*(1-SQUARE(g[i]))/(.001 + SQUARE(g[i])*(1-SQUARE(Exp[i])));
        r[i] = sqrt(MIN(1, MAX(0, r[i])));
#endif
        r[i] *= sqrt(Ex[i]/(1e-8+Ep[i]));
    }
    interp_band_gain(rf, r);
    for (i=0;i<FREQ_SIZE;i++) {
        X[i].r += rf[i]*P[i].r;
        X[i].i += rf[i]*P[i].i;
    }
    float newE[NB_BANDS];
    compute_band_energy(newE, X);
    float norm[NB_BANDS];
    float normf[FREQ_SIZE]={0};
    for (i=0;i<NB_BANDS;i++) {
        norm[i] = sqrt(Ex[i]/(1e-8+newE[i]));
    }
    interp_band_gain(normf, norm);
    for (i=0;i<FREQ_SIZE;i++) {
            X[i].r *= normf[i];
            X[i].i *= normf[i];
        }
    }







