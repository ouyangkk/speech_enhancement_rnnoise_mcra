//
// Created by 欧阳康的macbook on 2022/2/12.
//
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include "fft.h"
#include "define.h"
#include "Platform.h"
#include "lpc_pitch.h"
#include "rnn.h"
#include "rnn_weight.h"
#include "denoise_function.h"
static float uni_rand() {
    return rand()/(double)RAND_MAX-.5;
}

extern int eband5ms[];
static void rand_resp(float *a, float *b) {
    a[0] = .75*uni_rand();
    a[1] = .75*uni_rand();
    b[0] = .75*uni_rand();
    b[1] = .75*uni_rand();
}

int lowpass = FREQ_SIZE;
int band_lp = NB_BANDS;
int main(int argc, char **argv) {
    int i;
    int count=0;
    static const float a_hp[2] = {-1.99599, 0.99600};
    static const float b_hp[2] = {-2, 1};
    float a_noise[2] = {0};
    float b_noise[2] = {0};
    float a_sig[2] = {0};
    float b_sig[2] = {0};
    float mem_hp_x[2]={0};
    float mem_hp_n[2]={0};
    float mem_resp_x[2]={0};
    float mem_resp_n[2]={0};
    float x[FRAME_SIZE];
    float n[FRAME_SIZE];
    float xn[FRAME_SIZE];
    int vad_cnt=0;
    int gain_change_count=0;
    float speech_gain = 1, noise_gain = 1;
    FILE *f1, *f2;
    int maxCount;
    DenoiseState *st;
    DenoiseState *noise_state;
    DenoiseState *noisy;
    st = model_create();
    noise_state = model_create();
    noisy = model_create();
    if (argc!=4) {
        fprintf(stderr, "usage: %s <speech> <noise> <count>\n", argv[0]);
        return 1;
    }
    f1 = fopen(argv[1], "r");
    f2 = fopen(argv[2], "r");
    maxCount = 12000;
    for(i=0;i<150;i++) {
        short tmp[FRAME_SIZE];
        fread(tmp, sizeof(short), FRAME_SIZE, f2);
    }
    while (1) {
        complex X[FREQ_SIZE], Y[FREQ_SIZE], N[FREQ_SIZE], P[WINDOW_SIZE];
        float Ex[NB_BANDS], Ey[NB_BANDS], En[NB_BANDS], Ep[NB_BANDS];
        float Exp[NB_BANDS];
        float Ln[NB_BANDS];
        float features[NB_FEATURES];
        float g[NB_BANDS];
        short tmp[FRAME_SIZE];
        float vad=0;
        float E=0;
        if (count==maxCount) break;
        if ((count%1000)==0) fprintf(stderr, "%d\r", count);

#if DATA_RANDOM
        if (++gain_change_count > 2821) {
            speech_gain = pow(10., (-40+(rand()%60))/20.);
            noise_gain = pow(10., (-30+(rand()%50))/20.);
            if (rand()%10==0) noise_gain = 0;
            noise_gain *= speech_gain;
            if (rand()%10==0) speech_gain = 0;
            gain_change_count = 0;
            rand_resp(a_noise, b_noise);
            rand_resp(a_sig, b_sig);
            lowpass = FREQ_SIZE * 4000./8000. * pow(50., rand()/(double)RAND_MAX);
            for (i=0;i<NB_BANDS;i++) {
                if (eband5ms[i]<<SHIFT > lowpass) {
                    band_lp = i;
                    break;
                }
            }
        }
        if (speech_gain != 0) {
            fread(tmp, sizeof(short), FRAME_SIZE, f1);
            if (feof(f1)) {
                rewind(f1);
                fread(tmp, sizeof(short), FRAME_SIZE, f1);
            }
            for (i=0;i<FRAME_SIZE;i++) x[i] = speech_gain*tmp[i];
            for (i=0;i<FRAME_SIZE;i++) E += tmp[i]*(float)tmp[i];
        } else {
            for (i=0;i<FRAME_SIZE;i++) x[i] = 0;
            E = 0;
        }
        if (noise_gain!=0) {
            fread(tmp, sizeof(short), FRAME_SIZE, f2);
            if (feof(f2)) {
                rewind(f2);
                fread(tmp, sizeof(short), FRAME_SIZE, f2);
            }
            for (i=0;i<FRAME_SIZE;i++) n[i] = noise_gain*tmp[i];
        } else {
            for (i=0;i<FRAME_SIZE;i++) n[i] = 0;
        }
        biquad(x, mem_hp_x, x, b_hp, a_hp, FRAME_SIZE);
        biquad(x, mem_resp_x, x, b_sig, a_sig, FRAME_SIZE);
        biquad(n, mem_hp_n, n, b_hp, a_hp, FRAME_SIZE);
        biquad(n, mem_resp_n, n, b_noise, a_noise, FRAME_SIZE);

# endif

#if DATA_NOTRANDOM
        fread(tmp, sizeof(short), FRAME_SIZE, f1);
    if (feof(f1)) {
        rewind(f1);
        fread(tmp, sizeof(short), FRAME_SIZE, f1);
    }
    for (i=0;i<FRAME_SIZE;i++) x[i] = tmp[i];
    for (i=0;i<FRAME_SIZE;i++) E += tmp[i]*(float)tmp[i];

fread(tmp, sizeof(short), FRAME_SIZE, f2);
if (feof(f2)) {
    rewind(f2);
    fread(tmp, sizeof(short), FRAME_SIZE, f2);
}
for (i=0;i<FRAME_SIZE;i++) n[i] = tmp[i];

#endif

        for (i=0;i<FRAME_SIZE;i++) xn[i] = x[i] + n[i];
        if (E > 1e9f) {
            vad_cnt=0;
        } else if (E > 1e8f) {
            vad_cnt -= 5;
        } else if (E > 1e7f) {
            vad_cnt++;
        } else {
            vad_cnt+=2;
        }
        if (vad_cnt < 0) vad_cnt = 0;
        if (vad_cnt > 15) vad_cnt = 15;

        if (vad_cnt >= 10) vad = 0;
        else if (vad_cnt > 0) vad = 0.5f;
        else vad = 1.f;

        // if (E > 1e-1)
        // {
        //   vad = 1;
        // }
        // else
        // {
        //   vad = 0;
        // }


        frame_analysis(st, Y, Ey, x);
        frame_analysis(noise_state, N, En, n);
        for (i=0;i<NB_BANDS;i++) Ln[i] = log10(1e-2+En[i]);
        int silence = compute_frame_features(noisy, X, P, Ex, Ep, Exp, features, xn);
#if PITCH_F
        pitch_filter(X, P, Ex, Ep, Exp, g);
#endif
        for (i=0;i<NB_BANDS;i++) {
            g[i] = sqrt((Ey[i]+1e-6)/(Ex[i]+1e-6));
            if (g[i] > 1) g[i] = 1;
            if (silence || i > band_lp) g[i] = -1;
            if (Ey[i] < 5e-3 && Ex[i] < 5e-3) g[i] = -1;
            if (vad==0 && noise_gain==0) g[i] = -1;
        }
        count++;
#if 1
        fwrite(features, sizeof(float), NB_FEATURES, stdout);
        fwrite(g, sizeof(float), NB_BANDS, stdout);
        fwrite(Ln, sizeof(float), NB_BANDS, stdout);
        fwrite(&vad, sizeof(float), 1, stdout);
#endif
    }
    fprintf(stderr, "matrix size: %d x %d\n", count, NB_FEATURES + 2*NB_BANDS + 1);
    fclose(f1);
    fclose(f2);
    return 0;
}

