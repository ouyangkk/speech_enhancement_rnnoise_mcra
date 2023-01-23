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
#include "DSP_denoise.h"

/* */
complex DSP_X[FREQ_SIZE];

float denoise_process(DenoiseState *st, float *out, const float *in, int Block) {
    int i;
    complex X[FREQ_SIZE];
    complex P[WINDOW_SIZE];
    

    float x[FRAME_SIZE];
    float Ex[NB_BANDS], Ep[NB_BANDS];
    float Exp[NB_BANDS];
    float features[NB_FEATURES];
    float g[NB_BANDS];
    float gf[FREQ_SIZE]={1};
    float vad_prob = 0;
    int silence;
    static const float a_hp[2] = {-1.99599, 0.99600};
    static const float b_hp[2] = {-2, 1};
    biquad(x, st->mem_hp_x, in, b_hp, a_hp, FRAME_SIZE);

#if RNN_ONLY
    silence = compute_frame_features(st, X, P, Ex, Ep, Exp, features, x);
    if (!silence) {
        compute_rnn(&st->rnn, g, &vad_prob, features);
#if PITCH_F
        pitch_filter(X, P, Ex, Ep, Exp, g);
#endif
        for (i=0;i<NB_BANDS;i++) {
            float alpha = .6f;

            g[i] = MAX(g[i], alpha*st->lastg[i]);
            st->lastg[i] = g[i];

        }
        interp_band_gain(gf, g);
#if 1
        for (i=0;i<FREQ_SIZE;i++) {
            if (i<5 || i>C_F){
                X[i].r *= 0;
                X[i].i *= 0;
            }
            // else{
            X[i].r *= gf[i];
            X[i].i *= gf[i];
        }
#endif
    }

    frame_synthesis(st, out, X);
    return vad_prob;
#endif

#if DSP_ONLY  
    frame_analysis(st, X, Ex, x);
    G_calculate_process(st, X, DSP_X, Block);
    for (i=0;i<FREQ_SIZE;i++)
    {
                X[i].r = DSP_X[i].r;
                X[i].i = DSP_X[i].i;
                //X[i].r = X[i].r ; //for debug
                //X[i].r = X[i].r ;
            }
    frame_synthesis(st, out, X);
    return 0;

#endif

#if DSP_RNN

    silence = compute_frame_features(st, X, P, Ex, Ep, Exp, features, x);
    G_calculate_process(st, X, DSP_X, Block);

    if (!silence) {
        compute_rnn(&st->rnn, g, &vad_prob, features);
#if PITCH_F
        pitch_filter(X, P, Ex, Ep, Exp, g);
#endif
        for (i=0;i<NB_BANDS;i++) {
            float alpha = .6f;

            g[i] = MAX(g[i], (alpha*st->lastg[i]));
            st->lastg[i] = g[i];

        }
        interp_band_gain(gf, g);
#if 1
        for (i=0;i<FREQ_SIZE;i++) {
            if (i<5 || i>C_F){
                X[i].r *= 0;
                X[i].i *= 0;
            }
            else{
            X[i].r *= gf[i];
            X[i].i *= gf[i];
            X[i].r  = WEIGHT_FACTOR*X[i].r + (1-WEIGHT_FACTOR)*DSP_X[i].r;
            X[i].i  = WEIGHT_FACTOR*X[i].i + (1-WEIGHT_FACTOR)*DSP_X[i].i;
        }
  }
#endif


    }

    //X[i].r = DSP_X[i].r;
    //X[i].i = DSP_X[i].i;

    frame_synthesis(st, out, X);
    return vad_prob;
#endif
    }



int main(int argc, char **argv) {
  int Block = 0;
  int i;
  int first = 1;
  float x[FRAME_SIZE];
  FILE *f1, *fout;
  DenoiseState *st;
  st = model_create();
  if (argc!=3) {
    fprintf(stderr, "usage: %s <noisy speech> <output denoised>\n", argv[0]);
    return 1;
  }
  f1 = fopen(argv[1], "rb");
  fout = fopen(argv[2], "wb");
  while (1) {
    short tmp[FRAME_SIZE];
    fread(tmp, sizeof(short), FRAME_SIZE, f1);
    if (feof(f1)) break;
    for (i=0;i<FRAME_SIZE;i++) x[i] = tmp[i];
    denoise_process(st, x, x, Block);
    for (i=0;i<FRAME_SIZE;i++) tmp[i] = x[i];
    if (!first) fwrite(tmp, sizeof(short), FRAME_SIZE, fout);
    first = 0;
    Block = Block + 1 ; 
    }
  mem_destroy(st);
  fclose(f1);
  fclose(fout);
  return 0;
}