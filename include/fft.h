#ifndef FFT_H
#define FFT_H

#include "Platform.h"
#include <stdlib.h>
#include <math.h>

#if defined(ARM_x)
#include "arm/fft_arm.h"
#endif

#define FFT_COS(phase) (float) cos(phase)
#define FFT_SIN(phase) (float) sin(phase)

//fast Fourier transform gene
#define MAXFACTORS 8

typedef struct fft_state{
   int is_supported;
   void *priv;
}fft_state;

typedef struct {
    float r;
    float i;
}complex;


typedef struct _fft_state{
    int nfft;
    float scale;
    int shift;
    int16 factors[2*MAXFACTORS];
    const int16 *bitrev;
    const complex *twiddles;
    fft_state *fft_sup;
} _fft_state;

//The following are related operations for complex numbers

#define  C_ADD( res, a,b)\
    do { \
            CHECK_OVERFLOW_OP((a).r,+,(b).r)\
            CHECK_OVERFLOW_OP((a).i,+,(b).i)\
            (res).r=(a).r+(b).r;  (res).i=(a).i+(b).i; \
    }while(0)

#define  C_SUB( res, a,b)\
    do { \
            CHECK_OVERFLOW_OP((a).r,-,(b).r)\
            CHECK_OVERFLOW_OP((a).i,-,(b).i)\
            (res).r=(a).r-(b).r;  (res).i=(a).i-(b).i; \
    }while(0)

#define C_ADDTO( res , a)\
    do { \
            CHECK_OVERFLOW_OP((res).r,+,(a).r)\
            CHECK_OVERFLOW_OP((res).i,+,(a).i)\
            (res).r += (a).r;  (res).i += (a).i;\
    }while(0)

//#define C_SUBFROM( res , a)\
//    do {\
//            CHECK_OVERFLOW_OP((res).r,-,(a).r)\
//            CHECK_OVERFLOW_OP((res).i,-,(a).i)\
//            (res).r -= (a).r;  (res).i -= (a).i; \
//    }while(0)

#define C_MUL(m,a,b) \
    do{ (m).r = (a).r*(b).r - (a).i*(b).i;\
        (m).i = (a).r*(b).i + (a).i*(b).r; }while(0)

#define C_MULBYSCALAR( c, s ) \
    do{ (c).r *= (s);\
        (c).i *= (s); }while(0) 
                              
_fft_state *fft_alloc_twiddles(int nfft,void * mem,size_t * lenmem,\
                                        const _fft_state *base, int arch);


void fft_c(const _fft_state *st,const complex *fin,complex *fout);
//void ifft_c(const _fft_state *st,const complex *fin,complex *fout);

#define  kf_cexp(x,phase) \
        do{ \
                (x)->r = FFT_COS(phase);\
                (x)->i = FFT_SIN(phase);\
        }while(0)

int fft_alloc_arch(_fft_state *st);

_fft_state *fft_alloc_twiddles(int nfft,void * mem,size_t * lenmem,\
                               const _fft_state *base, int arch);
_fft_state *fft_alloc(int nfft,void * mem,size_t * lenmem, int arch);

void fft_free_arch(_fft_state *st);
void fft_free(const _fft_state *cfg, int arch);
void fft_main(const _fft_state *st,complex *fout);
void fft_c(const _fft_state *st,const complex *fin,complex *fout);



#endif