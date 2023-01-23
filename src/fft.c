#include "fft.h"
#include "Platform.h"

#define CUSTOM_MODES
static void kf_bfly2(
        complex * Fout,
        int m,
        int N
)
{
    complex * Fout2;
    int i;
    (void)m;
#ifdef CUSTOM_MODES
    if (m==1)
    {
        //assert(m==1);
        for (i=0;i<N;i++)
        {
            complex t;
            Fout2 = Fout + 1;
            t = *Fout2;
            C_SUB( *Fout2 ,  *Fout , t );
            C_ADDTO( *Fout ,  t );
            Fout += 2;
        }
    } else
#endif
    {
        float tw;
        //tw = QCONST16(0.7071067812f, 15);
        tw = 0.7071067812f;
        //assert(m==4);
        for (i=0;i<N;i++)
        {
            complex t;
            Fout2 = Fout + 4;
            t = Fout2[0];
            C_SUB( Fout2[0] ,  Fout[0] , t );
            C_ADDTO( Fout[0] ,  t );

            t.r = MUL(ADD(Fout2[1].r, Fout2[1].i), tw);
            t.i = MUL(SUB(Fout2[1].i, Fout2[1].r), tw);
            C_SUB( Fout2[1] ,  Fout[1] , t );
            C_ADDTO( Fout[1] ,  t );

            t.r = Fout2[2].i;
            t.i = -Fout2[2].r;
            C_SUB( Fout2[2] ,  Fout[2] , t );
            C_ADDTO( Fout[2] ,  t );

            t.r = MUL(SUB(Fout2[3].i, Fout2[3].r), tw);
            t.i = MUL(NEG(ADD(Fout2[3].i, Fout2[3].r)), tw);
            C_SUB( Fout2[3] ,  Fout[3] , t );
            C_ADDTO( Fout[3] ,  t );
            Fout += 8;
        }
    }
}

static void kf_bfly4(
        complex * Fout,
        const size_t fstride,
        const _fft_state *st,
        int m,
        int N,
        int mm
)
{
    int i;

    if (m==1)
    {
        /* Degenerate case where all the twiddles are 1. */
        for (i=0;i<N;i++)
        {
            complex scratch0, scratch1;

            C_SUB( scratch0 , *Fout, Fout[2] );
            C_ADDTO(*Fout, Fout[2]);
            C_ADD( scratch1 , Fout[1] , Fout[3] );
            C_SUB( Fout[2], *Fout, scratch1 );
            C_ADDTO( *Fout , scratch1 );
            C_SUB( scratch1 , Fout[1] , Fout[3] );

            Fout[1].r = ADD(scratch0.r, scratch1.i);
            Fout[1].i = SUB(scratch0.i, scratch1.r);
            Fout[3].r = SUB(scratch0.r, scratch1.i);
            Fout[3].i = ADD(scratch0.i, scratch1.r);
            Fout+=4;
        }
    } else {
        int j;
        complex scratch[6];
        const complex *tw1,*tw2,*tw3;
        const int m2=2*m;
        const int m3=3*m;
        complex * Fout_beg = Fout;
        for (i=0;i<N;i++)
        {
            Fout = Fout_beg + i*mm;
            tw3 = tw2 = tw1 = st->twiddles;
            /* m is guaranteed to be a multiple of 4. */
            for (j=0;j<m;j++)
            {
                C_MUL(scratch[0],Fout[m] , *tw1 );
                C_MUL(scratch[1],Fout[m2] , *tw2 );
                C_MUL(scratch[2],Fout[m3] , *tw3 );

                C_SUB( scratch[5] , *Fout, scratch[1] );
                C_ADDTO(*Fout, scratch[1]);
                C_ADD( scratch[3] , scratch[0] , scratch[2] );
                C_SUB( scratch[4] , scratch[0] , scratch[2] );
                C_SUB( Fout[m2], *Fout, scratch[3] );
                tw1 += fstride;
                tw2 += fstride*2;
                tw3 += fstride*3;
                C_ADDTO( *Fout , scratch[3] );

                Fout[m].r = ADD(scratch[5].r, scratch[4].i);
                Fout[m].i = SUB(scratch[5].i, scratch[4].r);
                Fout[m3].r = SUB(scratch[5].r, scratch[4].i);
                Fout[m3].i = ADD(scratch[5].i, scratch[4].r);
                ++Fout;
            }
        }
    }
}


#ifndef RADIX_TWO_ONLY

static void kf_bfly3(
        complex * Fout,
        const size_t fstride,
        const _fft_state *st,
        int m,
        int N,
        int mm
)
{
    int i;
    size_t k;
    const size_t m2 = 2*m;
    const complex *tw1,*tw2;
    complex scratch[5];
    complex epi3;

    complex * Fout_beg = Fout;
#ifdef FIXED_POINT
    /*epi3.r = -16384;*/ /* Unused */
   epi3.i = -28378;
#else
    epi3 = st->twiddles[fstride*m];
#endif
    for (i=0;i<N;i++)
    {
        Fout = Fout_beg + i*mm;
        tw1=tw2=st->twiddles;
        /* For non-custom modes, m is guaranteed to be a multiple of 4. */
        k=m;
        do {

            C_MUL(scratch[1],Fout[m] , *tw1);
            C_MUL(scratch[2],Fout[m2] , *tw2);

            C_ADD(scratch[3],scratch[1],scratch[2]);
            C_SUB(scratch[0],scratch[1],scratch[2]);
            tw1 += fstride;
            tw2 += fstride*2;

            Fout[m].r = SUB(Fout->r, HALF(scratch[3].r));
            Fout[m].i = SUB(Fout->i, HALF(scratch[3].i));

            C_MULBYSCALAR( scratch[0] , epi3.i );

            C_ADDTO(*Fout,scratch[3]);

            Fout[m2].r = ADD(Fout[m].r, scratch[0].i);
            Fout[m2].i = SUB(Fout[m].i, scratch[0].r);

            Fout[m].r = SUB(Fout[m].r, scratch[0].i);
            Fout[m].i = ADD(Fout[m].i, scratch[0].r);

            ++Fout;
        } while(--k);
    }
}


#ifndef OVERRIDE_kf_bfly5
static void kf_bfly5(
        complex * Fout,
        const size_t fstride,
        const _fft_state *st,
        int m,
        int N,
        int mm
)
{
    complex *Fout0,*Fout1,*Fout2,*Fout3,*Fout4;
    int i, u;
    complex scratch[13];
    const complex *tw;
    complex ya,yb;
    complex * Fout_beg = Fout;

#ifdef FIXED_POINT 
   ya.r = 10126;
   ya.i = -31164;
   yb.r = -26510;
   yb.i = -19261;
#else
    ya = st->twiddles[fstride*m];
    yb = st->twiddles[fstride*2*m];
#endif
    tw=st->twiddles;

    for (i=0;i<N;i++)
    {
        Fout = Fout_beg + i*mm;
        Fout0=Fout;
        Fout1=Fout0+m;
        Fout2=Fout0+2*m;
        Fout3=Fout0+3*m;
        Fout4=Fout0+4*m;

        /* For non-custom modes, m is guaranteed to be a multiple of 4. */
        for ( u=0; u<m; ++u ) {
            scratch[0] = *Fout0;

            C_MUL(scratch[1] ,*Fout1, tw[u*fstride]);
            C_MUL(scratch[2] ,*Fout2, tw[2*u*fstride]);
            C_MUL(scratch[3] ,*Fout3, tw[3*u*fstride]);
            C_MUL(scratch[4] ,*Fout4, tw[4*u*fstride]);

            C_ADD( scratch[7],scratch[1],scratch[4]);
            C_SUB( scratch[10],scratch[1],scratch[4]);
            C_ADD( scratch[8],scratch[2],scratch[3]);
            C_SUB( scratch[9],scratch[2],scratch[3]);

            Fout0->r = ADD(Fout0->r, ADD(scratch[7].r, scratch[8].r));
            Fout0->i = ADD(Fout0->i, ADD(scratch[7].i, scratch[8].i));

            scratch[5].r = ADD(scratch[0].r, ADD(MUL(scratch[7].r,ya.r), MUL(scratch[8].r,yb.r)));
            scratch[5].i = ADD(scratch[0].i, ADD(MUL(scratch[7].i,ya.r), MUL(scratch[8].i,yb.r)));

            scratch[6].r =  ADD(MUL(scratch[10].i,ya.i), MUL(scratch[9].i,yb.i));
            scratch[6].i = NEG(ADD(MUL(scratch[10].r,ya.i), MUL(scratch[9].r,yb.i)));

            C_SUB(*Fout1,scratch[5],scratch[6]);
            C_ADD(*Fout4,scratch[5],scratch[6]);

            scratch[11].r = ADD(scratch[0].r, ADD(MUL(scratch[7].r,yb.r), MUL(scratch[8].r,ya.r)));
            scratch[11].i = ADD(scratch[0].i, ADD(MUL(scratch[7].i,yb.r), MUL(scratch[8].i,ya.r)));
            scratch[12].r = SUB(MUL(scratch[9].i,ya.i), MUL(scratch[10].i,yb.i));
            scratch[12].i = SUB(MUL(scratch[10].r,yb.i), MUL(scratch[9].r,ya.i));

            C_ADD(*Fout2,scratch[11],scratch[12]);
            C_SUB(*Fout3,scratch[11],scratch[12]);

            ++Fout0;++Fout1;++Fout2;++Fout3;++Fout4;
        }
    }
}
#endif /* OVERRIDE_kf_bfly5 */


#endif


#ifdef CUSTOM_MODES

static
void compute_bitrev_table(
        int Fout,
        int16 *f,
        const size_t fstride,
        int in_stride,
        int16 * factors,
        const _fft_state *st
)
{
    const int p=*factors++; /* the radix  */
    const int m=*factors++; /* stage's fft length/p */

    /*printf ("fft %d %d %d %d %d %d\n", p*m, m, p, s2, fstride*in_stride, N);*/
    if (m==1)
    {
        int j;
        for (j=0;j<p;j++)
        {
            *f = Fout+j;
            f += fstride*in_stride;
        }
    } else {
        int j;
        for (j=0;j<p;j++)
        {
            compute_bitrev_table( Fout , f, fstride*p, in_stride, factors,st);
            f += fstride*in_stride;
            Fout += m;
        }
    }
}

/*  facbuf is populated by p1,m1,p2,m2, ...
    where
    p[i] * m[i] = m[i-1]
    m0 = n                  */
static
int kf_factor(int n,int16 * facbuf)
{
    int p=4;
    int i;
    int stages=0;
    int nbak = n;

    /*factor out powers of 4, powers of 2, then any remaining primes */
    do {
        while (n % p) {
            switch (p) {
                case 4: p = 2; break;
                case 2: p = 3; break;
                default: p += 2; break;
            }
            if (p>32000 || (int32)p*(int32)p > n)
                p = n;          /* no more factors, skip to end */
        }
        n /= p;
#ifdef RADIX_TWO_ONLY
        if (p!=2 && p != 4)
#else
        if (p>5)
#endif
        {
            return 0;
        }
        facbuf[2*stages] = p;
        if (p==2 && stages > 1)
        {
            facbuf[2*stages] = 4;
            facbuf[2] = 2;
        }
        stages++;
    } while (n > 1);
    n = nbak;
    /* Reverse the order to get the radix 4 at the end, so we can use the
       fast degenerate case. It turns out that reversing the order also
       improves the noise behaviour. */
    for (i=0;i<stages/2;i++)
    {
        int tmp;
        tmp = facbuf[2*i];
        facbuf[2*i] = facbuf[2*(stages-i-1)];
        facbuf[2*(stages-i-1)] = tmp;
    }
    for (i=0;i<stages;i++)
    {
        n /= facbuf[2*i];
        facbuf[2*i+1] = n;
    }
    return 1;
}

static void compute_twiddles(complex *twiddles, int nfft)
{
    int i;
#ifdef FIXED_POINT
    for (i=0;i<nfft;++i) {
      float phase = -i;
      kf_cexp2(twiddles+i, DIV(phase,nfft));
   }
#else
    for (i=0;i<nfft;++i) {
        const double pi=3.14159265358979323846264338327;
        double phase = ( -2*pi /nfft ) * i;
        kf_cexp(twiddles+i, phase );
    }
#endif
}

int fft_alloc_arch(_fft_state *st) {
    (void)st;
    return 0;
}

/*
 *
 * Allocates all necessary storage space for the fft and ifft.
 * The return value is a contiguous block of memory.  As such,
 * It can be freed with free().
 * */
_fft_state *fft_alloc_twiddles(int nfft,void * mem,size_t * lenmem,
                                        const _fft_state *base, int arch)
{
    _fft_state *st=NULL;
    size_t memneeded = sizeof(struct _fft_state);

    if ( lenmem==NULL ) {
        st = (_fft_state*)malloc(memneeded);
    }else{
        if (mem != NULL && *lenmem >= memneeded)
            st = (_fft_state*)mem;
        *lenmem = memneeded;
    }
    if (st) {
        int16 *bitrev;
        complex *twiddles;

        st->nfft=nfft;
//#ifdef FIXED_POINT
//        st->scale_shift = celt_ilog2(st->nfft);
//        if (st->nfft == 1<<st->scale_shift)
//           st->scale = 32767;
//        else
//           st->scale = (1073741824+st->nfft/2)/st->nfft>>(15-st->scale_shift);
//#else
        st->scale = 1.f/nfft;
//#endif
        if (base != NULL)
        {
            st->twiddles = base->twiddles;
            st->shift = 0;
            while (st->shift < 32 && nfft<<st->shift != base->nfft)
                st->shift++;
            if (st->shift>=32)
                goto fail;
        } else {
            st->twiddles = twiddles = (complex*)malloc(sizeof(complex)*nfft);
            compute_twiddles(twiddles, nfft);
            st->shift = -1;
        }
        if (!kf_factor(nfft,st->factors))
        {
            goto fail;
        }

        /* bitrev */
        st->bitrev = bitrev = (int16*)malloc(sizeof(int16)*nfft);
        if (st->bitrev==NULL)
            goto fail;
        compute_bitrev_table(0, bitrev, 1,1, st->factors,st);

        /* Initialize architecture specific fft parameters */
        if (fft_alloc_arch(st))
            goto fail;
    }
    return st;
    fail:
    fft_free(st, arch);
    return NULL;
}

_fft_state *fft_alloc(int nfft,void * mem,size_t * lenmem, int arch)
{
    return fft_alloc_twiddles(nfft, mem, lenmem, NULL, arch);
}

void fft_free_arch(_fft_state *st) {
    (void)st;
}

void fft_free(const _fft_state *cfg, int arch)
{
    if (cfg)
    {
        fft_free_arch((_fft_state *)cfg);
        free((int16*)cfg->bitrev);
        if (cfg->shift < 0)
            free((complex*)cfg->twiddles);
        free((_fft_state*)cfg);
    }
}

#endif /* CUSTOM_MODES */

void fft_main(const _fft_state *st,complex *fout)
{
    int m2, m;
    int p;
    int L;
    int fstride[MAXFACTORS];
    int i;
    int shift;

    /* st->shift can be -1 */
    shift = st->shift>0 ? st->shift : 0;

    fstride[0] = 1;
    L=0;
    do {
        p = st->factors[2*L];
        m = st->factors[2*L+1];
        fstride[L+1] = fstride[L]*p;
        L++;
    } while(m!=1);
    m = st->factors[2*L-1];
    for (i=L-1;i>=0;i--)
    {
        if (i!=0)
            m2 = st->factors[2*i-1];
        else
            m2 = 1;
        switch (st->factors[2*i])
        {
            case 2:
                kf_bfly2(fout, m, fstride[i]);
                break;
            case 4:
                kf_bfly4(fout,fstride[i]<<shift,st,m, fstride[i], m2);
                break;
#ifndef RADIX_TWO_ONLY
            case 3:
                kf_bfly3(fout,fstride[i]<<shift,st,m, fstride[i], m2);
                break;
            case 5:
                kf_bfly5(fout,fstride[i]<<shift,st,m, fstride[i], m2);
                break;
#endif
        }
        m = m2;
    }
}

void fft_c(const _fft_state *st,const complex *fin,complex *fout)
{
    int i;
    float scale;
#ifdef FIXED_POINT
    /* Allows us to scale with MULT16_32_Q16(), which is faster than
      MULT16_32_Q15() on ARM. */
   int scale_shift = st->scale_shift-1;
#endif
    scale = st->scale;

//    celt_assert2 (fin != fout, "In-place FFT not supported");
    /* Bit-reverse the input */
    for (i=0;i<st->nfft;i++)
    {
        complex x = fin[i];
        fout[st->bitrev[i]].r = MUL(scale, x.r);
        fout[st->bitrev[i]].i = MUL(scale, x.i);
    }
    fft_main(st, fout);
}


