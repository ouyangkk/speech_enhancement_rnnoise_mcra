#include "lpc_pitch.h"
#include "Platform.h"

static inline void xcorr_kernel(const float * x, const float * y, float sum[4], int len)
{
    int j;
    float y_0, y_1, y_2, y_3;
    //celt_assert(len>=3);
    y_3=0; /* gcc doesn't realize that y_3 can't be used uninitialized */
    y_0=*y++;
    y_1=*y++;
    y_2=*y++;
    for (j=0;j<len-3;j+=4)
    {
        float tmp;
        tmp = *x++;
        y_3=*y++;
        sum[0] = MAC(sum[0],tmp,y_0);
        sum[1] = MAC(sum[1],tmp,y_1);
        sum[2] = MAC(sum[2],tmp,y_2);
        sum[3] = MAC(sum[3],tmp,y_3);
        tmp=*x++;
        y_0=*y++;
        sum[0] = MAC(sum[0],tmp,y_1);
        sum[1] = MAC(sum[1],tmp,y_2);
        sum[2] = MAC(sum[2],tmp,y_3);
        sum[3] = MAC(sum[3],tmp,y_0);
        tmp=*x++;
        y_1=*y++;
        sum[0] = MAC(sum[0],tmp,y_2);
        sum[1] = MAC(sum[1],tmp,y_3);
        sum[2] = MAC(sum[2],tmp,y_0);
        sum[3] = MAC(sum[3],tmp,y_1);
        tmp=*x++;
        y_2=*y++;
        sum[0] = MAC(sum[0],tmp,y_3);
        sum[1] = MAC(sum[1],tmp,y_0);
        sum[2] = MAC(sum[2],tmp,y_1);
        sum[3] = MAC(sum[3],tmp,y_2);
    }
    if (j++<len)
    {
        float tmp = *x++;
        y_3=*y++;
        sum[0] = MAC(sum[0],tmp,y_0);
        sum[1] = MAC(sum[1],tmp,y_1);
        sum[2] = MAC(sum[2],tmp,y_2);
        sum[3] = MAC(sum[3],tmp,y_3);
    }
    if (j++<len)
    {
        float tmp=*x++;
        y_0=*y++;
        sum[0] = MAC(sum[0],tmp,y_1);
        sum[1] = MAC(sum[1],tmp,y_2);
        sum[2] = MAC(sum[2],tmp,y_3);
        sum[3] = MAC(sum[3],tmp,y_0);
    }
    if (j<len)
    {
        float tmp=*x++;
        y_1=*y++;
        sum[0] = MAC(sum[0],tmp,y_2);
        sum[1] = MAC(sum[1],tmp,y_3);
        sum[2] = MAC(sum[2],tmp,y_0);
        sum[3] = MAC(sum[3],tmp,y_1);
    }
}

static inline void dual_inner_prod(const float *x, const float *y01, const float *y02,
                                        int N, float *xy1, float *xy2)
{
    int i;
    float xy01=0;
    float xy02=0;
    for (i=0;i<N;i++)
    {
        xy01 = MAC(xy01, x[i], y01[i]);
        xy02 = MAC(xy02, x[i], y02[i]);
    }
    *xy1 = xy01;
    *xy2 = xy02;
}

static inline float single_inner_prod(const float *x,
                                              const float *y, int N)
{
    int i;
    float xy=0;
    for (i=0;i<N;i++)
        xy = MAC(xy, x[i], y[i]);
    return xy;
}
//Filters associated with autocorrelation
static void fir_filter(const float *x,
                      const float *num,
                      float *y,
                      int N,
                      float *mem)
{
    int i;
    float num0, num1, num2, num3, num4;
    float mem0, mem1, mem2, mem3, mem4;
    num0=num[0];
    num1=num[1];
    num2=num[2];
    num3=num[3];
    num4=num[4];
    mem0=mem[0];
    mem1=mem[1];
    mem2=mem[2];
    mem3=mem[3];
    mem4=mem[4];
    for (i=0;i<N;i++)
    {
        float sum = x[i];
        sum = MAC(sum,num0,mem0);
        sum = MAC(sum,num1,mem1);
        sum = MAC(sum,num2,mem2);
        sum = MAC(sum,num3,mem3);
        sum = MAC(sum,num4,mem4);
        mem4 = mem3;
        mem3 = mem2;
        mem2 = mem1;
        mem1 = mem0;
        mem0 = x[i];
        y[i] = sum;
    }
    mem[0]=mem0;
    mem[1]=mem1;
    mem[2]=mem2;
    mem[3]=mem3;
    mem[4]=mem4;
}
//
void pitch_xcorr(const float *_x, const float *_y,
                      float *xcorr, int len, int max_pitch)
{

#if 0 /* This is a simple version of the pitch correlation that should work
         well on DSPs like Blackfin and TI C5x/C6x */
   int i, j;
   for (i=0;i<max_pitch;i++)
   {
      float sum = 0;
      for (j=0;j<len;j++)
         sum = MAC(sum, _x[j], _y[i+j]);
      xcorr[i] = sum;
   }

#else /* Unrolled version of the pitch correlation -- runs faster on x86 and ARM */
    int i;
    /*The EDSP version requires that max_pitch is at least 1, and that _x is
       32-bit aligned.
      Since it's hard to put asserts in assembly, put them here.*/
    //celt_assert(max_pitch>0);
    //celt_assert((((unsigned char *)_x-(unsigned char *)NULL)&3)==0);
    for (i=0;i<max_pitch-3;i+=4)
    {
        float sum[4]={0,0,0,0};
        xcorr_kernel(_x, _y+i, sum, len);
        xcorr[i]=sum[0];
        xcorr[i+1]=sum[1];
        xcorr[i+2]=sum[2];
        xcorr[i+3]=sum[3];
    }
    /* In case max_pitch isn't a multiple of 4, do non-unrolled version. */
    for (;i<max_pitch;i++)
    {
        float sum;
        sum = single_inner_prod(_x, _y+i, len);
        xcorr[i] = sum;
    }
#endif
}

//alculate LPC coefficients
void lpc_analyze(
        float       *_lpc, /* out: [0...p-1] LPC coefficients      */
        const float *ac,  /* in:  [0...p] autocorrelation values  */
        int          p
)
{
    int i, j;
    float r;
    float error = ac[0];
    float *lpc = _lpc;

    RNN_CLEAR(lpc, p);
    if (ac[0] != 0)
    {
        for (i = 0; i < p; i++) {
            /* Sum up this iteration's reflection coefficient */
            float rr = 0;
            for (j = 0; j < i; j++)
                rr += MUL(lpc[j],ac[i - j]);
            rr += ac[i + 1];
            r = -rr/error;
            /*  Update LPC coefficients and total error */
            lpc[i] = r;
            for (j = 0; j < (i+1)>>1; j++)
            {
                float tmp1, tmp2;
                tmp1 = lpc[j];
                tmp2 = lpc[i-1-j];
                lpc[j]     = tmp1 + MUL(r,tmp2);
                lpc[i-1-j] = tmp2 + MUL(r,tmp1);
            }

            error = error - MUL(MUL(r,r),error);
            /* Bail out once we get 30 dB gain */
            if (error<.001f*ac[0])
                break;
        }
    }
}
//alculate LPC coefficients Autocorrelation function
int lpc_autocorr(
        const float *x,   /*  in: [0...n-1] samples x   */
        float       *ac,  /* out: [0...lag-1] ac values */
        const float      *window,
        int          overlap,
        int          lag,
        int          n)
{
    float d;
    int i, k;
    int fastN=n-lag;
    int shift;
    const float *xptr;
    float xx[n];
    //celt_assert(n>0);
    //celt_assert(overlap>=0);
    if (overlap == 0)
    {
        xptr = x;
    } else {
        for (i=0;i<n;i++)
            xx[i] = x[i];
        for (i=0;i<overlap;i++)
        {
            xx[i] = MUL(x[i],window[i]);
            xx[n-i-1] = MUL(x[n-i-1],window[i]);
        }
        xptr = xx;
    }
    shift=0;
    //pitch autocorrelation
    pitch_xcorr(xptr, xptr, ac, fastN, lag+1);
    for (k=0;k<=lag;k++)
    {
        for (i = k+fastN, d = 0; i < n; i++)
            d = MAC(d, xptr[i], xptr[i-k]);
        ac[k] += d;
    }

    return shift;
}

void pitch_residual(float *x[], float *x_lp,\
                      int len, int C)
{
    int i;
    float ac[5];
    float tmp=32767.f;
    float lpc[4], mem[5]={0,0,0,0,0};
    float lpc2[5];
    float c1 = 0.8f;
    for (i=1;i<len>>1;i++)
        x_lp[i] = HALF(HALF(x[0][(2*i-1)]+x[0][(2*i+1)])+x[0][2*i]);
    x_lp[0] = HALF(HALF(x[0][1])+x[0][0]);
    if (C==2)
    {
        for (i=1;i<len>>1;i++)
            x_lp[i] += HALF(HALF(x[1][(2*i-1)]+x[1][(2*i+1)])+x[1][2*i]);
        x_lp[0] += HALF(HALF(x[1][1])+x[1][0]);
    }

    lpc_autocorr(x_lp, ac, NULL, 0,
                   4, len>>1);

    /* Noise floor -40 dB */
    ac[0] *= 1.0001f;
    /* Lag windowing */
    for (i=1;i<=4;i++)
    {

        //According to the author's explanation, this is experience
        ac[i] -= ac[i]*(.008f*i)*(.008f*i);
    }

    lpc_analyze(lpc, ac, 4);
    for (i=0;i<4;i++)
    {
        tmp = MUL(0.9f,tmp);
        lpc[i] = MUL(lpc[i], tmp);
    }
    /* Add a zero */
    lpc2[0] = lpc[0] + c1;
    lpc2[1] = lpc[1] + MUL(c1,lpc[0]);
    lpc2[2] = lpc[2] + MUL(c1,lpc[1]);
    lpc2[3] = lpc[3] + MUL(c1,lpc[2]);
    lpc2[4] = MUL(c1,lpc[3]);
    fir_filter(x_lp, lpc2, x_lp, len>>1, mem);
}

static void find_best_pitch(float *xcorr, float *y, int len,
                            int max_pitch, int *best_pitch)
{
    int i, j;
    float Syy=1;
    float best_num[2];
    float best_den[2];
    best_num[0] = -1;
    best_num[1] = -1;
    best_den[0] = 0;
    best_den[1] = 0;
    best_pitch[0] = 0;
    best_pitch[1] = 1;
    for (j=0;j<len;j++)
        Syy = ADD(Syy, MUL(y[j],y[j]));
    for (i=0;i<max_pitch;i++)
    {
        if (xcorr[i]>0)
        {
            float num;
            float xcorr16;
            xcorr16 = xcorr[i];
            /* Considering the range of xcorr16, this should avoid both underflows
               and overflows (inf) when squaring xcorr16 */
            xcorr16 *= 1e-12f;
            num = MUL(xcorr16,xcorr16);
            if (MUL(num,best_den[1]) > MUL(best_num[1],Syy))
            {
                if (MUL(num,best_den[0]) > MUL(best_num[0],Syy))
                {
                    best_num[1] = best_num[0];
                    best_den[1] = best_den[0];
                    best_pitch[1] = best_pitch[0];
                    best_num[0] = num;
                    best_den[0] = Syy;
                    best_pitch[0] = i;
                } else {
                    best_num[1] = num;
                    best_den[1] = Syy;
                    best_pitch[1] = i;
                }
            }
        }
        Syy += MUL(y[i+len],y[i+len]) - MUL(y[i],y[i]);
        Syy = MAX(1, Syy);
    }
}

void pitch_search(const float *x_lp, float *y,
                  int len, int max_pitch, int *pitch)
{
    int i, j;
    int lag;
    int best_pitch[2]={0,0};
    int offset;

    //celt_assert(len>0);
    //celt_assert(max_pitch>0);
    lag = len+max_pitch;

    float x_lp4[len>>2];
    float y_lp4[lag>>2];
    float xcorr[max_pitch>>1];

    /* Downsample by 2 again */
    for (j=0;j<len>>2;j++)
        x_lp4[j] = x_lp[2*j];
    for (j=0;j<lag>>2;j++)
        y_lp4[j] = y[2*j];


    /* Coarse search with 4x decimation */

    pitch_xcorr(x_lp4, y_lp4, xcorr, len>>2, max_pitch>>2);

    find_best_pitch(xcorr, y_lp4, len>>2, max_pitch>>2, best_pitch);

    for (i=0;i<max_pitch>>1;i++)
    {
        float sum;
        xcorr[i] = 0;
        if (abs(i-2*best_pitch[0])>2 && abs(i-2*best_pitch[1])>2)
            continue;
        sum = single_inner_prod(x_lp, y+i, len>>1);
        xcorr[i] = MAX(-1, sum);
    }
    find_best_pitch(xcorr, y, len>>1, max_pitch>>1, best_pitch);

    /* Refine by pseudo-interpolation */
    if (best_pitch[0]>0 && best_pitch[0]<(max_pitch>>1)-1)
    {
        float a, b, c;
        a = xcorr[best_pitch[0]-1];
        b = xcorr[best_pitch[0]];
        c = xcorr[best_pitch[0]+1];
        if ((c-a) > MUL(0.7f,b-a))
            offset = 1;
        else if ((a-c) > MUL(0.7f,b-c))
            offset = -1;
        else
            offset = 0;
    } else {
        offset = 0;
    }
    *pitch = 2*best_pitch[0]-offset;
}

static float compute_pitch_gain(float xy, float xx, float yy)
{
    return xy/sqrt(1+xx*yy);
}

static const int second_check[16] = {0, 0, 3, 2, 3, 2, 5, 2, 3, 2, 3, 2, 5, 2, 3, 2};

//eliminate Octave Errors
float remove_doubling(float *x, int maxperiod, int minperiod,
                           int N, int *T0_, int prev_period, float prev_gain)
{
    int k, i, T, T0;
    float g, g0;
    float pg;
    float xy,xx,yy,xy2;
    float xcorr[3];
    float best_xy, best_yy;
    int offset;
    int minperiod0;

    minperiod0 = minperiod;
    maxperiod /= 2;
    minperiod /= 2;
    *T0_ /= 2;
    prev_period /= 2;
    N /= 2;
    x += maxperiod;
    if (*T0_>=maxperiod)
        *T0_=maxperiod-1;

    T = T0 = *T0_;
    float yy_lookup[maxperiod+1];
    dual_inner_prod(x, x, x-T0, N, &xx, &xy);
    yy_lookup[0] = xx;
    yy=xx;
    for (i=1;i<=maxperiod;i++)
    {
        yy = yy+MUL(x[-i],x[-i])-MUL(x[N-i],x[N-i]);
        yy_lookup[i] = MAX(0, yy);
    }
    yy = yy_lookup[T0];
    best_xy = xy;
    best_yy = yy;
    g = g0 = compute_pitch_gain(xy, xx, yy);
    /* Look for any pitch at T/k */
    for (k=2;k<=15;k++)
    {
        int T1, T1b;
        float g1;
        float cont=0;
        float thresh;
        T1 = (2*T0+k)/(2*k);
        if (T1 < minperiod)
            break;
        /* Look for another strong correlation at T1b */
        if (k==2)
        {
            if (T1+T0>maxperiod)
                T1b = T0;
            else
                T1b = T0+T1;
        } else
        {
            T1b = (2*second_check[k]*T0+k)/(2*k);
        }
        dual_inner_prod(x, &x[-T1], &x[-T1b], N, &xy, &xy2);
        xy = HALF(xy + xy2);
        yy = HALF(yy_lookup[T1] + yy_lookup[T1b]);
        g1 = compute_pitch_gain(xy, xx, yy);
        if (abs(T1-prev_period)<=1)
            cont = prev_gain;
        else if (abs(T1-prev_period)<=2 && 5*k*k < T0)
            cont = HALF(prev_gain);
        else
            cont = 0;
        thresh = MAX(0.3f, MUL(0.7f,g0)-cont);
        /* Bias against very high pitch (very short period) to avoid false-positives
           due to short-term correlation */
        if (T1<3*minperiod)
            thresh = MAX(0.4f, MUL(0.85f,g0)-cont);
        else if (T1<2*minperiod)
            thresh = MAX(0.5f, MUL(0.9f,g0)-cont);
        if (g1 > thresh)
        {
            best_xy = xy;
            best_yy = yy;
            T = T1;
            g = g1;
        }
    }
    best_xy = MAX(0, best_xy);
    if (best_yy <= best_xy)
        pg = 32767;
    else
        pg = best_xy/(best_yy+1);

    for (k=0;k<3;k++)
        xcorr[k] = single_inner_prod(x, x-(T+k-1), N);
    if ((xcorr[2]-xcorr[0]) > MUL(0.7f,xcorr[1]-xcorr[0]))
        offset = 1;
    else if ((xcorr[0]-xcorr[2]) > MUL(0.7f,xcorr[1]-xcorr[2]))
        offset = -1;
    else
        offset = 0;
    if (pg > g)
        pg = g;
    *T0_ = 2*T+offset;

    if (*T0_<minperiod0)
        *T0_=minperiod0;
    return pg;
}
//
// Created by 欧阳康的macbook on 2022/2/10.
//

