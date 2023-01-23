//
// Created by 欧阳康的macbook on 2022/2/11.
//
#ifndef PLATFORM_H
#define PLATFORM_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>




static inline void *rnnoise_alloc (size_t size)
{
    return malloc(size);
}

static inline void rnnoise_free (void *ptr)
{
    free(ptr);
}

#define RNN_COPY(dst, src, n) (memcpy((dst), (src), (n)*sizeof(*(dst)) + 0*((dst)-(src)) ))
#define RNN_MOVE(dst, src, n) (memmove((dst), (src), (n)*sizeof(*(dst)) + 0*((dst)-(src)) ))
#define RNN_CLEAR(dst, n) (memset((dst), 0, (n)*sizeof(*(dst))))

//#define __GNUC_PREREQ(_maj,_min) 0
#define CHECK_OVERFLOW_OP(a,op,b)
//#define assert(cond) {if (!(cond)) {celt_fatal("assertion failed: " #cond);}}
//C99 standard C
#if (defined(__STDC__) && __STDC__ && defined(__STDC_VERSION__) && __STDC_VERSION__ >= 199901L) || (defined(__GNUC__) && (defined(_STDINT_H) || defined(_STDINT_H_)) || defined (HAVE_STDINT_H))
#include <stdint.h>

   typedef int16_t int16;
   typedef uint16_t uint16;
   typedef int32_t int32;
   typedef uint32_t uint32;

// windows x86 simulation
// unix on windows
#  elif defined(__CYGWIN__)
#    include <_G_config.h>
     typedef _G_int32_t tws_int32;
     typedef _G_uint32_t tws_uint32;
     typedef _G_int16 tws_int16;
     typedef _G_uint16 tws_uint16;
// mingw such as clion simulation on windows
#  elif defined(__MINGW32__)
     typedef short tws_int16;
     typedef unsigned short tws_uint16;
     typedef int tws_int32;
     typedef unsigned int tws_uint32;

// MACOS system simulation
#elif defined(__MACOS__)
#  include <sys/types.h>
   typedef SInt16 int16;
   typedef UInt16 uint16;
   typedef SInt32 int32;
   typedef UInt32 uint32;

// TI  series
#elif defined(TI_C54X) || defined (TI_C55X)

   typedef short int16;
   typedef unsigned short uint16;
   typedef long int32;
   typedef unsigned long uint32;

#elif defined(TI_C6X)

   typedef short int16;
   typedef unsigned short uint16;
   typedef int int32;
   typedef unsigned int uint32;

#else
   typedef short int16;
   typedef unsigned short uint16;
   typedef int int32;
   typedef unsigned int uint32;

#endif

#define DIV(a,b)     (((float)(a))/(float)(b))
#define ADD(a,b) ((a)+(b))
#define SUB(a,b) ((a)-(b))
//#define ADD32(a,b) ((a)+(b))
//#define SUB32(a,b) ((a)-(b))

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

#define ABS(x) ((float)fabs(x))
//#define ABS32(x) ((float)fabs(x))

#define NEG(x) (-(x))
#define HALF(x)   (.5f*(x))

#define MUL(a,b) ((a)*(b))
//#define MULT16_16(a,b)  ((tws_val32)(a)*(tws_val32)(b))
#define MAC(c,a,b)     ((c)+(a)*(b))


#define SCALE 32768.f
#define SCALEIN(x)      ((x)*SCALE16)
#define SCALEOUT(x)     ((x)*(1/SCALE16))
#endif


