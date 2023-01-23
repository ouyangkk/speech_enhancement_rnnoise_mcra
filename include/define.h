//
// Created by 欧阳康的macbook on 2022/2/11.
//

#ifndef OM_DEBUG_DEFINE_H
#define OM_DEBUG_DEFINE_H
#define PI 3.14159265358979323846264338327950288
#define DATA_RANDOM 1
#define DATA_NOTRANDOM 0

#ifndef TRAINING
#define TRAINING 0
#endif

#define DSP_ONLY 0
#define RNN_ONLY 0
#define DSP_RNN 1
#define WEIGHT_FACTOR 0.5
#define PITCH_F 1
#define SHIFT 2

/*This will affect the number of fft points.
 If there is a delay, you can increase this value to 256, but you need to retrain the model */
#define FRAME_SIZE 160

#define NOISE_FACTOR 0.0001
#define  P_MIN 0.0001  //POW
#define  SNR_MAX 0.7//POW
#define  SNR_MIN 0.45 //pow
#define  SNR_P_MIN 0.0001
#define  SNR_P_MAX_MIN (log(SNR_P_MAX)-log(SNR_P_MIN))
#define  SNR_P_MAX 10
#define  WIN_HALF 4
#define  WIN_ALL 9
#define  PR_SNR 0.7
#define  C_F (FRAME_SIZE>>1)
#define WINDOW_SIZE (2*FRAME_SIZE)
#define FREQ_SIZE (FRAME_SIZE + 1)
#define SQUARE(x) ((x)*(x))


//When the audio sample rate is 16k, this value is fixed
#define PITCH_MIN_PERIOD 20
#define PITCH_MAX_PERIOD 256
#define PITCH_FRAME_SIZE 320
#define PITCH_BUF_SIZE (PITCH_MAX_PERIOD+PITCH_FRAME_SIZE)

//This will affect the size of the model
#define NB_BANDS 18

#define CEPS_MEM 8
#define NB_DELTA_CEPS 6

//number of features
#define NB_FEATURES (NB_BANDS+3*NB_DELTA_CEPS+2)

#endif //OM_DEBUG_DEFINE_H
