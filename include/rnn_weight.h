//
// Created by 欧阳康的macbook on 2022/2/11.
//

#ifndef OM_DEBUG_RNN_WEIGHT_H
#define OM_DEBUG_RNN_WEIGHT_H
#include "rnn.h"

#define INPUT_DENSE_SIZE 24
extern const DenseLayer input_dense;

#define VAD_GRU_SIZE 24
extern const GRULayer vad_gru;

#define NOISE_GRU_SIZE 48
extern const GRULayer noise_gru;

#define DENOISE_GRU_SIZE 96
extern const GRULayer denoise_gru;

#define DENOISE_OUTPUT_SIZE 18
extern const DenseLayer denoise_output;

#define VAD_OUTPUT_SIZE 1
extern const DenseLayer vad_output;

struct RNNState {
    float vad_gru_state[VAD_GRU_SIZE];
    float noise_gru_state[NOISE_GRU_SIZE];
    float denoise_gru_state[DENOISE_GRU_SIZE];
};


#endif //OM_DEBUG_RNN_WEIGHT_H
