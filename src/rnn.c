//
// Created by 欧阳康的macbook on 2022/2/11.
//
#include <math.h>
#include "Platform.h"
#include "table.h"
#include "rnn.h"
#include "rnn_weight.h"
#include <stdio.h>

static inline float tansig_approx(float x)
{
    int i;
    float y, dy;
    float sign=1;
    /* Tests are reversed to catch NaNs */
    if (!(x<8))
        return 1;
    if (!(x>-8))
        return -1;
//#ifndef FIXED_POINT
//    /* Another check in case of -ffast-math */
//    if (celt_isnan(x))
//        return 0;
//#endif
    if (x<0)
    {
        x=-x;
        sign=-1;
    }
    i = (int)floor(.5f+25*x);
    x -= .04f*i;
    y = tansig_table[i];
    dy = 1-y*y;
    y = y + x*dy*(1 - y*x);
    return sign*y;
}

static inline float sigmoid_approx(float x)
{
    return .5 + .5*tansig_approx(.5*x);
}

static inline float relu(float x)
{
    return x < 0 ? 0 : x;
}

void compute_dense(const DenseLayer *layer, float *output, const float *input)
{
    int i, j;
    int N, M;
    int stride;
    M = layer->nb_inputs;
    N = layer->nb_neurons;
    stride = N;
    for (i=0;i<N;i++)
    {
        /* Compute update gate. */
        float sum = layer->bias[i];
        for (j=0;j<M;j++)
            sum += layer->input_weights[j*stride + i]*input[j];
        output[i] = WEIGHTS_SCALE*sum;
    }
    if (layer->activation == ACTIVATION_SIGMOID) {
        for (i=0;i<N;i++)
            output[i] = sigmoid_approx(output[i]);
    } else if (layer->activation == ACTIVATION_TANH) {
        for (i=0;i<N;i++)
            output[i] = tansig_approx(output[i]);
    } else if (layer->activation == ACTIVATION_RELU) {
        for (i=0;i<N;i++)
            output[i] = relu(output[i]);
    } else {
        *(int*)0=0;
    }
}

void compute_gru(const GRULayer *gru, float *state, const float *input)
{
    int i, j;
    int N, M;
    int stride;
    float z[MAX_NEURONS];
    float r[MAX_NEURONS];
    float h[MAX_NEURONS];
    M = gru->nb_inputs;
    N = gru->nb_neurons;
    stride = 3*N;
    for (i=0;i<N;i++)
    {
        /* Compute update gate. */
        float sum = gru->bias[i];
        for (j=0;j<M;j++)
            sum += gru->input_weights[j*stride + i]*input[j];
        for (j=0;j<N;j++)
            sum += gru->recurrent_weights[j*stride + i]*state[j];
        z[i] = sigmoid_approx(WEIGHTS_SCALE*sum);
    }
    for (i=0;i<N;i++)
    {
        /* Compute reset gate. */
        float sum = gru->bias[N + i];
        for (j=0;j<M;j++)
            sum += gru->input_weights[N + j*stride + i]*input[j];
        for (j=0;j<N;j++)
            sum += gru->recurrent_weights[N + j*stride + i]*state[j];
        r[i] = sigmoid_approx(WEIGHTS_SCALE*sum);
    }
    for (i=0;i<N;i++)
    {
        /* Compute output. */
        float sum = gru->bias[2*N + i];
        for (j=0;j<M;j++)
            sum += gru->input_weights[2*N + j*stride + i]*input[j];
        for (j=0;j<N;j++)
            sum += gru->recurrent_weights[2*N + j*stride + i]*state[j]*r[j];
        if (gru->activation == ACTIVATION_SIGMOID) sum = sigmoid_approx(WEIGHTS_SCALE*sum);
        else if (gru->activation == ACTIVATION_TANH) sum = tansig_approx(WEIGHTS_SCALE*sum);
        else if (gru->activation == ACTIVATION_RELU) sum = relu(WEIGHTS_SCALE*sum);
        else *(int*)0=0;
        h[i] = z[i]*state[i] + (1-z[i])*sum;
    }
    for (i=0;i<N;i++)
        state[i] = h[i];
}

//这个是减少模型参数的关键因子
#define  INPUT_SIZE 38

void compute_rnn(RNNState *rnn, float *gains, float *vad, const float *input) {
    int i;
    float dense_out[MAX_NEURONS];
    float noise_input[MAX_NEURONS*3];
    float denoise_input[MAX_NEURONS*3];
    compute_dense(&input_dense, dense_out, input);
    compute_gru(&vad_gru, rnn->vad_gru_state, dense_out);
    compute_dense(&vad_output, vad, rnn->vad_gru_state);
    for (i=0;i<INPUT_DENSE_SIZE;i++) noise_input[i] = dense_out[i];
    for (i=0;i<VAD_GRU_SIZE;i++) noise_input[i+INPUT_DENSE_SIZE] = rnn->vad_gru_state[i];
    for (i=0;i<INPUT_SIZE;i++) noise_input[i+INPUT_DENSE_SIZE+VAD_GRU_SIZE] = input[i];
    compute_gru(&noise_gru, rnn->noise_gru_state, noise_input);

    for (i=0;i<VAD_GRU_SIZE;i++) denoise_input[i] = rnn->vad_gru_state[i];
    for (i=0;i<NOISE_GRU_SIZE;i++) denoise_input[i+VAD_GRU_SIZE] = rnn->noise_gru_state[i];
    for (i=0;i<INPUT_SIZE;i++) denoise_input[i+VAD_GRU_SIZE+NOISE_GRU_SIZE] = input[i];
    compute_gru(&denoise_gru, rnn->denoise_gru_state, denoise_input);
    compute_dense(&denoise_output, gains, rnn->denoise_gru_state);
}

