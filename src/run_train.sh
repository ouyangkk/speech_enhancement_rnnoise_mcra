
#! /bin/bash
#! /bin/bash


stage=0
if [ $stage -le 0 ]; then
gcc -DTRAINING=1 -Wall -W -g -I../include train_main.c denoise_function.c fft.c lpc_pitch.c rnn_weight.c rnn.c -o feature_extract -lm
fi

# Your folder should first contain the collected speech and noise signals, which may be pre-processed before, such as sample rate conversion and speech splicing operations
#if [ $stage -le 1 ]; then
./wav2pcm clean.wav clean.pcm
./wav2pcm noise.wav noise.pcm
fi

#This step is feature extraction
if [ $stage -le 2 ]; then
./feature_extract clean.pcm noise.pcm 12000 >feature.f32
fi

#Make sure your computer includes the python3 environment,This is the preprocessing step for training
if [ $stage -le 3 ]; then
python3 bin2hdf5.py feature.f32 12000 75 feature.h5
fi
#Make sure your computer contains tensorflow and keras environment,This step requires waiting. The waiting time is related to the amount of data and model design.
if [ $stage -le 4 ]; then
python3 train.py
fi

if [ $stage -le 5 ]; then
python3 python2c.py weights.hdf5 rnn_weight.c
fi


