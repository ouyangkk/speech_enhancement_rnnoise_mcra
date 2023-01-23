#! /bin/bash
#! /bin/bash


stage=0

if [ $stage -le 0 ]; then
./wav2pcm noisy.wav noisy.pcm
fi

if [ $stage -le 1 ]; then
./DENOISE noisy.pcm enhance.pcm
fi

if [ $stage -le 2 ]; then
./pcm2wav 1 16000 16 enhance.pcm enhance.wav
fi