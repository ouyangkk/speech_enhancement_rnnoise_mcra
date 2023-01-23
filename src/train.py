#!/usr/bin/python

from __future__ import print_function

import keras
from keras.models import Sequential
import tensorflow as tf
from keras.models import Model
from keras.layers import Input
from keras.layers import Dense
from keras.layers import LSTM
from keras.layers import GRU
from keras.layers import SimpleRNN
from keras.layers import Dropout
from keras.layers import concatenate
from keras import losses
from keras import regularizers
from keras.constraints import min_max_norm
import h5py
import keras.backend.tensorflow_backend as KTF

from keras.constraints import Constraint
from keras import backend as K
import numpy as np

KTF.set_session(tf.Session(config = tf.ConfigProto(device_count={'gpu':0})))
#import tensorflow as tf
#from keras.backend.tensorflow_backend import set_session
#config = tf.ConfigProto()
#config.gpu_options.per_process_gpu_memory_fraction = 0.42
#set_session(tf.Session(config=config))


def my_crossentropy(y_true, y_pred):
    return K.mean(2*K.abs(y_true-0.5) * K.binary_crossentropy(y_pred+1e-8, y_true), axis=-1)

def mymask(y_true):
    return K.minimum(y_true+1., 1.)

def msse(y_true, y_pred):
    return K.mean(mymask(y_true) * K.square(K.sqrt(y_pred+1e-8) - K.sqrt(y_true)), axis=-1)

def mycost(y_true, y_pred):
    return K.mean(mymask(y_true+1e-8) * (10*K.square(K.square((K.sqrt(y_pred+1e-5) - K.sqrt(y_true+1e-8))+(K.sqrt((
                                                                                                                     y_pred+1e-5)*(
            1-y_true+1e-8)))) + K.square(K.sqrt(y_pred+1e-8) - K.sqrt(y_true+1e-8)) + 0.01*K.binary_crossentropy((
            y_pred+1e-8),
                                                                                                  y_true+1e-8))), axis=-1)

def my_accuracy(y_true, y_pred):
    return K.mean(2*K.abs(y_true-0.5) * K.equal(y_true, K.round(y_pred+1e-5)), axis=-1)

class WeightClip(Constraint):
    '''Clips the weights incident to each hidden unit to be inside a range
    '''
    def __init__(self, c=2):
        self.c = c

    def __call__(self, p):
        return K.clip(p, -self.c, self.c)

    def get_config(self):
        return {'name': self.__class__.__name__,
            'c': self.c}

reg = 0.000001
constraint = WeightClip(0.499)

print('Build model...')
main_input = Input(shape=(None, 38), name='main_input')
tmp = Dense(24, activation='tanh', name='input_dense', kernel_constraint=constraint, bias_constraint=constraint)(main_input)
vad_gru = GRU(24, activation='tanh', recurrent_activation='sigmoid', return_sequences=True, name='vad_gru', kernel_regularizer=regularizers.l2(reg), recurrent_regularizer=regularizers.l2(reg), kernel_constraint=constraint, recurrent_constraint=constraint, bias_constraint=constraint)(tmp)
vad_output = Dense(1, activation='sigmoid', name='vad_output', kernel_constraint=None, bias_constraint=None)(vad_gru)
noise_input = keras.layers.concatenate([tmp, vad_gru, main_input])
noise_gru = GRU(48, activation='relu', recurrent_activation='sigmoid', return_sequences=True, name='noise_gru', kernel_regularizer=regularizers.l2(reg), recurrent_regularizer=regularizers.l2(reg), kernel_constraint=constraint, recurrent_constraint=constraint, bias_constraint=constraint)(noise_input)
denoise_input = keras.layers.concatenate([vad_gru, noise_gru, main_input])

denoise_gru = GRU(96, activation='tanh', recurrent_activation='sigmoid', return_sequences=True, name='denoise_gru', kernel_regularizer=regularizers.l2(reg), recurrent_regularizer=regularizers.l2(reg), kernel_constraint=constraint, recurrent_constraint=constraint, bias_constraint=constraint)(denoise_input)

denoise_output = Dense(18, activation='sigmoid', name='denoise_output', kernel_constraint=constraint, bias_constraint=constraint)(denoise_gru)

model = Model(inputs=main_input, outputs=[denoise_output, vad_output])

model.compile(loss=[mycost, my_crossentropy],
              metrics=[msse],
              optimizer='adam', loss_weights=[10, 0.5])


batch_size = 128

print('Loading data...')
with h5py.File('feature.h5', 'r') as hf:
    all_data = hf['data'][:]
print('done.')

window_size = 2000

nb_sequences = len(all_data)//window_size
print(nb_sequences, ' sequences')
x_train = all_data[:nb_sequences*window_size, :38]
print(x_train)
x_train = np.reshape(x_train, (nb_sequences, window_size, 38))

y_train = np.copy(all_data[:nb_sequences*window_size, 38:56])
print(y_train)
y_train = np.reshape(y_train, (nb_sequences, window_size, 18))

noise_train = np.copy(all_data[:nb_sequences*window_size, 56:74])
noise_train = np.reshape(noise_train, (nb_sequences, window_size, 18))

vad_train = np.copy(all_data[:nb_sequences*window_size, 74:75])
vad_train = np.reshape(vad_train, (nb_sequences, window_size, 1))

all_data = 0;
#x_train = x_train.astype('float32')
#y_train = y_train.astype('float32')

print(len(x_train), 'train sequences. x shape =', x_train.shape, 'y shape = ', y_train.shape)
# physical_devices = tf.config.experimental.list_physical_devices("GPU")
# if (physical_devices is not None):
#     tf.config.experimental.set_memory_growth(physical_devices[0], True)

print('Train...')
model.fit(x_train, [y_train, vad_train],
          batch_size=batch_size,
          epochs=60,
          validation_split=0.1)
model.save("weights.hdf5")
