import keras
from keras import layers
from keras import regularizers
import pandas as pd
import numpy as np
from scipy import stats
import sys,math

def model_training_parameter(data,encoding_dim,layer_dim,max_epoch,log_file):
  Standard_data=stats.zscore(data,axis=0)
  Standard_data[np.where(np.isnan(Standard_data))]=0
  
  x_train=Standard_data
  x_test=Standard_data
  
  input_img = keras.Input(shape=(x_train.shape[1],))
  encoded=input_img
  for i in range(len(layer_dim)):
    encoded = layers.Dense(int(layer_dim[i]), activation='relu',
                         activity_regularizer=regularizers.l1(10e-5))(encoded)
  encoded = layers.Dense(int(encoding_dim), activation='relu')(encoded)
  decoded=encoded
  for i in range(len(layer_dim)):
    decoded = layers.Dense(int(layer_dim[-1*i]), activation='relu',
                         activity_regularizer=regularizers.l1(10e-5))(decoded)
                         
  decoded = layers.Dense(x_train.shape[1],activation=None)(decoded)
  
  autoencoder =keras.Model(input_img, decoded)
  encoder = keras.Model(input_img,encoded)
  
  autoencoder.compile(optimizer='adam', loss='MSE')
  
  # creat CSVLogger callbacks
  csv_logger = CSVLogger(log_file, append=True)

  autoencoder.fit(x_train, x_train,
                  epochs=int(max_epoch),
                  batch_size=min(10000,x_train.shape[0]),
                  shuffle=True,
                  validation_data=(x_test, x_test),
                  callbacks=[csv_logger])

  encoded_imgs = encoder.predict(Standard_data)
  encoded_imgs=pd.DataFrame(encoded_imgs)

  return encoded_imgs
