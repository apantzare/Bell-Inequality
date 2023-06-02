# -*- coding: utf-8 -*-
"""
Created on Sat Apr 22 20:06:38 2023

@author: axel
"""

#import xgboost as xgb
import numpy as np
from sklearn.model_selection import train_test_split
import tensorflow as tf
from tensorflow import keras
#from keras.layers import Conv1D
from sklearn.decomposition import PCA
from tensorflow.keras.callbacks import ModelCheckpoint
import pandas as pd
#from sklearn import svm
#from sklearn import metrics
#from warnings import simplefilter
#from sklearn.exceptions import ConvergenceWarning
#import pickle

def custom_loss(y, y_pred):
    """
    Parameters
    # ----------
    y_pred   : predicted values for a batch if samples (should, are not be binary: 0 or 1)
    y   : correct values for the set of samples used (must be binary: 0 or 1)
    Returns
    -------
    out : the special loss
    """
    FP = tf.math.logical_and(tf.cast(y, tf.float32) == 0, tf.cast(y_pred, tf.float32) > 0.5)
    TP = tf.math.logical_and(tf.cast(y, tf.float32) == 1, tf.cast(y_pred, tf.float32) > 0.5)
    FPs = tf.math.reduce_sum(tf.cast(FP,tf.float32)*y_pred)
    #FN = tf.math.logical_and(tf.cast(y, tf.float32) == 1, tf.cast(y_pred, tf.float32) < 0.5)

    Precision = FPs /(tf.math.reduce_sum(tf.cast(TP,tf.float32))+1e-3)#tf.math.reduce_sum(tf.cast(FP,tf.float32))+
    return 0.05*Precision+tf.keras.losses.BinaryCrossentropy()(y, y_pred)#165


data=np.genfromtxt('data.dat')
target=np.genfromtxt('target.dat')

pca = PCA(n_components=40)
dataPCA = pca.fit_transform(data)

X_train, X_test, y_train, y_test = train_test_split(dataPCA, target, test_size=0.01, random_state=42)
train_dataset=tf.data.Dataset.from_tensor_slices((X_train, y_train))
test_dataset=tf.data.Dataset.from_tensor_slices((X_test, y_test))

train_dataset = train_dataset.batch(100)
test_dataset = test_dataset.batch(100)

Neuron=1000
drop=0.01
strategy=tf.distribute.MirroredStrategy()
with strategy.scope():
	model = tf.keras.models.Sequential([ 
    	tf.keras.layers.InputLayer(input_shape=(dataPCA.shape[1])), #lägga loop för att lägga till lager?
    	tf.keras.layers.BatchNormalization(),#bra för overfitting
    	tf.keras.layers.Dropout(drop),
		tf.keras.layers.Dense(Neuron,  activation='relu'),#1 hidden
		tf.keras.layers.Dropout(drop),
		tf.keras.layers.Dense(Neuron,  activation='relu'),#2 hidden
		tf.keras.layers.Dropout(drop),
		tf.keras.layers.Dense(Neuron,  activation='relu'),#3 hidden kernel_regularizer=tf.keras.regularizers.l2(0.0001),
		tf.keras.layers.Dropout(drop),
		tf.keras.layers.Dense(Neuron,  activation='relu'),#4 hidden
		tf.keras.layers.Dropout(drop),
		tf.keras.layers.Dense(Neuron,  activation='relu'),#5 hidden
		tf.keras.layers.Dropout(drop),
		tf.keras.layers.Dense(Neuron,  activation='relu'),#6 hidden
		tf.keras.layers.Dropout(drop),
		tf.keras.layers.Dense(Neuron,  activation='relu'),#7 hidden
		tf.keras.layers.Dropout(drop),
		tf.keras.layers.Dense(Neuron,  activation='relu'),#8 hidden    
		tf.keras.layers.Dropout(drop),
		tf.keras.layers.Dense(Neuron,  activation='relu'),#8 hidden    
		tf.keras.layers.Dropout(drop),
		tf.keras.layers.Dense(Neuron,  activation='relu'),#8 hidden    
		tf.keras.layers.Dropout(drop),
		tf.keras.layers.Dense(1, activation='hard_sigmoid')
		])

	model.compile(optimizer=keras.optimizers.experimental.Adam(learning_rate=0.0001),
		loss=tf.keras.losses.BinaryCrossentropy(), metrics=['BinaryAccuracy' ,'Precision'])


filepath = 'DNN_balanced_10layers_1000neurons_0001learning_0drop.hdf5'

checkpoint = ModelCheckpoint(filepath, monitor='val_loss', verbose=1, \
                              save_best_only=True, save_weights_only=False, \
                             mode='auto', save_frequency=1, steps_per_epoch=150)

history = model.fit(train_dataset, epochs=100, validation_data = test_dataset, callbacks=[checkpoint], verbose=2, steps_per_epoch=300)

result=[history.history['val_precision'],history.history['precision'],history.history['val_binary_accuracy'],history.history['binary_accuracy']]

model=keras.models.load_model(filepath)
np.savetxt('predDNN_0.dat',model.predict(X_test))
#np.savetxt('true.dat', y_test)
np.savetxt('result_DNN_0.dat',result)

headerList = ["mLL", "ptLL", "dPhi_LL", "dPhiLLmet", "MET", "mt", "jet_n", "goodjet_n", "goodbjet_n",
              "Pt1", "Eta1", "E1", "Phi1", "lep_charge1", "lep_type1", "lep_trackd0pvunbiased1", "lep_tracksigd0pvunbiased1", "lep_z01", "lep_ptcone301", "lep_etcone201", 
              "Pt2", "Eta2", "E2", "Phi2", "lep_charge2", "lep_type2", "lep_trackd0pvunbiased2", "lep_tracksigd0pvunbiased2", "lep_z02", "lep_ptcone302", "lep_etcone202", 
              "jeteta0", "jetMV2c100", "jetjvt0", "jetpt0", "jetphi0", "jetE0", "jeteta1", "jetMV2c101", "jetjvt1", "jetpt1", "jetphi1", "jetE1", "jeteta2", "jetMV2c102", "jetjvt2", "jetpt2", "jetphi2", "jetE2", "target"]
sel=np.arange(0,len(headerList)-1)

realdf=pd.read_csv('Data.csv', delimiter='\t', names=headerList)


realdf=realdf.loc[(realdf['ptLL']<400) & (realdf['MET']<600)& (realdf['Pt1']<150000) & (realdf['Pt2']<500000) & (realdf['MET']<250000) & (realdf['E1']<1000) & (realdf['lep_etcone201']>-10000) & (realdf['lep_tracksigd0pvunbiased1']<50) & (np.abs(realdf['lep_trackd0pvunbiased1'])<0.4) & (realdf['lep_tracksigd0pvunbiased1']<30) 
& (np.abs(realdf['lep_z01'])<2) & (realdf['lep_ptcone301']<10000) & (realdf['E2']<400) & (np.abs(realdf['lep_trackd0pvunbiased2'])<0.5) & (realdf['lep_tracksigd0pvunbiased2']<45) & (realdf['lep_ptcone302']<6000) & (np.abs(realdf['lep_etcone202'])<5000) & (realdf['jetpt1']<35000)]#ta bort sjuka outliars
realdf.to_numpy()
realdf=np.asarray(realdf)
realdata=realdf[:,sel]

for i in range(len(sel)): realdata[:,i]=(realdata[:,i]-np.mean(realdata[:,i]))/np.std(realdata[:,i]) 

pca = PCA(n_components=40)
realdataPCA = pca.fit_transform(realdata) #DNN
realoutputDNN=np.floor((model.predict(realdataPCA)).clip(min=0, max=1)).flatten()
np.savetxt('result_real_DNN_0.dat',realoutputDNN)
