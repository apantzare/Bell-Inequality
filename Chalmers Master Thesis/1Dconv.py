# -*- coding: utf-8 -*-
"""
Created on Sat Apr 22 21:09:21 2023

@author: axel
"""

#import xgboost as xgb
import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
import tensorflow as tf
from tensorflow import keras
from keras.layers import Conv1D
#from sklearn.decomposition import PCA
from tensorflow.keras.callbacks import ModelCheckpoint
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

newdata=data.copy()
newdataraw=np.c_[newdata, np.zeros((len(newdata[:,1]),1))]
newdata=np.c_[newdata, np.zeros((len(newdata[:,1]),1))]

div5plus0=[0, 5, 10, 15, 20, 25, 30, 35, 40, 45]
div5plus1=[1, 6, 11, 16, 21, 26, 31, 36, 41, 46]
div5plus2=[2, 7, 12, 17, 22, 27, 32, 37, 42, 47]
div5plus3=[3, 8, 13, 18, 23, 28, 33, 38, 43, 48]
div5plus4=[4, 9, 14, 19, 24, 29, 34, 39, 44, 49]


sel_chanal_1=[0, 1, 2, 3, 4, 5, 6, 7, 8, 30]
sel_chanal_2=[9, 10, 11, 12, 13, 14, 15, 16,17, 18]
sel_chanal_3=[20, 21, 22, 23, 24, 25, 26, 27, 28, 29]
sel_chanal_4=[31, 32, 33, 34, 35, 36, 37, 38, 39, 40]
sel_chanal_5=[41, 42, 43, 44, 45, 46, 47, 48, 49, 19]

newdata[:,div5plus0] = newdataraw[:,sel_chanal_1]
newdata[:,div5plus1] = newdataraw[:,sel_chanal_2]
newdata[:,div5plus2] = newdataraw[:,sel_chanal_3]
newdata[:,div5plus3] = newdataraw[:,sel_chanal_4]
newdata[:,div5plus4] = newdataraw[:,sel_chanal_5]

newdata=newdata.reshape(-1,10,5)
X_train, X_test, y_train, y_test = train_test_split(newdata, target, test_size=0.01, random_state=42)
train_dataset=tf.data.Dataset.from_tensor_slices((X_train, y_train))
test_dataset=tf.data.Dataset.from_tensor_slices((X_test, y_test))
train_dataset = train_dataset.batch(100)
test_dataset = test_dataset.batch(100)

drop=0.2
strategy=tf.distribute.MirroredStrategy()
with strategy.scope():
    model = keras.models.Sequential()
    model.add(Conv1D(128, kernel_size=3, activation='relu', input_shape=(10,5)))
    model.add(tf.keras.layers.MaxPooling1D(pool_size=2))
    model.add(Conv1D(256, kernel_size=3, activation='relu'))
    model.add(tf.keras.layers.MaxPooling1D(pool_size=2))
    model.add(tf.keras.layers.Flatten())
    model.add(tf.keras.layers.Dense(128, activation='relu'))
    model.add(tf.keras.layers.Dropout(drop))
    model.add(tf.keras.layers.Dense(1, activation='sigmoid'))
    
    model.compile(optimizer=keras.optimizers.Adam(learning_rate=0.00005),
                  loss=tf.keras.losses.BinaryCrossentropy(),
                  metrics=['BinaryAccuracy', 'Precision'])


filepath = '1Dconv_128+256+128_00005learning_.hdf5'

checkpoint = ModelCheckpoint(filepath, monitor='val_loss', verbose=1, \
                              save_best_only=True, save_weights_only=False, \
                             mode='auto', save_frequency=1)

history = model.fit(train_dataset, epochs=100, validation_data = test_dataset, callbacks=[checkpoint], verbose=2, steps_per_epoch=300)

result=[history.history['val_precision'],history.history['precision'],history.history['val_binary_accuracy'],history.history['binary_accuracy']]

model=keras.models.load_model(filepath)
np.savetxt('pred1Dconv.dat',model.predict(X_test))
#np.savetxt('true1Dconv.dat', y_test)
np.savetxt('result1Dconv.dat',result)

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

realdata1D=realdata.copy()
realdata1D=np.c_[realdata1D, np.zeros((len(realdata[:,1]),1))]
realdata1Draw=np.c_[realdata1D, np.zeros((len(realdata[:,1]),1))]

realdata1D[:,div5plus0] = realdata1Draw[:,sel_chanal_1]
realdata1D[:,div5plus1] = realdata1Draw[:,sel_chanal_2]
realdata1D[:,div5plus2] = realdata1Draw[:,sel_chanal_3]
realdata1D[:,div5plus3] = realdata1Draw[:,sel_chanal_4]
realdata1D[:,div5plus4] = realdata1Draw[:,sel_chanal_5]

realdata1D=realdata1D.reshape(-1,10,5)

realoutputDNN=model.predict(realdata1D)
np.savetxt('result_real_1Dconv.dat',realoutputDNN)
