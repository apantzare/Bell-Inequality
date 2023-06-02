# -*- coding: utf-8 -*-
"""
Created on Mon Apr 24 11:19:32 2023

@author: axel
"""

import numpy as np
from sklearn.decomposition import PCA
import pandas as pd

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
print(realdataPCA)
