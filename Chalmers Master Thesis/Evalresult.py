# -*- coding: utf-8 -*-
"""
Created on Sun Apr 23 17:42:10 2023

@author: axel
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn import metrics

def calcbell(phi1,phi2):
  Xp=np.cos(phi1)
  Yp=np.sin(phi1)
  Xm=np.cos(phi2)
  Ym=np.sin(phi2)
  
  term1=(np.multiply(Xp,Xm)+np.multiply(Yp,Ym)).flatten()
  term2=(np.multiply(np.multiply(Xp,Xp)-np.multiply(Yp,Yp),np.multiply(Xm,Xm)-np.multiply(Ym,Ym))).flatten()
  term3=np.multiply(np.multiply(Xp,Xm),np.multiply(Ym,Yp)).flatten()
  
  c1=8/np.sqrt(3)
  c2=25
  c3=100

  sum=c1*np.mean(term1)+c2*np.mean(term2)+c3*np.mean(term3)
  std=np.sqrt(c1**2*np.var(term1)+c2**2*np.var(term2)+c3**2*np.var(term3)+
  2*c1*c2*np.cov(term1,term2)[0,1]+2*c1*c3*np.cov(term1,term3)[0,1]+2*c2*c3*np.cov(term2,term3)[0,1])/np.sqrt(len(phi1))
  return [sum, std]

def resultfig(datafile, figfile):
    plt.rcParams.update({'font.size' : 18})
    fig,ax = plt.subplots(dpi=300)
    s=18
    ax2 = ax.twinx()
    ax.yaxis.label.set_color('m')
    ax.tick_params(axis='y', colors='m')
    ax2.yaxis.label.set_color('b')
    ax2.tick_params(axis='y', colors='b')

    result_DNN=np.genfromtxt(datafile)
    ax2.plot(result_DNN[0],label='Validationdata', linestyle='dotted' , color='b')
    ax2.plot(result_DNN[1], label='Trainingdata', linestyle='dashed', color='b')
    ax.plot(result_DNN[2],label='Validationdata', linestyle='dotted', color='m')
    ax.plot(result_DNN[3], label='Trainingdata', linestyle='dashed', color='m')
    
    ax2.spines['left'].set_color('m') 
    ax2.spines['right'].set_color('b') 
    ax.set_xlabel('Epochs', fontsize=s)
    ax.set_ylabel('Binary Accuracy', fontsize=s)
    ax2.set_ylabel('Precision', fontsize=s)
    ax.legend()

    plt.tight_layout()
    #fig.set_figwidth(5)
    #fig.set_figheight(4)
    fig.savefig(figfile)
    
def printresult(realoutput, true, pred, realdf):
    print("Totalt "+str(np.sum(realoutput))+" riktiga Higgs event")
    prec=round(100*metrics.precision_score(true,pred),2)
    acc=round(100*metrics.accuracy_score(true,pred),2)
    print("precision " + str(prec)+ "%, Accuracy "+ str(acc)+"%")
    phi1=np.extract(realoutput==1,realdf[:,12])
    phi2=np.extract(realoutput==1,realdf[:,23])
    bell=calcbell(phi1,phi2)
    print("For data Bell test shows "+str(round(bell[0],2))+"+-"+str(round(bell[1],2)))#+", with "+str(round((bell[0]-2)/bell[1],2))+" sigma certinty")
    return acc, prec, np.sum(realoutput), bell[0], bell[1]

resultfig("result1Dconv.dat", "train_1Dconv.png")
resultfig("result_DNN_0.dat", "train_DNN.png")

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

s=18
plt.rcParams.update({'font.size' : 18})
fig,ax = plt.subplots(dpi=300)
ax.hist(np.genfromtxt('result_real_DNN_0.dat'),bins=100)
ax.set_ylabel('Number of events',fontsize=s)
ax.set_xlabel('Probability of'+r' $H\rightarrow WW$',fontsize=s)
fig.tight_layout()
fig.savefig('real_DNN.png')
fig,ax = plt.subplots(dpi=300)
ax.hist(np.genfromtxt('result_real_1Dconv.dat'), bins=100)
ax.set_ylabel('Number of events',fontsize=s)
ax.set_xlabel('Probability of'+r' $H\rightarrow WW$',fontsize=s)
fig.tight_layout()
fig.savefig('real_1D.png')


rond=0.9
realoutputDNN = np.floor((np.genfromtxt('result_real_DNN_0.dat')/rond).clip(min=0, max=1)).flatten()
realoutput1Dconv = np.floor((np.genfromtxt('result_real_1Dconv.dat')/rond).clip(min=0, max=1)).flatten()

predDNN = np.floor((np.genfromtxt('predDNN_0.dat')/rond).clip(min=0, max=1)).flatten()
pred1D = np.floor((np.genfromtxt('pred1Dconv.dat')/rond).clip(min=0, max=1)).flatten()
true = np.genfromtxt('true.dat')


accuracy=[79,79, 59] #BCE
precision=[64,67, 33]
events=[4308,5192, 17394]
bell=[4.32,2.96, 6.44]
error=[0.27,0.23, 0.12]

accuracy[0], precision[0], events[0], bell[0], error[0] = printresult(realoutputDNN, true, predDNN, realdf)
accuracy[1], precision[1], events[1], bell[1], error[1] = printresult(realoutput1Dconv, true, pred1D, realdf)
accuracy[2], precision[2], events[2], bell[2], error[2] = printresult(realoutput1Dconv*realoutputDNN, true, predDNN*pred1D, realdf)

nummer=np.arange(1,len(accuracy)+1)
fig,ax = plt.subplots(dpi=300)
s=18
plt.rcParams.update({'font.size' : 18})

ax2 = ax.twinx()
ax.bar(nummer-0.2, accuracy,0.2, label='Accuracy')
ax.bar(nummer, precision,0.2, label='Precision')
ax2.bar(nummer+0.2, events,0.2, label='#Events', color='m')
ax.legend(loc='upper right')
ax.axis([0.5,len(accuracy)+0.5,80,100])
ax2.axis([0.5,len(accuracy)+0.5,0,11100])
X=['DNN', '1D', 'DNN+1D']
ax.set_xticks(nummer, X, rotation = 45, fontsize=s)

ax.set_ylabel('Accuracy/precision(%)',fontsize=s)
ax.set_xlabel('Model',fontsize=s)

ax2.yaxis.label.set_color('m')
ax2.set_ylabel('#Predicted events', color='m')
ax2.tick_params(axis='y', colors='m')
ax2.spines['right'].set_color('m') 
ax2.axhline(423,color='m', linestyle='dotted')
ax.grid(axis='y')
fig.tight_layout()
fig.savefig("acc0.9.png")

fig,ax = plt.subplots(dpi=300)
plt.rcParams.update({'font.size' : 18})
plt.errorbar(nummer, bell, error,fmt='o', alpha=0.5,color='blue', label='value $\pm \sigma$')
plt.hlines(2,nummer[0]-0.5, nummer[-1]+0.5, color='red', label='Bell limit', linestyle='dotted')
plt.fill_between(2*nummer-2,2.3,4.1,alpha=0.4, label='predicted value, Training data')
plt.fill_between(2*nummer-2,2.7,2.8, label='predicted value, Barr-Fabriecchi')
plt.axis([0.5,nummer[-1]+0.5,1,11])
plt.xlabel('Model', fontsize=s)
ax.set_xticks(nummer, X, rotation = 45, fontsize=s)
plt.ylabel('CGLMP value', fontsize=s)
plt.legend(loc='upper right', fontsize=s-3)
plt.tight_layout()
plt.grid(axis='y')
fig.savefig("bell09.png")

rond=0.5
realoutputDNN = np.floor((np.genfromtxt('result_real_DNN_0.dat')/rond).clip(min=0, max=1)).flatten()
realoutput1Dconv = np.floor((np.genfromtxt('result_real_1Dconv.dat')/rond).clip(min=0, max=1)).flatten()
predDNN = np.floor((np.genfromtxt('predDNN_0.dat')/rond).clip(min=0, max=1)).flatten()
pred1D = np.floor((np.genfromtxt('pred1Dconv.dat')/rond).clip(min=0, max=1)).flatten()

accuracy[0], precision[0], events[0], bell[0], error[0] = printresult(realoutputDNN, true, predDNN, realdf)
accuracy[1], precision[1], events[1], bell[1], error[1] = printresult(realoutput1Dconv, true, pred1D, realdf)
accuracy[2], precision[2], events[2], bell[2], error[2] = printresult(realoutput1Dconv*realoutputDNN, true, predDNN*pred1D, realdf)

nummer=np.arange(1,len(accuracy)+1)
fig,ax = plt.subplots(dpi=300)
s=18
plt.rcParams.update({'font.size' : 18})

ax2 = ax.twinx()
ax.bar(nummer-0.2, accuracy,0.2, label='Accuracy')
ax.bar(nummer, precision,0.2, label='Precision')
ax2.bar(nummer+0.2, events,0.2, label='#Events', color='m')
ax.legend(loc='upper right')
ax.axis([0.5,len(accuracy)+0.5,80,100])
ax2.axis([0.5,len(accuracy)+0.5,0,11100])
X=['DNN', '1D', 'DNN+1D']
ax.set_xticks(nummer, X, rotation = 45, fontsize=s)

ax.set_ylabel('Accuracy/precision(%)',fontsize=s)
ax.set_xlabel('Model',fontsize=s)

ax2.yaxis.label.set_color('m')
ax2.set_ylabel('#Predicted events', color='m')
ax2.tick_params(axis='y', colors='m')
ax2.spines['right'].set_color('m') 
ax2.axhline(423,color='m', linestyle='dotted')
ax.grid(axis='y')
fig.tight_layout()
fig.savefig("acc.png")

fig,ax = plt.subplots(dpi=300)
plt.rcParams.update({'font.size' : 18})
plt.errorbar(nummer, bell, error,fmt='o', alpha=0.5,color='blue', label='value $\pm \sigma$')
plt.hlines(2,nummer[0]-0.5, nummer[-1]+0.5, color='red', label='Bell limit', linestyle='dotted')
plt.fill_between(2*nummer-2,2.3,4.1,alpha=0.4, label='predicted value, Training data')
plt.fill_between(2*nummer-2,2.7,2.8, label='predicted value, Barr-Fabriecchi')
plt.axis([0.5,nummer[-1]+0.5,1,11])
plt.xlabel('Model', fontsize=s)
ax.set_xticks(nummer, X, rotation = 45, fontsize=s)
plt.ylabel('CGLMP value', fontsize=s)
plt.legend(loc='upper right', fontsize=s-3)
plt.tight_layout()
plt.grid(axis='y')
fig.savefig("bell.png")