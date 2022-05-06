# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from scipy.optimize import curve_fit

sigma=1

Dataterm1=[2.81781,3.10438,2.97804,2.77163]
Dataterm1error=[0.03202,0.03248,0.040248,0.040248]
Dataterm1error = [sigma*x for x in Dataterm1error]

Dataterm2=[0.124356,1.51799,0.667284,-0.386129]
Dataterm2error=[0.25436,0.28177,0.33922,0.38523]
Dataterm2error = [sigma*x for x in Dataterm2error]

Dataterm3=[-0.597703,0.889404,0.43801,-0.760452]
Dataterm3error=[0.25098,0.27646,0.32904,0.37496]
Dataterm3error = [sigma*x for x in Dataterm3error]

DataTot=[2.34447, 5.51177, 4.08333, 1.62504]
DataToterror=[0.38290,0.41069,0.49653,0.57182]
DataToterror = [sigma*x for x in DataToterror]


Higgsterm1=[2.95541, 3.20058, 3.21059, 3.10433]
Higgsterm1error=[0.00452, 0.004357, 0.00448, 0.00492]
Higgsterm1error = [sigma*x for x in Higgsterm1error]

Higgsterm2=[0.530123, 1.71271, 1.77822, 1.11786]
Higgsterm2error=[0.03618, 0.0382, 0.039267, 0.04218]
Higgsterm2error = [sigma*x for x in Higgsterm2error]

Higgsterm3=[0.553257, 1.74741, 1.8284, 1.16906]
Higgsterm3error=[0.03631, 0.03844, 0.03957, 0.042516]
Higgsterm3error = [sigma*x for x in Higgsterm3error]

HiggsTot=[4.03879, 6.6607, 6.81722, 5.39125]
HiggsToterror=[0.05516, 0.05663, 0.05823, 0.06305]
HiggsToterror = [sigma*x for x in HiggsToterror]


Bkgterm1=[2.79033, 3.14958, 3.02996, 2.8338]
Bkgterm1error=[0.0042, 0.00434, 0.00522, 0.00611]
Bkgterm1error = [sigma*x for x in Bkgterm1error]

Bkgterm2=[-0.451217, 1.36012, 0.730382, -0.372929]
Bkgterm2error=[0.033, 0.038, 0.0439, 0.04954]
Bkgterm2error = [sigma*x for x in Bkgterm2error]

Bkgterm3=[-0.432993, 1.39749, 0.781717, -0.36142]
Bkgterm3error=[0.033, 0.03787, 0.04376, 0.04932]
Bkgterm3error = [sigma*x for x in Bkgterm3error]

BkgTot=[1.90612, 5.90719, 4.54206, 2.09945]
BkgToterror=[0.05048, 0.05613, 0.06576, 0.07492]
BkgToterror = [sigma*x for x in BkgToterror]


Simterm1=[2.86501, 3.1748, 3.1289, 2.98837]
Simterm1error=[0.00308,0.00307,0.00341, 0.00386]
Simterm1error = [sigma*x for x in Simterm1error]

Simterm2=[-0.00724,1.53475, 1.30445, 0.47884]
Simterm2error=[0.02441, 0.02698, 0.02929, 0.03218]
Simterm2error = [sigma*x for x in Simterm2error]

Simterm3=[0.0132,1.5708, 1.35515, 0.51303]
Simterm3error=[0.0244, 0.02698, 0.02937, 0.032265]
Simterm3error = [sigma*x for x in Simterm3error]

SimTot=[2.87098, 6.2804, 5.78851, 3.98024]
SimToterror=[0.037295, 0.03987, 0.04369, 0.04845]
SimToterror = [sigma*x for x in SimToterror]

errorpred=[2.05091,2.01701, 1.79786,1.80836]

Sig=[117107/np.sqrt(117107+141739),100050/np.sqrt(100050+101952),94722/np.sqrt(94722+78174),83013/np.sqrt(83013+62278)]
N=[2463, 1951, 1402, 1094]
SB=[117107/(117107+141739),100050/(100050+101952),94722/(94722+78174),83013/(83013+62278)]
SBData=[0.07,0.088,0.12,0.15]
NHiggs=np.multiply(N,SBData)
Constraint=[0,1,2,3]

fig, ax1 = plt.subplots()

ax2 = ax1.twinx()
ax1.plot(Constraint, SB, 'r-')
ax2.plot(Constraint, Sig, 'b-')

ax1.xaxis.set_major_locator(MaxNLocator(integer=True))

ax1.set_xlabel('Number of constraints')
ax1.set_ylabel('S/(S+B)', color='r')
ax2.set_ylabel('Signal significance', color='b')

plt.show()

def specialplot(x,y1,y2,y3,y4,error1,error2,error3,error4, head):
    fig, ax1 = plt.subplots()

    plt.scatter(x,y1,c="g", alpha=0.8, marker='D',
            label="Data")
    plt.scatter(x,y2,c="b", alpha=0.8, marker='o',
            label="Higgs")
    plt.scatter(x,y3,c="r", alpha=0.8, marker='*',
            label="Background")
    plt.scatter(x,y4,c="0", alpha=0.8, marker='s',
            label="Complete simulation")
    plt.xlabel("Constraints")
    plt.ylabel("I value")
    plt.legend(loc='lower center')
    plt.title(head)
    ax1.xaxis.set_major_locator(MaxNLocator(integer=True))
    
    plt.errorbar(x, y1, yerr=error1, fmt="D", color="g")
    plt.errorbar(x, y2, yerr=error2, fmt="o", color="b")
    plt.errorbar(x, y3, yerr=error3, fmt="*", color="r")
    plt.errorbar(x, y4, yerr=error4, fmt="s", color="0")
    plt.show()
    
def func(x,a,c):    
    return a*x+c

def errorfunc(x,a,c):    
    return a*np.sqrt(x)+c

def Higgsfunc(x,a,b):    
    return a*np.exp(b*x)

def specialplot2(SB,xdata,ypure,ybkg,ysim,y2,error1pure,error1bkg,error1sim,error2,head):
    fig, ax1 = plt.subplots()
    x=np.linspace(0,1,100)
    k=[]
    for i in range(0,len(SB)):
        x1=[0,SB[i],1]
        y1=[ybkg[i],ysim[i],ypure[i]]
        error1=[error1bkg[i],error1sim[i],error1pure[i]]
        if i==1:
            col='b'
        elif i==2:
            col='y'
        elif i==3:
            col='m'
        elif i==0:
            col='k'
        plt.scatter(x1,y1,c=col, alpha=0.8, marker='s')#,
                    #label="Complete simulation "+str(i)+" constraints")
        plt.errorbar(x1, y1, yerr=error1, fmt="s", color=col)
        plt.errorbar(xdata[i], y2[i], yerr=error2[i], fmt="D", color=col)
        popt, pcov = curve_fit(func, x1, y1)
        yfit=func(x,*popt) 
        k.append(popt[0])
        style=col+'-'
        plt.plot(x, yfit, style,
         label=str(i)+" constraints")#'fit: a=%5.3f, c=%5.3f' % tuple(popt))
    
        plt.scatter(xdata[i],y2[i],c=col, alpha=0.8, marker='D')#,
                    #label="Data")
    

    plt.xlabel("S/(S+B)")
    plt.ylabel("I value")
    plt.legend(loc='lower right')
    plt.title(head)
    plt.show()   
    return k

specialplot(Constraint,DataTot,HiggsTot,BkgTot,SimTot,DataToterror,HiggsToterror,BkgToterror,SimToterror, "Total I")
specialplot(Constraint,Dataterm1,Higgsterm1,Bkgterm1,Simterm1,Dataterm1error,Higgsterm1error,Bkgterm1error, Simterm1error, "Term 1 I")
specialplot(Constraint,Dataterm2,Higgsterm2,Bkgterm2,Simterm2,Dataterm2error,Higgsterm2error,Bkgterm2error, Simterm2error, "Term 2 I")
specialplot(Constraint,Dataterm3,Higgsterm3,Bkgterm3,Simterm3,Dataterm3error,Higgsterm3error,Bkgterm3error, Simterm3error, "Term 3 I")

klist=specialplot2(SB,SBData,HiggsTot,BkgTot,SimTot, DataTot,HiggsToterror,BkgToterror,SimToterror, DataToterror, "Comparison simulation and data, Tot")
specialplot2(SB,SBData,Higgsterm1,Bkgterm1,Simterm1, Dataterm1,Higgsterm1error,Bkgterm1error,Simterm1error, Dataterm1error, "Comparison simulation and data, term 1")
specialplot2(SB,SBData,Higgsterm2,Bkgterm2,Simterm2, Dataterm2,Higgsterm2error,Bkgterm2error,Simterm2error, Dataterm2error, "Comparison simulation and data, term 2")
specialplot2(SB,SBData,Higgsterm3,Bkgterm3,Simterm3, Dataterm3,Higgsterm3error,Bkgterm3error,Simterm3error, Dataterm3error, "Comparison simulation and data, term 3")

plt.scatter(SBData,DataToterror,label="Data")
x=np.linspace(0,1,100)
popt, pcov = curve_fit(func, SBData, DataToterror)
yfit=func(x,*popt) 
plt.plot(x, yfit,label="Model")
plt.xlabel("S/(S+B)")
plt.ylabel("I error")
plt.legend(loc='lower right')
plt.title("Extrapolation of error")
plt.show()

Maxerror=func(1,*popt)
meanI=np.mean(DataTot)
meanX=np.mean(SBData)
meanError=np.mean(DataToterror)

plt.scatter(SBData,DataTot,c='b', alpha=0.8, marker='o',label="Data")
plt.errorbar(SBData, DataTot, yerr=DataToterror,fmt="o", color='b')
k=np.average(klist)
Y=[]
X=[1,1,1,1]
for i in range(0,4):
    Y.append(klist[i]*(1-SBData[i])+DataTot[i])
plt.scatter(X,Y,label="Predicted",c='y',marker='D')
plt.errorbar(X, Y, yerr=errorpred,c='y',fmt='D')
for i in range(0,4):
    plt.plot(x,klist[i]*x+(Y[i]-klist[i]*X[i]),label="Model with "+str(i) +" constraint")
plt.xlabel("S/(S+B)")
plt.ylabel("Predicted I value")
plt.axhline(y=2, color='r', linestyle='-', label="Bell limit")
plt.legend(loc='upper center')
plt.title("Extrapolation of Total I value")
plt.show()

plt.scatter(SBData, NHiggs,label="Data")
popt, pcov = curve_fit(func, SBData, NHiggs)
yfit=func(x,*popt)
plt.plot(x,yfit,label="model")
plt.xlabel("S/(S+B)")
plt.ylabel("Number of Higgs")
plt.legend(loc='upper right')
plt.title("Extrapolation of Number of Higgs events")
plt.show()

plt.scatter(SBData,N,label="Data")
x=np.linspace(0,1,100)
popt, pcov = curve_fit(Higgsfunc, SBData, N)
yfit=Higgsfunc(x,*popt) 
plt.plot(x, yfit,label="Model")
plt.xlabel("S/(S+B)")
plt.ylabel("Number of events")
plt.legend(loc='upper right')
plt.title("Extrapolation of number of events")
plt.ylim(0, 2500)
plt.show()