import ROOT
import numpy as np
import csv
#import os
#import heapq

def filterAndSave(tree, name, single, HWW):
    with open(name+'weights.csv','a') as f1:
        writer=csv.writer(f1, delimiter='\t',lineterminator='\n')
        for entryNum in range(0 , tree.GetEntries()):
            tree.GetEntry(entryNum)

            if(all(getattr(tree, "lep_isTightID")) and (getattr( tree ,"trigE")==1 or getattr( tree ,"trigM")==1) and  
               getattr( tree ,"lep_type")[0] != getattr( tree ,"lep_type")[1] and getattr( tree ,"lep_charge")[0] != getattr( tree ,"lep_charge")[1] and 
               getattr( tree ,"lep_n")==2):   #exakt 2??? 
                 
                if(getattr( tree ,"met_et")>30000 and all(i > 15000 for i in getattr( tree ,"lep_pt")) and 
                 0.1 > (getattr( tree ,"lep_ptcone30")[0]/getattr( tree ,"lep_pt")[0]) and #standard lepton selection
                 0.1 > (getattr( tree ,"lep_ptcone30")[1]/getattr( tree ,"lep_pt")[1]) and#          
                 getattr( tree ,"lep_etcone20")[0]/getattr( tree ,"lep_pt")[0]<0.1 and #
                 getattr( tree ,"lep_etcone20")[1]/getattr( tree ,"lep_pt")[1]<0.1 and #
                 (np.abs(getattr( tree ,"lep_eta")[0])<2.5) and np.abs(getattr( tree ,"lep_eta")[1])<2.5):# 

                     lepton0 = ROOT.TLorentzVector()
                     lepton1 = ROOT.TLorentzVector()
                     pt = getattr( tree ,"lep_pt")  
                     eta = getattr( tree ,"lep_eta")
                     phi = getattr( tree ,"lep_phi")
                     nrg = getattr( tree ,"lep_E")
                     #MET = getattr( tree ,"met_et")
                     METPhi = getattr( tree ,"met_phi")
                     
                     lepton0.SetPtEtaPhiE ( pt[0] , eta[0] , phi[0] , nrg[0])
                     lepton1.SetPtEtaPhiE ( pt[1] , eta[1] , phi[1] , nrg[1])
                     dilepton = lepton0 + lepton1
                     mll = dilepton.M()
                     
                     dphi=phi[0]-phi[1]
                     if (dphi>np.pi):
                         dphi=dphi-2*np.pi
                     elif (dphi<-np.pi):
                         dphi=dphi+2*np.pi
                         
                     if (dilepton.Phi()-METPhi>np.pi/2 and dilepton.Pt()>30000 and (mll>10000 and mll<55000) and np.abs(dphi)<1.8):
                        
                         goodjets=0
                         goodbjets=0
                         for j in range(getattr(tree, "jet_n")):
                             if getattr( tree ,"jet_pt")[j]>30000:
                                 goodjets+=1
                                 if goodjets==2:
                                     break
                                 
                             if (getattr(tree, "jet_MV2c10")[j]>0.1758475 and getattr( tree ,"jet_pt")[j]>20000):
                                 goodbjets+=1
                                 break
                         if (goodjets<=1 and goodbjets==0):
                            #mt = np.sqrt(2*sum(pt)*MET*(1-np.cos(METPhi)))/1000
                            #jet_pt = np.asarray(getattr( tree ,"jet_pt"))
                            #jet_E = np.asarray(getattr( tree ,"jet_E"))
                            #jet_phi = np.asarray(getattr( tree ,"jet_phi"))
                            #jet_jvt = np.asarray(getattr( tree ,"jet_jvt"))
                            #jet_mv2c10 = np.asarray(getattr( tree ,"jet_MV2c10"))

                            #while len(jet_pt)<=2:
                            #    jet_pt=np.append(jet_pt,0)
                            #    jet_E=np.append(jet_E,0)
                            #    jet_phi=np.append(jet_phi,0)
                            #    jet_jvt=np.append(jet_jvt,0)
                            #    jet_mv2c10=np.append(jet_mv2c10,0)
                            #lep_etcone20 = getattr( tree ,"lep_etcone20")
                            #lep_ptcone30 = getattr( tree ,"lep_ptcone30")
                            #lep_z0 = getattr( tree ,"lep_z0")
                            #lep_charge = getattr( tree ,"lep_z0")
                            #lep_type = getattr( tree ,"lep_z0")
                            #lep_trackd0pvunbiased = getattr( tree ,"lep_trackd0pvunbiased")
                            #lep_etcone20 = getattr( tree ,"lep_etcone20")
                            #lep_ptcone30 = getattr( tree ,"lep_ptcone30")
                            #row=[pt[0],pt[1], eta[0], eta[1], phi[0], phi[1], dphi, nrg[0], nrg[1], mll, MET, METPhi, mt, getattr( tree ,"jet_n"),
                            #     jet_pt[0], jet_pt[1], jet_pt[2], jet_phi[0], jet_phi[1], jet_phi[2], jet_E[0], jet_E[1], jet_E[2], jet_jvt[0], jet_jvt[1], jet_jvt[2],
                            #     jet_mv2c10[0], jet_mv2c10[1], jet_mv2c10[2], lep_z0[0], lep_z0[1], lep_charge[0], lep_charge[1], lep_type[0], lep_type[1], 
                            #     lep_etcone20[0], lep_etcone20[1], lep_ptcone30[0], lep_ptcone30[1], lep_trackd0pvunbiased[0], lep_trackd0pvunbiased[1],
                            #     target, getattr( tree ,"mcWeight")] #sista endast för sim
                            row=[getattr( tree ,"scaleFactor_ELE"),getattr( tree ,"scaleFactor_MUON"),getattr( tree ,"scaleFactor_LepTRIGGER"),getattr( tree ,"mcWeight"),getattr( tree ,"scaleFactor_PILEUP"),getattr( tree ,"SumWeights"), single, HWW]
                            writer.writerow(row)
                            #if target==1:
                            #    totweight+=(getattr( tree ,"scaleFactor_ELE")*getattr( tree ,"scaleFactor_MUON")*getattr( tree ,"scaleFactor_LepTRIGGER")*np.sign(getattr( tree ,"mcWeight"))*getattr( tree ,"scaleFactor_PILEUP"))#
                            #else:
                            #   totweight+=(getattr( tree ,"scaleFactor_ELE")*getattr( tree ,"scaleFactor_MUON")*getattr( tree ,"scaleFactor_LepTRIGGER")*getattr( tree ,"mcWeight")*getattr( tree ,"scaleFactor_PILEUP"))
        #return totweight/getattr( tree ,"SumWeights")
                            
#-------------------------MAIN-----------------------------------
name='MC'
#os.remove(name+'.csv')
with open(name+'/'+'files.txt') as f: 
    parts = f.read().splitlines()
    for part in parts:
        if ('single' in part):# or ('345324' in part):
            single=1
        else:
            single=0
        
        if ('H125_WW' in part):# or ('345324' in part):
            HWW=1
        else:
            HWW=0
        filename=str(name+'/'+str(part))
        print(filename)#+ " Target="+str(target))
        inFile = ROOT.TFile.Open(filename, "Read")
        tree = inFile.Get ("mini")
        filterAndSave(tree, name,single,HWW)
        #totweight=filterAndSave(tree, name,target)
        #totaldata.append(totweight)
        #print((np.sum(totaldata)))
headerList = ['pt0','pt1', 'eta0', 'eta1', 'phi0', 'phi1', 'dphi', 'nrg0', 'nrg1', 'mll', 'MET', 'METPhi', 'mt', 'jet_n', 
              'jet_pt0', 'jet_pt1', 'jet_pt2', 'jet_phi0', 'jet_phi1', 'jet_phi2', 'jet_E0', 'jet_E1', 'jet_E2', 'jet_jvt0', 'jet_jvt1', 'jet_jvt2',
              'jet_mv2c100', 'jet_mv2c101',  'jet_mv2c102', 'lep_z00', 'lep_z01', 'lep_charge0', 'lep_charge1', 'lep_type0', 'lep_type1', 
              'lep_etcone200', 'lep_etcone201', 'lep_ptcone300', 'lep_ptcone301', 'lep_trackd0pvunbiased0', 'lep_trackd0pvunbiased1',
               'target', 'mcWeight']