//////////////////////////////////////////////////////////////////////////////////////////
// Creates a TChain to be used by the Analysis.C class
#include "TROOT.h"
#include "TChain.h"
#include "TFile.h"
#include "TProof.h"

void main_HWWAnalysis(int proof = 0, int option= 0)
{
  // path to your local directory *or* URL, please change the default one!
  /* Local path example */
  //TString path = "/eos/project/a/atlas-outreach/projects/open-data/OpenDataTuples/renamedLargeRJets/2lep/";

  /* The URL to the ATLAS Open Data website repository */
  //TString path = "https://atlas-opendata.web.cern.ch/atlas-opendata/samples/2020/2lep/";
  
  /* The URL to the CERN Open Data portal repository */
  TString path = "~/GIT-REPOS/atlas-outreach-cpp-framework-13tev/Analysis/HWWAnalysis/";
  
  /* The XROOTD path to the CERN Open Data portal repository */
  //TString path = "root://eospublic.cern.ch//eos/opendata/atlas/OutreachDatasets/2020-01-22/2lep/";
  
  /* An external mirror */
  //TString path = "http://www.universidad.ch/ATLAS/outreach/open-data-2020/2lep/";

  //***************************************************************************************************//
  // adding chains of all MC and data samples
  // option==0 will run all one by one
  // Currently 5 options for MC (2,3,4,5,6) and 4 for data (11,12,13,14) which can be run in parallel
  // If using the options (11,12,13,14) of splitting data, 
  // you need to add the samples later with: hadd data.root dataA.root dataB.root dataC.root dataD.root
  //***************************************************************************************************//
  if (proof == 1)  TProof::Open("");
  
  if (option==10 || option==0){
  TChain* chain_data = new TChain("mini");
  chain_data->AddFile(path+"Data/data.2lep.root");
  if (proof == 1)  chain_data->SetProof();
  chain_data->Process("HWWAnalysis.C+","data");
  }
  
  if (option==1){
  TChain* chain_sim = new TChain("mini");
  chain_sim->AddFile(path+"MC/simulation.root");
  if (proof == 1)  chain_sim->SetProof();
  chain_sim->Process("HWWAnalysis.C+","simulation");
  }  
  
  if (option==32){
  TChain* chain_pure = new TChain("mini");
  chain_pure->AddFile(path+"MC/Pure.root");
  if (proof == 1)  chain_pure->SetProof();
  chain_pure->Process("HWWAnalysis.C+","Pure");
  }   
  
  if (option==31){
  TChain* chain_background = new TChain("mini");
  chain_background->AddFile(path+"MC/background.root");
  if (proof == 1)  chain_background->SetProof();
  chain_background->Process("HWWAnalysis.C+","background");
  }  
  
  if (option==30){
  TChain* chain_simulation = new TChain("mini");
  chain_simulation->AddFile(path+"MC/simulation.root");
  if (proof == 1)  chain_simulation->SetProof();
  chain_simulation->Process("HWWAnalysis.C+","simulation");
  }         
  
  //////////////////////  MC samples
  
  if (option==2 || option==0)   {
    
    // diboson
    TChain* chain_ZqqZll = new TChain("mini");
    chain_ZqqZll->AddFile(path+"MC/mc_363356.ZqqZll.2lep.root");
    if (proof == 1)  chain_ZqqZll->SetProof();
    chain_ZqqZll->Process("HWWAnalysis.C+","ZqqZll");
    
    TChain* chain_WqqZll = new TChain("mini");
    chain_WqqZll->AddFile(path+"MC/mc_363358.WqqZll.2lep.root");
    if (proof == 1)  chain_WqqZll->SetProof();
    chain_WqqZll->Process("HWWAnalysis.C+","WqqZll");
    
    TChain* chain_WpqqWmlv = new TChain("mini");
    chain_WpqqWmlv->AddFile(path+"MC/mc_363359.WpqqWmlv.2lep.root");
    if (proof == 1)  chain_WpqqWmlv->SetProof();
    chain_WpqqWmlv->Process("HWWAnalysis.C+","WpqqWmlv");
    
    TChain* chain_WplvWmqq = new TChain("mini");
    chain_WplvWmqq->AddFile(path+"MC/mc_363360.WplvWmqq.2lep.root");
    if (proof == 1)  chain_WplvWmqq->SetProof();
    chain_WplvWmqq->Process("HWWAnalysis.C+","WplvWmqq");
    
    TChain* chain_WlvZqq = new TChain("mini");
    chain_WlvZqq->AddFile(path+"MC/mc_363489.WlvZqq.2lep.root");
    if (proof == 1)  chain_WlvZqq->SetProof();
    chain_WlvZqq->Process("HWWAnalysis.C+","WlvZqq");
  }

  if (option==3 || option==0)   {

    TChain* chain_llll = new TChain("mini");
    chain_llll->AddFile(path+"MC/mc_363490.llll.2lep.root");
    if (proof == 1)  chain_llll->SetProof();
    chain_llll->Process("HWWAnalysis.C+","llll");
    
    TChain* chain_lllv = new TChain("mini");
    chain_lllv->AddFile(path+"MC/mc_363491.lllv.2lep.root");
    if (proof == 1)  chain_lllv->SetProof();
    chain_lllv->Process("HWWAnalysis.C+","lllv");
    
    TChain* chain_llvv = new TChain("mini");
    chain_llvv->AddFile(path+"MC/mc_363492.llvv.2lep.root");
    if (proof == 1)  chain_llvv->SetProof();
    chain_llvv->Process("HWWAnalysis.C+","llvv");
    
    TChain* chain_lvvv = new TChain("mini");
    chain_lvvv->AddFile(path+"MC/mc_363493.lvvv.2lep.root");
    if (proof == 1)  chain_lvvv->SetProof();
    chain_lvvv->Process("HWWAnalysis.C+","lvvv");
    
    //single top
    TChain* chain_single_top_tchan = new TChain("mini");
    chain_single_top_tchan->AddFile(path+"MC/mc_410011.single_top_tchan.2lep.root");
    if (proof == 1)  chain_single_top_tchan->SetProof();
    chain_single_top_tchan->Process("HWWAnalysis.C+","single_top_tchan");
    
    TChain* chain_single_antitop_tchan = new TChain("mini");
    chain_single_antitop_tchan->AddFile(path+"MC/mc_410012.single_antitop_tchan.2lep.root");
    if (proof == 1)  chain_single_antitop_tchan->SetProof();
    chain_single_antitop_tchan->Process("HWWAnalysis.C+","single_antitop_tchan");
    
    TChain* chain_single_top_schan = new TChain("mini");
    chain_single_top_schan->AddFile(path+"MC/mc_410025.single_top_schan.2lep.root");
    if (proof == 1)  chain_single_top_schan->SetProof();
    chain_single_top_schan->Process("HWWAnalysis.C+","single_top_schan");
    
    TChain* chain_single_antitop_schan = new TChain("mini");
    chain_single_antitop_schan->AddFile(path+"MC/mc_410026.single_antitop_schan.2lep.root");
    if (proof == 1)  chain_single_antitop_schan->SetProof();
    chain_single_antitop_schan->Process("HWWAnalysis.C+","single_antitop_schan");
    
    TChain* chain_single_top_wtchan = new TChain("mini");
    chain_single_top_wtchan->AddFile(path+"MC/mc_410013.single_top_wtchan.2lep.root");
    if (proof == 1)  chain_single_top_wtchan->SetProof();
    chain_single_top_wtchan->Process("HWWAnalysis.C+","single_top_wtchan");
    
    TChain* chain_single_antitop_wtchan = new TChain("mini");
    chain_single_antitop_wtchan->AddFile(path+"MC/mc_410014.single_antitop_wtchan.2lep.root");
    if (proof == 1)  chain_single_antitop_wtchan->SetProof();
    chain_single_antitop_wtchan->Process("HWWAnalysis.C+","single_antitop_wtchan");

 }


 if (option==4 || option==0)   {

    // Z+jets inclusive
 
    TChain* chain_Zee = new TChain("mini");
    chain_Zee->AddFile(path+"MC/mc_361106.Zee.2lep.root");
    if (proof == 1)  chain_Zee->SetProof();
    chain_Zee->Process("HWWAnalysis.C+","Zee");
 }

 if (option==5 || option==0)   {

    TChain* chain_Zmumu = new TChain("mini");
    chain_Zmumu->AddFile(path+"MC/mc_361107.Zmumu.2lep.root");
    if (proof == 1)  chain_Zmumu->SetProof();
    chain_Zmumu->Process("HWWAnalysis.C+","Zmumu");
 }


 if (option==6 || option==0)   {

    TChain* chain_Ztautau = new TChain("mini");
    chain_Ztautau->AddFile(path+"MC/mc_361108.Ztautau.2lep.root");
    if (proof == 1)  chain_Ztautau->SetProof();
    chain_Ztautau->Process("HWWAnalysis.C+","Ztautau");

    //ttbar
    TChain* chain_ttbar_lep = new TChain("mini");
    chain_ttbar_lep->AddFile(path+"MC/mc_410000.ttbar_lep.2lep.root");
    if (proof == 1)  chain_ttbar_lep->SetProof();
    chain_ttbar_lep->Process("HWWAnalysis.C+","ttbar_lep");
    
  // W+jets inclusive
  TChain* chain_Wplusenu = new TChain("mini");
  chain_Wplusenu->AddFile(path+"MC/mc_361100.Wplusenu.2lep.root");
  if (proof == 1)  chain_Wplusenu->SetProof();
  chain_Wplusenu->Process("HWWAnalysis.C+","Wplusenu");

  TChain* chain_Wplusmunu = new TChain("mini");
  chain_Wplusmunu->AddFile(path+"MC/mc_361101.Wplusmunu.2lep.root");
  if (proof == 1)  chain_Wplusmunu->SetProof();
  chain_Wplusmunu->Process("HWWAnalysis.C+","Wplusmunu");

  TChain* chain_Wplustaunu = new TChain("mini");
  chain_Wplustaunu->AddFile(path+"MC/mc_361102.Wplustaunu.2lep.root");
  if (proof == 1)  chain_Wplustaunu->SetProof();
  chain_Wplustaunu->Process("HWWAnalysis.C+","Wplustaunu");

  TChain* chain_Wminusenu = new TChain("mini");
  chain_Wminusenu->AddFile(path+"MC/mc_361103.Wminusenu.2lep.root");
  if (proof == 1)  chain_Wminusenu->SetProof();
  chain_Wminusenu->Process("HWWAnalysis.C+","Wminusenu");

  TChain* chain_Wminusmunu = new TChain("mini");
  chain_Wminusmunu->AddFile(path+"MC/mc_361104.Wminusmunu.2lep.root");
  if (proof == 1)  chain_Wminusmunu->SetProof();
  chain_Wminusmunu->Process("HWWAnalysis.C+","Wminusmunu");

  TChain* chain_Wminustaunu = new TChain("mini");
  chain_Wminustaunu->AddFile(path+"MC/mc_361105.Wminustaunu.2lep.root");
  if (proof == 1)  chain_Wminustaunu->SetProof();
  chain_Wminustaunu->Process("HWWAnalysis.C+","Wminustaunu");

  // ggH
  TChain* chain_ggH125_WW = new TChain("mini");
  chain_ggH125_WW->AddFile(path+"MC/mc_345324.ggH125_WW2lep.2lep.root");
  if (proof == 1)  chain_ggH125_WW->SetProof();
  chain_ggH125_WW->Process("HWWAnalysis.C+","ggH125_WW2lep");

  // VBF
  TChain* chain_VBFH125_WW = new TChain("mini");
  chain_VBFH125_WW->AddFile(path+"MC/mc_345323.VBFH125_WW2lep.2lep.root");
  if (proof == 1)  chain_VBFH125_WW->SetProof();
  chain_VBFH125_WW->Process("HWWAnalysis.C+","VBFH125_WW2lep");


 } // option 6

  
}
