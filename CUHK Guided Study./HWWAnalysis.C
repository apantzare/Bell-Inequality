/////////////////////////////////////////////////////////////
//// HWWAnalysis code
//// Author: ATLAS Collaboration (2019)
////
////
//// DISCLAIMER:
//// This Software is intended for educational use only!
//// Under no circumstances does it qualify to reproduce actual ATLAS analysis results or produce publishable results!
/////////////////////////////////////////////////////////////

#define HWWAnalysis_cxx
#include "TROOT.h"
#include "TPaveLabel.h"
#include "HWWAnalysis.h"
#include "HWWAnalysisHistograms.h"
#include <iostream>
#include <cstring>
#include <TAttMarker.h>
#include <string>
#include <vector>
#include <numeric>

#include <TH1.h>
#include <TH2.h>
#include <TStyle.h>
#include <TMath.h>
#include <TLorentzVector.h>

#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <algorithm>
#include <math.h>
#include <TF1.h>
#include <TH1.h>
#include <TLegend.h>
#include <TFile.h>
#include <TDirectory.h>
#include <TTree.h>
#include <TEnv.h>
#include "THStack.h"
#include "TFile.h"
#include "TTree.h"
#include "TKey.h"
#include "TCanvas.h"
#include <TStyle.h>
#include "TLatex.h"
#include "TImage.h"
#include "TLine.h"
#include "TColor.h"
#include "TROOT.h"
#include "TH2F.h"
#include "TMath.h"
#include "TPaveText.h" 
#include "TGraph.h"
string name;
// declaring variables for plot later
std::vector< float > Iarrterm1;
std::vector< float > Iarrterm2;
std::vector< float > Iarrterm3;
std::vector< float > xmLL;
std::vector< float > xpLLmet;
std::vector< float > xpt;
std::vector< float > xMET;
std::vector< float > xpLL;
std::vector< float > xiYminvec;
std::vector< float > xiYplusvec;
std::vector< float > xiXminvec;
std::vector< float > xiXplusvec;
void HWWAnalysis::Begin(TTree * )
{

  nEvents=0;

}

void HWWAnalysis::SlaveBegin(TTree * )
{
  TString option = GetOption();
  printf("Starting analysis with process option: %s \n", option.Data());
  
  name=option;
  
  define_histograms();
  
  FillOutputList();
}

Bool_t HWWAnalysis::Process(Long64_t entry)
{
  fChain->GetTree()->GetEntry(entry);
  nEvents++;
  if (nEvents % 1000000 == 0) std::cout << "Analysed a total of: " << nEvents << " events out of " << fChain->GetTree()->GetEntries() << " in this sample" << std::endl;
  
  if(fChain->GetTree()->GetEntries()>0)
    {
      // **********************************************************************************************//
      // Begin analysis selection, largely based on: ATLAS Collaboration, Phys. Lett. B 789 (2019) 508 //
      // **********************************************************************************************//

      //Scale factors
      Float_t scaleFactor = scaleFactor_ELE*scaleFactor_MUON*scaleFactor_LepTRIGGER*scaleFactor_PILEUP;

      //MC weight
      Float_t m_mcWeight = mcWeight;

      // read input option
      TString option = GetOption();
      if(option.Contains("single")) { m_mcWeight = (mcWeight/TMath::Abs(mcWeight)); } // set to 1 or -1 for single top MCs

      //Total weight
      Float_t weight = scaleFactor*m_mcWeight;

      // Make difference between data and MC
      if(option.Contains("data")) {  weight = 1.; }
      
      
      //Preselection cut for electron/muon trigger 
      if(trigE || trigM)
	{
	  
	  // Preselection of good leptons
	  int goodlep_index[2];
	  int goodlep_n = 0;
	  int lep_index =0;
	  
	  for(unsigned int i=0; i<lep_n; i++)
	    {

              TLorentzVector leptemp;  leptemp.SetPtEtaPhiE(lep_pt->at(i)/1000., lep_eta->at(i), lep_phi->at(i), lep_E->at(i)/1000.);
	      
	      // Lepton is Tight
	      if( lep_isTightID->at(i) )
	        {
		  
		  // standard lepton isolation requirement => strict isolation
		  if( lep_pt->at(i) >15000. && ( (lep_ptcone30->at(i)/lep_pt->at(i)) < 0.1) && ( (lep_etcone20->at(i) / lep_pt->at(i)) < 0.1 ) )
		    {
		      if ( lep_type->at(i)==11 && TMath::Abs(lep_eta->at(i)) < 2.47 && ( TMath::Abs(lep_eta->at(i)) < 1.37 || TMath::Abs(lep_eta->at(i)) > 1.52 ) ) {
			if( TMath::Abs(lep_trackd0pvunbiased->at(i))/lep_tracksigd0pvunbiased->at(i) < 5 && TMath::Abs(lep_z0->at(i)*TMath::Sin(leptemp.Theta())) < 0.5) {
			  goodlep_n = goodlep_n + 1;
			  goodlep_index[lep_index] = i;
			  lep_index++;
			}
		      }
		      // muon selection
		      if ( lep_type->at(i) ==13 && TMath::Abs(lep_eta->at(i)) < 2.5 ) {
			if( TMath::Abs(lep_trackd0pvunbiased->at(i))/lep_tracksigd0pvunbiased->at(i) < 3 && TMath::Abs(lep_z0->at(i)*TMath::Sin(leptemp.Theta())) < 0.5) {
			  goodlep_n = goodlep_n + 1;
			  goodlep_index[lep_index] = i;
			  lep_index++;
			}
		      }
		    }
		}// tight
            }
	  
	  
	  //Exactly two good leptons, leading lepton with pT > 22 GeV and the subleading lepton with pT > 15 GeV
	  if(goodlep_n==2)
	    {
	      
	      int goodlep1_index = goodlep_index[0];
	      int goodlep2_index = goodlep_index[1];
	      
              if(lep_pt->at(goodlep1_index) > 22000)
		{
		  
		  //two different-flavour opposite-sign leptons
		  if ( lep_charge->at(goodlep1_index) * lep_charge->at(goodlep2_index)  < 0 ) 
		    {
		      if ( lep_type->at(goodlep1_index) != lep_type->at(goodlep2_index) )
			{
			  
			  // TLorentzVector definitions
			  TLorentzVector Lepton_1  = TLorentzVector();
			  TLorentzVector Lepton_2  = TLorentzVector();
			  TLorentzVector      MeT  = TLorentzVector();
			  
			  Lepton_1.SetPtEtaPhiE(lep_pt->at(goodlep1_index), lep_eta->at(goodlep1_index), lep_phi->at(goodlep1_index),lep_E->at(goodlep1_index));
			  Lepton_2.SetPtEtaPhiE(lep_pt->at(goodlep2_index), lep_eta->at(goodlep2_index), lep_phi->at(goodlep2_index),lep_E->at(goodlep2_index));
			  MeT.SetPtEtaPhiE(met_et, 0, met_phi , met_et);
			  
			  TLorentzVector     Lepton_12 = TLorentzVector();
			  Lepton_12 = Lepton_1 + Lepton_2;
			  
			  float mLL       = Lepton_12.Mag()/1000.;
			  float ptLL      = Lepton_12.Pt()/1000.;
			  
			  float dPhi_LL  = TMath::Abs(lep_phi->at(goodlep1_index) - lep_phi->at(goodlep2_index) );
			  dPhi_LL        = dPhi_LL < TMath::Pi() ? dPhi_LL : 2*TMath::Pi() - dPhi_LL;
			  
			  Float_t MET = met_et/1000.;
			  
			  float dPhiLLmet = TMath::Abs( Lepton_12.Phi() - MeT.Phi() );
			  dPhiLLmet    = dPhiLLmet < TMath::Pi() ? dPhiLLmet : 2*TMath::Pi() - dPhiLLmet;
			  
			  float mt    = sqrt(2*Lepton_12.Pt()*MeT.Et()*(1-cos(Lepton_12.DeltaPhi(MeT))))/1000.;
			  
			  
			  //Preselection of good jets
			  int goodjet_n = 0;
			  int goodbjet_n = 0;
			  
			  int goodjet_index[jet_n];
			  int jet_index = 0;
			  
			  int goodbjet_index[jet_n];
			  int bjet_index = 0;
			  
			  for(unsigned int i=0; i<jet_n; i++)
			    {
			      if(jet_pt->at(i) > 20000. && TMath::Abs(jet_eta->at(i)) < 2.5)
				{
				  // JVT cleaning
				  bool jvt_pass=true;
				  if (jet_pt->at(i) < 60000. && TMath::Abs(jet_eta->at(i)) < 2.4 && jet_jvt->at(i) < 0.59) jvt_pass=false;
				  if (jvt_pass)
				    {
				      
				      // cut on 85% WP
				      if ( jet_MV2c10->at(i) > 0.1758475  && TMath::Abs(jet_eta->at(i)) < 2.5 )
					{
					  goodbjet_n++;
					  goodbjet_index[bjet_index] = i;
					  bjet_index++;
					}
				      
				      if (jet_pt->at(i)>30000.)
					{
					  goodjet_n++;
					  goodjet_index[jet_index] = i;
					  jet_index++;
					}
				      
				    }
				}
			    }
			  
			  // fill histograms before any cut

			  double names_of_global_variable_pre[]={mLL, ptLL, dPhi_LL, dPhiLLmet, MET, mt, (double)goodjet_n, (double)goodbjet_n};
			  TString histonames_of_global_variable_pre[]={"histI_mLL", "histI_ptLL", "histI_dPhi_LL", "histI_dPhiLLmet", "histI_etmiss", "histI_mt", "histI_n_jets", "histI_n_bjets"};
			  int length_global_pre = sizeof(names_of_global_variable_pre)/sizeof(names_of_global_variable_pre[0]);
			  for (int i=0; i<length_global_pre; i++)                                              {
			    FillHistogramsGlobal( names_of_global_variable_pre[i], weight, histonames_of_global_variable_pre[i]);
			  }
			  
			  //  remove low mass meson resonances and DY events; ggF regions, at least 1 jet
			  if ( mLL > 10 && goodjet_n <= 1 && MET > 20) //mt (&& mt < 120) not here initially MET originally over 20 but now under 70
			    {
			      if ( dPhiLLmet > TMath::Pi()/2 )
				{
				  
				  if ( ptLL > 30 )
				    {
				      
				      if ( mLL < 55 )	//45, 55 before, trying to Find pure higgs
					{
					  
					  if ( dPhi_LL < 1.8 ) 
					    {      
					      
					      if ( goodbjet_n ==0 ) 
						{
						double xiYplus=TMath::Sin(lep_phi->at(goodlep2_index)); //whis is - and +
						double xiXplus=TMath::Cos(lep_phi->at(goodlep2_index));
						double xiXminus=TMath::Cos(lep_phi->at(goodlep1_index));
						double xiYminus=TMath::Sin(lep_phi->at(goodlep1_index));
						//xiYminvec.push_back(xiYminus);
						//xiXplusvec.push_back(xiXplus);
						//xiYplusvec.push_back(xiYminus);
						//xiXminvec.push_back(xiXminus);
						xpLL.push_back(dPhi_LL);
						xmLL.push_back(mLL);
						xpLLmet.push_back(dPhiLLmet);
						xpt.push_back(mt);
						xMET.push_back(MET);
						
						Iarrterm1.push_back(xiXplus*xiXminus+xiYplus*xiYminus);
						Iarrterm2.push_back((xiXplus*xiXplus-xiYplus*xiYplus)*(xiXminus*xiXminus-	xiYminus*xiYminus));
						Iarrterm3.push_back(xiXminus*xiYminus*xiXplus*xiYplus);
						  //Start to fill histograms : definitions of x-axis variables
						  
						  double names_of_global_variable[]={mLL, ptLL, dPhi_LL, dPhiLLmet, MET, mt, (double)goodjet_n, (double)goodbjet_n};
						  
						  double names_of_leadlep_variable[]={Lepton_1.Pt()/1000., Lepton_1.Eta(), Lepton_1.E()/1000., Lepton_1.Phi(), (double)lep_charge->at(goodlep1_index), (double)lep_type->at(goodlep1_index) };
						  
						  double names_of_subleadlep_variable[]={Lepton_2.Pt()/1000., Lepton_2.Eta(), Lepton_2.E()/1000., Lepton_2.Phi(), (double)lep_charge->at(goodlep2_index), (double)lep_type->at(goodlep2_index)};
						  
						  
						  //Start to fill histograms : definitions of histogram names
						  
						  TString histonames_of_global_variable[]={"hist_mLL", "hist_ptLL", "hist_dPhi_LL", "hist_dPhiLLmet", "hist_etmiss", "hist_mt", "hist_n_jets", "hist_n_bjets"};
						  TString histonames_of_leadlep_variable[]={"hist_leadleptpt", "hist_leadlepteta", "hist_leadleptE", "hist_leadleptphi", "hist_leadleptch", "hist_leadleptID"};
						  TString histonames_of_subleadlep_variable[]={"hist_subleadleptpt", "hist_subleadlepteta", "hist_subleadleptE", "hist_subleadleptphi", "hist_subleadleptch", "hist_subleadleptID"};
						  
						  int length_global = sizeof(names_of_global_variable)/sizeof(names_of_global_variable[0]);
						  int length_leadlep = sizeof(names_of_leadlep_variable)/sizeof(names_of_leadlep_variable[0]);
						  int length_subleadlep = sizeof(names_of_subleadlep_variable)/sizeof(names_of_subleadlep_variable[0]);
						  
						  //Fill histograms
						  
						  for (int i=0; i<length_global; i++)
						    {
						      FillHistogramsGlobal( names_of_global_variable[i], weight, histonames_of_global_variable[i]);
						    }
						  
						  for (int i=0; i<length_leadlep; i++)
						    {
						      FillHistogramsLeadlept( names_of_leadlep_variable[i], weight, histonames_of_leadlep_variable[i]);
						    }
						  
						  for (int i=0; i<length_subleadlep; i++)
						    {
						      FillHistogramsSubleadlept( names_of_subleadlep_variable[i], weight, histonames_of_subleadlep_variable[i]);
						    }
						  
						  
						}
					    }
					}
				    } // selection			      
				} // jet cut
			    }
			}
		    }
		}
	    }
	}
    }

  return kTRUE;
}


void HWWAnalysis::SlaveTerminate()
{
}

void scatterplot(vector< float >& x, vector< float >& y, string option,string xlabel, string ylabel)
{
  std::string Folder="Plots "+option;
  string head = xlabel+ylabel;
  char* chary = &ylabel[0];
  char* charx = &xlabel[0];
  char* charhead = &head[0];    
  char* chardescription = &head[0];
  TH2F *hist = new TH2F(charhead,chardescription, x.size()/4500+30, *min_element(x.begin(),x.end()), *max_element(x.begin(),x.end()),y.size()/4500+30, *min_element(y.begin(),y.end()), *max_element(y.begin(),y.end()));
 
  for (int i=0;i<(int)x.size();i++)
    {
    hist->Fill(x[i],y[i]);
    }
  TCanvas *c = new TCanvas();
  hist->GetXaxis()->SetTitle(charx);
  hist->GetYaxis()->SetTitle(chary);
  hist->Draw("colz");

  std::string name = Folder+"/"+head+".png";
  c->SaveAs(name.c_str());
}

void overlay(string name)
{
	string filenameA="Pure"+name+".root";
	string filenameB="simulation"+name+".root";
	char* char_nameA = &filenameA[0];
	char* char_nameB = &filenameB[0];
	char* char_name = &name[0];
			
	TFile fA(char_nameA);
	TH1F *hA= (TH1F*)fA.Get(char_name);
	TFile fB(char_nameB);
	TH1F *hB= (TH1F*)fB.Get(char_name);
	TCanvas *c = new TCanvas();
	
	hA->Scale(1/hA->GetEntries());
	hB->Scale(1/hB->GetEntries());
	
	hA->SetLineStyle(1);
	hB->SetLineStyle(1);
	
	hB->SetLineColor(kRed);
	hA->SetLineColor(kBlue);	
	
	hA->Draw();
	hB->Draw("same");
	
	auto legend = new TLegend();
	legend->SetHeader("Legend","C");
	legend->AddEntry(hA,"Pure Higgs signal","lep");
	legend->AddEntry(hB,"Background signal","lep");
	legend->Draw();
	
	hB->GetYaxis()->SetTitle("Normalized Count");	
	hA->GetYaxis()->SetTitle("Normalized Count");	
	
	c->SaveAs((name+".png").c_str());
}

void hist(vector< float >& x, float coeff, string option, string head, string description)
{
  std::string Folder="Histograms "+option;

  char* charhead = &head[0];
  char* chardescription = &description[0];
  
  
  char* char_name;
  string str_name=option+head+".root";
  char_name = &str_name[0];
  
  TFile *f = new TFile(char_name,"recreate");
  f->cd();
  TH1F *hist = new TH1F(charhead,chardescription,x.size()/3000+30, *min_element(x.begin(),x.end()), *max_element(x.begin(),x.end()) );
  for (int i=0;i<(int)x.size();i++)
    {
    hist->Fill(coeff*x[i]);
    }
  hist->GetXaxis()->SetTitle(charhead);
  hist->GetYaxis()->SetTitle("Count");
  hist->SetDirectory(0);
  TH1::AddDirectory(kFALSE);
  
  hist->Write();
  TCanvas *c = new TCanvas();
  hist->Draw();
  delete f;
  std::string name = Folder+"/"+head+".png";
  c->SaveAs(name.c_str());
}

float var(vector< float >& x,float X)
{
float sum=0;
  for (int i=0;i<(int)x.size();i++)
    {
    sum+=pow(static_cast<double>(x[i]) - X,2);
    }
  double var=1.0/static_cast<double>(x.size()*(x.size()-1))*sum;
return var;
}


float cov(vector< float >& x,float X,vector< float >& y,float Y)
{
float sum=0;
  for (int i=0;i<(int)x.size();i++)
    {
    sum+=((x[i]) - X)*((y[i]) - Y);
    }
  double cov=1.0/static_cast<double>(x.size()*(x.size()-1))*sum;
return cov;
}

void HWWAnalysis::Terminate()
{
  
  TString filename_option = GetOption();
  printf("Writting with name option: %s \n", filename_option.Data());
  TString output_name="Output_HWWAnalysis/"+filename_option+".root";
  const char* filename = output_name;
  
//  Iarrterm1={Iarrterm1.begin(),Iarrterm1.begin()+80};
//  Iarrterm2={Iarrterm2.begin(),Iarrterm2.begin()+80};
//  Iarrterm3={Iarrterm3.begin(),Iarrterm3.begin()+80};
  
  double Iaverageterm1 = std::accumulate(Iarrterm1.begin(), Iarrterm1.end(), 0.0)/Iarrterm1.size();
  double Iaverageterm2 = std::accumulate(Iarrterm2.begin(), Iarrterm2.end(), 0.0)/Iarrterm2.size();
  double Iaverageterm3 = std::accumulate(Iarrterm3.begin(), Iarrterm3.end(), 0.0)/Iarrterm3.size();
  
  double Ivarterm1=var(Iarrterm1,Iaverageterm1);
  double Ivarterm2=var(Iarrterm2,Iaverageterm2);
  double Ivarterm3=var(Iarrterm3,Iaverageterm3); 
  
  double cov12=cov(Iarrterm1,Iaverageterm1,Iarrterm2,Iaverageterm2);
  double cov23=cov(Iarrterm2,Iaverageterm2,Iarrterm3,Iaverageterm3);
  double cov13=cov(Iarrterm1,Iaverageterm1,Iarrterm3,Iaverageterm3);
      
  float coeff1=8/sqrt(3);
  float coeff2=25;
  float coeff3=100;
  
//  hist(Iarrterm1, coeff1, GetOption(), "Term 1", "Distribution Term 1");
//  hist(Iarrterm2, coeff2, GetOption(), "Term 2", "Distribution Term 2");  
//  hist(Iarrterm3, coeff3, GetOption(), "Term 3", "Distribution Term 3");
  
//  hist(xiYplusvec, 1, GetOption(), "xiYPlus", "Distribution xiYPlus");
//  hist(xiXplusvec, 1, GetOption(), "xiXplus", "Distribution xiXplusvec");
//  hist(xiYminvec, 1, GetOption(), "xiYmin","Distribution xiYminvec");
//  hist(xiXminvec, 1, GetOption(), "xiXmin", "Distribution xiXminvec");
  
//  hist(xpt, 1, GetOption(), "pt", "Distribution pt");
//  hist(xpLL, 1, GetOption(), "pLL", "Distribution phiLL");
//  hist(xmLL, 1, GetOption(), "mLL", "Distribution mLL");
//  hist(xMET, 1, GetOption(), "MET", "Distribution MET");
//  hist(xpLLmet, 1, GetOption(), "pLLmet", "Distribution pLLmet");
  
  scatterplot(xpLL,Iarrterm1, GetOption(),"term1 ","dphiLL");
  scatterplot(xmLL,Iarrterm1, GetOption(),"term1"," mLL");
  scatterplot(xpt,Iarrterm1, GetOption(),"term1"," mt");
  scatterplot(xMET,Iarrterm1, GetOption(),"term1"," xMet");
  scatterplot(xpLLmet,Iarrterm1, GetOption(),"term1"," pLLMET");
  
  scatterplot(xpLL,Iarrterm2,GetOption(),"term2 ","dphiLL");
  scatterplot(xmLL,Iarrterm2,GetOption(),"term2 ","mLL");
  scatterplot(xpt,Iarrterm2,GetOption(),"term2 ","mt");
  scatterplot(xMET,Iarrterm2, GetOption(),"term2 ","MET");
  scatterplot(xpLLmet,Iarrterm2, GetOption(),"term2 ","pLLMET");
  
  scatterplot(xpLL,Iarrterm3,GetOption(),"term3"," dphiLL");
  scatterplot(xmLL,Iarrterm3,GetOption(),"term3"," mLL");
  scatterplot(xpt,Iarrterm3,GetOption(),"term3"," mt");
  scatterplot(xMET,Iarrterm3,GetOption(),"term3"," MET");
  scatterplot(xpLLmet,Iarrterm3,GetOption(),"term3"," pLLMET");
  
  scatterplot(xmLL,xpLL, GetOption(),"mLL"," dphiLL");  
  scatterplot(xmLL,xpt, GetOption(),"mLL"," pt");  

//  overlay("Term 1");
//  overlay("Term 2");
//  overlay("Term 3");
      
//  overlay("xiYPlus");
//  overlay("xiXplus");
//  overlay("xiYmin");
//  overlay("xiXmin");
  
//  overlay("pt");
//  overlay("pLL");
//  overlay("mLL");
//  overlay("MET");
//  overlay("pLLmet");  
      
  double I=coeff1*Iaverageterm1+coeff2*Iaverageterm2+coeff3*Iaverageterm3;
  double Istd=sqrt(pow(coeff1,2)*Ivarterm1+pow(coeff2,2)*Ivarterm2+pow(coeff3,2)*Ivarterm3+2*coeff1*coeff2*cov12+2*coeff2*coeff3*cov23+2*coeff1*coeff3*cov13);
  
  cout << "cov 1-2 " << 2*coeff2*coeff1*cov12 << ", cov 2-3 " << 2*coeff2*coeff3*cov23 << ", cov 1-3 " << 2*coeff1*coeff3*cov13 <<"\n";  
  cout << "I1= " << coeff1*Iaverageterm1 << " +- " << coeff1*sqrt(Ivarterm1) << "\n";  
  cout << "I2= " << coeff2*Iaverageterm2 << " +- " << coeff2*sqrt(Ivarterm2) << "\n";  
  cout << "I3= " << coeff3*Iaverageterm3 << " +- " << coeff3*sqrt(Ivarterm3) << "\n";
  cout << "I= " << I << " +- " << Istd << ", " << TMath::Abs((I-2)/Istd) << " sigmas from 2, after " << Iarrterm3.size() << " events." << "\n";

  TFile physicsoutput(filename,"recreate");
  WriteHistograms();
  physicsoutput.Close();
}
