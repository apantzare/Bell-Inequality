// main01.cc is a part of the PYTHIA event generator.
// Copyright (C) 2023 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Keywords: basic usage; charged multiplicity;

//g++ -I/home/axel/pythia8309/include `root-config --cflags` main91.cc -o main91 -lpythia8 -L/home/axel/pythia8309/lib `root-config --glibs`
//export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/home/axel/pythia8309/lib
// This is a simple test program. It fits on one slide in a talk.
// It studies the charged multiplicity distribution at the LHC.

#include "Pythia8/Pythia.h"
#include <iostream> 
#include <stdlib.h>
#include "TFile.h"
#include "TTree.h"
#include "TThread.h"
#include <cmath>
using namespace std;

double get_angle(double x, double y) {
    double angle = std::atan2(y, x);
    if (angle > M_PI) {
        angle -= 2 * M_PI;
    }
    return angle;
}
    
double fix_angle(double angle) {
    angle = fmod(angle + M_PI, 2*M_PI); // map angle to [0, 2*pi)
    if (angle < 0)
        angle += 2*M_PI; // map angle to [0, 2*pi]
    return angle - M_PI; // map angle to [-pi, pi]
}
    
using namespace Pythia8;

void *handle(void *ptr){
int ith = (long)ptr;
TFile *output = new TFile(Form("New_Data_%d.root ",ith), "recreate");

TTree *tree = new TTree("mini", "mini");

int njet;
int runnumber=1;
int eventnummer;
int channelnummer=1;
int mcWeight=1;

int nlep = 2;
int njetpre = 10;
std::vector<float> eTjetsS;
std::vector<float> etaJetsS;
std::vector<float> phiJetsS;
std::vector<float> Ejet;

std::vector<int> charge;  
std::vector<float> Eta;
std::vector<float> phi; 
std::vector<float> pT; 
std::vector<float> Elep; 
std::vector<int> type;    		
double met_et;
double met_phi;        


tree->Branch("runNumber", &runnumber);  
tree->Branch("eventNumber", &eventnummer);  
tree->Branch("channelNumber", &channelnummer);  
tree->Branch("mcWeight", &mcWeight);  

tree->Branch("jet_n", &njet);  
tree->Branch("lep_n", &nlep);   
tree->Branch("jet_pt", &eTjetsS);    
tree->Branch("jet_eta", &etaJetsS);    
tree->Branch("jet_phi", &phiJetsS);    
tree->Branch("jet_E", &Ejet);    
//tree->Branch("m", &m, "m/D");    
tree->Branch("lep_charge", &charge);   
tree->Branch("lep_eta", &Eta);    
tree->Branch("lep_phi", &phi);    
tree->Branch("lep_pt", &pT);    
tree->Branch("lep_E", &Elep);    
tree->Branch("lep_type", &type);    
tree->Branch("met_et", &met_et);    
tree->Branch("met_phi", &met_phi);  

 // Set up generation. Incoming pp beams is default.
Pythia pythia; // Declare Pythia object

int nevents = 1e3;

pythia.readString("print:quiet = on");
pythia.readString("Beams:eCM = 13000."); // 13 TeV CM energy.

pythia.readString("HiggsSM:gg2H = on"); // Switch on g g -> H.
pythia.readString("25:onMode = off");
pythia.readString("25:onIfAll = 24 24");

pythia.readString("24:onMode = off");
pythia.readString("24:onIfMatch = 11 12");
pythia.readString("24:onIfMatch = -11 12");
pythia.readString("24:onIfMatch = 13 14");
pythia.readString("24:onIfMatch = -13 14");

pythia.readString("Random:setSeed = on");
pythia.readString(Form("Random:seed = %d", ith));


pythia.init(); // Initialize.

// Common parameters for the two jet finders.
double etaMax   = 4.;
double radius   = 0.7;
double pTjetMin = 10.;
// Exclude neutrinos (and other invisible) from study.
int    nSel     = 2;
// Range and granularity of CellJet jet finder.
int    nEta     = 80;
int    nPhi     = 64;

// Set up SlowJet jet finder, with anti-kT clustering
// and pion mass assumed for non-photons..
SlowJet slowJet( -1, radius, pTjetMin, etaMax, nSel, 1);

// Set up CellJet jet finder.
CellJet cellJet( etaMax, nEta, nPhi, nSel);

// Generate event(s).
int Wcount=0;
nlep=0;
int ids[2] = {0,0};
int idx[2] = {0,0};

for(int i=0; i < nevents; i++){
	if(!pythia.next()) continue;
	int entries=pythia.event.size();
	Wcount=0;
	nlep=0;
	int id;

	for(int j=0; j < entries; j++){
    	id=pythia.event[j].id();
        if(id==24){
      		Wcount++;
      	}
      	if (Wcount==2){
        	if(id==11 || id==-11 || id==13 || id==-13){
        		//std::cout << i << std::endl;
        		//std::cout << id << std::endl;
          		nlep++;
          	}
         } 
      }                
      
	if(Wcount==2 && nlep==2){
    	nlep=0;
		for(int j=0; j < entries; j++){
    		int id = pythia.event[j].id();
    		if(id==11 || id==-11 || id==13 || id==-13){
    			ids[nlep]=abs(id);
       			idx[nlep]=j;
       			nlep++;

   			}
   		}

		if(ids[0]!=ids[1]){
		//std::cout << i << std::endl;
    		slowJet.analyze( pythia.event );
    		njet = slowJet.sizeJet();
    		double ptx=0.;
    		double pty=0.;  

    		for (int k = 0; k < njet; ++k) {
              eTjetsS.push_back(slowJet.pT(k)*1000) ;
              etaJetsS.push_back(slowJet.y(k)) ;
              phiJetsS.push_back(get_angle(cos(slowJet.phi(k)),sin(slowJet.phi(k))));              
              Ejet.push_back(slowJet.p(k)[0]*1000);
              ptx +=slowJet.p(k)[1];
              pty +=slowJet.p(k)[2];
            }
    		
            eventnummer=i;
    		//std::cout << "Event" << i << std::endl;

            for (int k = 0; k < nlep; ++k) {
                //m[k] = pythia.event[idx[k]].m();
        		charge.push_back(pythia.event[idx[k]].chargeType()/3);  
        		Eta.push_back(pythia.event[idx[k]].eta());
                phi.push_back(pythia.event[idx[k]].phi()); 
                pT.push_back(pythia.event[idx[k]].pT()*1000); 
                Elep.push_back(pythia.event[idx[k]].e()*1000); 
                type.push_back(ids[k]);
                ptx +=pythia.event[idx[k]].pT()*cos(pythia.event[idx[k]].phi());
                pty +=pythia.event[idx[k]].pT()*sin(pythia.event[idx[k]].phi());
            }
            met_et = 1000*sqrt(pow(ptx,2)+pow(pty,2));
            met_phi = get_angle(ptx,pty);
            tree->Fill();
    		
            eTjetsS.clear();
            etaJetsS.clear();
            phiJetsS.clear();
            Ejet.clear();

            charge.clear();
            Eta.clear();
            phi.clear();
            pT.clear();
            Elep.clear();
            type.clear();
    	}
    }
}

output->Write();
output->Close();
pythia.next(); // Generate an(other) event. Fill event record.

}

int main() {
    int ncores=128;
    TThread *th[ncores];
    for(int i=0; i<ncores; i++){
        th[i]=new TThread(Form("th%d",i), handle, (void*)i);
        th[i]->Run();
    }
    for(int i=0; i<ncores; i++){
        th[i]->Join();
    }
  return 0;
}
