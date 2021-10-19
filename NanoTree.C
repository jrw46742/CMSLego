#define NanoTree_cxx
#include "NanoTree.h"
#include <stdlib.h>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <algorithm>
#include <cmath>
#define RESET   "\033[0m"
#define BLACK   "\033[30m"      /* Black */
#define RED     "\033[31m"      /* Red */
#define GREEN   "\033[32m"      /* Green */
#define YELLOW  "\033[33m"      /* Yellow */
#define BLUE    "\033[34m"      /* Blue */
#define MAGENTA "\033[35m"      /* Magenta */
#define CYAN    "\033[36m"      /* Cyan */
#define WHITE   "\033[37m"      /* White */
#define BOLDBLACK   "\033[1m\033[30m"      /* Bold Black */
#define BOLDRED     "\033[1m\033[31m"      /* Bold Red */
#define BOLDGREEN   "\033[1m\033[32m"      /* Bold Green */
#define BOLDYELLOW  "\033[1m\033[33m"      /* Bold Yellow */
#define BOLDBLUE    "\033[1m\033[34m"      /* Bold Blue */
#define BOLDMAGENTA "\033[1m\033[35m"      /* Bold Magenta */
#define BOLDCYAN    "\033[1m\033[36m"      /* Bold Cyan */
#define BOLDWHITE   "\033[1m\033[37m"      /* Bold White */

#ifdef _WIN32
#include <Windows.h>
#else
#include <unistd.h>
#endif
#include <iostream>
#include <cstdlib>
#include <stdlib.h>     
using namespace std;

void NanoTree::Loop()
{
//   In a ROOT session, you can do:
//      root> .L NanoTree.C
//      root> NanoTree t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
/////////////////////////////////////////////////////////////////////////////////
//Smarter Nano:
//store muon_jetidx[0] and electron_jetidx[0]
//loop jets and don't print if ijet == those stored indices;   also don't print if fabs(eta) > 1.5
//maybe count N GOOD jets with pT > 30 && not outside our range && not the muon or electorn and print a "NOPE, MOVE ON" message if < 2   
   
   if (fChain == 0) return;    
   Long64_t nentries = fChain->GetEntriesFast();
   Long64_t nbytes = 0, nb = 0;

   int skip = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;

      if(fabs(Jet_eta[0]) < .7 && nJet > 2 && nMuon >= 2 && nElectron >= 1) {
	   cout << BOLDBLUE << "Event: " << event << RESET << endl;
	   for(int ijet = 0; ijet < nJet; ijet++){ //jet loop
	      if (fabs(Jet_eta[ijet]) < .7 && ijet != Muon_jetIdx[0] && ijet != Electron_jetIdx[0] && nJet > 2 && Jet_pt[ijet] > 30){
   	      std::cout << "Jet " << ijet << ": pt = " << Jet_pt[ijet] << ", eta = " << Jet_eta[ijet] << ", phi = " << Jet_phi[ijet] << std::endl;
	      }	
	   }	 
	   for(int i = 0; i < nMuon; i++){ //muon loop
	      cout << "Muon " << i << ": pt = " << Muon_pt[i] << ", eta = " << Muon_eta[i] << ", phi = " << Muon_phi[i] << ", Muon_jetIdx: " << Muon_jetIdx[i] << endl;
	   }
	   for(int i = 0; i < nElectron; i++){ //electron loop
    	      cout << "Electron " << i << ": pt = " << Electron_pt[i] << ", eta = " << Electron_eta[i] << ", phi = " << Electron_phi[i] << ", Electron_jetIdx: " << Electron_jetIdx[i] << endl;
	   }
	   for(int i = 0; i < nPhoton; i++){ //photon loop
	      cout << "Photon " << i << ": pt = " << Photon_pt[i] << ", eta = " << Photon_eta[i] << ", phi = " << Photon_phi[i] <<  endl;
	   }
	   cout << "\n" << endl;
      }
      else { 
	      cout << "########### SKIPPED EVENT: "<< event <<" #############\n"  << endl;
	      skip++;
	   }  	
   }
   cout << "Skipped " << skip << " events" << endl;
}

