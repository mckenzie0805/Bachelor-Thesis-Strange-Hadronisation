// Authors: Ioannis Sidiras, Paul Veen, Rik Spijkers

// C++ libraries
#include <iostream>
#include <vector>
//ROOT libraries
#include "TFile.h"
#include "TH1.h"
#include "TH1D.h"
#include "TTree.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TString.h"
#include "THStack.h"
#include "TAttMarker.h"
#include "TLegend.h"

#define PI 3.14159265
using namespace std;
using namespace TMath;

Double_t DeltaPhi(Double_t phi1, Double_t phi2){
	// Returns delta phi in range -pi/2 3pi/2
	return fmod(phi1-phi2+2.5*PI,2*PI)-0.5*PI;
}
	
bool IsStrange(Int_t particlepdg) {
	// Checks if given pdg code is that of a (anti)strange hadron, returns True/False
	// This also find ssbar mesons such as the phi
	Int_t pdg = std::abs(particlepdg);
	pdg /= 10; // get rid of the last digit, is not important for quark content
	if (pdg % 10 == 3) return true; // 3rd quark
	pdg /= 10;
	if (pdg % 10 == 3) return true; // 2nd quark
	pdg /= 10;
	if (pdg % 10 == 3) return true; // 1st quark
	return false;
}
	
void ssbar_CorrelationsAnalysis(TString inputFiles, TString outfilename){
	
	// Define the TChain
	TChain *ch1 = new TChain("tree");
	ch1->Add(inputFiles);

	TFile *output = new TFile(outfilename,"RECREATE");
	
	// Now we define vectors that carry the information at event level.
	vector<Int_t>* vID = 0;
	vector<Double_t>* vPt = 0;
	vector<Double_t>* vPhi = 0;
	vector<Double_t>* vStatus = 0;
	vector<Double_t>* vEta = 0;
	vector<Double_t>* vMother1 = 0;
	vector<Double_t>* vMotherID = 0;
	// Setting up chain branch addresses to the vectors defined above
	ch1->SetBranchAddress("ID",&vID);
	ch1->SetBranchAddress("PT",&vPt);
	ch1->SetBranchAddress("PHI",&vPhi);
	ch1->SetBranchAddress("ETA",&vEta);
	ch1->SetBranchAddress("STATUS",&vStatus);
	ch1->SetBranchAddress("MOTHER",&vMother1);
	ch1->SetBranchAddress("MOTHERID",&vMotherID);
	
	// Definition of variables
	Int_t aID, pID;
	Double_t pPt, pPhi, pStatus, pEta, pMotherID; //For trigger
	Double_t aPt, aPhi, aStatus, aEta, aMotherID; //For associate
	int nTrigger = 0;
	
	// Each vector is an event number of events analyzed is the total number of vectors
	int nEvents = ch1->GetEntries();
	
	cout<<"The number of events for this analysis is: "<<nEvents<<endl;
	
	// Definition of produced histograms
	// 3D histograms
	TH3D *hPttPtaDPhi = new TH3D("hDPtrPaDEta", "#Delta#eta p_{T} Trigger p_{T} Associate; p_{T} trigger (Gev/c); p_{T} associate (GeV/c); #Delta#eta", 100, 0, 50, 100, 0, 50, 80, -8, 8);
	
	// 2D histograms
	TH2D *hDPhiDEta = new TH2D("hDPhiDEta", "#Delta#phi and #Delta#eta; #Delta#phi (rad); #Delta#eta", 100, -PI/2, 3*PI/2, 80, -8, 8);
	
	// 1D histograms
	TH1D *hTrPt = new TH1D("hTrPt", "Trigger Transverse Momentum; p_{T} GeV/c; Counts", 100, 0, 50);
	TH1D *hDPhi = new TH1D("hDPhi", "#Delta#phi; #Delta#Phi (rad); Counts", 100, -PI/2, 3*PI/2);
	
	// Event Loop
	for(int iEvent = 0; iEvent < nEvents; iEvent++){
		ch1->GetEntry(iEvent);
		int nparticles = vID->size();
		for(int ipart = 0; ipart < nparticles; ipart++){
			pID = (*vID)[ipart];
			pPhi = (*vPhi)[ipart];
			pPt = (*vPt)[ipart];
			pStatus = (*vStatus)[ipart];
			pEta =(*vEta)[ipart];
			pMotherID = (*vMotherID)[ipart];
				
			if(abs(pID) == 3312){ // 3312 = pdg code for xi^-
				nTrigger++;
				hTrPt->Fill(pPt);
				
				for(int jpart = 0; jpart < nparticles; jpart++){
					if(jpart == ipart) continue; // Do not correlate with itself
					aID = (*vID)[jpart];
					aPhi = (*vPhi)[jpart];
					aPt = (*vPt)[jpart];
					aStatus = (*vStatus)[jpart];
					aEta = (*vEta)[jpart];
					aMotherID = (*vMotherID)[jpart];

					if(abs(aID) == 3312){// 3312 = pdg code for xi^-
						// if particles have opposite signed pdg codes then they count as +1 to our histogram (signal)
						// if they have the same sign we want to subtract the value from the histograms (background)
						int sign = 1;
						if (pID*aID > 0) sign = -1;
						hPttPtaDPhi->Fill(pPt, aPt, sign);
						hDPhiDEta->Fill(DeltaPhi(pPhi, aPhi), pEta-aEta, sign);
						hDPhi->Fill(DeltaPhi(pPhi, aPhi), sign);
					} // Associate Condition	
				} // Associate Loop
			} // Trigger Condition
		} // Trigger Loop
	} // Event loop

	output->Write();
	output->Close();
	cout << "The total number of triggers is: " << nTrigger << endl;
	
	cout << "File: " << outfilename << " has been created!" << endl;
}

