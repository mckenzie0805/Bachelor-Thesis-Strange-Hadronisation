// This script will provide trees with the pT phi eta PDG codes so they can be used to create more correlatiion plots. 
// In principle every vector is an event and every element of the vector is a particle produced at that event.  It also 
// creates a B+B+ correlation histogram so we can compare it with histogram created by the trees to check if everything works
// alright. Finally it produces a multiplicity histogram and Strange particles per event histogram.

//C++ libraries we are using
#include <iostream>
#include <cmath>
#include <cstring>
#include <chrono>
#include <vector>
//Root and pythia libraries
#include "Pythia8/Pythia.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TH1I.h"
#include "TTree.h"

#define PI 3.14159265
using namespace std;
using namespace Pythia8;

// Here we define some functions we are going to use
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

Double_t DeltaPhi(Double_t phi1, Double_t phi2){
	// Returns delta phi in range -pi/2 3pi/2
	return fmod(phi1-phi2+2.5*PI,2*PI)-0.5*PI;
	}


int main(int argc, char** argv){

	if(argc != 2){ // Main sould have two arguments one will be the ./CharmDeltaPhi and one other will be the filename or path	
		cout<<"Error in the number of arguments provided"<<endl;
		cout<<"Provide only the filepath/name."<<endl;
		cout<<"Terminating program"<<endl;
		return 1;
	}

	// Here we start keeping track of time
	auto start = chrono::high_resolution_clock::now();
	
	// Create output file
	TFile* output = new TFile(argv[1],"CREATE");
	if(!output->IsOpen()){
		cout<<"Error: File "<<argv[1]<<"already exists terminating program!"<<endl;
		return 1;
	}
	
	// Here I define the tree that the data will be stored.
	TTree *tree = new TTree("tree","bbbar correlations");

	
	// Here we define the variables we are going to need.
	Double_t pT, eta, phi;
	//Double_t pTTriger, pTAssociate;
	Double_t charge;
	Double_t status,mother,motherID;
	Int_t  id,size,nEvents,Strangeness;
	
	// Here I define the vectors
	vector<Int_t> vID;
	vector<Double_t> vPt;
	vector<Double_t> vEta;
	vector<Double_t> vPhi;
	vector<Double_t> vCharge;
	vector<Double_t> vStatus;
	vector<Double_t> vMother1;
	vector<Double_t> vMotherID;

	
	// Setting up tree Branches and histograms
	tree->Branch("ID",&vID);
	tree->Branch("PT",&vPt);
	tree->Branch("ETA",&vEta);
	tree->Branch("PHI",&vPhi);
	tree->Branch("CHARGE",&vCharge);
	tree->Branch("STATUS",&vStatus);
	tree->Branch("MOTHER",&vMother1);
	tree->Branch("MOTHERID",&vMotherID);
	tree->Branch("MULTIPLICITY",&size,"x/I");
	TH1D* hSize = new TH1D("hSize","Multiplicity",301,-0.5,300.5);
	TH1D* hidStrange = new TH1D("hidStrange","PDG Codes for Strange hadrons",12000,-6000,6000);
	//TH1D* hPtTriger = new TH1D("hPtTriger","p_{T} for triger B^{+} ",50,0,10);
	//TH1D* hPtAssociate = new TH1D("hPtAssociate", "p_{T} for associate B^{+}",50,0,10);
	TH1D* hStrangePart = new TH1D("hStrangePart", "Strange Particles Per Event",200,-0.5,200.5);

	
	// Kinematics constraints
	const Double_t pTmin = 0.15;//minimum pT
	const Double_t etamax = 4.;// maximum eta
	
	// The simulation is an object so we define it like
	Pythia pythia;
	
	// Simulation settings from a file
	//pythia.readFile("ssbar_monash.cmnd");
	pythia.readFile("ssbar_junctions.cmnd");
	//pythia.readFile("ssbar_ee.cmnd");
	//pythia.readFile("main07.cmnd");
	nEvents = pythia.mode("Main:numberOfEvents");
	
	// Here we create a radnom seed using the runtime so eachtime the simulation runs the outcome will be trully random.
	Int_t proccessid = getpid();
	string seedstr = "Random:seed = "+std::to_string((time(0)+proccessid)%900000000);
	pythia.readString("Random:setSeed = on");
	pythia.readString(seedstr);
	
	// Initializing simulation.
	pythia.init(); // Important must be at every simulation
	
	cout<<"Generating "<<nEvents<<" events!"<<endl;
	
	// Here we start the event loop
	for(int iEvent = 0; iEvent<nEvents; iEvent++){
		if(!pythia.next()) continue; // Generate the next event if there is an error in the current. This command must exist at every pythia script.
	
		int nPart = pythia.event.size(); // Particles produced in this event
		size = 0; // Initialiazing for multiplicity plot.
		Strangeness = 0; // Intialiazing for Strange production plot.
		// Initializing vectors
		vID.clear();
		vPt.clear();
		vEta.clear();
		vPhi.clear();
		vCharge.clear();
		vStatus.clear();
		vMother1.clear();
		vMotherID.clear();
		// Particle loop
		for(int iPart = 0; iPart<nPart; iPart++){
			const Particle &particle = pythia.event[iPart];
			if(!particle.isFinal()) continue; // Skip if the particle is not at its final state.
			
			id = particle.id();
			if((TMath::Abs(id) != 211)&&(TMath::Abs(id) != 2212)&&(!IsStrange(id))) continue;
			//if(!IsStrange(id)) continue; // check for strangeness, otherwise we don't care

			pT = particle.pT();
			eta = particle.eta();
			phi = particle.phi();
			charge = particle.charge();
			status = static_cast<Double_t> (particle.status());
			mother = static_cast<Double_t> (particle.mother1());
			motherID = static_cast<Double_t> (pythia.event[particle.mother1()].id());
			
			// Kinematics check
			if(pT < pTmin || abs(eta) > etamax ) continue;
			
			if(charge != 0) size++;
			
			//idStrange = id;
			hidStrange->Fill((Double_t) id);
			Strangeness++;
			// Filling vectors
			vID.push_back(id);
			vPt.push_back(pT);
			vEta.push_back(eta);
			vPhi.push_back(phi);
			vCharge.push_back(charge);
			vStatus.push_back(status);
			vMother1.push_back(mother);
			vMotherID.push_back(motherID);

		} // 1st particle loop
		hSize->Fill((Double_t) size);
		hStrangePart->Fill((Double_t) Strangeness);
		// In order not fill trees with empty vectors.
		// if(vID.empty() || vPt.empty() || vEta.empty() || vPhi.empty() || vCharge.empty() || vStatus.empty() || vMother1.empty() || vMotherID.empty() ) continue; 
		tree->Fill();
	} // End of event loop
	
	// write output and close it
	
	output->Write();
	cout<<"File has been created and its name is: "<<output->GetName()<<endl;
	output->Close();
	
	// Stop keeping time
	auto end = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::minutes>(end - start);
        cout << "This script took " << duration.count() << " minutes to run." << endl;
	
	return 0;
}




