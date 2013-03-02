#ifndef STOPTREELOOPER_H
#define STOPTREELOOPER_H

#include "TChain.h"
#include "TFile.h"
#include "TString.h"
#include "TH1F.h"

#include <iostream>
//#include "Math/LorentzVector.h"
/////#include "../Core/StopTree.h"
#include "../Core/STOPT.h"
 
#include <cmath>
#include <map>
#include <list>

#if not defined(__CINT__) || defined(__MAKECINT__)
// needs to be included when makecint runs (ACLIC)
#include "TMVA/Factory.h"
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#endif


using namespace std;

/* ------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------- */
const double PTMIN_J1   = 25;
const double PTMIN_J2   = 25;
const double PTMIN_BTAG = 30;
const double PTMIN_OTAG = 30; 
const double PTMIN_B    = 30;  //  This Two should be tigther than the  
const double PTMIN_O    = 30;  //  b-tagged versions.
const double BTAG_MIN = 0.679;

const bool __SORT = true;



class StopTree;

class StopTreeLooper {

    public:
   typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > LorentzVector;

        StopTreeLooper();
        ~StopTreeLooper();

        void setOutFileName(string filename); 
        void loop(TChain *chain, TString name);

	//	MT2struct Best_MT2Calculator_Ricardo(list<Candidate>, StopTree*, bool);
	//        list<Candidate> recoHadronicTop(StopTree*, bool, bool);
	//        list<Candidate> getBTaggedCands(list<Candidate> &candidates, StopTree* tree);
	void initBaby();
	void makeTree(const char*, TChain *chain);

        TTree  *outTree_;
        TFile  *outFile_;

	// which selections are passed
	Int_t sig_;
	Int_t cr1_;   
	Int_t cr4_;
	Int_t cr5_;   

	// cut and count selections
	Float_t t2ttLM_;
	Float_t t2ttHM_;

	// kinematic variables
	Float_t met_;
	Float_t chi2_;
	Float_t mt2w_;

	// event shapes
	Float_t htssl_;
	Float_t htosl_;
	Float_t htratiol_;
	Float_t htssm_;
	Float_t htosm_;
	Float_t htratiom_;

	Float_t dphimj1_;
	Float_t dphimj2_;
	Float_t dphimjmin_;

	// weights
	Float_t weight_;
	Float_t sltrigeff_;
	Float_t dltrigeff_;

	// jet counting
	Float_t nb_;
	Int_t njets_;

	// lepton variables
	Int_t passisotrk_;
	Int_t nlep_;

	Float_t lep1pt_;
	Float_t lep1eta_;
	Float_t dRleptB1_;

	Float_t lep2pt_;
	Float_t lep2eta_;
	Float_t dilmass_;

	// jet kinematics
	Float_t pt_b_;
	Float_t pt_J1_;
	Float_t pt_J2_;

	// susy variables
	Float_t mstop_;
	Float_t mlsp_;
	Float_t x_;

	Float_t rand_;
	vector<float> bdt_;

    private:

	static const float JET_PT = 30.;
	static const float JET_ETA = 2.4;

	static const int N_JETS_TO_CONSIDER = 6;
	static const int NJETS_CUT = 4;

	string m_outfilename_;

	//jet information
	vector<LorentzVector> jets;
	vector<float> jets_btag;
	vector<float> jets_sigma;
	float metphi;

        static const bool __apply_mva = false; 

};

#endif
