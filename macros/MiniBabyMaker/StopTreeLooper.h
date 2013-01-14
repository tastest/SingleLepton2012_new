#ifndef STOPTREELOOPER_H
#define STOPTREELOOPER_H

#include "TChain.h"
#include "TFile.h"
#include "TString.h"
#include "TH1F.h"

#include <iostream>
//#include "Math/LorentzVector.h"
#include "../Core/StopTree.h"
 
#include <cmath>
#include <map>
#include <list>

using namespace std;

/* ------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------- */
const double PTMIN_J1   = 25;
const double PTMIN_J2   = 25;
const double PTMIN_BTAG = 30;
const double PTMIN_OTAG = 30; 
const double PTMIN_B    = 30;  //  This Two should be tigther than the  
const double PTMIN_O    = 30;  //  b-tagged versions.
const double PDG_TOP_MASS = 173.5;
const double PDG_W_MASS = 80.385;
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

	MT2struct Best_MT2Calculator_Ricardo(list<Candidate>, StopTree*, bool);
        list<Candidate> recoHadronicTop(StopTree*, bool, bool);
        list<Candidate> getBTaggedCands(list<Candidate> &candidates, StopTree* tree);
	void initBaby();
	void makeTree(const char*);

        TTree  *outTree_;
        TFile  *outFile_;

	// which selections are passed
	Int_t sig_;
	Int_t cr1_;   
	Int_t cr4_;
	Int_t cr5_;   
	
	// kinematic variables
	Float_t met_;
	Float_t mt_;
	Float_t chi2_;
	Float_t mt2w_;
	Float_t mt2b_;
	Float_t mt2bl_;

	// "best" chi2 and MT2 variables
	Float_t chi2min_;
	Float_t chi2minprob_;
	Float_t chi2min_mt2b_;  
	Float_t chi2min_mt2bl_; 
	Float_t chi2min_mt2w_;  
	Float_t mt2bmin_;       
	Float_t mt2blmin_;      
	Float_t mt2wmin_;       
	Float_t mt2bmin_chi2_;  
	Float_t mt2blmin_chi2_; 
	Float_t mt2wmin_chi2_;  
	Float_t mt2bmin_chi2prob_;  
	Float_t mt2blmin_chi2prob_; 
	Float_t mt2wmin_chi2prob_;  
	Int_t   ncand_;

	// event shapes
	Float_t thrjet_;
	Float_t apljet_;
	Float_t sphjet_;
	Float_t cirjet_;

	Float_t thrjetl_;
	Float_t apljetl_;
	Float_t sphjetl_;
	Float_t cirjetl_;

	Float_t thrjetlm_;
	Float_t apljetlm_;
	Float_t sphjetlm_;
	Float_t cirjetlm_;

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

	// hadronic variables
	Int_t nb_;
	Int_t njets_;

	// lepton variables
	Int_t passisotrk_;
	Int_t passisotrkv2_;
	Int_t nlep_;

	Float_t lep1pt_;
	Float_t lep1eta_;

	Float_t lep2pt_;
	Float_t lep2eta_;
	Float_t dilmass_;

	// susy variables
	Float_t mstop_;
	Float_t mlsp_;
	Float_t x_;

	Float_t rand_;

    private:

	static const float JET_PT = 30.;
	static const float JET_ETA = 3.0;

	static const int N_JETS_TO_CONSIDER = 6;

	string m_outfilename_;

	int n_jets;
	vector<LorentzVector> jets;
	vector<float> btag;
	vector<int> mc;

};

#endif
