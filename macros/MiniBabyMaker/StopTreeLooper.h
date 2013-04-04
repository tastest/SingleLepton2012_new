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
        void doWHNtuple(); 
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

	Int_t whsig_;
	Int_t whcr1_;   

	// cut and count selections
	Float_t t2ttLM_;
	Float_t t2ttHM_;

	// kinematic variables
	Float_t mt_;
	Float_t mtup_;
	Float_t mtdown_;
	Float_t met_;
	Float_t metup_;
	Float_t metdown_;
	Float_t chi2_;
	Float_t chi2up_;
	Float_t chi2down_;
	Float_t mt2b_;
	Float_t mt2bl_;
	Float_t mt2w_;
	Float_t mt2wup_;
	Float_t mt2wdown_;

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
	Float_t nvtxweight_;
	Float_t sltrigeff_;
	Float_t dltrigeff_;
	Int_t   nsigevents_;

	// jet counting
	Int_t nb_;
	Int_t njets_;
	Int_t njets_up_;
	Int_t njets_down_;

	// lepton variables
	Int_t passisotrk_;
	Int_t passtauveto_;
	Int_t nlep_;

	Float_t lep1pt_;
	Float_t lep1eta_;
	Float_t dRleptB1_;

	Float_t lep2pt_;
	Float_t lep2eta_;
	Float_t dilmass_;

	// jet kinematics
	Float_t pt_b_;
	Float_t pt_b_up_;
	Float_t pt_b_down_;
	Float_t pt_J1_;
	Float_t pt_J2_;

        // wh event kinematics
	Float_t bbmass_;
	Float_t bbpt_;
	Float_t wpt_;
	Float_t bbwdphi_;

	// susy variables
	Float_t mstop_;
	Float_t mlsp_;
	Float_t x_;
	Float_t xsecsusy_;

	Float_t rand_;
	vector<float> bdt_;

    private:

	float JET_PT;
	float JET_ETA;
	float BJET_ETA;

	int N_JETS_TO_CONSIDER;
	int NJETS_CUT;

    float MET_CUT;

	string m_outfilename_;
    string m_minibabylabel_;

	//jet information
	vector<LorentzVector> jets;
	vector<LorentzVector> jets_up;
	vector<LorentzVector> jets_down;
	vector<LorentzVector> bjets;
	vector<float> jets_btag;
	vector<float> jets_up_btag;
	vector<float> jets_down_btag;
	vector<float> jets_sigma;
	vector<float> jets_up_sigma;
	vector<float> jets_down_sigma;
	float metphi;
	float metupphi;
	float metdownphi;

    static const bool __apply_mva = true; 
    static const bool __mini_branches = true;
    static const bool __add_babies = true; 

};

#endif
