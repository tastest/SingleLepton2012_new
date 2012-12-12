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

	// weights
	Float_t weight_;
	Float_t sltrigeff_;
	Float_t dltrigeff_;

	// hadronic variables
	Int_t nb_;
	Int_t njets_;

	// lepton variables
	Int_t passisotrk_;
	Int_t nlep_;

	Float_t lep1pt_;
	Float_t lep1eta_;

	Float_t lep2pt_;
	Float_t lep2eta_;
	Float_t dilmass_;

    private:

	string m_outfilename_;
	//for phi corrected met
	/* float t1metphicorr; */
	/* float t1metphicorrphi; */
	/* float t1metphicorrmt; */
	/* //for mt peak definition */
	/* float min_mtpeak; */
	/* float max_mtpeak;  */

};

#endif
