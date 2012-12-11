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

	//selection
	bool passEvtSelection(const StopTree *sTree, TString name);
	bool passOneLeptonSelection(const StopTree *sTree, bool isData);
	bool passTwoLeptonSelection(const StopTree *sTree, bool isData);
	bool passSingleLeptonSelection(const StopTree *sTree, bool isData);
	bool passDileptonSelection(const StopTree *sTree, bool isData);
	bool passLepPlusIsoTrkSelection(const StopTree *sTree, bool isData);
	bool passIsoTrkVeto(const StopTree *sTree);

	//helper
	float getdltrigweight(int id1, int id2);
	float getsltrigweight(int id1, float pt, float eta);
	float vtxweight_n( const int nvertices, TH1F *hist, bool isData );
	float dRbetweenVectors(LorentzVector vec1, LorentzVector vec2 );
	float getdphi( float phi1 , float phi2 );
	float getMT( float pt1 , float phi1 , float pt2 , float phi2 );
	pair<float,float> getPhiCorrMET( float met, float metphi, int nvtx, bool ismc);

	MT2struct Best_MT2Calculator_Ricardo(list<Candidate>, StopTree*, bool);
        list<Candidate> recoHadronicTop(StopTree*, bool, bool);
        list<Candidate> getBTaggedCands(list<Candidate> &candidates, StopTree* tree);

    private:

	string m_outfilename_;
	//for phi corrected met
	float t1metphicorr;
	float t1metphicorrphi;
	float t1metphicorrmt;
	//for mt peak definition
	float min_mtpeak;
	float max_mtpeak; 

};

#endif
