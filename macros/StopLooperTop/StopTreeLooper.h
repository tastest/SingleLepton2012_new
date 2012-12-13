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
const float PTMIN_J1 = 30;
const float PTMIN_J2 = 30;
const float PTMIN_B = 30;
const float PTMIN_O = 30;
const float JET_PT = 30.;
const float JET_ETA = 3.0;
const float PDG_TOP_MASS = 173.5;
const float PDG_W_MASS = 80.385;
const float BTAG_MED = 0.679;
const float BTAG_LOW = 0.244;

const int N_JETS_TO_CONSIDER = 6;
const int NUM_LEAD_JETS_0B = 3;
const int NUM_LEAD_JETS_1B = 3;

const bool __SORT = true;

class StopTree;

        typedef struct {
                float one_chi2;
                float two_mt2b, two_mt2bl, two_mt2w;
                float three_mt2b, three_mt2bl, three_mt2w;
                float four_chi2b, four_chi2bl, four_chi2w;
        } MT2CHI2;


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
	float vtxweight_n(const int nvertices, TH1F *hist, bool isData);
	float dRbetweenVectors(LorentzVector vec1, LorentzVector vec2);
	float getdphi(float phi1, float phi2);
	float getMT(float pt1, float phi1, float pt2, float phi2);
	pair<float, float> getPhiCorrMET(float met, float metphi, int nvtx,
			bool ismc);

	list<Candidate> recoHadronicTop(StopTree*, bool);
        list<Candidate> applyBConsistency(list<Candidate> &candidates, StopTree* tree, float btagcut);
        MT2CHI2 MT2CHI2Calculator(list<Candidate> candidates, StopTree* tree);


private:

	string m_outfilename_;
	//for phi corrected met
	float t1metphicorr;
	float t1metphicorrphi;
	float t1metphicorrmt;
	//for mt peak definition
	float min_mtpeak;
	float max_mtpeak;

	int n_jets;
	vector<LorentzVector> jets;
	vector<float> btag;
	vector<int> mc;

};

#endif
