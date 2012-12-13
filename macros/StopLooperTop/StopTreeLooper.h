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
