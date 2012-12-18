#ifndef stopUtils_h
#define stopUtils_h

#include <algorithm>
#include <iostream>
#include <map>
#include <vector>
#include <set>
#include <list>
#include <sstream>

#include "StopTree.h"

#include "Math/LorentzVector.h"
#include "Math/Point3D.h"
#include "TMath.h"
#include "TBranch.h"
#include "TTree.h"
#include "TH1F.h"
#include "TFile.h"
#include "TBits.h"

using namespace std;

list<Candidate> MT2Calculator(StopTree* tree, bool isData);
MT2struct Best_MT2Calculator(StopTree* tree, bool isData);

float getdltrigweight(int id1, int id2);
float getsltrigweight(int id1, float pt, float eta);

bool passEvtSelection(const StopTree *sTree, TString name);
bool passOneLeptonSelection(const StopTree *sTree, bool isData);
bool passTwoLeptonSelection(const StopTree *sTree, bool isData);
bool passIsoTrkVeto(const StopTree *sTree);
bool passIsoTrkVeto_v2(const StopTree *sTree);
bool passSingleLeptonSelection(const StopTree *sTree, bool isData);
bool passDileptonSelection(const StopTree *sTree, bool isData);
bool passLepPlusIsoTrkSelection(const StopTree *sTree, bool isData);
pair<float,float> getPhiCorrMET( float met, float metphi, int nvtx, bool ismc);

float vtxweight_n( const int nvertices, TH1F *hist, bool isData );
float getdphi( float phi1 , float phi2 );
float dRbetweenVectors(StopTree::LorentzVector vec1,StopTree::LorentzVector vec2 );

#endif
