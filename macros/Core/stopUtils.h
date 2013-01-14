#ifndef stopUtils_h
#define stopUtils_h

#include <algorithm>
#include <iostream>
#include <map>
#include <vector>
#include <set>
#include <list>
#include <sstream>

#include "STOPT.h"

#include "Math/LorentzVector.h"
#include "Math/Point3D.h"
#include "TMath.h"
#include "TBranch.h"
#include "TTree.h"
#include "TH1F.h"
#include "TFile.h"
#include "TBits.h"

using namespace std;

float getdltrigweight(int id1, int id2);
float getsltrigweight(int id1, float pt, float eta);

bool passEvtSelection(TString name);
bool passOneLeptonSelection(bool isData);
bool passTwoLeptonSelection(bool isData);
bool passIsoTrkVeto();
bool passIsoTrkVeto_v2();
bool passSingleLeptonSelection(bool isData);
bool passDileptonSelection(bool isData);
bool passLepPlusIsoTrkSelection(bool isData);
pair<float,float> getPhiCorrMET( float met, float metphi, int nvtx, bool ismc);

float getDataMCRatio(float eta);

float vtxweight_n( const int nvertices, TH1F *hist, bool isData );
float getdphi( float phi1 , float phi2 );
float dRbetweenVectors(LorentzVector vec1,LorentzVector vec2 );
float getMinDphi(float metPhi, LorentzVector vec1, LorentzVector vec2 );


#endif
