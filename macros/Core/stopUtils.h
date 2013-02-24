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

class Candidate : public TObject {
 public:
  float chi2, mt2w, mt2bl, mt2b;
  int j1, j2, bi, oi;
  float k1, k2;
  bool match;
  
  ClassDef(Candidate, 2)
    };

typedef vector<Candidate> CANDIDATES;


int leadingJetIndex(vector<LorentzVector> jets, int iskip1, int iskip2);
unsigned int getNJets();
vector<int> getBJetIndex(double discr, int iskip1, int iskip2);

float getdltrigweight(int id1, int id2);
float getsltrigweight(int id1, float pt, float eta);

bool passEvtSelection(TString name);
bool passOneLeptonSelection(bool isData);
bool passTwoLeptonSelection(bool isData);
bool passIsoTrkVeto();
bool passIsoTrkVeto_v2();
bool passIsoTrkVeto_v3();
bool passIsoTrkVeto_v4();
bool passSingleLeptonSelection(bool isData);
bool passDileptonSelection(bool isData);
bool passLepPlusIsoTrkSelection(bool isData);
bool passLepPlusIsoTrkSelection_noEMu(bool isData, bool isMu);

pair<float,float> getPhiCorrMET( float met, float metphi, int nvtx, bool ismc);

float getDataMCRatio(float eta);
float getDataMCRatioFix(float eta);

float vtxweight_n( const int nvertices, TH1F *hist, bool isData );
float getdphi( float phi1 , float phi2 );
float dRbetweenVectors(LorentzVector vec1,LorentzVector vec2 );
float getMinDphi(float metPhi, LorentzVector vec1, LorentzVector vec2 );
float getMT( float pt1 , float phi1 , float pt2 , float phi2 );

struct DorkyEventIdentifier {
  // this is a workaround for not having unique event id's in MC
  unsigned long int run, event,lumi;
  bool operator < (const DorkyEventIdentifier &) const;
  bool operator == (const DorkyEventIdentifier &) const;
};

bool is_duplicate (const DorkyEventIdentifier &id, std::set<DorkyEventIdentifier> &already_seen);
int load_badlaserevents  (char* filename, std::set<DorkyEventIdentifier> &events_lasercalib);
bool is_badLaserEvent (const DorkyEventIdentifier &id, std::set<DorkyEventIdentifier> &events_lasercalib);
double calculateMT2w(vector<LorentzVector> jets, vector<float> btag, LorentzVector lep, float met, float metphi);
double mt2wWrapper(LorentzVector lep, LorentzVector jet_o, LorentzVector jet_b, float met, float metphi);
double calculateChi2(vector<LorentzVector> jets, vector<float> sigma_jets);
double fc2 (double c1, double m12, double m22, double m02, bool verbose);
double fchi2 (double c1, double pt1, double sigma1, double pt2, double sigma2, double m12, double m22, double m02);
void minuitFunction(int&, double* , double &result, double par[], int);
double getChi2(LorentzVector jets_b, LorentzVector jets_j1, LorentzVector jets_j2, float sigma_b, float sigma_j1, float sigma_j2);

static const float BTAG_MED = 0.679;
static const float PDG_TOP_MASS = 173.5;
static const float PDG_W_MASS = 80.385;



#endif
