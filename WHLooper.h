#ifndef WHLOOPER_H
#define WHLOOPER_H

#include "TChain.h"
#include "TFile.h"
#include "TString.h"
#include "TH1F.h"
#include "TH2F.h"

#include <iostream>
#include "Math/LorentzVector.h"
 
#include <cmath>
#include <map>

using namespace std;

class WHLooper {

 public:
  typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > LorentzVector;

  enum csvpoint { CSVL, CSVM, CSVT };

  WHLooper();
  ~WHLooper();

  void setOutFileName(string filename); 
  void loop(TChain *chain, TString name);

 private:

  void fillHists1DWrapper(std::map<std::string, TH1F*>& h_1d, const float evtweight = 1., const std::string& dir = "");
  void fillHists1D(std::map<std::string, TH1F*>& h_1d, const float evtweight = 1., const std::string& dir = "", const std::string& suffix = "");
  std::vector<LorentzVector> getBJets(const csvpoint csv = WHLooper::CSVM);
  float getCSVCut(const csvpoint csv = WHLooper::CSVM);

  void MakeBabyNtuple (const char *);
  void FillBabyNtuple (const float evtweight);
  void CloseBabyNtuple ();

  string m_outfilename_;
  TFile* outfile_;
  //for phi corrected met
  float t1metphicorr;
  float t1metphicorrphi;
  float t1metphicorrmt;
  //for mt peak definition
  float min_mtpeak;
  float max_mtpeak; 
  // for storing mt2 values
  //  float mt2w_;

  std::vector<LorentzVector> myJets_;
  std::vector<LorentzVector> myBJets_;
  std::vector<int> myBJetsIdx_;

  // internal flags
  bool isWjets_;
  bool isTChiwh_;

  // for output minibaby
  TFile *babyFile_;
  TTree *babyTree_;

  UInt_t  run_;
  UInt_t  lumi_;
  UInt_t  event_;
  Int_t   leptype_;
  Float_t weight_;

  Float_t pfmet_;
  Float_t lep1mt_;
  Float_t mt2b_;
  Float_t mt2bl_;
  Float_t mt2w_;
  Float_t bbpt_;
  Float_t wpt_;
  Float_t bbwdphi_;

};

#endif
