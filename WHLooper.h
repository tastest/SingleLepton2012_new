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
  void fillHists2D(std::map<std::string, TH2F*>& h_2d, const float evtweight = 1., const std::string& dir = "", const std::string& suffix = "");
  void fillFlavorHists1D(std::map<std::string, TH1F*>& h_1d, const float evtweight = 1., const std::string& dir = "", const std::string& suffix = "");
  void fillJetAccHists(std::map<std::string, TH1F*>& h_1d, const float evtweight = 1., const std::string& dir = "", const std::string& suffix = "");

  float getCSVCut(const csvpoint csv = WHLooper::CSVM);

  void dumpEventInfo(const std::string& comment);

  string m_outfilename_;
  TFile* outfile_;
  //for phi corrected met
  float t1metphicorr;
  float t1metphicorrphi;
  float t1metphicorrmt;

  // store main cut values
  float CUT_BBMASS_LOW_;
  float CUT_BBMASS_HIGH_;
  float CUT_BBMASS_CR1_LOW_;
  float CUT_BBMASS_CR1_HIGH_;
  float CUT_BBMASS_CR8_LOW_;
  float CUT_BBMASS_CR8_HIGH_;
  float CUT_MET_PRESEL_;
  float CUT_MET_;
  float CUT_MT_PRESEL_;
  float CUT_MT_;
  float CUT_MT_CR13_LOW_;
  float CUT_MT_CR13_HIGH_;
  float CUT_MT2BL_;

  // internal vars to store
  std::vector<LorentzVector> jets_;
  std::vector<LorentzVector> bjets_;
  std::vector<LorentzVector> jets_fwd_;
  std::vector<float> jets_csv_;
  std::vector<int> jets_idx_;
  std::vector<int> jets_fwd_idx_;
  std::vector<int> bjets_idx_;
  std::vector<float> jets_smearcorrs_;
  int njets_;
  int njetsalleta_;
  int nbjets_;
  int nbjetst_;
  int nbjetsl_;
  float met_;
  float metphi_;
  float mt_;
  float mt2b_;
  float mt2bl_;
  float mt2w_;
  LorentzVector bb_;
  float met_soft_;
  float sumet_;
  float sumet_soft_;
  float ht_;
  float wpt_;
  float lepmetdphi_;
  float bbwdphi_;
  // for ttbar samples
  float genmt2bl_;
  bool tobtecveto_;

  // for CR3 (with two leptons)
  LorentzVector lep_;
  float pseudomet_lep_;
  float pseudometphi_lep_;
  float pseudomt_lep_;
  float dphi_pseudomet_lep_;
  float pseudomt2b_;
  float pseudomt2bl_;
  float pseudomt2w_;

  // internal flags
  bool isWjets_;
  bool isWNjets_;
  bool isWNjets_nobb_;
  bool isWNjets_onlybb_;
  bool isWHbb_;
  bool isWino_;
  bool isTChiwh_;
  bool isTChihhwwbb_;
  bool isScan_;
  bool isttmg_;
  bool isttsl_;
  bool isttdl_;
  bool istsl_;
  bool istdl_;
  bool isttvsl_;
  bool isttvdl_;
  bool isttvother_;

};

#endif
