#include "WHLooper.h"

//#include "../../CORE/jetSmearingTools.h"
//#include "../../CORE/Thrust.h"
//#include "../../CORE/EventShape.h"

#include "../Core/STOPT.h"
#include "../Core/stopUtils.h"
#include "../Plotting/PlotUtilities.h"
#include "../../Tools/BTagReshaping/BTagReshaping.h"
// #include "../Core/MT2Utility.h"
// #include "../Core/mt2bl_bisect.h"
// #include "../Core/mt2w_bisect.h"
// #include "../Core/PartonCombinatorics.h"

#include "TROOT.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TChainElement.h"
#include "TFile.h"
#include "TMath.h"
#include "TChain.h"
#include "TBenchmark.h"
#include "TVector2.h"
#include "Math/VectorUtil.h"
#include "Riostream.h"

#include <algorithm>
#include <utility>
#include <map>
#include <set>
#include <list>

using namespace Stop;

const bool doTauVeto = true;
const bool doLep2Veto = true;
const bool doCSVReshaping = false;
const bool doBtagSFs = true; // should only do one of reshaping and SFs at a time
const bool doWJetsOverlapRemoval = false; // obsolete, done automatically
const bool doISRReweight = false;
const bool doTopPtWeight = true;
const bool doJetSmearing = false;
const bool doWbbMtReweight = true;
const bool doTopPtReweight2 = false;

const bool doFlavorPlots = true;
const bool doNM1Plots = true;
const bool doWJetsPlots = false;
const bool doNvtxSplit = false;
const bool doJetAccPlots = false;

// regions to do
const bool blindSignal = true;
const bool doInclusive = false;
const bool doInclusiveMTTail = false;
const bool doSignal = false;
const bool doSignalMassLast = true;
const bool doSignalMETLast = true;
const bool doSignalSMWH = false;
const bool doCR1 = false; // high m(bb)
const bool doCR1METLast = true; // high m(bb), met cut last
const bool doCR2 = false; // lep + track
const bool doCR3 = false; // dilep
const bool doCR23 = true; // dilep + (lep+track)
const bool doCR4 = false; // dilep, high m(bb)
const bool doCR5 = false; // bveto
const bool doCR5METLast = true; // bveto
const bool doCR5InvMass = true; // bveto, inverted m(bb)
const bool doCR5MassLast = false; // bveto, mass last
const bool doCR6 = false; // 1 btag -- not used
const bool doCR6METLast = false; // 1 btag -- not used
const bool doCR7 = false; // high m(bb), 3 jets
const bool doCR8 = false; // low m(bb)
const bool doCR8METLast = true; // low m(bb)
const bool doCR9 = false; // high m(bb), 150-200
const bool doCR10 = false; // high m(bb), 200-250
const bool doCR11 = false; // bveto, high m(bb)
const bool doCR12 = false; // high m(bb), 4+ jets
const bool doCR13 = false; // mt bulk
const bool doCR14 = true; // inverted mbb region
const bool doStopSel = false;

std::set<DorkyEventIdentifier> already_seen; 
std::set<DorkyEventIdentifier> events_lasercalib; 
std::set<DorkyEventIdentifier> events_hcallasercalib; 

//--------------------------------------------------------------------
// 
// // This is meant to be passed as the third argument, the predicate, of the standard library sort algorithm
// inline bool sortByPt(const LorentzVector &vec1, const LorentzVector &vec2 ) {
//     return vec1.pt() > vec2.pt();
// 
// }
//--------------------------------------------------------------------

// from the internets: to replace a subset of a string
bool replace(std::string& str, const std::string& from, const std::string& to) {
    size_t start_pos = str.find(from);
    if(start_pos == std::string::npos)
        return false;
    str.replace(start_pos, from.length(), to);
    return true;
}

//--------------------------------------------------------------------

WHLooper::WHLooper()
{
  m_outfilename_ = "histos.root";
  t1metphicorr = -9999.;
  t1metphicorrphi = -9999.;
  t1metphicorrmt = -9999.;

  CUT_BBMASS_LOW_ = 100.0;
  //  CUT_BBMASS_HIGH_ = 140.0;
  CUT_BBMASS_HIGH_ = 150.0;
  CUT_BBMASS_CR1_LOW_ = 150.0;
  CUT_BBMASS_CR1_HIGH_ = 250.0;
  CUT_BBMASS_CR8_LOW_ = 50.0;
  CUT_BBMASS_CR8_HIGH_ = 100.0;
  CUT_MET_PRESEL_ = 50.0;
  // CUT_MET_PRESEL_ = 100.0;
  CUT_MET_ = 175.0;
  //  CUT_MT_PRESEL_ = 50.0;
  CUT_MT_PRESEL_ = 0.0;
  CUT_MT_ = 100.0;
  CUT_MT_CR13_LOW_ = 50.0;
  CUT_MT_CR13_HIGH_ = 100.0;
  CUT_MT2BL_ = 200.0;

}

//--------------------------------------------------------------------

WHLooper::~WHLooper()
{
}

//--------------------------------------------------------------------

void WHLooper::setOutFileName(string filename)
{
  m_outfilename_ = filename;

}

//--------------------------------------------------------------------

void WHLooper::loop(TChain *chain, TString name) {

  // Benchmark
  TBenchmark *bmark = new TBenchmark();
  bmark->Start("benchmark");

  if (name.Contains("wjets") || name.Contains("wbb")) {
    isWjets_ = true;
  } else {
    isWjets_ = false;
  }

  if (name.Contains("wjets")) {
    isWNjets_ = true;
  } else {
    isWNjets_ = false;
  }

  if (name.Contains("wjets_nobb")) {
    isWNjets_nobb_ = true;
    std::cout << "Recognize wjets_nobb sample" << std::endl;
  } else {
    isWNjets_nobb_ = false;
  }

  if (name.Contains("wjets_onlybb")) {
    isWNjets_onlybb_ = true;
    std::cout << "Recognize wjets_onlybb sample" << std::endl;
  } else {
    isWNjets_onlybb_ = false;
  }

  if (name.Contains("TChiwh")) {
    isTChiwh_ = true;
  } else {
    isTChiwh_ = false;
  }

  if (name.Contains("ttbar_mg")) {
    isttmg_ = true;
    std::cout << "Recognize ttbar madgraph sample" << std::endl;
  } else {
    isttmg_ = false;
  }

  if (name.Contains("ttbar_") && name.Contains("_1l")) {
    isttsl_ = true;
  } else {
    isttsl_ = false;
  }

  if (name.Contains("ttbar_") && name.Contains("_2l")) {
    isttdl_ = true;
  } else {
    isttdl_ = false;
  }

  if (name.Contains("single_top_1l")) {
    istsl_ = true;
  } else {
    istsl_ = false;
  }

  if (name.Contains("single_top_2l")) {
    istdl_ = true;
  } else {
    istdl_ = false;
  }

  if (name.Contains("ttv_1l")) {
    isttvsl_ = true;
    std::cout << "Recognize ttv_1l sample" << std::endl;
  } else {
    isttvsl_ = false;
  }

  if (name.Contains("ttv_2l")) {
    isttvdl_ = true;
    std::cout << "Recognize ttv_2l sample" << std::endl;
  } else {
    isttvdl_ = false;
  }

  if (name.Contains("ttv_other")) {
    isttvother_ = true;
    std::cout << "Recognize ttv_other sample" << std::endl;
  } else {
    isttvother_ = false;
  }

  //------------------------------
  // check for valid chain
  //------------------------------

  cout << "[WHLooper::loop] sample: " << name << endl;

  load_badlaserevents("../Core/badlaser_events.txt", events_lasercalib);
  load_badlaserevents("../Core/badhcallaser_events.txt", events_hcallasercalib);

  TObjArray *listOfFiles = chain->GetListOfFiles();
  TIter fileIter(listOfFiles);
  if (listOfFiles->GetEntries() == 0) {
    cout << "[WHLooper::loop] no files in chain" << endl;
    return;
  }

  BTagShapeInterface * nominalShape = new BTagShapeInterface("../../Tools/BTagReshaping/csvdiscr.root", 0.0, 0.0);

  //------------------------------
  // set up histograms
  //------------------------------


  gROOT->cd();

  cout << "[WHLooper::loop] creating output file: " << m_outfilename_ << endl;

  outfile_ = new TFile(m_outfilename_.c_str(),"RECREATE") ; 

  cout << "[WHLooper::loop] setting up histos" << endl;

  // inclusive region hists
  std::map<std::string, TH1F*> h_1d_inc_presel, h_1d_inc_2j, h_1d_inc_1b, h_1d_inc_2b, h_1d_inc_2j_mt, h_1d_inc_2j_mt_met100, h_1d_inc_2j_mt_met150, h_1d_inc_2j_mt_metcut;

  // signal region hists
  std::map<std::string, TH1F*> h_1d_sig_presel, h_1d_sig_final;
  // signal region nm1 hists
  std::map<std::string, TH1F*> h_1d_sig_met_nm1, h_1d_sig_mt_nm1, h_1d_sig_mt2bl_nm1;
  // signal mtpeak hists
  std::map<std::string, TH1F*> h_1d_sig_mtpeak_nomet, h_1d_sig_mtpeak_met;

  // signal region hists
  std::map<std::string, TH1F*> h_1d_sig_bbmasslast_presel, h_1d_sig_bbmasslast_final;
  // signal region nm1 hists
  std::map<std::string, TH1F*> h_1d_sig_bbmasslast_met_nm1, h_1d_sig_bbmasslast_mt_nm1, h_1d_sig_bbmasslast_mt2bl_nm1,h_1d_sig_bbmasslast_bbmass_nm1;
  // signal region nm1 hists
  std::map<std::string, TH1F*> h_1d_sig_bbmasslast_mtpeak, h_1d_sig_bbmasslast_mtcut, h_1d_sig_bbmasslast_mtcut_met100, h_1d_sig_bbmasslast_metcut, h_1d_sig_bbmasslast_met100, h_1d_sig_bbmasslast_met150;

  // signal region hists
  std::map<std::string, TH1F*> h_1d_sig_metlast_presel, h_1d_sig_metlast_final;
  // signal region nm1 hists
  std::map<std::string, TH1F*> h_1d_sig_metlast_met_nm1, h_1d_sig_metlast_mtfirst, h_1d_sig_metlast_mt_nm1, h_1d_sig_metlast_mt2bl_nm1, h_1d_sig_metlast_met100, h_1d_sig_metlast_met150;
  // 2d hists for correlations
  std::map<std::string, TH2F*> h_2d_sig_metlast_mt2bl_nm1;

  // signal SMWH region hists
  std::map<std::string, TH1F*> h_1d_sigsmwh_presel, h_1d_sigsmwh_final;
  // // signal region nm1 hists
  // std::map<std::string, TH1F*> h_1d_sig_met_nm1, h_1d_sig_mt_nm1, h_1d_sig_mt2bl_nm1;
  // // signal mtpeak hists
  // std::map<std::string, TH1F*> h_1d_sig_mtpeak_nomet, h_1d_sig_mtpeak_met;

  // cr1 hists
  std::map<std::string, TH1F*> h_1d_cr1_presel, h_1d_cr1_final;
  // cr1 nm1 hists
  std::map<std::string, TH1F*> h_1d_cr1_met_nm1, h_1d_cr1_mt_nm1, h_1d_cr1_mt2bl_nm1;
  // cr1 mtpeak hists
  std::map<std::string, TH1F*> h_1d_cr1_mtpeak_nomet, h_1d_cr1_mtpeak_met, h_1d_cr1_mtlast_mt2bl_nm1;

  // cr1, met last hists
  std::map<std::string, TH1F*> h_1d_cr1_metlast_presel, h_1d_cr1_metlast_final;
  // cr1, met last nm1 hists
  std::map<std::string, TH1F*> h_1d_cr1_metlast_mt_nm1, h_1d_cr1_metlast_mtfirst, h_1d_cr1_metlast_mt2bl_nm1, h_1d_cr1_metlast_met_nm1, h_1d_cr1_metlast_met100, h_1d_cr1_metlast_met150;

  // cr2 hists
  std::map<std::string, TH1F*> h_1d_cr2_presel, h_1d_cr2_final;
  // cr2 nm1 hists
  std::map<std::string, TH1F*> h_1d_cr2_bbmass_nm1, h_1d_cr2_met_nm1, h_1d_cr2_mt_nm1, h_1d_cr2_mt2blfirst, h_1d_cr2_mt2blcut, h_1d_cr2_mt2bl_nm1;

  // cr3 hists
  std::map<std::string, TH1F*> h_1d_cr3_presel, h_1d_cr3_final;
  // cr3 nm1 hists
  std::map<std::string, TH1F*> h_1d_cr3_bbmass_nm1, h_1d_cr3_met_nm1, h_1d_cr3_mt_nm1, h_1d_cr3_mt2blfirst, h_1d_cr3_mt2blcut, h_1d_cr3_mt2bl_nm1;

  // cr23 hists
  std::map<std::string, TH1F*> h_1d_cr23_presel, h_1d_cr23_final;
  // cr23 nm1 hists
  std::map<std::string, TH1F*> h_1d_cr23_bbmass_nm1, h_1d_cr23_met_nm1, h_1d_cr23_mt_nm1, h_1d_cr23_mt2blfirst, h_1d_cr23_mt2blcut, h_1d_cr23_mt2bl_nm1;

  // cr4 hists
  std::map<std::string, TH1F*> h_1d_cr4_presel, h_1d_cr4_final;
  // cr4 nm1 hists
  std::map<std::string, TH1F*> h_1d_cr4_met_nm1, h_1d_cr4_mt_nm1, h_1d_cr4_mt2bl_nm1;

  // cr5 hists
  std::map<std::string, TH1F*> h_1d_cr5_presel, h_1d_cr5_final;
  // cr5 nm1 hists
  std::map<std::string, TH1F*> h_1d_cr5_bbmass_nm1, h_1d_cr5_met_nm1, h_1d_cr5_met100, h_1d_cr5_mt_nm1, h_1d_cr5_mt2bl_nm1;

  // cr5, met last hists
  std::map<std::string, TH1F*> h_1d_cr5_metlast_presel, h_1d_cr5_metlast_final;
  // cr5, met last nm1 hists
  std::map<std::string, TH1F*> h_1d_cr5_metlast_bbmass_nm1, h_1d_cr5_metlast_mtfirst, h_1d_cr5_metlast_mt_nm1, h_1d_cr5_metlast_mt2bl_nm1, h_1d_cr5_metlast_met_nm1, h_1d_cr5_metlast_met100, h_1d_cr5_metlast_met150, h_1d_cr5_metlast_nomt_met100, h_1d_cr5_metlast_nomt_met150, h_1d_cr5_metlast_nomt_met175;

  // cr5, inv mass hists
  std::map<std::string, TH1F*> h_1d_cr5_invmass_presel, h_1d_cr5_invmass_final;
  // cr5, inv mass nm1 hists
  std::map<std::string, TH1F*> h_1d_cr5_invmass_bbmass_nm1, h_1d_cr5_invmass_mtfirst, h_1d_cr5_invmass_mt_nm1, h_1d_cr5_invmass_mt2bl_nm1, h_1d_cr5_invmass_met_nm1, h_1d_cr5_invmass_met100, h_1d_cr5_invmass_met150, h_1d_cr5_invmass_nomt_met100, h_1d_cr5_invmass_nomt_met150, h_1d_cr5_invmass_nomt_met175;

  // cr5 mass last region hists
  std::map<std::string, TH1F*> h_1d_cr5_bbmasslast_presel, h_1d_cr5_bbmasslast_final;
  // cr5 mass last region nm1 hists
  std::map<std::string, TH1F*> h_1d_cr5_bbmasslast_met_nm1, h_1d_cr5_bbmasslast_mt_nm1, h_1d_cr5_bbmasslast_mt2bl_nm1,h_1d_cr5_bbmasslast_bbmass_nm1;
  // cr5 mass last region nm1 hists
  std::map<std::string, TH1F*> h_1d_cr5_bbmasslast_mtpeak, h_1d_cr5_bbmasslast_mtcut, h_1d_cr5_bbmasslast_mtcut_met100, h_1d_cr5_bbmasslast_metcut, h_1d_cr5_bbmasslast_met100, h_1d_cr5_bbmasslast_met150;

  // cr6 hists
  std::map<std::string, TH1F*> h_1d_cr6_presel, h_1d_cr6_final;
  // cr6 nm1 hists
  std::map<std::string, TH1F*> h_1d_cr6_met_nm1, h_1d_cr6_mt_nm1, h_1d_cr6_mt2bl_nm1;

  // cr6, met last hists
  std::map<std::string, TH1F*> h_1d_cr6_metlast_presel, h_1d_cr6_metlast_final;
  // cr6, met last nm1 hists
  std::map<std::string, TH1F*> h_1d_cr6_metlast_bbmass_nm1, h_1d_cr6_metlast_mt_nm1, h_1d_cr6_metlast_mt2bl_nm1, h_1d_cr6_metlast_met_nm1, h_1d_cr6_metlast_met100, h_1d_cr6_metlast_met150;

  // cr7 hists
  std::map<std::string, TH1F*> h_1d_cr7_presel, h_1d_cr7_final;
  // cr7 nm1 hists
  std::map<std::string, TH1F*> h_1d_cr7_met_nm1, h_1d_cr7_mt_nm1, h_1d_cr7_mt2bl_nm1;

  // cr8 hists
  std::map<std::string, TH1F*> h_1d_cr8_presel, h_1d_cr8_final;
  // cr8 nm1 hists
  std::map<std::string, TH1F*> h_1d_cr8_met_nm1, h_1d_cr8_mt_nm1, h_1d_cr8_mt2bl_nm1;

  // cr8, met last hists
  std::map<std::string, TH1F*> h_1d_cr8_metlast_presel, h_1d_cr8_metlast_final;
  // cr8, met last nm1 hists
  std::map<std::string, TH1F*> h_1d_cr8_metlast_mt_nm1, h_1d_cr8_metlast_mtfirst, h_1d_cr8_metlast_mt2bl_nm1, h_1d_cr8_metlast_met_nm1, h_1d_cr8_metlast_met100, h_1d_cr8_metlast_met150;

  // cr9 hists
  std::map<std::string, TH1F*> h_1d_cr9_presel, h_1d_cr9_final;
  // cr9 nm1 hists
  std::map<std::string, TH1F*> h_1d_cr9_met_nm1, h_1d_cr9_mt_nm1, h_1d_cr9_mt2bl_nm1;
  // cr9 mtpeak hists
  std::map<std::string, TH1F*> h_1d_cr9_mtpeak_nomet, h_1d_cr9_mtpeak_met, h_1d_cr9_mtlast_mt2bl_nm1;

  // cr10 hists
  std::map<std::string, TH1F*> h_1d_cr10_presel, h_1d_cr10_final;
  // cr10 nm1 hists
  std::map<std::string, TH1F*> h_1d_cr10_met_nm1, h_1d_cr10_mt_nm1, h_1d_cr10_mt2bl_nm1;

  // cr11 hists
  std::map<std::string, TH1F*> h_1d_cr11_presel, h_1d_cr11_final;
  // cr11 nm1 hists
  std::map<std::string, TH1F*> h_1d_cr11_met_nm1, h_1d_cr11_mt_nm1, h_1d_cr11_mt2bl_nm1;

  // cr12 hists
  std::map<std::string, TH1F*> h_1d_cr12_presel, h_1d_cr12_final;
  // cr12 nm1 hists
  std::map<std::string, TH1F*> h_1d_cr12_bbmass_nm1, h_1d_cr12_met_nm1, h_1d_cr12_mt_nm1, h_1d_cr12_mt2bl_nm1;

  // cr13 hists
  std::map<std::string, TH1F*> h_1d_cr13_presel, h_1d_cr13_final;
  // cr13 nm1 hists
  std::map<std::string, TH1F*> h_1d_cr13_bbmass_nm1, h_1d_cr13_mt2bl_nm1, h_1d_cr13_met_nm1, h_1d_cr13_met100, h_1d_cr13_met150;

  // cr14 hists
  std::map<std::string, TH1F*> h_1d_cr14_presel, h_1d_cr14_final;
  // cr14 nm1 hists
  std::map<std::string, TH1F*> h_1d_cr14_mt_nm1, h_1d_cr14_mtfirst, h_1d_cr14_mt2bl_nm1, h_1d_cr14_met_nm1, h_1d_cr14_met100, h_1d_cr14_met150, h_1d_cr14_nomt_met100, h_1d_cr14_nomt_met150, h_1d_cr14_nomt_met175;

  // stop region hists
  std::map<std::string, TH1F*> h_1d_stop_presel, h_1d_stop_comp;
  std::map<std::string, TH1F*> h_1d_stop_met_nm1, h_1d_stop_mt_nm1, h_1d_stop_isotrk_nm1, h_1d_stop_tauveto_nm1 ;

  if (doInclusive) {
    outfile_->mkdir("inc_presel");
    outfile_->mkdir("inc_2j");
    if (doInclusiveMTTail) {
      outfile_->mkdir("inc_2j_mt");
      outfile_->mkdir("inc_2j_mt_met100");
      outfile_->mkdir("inc_2j_mt_met150");
      outfile_->mkdir("inc_2j_mt_metcut");
    }
    outfile_->mkdir("inc_1b");
    outfile_->mkdir("inc_2b");
  }

  if (doSignal) {
    outfile_->mkdir("sig_presel");
    if (doNM1Plots) {
      outfile_->mkdir("sig_met_nm1");
      outfile_->mkdir("sig_mtpeak_nomet");
      outfile_->mkdir("sig_mt_nm1");
      outfile_->mkdir("sig_mtpeak_met");
      outfile_->mkdir("sig_mt2bl_nm1");
    }
    outfile_->mkdir("sig_final");
  }

  if (doSignalMassLast) {
    outfile_->mkdir("sig_bbmasslast_presel");
    if (doNM1Plots) {
      outfile_->mkdir("sig_bbmasslast_mt2bl_nm1");
      outfile_->mkdir("sig_bbmasslast_mtpeak");
      outfile_->mkdir("sig_bbmasslast_mtcut");
      outfile_->mkdir("sig_bbmasslast_mtcut_met100");
      outfile_->mkdir("sig_bbmasslast_metcut");
      outfile_->mkdir("sig_bbmasslast_mt_nm1");
      outfile_->mkdir("sig_bbmasslast_met_nm1");
      outfile_->mkdir("sig_bbmasslast_met100");
      outfile_->mkdir("sig_bbmasslast_met150");
      outfile_->mkdir("sig_bbmasslast_bbmass_nm1");
    }
    outfile_->mkdir("sig_bbmasslast_final");
  }

  if (doSignalMETLast) {
    outfile_->mkdir("sig_metlast_presel");
    if (doNM1Plots) {
      outfile_->mkdir("sig_metlast_mt2bl_nm1");
      outfile_->mkdir("sig_metlast_mtfirst");
      outfile_->mkdir("sig_metlast_mt_nm1");
      outfile_->mkdir("sig_metlast_met_nm1");
      outfile_->mkdir("sig_metlast_met100");
      outfile_->mkdir("sig_metlast_met150");
    }
    outfile_->mkdir("sig_metlast_final");
  }

  if (doSignalSMWH) {
    outfile_->mkdir("sigsmwh_presel");
    // if (doNM1Plots) {
    //   outfile_->mkdir("sig_met_nm1");
    //   outfile_->mkdir("sig_mtpeak_nomet");
    //   outfile_->mkdir("sig_mt_nm1");
    //   outfile_->mkdir("sig_mtpeak_met");
    //   outfile_->mkdir("sig_mt2bl_nm1");
    // }
    outfile_->mkdir("sigsmwh_final");
  }

  if (doCR1) {
    outfile_->mkdir("cr1_presel");
    if (doNM1Plots) {
      outfile_->mkdir("cr1_met_nm1");
      outfile_->mkdir("cr1_mtpeak_nomet");
      outfile_->mkdir("cr1_mt_nm1");
      outfile_->mkdir("cr1_mtpeak_met");
      outfile_->mkdir("cr1_mt2bl_nm1");
      outfile_->mkdir("cr1_mtlast_mt2bl_nm1");
    }
    outfile_->mkdir("cr1_final");
  }

  if (doCR1METLast) {
    outfile_->mkdir("cr1_metlast_presel");
    if (doNM1Plots) {
      outfile_->mkdir("cr1_metlast_mt2bl_nm1");
      outfile_->mkdir("cr1_metlast_mtfirst");
      outfile_->mkdir("cr1_metlast_mt_nm1");
      outfile_->mkdir("cr1_metlast_met_nm1");
      outfile_->mkdir("cr1_metlast_met100");
      outfile_->mkdir("cr1_metlast_met150");
    }
    outfile_->mkdir("cr1_metlast_final");
  }

  if (doCR2) {
    outfile_->mkdir("cr2_presel");
    if (doNM1Plots) {
      outfile_->mkdir("cr2_bbmass_nm1");
      outfile_->mkdir("cr2_mt_nm1");
      outfile_->mkdir("cr2_mt2blfirst");
      outfile_->mkdir("cr2_met_nm1");
      outfile_->mkdir("cr2_mt2blcut");
      outfile_->mkdir("cr2_mt2bl_nm1");
    }
    outfile_->mkdir("cr2_final");
  }

  if (doCR3) {
    outfile_->mkdir("cr3_presel");
    if (doNM1Plots) {
      outfile_->mkdir("cr3_bbmass_nm1");
      outfile_->mkdir("cr3_mt_nm1");
      outfile_->mkdir("cr3_mt2blfirst");
      outfile_->mkdir("cr3_met_nm1");
      outfile_->mkdir("cr3_mt2blcut");
      outfile_->mkdir("cr3_mt2bl_nm1");
    }
    outfile_->mkdir("cr3_final");
  }

  if (doCR23) {
    outfile_->mkdir("cr23_presel");
    if (doNM1Plots) {
      outfile_->mkdir("cr23_bbmass_nm1");
      outfile_->mkdir("cr23_mt_nm1");
      outfile_->mkdir("cr23_mt2blfirst");
      outfile_->mkdir("cr23_met_nm1");
      outfile_->mkdir("cr23_mt2blcut");
      outfile_->mkdir("cr23_mt2bl_nm1");
    }
    outfile_->mkdir("cr23_final");
  }

  if (doCR4) {
    outfile_->mkdir("cr4_presel");
    if (doNM1Plots) {
      outfile_->mkdir("cr4_mt2bl_nm1");
      outfile_->mkdir("cr4_met_nm1");
      outfile_->mkdir("cr4_mt_nm1");
    }
    outfile_->mkdir("cr4_final");
  }

  if (doCR5) {
    outfile_->mkdir("cr5_presel");
    if (doNM1Plots) {
      outfile_->mkdir("cr5_bbmass_nm1");
      outfile_->mkdir("cr5_met_nm1");
      outfile_->mkdir("cr5_met100");
      outfile_->mkdir("cr5_mt_nm1");
      outfile_->mkdir("cr5_mt2bl_nm1");
    }
    outfile_->mkdir("cr5_final");
  }

  if (doCR5METLast) {
    outfile_->mkdir("cr5_metlast_presel");
    if (doNM1Plots) {
      outfile_->mkdir("cr5_metlast_bbmass_nm1");
      outfile_->mkdir("cr5_metlast_mt2bl_nm1");
      outfile_->mkdir("cr5_metlast_mtfirst");
      outfile_->mkdir("cr5_metlast_mt_nm1");
      outfile_->mkdir("cr5_metlast_nomt_met100");
      outfile_->mkdir("cr5_metlast_nomt_met150");
      outfile_->mkdir("cr5_metlast_nomt_met175");
      outfile_->mkdir("cr5_metlast_met_nm1");
      outfile_->mkdir("cr5_metlast_met100");
      outfile_->mkdir("cr5_metlast_met150");
    }
    outfile_->mkdir("cr5_metlast_final");
  }

  if (doCR5InvMass) {
    outfile_->mkdir("cr5_invmass_presel");
    if (doNM1Plots) {
      outfile_->mkdir("cr5_invmass_bbmass_nm1");
      outfile_->mkdir("cr5_invmass_mt2bl_nm1");
      outfile_->mkdir("cr5_invmass_mtfirst");
      outfile_->mkdir("cr5_invmass_mt_nm1");
      outfile_->mkdir("cr5_invmass_nomt_met100");
      outfile_->mkdir("cr5_invmass_nomt_met150");
      outfile_->mkdir("cr5_invmass_nomt_met175");
      outfile_->mkdir("cr5_invmass_met_nm1");
      outfile_->mkdir("cr5_invmass_met100");
      outfile_->mkdir("cr5_invmass_met150");
    }
    outfile_->mkdir("cr5_invmass_final");
  }

  if (doCR6) {
    outfile_->mkdir("cr6_presel");
    if (doNM1Plots) {
      outfile_->mkdir("cr6_met_nm1");
      outfile_->mkdir("cr6_mt_nm1");
      outfile_->mkdir("cr6_mt2bl_nm1");
    }
    outfile_->mkdir("cr6_final");
  }

  if (doCR6METLast) {
    outfile_->mkdir("cr6_metlast_presel");
    if (doNM1Plots) {
      outfile_->mkdir("cr6_metlast_bbmass_nm1");
      outfile_->mkdir("cr6_metlast_mt2bl_nm1");
      outfile_->mkdir("cr6_metlast_mt_nm1");
      outfile_->mkdir("cr6_metlast_met_nm1");
      outfile_->mkdir("cr6_metlast_met100");
      outfile_->mkdir("cr6_metlast_met150");
    }
    outfile_->mkdir("cr6_metlast_final");
  }

  if (doCR7) {
    outfile_->mkdir("cr7_presel");
    if (doNM1Plots) {
      outfile_->mkdir("cr7_met_nm1");
      outfile_->mkdir("cr7_mt_nm1");
      outfile_->mkdir("cr7_mt2bl_nm1");
    }
    outfile_->mkdir("cr7_final");
  }

  if (doCR8) {
    outfile_->mkdir("cr8_presel");
    if (doNM1Plots) {
      outfile_->mkdir("cr8_met_nm1");
      outfile_->mkdir("cr8_mt_nm1");
      outfile_->mkdir("cr8_mt2bl_nm1");
    }
    outfile_->mkdir("cr8_final");
  }

  if (doCR8METLast) {
    outfile_->mkdir("cr8_metlast_presel");
    if (doNM1Plots) {
      outfile_->mkdir("cr8_metlast_mt2bl_nm1");
      outfile_->mkdir("cr8_metlast_mtfirst");
      outfile_->mkdir("cr8_metlast_mt_nm1");
      outfile_->mkdir("cr8_metlast_met_nm1");
      outfile_->mkdir("cr8_metlast_met100");
      outfile_->mkdir("cr8_metlast_met150");
    }
    outfile_->mkdir("cr8_metlast_final");
  }

  if (doCR9) {
    outfile_->mkdir("cr9_presel");
    if (doNM1Plots) {
      outfile_->mkdir("cr9_met_nm1");
      outfile_->mkdir("cr9_mtpeak_nomet");
      outfile_->mkdir("cr9_mt_nm1");
      outfile_->mkdir("cr9_mtpeak_met");
      outfile_->mkdir("cr9_mt2bl_nm1");
      outfile_->mkdir("cr9_mtlast_mt2bl_nm1");
    }
    outfile_->mkdir("cr9_final");
  }

  if (doCR10) {
    outfile_->mkdir("cr10_presel");
    if (doNM1Plots) {
      outfile_->mkdir("cr10_met_nm1");
      outfile_->mkdir("cr10_mt_nm1");
      outfile_->mkdir("cr10_mt2bl_nm1");
    }
    outfile_->mkdir("cr10_final");
  }

  if (doCR11) {
    outfile_->mkdir("cr11_presel");
    if (doNM1Plots) {
      outfile_->mkdir("cr11_met_nm1");
      outfile_->mkdir("cr11_mt_nm1");
      outfile_->mkdir("cr11_mt2bl_nm1");
    }
    outfile_->mkdir("cr11_final");
  }

  if (doCR12) {
    outfile_->mkdir("cr12_presel");
    if (doNM1Plots) {
      outfile_->mkdir("cr12_met_nm1");
      outfile_->mkdir("cr12_mt_nm1");
      outfile_->mkdir("cr12_mt2bl_nm1");
    }
    outfile_->mkdir("cr12_final");
  }

  if (doCR13) {
    outfile_->mkdir("cr13_presel");
    if (doNM1Plots) {
      outfile_->mkdir("cr13_bbmass_nm1");
      outfile_->mkdir("cr13_mt2bl_nm1");
      outfile_->mkdir("cr13_met_nm1");
      outfile_->mkdir("cr13_met100");
      outfile_->mkdir("cr13_met150");
    }
    outfile_->mkdir("cr13_final");
  }

  if (doCR14) {
    outfile_->mkdir("cr14_presel");
    if (doNM1Plots) {
      outfile_->mkdir("cr14_mt2bl_nm1");
      outfile_->mkdir("cr14_mtfirst");
      outfile_->mkdir("cr14_mt_nm1");
      outfile_->mkdir("cr14_nomt_met100");
      outfile_->mkdir("cr14_nomt_met150");
      outfile_->mkdir("cr14_nomt_met175");
      outfile_->mkdir("cr14_met_nm1");
      outfile_->mkdir("cr14_met100");
      outfile_->mkdir("cr14_met150");
    }
    outfile_->mkdir("cr14_final");
  }

  if (doStopSel) {
    outfile_->mkdir("stop_presel");
    outfile_->mkdir("stop_met_nm1");
    outfile_->mkdir("stop_mt_nm1");
    outfile_->mkdir("stop_isotrk_nm1");
    outfile_->mkdir("stop_tauveto_nm1");
    outfile_->mkdir("stop_comp");
  }

  outfile_->cd();

  //------------------------------
  // vtx reweighting
  //------------------------------

  //  TFile* vtx_file = TFile::Open("../vtxreweight/vtxreweight_Summer12_DR53X-PU_S10_9p7ifb_Zselection.root");
  // TFile* vtx_file = TFile::Open("vtxreweight_Run2012D_Zselection.root");
  // if( vtx_file == 0 ){
  //   cout << "vtxreweight error, couldn't open vtx file. Quitting!"<< endl;
  //   exit(0);
  // }

  // TH1F* h_vtx_wgt = (TH1F*)vtx_file->Get("hratio");
  // h_vtx_wgt->SetName("h_vtx_wgt");

  //------------------------------
  // jet smearer object to do reco smearing
  //------------------------------

  std::vector<std::string> list_of_file_names;
  list_of_file_names.push_back("jetSmearData/Spring10_PtResolution_AK5PF.txt");
  list_of_file_names.push_back("jetSmearData/Spring10_PhiResolution_AK5PF.txt");
  list_of_file_names.push_back("jetSmearData/jet_resolutions.txt");
  //  JetSmearer *jetSmearer = makeJetSmearer(list_of_file_names);

  //------------------------------
  // file loop
  //------------------------------

  unsigned int nEventsChain=0;
  unsigned int nEvents = chain->GetEntries();
  nEventsChain = nEvents;
  ULong64_t nEventsTotal = 0;
  ULong64_t nEventsPass = 0;
  //  ULong64_t i_permille_old = 0;

  bool isData = name.Contains("data") ? true : false;

  while (TChainElement *currentFile = (TChainElement*)fileIter.Next()) {

    //----------------------------
    // load the stop baby tree
    //----------------------------

    cout << "[WHLooper::loop] running on file: " << currentFile->GetTitle() << endl;

    TFile *file = new TFile( currentFile->GetTitle() );
    TTree *tree = (TTree*)file->Get("t");
    stopt.Init(tree);

    //----------------------------
    // event loop
    //----------------------------

    ULong64_t nEvents = tree->GetEntriesFast();
    for(ULong64_t event = 0; event < nEvents; ++event) {
      stopt.GetEntry(event);

      //----------------------------
      // increment counters
      //----------------------------

      ++nEventsTotal;
      if (nEventsTotal%10000==0) {
	ULong64_t i_permille = (int)floor(1000 * nEventsTotal / float(nEventsChain));
	//if (i_permille != i_permille_old) {//this prints too often!
	// xterm magic from L. Vacavant and A. Cerri
	if (isatty(1)) {
	  printf("\015\033[32m ---> \033[1m\033[31m%4.1f%%"
		 "\033[0m\033[32m <---\033[0m\015", i_permille/10.);
	  fflush(stdout);
	}
	//	i_permille_old = i_permille;
      }

      //---------------------
      // skip duplicates
      //---------------------

      if( isData ) {
        DorkyEventIdentifier id = {stopt.run(), stopt.event(), stopt.lumi() };
        if (is_duplicate(id, already_seen) ){
          continue;
        }
	if (is_badLaserEvent(id,events_lasercalib) ){
	  //std::cout<< "Removed bad laser calibration event:" << run << "   " << event<<"\n";
	  continue;
	}
	if (is_badLaserEvent(id,events_hcallasercalib) ){
	  //std::cout<< "Removed bad hcal laser calibration event:" << run << "   " << event<<"\n";
	  continue;
	}
      }

      //---------------------
      // W+N partons samples: remove events with 2 gen b quarks
      //  to avoid overlap with wbb+jets
      //---------------------
      if (doWJetsOverlapRemoval && isWNjets_ && (stopt.nbs() == 2)) continue;
      else if (isWNjets_nobb_ && (stopt.nbs() == 2)) continue;
      else if (isWNjets_onlybb_ && (stopt.nbs() != 2)) continue;

      //---------------------------------------------------------------------------- 
      // ttV samples: split into 1,2 leptons, others
      //---------------------------------------------------------------------------- 

      if (isttvsl_ && (stopt.nleps() != 1) ) continue;
      else if (isttvdl_ && (stopt.nleps() != 2) ) continue;
      else if (isttvother_ && ((stopt.nleps() == 1) || (stopt.nleps() == 2)) ) continue;

      //---------------------------------------------------------------------------- 
      // determine event weight
      // make 2 example histograms of nvtx and corresponding weight
      //---------------------------------------------------------------------------- 

      const float lumi = 19.5; // full 2012
      //      const float lumi = 12.5; // 2012 periods A-C (approx..)
      //      const float lumi = 7.; // 2012 period D (approx)

      float evtweight_novtxweight = isData ? 1. : ( stopt.weight() * lumi * stopt.mgcor() );
      float evtweight = isData ? 1. : evtweight_novtxweight * stopt.nvtxweight();
      // float evtweight = evtweight_novtxweight;

      // to reweight from file - also need to comment stuff before
      // float vtxweight = vtxweight_n( stopt.nvtx(), h_vtx_wgt, isData );
      // evtweight *= vtxweight;

      // remove mgcor for ttbar mg samples
      if (isttmg_) {
	evtweight /= stopt.mgcor();
	// apply ISR weight to ttbar, if requested
	if (doISRReweight) {
	  if (stopt.ptttbar() > 250.) evtweight *= 0.8;
	  else if (stopt.ptttbar() > 150.) evtweight *= 0.9;
	  else if (stopt.ptttbar() > 120.) evtweight *= 0.95;
	}
      }
      // cross section weights for TChiwh samples:
      //  xsec (pb) * 1000 (pb to fb) * br(w->lv) 0.33 * br(h->bb) 0.58 / nevents (10000)
      float weight_lumi_br_nevents = 1.914E-02;
      if (isTChiwh_) {
	float xsecweight = 1.;
	if (name.Contains("TChiwh_125")) xsecweight = 4.8 * weight_lumi_br_nevents; 
	else if (name.Contains("TChiwh_150")) xsecweight = 2.4 * 0.5 * weight_lumi_br_nevents; // now 20k events
	else if (name.Contains("TChiwh_200")) xsecweight = 0.79 * 0.5 * weight_lumi_br_nevents; // now 20k events
	else if (name.Contains("TChiwh_250")) xsecweight = 0.32 * weight_lumi_br_nevents;
	else if (name.Contains("TChiwh_300")) xsecweight = 0.15 * weight_lumi_br_nevents;
	else if (name.Contains("TChiwh_350")) xsecweight = 0.074 * weight_lumi_br_nevents;
	else if (name.Contains("TChiwh_400")) xsecweight = 0.039 * 10 * weight_lumi_br_nevents; // only 1k events

	// reset weight here for signal MC. Ignore weight, nvtxweight, mgcor
	evtweight = xsecweight * lumi;
      }

      plot1D("h_nvtx_nosel",       stopt.nvtx(),       evtweight, h_1d_sig_presel, 40, 0, 40);
      plot1D("h_vtxweight_nosel", stopt.nvtxweight(), evtweight, h_1d_sig_presel, 41, -4., 4.);

      if (isttsl_ || isttdl_) {
	if (doTopPtReweight2) evtweight *= sqrt( TopPtWeight_v2(stopt.t().Pt()) * TopPtWeight_v2(stopt.tbar().Pt()) );
	else evtweight *= TopPtWeight(stopt.ptt());
      }

      // trigger effs
      float sltrigeff = isData ? 1. : 
	getsltrigweight(stopt.id1(), stopt.lep1().Pt(), stopt.lep1().Eta());
      float dltrigeff = isData ? 1. : 
	getdltrigweight(stopt.id1(), stopt.id2());

      float evtweight1l = evtweight * sltrigeff;
      float evtweight2l = evtweight * dltrigeff;

      //----------------------------------------------------------------------------
      // apply preselection:
      // rho 0-40 GeV, MET filters, >=1 good lepton, veto 2 leptons dR < 0.1
      //----------------------------------------------------------------------------

      if ( !passEvtSelection(name, false) ) continue;

      // data run selection: to run on certain periods
      // if (isData && (stopt.run() > 203755)) continue; // periods A-C
      // if (isData && (stopt.run() < 203768)) continue; // period D

      // require first vertex to be good (i.e. index == 0)
      //   to make sure the same vertex is used for btagging, leptons, etc
      if (stopt.indexfirstGoodVertex_()) continue;

      //----------------------------------------------------------------------------
      // Function to perform MET phi corrections on-the-fly
      // Default branches are: tree->t1metphicorr_ and tree->t1metphicorrmt_
      //----------------------------------------------------------------------------

      // pair<float, float> p_t1metphicorr = 
      // 	getPhiCorrMET( stopt.t1met10(), stopt.t1met10phi(), stopt.nvtx(), !isData);
      // t1metphicorr    = p_t1metphicorr.first;
      // t1metphicorrphi = p_t1metphicorr.second;
      // t1metphicorrmt  = getMT( stopt.lep1().Pt() , stopt.lep1().Phi() , t1metphicorr , t1metphicorrphi );  


      //----------------------------------------------------------------------------
      // ADD CODE BELOW THIS LINE
      //----------------------------------------------------------------------------

      // ----------------------------------------------------
      // gather all required variables to make all selections

      jets_.clear();
      jets_fwd_.clear();
      bjets_.clear();
      jets_csv_.clear();
      jets_idx_.clear();
      jets_fwd_idx_.clear();
      bjets_idx_.clear();
      jets_smearcorrs_.clear();
      njets_ = 0;
      njetsalleta_ = 0;
      nbjets_ = 0;
      nbjetst_ = 0;
      nbjetsl_ = 0;
      ht_ = 0.;

      for( unsigned int i = 0 ; i < stopt.pfjets().size() ; ++i ){

	LorentzVector thisjet;
	// if (doJetSmearing && !isData) thisjet = smearJet(stopt.pfjets().at(i), jetSmearer, true);
	// else thisjet = stopt.pfjets().at(i);	
	thisjet = stopt.pfjets().at(i);	

	jets_smearcorrs_.push_back(thisjet.pt()/stopt.pfjets().at(i).pt());

	// Basic jet selection
	if( thisjet.pt()<30 )  continue;
	if( fabs(thisjet.eta())>4.7 )  continue;
	//	if ( (fabs(thisjet.eta()) < 2.5) 
	//	     && (stopt.pfjets_beta2_0p5().at(i)<0.2) ) continue;
	// pileup MVA ID: use tight working point (= 0)
	//	if (!passMVAJetId(thisjet.pt(), thisjet.eta(), stopt.pfjets_mva5xPUid().at(i), 0)) continue;
	// pileup MVA ID: now use medium working point (= 1)
	if (!passMVAJetId(thisjet.pt(), thisjet.eta(), stopt.pfjets_mva5xPUid().at(i), 1)) continue;
	
	// count jets from all eta, save jets within |eta| < 2.4
	++njetsalleta_;
	ht_ += thisjet.pt();
	if (fabs(thisjet.eta()) > 2.4) {
	  jets_fwd_.push_back( thisjet );
	  jets_fwd_idx_.push_back(i);
	  continue;
	}
	++njets_;

	//	if (fabs(thisjet.eta()) < 2.4) ++njets_;

	jets_.push_back( thisjet );
	jets_idx_.push_back(i);

	// if (n_jets==1) 
	//   dphimj1 = getdphi(t1metphicorrphi, thisjet.phi() );
	// if (n_jets==2) {
	//   dphimj2 = getdphi(t1metphicorrphi, thisjet.phi() );
	//   dphimjmin = TMath::Min( dphimj1 , dphimj2 );
	// }

	float csv_nominal= stopt.pfjets_csv().at(i);

	//RESHAPING -- requires babies V20 or higher
	//  -- may need to validate for non-b-jets
	if (doCSVReshaping && !isData && !isTChiwh_) {
	  csv_nominal = nominalShape->reshape( thisjet.eta(),
					       thisjet.pt(),
					       stopt.pfjets_csv().at(i),
					       stopt.pfjets_mcflavorAlgo().at(i) ); 
	}

	jets_csv_.push_back( csv_nominal );

	//	if( (fabs(thisjet.eta()) <= 2.4) && (csv_nominal > 0.5) ) {
	if( (fabs(thisjet.eta()) <= 2.4) && (csv_nominal > getCSVCut(WHLooper::CSVM)) ) {
	  bjets_.push_back( thisjet );
	  bjets_idx_.push_back(i);
	  ++nbjets_;
	  if (csv_nominal > getCSVCut(WHLooper::CSVT)) {
	    ++nbjetst_;
	  }
	  // apply SF if requested
	  if (doBtagSFs && !isData && !isTChiwh_) {
	    float sf = getBtagSF(thisjet.pt(),thisjet.eta(),stopt.pfjets_mcflavorAlgo().at(i));
	    evtweight1l *= sf;
	    evtweight2l *= sf;
	  }
	}
	else if( (fabs(thisjet.eta()) <= 2.4) && (csv_nominal > getCSVCut(WHLooper::CSVL)) ) {
	  ++nbjetsl_;
	}

	//      	sigma_jets.push_back(stopt.pfjets_sigma().at(i));

        // //count jets that are not overlapping with second lepton
	// if (isData) continue;
	// if (stopt.nleps()!=2) continue;
	// if (stopt.mclep2().pt() < 30.) continue;
	// if (ROOT::Math::VectorUtil::DeltaR(stopt.mclep2(), jet) > 0.4 ) continue;
	// n_ljets--;

      } // loop over pfjets

      // -- set met/mt definitions for cuts
      met_ = stopt.t1metphicorr();
      metphi_ = stopt.t1metphicorrphi();
      mt_ = stopt.t1metphicorrmt();

      // -- alternate MET def: jets + leptons + remaining tracks
      // met_ = stopt.mettlj15();
      // metphi_ = stopt.mettlj15phi();
      // mt_ = getMT(stopt.lep1().pt(),stopt.lep1().phi(),stopt.mettlj15(),stopt.mettlj15phi());

      // sum Et: start from pfsumet, want type1
      //  want to add L2/L3 jet corrections: follow same procedure as type1met
      // get soft component of met - remove lepton and all passing jets 
      TVector2 lep(stopt.lep1().px(),stopt.lep1().py());
      TVector2 met;
      met.SetMagPhi(met_, metphi_);
      TVector2 w = lep+met;
      wpt_ = w.Mod();
      TVector2 pfmet_soft(met);
      pfmet_soft += lep;
      lepmetdphi_ = fabs(TVector2::Phi_mpi_pi(stopt.lep1().phi() - metphi_));
      sumet_ = stopt.pfsumet();
      sumet_soft_ = stopt.pfsumet() - stopt.lep1().pt();
      for (unsigned int i = 0; i < jets_.size(); ++i) {
	int jet_idx = jets_idx_[i];

	TVector2 jet;
	// jet pt in MET: jetpt_l2l3 = jetpt_uncor + jetpt_l1l2l3cor - jetpt_l1cor
	float jetpt_uncor = jets_[i].pt() / stopt.pfjets_corr()[jet_idx];
	float jetpt_l2l3 = jets_[i].pt() + jetpt_uncor - jetpt_uncor * stopt.pfjets_l1corr()[jet_idx]; 
	jet.SetMagPhi(jetpt_l2l3, jets_[i].phi());
	pfmet_soft += jet;

	sumet_ += jetpt_l2l3 - jetpt_uncor * stopt.pfjets_l1corr()[jet_idx];
	sumet_soft_ -= jetpt_uncor;

      } // loop over jets
      met_soft_ = pfmet_soft.Mod();

      // calculate mt2 vars after selecting jets -- require at least 2 bjets here to avoid wasting time..
      bb_ = LorentzVector();
      mt2b_ = -1.;
      mt2bl_ = -1.;
      mt2w_ = -1.;
      if (nbjets_ >= 2) {
	bb_ = bjets_.at(0) + bjets_.at(1);
	mt2b_ = calculateMT2w(jets_, jets_csv_, stopt.lep1(), met_, metphi_, MT2b);
	mt2bl_ = calculateMT2w(jets_, jets_csv_, stopt.lep1(), met_, metphi_, MT2bl);
	mt2w_ = calculateMT2w(jets_, jets_csv_, stopt.lep1(), met_, metphi_, MT2w);
      }
      else if ((nbjets_ == 1) && (njets_ >= 2)) {
	// 1 bjet: use bjet + highest pt other jet for dijet mass
	if (bjets_idx_.at(0) == jets_idx_.at(0)) bb_ = bjets_.at(0) + jets_.at(1);
	else bb_ = bjets_.at(0) + jets_.at(1);
      }
      else if (njets_ >= 2) {
	// 0 bjets: use two highest pt jets for invariant mass
	bb_ = jets_.at(0) + jets_.at(1);
      }

      // gen MT2bl for ttbar MC
      if (isttsl_ || isttdl_ || (isWjets_ && stopt.genbs().size() >= 2)) {
	std::vector<float> genbs_csv(stopt.genbs().size(), 0.99);
	genmt2bl_ = calculateMT2w(stopt.genbs(), genbs_csv, stopt.mclep1(), stopt.genmet(), stopt.genmetphi(), MT2bl);
      }
      bbwdphi_ = fabs(TVector2::Phi_mpi_pi(bb_.phi() - w.Phi()));

      // TVector2 lep(stopt.lep1().px(),stopt.lep1().py());
      // TVector2 met;
      // met.SetMagPhi(stopt.t1metphicorr(), metphi_);
      // TVector2 w = lep+met; 

      // end variables --------------------------------------
      // ----------------------------------------------------

      // weight for MT tail in Wbb sample to account for missing off-shell W contribution
      if (doWbbMtReweight && isWjets_ && !isWNjets_) {
	// weight up events with reco MT > 100 by 10%
	if (mt_ > 100. && met_ > 100.) {
	  evtweight1l *= 1.3;
	  evtweight2l *= 1.3;
	}
	else if (mt_ > 100. && met_ > 50.) {
	  evtweight1l *= 1.1;
	  evtweight2l *= 1.1;
	}
      }

      // ----------------------------------------------------
      // selections bits

      //      bool passisotrk = passIsoTrkVeto_v2();
      //      bool passisotrk = passIsoTrkVeto_v3() && (stopt.ngoodlep() == 1);
      bool passisotrk = passIsoTrkVeto_v4();
      if (doTauVeto) passisotrk &= passTauVeto();
      if (doLep2Veto) passisotrk &= (stopt.ngoodlep() == 1);

      // end selections bits --------------------------------
      // ----------------------------------------------------

      // ----------------------------------------------------
      // region selection and plots

      // always require at least 2 (central) jets
      if (njets_ < 2) continue;

      // require lead jet pt > 50 GeV
      if (jets_.at(0).pt() < 50.) continue;

      // try tightening isolation to see CR agreement
      //      if ( stopt.isopf1() > 0.1 ) continue;

      // -------------------------------------------
      // *** inclusive preselection:
      //   >= 1 lepton
      //   >= 2 (central) jets
      //   met > 50
      //
      //   then add btags

      if ( doInclusive
	   && passSingleLeptonSelection(isData) 
	   && (njets_ >= 2)
	   && (met_ > CUT_MET_PRESEL_) ) {

        fillHists1DWrapper(h_1d_inc_presel,evtweight1l,"inc_presel");

	if ( (njetsalleta_ == 2) ) {
	  fillHists1DWrapper(h_1d_inc_2j,evtweight1l,"inc_2j");
	}

	if (doInclusiveMTTail) {
	  if ( (njetsalleta_ == 2) && (mt_ > CUT_MT_) ) {
	    fillHists1DWrapper(h_1d_inc_2j_mt,evtweight1l,"inc_2j_mt");
	  }

	  if ( (njetsalleta_ == 2) && (mt_ > CUT_MT_) && (met_ > 100.) ) {
	    fillHists1DWrapper(h_1d_inc_2j_mt_met100,evtweight1l,"inc_2j_mt_met100");
	  }

	  if ( (njetsalleta_ == 2) && (mt_ > CUT_MT_) && (met_ > 150.) ) {
	    fillHists1DWrapper(h_1d_inc_2j_mt_met150,evtweight1l,"inc_2j_mt_met150");
	  }

	  if ( (njetsalleta_ == 2) && (mt_ > CUT_MT_) && (met_ > CUT_MET_) ) {
	    fillHists1DWrapper(h_1d_inc_2j_mt_metcut,evtweight1l,"inc_2j_mt_metcut");
	  }
	}

	bool fail = false;
	if ( !fail && (nbjets_ >= 1) ) {
	  fillHists1DWrapper(h_1d_inc_1b,evtweight1l,"inc_1b");
	}
	else fail = true;

	if (!fail && (nbjets_ >= 2 ) ) {
	  fillHists1DWrapper(h_1d_inc_2b,evtweight1l,"inc_2b");
	}
	else fail = true;
      }

      // -------------------------------------------
      // *** presel for signal region:
      //   == 1 lepton, iso track veto
      //   >= 2 bjets
      //   100 < m(bb) < 150
      //   met > 50

      // *** signal region:
      //   == 1 lepton, iso track veto
      //   == 2 jets, all eta
      //   == 2 bjets
      //   100 < m(bb) < 150
      //   met > 175 (cut at 50 for presel)
      //   mt > 100
      //   mt2bl > 200

      if ( doSignal
	   && passSingleLeptonSelection(isData) 
	   && passisotrk 
	   && (nbjets_ >= 2)
	   && (bb_.M() > CUT_BBMASS_LOW_) && (bb_.M() < CUT_BBMASS_HIGH_)
	   && (met_ > CUT_MET_PRESEL_) 
	   && (mt_ > CUT_MT_PRESEL_) 
	   && (!isData || !blindSignal) ) {

        fillHists1DWrapper(h_1d_sig_presel,evtweight1l,"sig_presel");

	bool fail = false;
	if ( !fail && (njetsalleta_ == 2) ) {
	  if (doNM1Plots) fillHists1DWrapper(h_1d_sig_met_nm1,evtweight1l,"sig_met_nm1");
	}
	else fail = true;

	// mt peak before cutting on met
	if (!fail && (mt_ > 50.) && (mt_ < 80.) ) {
	  if (doNM1Plots) fillHists1DWrapper(h_1d_sig_mtpeak_nomet,evtweight1l,"sig_mtpeak_nomet");
	}

	if (!fail && (met_ > CUT_MET_) ) {
	  if (doNM1Plots) fillHists1DWrapper(h_1d_sig_mt_nm1,evtweight1l,"sig_mt_nm1");
	}
	else fail = true;

	// mt peak after cutting on met
	if (!fail && (mt_ > 50.) && (mt_ < 80.) ) {
	  if (doNM1Plots) fillHists1DWrapper(h_1d_sig_mtpeak_met,evtweight1l,"sig_mtpeak_met");
	}

	if (!fail && (mt_ > CUT_MT_) ) {
	  if (doNM1Plots) fillHists1DWrapper(h_1d_sig_mt2bl_nm1,evtweight1l,"sig_mt2bl_nm1");
	}
	else fail = true;

	if (!fail && (mt2bl_ > CUT_MT2BL_) ) {
	  fillHists1DWrapper(h_1d_sig_final,evtweight1l,"sig_final");
	}

      } // signal region sel

      // -------------------------------------------
      // *** signal region, but applying m(bb) cut last
      //  - presel:
      //   == 1 lepton, iso track veto
      //   >= 2 bjets
      //   met > 50

      if ( doSignalMassLast
	   && passSingleLeptonSelection(isData) 
	   && passisotrk 
	   && (nbjets_ >= 2)
	   && (met_ > CUT_MET_PRESEL_) 
	   && (mt_ > CUT_MT_PRESEL_) ) {
	   //	   && (!isData || !blindSignal) ) {

	// for data: remove mass window, so we can compare bb mass shape outside signal region as cuts are applied
	bool fail = false;
	if (isData && blindSignal && (bb_.M() > CUT_BBMASS_LOW_) && (bb_.M() < CUT_BBMASS_CR1_LOW_) ) fail = true;

        if (!fail) fillHists1DWrapper(h_1d_sig_bbmasslast_presel,evtweight1l,"sig_bbmasslast_presel");

	if ( !fail && (njetsalleta_ == 2) ) {
	  if (doNM1Plots) fillHists1DWrapper(h_1d_sig_bbmasslast_mt2bl_nm1,evtweight1l,"sig_bbmasslast_mt2bl_nm1");
	}
	else fail = true;

	// mt peak before cutting on met/mt/mt2bl
	if (!fail && (mt_ > 50.) && (mt_ < 80.) ) {
	  if (doNM1Plots) fillHists1DWrapper(h_1d_sig_bbmasslast_mtpeak,evtweight1l,"sig_bbmasslast_mtpeak");
	}

	// plots for bbmass after each of the major cuts done separately
	if (!fail && (mt_ > CUT_MT_) ) {
	  if (doNM1Plots) fillHists1DWrapper(h_1d_sig_bbmasslast_mtcut,evtweight1l,"sig_bbmasslast_mtcut");
	  if (!fail && (met_ > 100.) ) {
	    if (doNM1Plots) fillHists1DWrapper(h_1d_sig_bbmasslast_mtcut_met100,evtweight1l,"sig_bbmasslast_mtcut_met100");
	  }
	}

	if (!fail && (met_ > CUT_MET_) ) {
	  if (doNM1Plots) fillHists1DWrapper(h_1d_sig_bbmasslast_metcut,evtweight1l,"sig_bbmasslast_metcut");
	}

	if (!fail && (mt2bl_ > CUT_MT2BL_) ) {
	  if (doNM1Plots) fillHists1DWrapper(h_1d_sig_bbmasslast_mt_nm1,evtweight1l,"sig_bbmasslast_mt_nm1");
	}
	else fail = true;

	if (!fail && (mt_ > CUT_MT_) ) {
	  if (doNM1Plots) fillHists1DWrapper(h_1d_sig_bbmasslast_met_nm1,evtweight1l,"sig_bbmasslast_met_nm1");
	}
	else fail = true;

	if (!fail && (met_ > 100.) ) {
	  fillHists1DWrapper(h_1d_sig_bbmasslast_met100,evtweight1l,"sig_bbmasslast_met100");
	}
	else fail = true;

	if (!fail && (met_ > 150.) ) {
	  fillHists1DWrapper(h_1d_sig_bbmasslast_met150,evtweight1l,"sig_bbmasslast_met150");
	}
	else fail = true;

	if (!fail && (met_ > CUT_MET_) ) {
	  fillHists1DWrapper(h_1d_sig_bbmasslast_bbmass_nm1,evtweight1l,"sig_bbmasslast_bbmass_nm1");
	}
	else fail = true;

	if (!fail && (bb_.M() > CUT_BBMASS_LOW_) && (bb_.M() < CUT_BBMASS_HIGH_) ) {
	  fillHists1DWrapper(h_1d_sig_bbmasslast_final,evtweight1l,"sig_bbmasslast_final");
	}

      } // signal region sel (bbmass last)

      // -------------------------------------------
      // *** signal region, but applying met cut last
      //  - presel:
      //   == 1 lepton, iso track veto
      //   >= 2 bjets
      //   100 < m(bb) < 150
      //   met > 50

      if ( doSignalMETLast
	   && passSingleLeptonSelection(isData) 
	   && passisotrk 
	   && (nbjets_ >= 2)
	   && (bb_.M() > CUT_BBMASS_LOW_) && (bb_.M() < CUT_BBMASS_HIGH_)
	   && (met_ > CUT_MET_PRESEL_) 
	   && (mt_ > CUT_MT_PRESEL_) 
	   && (!isData || !blindSignal) ) {

        fillHists1DWrapper(h_1d_sig_metlast_presel,evtweight1l,"sig_metlast_presel");

	bool fail = false;
	if ( !fail && (njetsalleta_ == 2) ) {
	  if (doNM1Plots) fillHists1DWrapper(h_1d_sig_metlast_mt2bl_nm1,evtweight1l,"sig_metlast_mt2bl_nm1");
	  if (doNM1Plots) fillHists2D(h_2d_sig_metlast_mt2bl_nm1,evtweight1l,"sig_metlast_mt2bl_nm1");
	}
	else fail = true;

	if (!fail && (mt_ > CUT_MT_) ) {
	  if (doNM1Plots) fillHists1DWrapper(h_1d_sig_metlast_mtfirst,evtweight1l,"sig_metlast_mtfirst");
	}

	if (!fail && (mt2bl_ > CUT_MT2BL_) ) {
	  if (doNM1Plots) fillHists1DWrapper(h_1d_sig_metlast_mt_nm1,evtweight1l,"sig_metlast_mt_nm1");
	}
	else fail = true;

	if (!fail && (mt_ > CUT_MT_) ) {
	  if (doNM1Plots) fillHists1DWrapper(h_1d_sig_metlast_met_nm1,evtweight1l,"sig_metlast_met_nm1");
	}
	else fail = true;

	if (!fail && (met_ > 100.) ) {
	  if (doNM1Plots) fillHists1DWrapper(h_1d_sig_metlast_met100,evtweight1l,"sig_metlast_met100");
	}
	else fail = true;

	if (!fail && (met_ > 150.) ) {
	  if (doNM1Plots) fillHists1DWrapper(h_1d_sig_metlast_met150,evtweight1l,"sig_metlast_met150");
	}
	else fail = true;

	if (!fail && (met_ > CUT_MET_) ) {
	  fillHists1DWrapper(h_1d_sig_metlast_final,evtweight1l,"sig_metlast_final");
	}
	else fail = true;

      } // signal region sel, met last

      // -------------------------------------------
      // *** signal region with SMWH cuts:
      //   == 1 lepton, iso track veto
      //   == 2 jets, all eta
      //   == 2 bjets
      //   100 < m(bb) < 150
      //   met > 45
      //   mt > 100
      //   mt2bl > 200

      if ( doSignalSMWH
	   && passSingleLeptonSelection(isData) 
	   && passisotrk 
	   && (nbjets_ >= 2)
	   && (nbjetst_ >= 1)
	   && (bb_.M() > CUT_BBMASS_LOW_) && (bb_.M() < CUT_BBMASS_HIGH_)
	   && (met_ > 45.) 
	   && (!isData || !blindSignal) ) {

        fillHists1DWrapper(h_1d_sigsmwh_presel,evtweight1l,"sigsmwh_presel");

	bool fail = false;
	if ( !fail && (njetsalleta_ == 2) ) {
	  //	  if (doNM1Plots) fillHists1DWrapper(h_1d_sig_met_nm1,evtweight1l,"sig_met_nm1");
	}
	else fail = true;

	if (!fail && (bb_.pt() > 100.) ) {
	  //	  if (doNM1Plots) fillHists1DWrapper(h_1d_sig_mt_nm1,evtweight1l,"sig_mt_nm1");
	}
	else fail = true;

	if (!fail && (bbwdphi_ > 2.95) ) {
	  //	  if (doNM1Plots) fillHists1DWrapper(h_1d_sig_mt_nm1,evtweight1l,"sig_mt_nm1");
	}
	else fail = true;

	if (!fail && (wpt_ > 180.) ) {
	  //	  if (doNM1Plots) fillHists1DWrapper(h_1d_sig_mt_nm1,evtweight1l,"sig_mt_nm1");
	}
	else fail = true;

	if (!fail && (lepmetdphi_ < TMath::Pi()/2.) ) {
	  fillHists1DWrapper(h_1d_sigsmwh_final,evtweight1l,"sigsmwh_final");
	}

      } // signal region sel

      // -------------------------------------------
      // *** CR1: 150 < m(bb) < 250 
      //  otherwise same as signal region

      if ( doCR1
	   && passSingleLeptonSelection(isData) 
	   && passisotrk 
	   && (nbjets_ >= 2)
	   && (bb_.M() > CUT_BBMASS_CR1_LOW_) && (bb_.M() < CUT_BBMASS_CR1_HIGH_) 
	   && (met_ > CUT_MET_PRESEL_) 
	   && (mt_ > CUT_MT_PRESEL_) ) {

        fillHists1DWrapper(h_1d_cr1_presel,evtweight1l,"cr1_presel");

	bool fail = false;
	if ( !fail && (njetsalleta_ == 2) ) {
	  if (doNM1Plots) fillHists1DWrapper(h_1d_cr1_met_nm1,evtweight1l,"cr1_met_nm1");
	}
	else fail = true;

	//	if (!fail && isData && (met_ > 420.)) dumpEventInfo("CR1 high MET event");

	// mt peak before cutting on met
	if (!fail && (mt_ > 50.) && (mt_ < 80.) ) {
	  if (doNM1Plots) fillHists1DWrapper(h_1d_cr1_mtpeak_nomet,evtweight1l,"cr1_mtpeak_nomet");
	}

	if (!fail && (met_ > CUT_MET_) ) {
	  if (doNM1Plots) fillHists1DWrapper(h_1d_cr1_mt_nm1,evtweight1l,"cr1_mt_nm1");
	}
	else fail = true;

	// mt peak after cutting on met
	if (!fail && (mt_ > 50.) && (mt_ < 80.) ) {
	  if (doNM1Plots) fillHists1DWrapper(h_1d_cr1_mtpeak_met,evtweight1l,"cr1_mtpeak_met");
	}

	if (!fail && (mt2bl_ > CUT_MT2BL_) ) {
	  if (doNM1Plots) fillHists1DWrapper(h_1d_cr1_mtlast_mt2bl_nm1,evtweight1l,"cr1_mtlast_mt2bl_nm1");
	}
	// else fail = true;

	if (!fail && (mt_ > CUT_MT_) ) {
	  if (doNM1Plots) fillHists1DWrapper(h_1d_cr1_mt2bl_nm1,evtweight1l,"cr1_mt2bl_nm1");
	}
	else fail = true;

	if (!fail && (mt2bl_ > CUT_MT2BL_) ) {
	  fillHists1DWrapper(h_1d_cr1_final,evtweight1l,"cr1_final");
	}

      } // CR1 region sel


      // -------------------------------------------
      // *** CR1, MET last: 150 < m(bb) < 250 
      //  otherwise same as signal region

      if ( doCR1METLast
	   && passSingleLeptonSelection(isData) 
	   && passisotrk 
	   && (nbjets_ >= 2)
	   && (bb_.M() > CUT_BBMASS_CR1_LOW_) && (bb_.M() < CUT_BBMASS_CR1_HIGH_) 
	   && (met_ > CUT_MET_PRESEL_) 
	   && (mt_ > CUT_MT_PRESEL_) ) {

        fillHists1DWrapper(h_1d_cr1_metlast_presel,evtweight1l,"cr1_metlast_presel");

	if (isData && (stopt.ngoodlep() >= 2)) dumpEventInfo("CR1 2lep event");

	bool fail = false;
	if ( !fail && (njetsalleta_ == 2) ) {
	  if (doNM1Plots) fillHists1DWrapper(h_1d_cr1_metlast_mt2bl_nm1,evtweight1l,"cr1_metlast_mt2bl_nm1");
	}
	else fail = true;

	if (!fail && (mt_ > CUT_MT_) ) {
	  fillHists1DWrapper(h_1d_cr1_metlast_mtfirst,evtweight1l,"cr1_metlast_mtfirst");
	}

	if (!fail && (mt2bl_ > CUT_MT2BL_) ) {
	  if (doNM1Plots) fillHists1DWrapper(h_1d_cr1_metlast_mt_nm1,evtweight1l,"cr1_metlast_mt_nm1");
	}
	else fail = true;

	if (!fail && (mt_ > CUT_MT_) ) {
	  fillHists1DWrapper(h_1d_cr1_metlast_met_nm1,evtweight1l,"cr1_metlast_met_nm1");
	}
	else fail = true;

	if (!fail && (met_ > 100.) ) {
	  if (doNM1Plots) fillHists1DWrapper(h_1d_cr1_metlast_met100,evtweight1l,"cr1_metlast_met100");
	}
	else fail = true;

	if (!fail && (met_ > 150.) ) {
	  if (doNM1Plots) fillHists1DWrapper(h_1d_cr1_metlast_met150,evtweight1l,"cr1_metlast_met150");
	}
	else fail = true;

	if (!fail && (met_ > CUT_MET_) ) {
	  if (doNM1Plots) fillHists1DWrapper(h_1d_cr1_metlast_final,evtweight1l,"cr1_metlast_final");
	}
	else fail = true;

      } // CR1 region sel


      // -------------------------------------------
      // *** CR2: 1 lepton + iso track, i.e. inverted iso track (or tau) veto, require exactly 1 lep
      //   otherwise same as signal region

      if ( doCR2
	   && passLepPlusIsoTrkSelectionWHMet(isData)
	   // && ( (passSingleLeptonSelection(isData) && !passisotrk)
	   // // && passLepPlusIsoTrkSelection(isData) 
	   // 	|| passLepPlusTauSelection(isData) )
	   //	   && (stopt.ngoodlep() == 1) 
	   && (nbjets_ >= 2)
	   //	   && (bb_.M() > CUT_BBMASS_LOW_) && (bb_.M() < CUT_BBMASS_HIGH_)
	   && (met_ > CUT_MET_PRESEL_) 
	   && (mt_ > CUT_MT_PRESEL_) ) {

        fillHists1DWrapper(h_1d_cr2_presel,evtweight1l,"cr2_presel");

	bool fail = false;
	if ( !fail && (njetsalleta_ == 2) ) {
	  if (doNM1Plots) fillHists1DWrapper(h_1d_cr2_bbmass_nm1,evtweight1l,"cr2_bbmass_nm1");
	}
	else fail = true;

	if (!fail && (bb_.M() > CUT_BBMASS_LOW_) && (bb_.M() < CUT_BBMASS_HIGH_) ) {
	  if (doNM1Plots) fillHists1DWrapper(h_1d_cr2_mt_nm1,evtweight1l,"cr2_mt_nm1");
	}
	else fail = true;

	if (!fail && (mt2bl_ > CUT_MT2BL_) ) {
	  if (doNM1Plots) fillHists1DWrapper(h_1d_cr2_mt2blfirst,evtweight1l,"cr2_mt2blfirst");
	}

	if (!fail && (mt_ > CUT_MT_) ) {
	  fillHists1DWrapper(h_1d_cr2_met_nm1,evtweight1l,"cr2_met_nm1");
	}
	else fail = true;

	//	if (!fail && isData && (met_ > 250. || mt2bl_ > 250.)) dumpEventInfo("CR2, high MET/MT2bl event");

	if (!fail && (mt2bl_ > CUT_MT2BL_) ) {
	  if (doNM1Plots) fillHists1DWrapper(h_1d_cr2_mt2blcut,evtweight1l,"cr2_mt2blcut");
	}

	if (!fail && (met_ > CUT_MET_) ) {
	  if (doNM1Plots) fillHists1DWrapper(h_1d_cr2_mt2bl_nm1,evtweight1l,"cr2_mt2bl_nm1");
	}
	else fail = true;

	if (!fail && (mt2bl_ > CUT_MT2BL_) ) {
	  if (doNM1Plots) fillHists1DWrapper(h_1d_cr2_final,evtweight1l,"cr2_final");
	}
	else fail = true;

      } // CR2 region sel


      // -------------------------------------------
      // *** CR3: 2 leptons, Z mass veto for same flavor, no iso track veto
      //   otherwise same as signal region
      //   !!! pesudo MET/MT etc defs to emulate losing 2nd lepton

      if ( doCR3
	   && passDileptonSelection(isData) 
	   //	   && (abs(stopt.id1()) != abs(stopt.id2()) || fabs( stopt.dilmass() - 91.) > 15. )
	   && (nbjets_ >= 2) ) {
	//	   && (bb_.M() > CUT_BBMASS_LOW_) && (bb_.M() < CUT_BBMASS_HIGH_) ) {

	// tight 3rd track veto used by stop people -- necessary?
	//	      if ( (stopt.trkpt10loose() <0.0001 || stopt.trkreliso10loose() > 0.1) 

	//calculate pseudo met and mt
	//find positive lepton - this is the one that is combined with the pseudomet to form the mT
	bool isfirstp = (stopt.id1() > 0) ? true : false;
	lep_ = isfirstp ? stopt.lep1() : stopt.lep2();
		
	//recalculate met
	float metx = met_ * cos( metphi_ );
	float mety = met_ * sin( metphi_ );
		
	//recalculate the MET with the positive lepton
	metx += isfirstp ? stopt.lep1().px() : stopt.lep2().px();
	mety += isfirstp ? stopt.lep1().py() : stopt.lep2().py();
		
	pseudomet_lep_    = sqrt(metx*metx + mety*mety);
	pseudometphi_lep_ = atan2( mety , metx );
		
	//recalculate the MT with the negative lepton
	pseudomt_lep_ = getMT( lep_.pt() , lep_.phi() , pseudomet_lep_ , pseudometphi_lep_ );
	//dphi between met and lepton
	dphi_pseudomet_lep_ = TVector2::Phi_mpi_pi( lep_.phi() - pseudometphi_lep_ );

	// recalculate mt2 vars also..
	pseudomt2b_ = calculateMT2w(jets_, jets_csv_, lep_, pseudomet_lep_, pseudometphi_lep_, MT2b);
	pseudomt2bl_ = calculateMT2w(jets_, jets_csv_, lep_, pseudomet_lep_, pseudometphi_lep_, MT2bl);
	pseudomt2w_ = calculateMT2w(jets_, jets_csv_, lep_, pseudomet_lep_, pseudometphi_lep_, MT2w);


	bool fail = false;
	//	if ( pseudomet_lep_ > CUT_MET_PRESEL_ ) {
	if ( met_ > CUT_MET_PRESEL_ ) {
          fillHists1DWrapper(h_1d_cr3_presel,evtweight2l,"cr3_presel");
	}
	else fail = true;

	if ( !fail && (njetsalleta_ == 2) ) {
	  if (doNM1Plots) fillHists1DWrapper(h_1d_cr3_bbmass_nm1,evtweight2l,"cr3_bbmass_nm1");
	}
	else fail = true;

	if ( !fail &&  (bb_.M() > CUT_BBMASS_LOW_) && (bb_.M() < CUT_BBMASS_HIGH_) ) {
	  if (doNM1Plots) fillHists1DWrapper(h_1d_cr3_mt_nm1,evtweight2l,"cr3_mt_nm1");
	}
	else fail = true;

	if (!fail && (mt2bl_ > CUT_MT2BL_) ) {
	  if (doNM1Plots) fillHists1DWrapper(h_1d_cr3_mt2blfirst,evtweight1l,"cr3_mt2blfirst");
	}

	//	if (!fail && (pseudomt_lep_ > CUT_MT_) ) {
	if (!fail && (mt_ > CUT_MT_) ) {
	  fillHists1DWrapper(h_1d_cr3_met_nm1,evtweight2l,"cr3_met_nm1");
	}
	else fail = true;

	if (!fail && (mt2bl_ > CUT_MT2BL_) ) {
	  if (doNM1Plots) fillHists1DWrapper(h_1d_cr3_mt2blcut,evtweight2l,"cr3_mt2blcut");
	}

	//	if (!fail && (pseudomet_lep_ > CUT_MET_) ) {
	if (!fail && (met_ > CUT_MET_) ) {
	  if (doNM1Plots) fillHists1DWrapper(h_1d_cr3_mt2bl_nm1,evtweight2l,"cr3_mt2bl_nm1");
	}
	else fail = true;

	//	if (!fail && isData && met_ > 270.) dumpEventInfo("CR3, high MET event");

	//	if (!fail && (pseudomt2bl_ > CUT_MT2BL_) ) {
	if (!fail && (mt2bl_ > CUT_MT2BL_) ) {
	  if (doNM1Plots) fillHists1DWrapper(h_1d_cr3_final,evtweight2l,"cr3_final");
	}
	else fail = true;

      } // CR3 region sel


      // -------------------------------------------
      // *** CR23: 1 lepton + iso track OR 2 lepton
      //   otherwise same as signal region

      if ( doCR23
	   && ( passLepPlusIsoTrkSelectionWHMet(isData)
	   //	   && ( ( passSingleLeptonSelection(isData) && !passisotrk && (stopt.ngoodlep() == 1) ) 
	   // && passLepPlusIsoTrkSelection(isData) 
		//		|| ( passLepPlusTauSelection(isData) && (stopt.ngoodlep() == 1) )
		|| ( passDileptonSelection(isData) ) )
		     //&& (abs(stopt.id1()) != abs(stopt.id2()) || fabs( stopt.dilmass() - 91.) > 15. ) ) )
	   && (nbjets_ >= 2)
	   //	   && (bb_.M() > CUT_BBMASS_LOW_) && (bb_.M() < CUT_BBMASS_HIGH_)
	   && (met_ > CUT_MET_PRESEL_) 
	   && (mt_ > CUT_MT_PRESEL_) ) {

        fillHists1DWrapper(h_1d_cr23_presel,evtweight1l,"cr23_presel");

	bool fail = false;
	if ( !fail && (njetsalleta_ == 2) ) {
	  if (doNM1Plots) fillHists1DWrapper(h_1d_cr23_bbmass_nm1,evtweight1l,"cr23_bbmass_nm1");
	}
	else fail = true;

	if (!fail && (bb_.M() > CUT_BBMASS_LOW_) && (bb_.M() < CUT_BBMASS_HIGH_) ) {
	  if (doNM1Plots) fillHists1DWrapper(h_1d_cr23_mt_nm1,evtweight1l,"cr23_mt_nm1");
	}
	else fail = true;

	if (!fail && (mt2bl_ > CUT_MT2BL_) ) {
	  if (doNM1Plots) fillHists1DWrapper(h_1d_cr23_mt2blfirst,evtweight1l,"cr23_mt2blfirst");
	}

	if (!fail && (mt_ > CUT_MT_) ) {
	  fillHists1DWrapper(h_1d_cr23_met_nm1,evtweight1l,"cr23_met_nm1");
	}
	else fail = true;

	//	if (!fail && isData && (met_ > 250. || mt2bl_ > 250.)) dumpEventInfo("CR23, high MET/MT2bl event");

	if (!fail && (mt2bl_ > CUT_MT2BL_) ) {
	  if (doNM1Plots) fillHists1DWrapper(h_1d_cr23_mt2blcut,evtweight1l,"cr23_mt2blcut");
	}

	if (!fail && (met_ > CUT_MET_) ) {
	  if (doNM1Plots) fillHists1DWrapper(h_1d_cr23_mt2bl_nm1,evtweight1l,"cr23_mt2bl_nm1");
	}
	else fail = true;

	if (!fail && (mt2bl_ > CUT_MT2BL_) ) {
	  if (doNM1Plots) fillHists1DWrapper(h_1d_cr23_final,evtweight1l,"cr23_final");
	}
	else fail = true;

      } // CR23 region sel


      // -------------------------------------------
      // *** CR4: 2 leptons, Z mass veto for same flavor, no iso track veto
      //   high mass: 150 < m(bb) < 250
      //   otherwise same as signal region
      //   !!! modified MET/MT etc defs to emulate losing 2nd lepton

      if ( doCR4
	   && passDileptonSelection(isData) 
	   && (abs(stopt.id1()) != abs(stopt.id2()) || fabs( stopt.dilmass() - 91.) > 15. )
	   && (nbjets_ >= 2)
	   && (bb_.M() > CUT_BBMASS_CR1_LOW_) && (bb_.M() < CUT_BBMASS_CR1_HIGH_) ) {

	// tight 3rd track veto used by stop people -- necessary?
	//	      if ( (stopt.trkpt10loose() <0.0001 || stopt.trkreliso10loose() > 0.1) 

	//calculate pseudo met and mt
	//find positive lepton - this is the one that is combined with the pseudomet to form the mT
	bool isfirstp = (stopt.id1() > 0) ? true : false;
	lep_ = isfirstp ? stopt.lep1() : stopt.lep2();
		
	//recalculate met
	float metx = met_ * cos( metphi_ );
	float mety = met_ * sin( metphi_ );
		
	//recalculate the MET with the positive lepton
	metx += isfirstp ? stopt.lep1().px() : stopt.lep2().px();
	mety += isfirstp ? stopt.lep1().py() : stopt.lep2().py();
		
	pseudomet_lep_    = sqrt(metx*metx + mety*mety);
	pseudometphi_lep_ = atan2( mety , metx );
		
	//recalculate the MT with the negative lepton
	pseudomt_lep_ = getMT( lep_.pt() , lep_.phi() , pseudomet_lep_ , pseudometphi_lep_ );
	//dphi between met and lepton
	dphi_pseudomet_lep_ = TVector2::Phi_mpi_pi( lep_.phi() - pseudometphi_lep_ );

	// recalculate mt2 vars also..
	pseudomt2b_ = calculateMT2w(jets_, jets_csv_, lep_, pseudomet_lep_, pseudometphi_lep_, MT2b);
	pseudomt2bl_ = calculateMT2w(jets_, jets_csv_, lep_, pseudomet_lep_, pseudometphi_lep_, MT2bl);
	pseudomt2w_ = calculateMT2w(jets_, jets_csv_, lep_, pseudomet_lep_, pseudometphi_lep_, MT2w);


	bool fail = false;
	//	if ( pseudomet_lep_ > CUT_MET_PRESEL_ ) {
	if ( met_ > CUT_MET_PRESEL_ ) {
          fillHists1DWrapper(h_1d_cr4_presel,evtweight2l,"cr4_presel");
	}
	else fail = true;

	if ( !fail && (njetsalleta_ == 2) ) {
	  if (doNM1Plots) fillHists1DWrapper(h_1d_cr4_mt2bl_nm1,evtweight2l,"cr4_mt2bl_nm1");
	}
	else fail = true;

	//	if (!fail && isData && met_ > 320.) dumpEventInfo("CR4, high MET event");

	//	if (!fail && (pseudomt2bl_ > CUT_MT2BL_) ) {
	if (!fail && (mt2bl_ > CUT_MT2BL_) ) {
	  fillHists1DWrapper(h_1d_cr4_met_nm1,evtweight2l,"cr4_met_nm1");
	}
	else fail = true;

	//	if (!fail && (pseudomet_lep_ > CUT_MET_) ) {
	if (!fail && (met_ > CUT_MET_) ) {
	  if (doNM1Plots) fillHists1DWrapper(h_1d_cr4_mt_nm1,evtweight2l,"cr4_mt_nm1");
	}
	else fail = true;

	//	if (!fail && (pseudomt_lep_ > CUT_MT_) ) {
	if (!fail && (mt_ > CUT_MT_) ) {
	  if (doNM1Plots) fillHists1DWrapper(h_1d_cr4_final,evtweight2l,"cr4_final");
	}
	else fail = true;

      } // CR4 region sel


      // -------------------------------------------
      // *** CR5: require exactly 0 btags
      //  otherwise same as signal

      if ( doCR5
	   && passSingleLeptonSelection(isData) 
	   && passisotrk 
	   && (nbjets_ == 0)
	   //	   && (bb_.M() > CUT_BBMASS_LOW_) && (bb_.M() < CUT_BBMASS_HIGH_)
	   && (met_ > CUT_MET_PRESEL_) 
	   && (mt_ > CUT_MT_PRESEL_) ) {

	// compute MT2 vars for 0 b events passing this presel
	mt2b_ = calculateMT2w(jets_, jets_csv_, stopt.lep1(), met_, metphi_, MT2b);
	mt2bl_ = calculateMT2w(jets_, jets_csv_, stopt.lep1(), met_, metphi_, MT2bl);
	mt2w_ = calculateMT2w(jets_, jets_csv_, stopt.lep1(), met_, metphi_, MT2w);

        fillHists1DWrapper(h_1d_cr5_presel,evtweight1l,"cr5_presel");

	bool fail = false;
	if ( !fail && (njetsalleta_ == 2) ) {
	  if (doNM1Plots) fillHists1DWrapper(h_1d_cr5_bbmass_nm1,evtweight1l,"cr5_bbmass_nm1");
	}
	else fail = true;

	if ( !fail && ((bb_.M() > CUT_BBMASS_LOW_) && (bb_.M() < CUT_BBMASS_HIGH_)) ) {
	  if (doNM1Plots) fillHists1DWrapper(h_1d_cr5_met_nm1,evtweight1l,"cr5_met_nm1");
	}
	else fail = true;

	if (!fail && (met_ > 100.) ) {
	  if (doNM1Plots) fillHists1DWrapper(h_1d_cr5_met100,evtweight1l,"cr5_met100");
	}
	else fail = true;

	if (!fail && (met_ > CUT_MET_) ) {
	  if (doNM1Plots) fillHists1DWrapper(h_1d_cr5_mt_nm1,evtweight1l,"cr5_mt_nm1");
	}
	else fail = true;

	if (!fail && (mt_ > CUT_MT_) ) {
	  if (doNM1Plots) fillHists1DWrapper(h_1d_cr5_mt2bl_nm1,evtweight1l,"cr5_mt2bl_nm1");
	}
	else fail = true;

	if (!fail && (mt2bl_ > CUT_MT2BL_) ) {
	  fillHists1DWrapper(h_1d_cr5_final,evtweight1l,"cr5_final");
	}

      } // CR5 region sel

      // -------------------------------------------
      // *** CR5, met last: require exactly 0 btags
      //  otherwise same as signal

      if ( doCR5METLast
	   && passSingleLeptonSelection(isData) 
	   && passisotrk 
	   && (nbjets_ == 0)
	   //	   && (bb_.M() > CUT_BBMASS_LOW_) && (bb_.M() < CUT_BBMASS_HIGH_)
	   && (met_ > CUT_MET_PRESEL_) 
	   && (mt_ > CUT_MT_PRESEL_) ) {

	// compute MT2 vars for 0 b events passing this presel
	mt2b_ = calculateMT2w(jets_, jets_csv_, stopt.lep1(), met_, metphi_, MT2b);
	mt2bl_ = calculateMT2w(jets_, jets_csv_, stopt.lep1(), met_, metphi_, MT2bl);
	mt2w_ = calculateMT2w(jets_, jets_csv_, stopt.lep1(), met_, metphi_, MT2w);

        fillHists1DWrapper(h_1d_cr5_metlast_presel,evtweight1l,"cr5_metlast_presel");

	bool fail = false;
	if ( !fail && (njetsalleta_ == 2) ) {
	  if (doNM1Plots) fillHists1DWrapper(h_1d_cr5_metlast_bbmass_nm1,evtweight1l,"cr5_metlast_bbmass_nm1");
	}
	else fail = true;

	if ( !fail && ((bb_.M() > CUT_BBMASS_LOW_) && (bb_.M() < CUT_BBMASS_HIGH_)) ) {
	  if (doNM1Plots) fillHists1DWrapper(h_1d_cr5_metlast_mt2bl_nm1,evtweight1l,"cr5_metlast_mt2bl_nm1");
	}
	else fail = true;

	if (!fail && (mt_ > CUT_MT_) ) {
	  if (doNM1Plots) fillHists1DWrapper(h_1d_cr5_metlast_mtfirst,evtweight1l,"cr5_metlast_mtfirst");
	}

	if (!fail && (mt2bl_ > CUT_MT2BL_) ) {
	  if (doNM1Plots) fillHists1DWrapper(h_1d_cr5_metlast_mt_nm1,evtweight1l,"cr5_metlast_mt_nm1");
	}
	else fail = true;

	if (!fail && (met_ > 100.) ) {
	  if (doNM1Plots) fillHists1DWrapper(h_1d_cr5_metlast_nomt_met100,evtweight1l,"cr5_metlast_nomt_met100");
	}

	if (!fail && (met_ > 150.) ) {
	  if (doNM1Plots) fillHists1DWrapper(h_1d_cr5_metlast_nomt_met150,evtweight1l,"cr5_metlast_nomt_met150");
	}

	if (!fail && (met_ > 175.) ) {
	  if (doNM1Plots) fillHists1DWrapper(h_1d_cr5_metlast_nomt_met175,evtweight1l,"cr5_metlast_nomt_met175");
	}

	if (!fail && (mt_ > CUT_MT_) ) {
	  if (doNM1Plots) fillHists1DWrapper(h_1d_cr5_metlast_met_nm1,evtweight1l,"cr5_metlast_met_nm1");
	}
	else fail = true;

	//	if (isData && !fail && (met_ > 500.)) dumpEventInfo("cr5 high met event");

	if (!fail && (met_ > 100.) ) {
	  if (doNM1Plots) fillHists1DWrapper(h_1d_cr5_metlast_met100,evtweight1l,"cr5_metlast_met100");
	}
	else fail = true;

	if (!fail && (met_ > 150.) ) {
	  if (doNM1Plots) fillHists1DWrapper(h_1d_cr5_metlast_met150,evtweight1l,"cr5_metlast_met150");
	}
	else fail = true;

	if (!fail && (met_ > CUT_MET_) ) {
	  fillHists1DWrapper(h_1d_cr5_metlast_final,evtweight1l,"cr5_metlast_final");
	}

      } // CR5 region sel

      // -------------------------------------------
      // *** CR5, inverted mass: require exactly 0 btags
      //  otherwise same as signal

      if ( doCR5InvMass
	   && passSingleLeptonSelection(isData) 
	   && passisotrk 
	   && (nbjets_ == 0)
	   //	   && (bb_.M() > CUT_BBMASS_LOW_) && (bb_.M() < CUT_BBMASS_HIGH_)
	   && (met_ > CUT_MET_PRESEL_) 
	   && (mt_ > CUT_MT_PRESEL_) ) {

	// compute MT2 vars for 0 b events passing this presel
	mt2b_ = calculateMT2w(jets_, jets_csv_, stopt.lep1(), met_, metphi_, MT2b);
	mt2bl_ = calculateMT2w(jets_, jets_csv_, stopt.lep1(), met_, metphi_, MT2bl);
	mt2w_ = calculateMT2w(jets_, jets_csv_, stopt.lep1(), met_, metphi_, MT2w);

        fillHists1DWrapper(h_1d_cr5_invmass_presel,evtweight1l,"cr5_invmass_presel");

	bool fail = false;
	if ( !fail && (njetsalleta_ == 2) ) {
	  if (doNM1Plots) fillHists1DWrapper(h_1d_cr5_invmass_bbmass_nm1,evtweight1l,"cr5_invmass_bbmass_nm1");
	}
	else fail = true;

	if ( !fail && ((bb_.M() < CUT_BBMASS_LOW_) || (bb_.M() > CUT_BBMASS_HIGH_)) ) {
	  if (doNM1Plots) fillHists1DWrapper(h_1d_cr5_invmass_mt2bl_nm1,evtweight1l,"cr5_invmass_mt2bl_nm1");
	}
	else fail = true;

	if (!fail && (mt_ > CUT_MT_) ) {
	  if (doNM1Plots) fillHists1DWrapper(h_1d_cr5_invmass_mtfirst,evtweight1l,"cr5_invmass_mtfirst");
	}

	if (!fail && (mt2bl_ > CUT_MT2BL_) ) {
	  if (doNM1Plots) fillHists1DWrapper(h_1d_cr5_invmass_mt_nm1,evtweight1l,"cr5_invmass_mt_nm1");
	}
	else fail = true;

	if (!fail && (met_ > 100.) ) {
	  if (doNM1Plots) fillHists1DWrapper(h_1d_cr5_invmass_nomt_met100,evtweight1l,"cr5_invmass_nomt_met100");
	}

	if (!fail && (met_ > 150.) ) {
	  if (doNM1Plots) fillHists1DWrapper(h_1d_cr5_invmass_nomt_met150,evtweight1l,"cr5_invmass_nomt_met150");
	}

	if (!fail && (met_ > 175.) ) {
	  if (doNM1Plots) fillHists1DWrapper(h_1d_cr5_invmass_nomt_met175,evtweight1l,"cr5_invmass_nomt_met175");
	}

	if (!fail && (mt_ > CUT_MT_) ) {
	  if (doNM1Plots) fillHists1DWrapper(h_1d_cr5_invmass_met_nm1,evtweight1l,"cr5_invmass_met_nm1");
	}
	else fail = true;

	//	if (isData && !fail && (met_ > 500.)) dumpEventInfo("cr5 high met event");

	if (!fail && (met_ > 100.) ) {
	  if (doNM1Plots) fillHists1DWrapper(h_1d_cr5_invmass_met100,evtweight1l,"cr5_invmass_met100");
	}
	else fail = true;

	if (!fail && (met_ > 150.) ) {
	  if (doNM1Plots) fillHists1DWrapper(h_1d_cr5_invmass_met150,evtweight1l,"cr5_invmass_met150");
	}
	else fail = true;

	if (!fail && (met_ > CUT_MET_) ) {
	  fillHists1DWrapper(h_1d_cr5_invmass_final,evtweight1l,"cr5_invmass_final");
	}

      } // CR5 region sel

      // -------------------------------------------
      // *** CR6: require exactly 1 btag
      //  veto on 2nd loose btag
      //  otherwise same as signal

      if ( doCR6
	   && passSingleLeptonSelection(isData) 
	   && passisotrk 
	   && (nbjets_ == 1)
	   && (nbjetsl_ == 1)
	   && (bb_.M() > CUT_BBMASS_LOW_) && (bb_.M() < CUT_BBMASS_HIGH_)
	   && (met_ > CUT_MET_PRESEL_) ) {

	// compute MT2 vars for 1 b events passing this presel
	mt2b_ = calculateMT2w(jets_, jets_csv_, stopt.lep1(), met_, metphi_, MT2b);
	mt2bl_ = calculateMT2w(jets_, jets_csv_, stopt.lep1(), met_, metphi_, MT2bl);
	mt2w_ = calculateMT2w(jets_, jets_csv_, stopt.lep1(), met_, metphi_, MT2w);

        fillHists1DWrapper(h_1d_cr6_presel,evtweight1l,"cr6_presel");

	bool fail = false;
	if ( !fail && (njetsalleta_ == 2) ) {
	  if (doNM1Plots) fillHists1DWrapper(h_1d_cr6_met_nm1,evtweight1l,"cr6_met_nm1");
	}
	else fail = true;

	if (!fail && (met_ > CUT_MET_) ) {
	  if (doNM1Plots) fillHists1DWrapper(h_1d_cr6_mt_nm1,evtweight1l,"cr6_mt_nm1");
	}
	else fail = true;

	if (!fail && (mt_ > CUT_MT_) ) {
	  if (doNM1Plots) fillHists1DWrapper(h_1d_cr6_mt2bl_nm1,evtweight1l,"cr6_mt2bl_nm1");
	}
	else fail = true;

	if (!fail && (mt2bl_ > CUT_MT2BL_) ) {
	  fillHists1DWrapper(h_1d_cr6_final,evtweight1l,"cr6_final");
	}

      } // CR6 region sel


      // -------------------------------------------
      // *** CR6, met last: require exactly 1 btag
      //  veto on 2nd loose btag
      //  otherwise same as signal

      if ( doCR6METLast
	   && passSingleLeptonSelection(isData) 
	   && passisotrk 
	   && (nbjets_ == 1)
	   && (nbjetsl_ == 1)
	   //	   && (bb_.M() > CUT_BBMASS_LOW_) && (bb_.M() < CUT_BBMASS_HIGH_)
	   && (met_ > CUT_MET_PRESEL_) 
	   && (mt_ > CUT_MT_PRESEL_) ) {

        fillHists1DWrapper(h_1d_cr6_metlast_presel,evtweight1l,"cr6_metlast_presel");

	bool fail = false;
	if ( !fail && (njetsalleta_ == 2) ) {
	  if (doNM1Plots) fillHists1DWrapper(h_1d_cr6_metlast_bbmass_nm1,evtweight1l,"cr6_metlast_bbmass_nm1");
	}
	else fail = true;

	if ( !fail && ((bb_.M() > CUT_BBMASS_LOW_) && (bb_.M() < CUT_BBMASS_HIGH_)) ) {
	  if (doNM1Plots) fillHists1DWrapper(h_1d_cr6_metlast_mt2bl_nm1,evtweight1l,"cr6_metlast_mt2bl_nm1");
	}
	else fail = true;

	if (!fail && (mt2bl_ > CUT_MT2BL_) ) {
	  if (doNM1Plots) fillHists1DWrapper(h_1d_cr6_metlast_mt_nm1,evtweight1l,"cr6_metlast_mt_nm1");
	}
	else fail = true;

	if (!fail && (mt_ > CUT_MT_) ) {
	  if (doNM1Plots) fillHists1DWrapper(h_1d_cr6_metlast_met_nm1,evtweight1l,"cr6_metlast_met_nm1");
	}
	else fail = true;

	if (!fail && (met_ > 100.) ) {
	  if (doNM1Plots) fillHists1DWrapper(h_1d_cr6_metlast_met100,evtweight1l,"cr6_metlast_met100");
	}
	else fail = true;

	if (!fail && (met_ > 150.) ) {
	  if (doNM1Plots) fillHists1DWrapper(h_1d_cr6_metlast_met150,evtweight1l,"cr6_metlast_met150");
	}
	else fail = true;

	if (!fail && (met_ > CUT_MET_) ) {
	  fillHists1DWrapper(h_1d_cr6_metlast_final,evtweight1l,"cr6_metlast_final");
	}

      } // CR6 region sel


      // -------------------------------------------
      // *** CR7: 150 < m(bb) < 250, exactly 3 jets
      //  otherwise same as signal region

      if ( doCR7
	   && passSingleLeptonSelection(isData) 
	   && passisotrk 
	   && (nbjets_ >= 2)
	   && (bb_.M() > CUT_BBMASS_CR1_LOW_) && (bb_.M() < CUT_BBMASS_CR1_HIGH_) 
	   && (met_ > CUT_MET_PRESEL_) 
	   && (mt_ > CUT_MT_PRESEL_) ) {

        fillHists1DWrapper(h_1d_cr7_presel,evtweight1l,"cr7_presel");

	bool fail = false;
	//	if ( !fail && (njetsalleta_ == 3) && (nbjets_ == 2) ) {
	if ( !fail && (njetsalleta_ == 3) ) {
	  if (doNM1Plots) fillHists1DWrapper(h_1d_cr7_met_nm1,evtweight1l,"cr7_met_nm1");
	}
	else fail = true;

	if (!fail && (met_ > CUT_MET_) ) {
	  if (doNM1Plots) fillHists1DWrapper(h_1d_cr7_mt_nm1,evtweight1l,"cr7_mt_nm1");
	}
	else fail = true;

	if (!fail && (mt_ > CUT_MT_) ) {
	  if (doNM1Plots) fillHists1DWrapper(h_1d_cr7_mt2bl_nm1,evtweight1l,"cr7_mt2bl_nm1");
	}
	else fail = true;

	if (!fail && (mt2bl_ > CUT_MT2BL_) ) {
	  fillHists1DWrapper(h_1d_cr7_final,evtweight1l,"cr7_final");
	}

      } // CR7 region sel


      // -------------------------------------------
      // *** CR8: 50 < m(bb) < 100 
      //  otherwise same as signal region

      if ( doCR8
	   && passSingleLeptonSelection(isData) 
	   && passisotrk 
	   && (nbjets_ >= 2)
	   && (bb_.M() > CUT_BBMASS_CR8_LOW_) && (bb_.M() < CUT_BBMASS_CR8_HIGH_) 
	   && (met_ > CUT_MET_PRESEL_) 
	   && (mt_ > CUT_MT_PRESEL_) ) {

        fillHists1DWrapper(h_1d_cr8_presel,evtweight1l,"cr8_presel");

	bool fail = false;
	if ( !fail && (njetsalleta_ == 2) ) {
	  if (doNM1Plots) fillHists1DWrapper(h_1d_cr8_met_nm1,evtweight1l,"cr8_met_nm1");
	}
	else fail = true;

	if (!fail && (met_ > CUT_MET_) ) {
	  if (doNM1Plots) fillHists1DWrapper(h_1d_cr8_mt_nm1,evtweight1l,"cr8_mt_nm1");
	}
	else fail = true;

	if (!fail && (mt_ > CUT_MT_) ) {
	  if (doNM1Plots) fillHists1DWrapper(h_1d_cr8_mt2bl_nm1,evtweight1l,"cr8_mt2bl_nm1");
	}
	else fail = true;

	if (!fail && (mt2bl_ > CUT_MT2BL_) ) {
	  fillHists1DWrapper(h_1d_cr8_final,evtweight1l,"cr8_final");
	}

      } // CR8 region sel


      // -------------------------------------------
      // *** CR8, MET last: 50 < m(bb) < 100 
      //  otherwise same as signal region

      if ( doCR8METLast
	   && passSingleLeptonSelection(isData) 
	   && passisotrk 
	   && (nbjets_ >= 2)
	   && (bb_.M() > CUT_BBMASS_CR8_LOW_) && (bb_.M() < CUT_BBMASS_CR8_HIGH_) 
	   && (met_ > CUT_MET_PRESEL_) 
	   && (mt_ > CUT_MT_PRESEL_) ) {

        fillHists1DWrapper(h_1d_cr8_metlast_presel,evtweight1l,"cr8_metlast_presel");

	bool fail = false;
	if ( !fail && (njetsalleta_ == 2) ) {
	  if (doNM1Plots) fillHists1DWrapper(h_1d_cr8_metlast_mt2bl_nm1,evtweight1l,"cr8_metlast_mt2bl_nm1");
	}
	else fail = true;

	if (!fail && (mt_ > CUT_MT_) ) {
	  fillHists1DWrapper(h_1d_cr8_metlast_mtfirst,evtweight1l,"cr8_metlast_mtfirst");
	}

	if (!fail && (mt2bl_ > CUT_MT2BL_) ) {
	  if (doNM1Plots) fillHists1DWrapper(h_1d_cr8_metlast_mt_nm1,evtweight1l,"cr8_metlast_mt_nm1");
	}
	else fail = true;

	if (!fail && (mt_ > CUT_MT_) ) {
	  fillHists1DWrapper(h_1d_cr8_metlast_met_nm1,evtweight1l,"cr8_metlast_met_nm1");
	}
	else fail = true;

	if (!fail && (met_ > 100.) ) {
	  if (doNM1Plots) fillHists1DWrapper(h_1d_cr8_metlast_met100,evtweight1l,"cr8_metlast_met100");
	}
	else fail = true;

	if (!fail && (met_ > 150.) ) {
	  if (doNM1Plots) fillHists1DWrapper(h_1d_cr8_metlast_met150,evtweight1l,"cr8_metlast_met150");
	}
	else fail = true;

	if (!fail && (met_ > CUT_MET_) ) {
	  if (doNM1Plots) fillHists1DWrapper(h_1d_cr8_metlast_final,evtweight1l,"cr8_metlast_final");
	}
	else fail = true;

      } // CR8 region sel


      // -------------------------------------------
      // *** CR idea: signal mass window, at least 3 jets, then exactly 3 jets
      //  otherwise same as signal region
      //  !! need to check to make sure signal contribution isn't too large..
      //  -- after final selection, bg ~3, sig ~2 => don't use this region!


      // -------------------------------------------
      // *** CR9: 150 < m(bb) < 200 
      //  otherwise same as signal region
      //  to study met shape disagreement in CR1

      if ( doCR9
	   && passSingleLeptonSelection(isData) 
	   && passisotrk 
	   && (nbjets_ >= 2)
	   && (bb_.M() > CUT_BBMASS_CR1_LOW_) && (bb_.M() < 200.) 
	   && (met_ > CUT_MET_PRESEL_) 
	   && (mt_ > CUT_MT_PRESEL_) ) {

        fillHists1DWrapper(h_1d_cr9_presel,evtweight1l,"cr9_presel");

	bool fail = false;
	if ( !fail && (njetsalleta_ == 2) ) {
	  if (doNM1Plots) fillHists1DWrapper(h_1d_cr9_met_nm1,evtweight1l,"cr9_met_nm1");
	}
	else fail = true;

	// mt peak before cutting on met
	if (!fail && (mt_ > 50.) && (mt_ < 80.) ) {
	  if (doNM1Plots) fillHists1DWrapper(h_1d_cr9_mtpeak_nomet,evtweight1l,"cr9_mtpeak_nomet");
	}

	if (!fail && (met_ > CUT_MET_) ) {
	  if (doNM1Plots) fillHists1DWrapper(h_1d_cr9_mt_nm1,evtweight1l,"cr9_mt_nm1");
	}
	else fail = true;

	// mt peak after cutting on met
	if (!fail && (mt_ > 50.) && (mt_ < 80.) ) {
	  if (doNM1Plots) fillHists1DWrapper(h_1d_cr9_mtpeak_met,evtweight1l,"cr9_mtpeak_met");
	}

	if (!fail && (mt2bl_ > CUT_MT2BL_) ) {
	  if (doNM1Plots) fillHists1DWrapper(h_1d_cr9_mtlast_mt2bl_nm1,evtweight1l,"cr9_mtlast_mt2bl_nm1");
	}
	// else fail = true;

	if (!fail && (mt_ > CUT_MT_) ) {
	  if (doNM1Plots) fillHists1DWrapper(h_1d_cr9_mt2bl_nm1,evtweight1l,"cr9_mt2bl_nm1");
	}
	else fail = true;

	if (!fail && (mt2bl_ > CUT_MT2BL_) ) {
	  fillHists1DWrapper(h_1d_cr9_final,evtweight1l,"cr9_final");
	}

      } // CR9 region sel


      // -------------------------------------------
      // *** CR10: 200 < m(bb) < 250 
      //  otherwise same as signal region
      //  to study met shape disagreement in CR1

      if ( doCR10
	   && passSingleLeptonSelection(isData) 
	   && passisotrk 
	   && (nbjets_ >= 2)
	   && (bb_.M() > 200.) && (bb_.M() < 250.) 
	   && (met_ > CUT_MET_PRESEL_) 
	   && (mt_ > CUT_MT_PRESEL_) ) {

        fillHists1DWrapper(h_1d_cr10_presel,evtweight1l,"cr10_presel");

	bool fail = false;
	if ( !fail && (njetsalleta_ == 2) ) {
	  if (doNM1Plots) fillHists1DWrapper(h_1d_cr10_met_nm1,evtweight1l,"cr10_met_nm1");
	}
	else fail = true;

	if (!fail && (met_ > CUT_MET_) ) {
	  if (doNM1Plots) fillHists1DWrapper(h_1d_cr10_mt_nm1,evtweight1l,"cr10_mt_nm1");
	}
	else fail = true;

	if (!fail && (mt_ > CUT_MT_) ) {
	  if (doNM1Plots) fillHists1DWrapper(h_1d_cr10_mt2bl_nm1,evtweight1l,"cr10_mt2bl_nm1");
	}
	else fail = true;

	if (!fail && (mt2bl_ > CUT_MT2BL_) ) {
	  fillHists1DWrapper(h_1d_cr10_final,evtweight1l,"cr10_final");
	}

      } // CR10 region sel


      // -------------------------------------------
      // *** CR11: require exactly 0 btags, 150 < m(bb) < 250
      //  otherwise same as signal

      if ( doCR11
	   && passSingleLeptonSelection(isData) 
	   && passisotrk 
	   && (nbjets_ == 0)
	   && (bb_.M() > CUT_BBMASS_CR1_LOW_) && (bb_.M() < CUT_BBMASS_CR1_HIGH_)
	   && (met_ > CUT_MET_PRESEL_) 
	   && (mt_ > CUT_MT_PRESEL_) ) {

	// compute MT2 vars for 0 b events passing this presel
	mt2b_ = calculateMT2w(jets_, jets_csv_, stopt.lep1(), met_, metphi_, MT2b);
	mt2bl_ = calculateMT2w(jets_, jets_csv_, stopt.lep1(), met_, metphi_, MT2bl);
	mt2w_ = calculateMT2w(jets_, jets_csv_, stopt.lep1(), met_, metphi_, MT2w);

        fillHists1DWrapper(h_1d_cr11_presel,evtweight1l,"cr11_presel");

	bool fail = false;
	if ( !fail && (njetsalleta_ == 2) ) {
	  if (doNM1Plots) fillHists1DWrapper(h_1d_cr11_met_nm1,evtweight1l,"cr11_met_nm1");
	}
	else fail = true;

	if (!fail && (met_ > CUT_MET_) ) {
	  if (doNM1Plots) fillHists1DWrapper(h_1d_cr11_mt_nm1,evtweight1l,"cr11_mt_nm1");
	}
	else fail = true;

	if (!fail && (mt_ > CUT_MT_) ) {
	  if (doNM1Plots) fillHists1DWrapper(h_1d_cr11_mt2bl_nm1,evtweight1l,"cr11_mt2bl_nm1");
	}
	else fail = true;

	if (!fail && (mt2bl_ > CUT_MT2BL_) ) {
	  fillHists1DWrapper(h_1d_cr11_final,evtweight1l,"cr11_final");
	}

      } // CR11 region sel

      // -------------------------------------------
      // *** CR12: 150 < m(bb) < 250, at least 4 jets
      //  otherwise same as signal region

      if ( doCR12
	   && passSingleLeptonSelection(isData) 
	   && passisotrk 
	   && (njetsalleta_ >= 4)
	   && (nbjets_ >= 2)
	   //	   && (bb_.M() > CUT_BBMASS_CR1_LOW_) && (bb_.M() < CUT_BBMASS_CR1_HIGH_) 
	   && (met_ > CUT_MET_PRESEL_) 
	   && (mt_ > CUT_MT_PRESEL_) ) {

        fillHists1DWrapper(h_1d_cr12_presel,evtweight1l,"cr12_presel");

	bool fail = false;
	if ( !fail && (bb_.M() > CUT_BBMASS_CR1_LOW_) && (bb_.M() < CUT_BBMASS_CR1_HIGH_) ) {
	  if (doNM1Plots) fillHists1DWrapper(h_1d_cr12_met_nm1,evtweight1l,"cr12_met_nm1");
	}
	else fail = true;

	if (!fail && (met_ > CUT_MET_) ) {
	  if (doNM1Plots) fillHists1DWrapper(h_1d_cr12_mt_nm1,evtweight1l,"cr12_mt_nm1");
	}
	else fail = true;

	if (!fail && (mt_ > CUT_MT_) ) {
	  if (doNM1Plots) fillHists1DWrapper(h_1d_cr12_mt2bl_nm1,evtweight1l,"cr12_mt2bl_nm1");
	}
	else fail = true;

	if (!fail && (mt2bl_ > CUT_MT2BL_) ) {
	  fillHists1DWrapper(h_1d_cr12_final,evtweight1l,"cr12_final");
	}

      } // CR12 region sel


      // -------------------------------------------
      // *** CR13: MT peak region
      //  - presel:
      //   == 1 lepton, iso track veto
      //   >= 2 bjets
      //   100 < m(bb) < 150
      //   met > 50

      if ( doCR13
	   && passSingleLeptonSelection(isData) 
	   && passisotrk 
	   && (nbjets_ >= 2)
	   //	   && (bb_.M() > CUT_BBMASS_LOW_) && (bb_.M() < CUT_BBMASS_HIGH_)
	   && (met_ > CUT_MET_PRESEL_) 
	   && (mt_ > CUT_MT_CR13_LOW_) && (mt_ <= CUT_MT_CR13_HIGH_) ) {

        fillHists1DWrapper(h_1d_cr13_presel,evtweight1l,"cr13_presel");

	bool fail = false;
	if ( !fail && (njetsalleta_ == 2) ) {
	  if (doNM1Plots) fillHists1DWrapper(h_1d_cr13_bbmass_nm1,evtweight1l,"cr13_bbmass_nm1");
	}
	else fail = true;

	if ( !fail && (bb_.M() > CUT_BBMASS_LOW_) && (bb_.M() < CUT_BBMASS_HIGH_) ) {
	  if (doNM1Plots) fillHists1DWrapper(h_1d_cr13_mt2bl_nm1,evtweight1l,"cr13_mt2bl_nm1");
	}
	else fail = true;

	if (!fail && (mt2bl_ > CUT_MT2BL_) ) {
	  if (doNM1Plots) fillHists1DWrapper(h_1d_cr13_met_nm1,evtweight1l,"cr13_met_nm1");
	}
	else fail = true;

	if (!fail && (met_ > 100.) ) {
	  if (doNM1Plots) fillHists1DWrapper(h_1d_cr13_met100,evtweight1l,"cr13_met100");
	}
	else fail = true;

	if (!fail && (met_ > 150.) ) {
	  if (doNM1Plots) fillHists1DWrapper(h_1d_cr13_met150,evtweight1l,"cr13_met150");
	}
	else fail = true;

	if (!fail && (met_ > CUT_MET_) ) {
	  fillHists1DWrapper(h_1d_cr13_final,evtweight1l,"cr13_final");
	}
	else fail = true;

      } // CR13 region sel


      // -------------------------------------------
      // *** CR13: MT peak region, inverted mbb
      //  - presel:
      //   == 1 lepton, iso track veto
      //   >= 2 bjets
      //   100 < m(bb) < 150
      //   met > 50

      if ( doCR13
	   && passSingleLeptonSelection(isData) 
	   && passisotrk 
	   && (nbjets_ >= 2)
	   //	   && (bb_.M() > CUT_BBMASS_LOW_) && (bb_.M() < CUT_BBMASS_HIGH_)
	   && (met_ > CUT_MET_PRESEL_) 
	   && (mt_ > CUT_MT_CR13_LOW_) && (mt_ <= CUT_MT_CR13_HIGH_) ) {

        fillHists1DWrapper(h_1d_cr13_presel,evtweight1l,"cr13_presel");

	bool fail = false;
	if ( !fail && (njetsalleta_ == 2) ) {
	  if (doNM1Plots) fillHists1DWrapper(h_1d_cr13_bbmass_nm1,evtweight1l,"cr13_bbmass_nm1");
	}
	else fail = true;

	if ( !fail && (bb_.M() > CUT_BBMASS_LOW_) && (bb_.M() < CUT_BBMASS_HIGH_) ) {
	  if (doNM1Plots) fillHists1DWrapper(h_1d_cr13_mt2bl_nm1,evtweight1l,"cr13_mt2bl_nm1");
	}
	else fail = true;

	if (!fail && (mt2bl_ > CUT_MT2BL_) ) {
	  if (doNM1Plots) fillHists1DWrapper(h_1d_cr13_met_nm1,evtweight1l,"cr13_met_nm1");
	}
	else fail = true;

	if (!fail && (met_ > 100.) ) {
	  if (doNM1Plots) fillHists1DWrapper(h_1d_cr13_met100,evtweight1l,"cr13_met100");
	}
	else fail = true;

	if (!fail && (met_ > 150.) ) {
	  if (doNM1Plots) fillHists1DWrapper(h_1d_cr13_met150,evtweight1l,"cr13_met150");
	}
	else fail = true;

	if (!fail && (met_ > CUT_MET_) ) {
	  fillHists1DWrapper(h_1d_cr13_final,evtweight1l,"cr13_final");
	}
	else fail = true;

      } // CR13 region sel


      // -------------------------------------------
      // *** CR14: inverted m(bb) cut
      //  otherwise same as signal region

      if ( doCR14
	   && passSingleLeptonSelection(isData) 
	   && passisotrk 
	   && (nbjets_ >= 2)
	   && ((bb_.M() < CUT_BBMASS_LOW_) || (bb_.M() > CUT_BBMASS_HIGH_)) 
	   && (met_ > CUT_MET_PRESEL_) 
	   && (mt_ > CUT_MT_PRESEL_) ) {

        fillHists1DWrapper(h_1d_cr14_presel,evtweight1l,"cr14_presel");

	bool fail = false;
	if ( !fail && (njetsalleta_ == 2) ) {
	  if (doNM1Plots) fillHists1DWrapper(h_1d_cr14_mt2bl_nm1,evtweight1l,"cr14_mt2bl_nm1");
	}
	else fail = true;

	if (!fail && (mt_ > CUT_MT_) ) {
	  fillHists1DWrapper(h_1d_cr14_mtfirst,evtweight1l,"cr14_mtfirst");
	}

	if (!fail && (mt2bl_ > CUT_MT2BL_) ) {
	  if (doNM1Plots) fillHists1DWrapper(h_1d_cr14_mt_nm1,evtweight1l,"cr14_mt_nm1");
	}
	else fail = true;

	if (!fail && (met_ > 100.) ) {
	  if (doNM1Plots) fillHists1DWrapper(h_1d_cr14_nomt_met100,evtweight1l,"cr14_nomt_met100");
	}

	if (!fail && (met_ > 150.) ) {
	  if (doNM1Plots) fillHists1DWrapper(h_1d_cr14_nomt_met150,evtweight1l,"cr14_nomt_met150");
	}

	if (!fail && (met_ > 175.) ) {
	  if (doNM1Plots) fillHists1DWrapper(h_1d_cr14_nomt_met175,evtweight1l,"cr14_nomt_met175");
	}

	if (!fail && (mt_ > CUT_MT_) ) {
	  fillHists1DWrapper(h_1d_cr14_met_nm1,evtweight1l,"cr14_met_nm1");
	}
	else fail = true;

	if (!fail && (met_ > 100.) ) {
	  if (doNM1Plots) fillHists1DWrapper(h_1d_cr14_met100,evtweight1l,"cr14_met100");
	}
	else fail = true;

	if (!fail && (met_ > 150.) ) {
	  if (doNM1Plots) fillHists1DWrapper(h_1d_cr14_met150,evtweight1l,"cr14_met150");
	}
	else fail = true;

	if (!fail && (met_ > CUT_MET_) ) {
	  if (doNM1Plots) fillHists1DWrapper(h_1d_cr14_final,evtweight1l,"cr14_final");
	}
	else fail = true;

      } // CR14 region sel


      // -------------------------------------------
      // *** Stop Presel region
      //   >= 1 lepton
      //   >= 4 jets
      //   blind data (includes stop signal region)

      // *** Stop Presel comparison region
      //   >= 1 lepton
      //   iso track veto v4
      //   >= 4 jets
      //   >= 1 bjet
      //   met > 150
      //   mt > 120
      //   tau veto
      //   blind data (includes stop signal region)

      if ( doStopSel
	   && passSingleLeptonSelection(isData) 
	   && (njets_ >= 4)
	   && !isData ) {

	float weight = evtweight_novtxweight;
	//	float weight = evtweight1l;

        fillHists1DWrapper(h_1d_stop_presel,weight,"stop_presel");

	bool fail = false;
	if ( !fail && (nbjets_ >= 1) ) {
	  fillHists1DWrapper(h_1d_stop_met_nm1,weight,"stop_met_nm1");
	} 
	else fail = true;

	if ( !fail && (met_ >= 150.) ) {
	  fillHists1DWrapper(h_1d_stop_mt_nm1,weight,"stop_mt_nm1");
	} 
	else fail = true;

	if ( !fail && (mt_ > 120.) ) {
	  fillHists1DWrapper(h_1d_stop_isotrk_nm1,weight,"stop_isotrk_nm1");
	} 
	else fail = true;

	if ( !fail && passIsoTrkVeto_v4() ) {
	  fillHists1DWrapper(h_1d_stop_tauveto_nm1,weight,"stop_tauveto_nm1");
	} 
	else fail = true;

	if ( !fail && passTauVeto() ) {
	  fillHists1DWrapper(h_1d_stop_comp,weight,"stop_comp");
	} 
	else fail = true;

      } // stop region sel

      // end regions ----------------------------------------
      // ----------------------------------------------------

      ++nEventsPass;

    } // end event loop

    // delete tree;
    
  } // end file loop
  
    //
    // finish
    //

  if (doInclusive) {
    savePlotsDir(h_1d_inc_presel,outfile_,"inc_presel");
    savePlotsDir(h_1d_inc_2j,outfile_,"inc_2j");
    if (doInclusiveMTTail) {
      savePlotsDir(h_1d_inc_2j_mt,outfile_,"inc_2j_mt");
      savePlotsDir(h_1d_inc_2j_mt_met100,outfile_,"inc_2j_mt_met100");
      savePlotsDir(h_1d_inc_2j_mt_met150,outfile_,"inc_2j_mt_met150");
      savePlotsDir(h_1d_inc_2j_mt_metcut,outfile_,"inc_2j_mt_metcut");
    }
    savePlotsDir(h_1d_inc_1b,outfile_,"inc_1b");
    savePlotsDir(h_1d_inc_2b,outfile_,"inc_2b");
  }

  if (doSignal) {
    savePlotsDir(h_1d_sig_presel,outfile_,"sig_presel");
    if (doNM1Plots) {
      savePlotsDir(h_1d_sig_met_nm1,outfile_,"sig_met_nm1");
      savePlotsDir(h_1d_sig_mtpeak_nomet,outfile_,"sig_mtpeak_nomet");
      savePlotsDir(h_1d_sig_mt_nm1,outfile_,"sig_mt_nm1");
      savePlotsDir(h_1d_sig_mtpeak_met,outfile_,"sig_mtpeak_met");
      savePlotsDir(h_1d_sig_mt2bl_nm1,outfile_,"sig_mt2bl_nm1");
    }
    savePlotsDir(h_1d_sig_final,outfile_,"sig_final");
  }

  if (doSignalMassLast) {
    savePlotsDir(h_1d_sig_bbmasslast_presel,outfile_,"sig_bbmasslast_presel");
    if (doNM1Plots) {
      savePlotsDir(h_1d_sig_bbmasslast_mt2bl_nm1,outfile_,"sig_bbmasslast_mt2bl_nm1");
      savePlotsDir(h_1d_sig_bbmasslast_mtpeak,outfile_,"sig_bbmasslast_mtpeak");
      savePlotsDir(h_1d_sig_bbmasslast_mtcut,outfile_,"sig_bbmasslast_mtcut");
      savePlotsDir(h_1d_sig_bbmasslast_mtcut_met100,outfile_,"sig_bbmasslast_mtcut_met100");
      savePlotsDir(h_1d_sig_bbmasslast_metcut,outfile_,"sig_bbmasslast_metcut");
      savePlotsDir(h_1d_sig_bbmasslast_mt_nm1,outfile_,"sig_bbmasslast_mt_nm1");
      savePlotsDir(h_1d_sig_bbmasslast_met_nm1,outfile_,"sig_bbmasslast_met_nm1");
      savePlotsDir(h_1d_sig_bbmasslast_met100,outfile_,"sig_bbmasslast_met100");
      savePlotsDir(h_1d_sig_bbmasslast_met150,outfile_,"sig_bbmasslast_met150");
      savePlotsDir(h_1d_sig_bbmasslast_bbmass_nm1,outfile_,"sig_bbmasslast_bbmass_nm1");
    }
    savePlotsDir(h_1d_sig_bbmasslast_final,outfile_,"sig_bbmasslast_final");
  }

  if (doSignalMETLast) {
    savePlotsDir(h_1d_sig_metlast_presel,outfile_,"sig_metlast_presel");
    if (doNM1Plots) {
      savePlotsDir(h_1d_sig_metlast_mt2bl_nm1,outfile_,"sig_metlast_mt2bl_nm1");
      savePlots2Dir(h_2d_sig_metlast_mt2bl_nm1,outfile_,"sig_metlast_mt2bl_nm1");
      savePlotsDir(h_1d_sig_metlast_mtfirst,outfile_,"sig_metlast_mtfirst");
      savePlotsDir(h_1d_sig_metlast_mt_nm1,outfile_,"sig_metlast_mt_nm1");
      savePlotsDir(h_1d_sig_metlast_met_nm1,outfile_,"sig_metlast_met_nm1");
      savePlotsDir(h_1d_sig_metlast_met100,outfile_,"sig_metlast_met100");
      savePlotsDir(h_1d_sig_metlast_met150,outfile_,"sig_metlast_met150");
    }
    savePlotsDir(h_1d_sig_metlast_final,outfile_,"sig_metlast_final");
  }

  if (doSignalSMWH) {
    savePlotsDir(h_1d_sigsmwh_presel,outfile_,"sigsmwh_presel");
    // if (doNM1Plots) {
    //   savePlotsDir(h_1d_sig_met_nm1,outfile_,"sig_met_nm1");
    //   savePlotsDir(h_1d_sig_mtpeak_nomet,outfile_,"sig_mtpeak_nomet");
    //   savePlotsDir(h_1d_sig_mt_nm1,outfile_,"sig_mt_nm1");
    //   savePlotsDir(h_1d_sig_mtpeak_met,outfile_,"sig_mtpeak_met");
    //   savePlotsDir(h_1d_sig_mt2bl_nm1,outfile_,"sig_mt2bl_nm1");
    // }
    savePlotsDir(h_1d_sigsmwh_final,outfile_,"sigsmwh_final");
  }

  if (doCR1) {
    savePlotsDir(h_1d_cr1_presel,outfile_,"cr1_presel");
    if (doNM1Plots) {
      savePlotsDir(h_1d_cr1_met_nm1,outfile_,"cr1_met_nm1");
      savePlotsDir(h_1d_cr1_mtpeak_nomet,outfile_,"cr1_mtpeak_nomet");
      savePlotsDir(h_1d_cr1_mt_nm1,outfile_,"cr1_mt_nm1");
      savePlotsDir(h_1d_cr1_mtpeak_met,outfile_,"cr1_mtpeak_met");
      savePlotsDir(h_1d_cr1_mt2bl_nm1,outfile_,"cr1_mt2bl_nm1");
      savePlotsDir(h_1d_cr1_mtlast_mt2bl_nm1,outfile_,"cr1_mtlast_mt2bl_nm1");
    }
    savePlotsDir(h_1d_cr1_final,outfile_,"cr1_final");
  }

  if (doCR1METLast) {
    savePlotsDir(h_1d_cr1_metlast_presel,outfile_,"cr1_metlast_presel");
    if (doNM1Plots) {
      savePlotsDir(h_1d_cr1_metlast_mt2bl_nm1,outfile_,"cr1_metlast_mt2bl_nm1");
      savePlotsDir(h_1d_cr1_metlast_mtfirst,outfile_,"cr1_metlast_mtfirst");
      savePlotsDir(h_1d_cr1_metlast_mt_nm1,outfile_,"cr1_metlast_mt_nm1");
      savePlotsDir(h_1d_cr1_metlast_met_nm1,outfile_,"cr1_metlast_met_nm1");
      savePlotsDir(h_1d_cr1_metlast_met100,outfile_,"cr1_metlast_met100");
      savePlotsDir(h_1d_cr1_metlast_met150,outfile_,"cr1_metlast_met150");
    }
    savePlotsDir(h_1d_cr1_metlast_final,outfile_,"cr1_metlast_final");
  }

  if (doCR2) {
    savePlotsDir(h_1d_cr2_presel,outfile_,"cr2_presel");
    if (doNM1Plots) {
      savePlotsDir(h_1d_cr2_bbmass_nm1,outfile_,"cr2_bbmass_nm1");
      savePlotsDir(h_1d_cr2_mt_nm1,outfile_,"cr2_mt_nm1");
      savePlotsDir(h_1d_cr2_mt2blfirst,outfile_,"cr2_mt2blfirst");
      savePlotsDir(h_1d_cr2_met_nm1,outfile_,"cr2_met_nm1");
      savePlotsDir(h_1d_cr2_mt2blcut,outfile_,"cr2_mt2blcut");
      savePlotsDir(h_1d_cr2_mt2bl_nm1,outfile_,"cr2_mt2bl_nm1");
    }
    savePlotsDir(h_1d_cr2_final,outfile_,"cr2_final");
  }

  if (doCR3) {
    savePlotsDir(h_1d_cr3_presel,outfile_,"cr3_presel");
    if (doNM1Plots) {
      savePlotsDir(h_1d_cr3_bbmass_nm1,outfile_,"cr3_bbmass_nm1");
      savePlotsDir(h_1d_cr3_mt_nm1,outfile_,"cr3_mt_nm1");
      savePlotsDir(h_1d_cr3_mt2blfirst,outfile_,"cr3_mt2blfirst");
      savePlotsDir(h_1d_cr3_met_nm1,outfile_,"cr3_met_nm1");
      savePlotsDir(h_1d_cr3_mt2blcut,outfile_,"cr3_mt2blcut");
      savePlotsDir(h_1d_cr3_mt2bl_nm1,outfile_,"cr3_mt2bl_nm1");
    }
    savePlotsDir(h_1d_cr3_final,outfile_,"cr3_final");
  }

  if (doCR23) {
    savePlotsDir(h_1d_cr23_presel,outfile_,"cr23_presel");
    if (doNM1Plots) {
      savePlotsDir(h_1d_cr23_bbmass_nm1,outfile_,"cr23_bbmass_nm1");
      savePlotsDir(h_1d_cr23_mt_nm1,outfile_,"cr23_mt_nm1");
      savePlotsDir(h_1d_cr23_mt2blfirst,outfile_,"cr23_mt2blfirst");
      savePlotsDir(h_1d_cr23_met_nm1,outfile_,"cr23_met_nm1");
      savePlotsDir(h_1d_cr23_mt2blcut,outfile_,"cr23_mt2blcut");
      savePlotsDir(h_1d_cr23_mt2bl_nm1,outfile_,"cr23_mt2bl_nm1");
    }
    savePlotsDir(h_1d_cr23_final,outfile_,"cr23_final");
  }

  if (doCR4) {
    savePlotsDir(h_1d_cr4_presel,outfile_,"cr4_presel");
    if (doNM1Plots) {
      savePlotsDir(h_1d_cr4_met_nm1,outfile_,"cr4_met_nm1");
      savePlotsDir(h_1d_cr4_mt_nm1,outfile_,"cr4_mt_nm1");
      savePlotsDir(h_1d_cr4_mt2bl_nm1,outfile_,"cr4_mt2bl_nm1");
    }
    savePlotsDir(h_1d_cr4_final,outfile_,"cr4_final");
  }

  if (doCR5) {
    savePlotsDir(h_1d_cr5_presel,outfile_,"cr5_presel");
    if (doNM1Plots) {
      savePlotsDir(h_1d_cr5_bbmass_nm1,outfile_,"cr5_bbmass_nm1");
      savePlotsDir(h_1d_cr5_met_nm1,outfile_,"cr5_met_nm1");
      savePlotsDir(h_1d_cr5_met100,outfile_,"cr5_met100");
      savePlotsDir(h_1d_cr5_mt_nm1,outfile_,"cr5_mt_nm1");
      savePlotsDir(h_1d_cr5_mt2bl_nm1,outfile_,"cr5_mt2bl_nm1");
    }
    savePlotsDir(h_1d_cr5_final,outfile_,"cr5_final");
  }

  if (doCR5METLast) {
    savePlotsDir(h_1d_cr5_metlast_presel,outfile_,"cr5_metlast_presel");
    if (doNM1Plots) {
      savePlotsDir(h_1d_cr5_metlast_mt2bl_nm1,outfile_,"cr5_metlast_mt2bl_nm1");
      savePlotsDir(h_1d_cr5_metlast_mtfirst,outfile_,"cr5_metlast_mtfirst");
      savePlotsDir(h_1d_cr5_metlast_mt_nm1,outfile_,"cr5_metlast_mt_nm1");
      savePlotsDir(h_1d_cr5_metlast_nomt_met100,outfile_,"cr5_metlast_nomt_met100");
      savePlotsDir(h_1d_cr5_metlast_nomt_met150,outfile_,"cr5_metlast_nomt_met150");
      savePlotsDir(h_1d_cr5_metlast_nomt_met175,outfile_,"cr5_metlast_nomt_met175");
      savePlotsDir(h_1d_cr5_metlast_met_nm1,outfile_,"cr5_metlast_met_nm1");
      savePlotsDir(h_1d_cr5_metlast_met100,outfile_,"cr5_metlast_met100");
      savePlotsDir(h_1d_cr5_metlast_met150,outfile_,"cr5_metlast_met150");
    }
    savePlotsDir(h_1d_cr5_metlast_final,outfile_,"cr5_metlast_final");
  }

  if (doCR5InvMass) {
    savePlotsDir(h_1d_cr5_invmass_presel,outfile_,"cr5_invmass_presel");
    if (doNM1Plots) {
      savePlotsDir(h_1d_cr5_invmass_mt2bl_nm1,outfile_,"cr5_invmass_mt2bl_nm1");
      savePlotsDir(h_1d_cr5_invmass_mtfirst,outfile_,"cr5_invmass_mtfirst");
      savePlotsDir(h_1d_cr5_invmass_mt_nm1,outfile_,"cr5_invmass_mt_nm1");
      savePlotsDir(h_1d_cr5_invmass_nomt_met100,outfile_,"cr5_invmass_nomt_met100");
      savePlotsDir(h_1d_cr5_invmass_nomt_met150,outfile_,"cr5_invmass_nomt_met150");
      savePlotsDir(h_1d_cr5_invmass_nomt_met175,outfile_,"cr5_invmass_nomt_met175");
      savePlotsDir(h_1d_cr5_invmass_met_nm1,outfile_,"cr5_invmass_met_nm1");
      savePlotsDir(h_1d_cr5_invmass_met100,outfile_,"cr5_invmass_met100");
      savePlotsDir(h_1d_cr5_invmass_met150,outfile_,"cr5_invmass_met150");
    }
    savePlotsDir(h_1d_cr5_invmass_final,outfile_,"cr5_invmass_final");
  }

  if (doCR6) {
    savePlotsDir(h_1d_cr6_presel,outfile_,"cr6_presel");
    if (doNM1Plots) {
      savePlotsDir(h_1d_cr6_met_nm1,outfile_,"cr6_met_nm1");
      savePlotsDir(h_1d_cr6_mt_nm1,outfile_,"cr6_mt_nm1");
      savePlotsDir(h_1d_cr6_mt2bl_nm1,outfile_,"cr6_mt2bl_nm1");
    }
    savePlotsDir(h_1d_cr6_final,outfile_,"cr6_final");
  }

  if (doCR6METLast) {
    savePlotsDir(h_1d_cr6_metlast_presel,outfile_,"cr6_metlast_presel");
    if (doNM1Plots) {
      savePlotsDir(h_1d_cr6_metlast_mt2bl_nm1,outfile_,"cr6_metlast_mt2bl_nm1");
      savePlotsDir(h_1d_cr6_metlast_mt_nm1,outfile_,"cr6_metlast_mt_nm1");
      savePlotsDir(h_1d_cr6_metlast_met_nm1,outfile_,"cr6_metlast_met_nm1");
      savePlotsDir(h_1d_cr6_metlast_met100,outfile_,"cr6_metlast_met100");
      savePlotsDir(h_1d_cr6_metlast_met150,outfile_,"cr6_metlast_met150");
    }
    savePlotsDir(h_1d_cr6_metlast_final,outfile_,"cr6_metlast_final");
  }

  if (doCR7) {
    savePlotsDir(h_1d_cr7_presel,outfile_,"cr7_presel");
    if (doNM1Plots) {
      savePlotsDir(h_1d_cr7_met_nm1,outfile_,"cr7_met_nm1");
      savePlotsDir(h_1d_cr7_mt_nm1,outfile_,"cr7_mt_nm1");
      savePlotsDir(h_1d_cr7_mt2bl_nm1,outfile_,"cr7_mt2bl_nm1");
    }
    savePlotsDir(h_1d_cr7_final,outfile_,"cr7_final");
  }

  if (doCR8) {
    savePlotsDir(h_1d_cr8_presel,outfile_,"cr8_presel");
    if (doNM1Plots) {
      savePlotsDir(h_1d_cr8_met_nm1,outfile_,"cr8_met_nm1");
      savePlotsDir(h_1d_cr8_mt_nm1,outfile_,"cr8_mt_nm1");
      savePlotsDir(h_1d_cr8_mt2bl_nm1,outfile_,"cr8_mt2bl_nm1");
    }
    savePlotsDir(h_1d_cr8_final,outfile_,"cr8_final");
  }

  if (doCR8METLast) {
    savePlotsDir(h_1d_cr8_metlast_presel,outfile_,"cr8_metlast_presel");
    if (doNM1Plots) {
      savePlotsDir(h_1d_cr8_metlast_mt2bl_nm1,outfile_,"cr8_metlast_mt2bl_nm1");
      savePlotsDir(h_1d_cr8_metlast_mtfirst,outfile_,"cr8_metlast_mtfirst");
      savePlotsDir(h_1d_cr8_metlast_mt_nm1,outfile_,"cr8_metlast_mt_nm1");
      savePlotsDir(h_1d_cr8_metlast_met_nm1,outfile_,"cr8_metlast_met_nm1");
      savePlotsDir(h_1d_cr8_metlast_met100,outfile_,"cr8_metlast_met100");
      savePlotsDir(h_1d_cr8_metlast_met150,outfile_,"cr8_metlast_met150");
    }
    savePlotsDir(h_1d_cr8_metlast_final,outfile_,"cr8_metlast_final");
  }

  if (doCR9) {
    savePlotsDir(h_1d_cr9_presel,outfile_,"cr9_presel");
    if (doNM1Plots) {
      savePlotsDir(h_1d_cr9_met_nm1,outfile_,"cr9_met_nm1");
      savePlotsDir(h_1d_cr9_mtpeak_nomet,outfile_,"cr9_mtpeak_nomet");
      savePlotsDir(h_1d_cr9_mt_nm1,outfile_,"cr9_mt_nm1");
      savePlotsDir(h_1d_cr9_mtpeak_met,outfile_,"cr9_mtpeak_met");
      savePlotsDir(h_1d_cr9_mt2bl_nm1,outfile_,"cr9_mt2bl_nm1");
      savePlotsDir(h_1d_cr9_mtlast_mt2bl_nm1,outfile_,"cr9_mtlast_mt2bl_nm1");
    }
    savePlotsDir(h_1d_cr9_final,outfile_,"cr9_final");
  }

  if (doCR10) {
    savePlotsDir(h_1d_cr10_presel,outfile_,"cr10_presel");
    if (doNM1Plots) {
      savePlotsDir(h_1d_cr10_met_nm1,outfile_,"cr10_met_nm1");
      savePlotsDir(h_1d_cr10_mt_nm1,outfile_,"cr10_mt_nm1");
      savePlotsDir(h_1d_cr10_mt2bl_nm1,outfile_,"cr10_mt2bl_nm1");
    }
    savePlotsDir(h_1d_cr10_final,outfile_,"cr10_final");
  }

  if (doCR11) {
    savePlotsDir(h_1d_cr11_presel,outfile_,"cr11_presel");
    if (doNM1Plots) {
      savePlotsDir(h_1d_cr11_met_nm1,outfile_,"cr11_met_nm1");
      savePlotsDir(h_1d_cr11_mt_nm1,outfile_,"cr11_mt_nm1");
      savePlotsDir(h_1d_cr11_mt2bl_nm1,outfile_,"cr11_mt2bl_nm1");
    }
    savePlotsDir(h_1d_cr11_final,outfile_,"cr11_final");
  }

  if (doCR12) {
    savePlotsDir(h_1d_cr12_presel,outfile_,"cr12_presel");
    if (doNM1Plots) {
      savePlotsDir(h_1d_cr12_met_nm1,outfile_,"cr12_met_nm1");
      savePlotsDir(h_1d_cr12_mt_nm1,outfile_,"cr12_mt_nm1");
      savePlotsDir(h_1d_cr12_mt2bl_nm1,outfile_,"cr12_mt2bl_nm1");
    }
    savePlotsDir(h_1d_cr12_final,outfile_,"cr12_final");
  }

  if (doCR13) {
    savePlotsDir(h_1d_cr13_presel,outfile_,"cr13_presel");
    if (doNM1Plots) {
      savePlotsDir(h_1d_cr13_mt2bl_nm1,outfile_,"cr13_mt2bl_nm1");
      savePlotsDir(h_1d_cr13_met_nm1,outfile_,"cr13_met_nm1");
      savePlotsDir(h_1d_cr13_met100,outfile_,"cr13_met100");
      savePlotsDir(h_1d_cr13_met150,outfile_,"cr13_met150");
    }
    savePlotsDir(h_1d_cr13_final,outfile_,"cr13_final");
  }

  if (doCR14) {
    savePlotsDir(h_1d_cr14_presel,outfile_,"cr14_presel");
    if (doNM1Plots) {
      savePlotsDir(h_1d_cr14_mt2bl_nm1,outfile_,"cr14_mt2bl_nm1");
      savePlotsDir(h_1d_cr14_mtfirst,outfile_,"cr14_mtfirst");
      savePlotsDir(h_1d_cr14_mt_nm1,outfile_,"cr14_mt_nm1");
      savePlotsDir(h_1d_cr14_nomt_met100,outfile_,"cr14_nomt_met100");
      savePlotsDir(h_1d_cr14_nomt_met150,outfile_,"cr14_nomt_met150");
      savePlotsDir(h_1d_cr14_nomt_met175,outfile_,"cr14_nomt_met175");
      savePlotsDir(h_1d_cr14_met_nm1,outfile_,"cr14_met_nm1");
      savePlotsDir(h_1d_cr14_met100,outfile_,"cr14_met100");
      savePlotsDir(h_1d_cr14_met150,outfile_,"cr14_met150");
    }
    savePlotsDir(h_1d_cr14_final,outfile_,"cr14_final");
  }

  if (doStopSel) {
    savePlotsDir(h_1d_stop_presel,outfile_,"stop_presel");
    savePlotsDir(h_1d_stop_met_nm1,outfile_,"stop_met_nm1");
    savePlotsDir(h_1d_stop_mt_nm1,outfile_,"stop_mt_nm1");
    savePlotsDir(h_1d_stop_isotrk_nm1,outfile_,"stop_isotrk_nm1");
    savePlotsDir(h_1d_stop_tauveto_nm1,outfile_,"stop_tauveto_nm1");
    savePlotsDir(h_1d_stop_comp,outfile_,"stop_comp");
  }

  outfile_->Write();
  outfile_->Close();
  delete outfile_;

  already_seen.clear();

  gROOT->cd();

  bmark->Stop("benchmark");
  cout << endl;
  cout << nEventsTotal << " Events Processed" << endl;
  if (!isData)  cout << nEventsPass << " Events Passed" << endl;
  cout << "------------------------------" << endl;
  cout << "CPU  Time:	" << Form( "%.01f s", bmark->GetCpuTime("benchmark")  ) << ", Rate: " << Form( "%.1f Hz", float(nEventsTotal)/bmark->GetCpuTime("benchmark")) << endl;
  cout << "Real Time:	" << Form( "%.01f s", bmark->GetRealTime("benchmark") ) << ", Rate: " << Form( "%.1f Hz", float(nEventsTotal)/bmark->GetRealTime("benchmark")) << endl;
  cout << endl;
  delete bmark;

}

//--------------------------------------------------------------------

float WHLooper::getCSVCut(const csvpoint csv) {
  float csvcut = 10.;
  if (csv == WHLooper::CSVM) csvcut = 0.679;
  else if (csv == WHLooper::CSVL) csvcut = 0.244;
  else if (csv == WHLooper::CSVT) csvcut = 0.898;

  return csvcut;
}

//--------------------------------------------------------------------

void WHLooper::fillHists1DWrapper(std::map<std::string, TH1F*>& h_1d, const float evtweight, const std::string& dir) {

  fillHists1D(h_1d, evtweight, dir);
  fillFlavorHists1D(h_1d, evtweight, dir);
  if (doFlavorPlots) {
    // single lepton regions: separate into e, m
    if ((dir.find("cr3_") == std::string::npos) && (dir.find("cr4_") == std::string::npos)) {
      if (stopt.leptype() == 0) fillFlavorHists1D(h_1d,evtweight,dir,"_e");
      else if (stopt.leptype() == 1) fillFlavorHists1D(h_1d,evtweight,dir,"_m");
    }

    // for dilepton regions (cr3/4), separate into ee, mm, em
    else {
      int id1 = fabs(stopt.id1());
      int id2 = fabs(stopt.id2());
      if (id1 == 11 && id2 == 11) fillFlavorHists1D(h_1d,evtweight,dir,"_ee");
      else if (id1 == 11 && id2 == 13) fillFlavorHists1D(h_1d,evtweight,dir,"_mm");
      else fillFlavorHists1D(h_1d,evtweight,dir,"_em");
    }
  }

  if (doWJetsPlots && isWjets_ && (dir.find("inc_") != std::string::npos)) {
    if (stopt.nbs() == 0) fillHists1D(h_1d,evtweight,dir,"_0genb");
    else if (stopt.nbs() == 1) fillHists1D(h_1d,evtweight,dir,"_1genb");
    else if (stopt.nbs() == 2) fillHists1D(h_1d,evtweight,dir,"_2genb");
  }

}

//--------------------------------------------------------------------

void WHLooper::fillHists1D(std::map<std::string, TH1F*>& h_1d, const float evtweight, const std::string& dir, const std::string& suffix) {

  outfile_->cd(dir.c_str());

  TVector2 lep(stopt.lep1().px(),stopt.lep1().py());
  TVector2 met;
  met.SetMagPhi(met_, metphi_);
  TVector2 trkmet;
  met.SetMagPhi(stopt.trkmet(), stopt.trkmetphi());
  TVector2 w = lep+met; 
  TVector2 trkmet_nolep = trkmet + lep;

  plot1D("h_njets"+suffix,        njets_,              evtweight, h_1d, 10, 0., 10.);
  plot1D("h_njetsalleta"+suffix,  njetsalleta_,        evtweight, h_1d, 10, 0., 10.);
  plot1D("h_nbjets"+suffix,       nbjets_,    evtweight, h_1d, 5, 0., 5.);

  if (met_ > 100.) {
    plot1D("h_pfmetcalometdphi_met100"+suffix,getdphi(stopt.t1metphicorrphi(), stopt.calometphi()),  evtweight, h_1d, 50, 0., TMath::Pi());
  }
  if (met_ > 150.) {
    plot1D("h_pfmetcalometdphi_met150"+suffix,getdphi(stopt.t1metphicorrphi(), stopt.calometphi()),  evtweight, h_1d, 50, 0., TMath::Pi());
  }
  if (met_ > 175.) {
    plot1D("h_pfmetcalometdphi_met175"+suffix,getdphi(stopt.t1metphicorrphi(), stopt.calometphi()),  evtweight, h_1d, 50, 0., TMath::Pi());
  }

  // phi cor met validation
  // plot1D("h_metdiff"+suffix,        met_ - stopt.pfmet(),    evtweight, h_1d, 500, -250., 250.);
  // plot1D("h_metphidiff"+suffix,  fabs(TVector2::Phi_mpi_pi(metphi_ - stopt.pfmetphi())),    evtweight, h_1d,  50, 0., TMath::Pi());

  plot1D("h_bbmass"+suffix,       bb_.M(),       evtweight, h_1d, 1000, 0., 1000.);
  plot1D("h_bbpt"+suffix,       bb_.pt(),       evtweight, h_1d, 500, 0., 500.);
  // plot1D("h_bblep1dr"+suffix,  ROOT::Math::VectorUtil::DeltaR( bb_ , stopt.lep1() ), evtweight, h_1d, 100, 0., 2.*TMath::Pi());
  // plot1D("h_bblep1dphi"+suffix,  fabs(TVector2::Phi_mpi_pi(bb_.phi() - stopt.lep1().phi())), evtweight, h_1d, 50, 0., TMath::Pi());

  //  plot1D("h_bbwdphi"+suffix,  fabs(TVector2::Phi_mpi_pi(bb_.phi() - w.Phi())), evtweight, h_1d, 50, 0., TMath::Pi());
  plot1D("h_bbwdphi"+suffix,  bbwdphi_, evtweight, h_1d, 50, 0., TMath::Pi());
  // plot1D("h_bbwdpt"+suffix,   bb_.pt() - w.Mod(),       evtweight, h_1d, 500, -250., 250.);

  // jet smearing corrs
  for (unsigned int i=0; i < jets_smearcorrs_.size(); ++i) {
    plot1D("h_jetsmearcorrs"+suffix,   jets_smearcorrs_.at(i),  evtweight, h_1d, 100, 0., 2.);
  }

  if (isWjets_ && (stopt.nbs() == 2)) {
    plot1D("h_genb1pt"+suffix,       stopt.genbs().at(0).pt(),       evtweight, h_1d, 1000, 0., 1000.);
    plot1D("h_genb2pt"+suffix,       stopt.genbs().at(1).pt(),       evtweight, h_1d, 1000, 0., 1000.);
  }

  // plots split by nvtx
  if (doNvtxSplit) {
    if (stopt.nvtx() < 15.) {
      plot1D("h_lep1mt_lowpu"+suffix,       mt_,       evtweight, h_1d, 1000, 0., 1000.);
      plot1D("h_met_lowpu"+suffix,          met_,    evtweight, h_1d, 500, 0., 500.);
      plot1D("h_mt2bl_lowpu"+suffix,  mt2bl_, evtweight, h_1d, 1000, 0., 1000.);
    } else {
      plot1D("h_lep1mt_highpu"+suffix,       mt_,       evtweight, h_1d, 1000, 0., 1000.);
      plot1D("h_met_highpu"+suffix,          met_,    evtweight, h_1d, 500, 0., 500.);
      plot1D("h_mt2bl_highpu"+suffix,  mt2bl_, evtweight, h_1d, 1000, 0., 1000.);
    }
  } // if doNvtxSplit

  // plot tau veto result
  plot1D("h_passtauveto"+suffix,  (int)passTauVeto(),  evtweight, h_1d, 2, 0., 2.);

  if (doWJetsPlots && isWjets_) {
    plot1D("h_nbs",       stopt.nbs(),       evtweight, h_1d, 5, 0, 5);

    // requires babies V20 or higher
    if (stopt.nbs() == 0) {
      plot1D("h_jet1flavor_0b"+suffix, abs(stopt.pfjets_mcflavorAlgo().at(jets_idx_.at(0))) , evtweight, h_1d, 23, -1., 22.);
      plot1D("h_jet2flavor_0b"+suffix, abs(stopt.pfjets_mcflavorAlgo().at(jets_idx_.at(1))) , evtweight, h_1d, 23, -1., 22.);
    } else if (stopt.nbs() == 1) {
      plot1D("h_jet1flavor_1b"+suffix, abs(stopt.pfjets_mcflavorAlgo().at(jets_idx_.at(0))) , evtweight, h_1d, 23, -1., 22.);
      plot1D("h_jet2flavor_1b"+suffix, abs(stopt.pfjets_mcflavorAlgo().at(jets_idx_.at(1))) , evtweight, h_1d, 23, -1., 22.);
    } else if (stopt.nbs() == 2) {
      plot1D("h_jet1flavor_2b"+suffix, abs(stopt.pfjets_mcflavorAlgo().at(jets_idx_.at(0))) , evtweight, h_1d, 23, -1., 22.);
      plot1D("h_jet2flavor_2b"+suffix, abs(stopt.pfjets_mcflavorAlgo().at(jets_idx_.at(1))) , evtweight, h_1d, 23, -1., 22.);
    }

    // requires babies V21 or higher
    // if (stopt.nbs() == 1) plot1D("h_genbjet1pt_1b"+suffix,       stopt.genbs().at(0).pt(),       evtweight, h_1d, 1000, 0., 1000.);
    // if (stopt.nbs() == 2) {
    //   // the genbs vector in the babies isn't presorted by pt..
    //   float genb1pt = stopt.genbs().at(0).pt();
    //   float genb2pt = 0.;
    //   if (stopt.genbs().at(1).pt() > genb1pt) {
    // 	genb1pt = stopt.genbs().at(1).pt();
    // 	genb2pt = stopt.genbs().at(0).pt();
    //   } else {
    // 	genb2pt = stopt.genbs().at(1).pt();
    //   }
    //   plot1D("h_genbjet1pt"+suffix,       genb1pt,       evtweight, h_1d, 1000, 0., 1000.);
    //   plot1D("h_genbjet2pt"+suffix,       genb2pt,       evtweight, h_1d, 1000, 0., 1000.);
    // }
  } // wjets plots

  // bjets and bbbar plots
  if (nbjets_ >= 2) {
    plot1D("h_bjet1pt"+suffix,       bjets_[0].pt(),       evtweight, h_1d, 1000, 0., 1000.);
    plot1D("h_bjet2pt"+suffix,       bjets_[1].pt(),       evtweight, h_1d, 1000, 0., 1000.);
    plot1D("h_bjet1eta"+suffix,       bjets_[0].eta(),       evtweight, h_1d, 100, -3., 3.);
    plot1D("h_bjet2eta"+suffix,       bjets_[1].eta(),       evtweight, h_1d, 100, -3., 3.);

    // // need V00-02-21 or higher babies for these vars
    // plot1D("h_bjet1pumva"+suffix, stopt.pfjets_mva5xPUid().at(bjets_idx_.at(0)) , evtweight, h_1d, 100, -1., 1.);
    // plot1D("h_bjet2pumva"+suffix, stopt.pfjets_mva5xPUid().at(bjets_idx_.at(1)) , evtweight, h_1d, 100, -1., 1.);

    plot1D("h_bjet1metdphi"+suffix,  fabs(TVector2::Phi_mpi_pi(bjets_[0].phi() - metphi_)),  evtweight, h_1d, 50, 0., TMath::Pi());
    plot1D("h_bjet2metdphi"+suffix,  fabs(TVector2::Phi_mpi_pi(bjets_[1].phi() - metphi_)),  evtweight, h_1d, 50, 0., TMath::Pi());

    float bbdr = ROOT::Math::VectorUtil::DeltaR( bjets_.at(0) , bjets_.at(1) );
    plot1D("h_bbdr"+suffix,  bbdr, evtweight, h_1d, 100, 0., 2.*TMath::Pi());
    plot1D("h_bbdphi"+suffix,  fabs(TVector2::Phi_mpi_pi(bjets_[0].phi() - bjets_[1].phi())), evtweight, h_1d, 50, 0., TMath::Pi());
    plot1D("h_bbdeta"+suffix,  bjets_.at(0).eta() - bjets_.at(1).eta() , evtweight, h_1d, 100, -6., 6.);

    float lep1bjet1dr = ROOT::Math::VectorUtil::DeltaR( stopt.lep1(), bjets_.at(0) );
    float lep1bjet2dr = ROOT::Math::VectorUtil::DeltaR( stopt.lep1(), bjets_.at(1) );
    plot1D("h_lep1bjet1dr"+suffix, lep1bjet1dr , evtweight, h_1d, 100, 0., 2.*TMath::Pi());
    plot1D("h_lep1bjet2dr"+suffix, lep1bjet2dr , evtweight, h_1d, 100, 0., 2.*TMath::Pi());

    plot1D("h_lep1bjetmindr"+suffix, TMath::Min(lep1bjet1dr,lep1bjet2dr) , evtweight, h_1d, 100, 0., 2.*TMath::Pi());

    LorentzVector lbb = stopt.lep1() + bjets_.at(0) + bjets_.at(1);
    plot1D("h_lbbpt"+suffix,          lbb.pt(),    evtweight, h_1d, 500, 0., 500.);

    // maria variable: M(bb) * dR(bb) / pt(bb)
    plot1D("h_bbmdrpt"+suffix, bb_.M() * bbdr / bb_.pt(), evtweight, h_1d, 200, 0., 2.*TMath::Pi());

    if (doNvtxSplit) {
      if (stopt.nvtx() < 15.) {
	plot1D("h_lbbpt_lowpu"+suffix,          lbb.pt(),    evtweight, h_1d, 500, 0., 500.);
      } else {
	plot1D("h_lbbpt_highpu"+suffix,          lbb.pt(),    evtweight, h_1d, 500, 0., 500.);
      }
    } // if doNvtxSplit


    // // plots for low pt jets: leading pt < 100
    // //  to investigate data/MC disagreement
    // if (bjets_[0].pt() < 100.) {
    //   plot1D("h_bjet1eta_lowpt"+suffix,       bjets_[0].eta(),       evtweight, h_1d, 100, -3., 3.);
    //   plot1D("h_bjet2eta_lowpt"+suffix,       bjets_[1].eta(),       evtweight, h_1d, 100, -3., 3.);

    //   plot1D("h_bbdr_lowpt"+suffix,  ROOT::Math::VectorUtil::DeltaR( bjets_.at(0) , bjets_.at(1) ), evtweight, h_1d, 100, 0., 2.*TMath::Pi());
    //   plot1D("h_bbdeta_lowpt"+suffix,  bjets_.at(0).eta() - bjets_.at(1).eta() , evtweight, h_1d, 100, -6., 6.);
    // }

    // else {
    //   plot1D("h_bjet1eta_highpt"+suffix,       bjets_[0].eta(),       evtweight, h_1d, 100, -3., 3.);
    //   plot1D("h_bjet2eta_highpt"+suffix,       bjets_[1].eta(),       evtweight, h_1d, 100, -3., 3.);

    //   plot1D("h_bbdr_highpt"+suffix,  ROOT::Math::VectorUtil::DeltaR( bjets_.at(0) , bjets_.at(1) ), evtweight, h_1d, 100, 0., 2.*TMath::Pi());
    //   plot1D("h_bbdeta_highpt"+suffix,  bjets_.at(0).eta() - bjets_.at(1).eta() , evtweight, h_1d, 100, -6., 6.);
    // }

    // LorentzVector b1lep1 = bjets_.at(0) + stopt.lep1();
    // plot1D("h_bjet1lep1mass"+suffix,       b1lep1.M(),       evtweight, h_1d, 1000, 0., 1000.);
    // float bjet1lep1dphi = fabs(TVector2::Phi_mpi_pi(bjets_[0].phi() - stopt.lep1().phi()));
    // plot1D("h_bjet1lep1dphi"+suffix,  bjet1lep1dphi, evtweight, h_1d, 50, 0., TMath::Pi());
    // LorentzVector b2lep1 = bjets_.at(1) + stopt.lep1();
    // plot1D("h_bjet2lep1mass"+suffix,       b2lep1.M(),       evtweight, h_1d, 1000, 0., 1000.);
    // float bjet2lep1dphi = fabs(TVector2::Phi_mpi_pi(bjets_[1].phi() - stopt.lep1().phi()));
    // plot1D("h_bjet2lep1dphi"+suffix,  bjet2lep1dphi, evtweight, h_1d, 50, 0., TMath::Pi());

    // plot1D("h_bjetlep1mindphi"+suffix,  TMath::Min(bjet1lep1dphi,bjet2lep1dphi), evtweight, h_1d, 50, 0., TMath::Pi());

    // use loose btags here in case i plot before requiring 2 med
    // std::vector<int> bjetIdx = getBJetIndex(WHLooper::CSVL,-1,-1);
    // plot1D("h_bjet1mc3"+suffix, stopt.pfjets_mc3().at(bjetIdx.at(0)) , evtweight, h_1d, 40, -20., 20.);
    // plot1D("h_bjet2mc3"+suffix, stopt.pfjets_mc3().at(bjetIdx.at(1)) , evtweight, h_1d, 40, -20., 20.);

    // need V00-02-20 or higher babies for these vars
    // plot1D("h_bjet1flavor"+suffix, abs(stopt.pfjets_mcflavorAlgo().at(bjetIdx.at(0))) , evtweight, h_1d, 23, -1., 22.);
    // plot1D("h_bjet2flavor"+suffix, abs(stopt.pfjets_mcflavorAlgo().at(bjetIdx.at(1))) , evtweight, h_1d, 23, -1., 22.);
  } // if nbjets >= 2

  if (njets_ >= 2) {
    plot1D("h_jet1pt"+suffix,       jets_[0].pt(),       evtweight, h_1d, 1000, 0., 1000.);
    plot1D("h_jet2pt"+suffix,       jets_[1].pt(),       evtweight, h_1d, 1000, 0., 1000.);
    plot1D("h_jet1eta"+suffix,      jets_[0].eta(),      evtweight, h_1d, 100, -3., 3.);
    plot1D("h_jet2eta"+suffix,      jets_[1].eta(),      evtweight, h_1d, 100, -3., 3.);

    plot1D("h_jet1metdphi"+suffix,  fabs(TVector2::Phi_mpi_pi(jets_[0].phi() - metphi_)),  evtweight, h_1d, 50, 0., TMath::Pi());
    plot1D("h_jet2metdphi"+suffix,  fabs(TVector2::Phi_mpi_pi(jets_[1].phi() - metphi_)),  evtweight, h_1d, 50, 0., TMath::Pi());

    float jjdr = ROOT::Math::VectorUtil::DeltaR( jets_.at(0) , jets_.at(1) );
    plot1D("h_jjdr"+suffix,  jjdr, evtweight, h_1d, 100, 0., 2.*TMath::Pi());
    plot1D("h_jjdphi"+suffix,  fabs(TVector2::Phi_mpi_pi(jets_[0].phi() - jets_[1].phi())), evtweight, h_1d, 50, 0., TMath::Pi());

    // maria variable: M(jj) * dR(jj) / pt(jj)
    LorentzVector jj = jets_.at(0) + jets_.at(1);
    plot1D("h_jjmdrpt"+suffix, jj.M() * jjdr / jj.pt(), evtweight, h_1d, 200, 0., 2.*TMath::Pi());

    plot1D("h_jet1csv"+suffix,      jets_csv_.at(0),      evtweight, h_1d, 100, 0., 1.);
    plot1D("h_jet2csv"+suffix,      jets_csv_.at(1),      evtweight, h_1d, 100, 0., 1.);

    // // plot csv in bins of jet pt
    // for (unsigned int ijet = 0; ijet < jets_.size(); ++ijet) {
    //   if (jets_[ijet].pt() < 50.)          plot1D("h_jetcsv_pt0"+suffix,      jets_csv_.at(ijet),      evtweight, h_1d, 100, 0., 1.);
    //   else if (jets_[ijet].pt() < 60.)     plot1D("h_jetcsv_pt1"+suffix,      jets_csv_.at(ijet),      evtweight, h_1d, 100, 0., 1.);
    //   else if (jets_[ijet].pt() < 70.)     plot1D("h_jetcsv_pt2"+suffix,      jets_csv_.at(ijet),      evtweight, h_1d, 100, 0., 1.);
    //   else if (jets_[ijet].pt() < 80.)     plot1D("h_jetcsv_pt3"+suffix,      jets_csv_.at(ijet),      evtweight, h_1d, 100, 0., 1.);
    //   else if (jets_[ijet].pt() < 90.)     plot1D("h_jetcsv_pt4"+suffix,      jets_csv_.at(ijet),      evtweight, h_1d, 100, 0., 1.);
    //   else if (jets_[ijet].pt() < 100.)    plot1D("h_jetcsv_pt5"+suffix,      jets_csv_.at(ijet),      evtweight, h_1d, 100, 0., 1.);
    //   else if (jets_[ijet].pt() < 120.)    plot1D("h_jetcsv_pt6"+suffix,      jets_csv_.at(ijet),      evtweight, h_1d, 100, 0., 1.);
    //   else if (jets_[ijet].pt() < 150.)    plot1D("h_jetcsv_pt7"+suffix,      jets_csv_.at(ijet),      evtweight, h_1d, 100, 0., 1.);
    //   else if (jets_[ijet].pt() < 200.)    plot1D("h_jetcsv_pt8"+suffix,      jets_csv_.at(ijet),      evtweight, h_1d, 100, 0., 1.);
    //   else                                 plot1D("h_jetcsv_pt9"+suffix,      jets_csv_.at(ijet),      evtweight, h_1d, 100, 0., 1.);
    // }

    // need V00-02-21 or higher babies for these vars
    plot1D("h_jet1pumva"+suffix, stopt.pfjets_mva5xPUid().at(jets_idx_.at(0)) , evtweight, h_1d, 100, -1., 1.);
    plot1D("h_jet2pumva"+suffix, stopt.pfjets_mva5xPUid().at(jets_idx_.at(1)) , evtweight, h_1d, 100, -1., 1.);

    plot1D("h_jet1chmneudiff"+suffix, stopt.pfjets_chm().at(jets_idx_.at(0)) - stopt.pfjets_neu().at(jets_idx_.at(0)) , evtweight, h_1d, 200, -50, 150);
    plot1D("h_jet2chmneudiff"+suffix, stopt.pfjets_chm().at(jets_idx_.at(1)) - stopt.pfjets_neu().at(jets_idx_.at(1)) , evtweight, h_1d, 200, -50, 150);

    // // requires babies V28 or higher
    // plot1D("h_jet1tobtecmult"+suffix, stopt.pfjets_tobtecmult().at(jets_idx_.at(0)) , evtweight, h_1d, 100, 0, 100);
    // plot1D("h_jet1tobtecmultfrac"+suffix, stopt.pfjets_tobtecmult().at(jets_idx_.at(0)) / stopt.pfjets_chm().at(jets_idx_.at(0)) , evtweight, h_1d, 50, 0., 1.);
    // plot1D("h_jet1tobtecfrac"+suffix, stopt.pfjets_tobtecfrac().at(jets_idx_.at(0)) , evtweight, h_1d, 50, 0., 1.);
    // plot1D("h_jet2tobtecmult"+suffix, stopt.pfjets_tobtecmult().at(jets_idx_.at(1)) , evtweight, h_1d, 100, 0, 100);
    // plot1D("h_jet2tobtecmultfrac"+suffix, stopt.pfjets_tobtecmult().at(jets_idx_.at(1)) / stopt.pfjets_chm().at(jets_idx_.at(1)) , evtweight, h_1d, 50, 0., 1.);
    // plot1D("h_jet2tobtecfrac"+suffix, stopt.pfjets_tobtecfrac().at(jets_idx_.at(1)) , evtweight, h_1d, 50, 0., 1.);


    // if ( (fabs(jets_[0].eta()) > 0.9) && (fabs(jets_[0].eta()) < 1.9) ) {
    //   plot1D("h_jet1chmneudiff_badeta"+suffix, stopt.pfjets_chm().at(jets_idx_.at(0)) - stopt.pfjets_neu().at(jets_idx_.at(0)) , evtweight, h_1d, 200, -50, 150);
    //   plot1D("h_jet1tobtecmult_badeta"+suffix, stopt.pfjets_tobtecmult().at(jets_idx_.at(0)) , evtweight, h_1d, 100, 0, 100);
    //   plot1D("h_jet1tobtecmultfrac_badeta"+suffix, stopt.pfjets_tobtecmult().at(jets_idx_.at(0)) / stopt.pfjets_chm().at(jets_idx_.at(0)) , evtweight, h_1d, 50, 0., 1.);
    //   plot1D("h_jet1tobtecfrac_badeta"+suffix, stopt.pfjets_tobtecfrac().at(jets_idx_.at(0)) , evtweight, h_1d, 50, 0., 1.);
    // } else {
    //   plot1D("h_jet1chmneudiff_goodeta"+suffix, stopt.pfjets_chm().at(jets_idx_.at(0)) - stopt.pfjets_neu().at(jets_idx_.at(0)) , evtweight, h_1d, 200, -50, 150);
    //   plot1D("h_jet1tobtecmult_goodeta"+suffix, stopt.pfjets_tobtecmult().at(jets_idx_.at(0)) , evtweight, h_1d, 100, 0, 100);
    //   plot1D("h_jet1tobtecmultfrac_goodeta"+suffix, stopt.pfjets_tobtecmult().at(jets_idx_.at(0)) / stopt.pfjets_chm().at(jets_idx_.at(0)) , evtweight, h_1d, 50, 0., 1.);
    //   plot1D("h_jet1tobtecfrac_goodeta"+suffix, stopt.pfjets_tobtecfrac().at(jets_idx_.at(0)) , evtweight, h_1d, 50, 0., 1.);
    // }

    // if ( (fabs(jets_[1].eta()) > 0.9) && (fabs(jets_[1].eta()) < 1.9) ) {
    //   plot1D("h_jet2chmneudiff_badeta"+suffix, stopt.pfjets_chm().at(jets_idx_.at(1)) - stopt.pfjets_neu().at(jets_idx_.at(1)) , evtweight, h_1d, 200, -50, 150);
    //   plot1D("h_jet2tobtecmult_badeta"+suffix, stopt.pfjets_tobtecmult().at(jets_idx_.at(1)) , evtweight, h_1d, 100, 0, 100);
    //   plot1D("h_jet2tobtecmultfrac_badeta"+suffix, stopt.pfjets_tobtecmult().at(jets_idx_.at(1)) / stopt.pfjets_chm().at(jets_idx_.at(1)) , evtweight, h_1d, 50, 0., 1.);
    //   plot1D("h_jet2tobtecfrac_badeta"+suffix, stopt.pfjets_tobtecfrac().at(jets_idx_.at(1)) , evtweight, h_1d, 50, 0., 1.);
    // } else {
    //   plot1D("h_jet2chmneudiff_goodeta"+suffix, stopt.pfjets_chm().at(jets_idx_.at(1)) - stopt.pfjets_neu().at(jets_idx_.at(1)) , evtweight, h_1d, 200, -50, 150);
    //   plot1D("h_jet2tobtecmult_goodeta"+suffix, stopt.pfjets_tobtecmult().at(jets_idx_.at(1)) , evtweight, h_1d, 100, 0, 100);
    //   plot1D("h_jet2tobtecmultfrac_goodeta"+suffix, stopt.pfjets_tobtecmult().at(jets_idx_.at(1)) / stopt.pfjets_chm().at(jets_idx_.at(1)) , evtweight, h_1d, 50, 0., 1.);
    //   plot1D("h_jet2tobtecfrac_goodeta"+suffix, stopt.pfjets_tobtecfrac().at(jets_idx_.at(1)) , evtweight, h_1d, 50, 0., 1.);
    // }

    // if ( stopt.pfjets_beta2_0p5().at(jets_idx_.at(0)) < 0.2 ) {
    //   plot1D("h_jet1pumva_lowbeta"+suffix, stopt.pfjets_mva5xPUid().at(jets_idx_.at(0)) , evtweight, h_1d, 100, -1., 1.);
    // } else {
    //   plot1D("h_jet1pumva_highbeta"+suffix, stopt.pfjets_mva5xPUid().at(jets_idx_.at(0)) , evtweight, h_1d, 100, -1., 1.);
    // }

    // if ( stopt.pfjets_beta2_0p5().at(jets_idx_.at(1)) < 0.2 ) {
    //   plot1D("h_jet2pumva_lowbeta"+suffix, stopt.pfjets_mva5xPUid().at(jets_idx_.at(1)) , evtweight, h_1d, 100, -1., 1.);
    // } else {
    //   plot1D("h_jet2pumva_highbeta"+suffix, stopt.pfjets_mva5xPUid().at(jets_idx_.at(1)) , evtweight, h_1d, 100, -1., 1.);
    // }

    // if ( fabs(jets_[0].eta()) < 2.5 ) {
    //   plot1D("h_jet1pumva_central"+suffix, stopt.pfjets_mva5xPUid().at(jets_idx_.at(0)) , evtweight, h_1d, 100, -1., 1.);
    // } else {
    //   plot1D("h_jet1pumva_forward"+suffix, stopt.pfjets_mva5xPUid().at(jets_idx_.at(0)) , evtweight, h_1d, 100, -1., 1.);
    // }

    // if ( fabs(jets_[1].eta()) < 2.5 ) {
    //   plot1D("h_jet2pumva_central"+suffix, stopt.pfjets_mva5xPUid().at(jets_idx_.at(1)) , evtweight, h_1d, 100, -1., 1.);
    // } else {
    //   plot1D("h_jet2pumva_forward"+suffix, stopt.pfjets_mva5xPUid().at(jets_idx_.at(1)) , evtweight, h_1d, 100, -1., 1.);
    // }

    // need V00-02-20 or higher babies for these vars
    plot1D("h_jet1flavor"+suffix, abs(stopt.pfjets_mcflavorAlgo().at(jets_idx_.at(0))) , evtweight, h_1d, 23, -1., 22.);
    plot1D("h_jet2flavor"+suffix, abs(stopt.pfjets_mcflavorAlgo().at(jets_idx_.at(1))) , evtweight, h_1d, 23, -1., 22.);

    // plots for 3rd, 4th jet
    if (njets_ >= 3) {
      plot1D("h_jet3pt"+suffix,       jets_[2].pt(),       evtweight, h_1d, 1000, 0., 1000.);
      plot1D("h_jet3eta"+suffix,      jets_[2].eta(),      evtweight, h_1d, 100, -3., 3.);
      if (njets_ >= 4) {
	plot1D("h_jet4pt"+suffix,       jets_[3].pt(),       evtweight, h_1d, 1000, 0., 1000.);
	plot1D("h_jet4eta"+suffix,      jets_[3].eta(),      evtweight, h_1d, 100, -3., 3.);
      }
    }


  } // central jets

  if (jets_fwd_.size() > 0) {
    plot1D("h_fwdjet1pt"+suffix,       jets_fwd_[0].pt(),       evtweight, h_1d, 1000, 0., 1000.);
    plot1D("h_fwdjet1eta"+suffix,      jets_fwd_[0].eta(),      evtweight, h_1d, 200, -6., 6.);
    plot1D("h_fwdjet1pumva"+suffix, stopt.pfjets_mva5xPUid().at(jets_fwd_idx_.at(0)) , evtweight, h_1d, 100, -1., 1.);

    if (jets_fwd_.size() > 1) {
      plot1D("h_fwdjet2pt"+suffix,       jets_fwd_[1].pt(),       evtweight, h_1d, 1000, 0., 1000.);
      plot1D("h_fwdjet2eta"+suffix,      jets_fwd_[1].eta(),      evtweight, h_1d, 200, -6., 6.);
      plot1D("h_fwdjet2pumva"+suffix, stopt.pfjets_mva5xPUid().at(jets_fwd_idx_.at(1)) , evtweight, h_1d, 100, -1., 1.);
    }
  } // fwd jets

  //   }

  // }

  plot1D("h_nvtx",      stopt.nvtx(),       evtweight, h_1d, 40, 0, 40);
  plot1D("h_vtxweight", stopt.nvtxweight(), evtweight, h_1d, 41, -4., 4.);
  plot1D("h_rhovor",    stopt.rhovor(),     evtweight, h_1d, 500, 0, 50);


  // plots for CR7 (high mass + 3 jets)
  if (dir.find("cr7") != std::string::npos) {
    if (njets_ >= 3) {
      plot1D("h_jet3pt"+suffix,       jets_[2].pt(),       evtweight, h_1d, 1000, 0., 1000.);
      plot1D("h_jet3eta"+suffix,      jets_[2].eta(),      evtweight, h_1d, 100, -3., 3.);
      plot1D("h_jet3csv"+suffix,      jets_csv_.at(2),      evtweight, h_1d, 100, 0., 1.);
      plot1D("h_jet3flavor"+suffix, abs(stopt.pfjets_mcflavorAlgo().at(jets_idx_.at(2))) , evtweight, h_1d, 23, -1., 22.);
      if (nbjets_ >= 3) {
	plot1D("h_bjet3pt"+suffix,       bjets_[2].pt(),       evtweight, h_1d, 1000, 0., 1000.);
	plot1D("h_bjet3eta"+suffix,       bjets_[2].eta(),       evtweight, h_1d, 100, -3., 3.);
        plot1D("h_bjet3flavor"+suffix, abs(stopt.pfjets_mcflavorAlgo().at(bjets_idx_.at(2))) , evtweight, h_1d, 23, -1., 22.);
      }
    }
  }

  // ttbar plots
  if (isttsl_ || isttdl_) {
    plot1D("h_gentpt"+suffix,     stopt.ptt(),       evtweight, h_1d, 1000, 0., 1000.);
    plot1D("h_gentbarpt"+suffix,     stopt.pttbar(),       evtweight, h_1d, 1000, 0., 1000.);
    plot1D("h_genttbarpt"+suffix,     stopt.ptttbar(),       evtweight, h_1d, 1000, 0., 1000.);
    plot1D("h_topptweight"+suffix,     TopPtWeight(stopt.ptt()),       1., h_1d, 100, 0., 2.);
  }

  // gen mt2bl plot
  if (isttsl_ || isttdl_ || (isWjets_ && stopt.genbs().size() >= 2)) {
    // use lep1 for genmt (and MT2bl) for tt2l..
    plot1D("h_genmt2bl"+suffix,     genmt2bl_,       evtweight, h_1d, 1000, 0., 1000.);
  }

  // gen level mt plots
  if (isttsl_ || isttdl_ || istsl_ || istdl_ || isWjets_) {
    plot1D("h_genmtgenmet"+suffix,     getMT( stopt.mclep1().pt() , stopt.mclep1().phi() , stopt.genmet(), stopt.genmetphi() ),  evtweight, h_1d, 1000, 0., 1000.);
    // single lepton backgrounds: mt(lep,nu)
    if (isttsl_ || istsl_ || isWjets_) {
      plot1D("h_genmtln"+suffix,     stopt.mcmtln(),       evtweight, h_1d, 1000, 0., 1000.);
      plot1D("h_genmln"+suffix,     stopt.mcmln(),       evtweight, h_1d, 1000, 0., 1000.);
    }
  }

  // // tt single lepton plots: find out where missing jets go
  // if (isttsl_ && (njetsalleta_ == 2) && (dir.find("cr1") != std::string::npos)) {

  //   // loop over genjets and genqgs to see which are outside pt, eta acceptance
  //   int ngenjets = 0;
  //   int njets_lepolap = 0;
  //   int njets_outsidept = 0;
  //   int njets_outsideeta = 0;
  //   int njets_outsidepteta = 0;
  //   for ( unsigned int i = 0; i < stopt.genjets().size(); ++i ) {

  //     if (stopt.genjets()[i].pt() < 25.) continue;
  //     // try to match to reco jets
  //     bool matched = false;
  //     for (unsigned int j = 0; j < jets_.size(); ++j) {
  // 	if (ROOT::Math::VectorUtil::DeltaR(stopt.genjets().at(i),jets_.at(j)) < 0.4) {
  // 	  matched = true;
  // 	  break;
  // 	}
  //     }
  //     if (matched) continue;

  //     ++ngenjets;

  //     // try matching to reco leptons, use cone of 0.5 to account for some spread between gen/reco jets..
  //     if (ROOT::Math::VectorUtil::DeltaR(stopt.genjets().at(i),stopt.lep1()) < 0.5) ++njets_lepolap;

  //     // draw pt, eta distributions for gen jets with pt > 25 that aren't selected at reco level

  //     plot1D("h_nonreco_genjetpt"+suffix, stopt.genjets()[i].pt(), evtweight, h_1d, 500, 0., 500.);
  //     plot1D("h_nonreco_genjeteta"+suffix,      stopt.genjets()[i].eta(),      evtweight, h_1d, 100, -3., 3.);

  //     if (stopt.genjets()[i].pt() < 30.) ++njets_outsidept;
  //     if (fabs(stopt.genjets()[i].eta()) > 4.7) ++njets_outsideeta;
  //     if ((stopt.genjets()[i].pt() < 30.) || (fabs(stopt.genjets()[i].eta()) > 4.7) )  ++njets_outsidepteta;

  //     if (stopt.genjets()[i].pt() > 250.) {
  // 	dumpEventInfo("High pt missed gen jet");
  // 	std::cout << " -- missed genjet pt: " << stopt.genjets()[i].pt() << ", eta: " << stopt.genjets()[i].eta() << std::endl;
  //     }

  //   }
  //   plot1D("h_ngenjets_outsidept"+suffix, njets_outsidept, evtweight, h_1d, 5, 0., 5.);
  //   plot1D("h_ngenjets_outsideeta"+suffix, njets_outsideeta, evtweight, h_1d, 5, 0., 5.);
  //   plot1D("h_ngenjets_outsidepteta"+suffix, njets_outsidepteta, evtweight, h_1d, 5, 0., 5.);
  //   plot1D("h_nmissedgenjets"+suffix, ngenjets - njets_lepolap - njets_outsidepteta, evtweight, h_1d, 5, 0., 5.);


  //   int nqgs_outsidept = 0;
  //   int nqgs_outsideeta = 0;
  //   int nqgs_outsidepteta = 0;
  //   for ( unsigned int i = 0; i < stopt.genqgs().size(); ++i ) {
  //     if (stopt.genqgs()[i].pt() < 30.) ++nqgs_outsidept;
  //     if (fabs(stopt.genqgs()[i].eta()) > 4.7) ++nqgs_outsideeta;
  //     if ((stopt.genqgs()[i].pt() < 30.) || (fabs(stopt.genqgs()[i].eta()) > 4.7) )  ++nqgs_outsidepteta;
  //   }
  //   plot1D("h_ngenqgs_outsidept"+suffix, nqgs_outsidept, evtweight, h_1d, 5, 0., 5.);
  //   plot1D("h_ngenqgs_outsideeta"+suffix, nqgs_outsideeta, evtweight, h_1d, 5, 0., 5.);
  //   plot1D("h_ngenqgs_outsidepteta"+suffix, nqgs_outsidepteta, evtweight, h_1d, 5, 0., 5.);

  // }

  // // require babies V27 or higher
  // plot1D("h_mht15"+suffix,          stopt.mht15(),    evtweight, h_1d, 500, 0., 500.);
  // plot1D("h_trkmet_mht15"+suffix,          stopt.trkmet_mht15(),    evtweight, h_1d, 500, 0., 500.);
  // plot1D("h_mettlj15"+suffix,          stopt.mettlj15(),    evtweight, h_1d, 500, 0., 500.);
  // plot1D("h_mttlj15"+suffix, getMT(stopt.lep1().pt(),stopt.lep1().phi(),stopt.mettlj15(),stopt.mettlj15phi()),    evtweight, h_1d, 1000, 0., 1000.);

  // plot true pt of W, Higgs
  if (isTChiwh_) {
    LorentzVector genw = stopt.mclep() + stopt.mcnu();
    plot1D("h_genwpt"+suffix,       genw.pt(),       evtweight, h_1d, 1000, 0., 1000.);
    if (stopt.genbs().size() == 2) {
      LorentzVector genh = stopt.genbs().at(0) + stopt.genbs().at(1);
      plot1D("h_genhpt"+suffix,       genh.pt(),       evtweight, h_1d, 1000, 0., 1000.);
    }
  }

  if (doJetAccPlots && (njetsalleta_ == 2)) fillJetAccHists(h_1d,evtweight,dir,suffix);

  return;
}

//--------------------------------------------------------------------
// hists to fill separate versions based on lepton flavor

void WHLooper::fillFlavorHists1D(std::map<std::string, TH1F*>& h_1d, const float evtweight, const std::string& dir, const std::string& suffix) {

  outfile_->cd(dir.c_str());

  // events histogram: bin 1: raw, bin 2: weighted
  plot1D("h_events"+suffix,       0.5,       1., h_1d, 2, 0., 2.);
  plot1D("h_events"+suffix,       1.5,       evtweight, h_1d, 2, 0., 2.);

  TVector2 lep(stopt.lep1().px(),stopt.lep1().py());
  TVector2 met;
  met.SetMagPhi(met_, metphi_);
  TVector2 trkmet;
  met.SetMagPhi(stopt.trkmet(), stopt.trkmetphi());
  TVector2 w = lep+met; 
  TVector2 trkmet_nolep = trkmet + lep;

  LorentzVector ljets = stopt.lep1();
  for (unsigned int i = 0; i < jets_.size(); ++i) {
    ljets += jets_.at(i);
  }
  plot1D("h_ljetspt"+suffix,  ljets.pt(),    evtweight, h_1d, 500, 0., 500.);

  plot1D("h_lep1pt"+suffix,       stopt.lep1().pt(),       evtweight, h_1d, 1000, 0., 1000.);
  plot1D("h_lep1eta"+suffix,      stopt.lep1().eta(),       evtweight, h_1d, 100, -3., 3.);
  plot1D("h_lep1mt"+suffix,       mt_,       evtweight, h_1d, 1000, 0., 1000.);
  plot1D("h_lep1isopf"+suffix,    stopt.isopf1(),       evtweight, h_1d, 100, 0., 0.5);
  plot1D("h_met"+suffix,          met_,    evtweight, h_1d, 500, 0., 500.);
  plot1D("h_sumet"+suffix,        sumet_,   evtweight, h_1d, 3000, 0., 3000.);
  plot1D("h_sumet_soft"+suffix,   sumet_soft_,   evtweight, h_1d, 3000, 0., 3000.);
  plot1D("h_ht"+suffix,           ht_,   evtweight, h_1d, 1000, 0., 1000.);
  plot1D("h_htlep"+suffix,        ht_ + stopt.lep1().pt(),   evtweight, h_1d, 1000, 0., 1000.);
  plot1D("h_metsig"+suffix,       met_/sqrt(sumet_),   evtweight, h_1d, 500, 0., 25.);
  plot1D("h_pfmet"+suffix,        stopt.pfmet(),    evtweight, h_1d, 500, 0., 500.);
  plot1D("h_t1met10"+suffix,      stopt.t1met10(),    evtweight, h_1d, 500, 0., 500.);
  plot1D("h_trkmet"+suffix,       stopt.trkmet(),    evtweight, h_1d, 500, 0., 500.);
  plot1D("h_met_soft"+suffix,     met_soft_,    evtweight, h_1d, 500, 0., 500.);
  plot1D("h_genmet"+suffix,       stopt.genmet(),    evtweight, h_1d, 500, 0., 500.);
  // plot1D("h_genmet_minus_pfmet"+suffix,       stopt.genmet() - met_,    evtweight, h_1d, 500, -250., 250.);
  // plot1D("h_genmet_minus_pfmet_div_genmet"+suffix,       (stopt.genmet() - met_)/stopt.genmet(),    evtweight, h_1d, 500, -5.0, 5.0);
  // plot1D("h_pfsumet"+suffix,      stopt.pfsumet(),    evtweight, h_1d, 1500, 0., 1500.);
  // plot1D("h_pfmetsig"+suffix,     stopt.pfmet()/sqrt(stopt.pfsumet()),   evtweight, h_1d, 500, 0., 20.);
  plot1D("h_ngoodlep"+suffix,       stopt.ngoodlep(),    evtweight, h_1d, 5, 0., 5.);
  plot1D("h_wpt"+suffix,          wpt_,       evtweight, h_1d, 1000, 0., 1000.);
  plot1D("h_lep1metdphi"+suffix,  lepmetdphi_,  evtweight, h_1d, 50, 0., TMath::Pi());

  plot1D("h_mt2b"+suffix,   mt2b_,  evtweight, h_1d, 1000, 0., 1000.);
  plot1D("h_mt2bl"+suffix,  mt2bl_, evtweight, h_1d, 1000, 0., 1000.);
  plot1D("h_mt2w"+suffix,   mt2w_,  evtweight, h_1d, 1000, 0., 1000.);

  if (stopt.ngoodlep() >= 2) {
    plot1D("h_lep2id"+suffix,   abs(stopt.id2()),  evtweight, h_1d, 3, 11, 14);
  }

  // // plots for CR3/CR4 (2 leptons)
  // if ((dir.find("cr3") != std::string::npos) || (dir.find("cr4") != std::string::npos)) {
  //   plot1D("h_leppt"+suffix,       lep_.pt(),       evtweight, h_1d, 1000, 0., 1000.);
  //   plot1D("h_pseudomt_lep"+suffix,       pseudomt_lep_,       evtweight, h_1d, 1000, 0., 1000.);
  //   plot1D("h_pseudomet_lep"+suffix,        pseudomet_lep_,    evtweight, h_1d, 500, 0., 500.);
  //   plot1D("h_dphi_pseudomet_lep"+suffix,  fabs(dphi_pseudomet_lep_),  evtweight, h_1d, 50, 0., TMath::Pi());
  //   plot1D("h_pseudomt2b"+suffix,   pseudomt2b_,  evtweight, h_1d, 1000, 0., 1000.);
  //   plot1D("h_pseudomt2bl"+suffix,  pseudomt2bl_, evtweight, h_1d, 1000, 0., 1000.);
  //   plot1D("h_pseudomt2w"+suffix,   pseudomt2w_,  evtweight, h_1d, 1000, 0., 1000.);
  // }

  // plots for CR2 (1 lepton + iso track/pfcand)
  if (dir.find("cr2") != std::string::npos) {
    plot1D("h_isotrkpt"+suffix,       stopt.pfcandOS10looseZ().pt(),       evtweight, h_1d, 1000, 0., 1000.);
    plot1D("h_isotrketa"+suffix,      stopt.pfcandOS10looseZ().eta(),       evtweight, h_1d, 100, -3., 3.);
    plot1D("h_lep1isotrkdphi"+suffix,  fabs(TVector2::Phi_mpi_pi(stopt.lep1().phi() - stopt.pfcandOS10looseZ().phi())),  evtweight, h_1d, 50, 0., TMath::Pi());
    if (stopt.pfcandpt5looseZ()  <9998.) {
      plot1D("h_lepcandpt"+suffix,       stopt.pfcand5looseZ().pt(),       evtweight, h_1d, 1000, 0., 1000.);
      plot1D("h_lepcandeta"+suffix,      stopt.pfcand5looseZ().eta(),       evtweight, h_1d, 100, -3., 3.);
      plot1D("h_lep1lepcanddphi"+suffix,  fabs(TVector2::Phi_mpi_pi(stopt.lep1().phi() - stopt.pfcand5looseZ().phi())),  evtweight, h_1d, 50, 0., TMath::Pi());
    }
  }

  // plots for 2 lep events (large overlap with cr3 above, obviously)
  if (stopt.ngoodlep() >= 2) {
    plot1D("h_lep2pt"+suffix,       stopt.lep2().pt(),       evtweight, h_1d, 1000, 0., 1000.);
    plot1D("h_lep2eta"+suffix,      stopt.lep2().eta(),       evtweight, h_1d, 100, -3., 3.);
    plot1D("h_dildphi"+suffix,  fabs(TVector2::Phi_mpi_pi(stopt.lep1().phi() - stopt.lep2().phi())),  evtweight, h_1d, 50, 0., TMath::Pi());
    plot1D("h_dilmass"+suffix,       stopt.dilmass(),       evtweight, h_1d, 1000, 0., 1000.);
    plot1D("h_dilpt"+suffix,       stopt.dilpt(),       evtweight, h_1d, 1000, 0., 1000.);
  }

  return;
}

//--------------------------------------------------------------------

void WHLooper::fillHists2D(std::map<std::string, TH2F*>& h_2d, const float evtweight, const std::string& dir, const std::string& suffix) {

  outfile_->cd(dir.c_str());

  plot2D("h_lep1mt_vs_met"+suffix, met_, mt_, evtweight, h_2d, 1000, 0., 1000., 500, 0., 500.);
  plot2D("h_mt2bl_vs_met"+suffix, met_, mt2bl_, evtweight, h_2d, 1000, 0., 1000., 1000, 0., 1000.);
  plot2D("h_mt2bl_vs_lep1mt"+suffix, mt_, mt2bl_, evtweight, h_2d, 500, 0., 500., 1000, 0., 1000.);

  return;
}

//--------------------------------------------------------------------

void WHLooper::fillJetAccHists(std::map<std::string, TH1F*>& h_1d, const float evtweight, const std::string& dir, const std::string& suffix) {

  int nptacc = 0;
  int netaacc = 0;
  int nlepolap = 0;
  int njetid = 0;
  int npumedid = 0;
  int nputightid = 0;
  int nrestails = 0;
  int nrescore = 0;
  int nmerged = 0;
  int nnoreco = 0;
  int nnorecounder30 = 0;
  int ngood = 0;

  float lepolapdr = 99.;

  // loop over gen jets
  for ( unsigned int igen = 0; igen < stopt.genjets().size(); ++igen ) {
    int qgmatched = 0;

    // loop over status 3 quarks/gluons and look for matches
    for ( unsigned int iqg = 0; iqg < stopt.genqgs().size(); ++iqg ) {
      if (ROOT::Math::VectorUtil::DeltaR(stopt.genjets().at(igen),stopt.genqgs().at(iqg)) < 0.4) ++qgmatched;
    }

    // matched to multiple quarks/gluons: likely jet has merged
    if (qgmatched > 1) {
      ++nmerged;
      //      nmerged += qgmatched - 1;
      // if (qgmatched > 2) {
      // 	std::cout << "found genjet matched to multiple partons: " << qgmatched << std::endl;
      // 	dumpEventInfo("event with merged genjet");
      // }
    } 

    // should do other checks here before dumping genjet?
    if (qgmatched == 0) continue;

    // check to see if genjet is obviously outside acceptance
    if (fabs(stopt.genjets().at(igen).eta()) > 4.7) {
      ++netaacc;
      continue;
    }
    if (stopt.genjets().at(igen).pt() < 20.) {
      ++nptacc;
      continue;
    }

    // may need to do better check to make sure jet isn't actually a lepton..

    // now check to see if we have any reco matches
    //  first check lepton overlap
    bool lepolap = false;
    for (unsigned int ifail = 0; ifail < stopt.pfjets_faillepolap().size(); ++ifail) {
      if (ROOT::Math::VectorUtil::DeltaR(stopt.genjets().at(igen),stopt.pfjets_faillepolap().at(ifail)) < 0.4) {
	float lepdr = ROOT::Math::VectorUtil::DeltaR(stopt.lep1(),stopt.pfjets_faillepolap().at(ifail));
	++nlepolap;
	lepolapdr = lepdr;
	lepolap = true;
      }
    } // reco jets overlapping leptons
    if (lepolap) continue;

    // check for reco jets that failed pf jet id
    bool failid = false;
    for (unsigned int ifail = 0; ifail < stopt.pfjets_failjetid().size(); ++ifail) {
      if (ROOT::Math::VectorUtil::DeltaR(stopt.genjets().at(igen),stopt.pfjets_failjetid().at(ifail)) < 0.4) {
	++njetid;
	failid = true;
      }
    } // reco jets failing jet id
    if (failid) continue;

    // now loop over all reco jets that didn't fail lepton olap or jet id
    bool recomatch = false;
    for (unsigned int ijet = 0; ijet < stopt.pfjets().size(); ++ijet) {
      if (ROOT::Math::VectorUtil::DeltaR(stopt.genjets().at(igen),stopt.pfjets().at(ijet)) < 0.4) {
	recomatch = true;
	// pt acceptance: check for reco pt < 30
	if (stopt.pfjets().at(ijet).pt() < 30.) {
	  // first check for low genjet pt, chalk up to pt acceptance
	  if (stopt.genjets().at(igen).pt() < 25.) {
	    ++nptacc;
	    break;
	  }

	  // check for resolution effects: (genjet pt - reco pt)/sigma, using jet sigma
	  float pterr = stopt.pfjets_sigma().at(ijet) * stopt.pfjets().at(ijet).pt();
	  float sigmadiff = (stopt.genjets().at(igen).pt() - stopt.pfjets().at(ijet).pt()) / pterr;
	  plot1D("h_jetsigmadiff"+suffix,      sigmadiff,    evtweight, h_1d, 100, -6., 6.);

	  if ( fabs(sigmadiff) > 2.0 ) {
	    ++nrestails;
	  } else {
	    ++nrescore;
	  }
	  break;
	} // below pt thresh
	// outside eta (shouldn't happen here..)
	if (fabs(stopt.pfjets().at(ijet).eta()) > 4.7) {
	  ++netaacc;
	  break;
	}
	// fail med pileup id
	if (!passMVAJetId(stopt.pfjets().at(ijet).pt(), stopt.pfjets().at(ijet).eta(), stopt.pfjets_mva5xPUid().at(ijet), 1)) {
	  ++npumedid;
	  break;
	}
	// fail tight pileup id
	if (!passMVAJetId(stopt.pfjets().at(ijet).pt(), stopt.pfjets().at(ijet).eta(), stopt.pfjets_mva5xPUid().at(ijet), 0)) {
	  ++nputightid;
	  break;
	}
	// if we get here, we're matched to a passing reco jet - yay!
	++ngood;

      }
    } // loop over good reco jets
    if (recomatch) continue;

    // no reco match: check if genjet pt is under 30, and chalk up to pt acceptance..
    if (stopt.genjets().at(igen).pt() < 30.) {
      ++nnorecounder30;
      continue;
    }

    // if we get here, we don't know why we lost the genjet and we don't have a reco match..
    ++nnoreco;
    // std::cout << "Didn't match gen jet with pt: " << stopt.genjets().at(igen).pt() << ", eta: " 
    // 	      << stopt.genjets().at(igen).eta() << ", phi: " << stopt.genjets().at(igen).phi() << std::endl;
    // dumpEventInfo("missed gen jet");

  } // loop over genjets

  enum lostjet { PTACC = 0, ETAACC = 1, LEPOLAP = 2, JETID = 3, PUMEDID = 4, PUTIGHTID = 5, RESCORE = 6, RESTAILS = 7, NORECOUNDER30 = 8, NORECO = 9, MERGED = 10 };

  TH1F* h_lostjets = getHist1D("h_lostjets"+suffix, h_1d, 11, 0., 11.);
  if (nptacc) h_lostjets->Fill(PTACC,evtweight*nptacc);
  if (netaacc) h_lostjets->Fill(ETAACC,evtweight*netaacc);
  if (nlepolap) h_lostjets->Fill(LEPOLAP,evtweight*nlepolap);
  if (njetid) h_lostjets->Fill(JETID,evtweight*njetid);
  if (npumedid) h_lostjets->Fill(PUMEDID,evtweight*npumedid);
  if (nputightid) h_lostjets->Fill(PUTIGHTID,evtweight*nputightid);
  if (nrescore) h_lostjets->Fill(RESCORE,evtweight*nrescore);
  if (nrestails) h_lostjets->Fill(RESTAILS,evtweight*nrestails);
  if (nnorecounder30) h_lostjets->Fill(NORECOUNDER30,evtweight*nnorecounder30);
  if (nnoreco) h_lostjets->Fill(NORECO,evtweight*nnoreco);
  if (nmerged) h_lostjets->Fill(MERGED,evtweight*nmerged);

  if (nlepolap) {
    plot1D("h_lepolapdr"+suffix,          lepolapdr,    evtweight, h_1d, 100, 0., 2.*TMath::Pi());
  }

  // this can happen from ISR jets added by pythia, etc. not necessarily useful
  // if (ngood != njetsalleta_) {
  //   std::cout << "WARNING: didn't match all reco jets: found " << ngood 
  // 	      << " matched to genjets, " << njetsalleta_ << " at reco level" << std::endl;
  //   dumpEventInfo("event with not all reco jets matched to gen");
  // }

  return;
}

//--------------------------------------------------------------------

void WHLooper::dumpEventInfo(const std::string& comment) {

  std::cout << "- Dumping info for event: " << comment << std::endl 
	    << "-- run: " << stopt.run() << ", lumi: " <<  stopt.lumi() << ", event: " <<  stopt.event() << std::endl
	    << "-- leptype: " << stopt.leptype() << ", lep pt: " << stopt.lep1().pt() << ", lep eta: " << stopt.lep1().eta() << std::endl
	    << "-- lep2id: " << stopt.id2() << ", lep2 pt: " << stopt.lep2().pt() << ", lep2 eta: " << stopt.lep2().eta() << std::endl
            << "-- met: " << met_ << ", mt: " << mt_ <<  ", mt2bl: " << mt2bl_ <<", njetsalleta: " << njetsalleta_ << std::endl
	    << "-- jet1 pt: " << jets_.at(0).pt() << ", jet1 eta: " << jets_.at(0).eta() << std::endl
    //	    << "-- jet1 chm: " << stopt.pfjets_chm().at(jets_idx_.at(0)) << ", jet1 neu: " << stopt.pfjets_neu().at(jets_idx_.at(0)) << std::endl
	    << "-- jet2 pt: " << jets_.at(1).pt() << ", jet2 eta: " << jets_.at(1).eta() << std::endl;
  //	    << "-- jet2 chm: " << stopt.pfjets_chm().at(jets_idx_.at(1)) << ", jet2 neu: " << stopt.pfjets_neu().at(jets_idx_.at(1)) << std::endl;

  return;
}
