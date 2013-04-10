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

const bool doReshaping = true;
const bool doWJetsOverlapRemoval = true;

const bool doFlavorPlots = true;
const bool doNM1Plots = true;
const bool doWJetsPlots = false;

// regions to do
const bool blindSignal = true;
const bool doInclusive = true;
const bool doSignal = true;
const bool doSignalMassLast = false;
const bool doCR1 = true; // high m(bb)
const bool doCR2 = true; // lep + track
const bool doCR3 = true; // dilep
const bool doCR4 = true; // dilep, high m(bb)
const bool doCR5 = true; // bveto
const bool doCR6 = false; // 1 btag -- not used
const bool doCR7 = true; // high m(bb), 3 jets
const bool doCR8 = true; // low m(bb)
const bool doCR9 = false; // high m(bb), 150-200
const bool doCR10 = false; // high m(bb), 200-250
const bool doStopSel = false;

std::set<DorkyEventIdentifier> already_seen; 
std::set<DorkyEventIdentifier> events_lasercalib; 
std::set<DorkyEventIdentifier> events_hcallasercalib; 

//--------------------------------------------------------------------

// This is meant to be passed as the third argument, the predicate, of the standard library sort algorithm
inline bool sortByPt(const LorentzVector &vec1, const LorentzVector &vec2 ) {
    return vec1.pt() > vec2.pt();
}

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
  CUT_BBMASS_HIGH_ = 140.0;
  CUT_BBMASS_CR1_LOW_ = 150.0;
  CUT_BBMASS_CR1_HIGH_ = 250.0;
  CUT_BBMASS_CR8_LOW_ = 50.0;
  CUT_BBMASS_CR8_HIGH_ = 100.0;
  CUT_MET_PRESEL_ = 50.0;
  CUT_MET_ = 175.0;
  //  CUT_MT_PRESEL_ = 50.0;
  CUT_MT_PRESEL_ = 0.0;
  CUT_MT_ = 100.0;
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

  if (name.Contains("TChiwh")) {
    isTChiwh_ = true;
  } else {
    isTChiwh_ = false;
  }

  if (name.Contains("ttbar_mg")) {
    isttmg_ = true;
  } else {
    isttmg_ = false;
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
  std::map<std::string, TH1F*> h_1d_inc_presel, h_1d_inc_1b, h_1d_inc_2b;

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

  // cr1 hists
  std::map<std::string, TH1F*> h_1d_cr1_presel, h_1d_cr1_final;
  // cr1 nm1 hists
  std::map<std::string, TH1F*> h_1d_cr1_met_nm1, h_1d_cr1_mt_nm1, h_1d_cr1_mt2bl_nm1;
  // cr1 mtpeak hists
  std::map<std::string, TH1F*> h_1d_cr1_mtpeak_nomet, h_1d_cr1_mtpeak_met;

  // cr2 hists
  std::map<std::string, TH1F*> h_1d_cr2_presel, h_1d_cr2_final;
  // cr2 nm1 hists
  std::map<std::string, TH1F*> h_1d_cr2_met_nm1, h_1d_cr2_mt_nm1, h_1d_cr2_mt2bl_nm1;

  // cr3 hists
  std::map<std::string, TH1F*> h_1d_cr3_presel, h_1d_cr3_final;
  // cr3 nm1 hists
  std::map<std::string, TH1F*> h_1d_cr3_met_nm1, h_1d_cr3_mt_nm1, h_1d_cr3_mt2bl_nm1;

  // cr4 hists
  std::map<std::string, TH1F*> h_1d_cr4_presel, h_1d_cr4_final;
  // cr4 nm1 hists
  std::map<std::string, TH1F*> h_1d_cr4_met_nm1, h_1d_cr4_mt_nm1, h_1d_cr4_mt2bl_nm1;

  // cr5 hists
  std::map<std::string, TH1F*> h_1d_cr5_presel, h_1d_cr5_final;
  // cr5 nm1 hists
  std::map<std::string, TH1F*> h_1d_cr5_njetsalleta_nm1, h_1d_cr5_met_nm1, h_1d_cr5_mt_nm1, h_1d_cr5_mt2bl_nm1;

  // cr6 hists
  std::map<std::string, TH1F*> h_1d_cr6_presel, h_1d_cr6_final;
  // cr6 nm1 hists
  std::map<std::string, TH1F*> h_1d_cr6_met_nm1, h_1d_cr6_mt_nm1, h_1d_cr6_mt2bl_nm1;

  // cr7 hists
  std::map<std::string, TH1F*> h_1d_cr7_presel, h_1d_cr7_final;
  // cr7 nm1 hists
  std::map<std::string, TH1F*> h_1d_cr7_met_nm1, h_1d_cr7_mt_nm1, h_1d_cr7_mt2bl_nm1;

  // cr8 hists
  std::map<std::string, TH1F*> h_1d_cr8_presel, h_1d_cr8_final;
  // cr8 nm1 hists
  std::map<std::string, TH1F*> h_1d_cr8_met_nm1, h_1d_cr8_mt_nm1, h_1d_cr8_mt2bl_nm1;

  // cr9 hists
  std::map<std::string, TH1F*> h_1d_cr9_presel, h_1d_cr9_final;
  // cr9 nm1 hists
  std::map<std::string, TH1F*> h_1d_cr9_met_nm1, h_1d_cr9_mt_nm1, h_1d_cr9_mt2bl_nm1;

  // cr10 hists
  std::map<std::string, TH1F*> h_1d_cr10_presel, h_1d_cr10_final;
  // cr10 nm1 hists
  std::map<std::string, TH1F*> h_1d_cr10_met_nm1, h_1d_cr10_mt_nm1, h_1d_cr10_mt2bl_nm1;

  // stop region hists
  std::map<std::string, TH1F*> h_1d_stop_presel, h_1d_stop_comp;
  std::map<std::string, TH1F*> h_1d_stop_met_nm1, h_1d_stop_mt_nm1, h_1d_stop_isotrk_nm1, h_1d_stop_tauveto_nm1 ;


  if (doInclusive) {
    outfile_->mkdir("inc_presel");
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
      outfile_->mkdir("sig_bbmasslast_met_nm1");
      outfile_->mkdir("sig_bbmasslast_mt_nm1");
      outfile_->mkdir("sig_bbmasslast_mt2bl_nm1");
      outfile_->mkdir("sig_bbmasslast_bbmass_nm1");
    }
    outfile_->mkdir("sig_bbmasslast_final");
  }

  if (doCR1) {
    outfile_->mkdir("cr1_presel");
    if (doNM1Plots) {
      outfile_->mkdir("cr1_met_nm1");
      outfile_->mkdir("cr1_mtpeak_nomet");
      outfile_->mkdir("cr1_mt_nm1");
      outfile_->mkdir("cr1_mtpeak_met");
      outfile_->mkdir("cr1_mt2bl_nm1");
    }
    outfile_->mkdir("cr1_final");
  }

  if (doCR2) {
    outfile_->mkdir("cr2_presel");
    if (doNM1Plots) {
      outfile_->mkdir("cr2_met_nm1");
      outfile_->mkdir("cr2_mt_nm1");
      outfile_->mkdir("cr2_mt2bl_nm1");
    }
    outfile_->mkdir("cr2_final");
  }

  if (doCR3) {
    outfile_->mkdir("cr3_presel");
    if (doNM1Plots) {
      outfile_->mkdir("cr3_met_nm1");
      outfile_->mkdir("cr3_mt_nm1");
      outfile_->mkdir("cr3_mt2bl_nm1");
    }
    outfile_->mkdir("cr3_final");
  }

  if (doCR4) {
    outfile_->mkdir("cr4_presel");
    if (doNM1Plots) {
      outfile_->mkdir("cr4_met_nm1");
      outfile_->mkdir("cr4_mt_nm1");
      outfile_->mkdir("cr4_mt2bl_nm1");
    }
    outfile_->mkdir("cr4_final");
  }

  if (doCR5) {
    outfile_->mkdir("cr5_presel");
    if (doNM1Plots) {
      outfile_->mkdir("cr5_njetsalleta_nm1");
      outfile_->mkdir("cr5_met_nm1");
      outfile_->mkdir("cr5_mt_nm1");
      outfile_->mkdir("cr5_mt2bl_nm1");
    }
    outfile_->mkdir("cr5_final");
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

  if (doCR9) {
    outfile_->mkdir("cr9_presel");
    if (doNM1Plots) {
      outfile_->mkdir("cr9_met_nm1");
      outfile_->mkdir("cr9_mt_nm1");
      outfile_->mkdir("cr9_mt2bl_nm1");
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
      if (isttmg_) evtweight /= stopt.mgcor();
      // cross section weights for TChiwh samples:
      //  xsec (pb) * 1000 (pb to fb) * br(w->lv) 0.33 * br(h->bb) 0.58 / nevents (10000)
      float weight_lumi_br_nevents = 1.914E-02;
      if (isTChiwh_) {
	float xsecweight = 1.;
	if (name.Contains("TChiwh_150")) xsecweight = 2.4 * 0.5 * weight_lumi_br_nevents; // now 20k events
	else if (name.Contains("TChiwh_200")) xsecweight = 0.79 * 0.5 * weight_lumi_br_nevents; // now 20k events
	else if (name.Contains("TChiwh_250")) xsecweight = 0.32 * weight_lumi_br_nevents;
	else if (name.Contains("TChiwh_300")) xsecweight = 0.15 * weight_lumi_br_nevents;
	else if (name.Contains("TChiwh_350")) xsecweight = 0.074 * weight_lumi_br_nevents;

	// reset weight here for signal MC. Ignore weight, nvtxweight, mgcor
	evtweight = xsecweight * lumi;
      }

      plot1D("h_nvtx_nosel",       stopt.nvtx(),       evtweight, h_1d_sig_presel, 40, 0, 40);
      plot1D("h_vtxweight_nosel", stopt.nvtxweight(), evtweight, h_1d_sig_presel, 41, -4., 4.);

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

      if ( !passEvtSelection(name) ) continue;

      // data run selection: to run on certain periods
      // if (isData && (stopt.run() > 203755)) continue; // periods A-C
      // if (isData && (stopt.run() < 203768)) continue; // period D

      // require first vertex to be good (i.e. index == 0)
      //   to make sure the same vertex is used for btagging, leptons, etc
      //      if (stopt.indexfirstGoodVertex_()) continue;

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
      njets_ = 0;
      njetsalleta_ = 0;
      nbjets_ = 0;
      nbjetsl_ = 0;

      for( unsigned int i = 0 ; i < stopt.pfjets().size() ; ++i ){
	
	// basic jet selection
	if( stopt.pfjets().at(i).pt()<30 )  continue;
	if( fabs(stopt.pfjets().at(i).eta())>4.7 )  continue;
	//	if ( (fabs(stopt.pfjets().at(i).eta()) < 2.5) 
	//	     && (stopt.pfjets_beta2_0p5().at(i)<0.2) ) continue;
	// pileup MVA ID: use tight working point (= 0)
	if (!passMVAJetId(stopt.pfjets().at(i).pt(), stopt.pfjets().at(i).eta(), stopt.pfjets_mva5xPUid().at(i), 0)) continue;
	
	// count jets from all eta, save jets within |eta| < 2.4
	++njetsalleta_;
	if (fabs(stopt.pfjets().at(i).eta()) > 2.4) {
	  jets_fwd_.push_back( stopt.pfjets().at(i) );
	  jets_fwd_idx_.push_back(i);
	  continue;
	}
	++njets_;

	//	if (fabs(stopt.pfjets().at(i).eta()) < 2.4) ++njets_;

	jets_.push_back( stopt.pfjets().at(i) );
	jets_idx_.push_back(i);

	// if (n_jets==1) 
	//   dphimj1 = getdphi(t1metphicorrphi, stopt.pfjets().at(i).phi() );
	// if (n_jets==2) {
	//   dphimj2 = getdphi(t1metphicorrphi, stopt.pfjets().at(i).phi() );
	//   dphimjmin = TMath::Min( dphimj1 , dphimj2 );
	// }

	float csv_nominal= stopt.pfjets_csv().at(i);

	//RESHAPING -- requires babies V20 or higher
	//  -- may need to validate for non-b-jets
	if (doReshaping && !isData && !isTChiwh_) {
	  csv_nominal = nominalShape->reshape( stopt.pfjets().at(i).eta(),
					       stopt.pfjets().at(i).pt(),
					       stopt.pfjets_csv().at(i),
					       stopt.pfjets_mcflavorAlgo().at(i) ); 
	}

	jets_csv_.push_back( csv_nominal );

	if( (fabs(stopt.pfjets().at(i).eta()) <= 2.4) && (csv_nominal > getCSVCut(WHLooper::CSVM)) ) {
	  bjets_.push_back( stopt.pfjets().at(i) );
	  bjets_idx_.push_back(i);
	  ++nbjets_;
	}
	else if( (fabs(stopt.pfjets().at(i).eta()) <= 2.4) && (csv_nominal > getCSVCut(WHLooper::CSVL)) ) {
	  ++nbjetsl_;
	}

	//      	sigma_jets.push_back(stopt.pfjets_sigma().at(i));

        // //count jets that are not overlapping with second lepton
	// if (isData) continue;
	// if (stopt.nleps()!=2) continue;
	// if (stopt.mclep2().pt() < 30.) continue;
	// if (ROOT::Math::VectorUtil::DeltaR(stopt.mclep2(), stopt.pfjets().at(i)) > 0.4 ) continue;
	// n_ljets--;

      } // loop over pfjets

      met_ = stopt.t1metphicorr();
      metphi_ = stopt.t1metphicorrphi();
      mt_ = stopt.t1metphicorrmt();

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

      // TVector2 lep(stopt.lep1().px(),stopt.lep1().py());
      // TVector2 met;
      // met.SetMagPhi(stopt.t1metphicorr(), metphi_);
      // TVector2 w = lep+met; 

      // end variables --------------------------------------
      // ----------------------------------------------------

      // ----------------------------------------------------
      // selections bits

      //      bool passisotrk = passIsoTrkVeto_v2();
      //      bool passisotrk = passIsoTrkVeto_v3() && (stopt.ngoodlep() == 1);
      bool passisotrk = passIsoTrkVeto_v4();
      //      bool passisotrk = passIsoTrkVeto_v4() && passTauVeto();

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
      //   100 < m(bb) < 140
      //   met > 50

      // *** signal region:
      //   == 1 lepton, iso track veto
      //   == 2 jets, all eta
      //   == 2 bjets
      //   100 < m(bb) < 140
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
	if (!fail && (mt_ > 50.) && (mt_ < 100.) ) {
	  if (doNM1Plots) fillHists1DWrapper(h_1d_sig_mtpeak_nomet,evtweight1l,"sig_mtpeak_nomet");
	}

	if (!fail && (met_ > CUT_MET_) ) {
	  if (doNM1Plots) fillHists1DWrapper(h_1d_sig_mt_nm1,evtweight1l,"sig_mt_nm1");
	}
	else fail = true;

	// mt peak after cutting on met
	if (!fail && (mt_ > 50.) && (mt_ < 100.) ) {
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
	   && (mt_ > CUT_MT_PRESEL_) 
	   && (!isData || !blindSignal) ) {

        fillHists1DWrapper(h_1d_sig_bbmasslast_presel,evtweight1l,"sig_bbmasslast_presel");

	bool fail = false;
	if ( !fail && (njetsalleta_ == 2) ) {
	  if (doNM1Plots) fillHists1DWrapper(h_1d_sig_bbmasslast_met_nm1,evtweight1l,"sig_bbmasslast_met_nm1");
	}
	else fail = true;

	if (!fail && (met_ > CUT_MET_) ) {
	  if (doNM1Plots) fillHists1DWrapper(h_1d_sig_bbmasslast_mt_nm1,evtweight1l,"sig_bbmasslast_mt_nm1");
	}
	else fail = true;

	if (!fail && (mt_ > CUT_MT_) ) {
	  if (doNM1Plots) fillHists1DWrapper(h_1d_sig_bbmasslast_mt2bl_nm1,evtweight1l,"sig_bbmasslast_mt2bl_nm1");
	}
	else fail = true;

	if (!fail && (mt2bl_ > CUT_MT2BL_) ) {
	  fillHists1DWrapper(h_1d_sig_bbmasslast_bbmass_nm1,evtweight1l,"sig_bbmasslast_bbmass_nm1");
	}
	else fail = true;

	if (!fail && (bb_.M() > CUT_BBMASS_LOW_) && (bb_.M() < CUT_BBMASS_HIGH_) ) {
	  fillHists1DWrapper(h_1d_sig_bbmasslast_final,evtweight1l,"sig_bbmasslast_final");
	}

      } // signal region sel (bbmass last)

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

	// mt peak before cutting on met
	if (!fail && (mt_ > 50.) && (mt_ < 100.) ) {
	  if (doNM1Plots) fillHists1DWrapper(h_1d_cr1_mtpeak_nomet,evtweight1l,"cr1_mtpeak_nomet");
	}

	if (!fail && (met_ > CUT_MET_) ) {
	  if (doNM1Plots) fillHists1DWrapper(h_1d_cr1_mt_nm1,evtweight1l,"cr1_mt_nm1");
	}
	else fail = true;

	// mt peak after cutting on met
	if (!fail && (mt_ > 50.) && (mt_ < 100.) ) {
	  if (doNM1Plots) fillHists1DWrapper(h_1d_cr1_mtpeak_met,evtweight1l,"cr1_mtpeak_met");
	}

	if (!fail && (mt_ > CUT_MT_) ) {
	  if (doNM1Plots) fillHists1DWrapper(h_1d_cr1_mt2bl_nm1,evtweight1l,"cr1_mt2bl_nm1");
	}
	else fail = true;

	if (!fail && (mt2bl_ > CUT_MT2BL_) ) {
	  fillHists1DWrapper(h_1d_cr1_final,evtweight1l,"cr1_final");
	}

      } // CR1 region sel


      // -------------------------------------------
      // *** CR2: 1 lepton + iso track, i.e. inverted iso track veto, require exactly 1 lep
      //   otherwise same as signal region

      if ( doCR2
	   && passSingleLeptonSelection(isData) 
	   && !passisotrk && (stopt.ngoodlep() == 1) 
	   && (nbjets_ >= 2)
	   && (bb_.M() > CUT_BBMASS_LOW_) && (bb_.M() < CUT_BBMASS_HIGH_)
	   && (met_ > CUT_MET_PRESEL_) 
	   && (mt_ > CUT_MT_PRESEL_) ) {

        fillHists1DWrapper(h_1d_cr2_presel,evtweight1l,"cr2_presel");

	bool fail = false;
	if ( !fail && (njetsalleta_ == 2) ) {
	  if (doNM1Plots) fillHists1DWrapper(h_1d_cr2_met_nm1,evtweight1l,"cr2_met_nm1");
	}
	else fail = true;

	if (!fail && (met_ > CUT_MET_) ) {
	  if (doNM1Plots) fillHists1DWrapper(h_1d_cr2_mt_nm1,evtweight1l,"cr2_mt_nm1");
	}
	else fail = true;

	if (!fail && (mt_ > CUT_MT_) ) {
	  if (doNM1Plots) fillHists1DWrapper(h_1d_cr2_mt2bl_nm1,evtweight1l,"cr2_mt2bl_nm1");
	}
	else fail = true;

	if (!fail && (mt2bl_ > CUT_MT2BL_) ) {
	  fillHists1DWrapper(h_1d_cr2_final,evtweight1l,"cr2_final");
	}

      } // CR2 region sel


      // -------------------------------------------
      // *** CR3: 2 leptons, Z mass veto for same flavor, no iso track veto
      //   otherwise same as signal region
      //   !!! pesudo MET/MT etc defs to emulate losing 2nd lepton

      if ( doCR3
	   && passDileptonSelection(isData) 
	   && (abs(stopt.id1()) != abs(stopt.id2()) || fabs( stopt.dilmass() - 91.) > 15. )
	   && (nbjets_ >= 2)
	   && (bb_.M() > CUT_BBMASS_LOW_) && (bb_.M() < CUT_BBMASS_HIGH_) ) {

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
	  if (doNM1Plots) fillHists1DWrapper(h_1d_cr3_met_nm1,evtweight2l,"cr3_met_nm1");
	}
	else fail = true;

	//	if (!fail && isData && met_ > 270.) dumpEventInfo("CR3, high MET event");

	//	if (!fail && (pseudomet_lep_ > CUT_MET_) ) {
	if (!fail && (met_ > CUT_MET_) ) {
	  if (doNM1Plots) fillHists1DWrapper(h_1d_cr3_mt_nm1,evtweight2l,"cr3_mt_nm1");
	}
	else fail = true;

	//	if (!fail && (pseudomt_lep_ > CUT_MT_) ) {
	if (!fail && (mt_ > CUT_MT_) ) {
	  if (doNM1Plots) fillHists1DWrapper(h_1d_cr3_mt2bl_nm1,evtweight2l,"cr3_mt2bl_nm1");
	}
	else fail = true;

	//	if (!fail && (pseudomt2bl_ > CUT_MT2BL_) ) {
	if (!fail && (mt2bl_ > CUT_MT2BL_) ) {
	  fillHists1DWrapper(h_1d_cr3_final,evtweight2l,"cr3_final");
	}

      } // CR3 region sel


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
	  if (doNM1Plots) fillHists1DWrapper(h_1d_cr4_met_nm1,evtweight2l,"cr4_met_nm1");
	}
	else fail = true;

	//	if (!fail && isData && met_ > 320.) dumpEventInfo("CR4, high MET event");

	//	if (!fail && (pseudomet_lep_ > CUT_MET_) ) {
	if (!fail && (met_ > CUT_MET_) ) {
	  if (doNM1Plots) fillHists1DWrapper(h_1d_cr4_mt_nm1,evtweight2l,"cr4_mt_nm1");
	}
	else fail = true;

	//	if (!fail && (pseudomt_lep_ > CUT_MT_) ) {
	if (!fail && (mt_ > CUT_MT_) ) {
	  if (doNM1Plots) fillHists1DWrapper(h_1d_cr4_mt2bl_nm1,evtweight2l,"cr4_mt2bl_nm1");
	}
	else fail = true;

	//	if (!fail && (pseudomt2bl_ > CUT_MT2BL_) ) {
	if (!fail && (mt2bl_ > CUT_MT2BL_) ) {
	  fillHists1DWrapper(h_1d_cr4_final,evtweight2l,"cr4_final");
	}

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
	if ( !fail && ((bb_.M() > CUT_BBMASS_LOW_) && (bb_.M() < CUT_BBMASS_HIGH_)) ) {
	  if (doNM1Plots) fillHists1DWrapper(h_1d_cr5_njetsalleta_nm1,evtweight1l,"cr5_njetsalleta_nm1");
	}
	else fail = true;

	if ( !fail && (njetsalleta_ == 2) ) {
	  if (doNM1Plots) fillHists1DWrapper(h_1d_cr5_met_nm1,evtweight1l,"cr5_met_nm1");
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
	   && (bb_.M() > 150.) && (bb_.M() < 200.) 
	   && (met_ > CUT_MET_PRESEL_) 
	   && (mt_ > CUT_MT_PRESEL_) ) {

        fillHists1DWrapper(h_1d_cr9_presel,evtweight1l,"cr9_presel");

	bool fail = false;
	if ( !fail && (njetsalleta_ == 2) ) {
	  if (doNM1Plots) fillHists1DWrapper(h_1d_cr9_met_nm1,evtweight1l,"cr9_met_nm1");
	}
	else fail = true;

	if (!fail && (met_ > CUT_MET_) ) {
	  if (doNM1Plots) fillHists1DWrapper(h_1d_cr9_mt_nm1,evtweight1l,"cr9_mt_nm1");
	}
	else fail = true;

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

  if (doSignal) {
    savePlotsDir(h_1d_inc_presel,outfile_,"inc_presel");
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
      savePlotsDir(h_1d_sig_bbmasslast_met_nm1,outfile_,"sig_bbmasslast_met_nm1");
      savePlotsDir(h_1d_sig_bbmasslast_mt_nm1,outfile_,"sig_bbmasslast_mt_nm1");
      savePlotsDir(h_1d_sig_bbmasslast_mt2bl_nm1,outfile_,"sig_bbmasslast_mt2bl_nm1");
      savePlotsDir(h_1d_sig_bbmasslast_bbmass_nm1,outfile_,"sig_bbmasslast_bbmass_nm1");
    }
    savePlotsDir(h_1d_sig_bbmasslast_final,outfile_,"sig_bbmasslast_final");
  }

  if (doCR1) {
    savePlotsDir(h_1d_cr1_presel,outfile_,"cr1_presel");
    if (doNM1Plots) {
      savePlotsDir(h_1d_cr1_met_nm1,outfile_,"cr1_met_nm1");
      savePlotsDir(h_1d_cr1_mtpeak_nomet,outfile_,"cr1_mtpeak_nomet");
      savePlotsDir(h_1d_cr1_mt_nm1,outfile_,"cr1_mt_nm1");
      savePlotsDir(h_1d_cr1_mtpeak_met,outfile_,"cr1_mtpeak_met");
      savePlotsDir(h_1d_cr1_mt2bl_nm1,outfile_,"cr1_mt2bl_nm1");
    }
    savePlotsDir(h_1d_cr1_final,outfile_,"cr1_final");
  }

  if (doCR2) {
    savePlotsDir(h_1d_cr2_presel,outfile_,"cr2_presel");
    if (doNM1Plots) {
      savePlotsDir(h_1d_cr2_met_nm1,outfile_,"cr2_met_nm1");
      savePlotsDir(h_1d_cr2_mt_nm1,outfile_,"cr2_mt_nm1");
      savePlotsDir(h_1d_cr2_mt2bl_nm1,outfile_,"cr2_mt2bl_nm1");
    }
    savePlotsDir(h_1d_cr2_final,outfile_,"cr2_final");
  }

  if (doCR3) {
    savePlotsDir(h_1d_cr3_presel,outfile_,"cr3_presel");
    if (doNM1Plots) {
      savePlotsDir(h_1d_cr3_met_nm1,outfile_,"cr3_met_nm1");
      savePlotsDir(h_1d_cr3_mt_nm1,outfile_,"cr3_mt_nm1");
      savePlotsDir(h_1d_cr3_mt2bl_nm1,outfile_,"cr3_mt2bl_nm1");
    }
    savePlotsDir(h_1d_cr3_final,outfile_,"cr3_final");
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
      savePlotsDir(h_1d_cr5_njetsalleta_nm1,outfile_,"cr5_njetsalleta_nm1");
      savePlotsDir(h_1d_cr5_met_nm1,outfile_,"cr5_met_nm1");
      savePlotsDir(h_1d_cr5_mt_nm1,outfile_,"cr5_mt_nm1");
      savePlotsDir(h_1d_cr5_mt2bl_nm1,outfile_,"cr5_mt2bl_nm1");
    }
    savePlotsDir(h_1d_cr5_final,outfile_,"cr5_final");
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

  if (doCR9) {
    savePlotsDir(h_1d_cr9_presel,outfile_,"cr9_presel");
    if (doNM1Plots) {
      savePlotsDir(h_1d_cr9_met_nm1,outfile_,"cr9_met_nm1");
      savePlotsDir(h_1d_cr9_mt_nm1,outfile_,"cr9_mt_nm1");
      savePlotsDir(h_1d_cr9_mt2bl_nm1,outfile_,"cr9_mt2bl_nm1");
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
  if (doFlavorPlots) {
    // single lepton regions: separate into e, m
    if ((dir.find("cr3_") == std::string::npos) && (dir.find("cr4_") == std::string::npos)) {
      if (stopt.leptype() == 0) fillHists1D(h_1d,evtweight,dir,"_e");
      else if (stopt.leptype() == 1) fillHists1D(h_1d,evtweight,dir,"_m");
    }

    // for dilepton regions (cr3/4), separate into ee, mm, em
    else {
      int id1 = fabs(stopt.id1());
      int id2 = fabs(stopt.id2());
      if (id1 == 11 && id2 == 11) fillHists1D(h_1d,evtweight,dir,"_ee");
      else if (id1 == 11 && id2 == 13) fillHists1D(h_1d,evtweight,dir,"_mm");
      else fillHists1D(h_1d,evtweight,dir,"_em");
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

  // events histogram: bin 1: raw, bin 2: weighted
  plot1D("h_events"+suffix,       0.5,       1., h_1d, 2, 0., 2.);
  plot1D("h_events"+suffix,       1.5,       evtweight, h_1d, 2, 0., 2.);

  TVector2 lep(stopt.lep1().px(),stopt.lep1().py());
  TVector2 met;
  met.SetMagPhi(met_, metphi_);
  TVector2 w = lep+met; 

  plot1D("h_lep1pt"+suffix,       stopt.lep1().pt(),       evtweight, h_1d, 1000, 0., 1000.);
  plot1D("h_lep1eta"+suffix,      stopt.lep1().eta(),       evtweight, h_1d, 100, -3., 3.);
  plot1D("h_lep1mt"+suffix,       mt_,       evtweight, h_1d, 1000, 0., 1000.);
  plot1D("h_lep1isopf"+suffix,    stopt.isopf1(),       evtweight, h_1d, 100, 0., 0.5);
  plot1D("h_met"+suffix,          met_,    evtweight, h_1d, 500, 0., 500.);
  plot1D("h_pfmet"+suffix,        stopt.pfmet(),    evtweight, h_1d, 500, 0., 500.);
  plot1D("h_t1met10"+suffix,      stopt.t1met10(),    evtweight, h_1d, 500, 0., 500.);
  plot1D("h_genmet"+suffix,       stopt.genmet(),    evtweight, h_1d, 500, 0., 500.);
  // plot1D("h_genmet_minus_pfmet"+suffix,       stopt.genmet() - met_,    evtweight, h_1d, 500, -250., 250.);
  // plot1D("h_genmet_minus_pfmet_div_genmet"+suffix,       (stopt.genmet() - met_)/stopt.genmet(),    evtweight, h_1d, 500, -5.0, 5.0);
  // plot1D("h_pfsumet"+suffix,      stopt.pfsumet(),    evtweight, h_1d, 1500, 0., 1500.);
  // plot1D("h_pfmetsig"+suffix,     stopt.pfmet()/sqrt(stopt.pfsumet()),   evtweight, h_1d, 500, 0., 20.);
  plot1D("h_njets"+suffix,        njets_,              evtweight, h_1d, 10, 0., 10.);
  plot1D("h_njetsalleta"+suffix,  njetsalleta_,        evtweight, h_1d, 10, 0., 10.);
  plot1D("h_nbjets"+suffix,       nbjets_,    evtweight, h_1d, 5, 0., 5.);
  plot1D("h_ngoodlep"+suffix,       stopt.ngoodlep(),    evtweight, h_1d, 5, 0., 5.);
  plot1D("h_wpt"+suffix,          w.Mod(),       evtweight, h_1d, 1000, 0., 1000.);
  plot1D("h_lep1metdphi"+suffix,  fabs(TVector2::Phi_mpi_pi(stopt.lep1().phi() - metphi_)),  evtweight, h_1d, 50, 0., TMath::Pi());

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

  plot1D("h_bbwdphi"+suffix,  fabs(TVector2::Phi_mpi_pi(bb_.phi() - w.Phi())), evtweight, h_1d, 50, 0., TMath::Pi());
  // plot1D("h_bbwdpt"+suffix,   bb_.pt() - w.Mod(),       evtweight, h_1d, 500, -250., 250.);

  plot1D("h_mt2b"+suffix,   mt2b_,  evtweight, h_1d, 1000, 0., 1000.);
  plot1D("h_mt2bl"+suffix,  mt2bl_, evtweight, h_1d, 1000, 0., 1000.);
  plot1D("h_mt2w"+suffix,   mt2w_,  evtweight, h_1d, 1000, 0., 1000.);

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

    plot1D("h_bbdr"+suffix,  ROOT::Math::VectorUtil::DeltaR( bjets_.at(0) , bjets_.at(1) ), evtweight, h_1d, 100, 0., 2.*TMath::Pi());
    plot1D("h_bbdphi"+suffix,  fabs(TVector2::Phi_mpi_pi(bjets_[0].phi() - bjets_[1].phi())), evtweight, h_1d, 50, 0., TMath::Pi());
    plot1D("h_bbdeta"+suffix,  bjets_.at(0).eta() - bjets_.at(1).eta() , evtweight, h_1d, 100, -6., 6.);

    float lep1bjet1dr = ROOT::Math::VectorUtil::DeltaR( stopt.lep1(), bjets_.at(0) );
    float lep1bjet2dr = ROOT::Math::VectorUtil::DeltaR( stopt.lep1(), bjets_.at(1) );
    plot1D("h_lep1bjet1dr"+suffix, lep1bjet1dr , evtweight, h_1d, 100, 0., 2.*TMath::Pi());
    plot1D("h_lep1bjet2dr"+suffix, lep1bjet2dr , evtweight, h_1d, 100, 0., 2.*TMath::Pi());

    plot1D("h_lep1bjetmindr"+suffix, TMath::Min(lep1bjet1dr,lep1bjet2dr) , evtweight, h_1d, 100, 0., 2.*TMath::Pi());

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

    plot1D("h_jjdr"+suffix,  ROOT::Math::VectorUtil::DeltaR( jets_.at(0) , jets_.at(1) ), evtweight, h_1d, 100, 0., 2.*TMath::Pi());
    plot1D("h_jjdphi"+suffix,  fabs(TVector2::Phi_mpi_pi(jets_[0].phi() - jets_[1].phi())), evtweight, h_1d, 50, 0., TMath::Pi());

    plot1D("h_jet1csv"+suffix,      jets_csv_.at(0),      evtweight, h_1d, 100, 0., 1.);
    plot1D("h_jet2csv"+suffix,      jets_csv_.at(1),      evtweight, h_1d, 100, 0., 1.);

    // need V00-02-21 or higher babies for these vars
    plot1D("h_jet1pumva"+suffix, stopt.pfjets_mva5xPUid().at(jets_idx_.at(0)) , evtweight, h_1d, 100, -1., 1.);
    plot1D("h_jet2pumva"+suffix, stopt.pfjets_mva5xPUid().at(jets_idx_.at(1)) , evtweight, h_1d, 100, -1., 1.);

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

  // plots for CR2 (1 lepton + iso track/pfcand)
  if (dir.find("cr2") != std::string::npos) {
    plot1D("h_isotrkpt"+suffix,       stopt.pfcandOS10looseZ().pt(),       evtweight, h_1d, 1000, 0., 1000.);
    plot1D("h_isotrketa"+suffix,      stopt.pfcandOS10looseZ().eta(),       evtweight, h_1d, 100, -3., 3.);
    plot1D("h_lep1isotrkdphi"+suffix,  fabs(TVector2::Phi_mpi_pi(stopt.lep1().phi() - stopt.pfcandOS10looseZ().phi())),  evtweight, h_1d, 50, 0., TMath::Pi());
  }

  // plots for CR3/CR4 (2 leptons)
  if ((dir.find("cr3") != std::string::npos) || (dir.find("cr4") != std::string::npos)) {
    plot1D("h_leppt"+suffix,       lep_.pt(),       evtweight, h_1d, 1000, 0., 1000.);
    plot1D("h_pseudomt_lep"+suffix,       pseudomt_lep_,       evtweight, h_1d, 1000, 0., 1000.);
    plot1D("h_pseudomet_lep"+suffix,        pseudomet_lep_,    evtweight, h_1d, 500, 0., 500.);
    plot1D("h_dphi_pseudomet_lep"+suffix,  fabs(dphi_pseudomet_lep_),  evtweight, h_1d, 50, 0., TMath::Pi());
    plot1D("h_pseudomt2b"+suffix,   pseudomt2b_,  evtweight, h_1d, 1000, 0., 1000.);
    plot1D("h_pseudomt2bl"+suffix,  pseudomt2bl_, evtweight, h_1d, 1000, 0., 1000.);
    plot1D("h_pseudomt2w"+suffix,   pseudomt2w_,  evtweight, h_1d, 1000, 0., 1000.);
  }

  // plots for 2 lep events (large overlap with cr3 above, obviously)
  if (stopt.ngoodlep() >= 2) {
    plot1D("h_lep2pt"+suffix,       stopt.lep2().pt(),       evtweight, h_1d, 1000, 0., 1000.);
    plot1D("h_lep2eta"+suffix,      stopt.lep2().eta(),       evtweight, h_1d, 100, -3., 3.);
    plot1D("h_dildphi"+suffix,  fabs(TVector2::Phi_mpi_pi(stopt.lep1().phi() - stopt.lep2().phi())),  evtweight, h_1d, 50, 0., TMath::Pi());
    plot1D("h_dilmass"+suffix,       stopt.dilmass(),       evtweight, h_1d, 1000, 0., 1000.);
    plot1D("h_dilpt"+suffix,       stopt.dilpt(),       evtweight, h_1d, 1000, 0., 1000.);
  }

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

  return;
}

//--------------------------------------------------------------------

void WHLooper::dumpEventInfo(const std::string& comment) {

  std::cout << "- Dumping info for event: " << comment << std::endl 
	    << "-- run: " << stopt.run() << ", event: " <<  stopt.event() << ", lumi: " <<  stopt.lumi() << std::endl
	    << "-- leptype: " << stopt.leptype() << ", lep pt: " << stopt.lep1().pt() << ", lep eta: " << stopt.lep1().eta() << std::endl
            << "-- met: " << met_ << ", mt: " << mt_ << ", njetsalleta: " << njetsalleta_ << std::endl
	    << "-- jet1 pt: " << jets_.at(0).pt() << ", jet1 eta: " << jets_.at(0).eta() << std::endl
	    << "-- jet2 pt: " << jets_.at(1).pt() << ", jet2 eta: " << jets_.at(1).eta() << std::endl;

  return;
}
