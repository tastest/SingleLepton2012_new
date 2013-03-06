#include "WHLooper.h"

//#include "../../CORE/jetSmearingTools.h"
//#include "../../CORE/Thrust.h"
//#include "../../CORE/EventShape.h"

#include "../Core/STOPT.h"
#include "../Core/stopUtils.h"
#include "../Plotting/PlotUtilities.h"
#include "../Core/MT2Utility.h"
#include "../Core/mt2bl_bisect.h"
#include "../Core/mt2w_bisect.h"
#include "../Core/PartonCombinatorics.h"

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

const bool doFlavorPlots = false;
const bool doNM1Plots = true;
const bool doMiniBaby = false;

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
  min_mtpeak = -9999.;
  max_mtpeak = -9999.; 
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

  if (name.Contains("TChiwh")) {
    isTChiwh_ = true;
  } else {
    isTChiwh_ = false;
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

  //------------------------------
  // set up histograms
  //------------------------------

  if (doMiniBaby) {
    string outbabyfile = m_outfilename_;
    replace(outbabyfile,"histos","minibaby");
    cout << "[WHLooper::loop] creating minibaby file: " << outbabyfile << endl;
    MakeBabyNtuple(outbabyfile.c_str());
  }

  gROOT->cd();

  cout << "[WHLooper::loop] creating output file: " << m_outfilename_ << endl;

  outfile_ = new TFile(m_outfilename_.c_str(),"RECREATE") ; 

  cout << "[WHLooper::loop] setting up histos" << endl;

  std::map<std::string, TH1F*> h_1d_presel;
  std::map<std::string, TH1F*> h_1d_presel_nobs;
  std::map<std::string, TH1F*> h_1d_presel_bs;
  std::map<std::string, TH1F*> h_1d_njets_nm1;
  std::map<std::string, TH1F*> h_1d_njetsalleta_nm1;
  std::map<std::string, TH1F*> h_1d_bbmass_nm1;
  std::map<std::string, TH1F*> h_1d_bbmass_nobs_nm1;
  std::map<std::string, TH1F*> h_1d_bbmass_bs_nm1;
  std::map<std::string, TH1F*> h_1d_pfmet_nm1;
  std::map<std::string, TH1F*> h_1d_lep1mt_nm1;
  std::map<std::string, TH1F*> h_1d_mt2bl_nm1;
  // std::map<std::string, TH1F*> h_1d_mt2w_nm1;
  // std::map<std::string, TH1F*> h_1d_bbpt_nm1;
  // std::map<std::string, TH1F*> h_1d_wpt_nm1;
  // std::map<std::string, TH1F*> h_1d_bbwdphi_nm1;
  std::map<std::string, TH1F*> h_1d_final;
  std::map<std::string, TH1F*> h_1d_final_nobs;
  std::map<std::string, TH1F*> h_1d_final_bs;
  std::map<std::string, TH2F*> h_2d_final;

  outfile_->mkdir("presel");
  if (isWjets_) {
    outfile_->mkdir("presel_nobs");
    outfile_->mkdir("presel_bs");
  }
  if (doNM1Plots) {
    outfile_->mkdir("njets_nm1");
    outfile_->mkdir("njetsalleta_nm1");
    outfile_->mkdir("bbmass_nm1");
    if (isWjets_) {
      outfile_->mkdir("bbmass_nobs_nm1");
      outfile_->mkdir("bbmass_bs_nm1");
    }
    outfile_->mkdir("pfmet_nm1");
    outfile_->mkdir("lep1mt_nm1");
    outfile_->mkdir("mt2bl_nm1");
    // outfile_->mkdir("mt2w_nm1");
    // outfile_->mkdir("bbpt_nm1");
    // outfile_->mkdir("wpt_nm1");
    // outfile_->mkdir("bbwdphi_nm1");
  }
  outfile_->mkdir("final");
  if (isWjets_) {
    outfile_->mkdir("final_nobs");
    outfile_->mkdir("final_bs");
  }
  outfile_->cd("presel");

  //------------------------------
  // vtx reweighting
  //------------------------------

  // TFile* vtx_file = TFile::Open("../vtxreweight/vtxreweight_Summer12_DR53X-PU_S10_9p7ifb_Zselection.root");
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

      //---------------------------------------------------------------------------- 
      // determine event weight
      // make 2 example histograms of nvtx and corresponding weight
      //---------------------------------------------------------------------------- 

      float evtweight = isData ? 1. : ( stopt.weight() * 19.5 * stopt.nvtxweight() * stopt.mgcor() );
      // cross section weights for TChiwh samples:
      //  xsec (pb) * 1000 (pb to fb) * br(w->lv) 0.33 * br(h->bb) 0.58 / nevents (10000)
      float weight_lumi_br_nevents = 1.914E-02;
      if (isTChiwh_) {
	if (name.Contains("TChiwh_150_1")) evtweight *= 2.4 * 0.5 * weight_lumi_br_nevents; // now 20k events
	else if (name.Contains("TChiwh_200_1")) evtweight *= 0.79 * 0.5 * weight_lumi_br_nevents; // now 20k events
	else if (name.Contains("TChiwh_250_1")) evtweight *= 0.32 * weight_lumi_br_nevents;
	else if (name.Contains("TChiwh_300_1")) evtweight *= 0.15 * weight_lumi_br_nevents;
	else if (name.Contains("TChiwh_350_1")) evtweight *= 0.074 * weight_lumi_br_nevents;

	// divide back out nvtxweight -- shouldn't be applied to signal
	evtweight /=  stopt.nvtxweight();
      }
      // to reweight from file - also need to comment stuff before
      //      float vtxweight = vtxweight_n( nvtx, h_vtx_wgt, isData );

      plot1D("h_nvtx_nosel",       stopt.nvtx(),       evtweight, h_1d_presel, 40, 0, 40);
      plot1D("h_vtxweight_nosel", stopt.nvtxweight(), evtweight, h_1d_presel, 41, -4., 4.);

      //----------------------------------------------------------------------------
      // apply preselection:
      // rho 0-40 GeV, MET filters, >=1 good lepton, veto 2 leptons dR < 0.1
      //----------------------------------------------------------------------------

      if ( !passEvtSelection(name) ) continue;

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

      // require 1 lepton sel + iso track veto
      if (!passSingleLeptonSelection(isData)) continue;
      if (stopt.ngoodlep() != 1) continue;
      if (!passIsoTrkVeto_v3()) continue;

      // set mt2 vars to dummy value
      mt2b_ = -1.;
      mt2bl_ = -1.;
      mt2w_ = -1.;

      // require 2 bjets
      myBJets_ = getBJets(WHLooper::CSVM);
      //      myBJets_ = getBJets(WHLooper::CSVL);
      int njets = getNJets();
      int njetsalleta = getNJets(4.7);
      int nbjets = myBJets_.size();
      //      if (njetsalleta != 2) continue;
      //      if (stopt.pfmet() < 150.) continue;
      if (nbjets < 2) continue;

      // specific cuts: dummy signal region
      //      if (nbjets < 2) continue;
      std::vector<LorentzVector> bjets_csvm = getBJets(WHLooper::CSVM);
      int nbjets_csvm = bjets_csvm.size();
      LorentzVector bb = myBJets_.at(0) + myBJets_.at(1);
      float lep1mt = getMT(stopt.lep1().pt(), stopt.lep1().phi(), stopt.t1metphicorr(), stopt.t1metphicorrphi() );
      TVector2 lep(stopt.lep1().px(),stopt.lep1().py());
      TVector2 met;
      met.SetMagPhi(stopt.t1metphicorr(), stopt.t1metphicorrphi());
      TVector2 w = lep+met; 

      if (doNM1Plots) fillHists1DWrapper(h_1d_njets_nm1,evtweight,"njets_nm1");
      if (njets > 2) continue;
      if (doNM1Plots) fillHists1DWrapper(h_1d_njetsalleta_nm1,evtweight,"njetsalleta_nm1");
      if (njetsalleta > 2) continue;

      //      if (doNM1Plots) fillHists1DWrapper(h_1d_pfmet_nm1,evtweight,"pfmet_nm1");
      if (stopt.t1metphicorr() < 50.) continue;

      // calculate mt2w, after requiring exactly two jets
      // make dummy vector of csv values for bjets, for mt2w calc
      std::vector<float> bjets_csv(2, 0.99);
      mt2b_ = calculateMT2w(bjets_csvm, bjets_csv, stopt.lep1(), stopt.t1metphicorr(), stopt.t1metphicorrphi(), MT2b);
      mt2bl_ = calculateMT2w(bjets_csvm, bjets_csv, stopt.lep1(), stopt.t1metphicorr(), stopt.t1metphicorrphi(), MT2bl);
      mt2w_ = calculateMT2w(bjets_csvm, bjets_csv, stopt.lep1(), stopt.t1metphicorr(), stopt.t1metphicorrphi(), MT2w);

      // consider tightening bbmass to 110,140 for low mass sel??
      if (doNM1Plots) fillHists1DWrapper(h_1d_bbmass_nm1,evtweight,"bbmass_nm1");
      if (doNM1Plots && isWjets_) {
	if (stopt.nbs() == 0) fillHists1DWrapper(h_1d_bbmass_nobs_nm1,evtweight,"bbmass_nobs_nm1");
	else fillHists1DWrapper(h_1d_bbmass_bs_nm1,evtweight,"bbmass_bs_nm1");
      }
      //           if (bb.M() < 95. || bb.M() > 150.) continue;
      //      if (bb.M() < 100. || bb.M() > 140.) continue;
      if (bb.M() < 150.) continue;

      // --- consider above cuts preselection, fill presel histos here
      // also fill custom minibaby if using
      if (doMiniBaby) FillBabyNtuple(evtweight);
      fillHists1DWrapper(h_1d_presel,evtweight,"presel");
      if (isWjets_) {
      	if (stopt.nbs() == 0) fillHists1DWrapper(h_1d_presel_nobs,evtweight,"presel_nobs");
      	else fillHists1DWrapper(h_1d_presel_bs,evtweight,"presel_bs");
      }

      if (doNM1Plots) fillHists1DWrapper(h_1d_pfmet_nm1,evtweight,"pfmet_nm1");
      if (stopt.t1metphicorr() < 175.) continue;
      //      if (stopt.t1metphicorr() < 150.) continue;
      // if (stopt.pfmet() < 80.) continue;
      if (doNM1Plots) fillHists1DWrapper(h_1d_lep1mt_nm1,evtweight,"lep1mt_nm1");
      //      if (lep1mt < 120.) continue;
      if (lep1mt < 100.) continue;
      if (doNM1Plots) fillHists1DWrapper(h_1d_mt2bl_nm1,evtweight,"mt2bl_nm1");
      //      if (mt2bl_ < 175.) continue;
      if (mt2bl_ < 200.) continue;
      // if (doNM1Plots) fillHists1DWrapper(h_1d_mt2w_nm1,evtweight,"mt2w_nm1");
      // if (mt2w_ < 175.) continue;
      // //      if (mt2w_ < 200.) continue;
      // if (doNM1Plots) fillHists1DWrapper(h_1d_bbpt_nm1,evtweight,"bbpt_nm1");
      // if (bb.pt() < 150.) continue;
      // //      if (bb.pt() < 175.) continue;
      // //      if (bb.pt() < 220.) continue;
      // if (doNM1Plots) fillHists1DWrapper(h_1d_wpt_nm1,evtweight,"wpt_nm1");
      // if (w.Mod() < 150.) continue;
      // //      if (w.Mod() < 220.) continue;
      // if (doNM1Plots) fillHists1DWrapper(h_1d_bbwdphi_nm1,evtweight,"bbwdphi_nm1");
      // if (fabs(TVector2::Phi_mpi_pi(bb.phi() - w.Phi())) < 2.95) continue;

      ++nEventsPass;

      // fill hists
      plot1D("h_nvtx",       stopt.nvtx(),       evtweight, h_1d_final, 40, 0, 40);
      plot1D("h_vtxweight", stopt.nvtxweight(), evtweight, h_1d_final, 41, -4., 4.);

      fillHists1DWrapper(h_1d_final,evtweight,"final");
      if (isWjets_) {
	if (stopt.nbs() == 0) fillHists1DWrapper(h_1d_final_nobs,evtweight,"final_nobs");
	else fillHists1DWrapper(h_1d_final_bs,evtweight,"final_bs");
      }

      plot2D("h_wpt_vs_bbpt", w.Mod(), bb.pt(), evtweight, h_2d_final, 500, 0, 1000, 500, 0, 1000);

    } // end event loop

    // delete tree;
    
  } // end file loop
  
    //
    // finish
    //

  //  savePlots12(h_1d, h_2d, (char*)m_outfilename_.c_str());

  savePlotsDir(h_1d_presel,outfile_,"presel");
  if (isWjets_) {
    savePlotsDir(h_1d_presel_nobs,outfile_,"presel_nobs");
    savePlotsDir(h_1d_presel_bs,outfile_,"presel_bs");
  }
  if (doNM1Plots) {
    savePlotsDir(h_1d_njets_nm1,outfile_,"njets_nm1");
    savePlotsDir(h_1d_njetsalleta_nm1,outfile_,"njetsalleta_nm1");
    savePlotsDir(h_1d_bbmass_nm1,outfile_,"bbmass_nm1");
    if (isWjets_) {
      savePlotsDir(h_1d_bbmass_nobs_nm1,outfile_,"bbmass_nobs_nm1");
      savePlotsDir(h_1d_bbmass_bs_nm1,outfile_,"bbmass_bs_nm1");
    }
    savePlotsDir(h_1d_pfmet_nm1,outfile_,"pfmet_nm1");
    savePlotsDir(h_1d_lep1mt_nm1,outfile_,"lep1mt_nm1");
    savePlotsDir(h_1d_mt2bl_nm1,outfile_,"mt2bl_nm1");
    // savePlotsDir(h_1d_mt2w_nm1,outfile_,"mt2w_nm1");
    // savePlotsDir(h_1d_bbpt_nm1,outfile_,"bbpt_nm1");
    // savePlotsDir(h_1d_wpt_nm1,outfile_,"wpt_nm1");
    // savePlotsDir(h_1d_bbwdphi_nm1,outfile_,"bbwdphi_nm1");
  }

  savePlotsDir(h_1d_final,outfile_,"final");
  savePlots2Dir(h_2d_final,outfile_,"final");
  if (isWjets_) {
    savePlotsDir(h_1d_final_nobs,outfile_,"final_nobs");
    savePlotsDir(h_1d_final_bs,outfile_,"final_bs");
  }  

  outfile_->Write();
  outfile_->Close();
  delete outfile_;

  already_seen.clear();

  gROOT->cd();

  if (doMiniBaby) CloseBabyNtuple();

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

std::vector<LorentzVector> WHLooper::getBJets(const csvpoint csv) {

  std::vector<LorentzVector> bjets;
  float csvcut = getCSVCut(csv);

  std::vector<int> bjetIdx = getBJetIndex(csvcut,-1,-1);

  for (unsigned int i=0; i<bjetIdx.size(); ++i) {
    bjets.push_back(stopt.pfjets().at(i));
  }

  // sort by pt -- shouldn't be needed
  //  sort(bjets.begin()  , bjets.end()  , sortByPt);

  return bjets;

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
    if (stopt.leptype() == 0) fillHists1D(h_1d,evtweight,dir,"_e");
    else if (stopt.leptype() == 1) fillHists1D(h_1d,evtweight,dir,"_m");
  }
}

//--------------------------------------------------------------------

void WHLooper::fillHists1D(std::map<std::string, TH1F*>& h_1d, const float evtweight, const std::string& dir, const std::string& suffix) {

  outfile_->cd(dir.c_str());

  float lep1mt = getMT(stopt.lep1().pt(), stopt.lep1().phi(), stopt.t1metphicorr(), stopt.t1metphicorrphi() );
  int njets = getNJets();
  int njetsalleta = getNJets(4.7);
  TVector2 lep(stopt.lep1().px(),stopt.lep1().py());
  TVector2 met;
  met.SetMagPhi(stopt.t1metphicorr(), stopt.t1metphicorrphi());
  TVector2 w = lep+met; 

  plot1D("h_lep1pt"+suffix,       stopt.lep1().pt(),       evtweight, h_1d, 1000, 0., 1000.);
  plot1D("h_lep1eta"+suffix,      stopt.lep1().eta(),       evtweight, h_1d, 100, -3., 3.);
  plot1D("h_lep1mt"+suffix,       lep1mt,       evtweight, h_1d, 1000, 0., 1000.);
  plot1D("h_pfmet"+suffix,        stopt.t1metphicorr(),    evtweight, h_1d, 500, 0., 500.);
  plot1D("h_pfsumet"+suffix,      stopt.pfsumet(),    evtweight, h_1d, 1500, 0., 1500.);
  plot1D("h_pfmetsig"+suffix,     stopt.pfmet()/sqrt(stopt.pfsumet()),   evtweight, h_1d, 500, 0., 20.);
  plot1D("h_njets"+suffix,        njets,              evtweight, h_1d, 10, 0., 10.);
  plot1D("h_njetsalleta"+suffix,  njetsalleta,        evtweight, h_1d, 10, 0., 10.);
  plot1D("h_nbjets"+suffix,       myBJets_.size(),    evtweight, h_1d, 5, 0., 5.);
  plot1D("h_wpt"+suffix,          w.Mod(),       evtweight, h_1d, 1000, 0., 1000.);
  plot1D("h_lep1metdphi"+suffix,  fabs(TVector2::Phi_mpi_pi(stopt.lep1().phi() - stopt.t1metphicorrphi())),  evtweight, h_1d, 50, 0., TMath::Pi());

  // phi cor met validation
  plot1D("h_metdiff"+suffix,        stopt.t1metphicorr() - stopt.pfmet(),    evtweight, h_1d, 500, -250., 250.);
  plot1D("h_metphidiff"+suffix,  fabs(TVector2::Phi_mpi_pi(stopt.t1metphicorrphi() - stopt.pfmetphi())),    evtweight, h_1d,  50, 0., TMath::Pi());


  if (isWjets_) {
    plot1D("h_nbs",       stopt.nbs(),       evtweight, h_1d, 5, 0, 5);
  }

  // bjets and bbbar plots
  if (myBJets_.size() >= 2) {
    plot1D("h_bjet1pt"+suffix,       myBJets_[0].pt(),       evtweight, h_1d, 500, 0., 500.);
    plot1D("h_bjet2pt"+suffix,       myBJets_[1].pt(),       evtweight, h_1d, 500, 0., 500.);
    plot1D("h_bjet1eta"+suffix,       myBJets_[0].eta(),       evtweight, h_1d, 100, -3., 3.);
    plot1D("h_bjet2eta"+suffix,       myBJets_[1].eta(),       evtweight, h_1d, 100, -3., 3.);

    LorentzVector bb = myBJets_.at(0) + myBJets_.at(1);
    plot1D("h_bbmass"+suffix,       bb.M(),       evtweight, h_1d, 1000, 0., 1000.);
    plot1D("h_bbpt"+suffix,       bb.pt(),       evtweight, h_1d, 500, 0., 500.);
    plot1D("h_bbdr"+suffix,  ROOT::Math::VectorUtil::DeltaR( myBJets_.at(0) , myBJets_.at(1) ), evtweight, h_1d, 100, 0., 2.*TMath::Pi());
    plot1D("h_bbdphi"+suffix,  fabs(TVector2::Phi_mpi_pi(myBJets_[0].phi() - myBJets_[1].phi())), evtweight, h_1d, 50, 0., TMath::Pi());

    plot1D("h_bblep1dr"+suffix,  ROOT::Math::VectorUtil::DeltaR( bb , stopt.lep1() ), evtweight, h_1d, 100, 0., 2.*TMath::Pi());
    plot1D("h_bblep1dphi"+suffix,  fabs(TVector2::Phi_mpi_pi(bb.phi() - stopt.lep1().phi())), evtweight, h_1d, 50, 0., TMath::Pi());

    plot1D("h_bbwdphi"+suffix,  fabs(TVector2::Phi_mpi_pi(bb.phi() - w.Phi())), evtweight, h_1d, 50, 0., TMath::Pi());
    plot1D("h_bbwdpt"+suffix,   bb.pt() - w.Mod(),       evtweight, h_1d, 500, -250., 250.);

    plot1D("h_bbwsumpt"+suffix,       bb.pt()+w.Mod(),       evtweight, h_1d, 1000, 0., 1000.);

    plot1D("h_allsumpt"+suffix,      myBJets_[0].pt()+ myBJets_[1].pt()+stopt.lep1().pt()+stopt.t1metphicorr() , evtweight, h_1d, 1500, 0., 1500.);

    LorentzVector b1lep1 = myBJets_.at(0) + stopt.lep1();
    plot1D("h_bjet1lep1mass"+suffix,       b1lep1.M(),       evtweight, h_1d, 1000, 0., 1000.);
    float bjet1lep1dphi = fabs(TVector2::Phi_mpi_pi(myBJets_[0].phi() - stopt.lep1().phi()));
    plot1D("h_bjet1lep1dphi"+suffix,  bjet1lep1dphi, evtweight, h_1d, 50, 0., TMath::Pi());
    LorentzVector b2lep1 = myBJets_.at(1) + stopt.lep1();
    plot1D("h_bjet2lep1mass"+suffix,       b2lep1.M(),       evtweight, h_1d, 1000, 0., 1000.);
    float bjet2lep1dphi = fabs(TVector2::Phi_mpi_pi(myBJets_[1].phi() - stopt.lep1().phi()));
    plot1D("h_bjet2lep1dphi"+suffix,  bjet2lep1dphi, evtweight, h_1d, 50, 0., TMath::Pi());

    plot1D("h_bjetlep1mindphi"+suffix,  TMath::Min(bjet1lep1dphi,bjet2lep1dphi), evtweight, h_1d, 50, 0., TMath::Pi());

    plot1D("h_mt2b"+suffix,   mt2b_,  evtweight, h_1d, 1000, 0., 1000.);
    plot1D("h_mt2bl"+suffix,  mt2bl_, evtweight, h_1d, 1000, 0., 1000.);
    plot1D("h_mt2w"+suffix,   mt2w_,  evtweight, h_1d, 1000, 0., 1000.);

    // use loose btags here in case i plot before requiring 2 med
    std::vector<int> bjetIdx = getBJetIndex(WHLooper::CSVL,-1,-1);
    plot1D("h_bjet1mc3"+suffix, stopt.pfjets_mc3().at(bjetIdx.at(0)) , evtweight, h_1d, 40, -20., 20.);
    plot1D("h_bjet2mc3"+suffix, stopt.pfjets_mc3().at(bjetIdx.at(1)) , evtweight, h_1d, 40, -20., 20.);

    // need V00-02-20 or higher babies for these vars
    //    plot1D("h_bjet1flavor"+suffix, abs(stopt.pfjets_mcflavorAlgo().at(bjetIdx.at(0))) , evtweight, h_1d, 22, 0., 22.);
    //    plot1D("h_bjet2flavor"+suffix, abs(stopt.pfjets_mcflavorAlgo().at(bjetIdx.at(1))) , evtweight, h_1d, 22, 0., 22.);
  }

  return;
}

//______________________________________________________________________________________
void WHLooper::MakeBabyNtuple (const char* babyFileName)
{

  TDirectory *rootdir = gDirectory->GetDirectory("Rint:");
  rootdir->cd();

  babyFile_ = new TFile(Form("%s", babyFileName), "RECREATE");
  babyFile_->cd();
  babyTree_ = new TTree("T1", "A Baby Ntuple");

  //event stuff
  babyTree_->Branch("run",          &run_,          "run/I"          );
  babyTree_->Branch("lumi",         &lumi_,         "lumi/I"         );
  babyTree_->Branch("event",        &event_,        "event/I"        );
  babyTree_->Branch("leptype",               &leptype_,               "leptype/I");
  babyTree_->Branch("weight",       &weight_,       "weight/F"       );

  //vars
  babyTree_->Branch("pfmet",        &pfmet_,        "pfmet/F"      );
  babyTree_->Branch("lep1mt",       &lep1mt_,       "lep1mt/F"      );
  babyTree_->Branch("mt2w",         &mt2w_,         "mt2w/F"      );
  babyTree_->Branch("bbpt",         &bbpt_,         "bbpt/F"      );
  babyTree_->Branch("wpt",          &wpt_,          "wpt/F"      );
  babyTree_->Branch("bbwdphi",      &bbwdphi_,      "bbwdphi/F"      );
}

//______________________________________________________________________________________
void WHLooper::FillBabyNtuple (const float evtweight)
{

  run_ = stopt.run();
  lumi_ = stopt.lumi();
  event_ = stopt.event();
  leptype_ = stopt.leptype();
  weight_ = evtweight;

  float lep1mt = getMT(stopt.lep1().pt(), stopt.lep1().phi(), stopt.t1metphicorr(), stopt.t1metphicorrphi() );
  LorentzVector bb = myBJets_.at(0) + myBJets_.at(1);
  TVector2 lep(stopt.lep1().px(),stopt.lep1().py());
  TVector2 met;
  met.SetMagPhi(stopt.t1metphicorr(), stopt.t1metphicorrphi());
  TVector2 w = lep+met; 

  pfmet_ = stopt.t1metphicorr();
  lep1mt_ = lep1mt;
  //   mt2w_; // filled in main loop
  bbpt_ = bb.pt();
  wpt_ = w.Mod();
  bbwdphi_ = fabs(TVector2::Phi_mpi_pi(bb.phi() - w.Phi()));

  babyTree_->Fill();
}

//______________________________________________________________________________________
void WHLooper::CloseBabyNtuple ()
{
  babyFile_->cd();
  babyTree_->Write();
  babyFile_->Close();
}

