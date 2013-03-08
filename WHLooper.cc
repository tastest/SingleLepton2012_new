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

// regions to do
const bool blindSignal = true;
const bool doSignal = true;
const bool doCR1 = true;
const bool doStopSel = true;

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

  gROOT->cd();

  cout << "[WHLooper::loop] creating output file: " << m_outfilename_ << endl;

  outfile_ = new TFile(m_outfilename_.c_str(),"RECREATE") ; 

  cout << "[WHLooper::loop] setting up histos" << endl;

  // signal region hists
  std::map<std::string, TH1F*> h_1d_sig_presel, h_1d_sig_final;
  // signal region nm1 hists
  std::map<std::string, TH1F*> h_1d_sig_nbjets_nm1, h_1d_sig_met_nm1, h_1d_sig_mt_nm1, h_1d_sig_mt2bl_nm1;

  // cr1 hists
  std::map<std::string, TH1F*> h_1d_cr1_presel, h_1d_cr1_final;
  // cr1 nm1 hists
  std::map<std::string, TH1F*> h_1d_cr1_nbjets_nm1, h_1d_cr1_met_nm1, h_1d_cr1_mt_nm1, h_1d_cr1_mt2bl_nm1;

  // stop region hists
  std::map<std::string, TH1F*> h_1d_stop_presel, h_1d_stop_comp;


  if (doSignal) {
    outfile_->mkdir("sig_presel");
    if (doNM1Plots) {
      outfile_->mkdir("sig_nbjets_nm1");
      outfile_->mkdir("sig_met_nm1");
      outfile_->mkdir("sig_mt_nm1");
      outfile_->mkdir("sig_mt2bl_nm1");
    }
    outfile_->mkdir("sig_final");
  }

  if (doCR1) {
    outfile_->mkdir("cr1_presel");
    if (doNM1Plots) {
      outfile_->mkdir("cr1_nbjets_nm1");
      outfile_->mkdir("cr1_met_nm1");
      outfile_->mkdir("cr1_mt_nm1");
      outfile_->mkdir("cr1_mt2bl_nm1");
    }
    outfile_->mkdir("cr1_final");
  }

  if (doStopSel) {
    outfile_->mkdir("stop_presel");
    outfile_->mkdir("stop_comp");
  }

  outfile_->cd();

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
	float xsecweight = 1.;
	if (name.Contains("TChiwh_150")) xsecweight = 2.4 * 0.5 * weight_lumi_br_nevents; // now 20k events
	else if (name.Contains("TChiwh_200")) xsecweight = 0.79 * 0.5 * weight_lumi_br_nevents; // now 20k events
	else if (name.Contains("TChiwh_250")) xsecweight = 0.32 * weight_lumi_br_nevents;
	else if (name.Contains("TChiwh_300")) xsecweight = 0.15 * weight_lumi_br_nevents;
	else if (name.Contains("TChiwh_350")) xsecweight = 0.074 * weight_lumi_br_nevents;

	// reset weight here for signal MC. Ignore weight, nvtxweight, mgcor
	evtweight = xsecweight * 19.5;
      }
      // to reweight from file - also need to comment stuff before
      //      float vtxweight = vtxweight_n( nvtx, h_vtx_wgt, isData );

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
      bjets_.clear();
      jets_csv_.clear();
      jets_idx_.clear();
      bjets_idx_.clear();
      njets_ = 0;
      njetsalleta_ = 0;
      nbjets_ = 0;

      for( unsigned int i = 0 ; i < stopt.pfjets().size() ; ++i ){
	
	// basic jet selection
	if( stopt.pfjets().at(i).pt()<30 )  continue;
	if( fabs(stopt.pfjets().at(i).eta())>4.7 )  continue;
	if ( (fabs(stopt.pfjets().at(i).eta()) < 2.5) 
	     && (stopt.pfjets_beta2_0p5().at(i)<0.2) ) continue;
	
	++njetsalleta_;
	if (fabs(stopt.pfjets().at(i).eta()) <= 2.4) ++njets_;

      	jets_.push_back( stopt.pfjets().at(i) );
	jets_idx_.push_back(i);
	// if (n_jets==1) 
	//   dphimj1 = getdphi(t1metphicorrphi, stopt.pfjets().at(i).phi() );
	// if (n_jets==2) {
	//   dphimj2 = getdphi(t1metphicorrphi, stopt.pfjets().at(i).phi() );
	//   dphimjmin = TMath::Min( dphimj1 , dphimj2 );
	// }

	float csv_nominal= stopt.pfjets_csv().at(i);

	//RESHAPING -- TO DO WITH UPDATED BABIES
	//only reshape for b jets ---> use status 3 matching information
	//float csv_nominal=nominalShape->reshape(stopt.pfjets().at(i).Eta(),stopt.pfjets().at(i).Pt(),stopt.pfjets_csv().at(i),(stopt.pfjets_flav().at(i)==5?5:0)); 
	//treat anything not matched to b or c as light
	// float csv_nominal=nominalShape->reshape( stopt.pfjets().at(i).eta(),
	// 					 stopt.pfjets().at(i).pt(),
	// 					 stopt.pfjets_csv().at(i),
	// 					( stopt.pfjets_mcflavorAlgo().at(i)>3 ? stopt.pfjets_mcflavorAlgo().at(i) : 1 ) ); 

	jets_csv_.push_back( csv_nominal );

	if( (fabs(stopt.pfjets().at(i).eta()) <= 2.4) && (csv_nominal > getCSVCut(WHLooper::CSVM)) ) {
	  bjets_.push_back( stopt.pfjets().at(i) );
	  bjets_idx_.push_back(i);
	  ++nbjets_;
	}

	//      	sigma_jets.push_back(stopt.pfjets_sigma().at(i));

        // //count jets that are not overlapping with second lepton
	// if (isData) continue;
	// if (stopt.nleps()!=2) continue;
	// if (stopt.mclep2().pt() < 30.) continue;
	// if (ROOT::Math::VectorUtil::DeltaR(stopt.mclep2(), stopt.pfjets().at(i)) > 0.4 ) continue;
	// n_ljets--;

      } // loop over pfjets

      // calculate mt2 vars after selecting jets -- require at least 2 bjets
      LorentzVector bb;
      mt2b_ = -1.;
      mt2bl_ = -1.;
      mt2w_ = -1.;
      if (nbjets_ >= 2) {
	bb = bjets_.at(0) + bjets_.at(1);
	mt2b_ = calculateMT2w(jets_, jets_csv_, stopt.lep1(), stopt.t1metphicorr(), stopt.t1metphicorrphi(), MT2b);
	mt2bl_ = calculateMT2w(jets_, jets_csv_, stopt.lep1(), stopt.t1metphicorr(), stopt.t1metphicorrphi(), MT2bl);
	mt2w_ = calculateMT2w(jets_, jets_csv_, stopt.lep1(), stopt.t1metphicorr(), stopt.t1metphicorrphi(), MT2w);
      }

      // TVector2 lep(stopt.lep1().px(),stopt.lep1().py());
      // TVector2 met;
      // met.SetMagPhi(stopt.t1metphicorr(), stopt.t1metphicorrphi());
      // TVector2 w = lep+met; 

      // end variables --------------------------------------
      // ----------------------------------------------------

      // ----------------------------------------------------
      // selections bits

      // should replace string comps with enums..

      bool dataset_1l=false;

      if((isData) && name.Contains("muo") 
	 && (abs(stopt.id1()) == 13 ))  dataset_1l=true;
      if((isData) && name.Contains("ele") 
	 && (abs(stopt.id1()) == 11 ))  dataset_1l=true;

      if(!isData) dataset_1l=true;

      bool dataset_2l=false;

      if((isData) && name.Contains("dimu") 
	 && (abs(stopt.id1()) == 13 ) 
	 && (abs(stopt.id2())==13)) dataset_2l=true;
      if((isData) && name.Contains("diel") 
	 && (abs(stopt.id1()) == 11 ) 
	 && (abs(stopt.id2())==11)) dataset_2l=true;
      if((isData) && name.Contains("mueg") 
	 && abs(stopt.id1()) != abs(stopt.id2())) 
	dataset_2l=true;

      if(!isData) dataset_2l=true;

      //      bool passisotrk = passIsoTrkVeto_v2();
      bool passisotrk = passIsoTrkVeto_v3() && (stopt.ngoodlep() == 1);
      //      bool passisotrk = passIsoTrkVeto_v4() && passTauVeto();

      // end selections bits --------------------------------
      // ----------------------------------------------------

      // ----------------------------------------------------
      // region selection and plots

      // always require at least 2 (central) jets
      if (njets_ < 2) continue;

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
	   && (bb.M() > 100.) && (bb.M() < 140.)
	   && (stopt.t1metphicorr() > 50.) 
	   && (!isData || !blindSignal) ) {

        fillHists1DWrapper(h_1d_sig_presel,evtweight1l,"sig_presel");

	bool fail = false;
	if ( !fail && (njetsalleta_ == 2) ) {
	  if (doNM1Plots) fillHists1DWrapper(h_1d_sig_nbjets_nm1,evtweight1l,"sig_nbjets_nm1");
	} 
	else fail = true;

	if (!fail && (nbjets_ == 2) ) {
	  if (doNM1Plots) fillHists1DWrapper(h_1d_sig_met_nm1,evtweight1l,"sig_met_nm1");
	}
	else fail = true;

	if (!fail && (stopt.t1metphicorr() > 175.) ) {
	  if (doNM1Plots) fillHists1DWrapper(h_1d_sig_mt_nm1,evtweight1l,"sig_mt_nm1");
	}
	else fail = true;

	if (!fail && (stopt.t1metphicorrmt() > 120.) ) {
	  if (doNM1Plots) fillHists1DWrapper(h_1d_sig_mt2bl_nm1,evtweight1l,"sig_mt2bl_nm1");
	}
	else fail = true;

	if (!fail && (mt2bl_ > 200.) ) {
	  fillHists1DWrapper(h_1d_sig_final,evtweight1l,"sig_final");
	}

      } // signal region sel

      // -------------------------------------------
      // *** CR1: m(bb) > 150 (can also put upper bound)
      //  otherwise same as signal region

      if ( doStopSel
	   && passSingleLeptonSelection(isData) 
	   && passisotrk 
	   && (nbjets_ >= 2)
	   && (bb.M() > 150.)
	   && (stopt.t1metphicorr() > 50.) ) {

        fillHists1DWrapper(h_1d_cr1_presel,evtweight1l,"cr1_presel");

	bool fail = false;
	if ( !fail && (njetsalleta_ == 2) ) {
	  if (doNM1Plots) fillHists1DWrapper(h_1d_cr1_nbjets_nm1,evtweight1l,"cr1_nbjets_nm1");
	} 
	else fail = true;

	if (!fail && (nbjets_ == 2) ) {
	  if (doNM1Plots) fillHists1DWrapper(h_1d_cr1_met_nm1,evtweight1l,"cr1_met_nm1");
	}
	else fail = true;

	if (!fail && (stopt.t1metphicorr() > 175.) ) {
	  if (doNM1Plots) fillHists1DWrapper(h_1d_cr1_mt_nm1,evtweight1l,"cr1_mt_nm1");
	}
	else fail = true;

	if (!fail && (stopt.t1metphicorrmt() > 120.) ) {
	  if (doNM1Plots) fillHists1DWrapper(h_1d_cr1_mt2bl_nm1,evtweight1l,"cr1_mt2bl_nm1");
	}
	else fail = true;

	if (!fail && (mt2bl_ > 200.) ) {
	  fillHists1DWrapper(h_1d_cr1_final,evtweight1l,"cr1_final");
	}

      } // cr1 region sel


      // -------------------------------------------
      // *** Stop Presel region
      //   >= 1 lepton
      //   iso track veto v2
      //   >= 4 jets
      //   >= 1 bjet
      //   met > 100
      //   blind data (includes stop signal region)

      // *** Stop Presel comparison region
      //   >= 1 lepton
      //   iso track veto v2
      //   >= 4 jets
      //   >= 1 bjet
      //   met > 150
      //   mt > 120
      //   blind data (includes stop signal region)

      if ( doStopSel
	   && passSingleLeptonSelection(isData) 
	   && passIsoTrkVeto_v2()
	   && (njets_ >= 4)
	   && (nbjets_ >= 1)
	   && (stopt.t1metphicorr() > 100.) 
	   && !isData ) {

        fillHists1DWrapper(h_1d_stop_presel,evtweight1l,"stop_presel");

	if ( (stopt.t1metphicorr() > 150.) && (stopt.t1metphicorrmt() > 120.) ) {
	  fillHists1DWrapper(h_1d_stop_comp,evtweight1l,"stop_comp");
	} 

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
    savePlotsDir(h_1d_sig_presel,outfile_,"sig_presel");
    if (doNM1Plots) {
      savePlotsDir(h_1d_sig_nbjets_nm1,outfile_,"sig_nbjets_nm1");
      savePlotsDir(h_1d_sig_met_nm1,outfile_,"sig_met_nm1");
      savePlotsDir(h_1d_sig_mt_nm1,outfile_,"sig_mt_nm1");
      savePlotsDir(h_1d_sig_mt2bl_nm1,outfile_,"sig_mt2bl_nm1");
    }
    savePlotsDir(h_1d_sig_final,outfile_,"sig_final");
  }

  if (doCR1) {
    savePlotsDir(h_1d_cr1_presel,outfile_,"cr1_presel");
    if (doNM1Plots) {
      savePlotsDir(h_1d_cr1_nbjets_nm1,outfile_,"cr1_nbjets_nm1");
      savePlotsDir(h_1d_cr1_met_nm1,outfile_,"cr1_met_nm1");
      savePlotsDir(h_1d_cr1_mt_nm1,outfile_,"cr1_mt_nm1");
      savePlotsDir(h_1d_cr1_mt2bl_nm1,outfile_,"cr1_mt2bl_nm1");
    }
    savePlotsDir(h_1d_cr1_final,outfile_,"cr1_final");
  }

  if (doStopSel) {
    savePlotsDir(h_1d_stop_presel,outfile_,"stop_presel");
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
    if (stopt.leptype() == 0) fillHists1D(h_1d,evtweight,dir,"_e");
    else if (stopt.leptype() == 1) fillHists1D(h_1d,evtweight,dir,"_m");
  }
}

//--------------------------------------------------------------------

void WHLooper::fillHists1D(std::map<std::string, TH1F*>& h_1d, const float evtweight, const std::string& dir, const std::string& suffix) {

  outfile_->cd(dir.c_str());

  TVector2 lep(stopt.lep1().px(),stopt.lep1().py());
  TVector2 met;
  met.SetMagPhi(stopt.t1metphicorr(), stopt.t1metphicorrphi());
  TVector2 w = lep+met; 

  plot1D("h_lep1pt"+suffix,       stopt.lep1().pt(),       evtweight, h_1d, 1000, 0., 1000.);
  plot1D("h_lep1eta"+suffix,      stopt.lep1().eta(),       evtweight, h_1d, 100, -3., 3.);
  plot1D("h_lep1mt"+suffix,       stopt.t1metphicorrmt(),       evtweight, h_1d, 1000, 0., 1000.);
  plot1D("h_pfmet"+suffix,        stopt.t1metphicorr(),    evtweight, h_1d, 500, 0., 500.);
  // plot1D("h_pfsumet"+suffix,      stopt.pfsumet(),    evtweight, h_1d, 1500, 0., 1500.);
  // plot1D("h_pfmetsig"+suffix,     stopt.pfmet()/sqrt(stopt.pfsumet()),   evtweight, h_1d, 500, 0., 20.);
  plot1D("h_njets"+suffix,        njets_,              evtweight, h_1d, 10, 0., 10.);
  plot1D("h_njetsalleta"+suffix,  njetsalleta_,        evtweight, h_1d, 10, 0., 10.);
  plot1D("h_nbjets"+suffix,       nbjets_,    evtweight, h_1d, 5, 0., 5.);
  plot1D("h_wpt"+suffix,          w.Mod(),       evtweight, h_1d, 1000, 0., 1000.);
  plot1D("h_lep1metdphi"+suffix,  fabs(TVector2::Phi_mpi_pi(stopt.lep1().phi() - stopt.t1metphicorrphi())),  evtweight, h_1d, 50, 0., TMath::Pi());

  // phi cor met validation
  // plot1D("h_metdiff"+suffix,        stopt.t1metphicorr() - stopt.pfmet(),    evtweight, h_1d, 500, -250., 250.);
  // plot1D("h_metphidiff"+suffix,  fabs(TVector2::Phi_mpi_pi(stopt.t1metphicorrphi() - stopt.pfmetphi())),    evtweight, h_1d,  50, 0., TMath::Pi());


  if (isWjets_) {
    plot1D("h_nbs",       stopt.nbs(),       evtweight, h_1d, 5, 0, 5);
  }

  // bjets and bbbar plots
  if (nbjets_ >= 2) {
    plot1D("h_bjet1pt"+suffix,       bjets_[0].pt(),       evtweight, h_1d, 500, 0., 500.);
    plot1D("h_bjet2pt"+suffix,       bjets_[1].pt(),       evtweight, h_1d, 500, 0., 500.);
    plot1D("h_bjet1eta"+suffix,       bjets_[0].eta(),       evtweight, h_1d, 100, -3., 3.);
    plot1D("h_bjet2eta"+suffix,       bjets_[1].eta(),       evtweight, h_1d, 100, -3., 3.);

    LorentzVector bb = bjets_.at(0) + bjets_.at(1);
    plot1D("h_bbmass"+suffix,       bb.M(),       evtweight, h_1d, 1000, 0., 1000.);
    plot1D("h_bbpt"+suffix,       bb.pt(),       evtweight, h_1d, 500, 0., 500.);
    // plot1D("h_bbdr"+suffix,  ROOT::Math::VectorUtil::DeltaR( bjets_.at(0) , bjets_.at(1) ), evtweight, h_1d, 100, 0., 2.*TMath::Pi());
    // plot1D("h_bbdphi"+suffix,  fabs(TVector2::Phi_mpi_pi(bjets_[0].phi() - bjets_[1].phi())), evtweight, h_1d, 50, 0., TMath::Pi());

    // plot1D("h_bblep1dr"+suffix,  ROOT::Math::VectorUtil::DeltaR( bb , stopt.lep1() ), evtweight, h_1d, 100, 0., 2.*TMath::Pi());
    // plot1D("h_bblep1dphi"+suffix,  fabs(TVector2::Phi_mpi_pi(bb.phi() - stopt.lep1().phi())), evtweight, h_1d, 50, 0., TMath::Pi());

    plot1D("h_bbwdphi"+suffix,  fabs(TVector2::Phi_mpi_pi(bb.phi() - w.Phi())), evtweight, h_1d, 50, 0., TMath::Pi());
    // plot1D("h_bbwdpt"+suffix,   bb.pt() - w.Mod(),       evtweight, h_1d, 500, -250., 250.);

    // plot1D("h_bbwsumpt"+suffix,       bb.pt()+w.Mod(),       evtweight, h_1d, 1000, 0., 1000.);

    // plot1D("h_allsumpt"+suffix,      bjets_[0].pt()+ bjets_[1].pt()+stopt.lep1().pt()+stopt.t1metphicorr() , evtweight, h_1d, 1500, 0., 1500.);

    // LorentzVector b1lep1 = bjets_.at(0) + stopt.lep1();
    // plot1D("h_bjet1lep1mass"+suffix,       b1lep1.M(),       evtweight, h_1d, 1000, 0., 1000.);
    // float bjet1lep1dphi = fabs(TVector2::Phi_mpi_pi(bjets_[0].phi() - stopt.lep1().phi()));
    // plot1D("h_bjet1lep1dphi"+suffix,  bjet1lep1dphi, evtweight, h_1d, 50, 0., TMath::Pi());
    // LorentzVector b2lep1 = bjets_.at(1) + stopt.lep1();
    // plot1D("h_bjet2lep1mass"+suffix,       b2lep1.M(),       evtweight, h_1d, 1000, 0., 1000.);
    // float bjet2lep1dphi = fabs(TVector2::Phi_mpi_pi(bjets_[1].phi() - stopt.lep1().phi()));
    // plot1D("h_bjet2lep1dphi"+suffix,  bjet2lep1dphi, evtweight, h_1d, 50, 0., TMath::Pi());

    // plot1D("h_bjetlep1mindphi"+suffix,  TMath::Min(bjet1lep1dphi,bjet2lep1dphi), evtweight, h_1d, 50, 0., TMath::Pi());

    plot1D("h_mt2b"+suffix,   mt2b_,  evtweight, h_1d, 1000, 0., 1000.);
    plot1D("h_mt2bl"+suffix,  mt2bl_, evtweight, h_1d, 1000, 0., 1000.);
    plot1D("h_mt2w"+suffix,   mt2w_,  evtweight, h_1d, 1000, 0., 1000.);

    // use loose btags here in case i plot before requiring 2 med
    // std::vector<int> bjetIdx = getBJetIndex(WHLooper::CSVL,-1,-1);
    // plot1D("h_bjet1mc3"+suffix, stopt.pfjets_mc3().at(bjetIdx.at(0)) , evtweight, h_1d, 40, -20., 20.);
    // plot1D("h_bjet2mc3"+suffix, stopt.pfjets_mc3().at(bjetIdx.at(1)) , evtweight, h_1d, 40, -20., 20.);

    // need V00-02-20 or higher babies for these vars
    // plot1D("h_bjet1flavor"+suffix, abs(stopt.pfjets_mcflavorAlgo().at(bjetIdx.at(0))) , evtweight, h_1d, 23, -1., 22.);
    // plot1D("h_bjet2flavor"+suffix, abs(stopt.pfjets_mcflavorAlgo().at(bjetIdx.at(1))) , evtweight, h_1d, 23, -1., 22.);
  } // if nbjets >= 2

  // if (njets >= 2) {
  //   std::vector<int> jetIdx = getJetIndex();
  //   // need V00-02-20 or higher babies for these vars
  //   plot1D("h_jet1flavor"+suffix, abs(stopt.pfjets_mcflavorAlgo().at(jetIdx.at(0))) , evtweight, h_1d, 23, -1., 22.);
  //   plot1D("h_jet2flavor"+suffix, abs(stopt.pfjets_mcflavorAlgo().at(jetIdx.at(1))) , evtweight, h_1d, 23, -1., 22.);

  //   if (stopt.nbs() == 0) {
  //     plot1D("h_jet1flavor_nobs"+suffix, abs(stopt.pfjets_mcflavorAlgo().at(jetIdx.at(0))) , evtweight, h_1d, 23, -1., 22.);
  //     plot1D("h_jet2flavor_nobs"+suffix, abs(stopt.pfjets_mcflavorAlgo().at(jetIdx.at(1))) , evtweight, h_1d, 23, -1., 22.);
  //   } else {
  //     plot1D("h_jet1flavor_bs"+suffix, abs(stopt.pfjets_mcflavorAlgo().at(jetIdx.at(0))) , evtweight, h_1d, 23, -1., 22.);
  //     plot1D("h_jet2flavor_bs"+suffix, abs(stopt.pfjets_mcflavorAlgo().at(jetIdx.at(1))) , evtweight, h_1d, 23, -1., 22.);
  //   }

  // }

  plot1D("h_nvtx",      stopt.nvtx(),       evtweight, h_1d, 40, 0, 40);
  plot1D("h_vtxweight", stopt.nvtxweight(), evtweight, h_1d, 41, -4., 4.);

  return;
}
