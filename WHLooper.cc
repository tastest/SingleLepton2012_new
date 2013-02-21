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

std::set<DorkyEventIdentifier> already_seen; 
std::set<DorkyEventIdentifier> events_lasercalib; 
std::set<DorkyEventIdentifier> events_hcallasercalib; 

//--------------------------------------------------------------------

// This is meant to be passed as the third argument, the predicate, of the standard library sort algorithm
inline bool sortByPt(const LorentzVector &vec1, const LorentzVector &vec2 ) {
    return vec1.pt() > vec2.pt();
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

  cout << "[WHLooper::loop] setting up histos" << endl;

  std::map<std::string, TH1F*> h_1d;

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
      if (name.Contains("TChiwh_250_1")) evtweight *= 6.1E-03;
      // to reweight from file - also need to comment stuff before
      //      float vtxweight = vtxweight_n( nvtx, h_vtx_wgt, isData );

      plot1D("h_nvtx_nosel",       stopt.nvtx(),       evtweight, h_1d, 40, 0, 40);
      plot1D("h_vtxweight_nosel", stopt.nvtxweight(), evtweight, h_1d, 41, -4., 4.);

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
      if (!passIsoTrkVeto_v2()) continue;

      // require 2 bjets
      myBJets_ = getBJets(WHLooper::CSVM);
      int nbjets = myBJets_.size();
      if (nbjets < 2) continue;

      // specific cuts: dummy signal region
      LorentzVector bb = myBJets_.at(0) + myBJets_.at(1);
      float lep1mt = getMT(stopt.lep1().pt(), stopt.lep1().phi(), stopt.pfmet(), stopt.pfmetphi() );
      int njets = getNJets();
      if (bb.M() < 95. || bb.M() > 150.) continue;
      if (lep1mt < 120.) continue;
      //      if (nbjets > 2 || njets > 3) continue;
      if (njets > 2) continue;
      if (stopt.pfmet() < 150.) continue;
      if (bb.pt() < 150.) continue;

      ++nEventsPass;

      // fill hists
      plot1D("h_nvtx",       stopt.nvtx(),       evtweight, h_1d, 40, 0, 40);
      plot1D("h_vtxweight", stopt.nvtxweight(), evtweight, h_1d, 41, -4., 4.);

      fillHists1D(h_1d,evtweight);
      if (stopt.leptype() == 0) fillHists1D(h_1d,evtweight,"_e");
      else if (stopt.leptype() == 1) fillHists1D(h_1d,evtweight,"_m");

    } // end event loop

    // delete tree;
    
  } // end file loop
  
    //
    // finish
    //

  savePlots(h_1d, (char*)m_outfilename_.c_str());

  // TFile outfile(m_outfilename_.c_str(),"RECREATE") ; 
  // printf("[WHLooper::loop] Saving histograms to %s\n", m_outfilename_.c_str());
  
  // std::map<std::string, TH1F*>::iterator it1d;
  // for(it1d=h_1d.begin(); it1d!=h_1d.end(); it1d++) {
  //   it1d->second->Write(); 
  //   delete it1d->second;
  // }
  
  // outfile.Write();
  // outfile.Close();

  already_seen.clear();

  gROOT->cd();

  bmark->Stop("benchmark");
  cout << endl;
  cout << nEventsTotal << " Events Processed" << endl;
  cout << nEventsPass << " Events Passed" << endl;
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

  // sort by pt
  sort(bjets.begin()  , bjets.end()  , sortByPt);

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

void WHLooper::fillHists1D(std::map<std::string, TH1F*>& h_1d, const float evtweight, const std::string& suffix) {

  float lep1mt = getMT(stopt.lep1().pt(), stopt.lep1().phi(), stopt.pfmet(), stopt.pfmetphi() );
  int njets = getNJets();
  TVector2 lep(stopt.lep1().px(),stopt.lep1().py());
  TVector2 met;
  met.SetMagPhi(stopt.pfmet(),stopt.pfmetphi());
  TVector2 w = lep+met; 

  plot1D(string("h_lep1pt")+suffix,       stopt.lep1().pt(),       evtweight, h_1d, 1000, 0., 1000.);
  plot1D(string("h_lep1eta")+suffix,      stopt.lep1().eta(),       evtweight, h_1d, 100, -3., 3.);
  plot1D(string("h_lep1mt")+suffix,       lep1mt,       evtweight, h_1d, 1000, 0., 1000.);
  plot1D(string("h_pfmet")+suffix,        stopt.pfmet(),    evtweight, h_1d, 500, 0., 500.);
  plot1D(string("h_njets")+suffix,        njets,              evtweight, h_1d, 10, 0., 10.);
  plot1D(string("h_nbjets")+suffix,       myBJets_.size(),    evtweight, h_1d, 5, 0., 5.);
  plot1D(string("h_wpt")+suffix,          w.Mod(),       evtweight, h_1d, 1000, 0., 1000.);

  // bjets and bbbar plots
  if (myBJets_.size() >= 2) {
    plot1D(string("h_bjet1pt")+suffix,       myBJets_[0].pt(),       evtweight, h_1d, 500, 0., 500.);
    plot1D(string("h_bjet2pt")+suffix,       myBJets_[1].pt(),       evtweight, h_1d, 500, 0., 500.);
    plot1D(string("h_bjet1eta")+suffix,       myBJets_[0].eta(),       evtweight, h_1d, 100, -3., 3.);
    plot1D(string("h_bjet2eta")+suffix,       myBJets_[1].eta(),       evtweight, h_1d, 100, -3., 3.);

    LorentzVector bb = myBJets_.at(0) + myBJets_.at(1);
    plot1D(string("h_bbmass")+suffix,       bb.M(),       evtweight, h_1d, 1000, 0., 1000.);
    plot1D(string("h_bbpt")+suffix,       bb.pt(),       evtweight, h_1d, 500, 0., 500.);
    plot1D(string("h_bbdphi")+suffix,  TVector2::Phi_0_2pi(myBJets_[0].phi() - myBJets_[1].phi()), evtweight, h_1d, 100, 0., 2.*TMath::Pi());
    plot1D(string("h_bbdr")+suffix,  ROOT::Math::VectorUtil::DeltaR( myBJets_.at(0) , myBJets_.at(1) ), evtweight, h_1d, 100, 0., 2.*TMath::Pi());

    plot1D(string("h_bblep1dr")+suffix,  ROOT::Math::VectorUtil::DeltaR( bb , stopt.lep1() ), evtweight, h_1d, 100, 0., 2.*TMath::Pi());

  }

  return;
}
