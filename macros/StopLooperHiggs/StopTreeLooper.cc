#include "StopTreeLooper.h"

//#include "../../CORE/jetSmearingTools.h"
//#include "../../CORE/Thrust.h"
//#include "../../CORE/EventShape.h"
#include "TStopwatch.h"

#include "Math/VectorUtil.h"
#include "../Core/STOPT.h"
#include "../Core/stopUtils.h"
#include "../Plotting/PlotUtilities.h"
#include "../Core/MT2Utility.h"
#include "../Core/mt2bl_bisect.h"
#include "../Core/mt2w_bisect.h"
#include "../Core/PartonCombinatorics.h"
#include "../../Tools/BTagReshaping/BTagReshaping.h"

#include "TROOT.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TChainElement.h"
#include "TFile.h"
#include "TMath.h"
#include "TChain.h"
#include "Riostream.h"
#include "TFitter.h"
#include "TRandom.h"

#include <algorithm>
#include <utility>
#include <map>
#include <set>
#include <list>

using namespace Stop;

std::set<DorkyEventIdentifier> already_seen;
std::set<DorkyEventIdentifier> events_lasercalib;
std::set<DorkyEventIdentifier> events_hcallasercalib;

StopTreeLooper::StopTreeLooper()
{
  m_outfilename_ = "histos.root";
  min_njets = -9999;
  t1metphicorr = -9999.;
  t1metphicorrphi = -9999.;
  t1metphicorrmt = -9999.;
  t1metphicorr_lep = -9999.;
  t1metphicorrphi_lep = -9999.;
  t1metphicorrmt_lep = -9999.;
  dphi_pseudometlep = -9999.;	
  leppt = -9999.;	
  dphimjmin = -9999.;
  dphimj1 = -9999.;
  dphimj2 = -9999.;
  pt_b = -9999.;
  mbb = -9999.;
  dRleptB1 = -9999.;
  htssl = -9999.;
  htosl = -9999.;
  htratiol = -9999.;
  htssm = -9999.;
  htosm = -9999.;
  htratiom = -9999.;
  min_mtpeak = -9999.;
  max_mtpeak = -9999.; 
  n_jets  = -9999;
  n_bjets = -9999;
  n_ljets = -9999;
  chi2min = -9999.;
  chi2minprob = -9999.;
  mt2bmin = -9999.;
  mt2blmin = -9999.;
  mt2wmin = -9999.;
  pfcalo_metratio = -9999.;
  pfcalo_metdphi  = -9999.;
}

StopTreeLooper::~StopTreeLooper()
{
}

void StopTreeLooper::setOutFileName(string filename)
{
  m_outfilename_ = filename;

}

void StopTreeLooper::loop(TChain *chain, TString name)
{

  TStopwatch stwatch;
  //------------------------------
  // check for valid chain
  //------------------------------

  printf("[StopTreeLooper::loop] %s\n", name.Data());

  load_badlaserevents("../Core/badlaser_events.txt", events_lasercalib);
  load_badlaserevents("../Core/badhcallaser_events.txt", events_hcallasercalib);

  TObjArray *listOfFiles = chain->GetListOfFiles();
  TIter fileIter(listOfFiles);
  if (listOfFiles->GetEntries() == 0) {
    cout << "[StopTreeLooper::loop] no files in chain" << endl;
    return;
  }

  //------------------------------------------------------------------------------------------------------
  // set csv discriminator reshaping
  //------------------------------------------------------------------------------------------------------
  
  BTagShapeInterface * nominalShape = new BTagShapeInterface("../../Tools/BTagReshaping/csvdiscr.root", 0.0, 0.0);

  //------------------------------
  // set up histograms
  //------------------------------

  gROOT->cd();

  cout << "[StopTreeLooper::loop] setting up histos" << endl;

  //for signal region 
  std::map<std::string, TH1F*> h_1d_sig;

  //-----------------------------------
  // PU reweighting based on true PU
  //-----------------------------------

  TFile* pu_file = TFile::Open("../vtxreweight/puWeights_Summer12_53x_True_19p5ifb.root");
  if( pu_file == 0 ){
    cout << "vtxreweight error, couldn't open vtx file. Quitting!"<< endl;
    exit(0);
  }

  TH1F* h_pu_wgt = (TH1F*)pu_file->Get("puWeights");
  h_pu_wgt->SetName("h_pu_wgt");

  //------------------------------
  // file loop
  //------------------------------

  unsigned int nEventsChain=0;
  unsigned int nEvents = chain->GetEntries();
  nEventsChain = nEvents;
  ULong64_t nEventsTotal = 0;
  int nevt_check = 0;

  bool isData = name.Contains("data") ? true : false;

  //Define jet multiplicity requirement
  min_njets = 4;
  printf("[StopTreeLooper::loop] N JET min. requirement for signal %i \n", min_njets);
  min_nbjets = 2;
  printf("[StopTreeLooper::loop] N B-JET min. requirement for signal %i \n", min_nbjets);

  //Define peak region 
  min_mtpeak = 50.; 
  max_mtpeak = 80.; 
  printf("[StopTreeLooper::loop] MT PEAK definition %.0f - %.0f GeV \n", min_mtpeak, max_mtpeak);
  min_mbb = 100.; 
  max_mbb = 150.; 
  printf("[StopTreeLooper::loop] MBB PEAK definition %.0f - %.0f GeV \n", min_mbb, max_mbb);

  cout << "[StopTreeLooper::loop] running over chain with total entries " << nEvents << endl;

  stwatch.Start();

  while (TChainElement *currentFile = (TChainElement*)fileIter.Next()) {
  
    //----------------------------
    // load the stop baby tree
    //----------------------------

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
	  stwatch.Stop();
	  if (i_permille%100<0.0001)
	    cout<<"At "<<i_permille/10.<<"% time is "<<stwatch.RealTime()<<" cpu: "<<stwatch.CpuTime()<<endl;
	  stwatch.Start();//stwatch.Continue();
	  
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

      nevt_check++;
      //---------------------------------------------------------------------------- 
      // determine event weight
      // make 2 example histograms of nvtx and corresponding weight
      //---------------------------------------------------------------------------- 

      // to reweight from the nvtx distribution
      // float evtweight = isData ? 1. : 
      // 	( stopt.weight() * 19.5 * stopt.nvtxweight() * stopt.mgcor() );
      float puweight = vtxweight_n( stopt.ntruepu(), h_pu_wgt, isData );
      float evtweight = isData ? 1. : 
	( stopt.weight() * 19.5 * puweight );
      if (!name.Contains("lmg")) evtweight *= stopt.mgcor();

      float nEts=100000;
      float BR=0.60*0.60;
      float lumi= 19.5;
      if(name.Contains("T6tthh_450")) evtweight = ( ( BR * 0.169668 * 1000.0 * lumi ) / nEts);
      if(name.Contains("T6tthh_350")) evtweight = ( ( BR * 0.807323 * 1000.0 * lumi ) / nEts);

      if(name.Contains("T6ttzz_450")) evtweight = ( ( 0.169668 * 1000.0 * lumi ) / (5*nEts));
      if(name.Contains("T6ttzz_350")) evtweight = ( ( 0.807323 * 1000.0 * lumi ) / (5*nEts));

      // plot1D("h_vtx",       stopt.nvtx(), evtweight, h_1d, 40, 0, 40);
      // plot1D("h_vtxweight",     puweight, evtweight, h_1d, 41, -4., 4.);

      //----------------------------------------------------------------------------
      // apply preselection:
      // rho 0-40 GeV, MET filters, >=1 good lepton, veto 2 leptons dR < 0.1
      //----------------------------------------------------------------------------

      if(stopt.indexfirstGoodVertex_()) continue;
      if ( !passEvtSelection(name), false ) continue;

      //----------------------------------------------------------------------------
      // Function to perform MET phi corrections on-the-fly
      // Default branches are: tree->t1metphicorr_ and tree->t1metphicorrmt_
      //----------------------------------------------------------------------------

      pair<float, float> p_t1metphicorr = 
      	getPhiCorrMET( stopt.t1met10(), stopt.t1met10phi(), stopt.nvtx(), !isData);
      t1metphicorr    = p_t1metphicorr.first;
      t1metphicorrphi = p_t1metphicorr.second;
      t1metphicorrmt  = getMT( stopt.lep1().Pt() , stopt.lep1().Phi() , t1metphicorr , t1metphicorrphi );  

      pfcalo_metratio = t1metphicorr/stopt.calomet();
      pfcalo_metdphi  = getdphi(t1metphicorrphi, stopt.calometphi());

      if (t1metphicorr<50) continue;

      //----------------------------------------------------------------------------
      // get jet information
      //----------------------------------------------------------------------------

      jets.clear();
      btag.clear();
      n_jets  = 0;
      n_ljets  = 0;
      n_bjets = 0;

      bool rejectTOBTEC = false;

      for( unsigned int i = 0 ; i < stopt.pfjets().size() ; ++i ){
	
	if( stopt.pfjets().at(i).pt()<30 )  continue;
	if( fabs(stopt.pfjets().at(i).eta())>2.4 )  continue;
	//	if( stopt.pfjets_beta2_0p5().at(i)<0.2 )  continue;
	//	bool passMediumPUid = passMVAJetId(stopt.pfjets().at(i).pt(), stopt.pfjets().at(i).eta(),stopt.pfjets_mvaPUid().at(i),1);
	bool passTightPUid = passMVAJetId(stopt.pfjets().at(i).pt(), stopt.pfjets().at(i).eta(),stopt.pfjets_mva5xPUid().at(i),0);

	if(!passTightPUid) continue;

	if(abs(stopt.pfjets().at(i).eta())>1 && (stopt.pfjets_chm().at(i) - stopt.pfjets_neu().at(i)) > 40 ) rejectTOBTEC = true;

	n_jets++;

	//to not use reshaped discriminator
	//float csv_nominal= stopt.pfjets_csv().at(i);
	float csv_nominal=isData ? stopt.pfjets_csv().at(i)
	  : nominalShape->reshape( stopt.pfjets().at(i).eta(),
				   stopt.pfjets().at(i).pt(),
				   stopt.pfjets_csv().at(i),
				   stopt.pfjets_mcflavorAlgo().at(i) ); 

	btag.push_back( csv_nominal );

	//tighten kinematic requirements for bjets
	if( fabs(stopt.pfjets().at(i).eta())>2.1 )  continue;
	if( stopt.pfjets().at(i).pt()<40 )  continue;

	if (csv_nominal > 0.240 && csv_nominal < 0.679) n_ljets++;

	if (csv_nominal < 0.679) continue;
	  n_bjets++;
	  jets.push_back( stopt.pfjets().at(i) );

      } 

      //------------------------------------------ 
      // baseline selection
      //------------------------------------------ 
      
      if(rejectTOBTEC) continue;
      if (n_jets<min_njets) continue;
      if (n_bjets<min_nbjets) continue;

      //------------------------------------------ 
      // datasets bit
      //------------------------------------------ 
      
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

      float trigweight = isData ? 1. : getsltrigweight(stopt.id1(), stopt.lep1().Pt(), stopt.lep1().Eta());
      float trigweight_dl = isData ? 1. : getdltrigweight(stopt.id1(), stopt.id2());

      //----------------------------------------------------------------------------
      // calculate mbb
      //----------------------------------------------------------------------------      

      mbb = -9999.;
      for (int i=0; i<jets.size(); ++i) 
	for (int j=i+1; j<jets.size(); ++j) {
	  float dR_bbpair = dRbetweenVectors( jets.at(i), jets.at(j) );
	  if ( dR_bbpair > 2./3.*TMath::Pi() ) continue;
	  LorentzVector bbpair = jets.at(i) + jets.at(j);
	  if ( bbpair.M()/bbpair.pt()/dR_bbpair > 0.65 ) continue;
	  if ( fabs(jets.at(i).Rapidity()-jets.at(j).Rapidity()) > 1.2 ) continue;
	  if ( (mbb<0.) || ( fabs(mbb-125.) > fabs(bbpair.M()-125.) ) ) 
	    mbb = bbpair.M();
	}

      //----------------------------------------------------------------------------
      // histogram tags
      //----------------------------------------------------------------------------      
 
      //iso-trk-veto & tau veto
      bool passisotrk = passIsoTrkVeto_v4() && passTauVeto();
 
      //flavor types
      string flav_tag_sl;
      if ( abs(stopt.id1())==13 ) flav_tag_sl = "_muo";
      else if ( abs(stopt.id1())==11 ) flav_tag_sl = "_ele";
      else flav_tag_sl = "_mysterysl";
      string flav_tag_dl;
      if      ( abs(stopt.id1()) == abs(stopt.id2()) && abs(stopt.id1()) == 13 ) 
	flav_tag_dl = "_dimu";
      else if ( abs(stopt.id1()) == abs(stopt.id2()) && abs(stopt.id1()) == 11 ) 
	flav_tag_dl = "_diel";
      else if ( abs(stopt.id1()) != abs(stopt.id2()) && abs(stopt.id1()) == 13 ) 
	flav_tag_dl = "_muel";
      else if ( abs(stopt.id1()) != abs(stopt.id2()) && abs(stopt.id1()) == 11 ) 
	flav_tag_dl = "_elmu";
      else flav_tag_dl = "_mysterydl";
      string basic_flav_tag_dl = flav_tag_dl;
      if ( abs(stopt.id1()) != abs(stopt.id2()) ) basic_flav_tag_dl = "_mueg";

      //
      // SIGNAL REGIONS - single lepton + 3 or >=4 b-tags
      //

      if ( dataset_1l && passSingleLeptonSelection(isData) 
	   && stopt.ngoodlep()==1
	   && passisotrk 
	   //	   && t1metphicorrmt> 120. 
	   && n_bjets>=4 )
	{
	  makeSIGPlots( evtweight*trigweight, h_1d_sig, "_1l4b", flav_tag_sl, 120. );
	  makeSIGPlots( evtweight*trigweight, h_1d_sig, "_1l4b", "", 120. );
	}

      if ( dataset_1l && passSingleLeptonSelection(isData) 
	   && stopt.ngoodlep()==1
	   && passisotrk 
	   && n_jets>=5 
	   //	   && t1metphicorrmt> 150.  
	   && n_bjets==3 
	   && n_ljets>0 )  
	{
	  makeSIGPlots( evtweight*trigweight, h_1d_sig, "_1l3b", flav_tag_sl, 150. );
	  makeSIGPlots( evtweight*trigweight, h_1d_sig, "_1l3b", "", 150. );
	}

      if ( dataset_1l && passSingleLeptonSelection(isData) 
	   && stopt.ngoodlep()==1
	   && passisotrk 
	   && n_jets>=5 
	   //	   && t1metphicorrmt> 150.  
	   && n_bjets==2 )
	{
	  makeSIGPlots( evtweight*trigweight, h_1d_sig, "_1l2b", flav_tag_sl, 150. );
	  makeSIGPlots( evtweight*trigweight, h_1d_sig, "_1l2b", "", 150. );
	}

      //
      // SIGNAL REGIONS - dilepton + 3 or >=4 b-tags
      //

      if ( dataset_2l && passDileptonSelection(isData) 
	   //	   && mbb > min_mbb && mbb < max_mbb 
	   && n_jets>=4 
	   && n_bjets>=4 )
	{
	  makeSIGPlots( evtweight*trigweight_dl, h_1d_sig, "_2l4b", flav_tag_dl, -1. );
	  makeSIGPlots( evtweight*trigweight_dl, h_1d_sig, "_2l4b", "", -1. );
	}

      if ( dataset_2l && passDileptonSelection(isData) 
	   //	   && mbb > min_mbb && mbb < max_mbb
	   && n_jets>=5 
	   && n_bjets==3 )
	{
	  makeSIGPlots( evtweight*trigweight_dl, h_1d_sig, "_2l3b", flav_tag_dl, -1. );
	  makeSIGPlots( evtweight*trigweight_dl, h_1d_sig, "_2l3b", "", -1. );
	}

      if ( dataset_2l && passDileptonSelection(isData) 
	   //	   && mbb > min_mbb && mbb < max_mbb
	   && n_jets>=4 
	   && n_bjets==2 )
	{
	  makeSIGPlots( evtweight*trigweight_dl, h_1d_sig, "_2l2b", flav_tag_dl, -1. );
	  makeSIGPlots( evtweight*trigweight_dl, h_1d_sig, "_2l2b", "", -1. );
	}

    } // end event loop
    
   // delete tree;
    
    stwatch.Stop();
    cout<<"time is "<<stwatch.CpuTime()<<endl;
    stwatch.Start();
    
  } // end file loop
  
    //
    // finish
    //
  
  cout<<"N EVENT CHECK "<<nevt_check<<endl;

  TFile outfile_sig(Form("SIG%s",m_outfilename_.c_str()),"RECREATE") ; 
  printf("[StopTreeLooper::loop] Saving SIG histograms to %s\n", m_outfilename_.c_str());

  std::map<std::string, TH1F*>::iterator it1d_sig;
  for(it1d_sig=h_1d_sig.begin(); it1d_sig!=h_1d_sig.end(); it1d_sig++) {
    it1d_sig->second->Write(); 
    delete it1d_sig->second;
  }

  outfile_sig.Write();
  outfile_sig.Close();

  already_seen.clear();

  gROOT->cd();

}

void StopTreeLooper::makeSIGPlots( float evtweight, std::map<std::string, TH1F*> &h_1d, 
				   string tag_selection, string flav_tag, float mtcut ) 
{

  int nbins = 50;
  float h_xmin = 0.;
  float h_xmax = 500.;
  float x_ovflw = h_xmax-0.001;
  
  //default met
  plot1D("h_sig_met"+tag_selection+flav_tag, min(t1metphicorr, x_ovflw), evtweight, h_1d, nbins-5, 50, h_xmax);
  //lepton pt - enters mT calculation
  plot1D("h_sig_leppt"+tag_selection+flav_tag, min(stopt.lep1().Pt(), x_ovflw), evtweight, h_1d, nbins, h_xmin, h_xmax);
  //angle between lepton and met
  float dphi_metlep = getdphi( stopt.lep1().Phi() , t1metphicorrphi );
  plot1D("h_sig_dphi_metlep"+tag_selection+flav_tag, dphi_metlep, evtweight, h_1d, 15, 0., TMath::Pi());
  //b-pT
  if (jets.size()>0) 
    plot1D("h_sig_bpt1"+tag_selection+flav_tag, jets.at(0).pt(), evtweight, h_1d, 50, 30., 400.);
  if (jets.size()>1) 
    plot1D("h_sig_bpt2"+tag_selection+flav_tag, jets.at(1).pt(), evtweight, h_1d, 50, 30., 400.);

  //mbb distribution
  float mbb_count = -1.;
  if ( mbb <= min_mbb ) mbb_count = 0.5;
  else if ( mbb > min_mbb && mbb < max_mbb ) mbb_count = 1.5;
  else if ( mbb >= max_mbb ) mbb_count = 2.5;
  else mbb_count = -1.;

  float mbbplot = mbb<h_xmin ? h_xmin+0.001 : mbb;
  plot1D("h_sig_mbb"+tag_selection+flav_tag, min(mbbplot, x_ovflw), evtweight, h_1d, nbins, h_xmin, h_xmax);
  plot1D("h_sig_mbb_count"+tag_selection+flav_tag, mbb_count, evtweight, h_1d, 3, 0, 3);
  
  //MT
  float mt_count = -1.;
  if ( t1metphicorrmt > min_mtpeak 
       && t1metphicorrmt < max_mtpeak )    mt_count = 0.5;
  else if ( t1metphicorrmt > mtcut ) mt_count = 1.5;
  if (mtcut<0.) mt_count = 1.5; 
 
  //binning for mT plots
  nbins = 30;
  h_xmin = 0.;
  h_xmax = 300.;
  x_ovflw = h_xmax-0.001;
  plot1D("h_sig_mt"      +tag_selection+flav_tag, min(t1metphicorrmt, x_ovflw), evtweight, h_1d, nbins, h_xmin, h_xmax);
  plot1D("h_sig_mt_count"+tag_selection+flav_tag, mt_count, evtweight, h_1d, 2, 0, 2);
  
  plot1D("h_sig_njets"+tag_selection+flav_tag, min(n_jets,7),  evtweight, h_1d, 7,0,7);
  plot1D("h_sig_nbjets"+tag_selection+flav_tag, min(n_bjets,5), evtweight, h_1d, 5, 0, 5);

}
