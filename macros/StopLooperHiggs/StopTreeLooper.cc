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
  mbb40 = -9999.;
  n_mbb = -9999;
  n_mbb40 = -9999;
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
  n_bjets40 = -9999;
  n_bjets30 = -9999;
  n_ljets = -9999;
  //systematics
  n_bjets40_upBCShape = -9999;
  n_bjets30_upBCShape = -9999;
  n_ljets_upBCShape = -9999;
  n_bjets40_downBCShape = -9999;
  n_bjets30_downBCShape = -9999;
  n_ljets_downBCShape = -9999;
  n_bjets40_upLShape = -9999;
  n_bjets30_upLShape = -9999;
  n_ljets_upLShape = -9999;
  n_bjets40_downLShape = -9999;
  n_bjets30_downLShape = -9999;
  n_ljets_downLShape = -9999;
  chi2min = -9999.;
  chi2minprob = -9999.;
  mt2bmin = -9999.;
  mt2blmin = -9999.;
  mt2wmin = -9999.;
  pfcalo_metratio = -9999.;
  pfcalo_metdphi  = -9999.;
  issigmbb = false;
  issigmbb40 = false;
  issigmt120 = false;
  issigmt150 = false;
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
  //systematic variations for payloads
  BTagShapeInterface * upBCShape    = new BTagShapeInterface("../../Tools/BTagReshaping/csvdiscr.root",  1.0 ,  0.0);
  BTagShapeInterface * downBCShape  = new BTagShapeInterface("../../Tools/BTagReshaping/csvdiscr.root", -1.0 ,  0.0);
  BTagShapeInterface * upLShape     = new BTagShapeInterface("../../Tools/BTagReshaping/csvdiscr.root",  0.0 ,  1.0);
  BTagShapeInterface * downLShape   = new BTagShapeInterface("../../Tools/BTagReshaping/csvdiscr.root",  0.0 , -1.0);

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
  bool isttbar = (name.Contains("ttsl") || name.Contains("ttdl") ||
		  name.Contains("ttall") || name.Contains("tt_")) ? true : false;

  cout<<"[StopTreeLooper::loop] Running on sample "<<name.Data()<<" is ttbar "<<isttbar<<endl;

  //Define jet multiplicity requirement
  min_njets = 4;
  printf("[StopTreeLooper::loop] N JET min. requirement for signal %i \n", min_njets);
  min_nbjets = 0;
  printf("[StopTreeLooper::loop] N B-JET min. requirement for signal %i \n", min_nbjets);

  //Define peak region 
  min_mtpeak = 50.; 
  max_mtpeak = 100.; 
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
      if ( !passEvtSelection(name, false) ) continue;

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

      if (t1metphicorr<50.) continue;

      //----------------------------------------------------------------------------
      // get jet information
      //----------------------------------------------------------------------------

      bjets40.clear();
      bjets30.clear();
      btag.clear();
      n_jets  = 0;
      n_ljets  = 0;
      n_bjets40 = 0;
      n_bjets30 = 0;
      //systematics
      bjets40_upBCShape.clear();
      bjets30_upBCShape.clear();
      n_ljets_upBCShape  = 0;
      n_bjets40_upBCShape = 0;
      n_bjets30_upBCShape = 0;
      bjets40_downBCShape.clear();
      bjets30_downBCShape.clear();
      n_ljets_downBCShape  = 0;
      n_bjets40_downBCShape = 0;
      n_bjets30_downBCShape = 0;
      bjets40_upLShape.clear();
      bjets30_upLShape.clear();
      n_ljets_upLShape  = 0;
      n_bjets40_upLShape = 0;
      n_bjets30_upLShape = 0;
      bjets40_downLShape.clear();
      bjets30_downLShape.clear();
      n_ljets_downLShape  = 0;
      n_bjets40_downLShape = 0;
      n_bjets30_downLShape = 0;

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

	//float csv_upBCShape=(isData || name.Contains("T2")) ? stopt.pfjets_csv().at(i)
	float csv_upBCShape=(isData) ? stopt.pfjets_csv().at(i)
	  : upBCShape->reshape( stopt.pfjets().at(i).eta(),
				stopt.pfjets().at(i).pt(),
				stopt.pfjets_csv().at(i),
				stopt.pfjets_mcflavorAlgo().at(i) );
	
	//float csv_downBCShape=(isData || name.Contains("T2")) ? stopt.pfjets_csv().at(i)
	float csv_downBCShape=(isData) ? stopt.pfjets_csv().at(i)
	  : downBCShape->reshape( stopt.pfjets().at(i).eta(),
				  stopt.pfjets().at(i).pt(),
				  stopt.pfjets_csv().at(i),
				  stopt.pfjets_mcflavorAlgo().at(i) );
	
	//float csv_upLShape=(isData || name.Contains("T2")) ? stopt.pfjets_csv().at(i)
	float csv_upLShape=(isData) ? stopt.pfjets_csv().at(i)
	  : upLShape->reshape( stopt.pfjets().at(i).eta(),
			       stopt.pfjets().at(i).pt(),
			       stopt.pfjets_csv().at(i),
			       stopt.pfjets_mcflavorAlgo().at(i) );
	
	//float csv_downLShape=(isData || name.Contains("T2")) ? stopt.pfjets_csv().at(i)
	float csv_downLShape=(isData) ? stopt.pfjets_csv().at(i)
	  : downLShape->reshape( stopt.pfjets().at(i).eta(),
				 stopt.pfjets().at(i).pt(),
				 stopt.pfjets_csv().at(i),
				 stopt.pfjets_mcflavorAlgo().at(i) );
	
	//tighten kinematic requirements for bjets
	//if( fabs(stopt.pfjets().at(i).eta())>2.1 )  continue;

	if (stopt.pfjets().at(i).pt()>40. && csv_nominal > 0.240 && csv_nominal < 0.679) n_ljets++;
	if (stopt.pfjets().at(i).pt()>40. && csv_upBCShape > 0.240 && csv_upBCShape < 0.679) n_ljets_upBCShape++;
	if (stopt.pfjets().at(i).pt()>40. && csv_downBCShape > 0.240 && csv_downBCShape < 0.679) n_ljets_downBCShape++;
	if (stopt.pfjets().at(i).pt()>40. && csv_upLShape > 0.240 && csv_upLShape < 0.679) n_ljets_upLShape++;
	if (stopt.pfjets().at(i).pt()>40. && csv_downLShape > 0.240 && csv_downLShape < 0.679) n_ljets_downLShape++;

	if (csv_nominal >= 0.679) {
	  n_bjets30++;
	  bjets30.push_back( stopt.pfjets().at(i) );
	  if( stopt.pfjets().at(i).pt()>=40 ) {
	    n_bjets40++;
	    bjets40.push_back( stopt.pfjets().at(i) );
	  }
	}


	if (csv_upBCShape >= 0.679) {
	  n_bjets30_upBCShape++;
	  bjets30_upBCShape.push_back( stopt.pfjets().at(i) );
	  if( stopt.pfjets().at(i).pt()>=40 ) {
	    n_bjets40_upBCShape++;
	    bjets40_upBCShape.push_back( stopt.pfjets().at(i) );
	  }
	}
	if (csv_downBCShape >= 0.679) {
	  n_bjets30_downBCShape++;
	  bjets30_downBCShape.push_back( stopt.pfjets().at(i) );
	  if( stopt.pfjets().at(i).pt()>=40 ) {
	    n_bjets40_downBCShape++;
	    bjets40_downBCShape.push_back( stopt.pfjets().at(i) );
	  }
	}
	if (csv_upLShape >= 0.679) {
	  n_bjets30_upLShape++;
	  bjets30_upLShape.push_back( stopt.pfjets().at(i) );
	  if( stopt.pfjets().at(i).pt()>=40 ) {
	    n_bjets40_upLShape++;
	    bjets40_upLShape.push_back( stopt.pfjets().at(i) );
	  }
	}
	if (csv_downLShape >= 0.679) {
	  n_bjets30_downLShape++;
	  bjets30_downLShape.push_back( stopt.pfjets().at(i) );
	  if( stopt.pfjets().at(i).pt()>=40 ) {
	    n_bjets40_downLShape++;
	    bjets40_downLShape.push_back( stopt.pfjets().at(i) );
	  }
	}

      } 

      //------------------------------------------ 
      // baseline selection
      //------------------------------------------ 
      
      if(rejectTOBTEC) continue;
      if (n_jets<min_njets) continue;
      if (n_bjets30<min_nbjets) continue;

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

      //top pt weight for systematics studies
      float topptweight = isttbar ? sqrt( TopPtWeight_v2(stopt.t().Pt()) * TopPtWeight_v2(stopt.tbar().Pt()) ) : 1.;
      
      //----------------------------------------------------------------------------
      // calculate mbb
      //----------------------------------------------------------------------------      

      n_mbb = 0;
      n_mbb40 = 0;
      mbb = getMbbWithCount(bjets30, n_mbb);
      mbb40 = getMbbWithCount(bjets40, n_mbb40);
      //systematics
      mbb_upBCShape = getMbb(bjets30_upBCShape);
      mbb40_upBCShape = getMbb(bjets40_upBCShape);
      mbb_downBCShape = getMbb(bjets30_downBCShape);
      mbb40_downBCShape = getMbb(bjets40_downBCShape);
      mbb_upLShape = getMbb(bjets30_upLShape);
      mbb40_upLShape = getMbb(bjets40_upLShape);
      mbb_downLShape = getMbb(bjets30_downLShape);
      mbb40_downLShape = getMbb(bjets40_downLShape);

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
      string njets_tag = Form("_%ij", n_jets<7 ? n_jets : 6);

      //
      // SIGNAL REGIONS - dilepton + 3 or >=4 b-tags
      //

      //blind signal regions
      issigmbb = (isData && ( ( n_mbb==1 && mbb >= min_mbb && mbb < max_mbb ) || n_mbb>1 ) ) ? true : false;
      issigmbb40 = (isData && ( ( n_mbb40==1 && mbb40 >= min_mbb && mbb40 < max_mbb ) || n_mbb40>1 ) ) ? true : false;

      if ( dataset_2l && passDileptonSelection(isData) 
	   && n_jets>=4 
	   && n_bjets30>=4 
	   && !issigmbb )
	{
	  makeSIGPlots( evtweight*trigweight_dl, h_1d_sig, "_2l4b", flav_tag_dl, -1. );
	  makeSIGPlots( evtweight*trigweight_dl, h_1d_sig, "_2l4b", "", -1. );
	  makeSIGPlots( evtweight*trigweight_dl*topptweight, h_1d_sig, "_topptwgt_2l4b", "", -1. );
	}

      if ( dataset_2l && passDileptonSelection(isData) 
	   && n_jets>=5 
	   && n_bjets40==3 
	   && n_bjets30==3 
	   && !issigmbb )
	{
	  makeSIGPlots( evtweight*trigweight_dl, h_1d_sig, "_2l3b", flav_tag_dl, -1. );
	  makeSIGPlots( evtweight*trigweight_dl, h_1d_sig, "_2l3b", "", -1. );
	  makeSIGPlots( evtweight*trigweight_dl*topptweight, h_1d_sig, "_topptwgt_2l3b", "", -1. );
	}

      if ( dataset_2l && passDileptonSelection(isData) 
	   && n_jets>=4 
	   && n_bjets30==2) 
	{
	  makeSIGPlots( evtweight*trigweight_dl, h_1d_sig, "_2l2b", flav_tag_dl, -1. );
	  makeSIGPlots( evtweight*trigweight_dl, h_1d_sig, "_2l2b", "", -1. );
	  //jet dependent distributions for systematics studies
	  makeSIGPlots( evtweight*trigweight_dl, h_1d_sig, "_2l2b"+njets_tag, "", -1. );
	  makeSIGPlots( evtweight*trigweight_dl*topptweight, h_1d_sig, "_topptwgt_2l2b", "", -1. );
	  if (n_jets>=5) {
	    makeSIGPlots( evtweight*trigweight_dl, h_1d_sig, "_2l2b_g5j", "", -1. );
	    makeSIGPlots( evtweight*trigweight_dl*topptweight, h_1d_sig, "_topptwgt_2l2b_g5j", "", -1. );
	  }
	}

      //
      // SIGNAL REGIONS - single lepton + 3 or >=4 b-tags
      //

      issigmt120 = (isData && t1metphicorrmt > 120.) ? true : false; 
      issigmt150 = (isData && t1metphicorrmt > 150.) ? true : false; 

      if ( dataset_1l && passSingleLeptonSelection(isData) 
	   && stopt.ngoodlep()==1
	   && passisotrk 
	   && n_jets>=4 )
	{
	  if ( n_bjets30>=4 ) {
	    if ( !issigmt120 ) {
	      makeSIGPlots( evtweight*trigweight, h_1d_sig, "_1l4b", flav_tag_sl, 120. );
	      makeSIGPlots( evtweight*trigweight, h_1d_sig, "_1l4b", "", 120. );
	      if (n_jets>=6) makeSIGPlots( evtweight*trigweight, h_1d_sig, "_1l4b_hinj", "", 120. );
	      makeSIGPlots( evtweight*trigweight*topptweight, h_1d_sig, "_topptwgt_1l4b", "", 120. );
	    }
	    //distributions in MT peak region --> here don't have to worry about signal contamination
	    if ( t1metphicorrmt > min_mtpeak && t1metphicorrmt < max_mtpeak ) {
	      makeSIGPlots( evtweight*trigweight, h_1d_sig, "_mtpeak_1l4b_g4j", "", 120. );
	      if ( n_jets>=5 ) makeSIGPlots( evtweight*trigweight, h_1d_sig, "_mtpeak_1l4b_g5j", "", 150. );
	    }
	  }
	  //systematics
	  if ( n_bjets30_upBCShape>=4 ) {
	    makeMinSIGPlots( evtweight*trigweight, h_1d_sig, "_1l4b_upBCShape", "", 120., mbb_upBCShape );
	    if ( n_jets>=5 ) 
	      makeMinSIGPlots( evtweight*trigweight, h_1d_sig, "_mtpeak_1l4b_g5j_upBCShape", "", 150., mbb_upBCShape );
	  }
	  if ( n_bjets30_downBCShape>=4 ) {
	    makeMinSIGPlots( evtweight*trigweight, h_1d_sig, "_1l4b_downBCShape", "", 120., mbb_downBCShape );
	    if ( n_jets>=5 ) 
	      makeMinSIGPlots( evtweight*trigweight, h_1d_sig, "_mtpeak_1l4b_g5j_downBCShape", "", 150., mbb_downBCShape );
	  }
	  if ( n_bjets30_upLShape>=4 ) {
	    makeMinSIGPlots( evtweight*trigweight, h_1d_sig, "_1l4b_upLShape", "", 120., mbb_upLShape );
	    if ( n_jets>=5 ) 
	      makeMinSIGPlots( evtweight*trigweight, h_1d_sig, "_mtpeak_1l4b_g5j_upLShape", "", 150., mbb_upLShape );
	  }
	  if ( n_bjets30_downLShape>=4 ) {
	    makeMinSIGPlots( evtweight*trigweight, h_1d_sig, "_1l4b_downLShape", "", 120., mbb_downLShape );
	    if ( n_jets>=5 ) 
	      makeMinSIGPlots( evtweight*trigweight, h_1d_sig, "_mtpeak_1l4b_g5j_downLShape", "", 150., mbb_downLShape );
	  }
	}

      if ( dataset_1l && passSingleLeptonSelection(isData) 
	   && stopt.ngoodlep()==1
	   && passisotrk  
	   && n_bjets40==3 
	   && n_bjets30==3 
	   && n_ljets>0 )  
	{
	  if ( !issigmt150 && n_jets>=5 ) {
	    makeSIGPlots( evtweight*trigweight, h_1d_sig, "_1l3b", flav_tag_sl, 150. );
	    makeSIGPlots( evtweight*trigweight, h_1d_sig, "_1l3b", "", 150. );
	    if (n_jets>=7) makeSIGPlots( evtweight*trigweight, h_1d_sig, "_1l3b_hinj", "", 150. );
	    makeSIGPlots( evtweight*trigweight*topptweight, h_1d_sig, "_topptwgt_1l3b", "", 150. );
	  }
	  //distributions in MT peak region --> here don't have to worry about signal contamination
	  if ( t1metphicorrmt > min_mtpeak && t1metphicorrmt < max_mtpeak ) {
	    makeSIGPlots( evtweight*trigweight, h_1d_sig, "_mtpeak_1l3b_g4j", "", 120. );
	    if ( n_jets>=5 ) makeSIGPlots( evtweight*trigweight, h_1d_sig, "_mtpeak_1l3b_g5j", "", 150. );
	  }
	}

      if ( dataset_1l && passSingleLeptonSelection(isData) 
	   && stopt.ngoodlep()==1
	   && passisotrk 
	   && n_jets>=4 
	   && n_bjets30==2)
	{
	  makeSIGPlots( evtweight*trigweight, h_1d_sig, "_1l2b_mt150_g4j", flav_tag_sl, 150. );
	  makeSIGPlots( evtweight*trigweight, h_1d_sig, "_1l2b_mt150_g4j", "", 150. );
	  makeSIGPlots( evtweight*trigweight, h_1d_sig, "_1l2b_mt120_g4j", "", 120. );
	  makeSIGPlots( evtweight*trigweight*topptweight, h_1d_sig, "_topptwgt_1l2b_mt150_g4j", "", 150. );
	  makeSIGPlots( evtweight*trigweight*topptweight, h_1d_sig, "_topptwgt_1l2b_mt120_g4j", "", 120. );
	  //studies of jet dependence
	  makeSIGPlots( evtweight*trigweight, h_1d_sig, "_1l2b_mt150"+njets_tag, "", 150. );
	  makeSIGPlots( evtweight*trigweight, h_1d_sig, "_1l2b_mt120"+njets_tag, "", 120. );
	  //inclusive samples
	  if (n_jets>=5) {
	    makeSIGPlots( evtweight*trigweight, h_1d_sig, "_1l2b_mt150_g5j", "", 150. );
	    makeSIGPlots( evtweight*trigweight, h_1d_sig, "_1l2b_mt120_g5j", "", 120. );
	    makeSIGPlots( evtweight*trigweight*topptweight, h_1d_sig, "_topptwgt_1l2b_mt150_g5j", "", 150. );
	    makeSIGPlots( evtweight*trigweight*topptweight, h_1d_sig, "_topptwgt_1l2b_mt120_g5j", "", 120. );
	    if (n_jets>=6) {
	      makeSIGPlots( evtweight*trigweight, h_1d_sig, "_1l2b_mt150_g6j", "", 150. );
	      makeSIGPlots( evtweight*trigweight, h_1d_sig, "_1l2b_mt120_g6j", "", 120. );
	      makeSIGPlots( evtweight*trigweight*topptweight, h_1d_sig, "_topptwgt_1l2b_mt150_g6j", "", 150. );
	      makeSIGPlots( evtweight*trigweight*topptweight, h_1d_sig, "_topptwgt_1l2b_mt120_g6j", "", 120. );
	    }
	  }
	  //distributions in MT peak region --> here don't have to worry about signal contamination
	  if ( t1metphicorrmt > min_mtpeak && t1metphicorrmt < max_mtpeak ) {
	    makeSIGPlots( evtweight*trigweight, h_1d_sig, "_mtpeak_1l2b_g4j", "", 120. );
	    if (n_jets>=5) makeSIGPlots( evtweight*trigweight, h_1d_sig, "_mtpeak_1l2b_g5j", "", 150. );
	  }
	}

      if ( dataset_1l && passSingleLeptonSelection(isData) 
	   && stopt.ngoodlep()==1
	   && passisotrk 
	   && n_jets>=4 
	   && (n_bjets30==1 || 
	       n_bjets30==2) )
	{
	  makeSIGPlots( evtweight*trigweight, h_1d_sig, "_1l1or2b_mt150_g4j", flav_tag_sl, 150. );
	  makeSIGPlots( evtweight*trigweight, h_1d_sig, "_1l1or2b_mt150_g4j", "", 150. );
	  makeSIGPlots( evtweight*trigweight, h_1d_sig, "_1l1or2b_mt120_g4j", "", 120. );
	  makeSIGPlots( evtweight*trigweight*topptweight, h_1d_sig, "_topptwgt_1l1or2b_mt150_g4j", "", 150. );
	  makeSIGPlots( evtweight*trigweight*topptweight, h_1d_sig, "_topptwgt_1l1or2b_mt120_g4j", "", 120. );
	  //studies of jet dependence
	  makeSIGPlots( evtweight*trigweight, h_1d_sig, "_1l1or2b_mt150"+njets_tag, "", 150. );
	  makeSIGPlots( evtweight*trigweight, h_1d_sig, "_1l1or2b_mt120"+njets_tag, "", 120. );
	  //inclusive samples
	  if (n_jets>=5) {
	    makeSIGPlots( evtweight*trigweight, h_1d_sig, "_1l1or2b_mt150_g5j", "", 150. );
	    makeSIGPlots( evtweight*trigweight, h_1d_sig, "_1l1or2b_mt120_g5j", "", 120. );
	    makeSIGPlots( evtweight*trigweight*topptweight, h_1d_sig, "_topptwgt_1l1or2b_mt150_g5j", "", 150. );
	    makeSIGPlots( evtweight*trigweight*topptweight, h_1d_sig, "_topptwgt_1l1or2b_mt120_g5j", "", 120. );
	    if (n_jets>=6) {
	      makeSIGPlots( evtweight*trigweight, h_1d_sig, "_1l1or2b_mt150_g6j", "", 150. );
	      makeSIGPlots( evtweight*trigweight, h_1d_sig, "_1l1or2b_mt120_g6j", "", 120. );
	      makeSIGPlots( evtweight*trigweight*topptweight, h_1d_sig, "_topptwgt_1l1or2b_mt150_g6j", "", 150. );
	      makeSIGPlots( evtweight*trigweight*topptweight, h_1d_sig, "_topptwgt_1l1or2b_mt120_g6j", "", 120. );
	    }
	  }
	  //distributions in MT peak region --> here don't have to worry about signal contamination
	  if ( t1metphicorrmt > min_mtpeak && t1metphicorrmt < max_mtpeak ) {
	    makeSIGPlots( evtweight*trigweight, h_1d_sig, "_mtpeak_1l1or2b_g4j", "", 120. );
	    if (n_jets>=5) makeSIGPlots( evtweight*trigweight, h_1d_sig, "_mtpeak_1l1or2b_g5j", "", 150. );
	  }
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
				   string tag_selection, string flav_tag, 
				   float mtcut) 
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
  if ( bjets30.size()>0 )
    plot1D("h_sig_bpt1"+tag_selection+flav_tag, bjets30.at(0).pt(), evtweight, h_1d, 50, 30., 400.);
  if ( bjets30.size()>1 ) 
    plot1D("h_sig_bpt2"+tag_selection+flav_tag, bjets30.at(1).pt(), evtweight, h_1d, 50, 30., 400.);
  
  //binning for mbb plots
  nbins = 50;
  h_xmin = 0.;
  h_xmax = 250.;
  x_ovflw = h_xmax-0.001;

  //mbb distribution
  //new definition first bin is CR, second is SR and third is bad pairs
  float mbb_count = -1.;
  if ( n_mbb==1 && ( mbb < min_mbb || mbb >= max_mbb ) ) mbb_count = 0.5;
  else if ( n_mbb==1 && ( mbb >= min_mbb && mbb < max_mbb ) ) mbb_count = 1.5;
  else if ( n_mbb>1 ) mbb_count = 1.5;
  else if ( n_mbb==0 ) mbb_count = 2.5;
  else mbb_count = -1.;

  float mbbplot = mbb<h_xmin ? h_xmin+0.001 : mbb;
  plot1D("h_sig_mbb"+tag_selection+flav_tag, min(mbbplot, x_ovflw), evtweight, h_1d, nbins, h_xmin, h_xmax);
  plot1D("h_sig_mbb_count"+tag_selection+flav_tag, mbb_count, evtweight, h_1d, 3, 0, 3);

  //alternative Mbb sideband region for systematic uncertainty
  //using same distribution, not including combinations failing kinematic requirements
  plot1D("h_sig_mbb_count_alt"+tag_selection+flav_tag, mbb_count, evtweight, h_1d, 3, 0, 3);
  
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

  //alternative MT peak region for systematic uncertainty
  float mt_count_alt = -1.;
  if ( t1metphicorrmt < 80. )        mt_count_alt = 0.5;
  else if ( t1metphicorrmt > mtcut ) mt_count_alt = 1.5;
  if (mtcut<0.) mt_count_alt = 1.5; 
  plot1D("h_sig_mt_count_alt"+tag_selection+flav_tag, mt_count_alt, evtweight, h_1d, 2, 0, 2);

  plot1D("h_sig_njets"+tag_selection+flav_tag, min(n_jets,7),  evtweight, h_1d, 7,0,7);
  plot1D("h_sig_nbjets"+tag_selection+flav_tag, min(n_bjets30,5), evtweight, h_1d, 5, 0, 5);

}



void StopTreeLooper::makeMinSIGPlots( float evtweight, std::map<std::string, TH1F*> &h_1d, 
				      string tag_selection, string flav_tag, 
				      float mtcut, float mbbval) 
{

  //binning for mbb plots
  int nbins = 50;
  float h_xmin = 0.;
  float h_xmax = 250.;
  float x_ovflw = h_xmax-0.001;

  //mbb distribution
  float mbbval_count = -1.;
  if ( mbbval < min_mbb ) mbbval_count = 0.5;
  else if ( mbbval >= min_mbb && mbbval < max_mbb ) mbbval_count = 1.5;
  else if ( mbbval >= max_mbb ) mbbval_count = 2.5;
  else mbbval_count = -1.;

  float mbbvalplot = mbbval<h_xmin ? h_xmin+0.001 : mbbval;
  plot1D("h_sig_mbb"+tag_selection+flav_tag, min(mbbvalplot, x_ovflw), evtweight, h_1d, nbins, h_xmin, h_xmax);
  plot1D("h_sig_mbb_count"+tag_selection+flav_tag, mbbval_count, evtweight, h_1d, 3, 0, 3);

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

}


float StopTreeLooper::getMbb( vector<LorentzVector> &bjets ) 
{
  float mbbval = -9999.;
  for (int i=0; i<bjets.size(); ++i) 
    for (int j=i+1; j<bjets.size(); ++j) {
      float dR_bbpair = dRbetweenVectors( bjets.at(i), bjets.at(j) );
      if ( dR_bbpair > 2./3.*TMath::Pi() ) continue;
      LorentzVector bbpair = bjets.at(i) + bjets.at(j);
      if ( bbpair.M()/bbpair.pt()/dR_bbpair > 0.65 ) continue;
      if ( fabs(bjets.at(i).Rapidity()-bjets.at(j).Rapidity()) > 1.2 ) continue;
      if ( (mbbval<0.) || ( fabs(mbbval-125.) > fabs(bbpair.M()-125.) ) ) 
	mbbval = bbpair.M();
    }
  return mbbval;
}


// updated version that keeps track of the number of pairs
float StopTreeLooper::getMbbWithCount( vector<LorentzVector> &bjets , int &npairs ) 
{
  vector<float> mbb;
  vector<int> i_mbb;
  vector<int> j_mbb;
  npairs = 0;
  for (int i=0; i<(int)bjets.size(); ++i) 
    for (int j=i+1; j<(int)bjets.size(); ++j) {
      float dR_bbpair = dRbetweenVectors( bjets.at(i), bjets.at(j) );
      if ( dR_bbpair > 2./3.*TMath::Pi() ) continue;
      LorentzVector bbpair = bjets.at(i) + bjets.at(j);
      if ( bbpair.M()/bbpair.pt()/dR_bbpair > 0.65 ) continue;
      if ( fabs(bjets.at(i).Rapidity()-bjets.at(j).Rapidity()) > 1.2 ) continue;
      //store all pairs
      mbb.push_back(bbpair.M());
      i_mbb.push_back(i);
      j_mbb.push_back(j);
    }

  float mbbval = -9999.;
  int i_mbbval = -9;
  int j_mbbval = -9;
  //find the pair closest to 125
  for (int k=0; k<(int)mbb.size(); ++k) {
    if ( (mbbval<0.) || ( fabs(mbbval-125.) > fabs(mbb.at(k)-125.) ) ) {
      mbbval = mbb.at(k);
      i_mbbval = i_mbb.at(k);
      j_mbbval = j_mbb.at(k);
    }
  } 
  if (mbbval>0.) npairs++;
 
  //now count additional pairs
  for (int k=0; k<(int)mbb.size(); ++k) {
    //do not use jets in existing pairing
    if ( i_mbb.at(k)==i_mbbval || 
	 j_mbb.at(k)==i_mbbval ||
	 i_mbb.at(k)==j_mbbval || 
	 j_mbb.at(k)==j_mbbval ) continue;
    npairs++;
  } 
  if (npairs>2) npairs = 2;

  return mbbval;
}

