#include "StopTreeLooper.h"

//#include "../../CORE/jetSmearingTools.h"
//#include "../../CORE/Thrust.h"
//#include "../../CORE/EventShape.h"

#include "../Core/STOPT.h"
#include "../Core/stopUtils.h"
#include "../Plotting/PlotUtilities.h"
//#include "../Core/MT2Utility.h"
#include "../Core/MT2.h"
#include "../Core/mt2bl_bisect.h"
#include "../Core/mt2w_bisect.h"
//#include "../Core/PartonCombinatorics.h"
#include "../../Tools/BTagReshaping/BTagReshaping.h"
#include "../../CORE/utilities.h"

#include "TROOT.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TChainElement.h"
#include "TFile.h"
#include "TMath.h"
#include "TChain.h"
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

StopTreeLooper::StopTreeLooper()
{
  m_outfilename_ = "histos.root";
  t1metphicorr = -9999.;
  t1metphicorrphi = -9999.;
  t1metphicorrmt = -9999.;
  min_mtpeak = -9999.;
  max_mtpeak = -9999.; 
}

StopTreeLooper::~StopTreeLooper()
{
}

void StopTreeLooper::setOutFileName(string filename)
{
  m_outfilename_ = filename;

}

//--------------------------------------------------------------------

void StopTreeLooper::loop(TChain *chain, TString name)
{

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

  //------------------------------
  // set up histograms
  //------------------------------

  gROOT->cd();

  cout << "[StopTreeLooper::loop] setting up histos" << endl;

  std::map<std::string, TH1F*> h_1d;
  std::map<std::string, TH2F*> h_2d;
  std::map<std::string, TH1F*> h_1d_cr1;
  std::map<std::string, TH1F*> h_1d_cr2;
  std::map<std::string, TH1F*> h_1d_cr4;
  std::map<std::string, TH1F*> h_1d_cr3;
  std::map<std::string, TH1F*> h_1d_cr5;
  std::map<std::string, TH1F*> h_1d_qg;
  std::map<std::string, TH2F*> h_2d_qg;

  //------------------------------------------------------------------------------------------------------
  // set csv discriminator reshaping
  //------------------------------------------------------------------------------------------------------
  
  BTagShapeInterface * nominalShape = new BTagShapeInterface("../../Tools/BTagReshaping/csvdiscr.root", 0.0, 0.0);

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
  //  ULong64_t i_permille_old = 0;

  bool isData = name.Contains("data") ? true : false;

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


      //------------------------------------------                                                                                                                             
      // datasets bit                                                                                                                                                          
      //------------------------------------------                                                                                                                             

      bool dataset_1l=false;

      if((isData) && name.Contains("muo") && (abs(stopt.id1()) == 13 ))  dataset_1l=true;
      if((isData) && name.Contains("ele") && (abs(stopt.id1()) == 11 ))  dataset_1l=true;

      if(!isData) dataset_1l=true;

      bool dataset_CR4=false;

      if((isData) && name.Contains("dimu") && (abs(stopt.id1()) == 13 ) && (abs(stopt.id2())==13) && (fabs( stopt.dilmass() - 91.) > 15.)) dataset_CR4=true;
      if((isData) && name.Contains("diel") && (abs(stopt.id1()) == 11 ) && (abs(stopt.id2())==11) && (fabs( stopt.dilmass() - 91.) > 15.)) dataset_CR4=true;
      if((isData) && name.Contains("mueg") && abs(stopt.id1()) != abs(stopt.id2())) dataset_CR4=true;

      if(!isData) dataset_CR4=true;

      // remove the runs that are bad for the pixel alignement
      //      if(stopt.run()>=207883 && stopt.run()<=208307) continue;

      //---------------------------------------------------------------------------- 
      // determine event weight
      // make 2 example histograms of nvtx and corresponding weight
      //---------------------------------------------------------------------------- 

      float puweight = vtxweight_n( stopt.ntruepu(), h_pu_wgt, isData );
      float evtweight = isData ? 1. : ( stopt.weight() * 19.5 * puweight);

      // to reweight from file - also need to comment stuff before
      //      float vtxweight = vtxweight_n( nvtx, h_vtx_wgt, isData );

      plot1D("h_vtx",       stopt.nvtx(),       evtweight, h_1d, 40, 0, 40);
      plot1D("h_vtxweight", stopt.nvtxweight(), evtweight, h_1d, 41, -4., 4.);

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

      jetsP4.clear();
      csv.clear();
      sigma.clear();
      mc3.clear();
      mc1.clear();
      lrm.clear();

      int n_ljets=0;
      int n_taujets=0;

      for( unsigned int i = 0 ; i < stopt.pfjets().size() ; ++i ){

        if( stopt.pfjets().at(i).pt()<30 )  continue;
        if( fabs(stopt.pfjets().at(i).eta())>2.4 )  continue;
	if( stopt.pfjets_beta2_0p5().at(i) < 0.2) continue;

	
	float csv_nominal=isData ? stopt.pfjets_csv().at(i)
	  : nominalShape->reshape( stopt.pfjets().at(i).eta(),
				   stopt.pfjets().at(i).pt(),
				   stopt.pfjets_csv().at(i),
				   stopt.pfjets_mcflavorAlgo().at(i) ); 

        // fill again because of the beta and csv                                                                                                                      
        jetsP4.push_back( stopt.pfjets().at(i) );
        csv.push_back( csv_nominal );
        mc3.push_back( stopt.pfjets_mc3().at(i) );
        mc1.push_back( stopt.pfjets_mcflavorAlgo().at(i) );
        sigma.push_back( stopt.pfjets_sigma().at(i) );
        lrm.push_back( stopt.pfjets_lrm().at(i));

	n_ljets++;
	n_taujets++;
	
	if(stopt.id2()!=(-1)) {
	  if(deltaR(stopt.pfjets().at(i).eta() , stopt.pfjets().at(i).phi() , stopt.lep2().eta(), stopt.lep2().phi())>0.4) continue; 
	  n_ljets--;
	}

	if(stopt.pfTau_leadPtcandID()!=(-1)) {
	  if(deltaR(stopt.pfjets().at(i).eta() , stopt.pfjets().at(i).phi() , stopt.pfTau().eta(), stopt.pfTau().phi())>0.4) continue; 
	  n_taujets--;
	}
	
      }


      ///////                                                                                                                                                     
      /// SIGNAL REGION 1l + 3b                                                                                                                                   
      ///                                                                                                                                                         

      if(passSingleLeptonSelection(isData) && stopt.t1metphicorr()>50 && stopt.t1metphicorrmt()>175 && jetsP4.size()>=4 && passIsoTrkVeto_v4() && passTauVeto()) {

        plot1D("h_TABLE_1l",  0  ,       evtweight, h_1d, 5, 0, 5);

	vector<int> indexLooseB=getBJetIndex(0.240, -1., -1., jetsP4, csv, lrm , 40.,2.1, false);

	if(indexLooseB.size()>2) {
	  plot1D("h_TABLE_1l",  1  ,       evtweight, h_1d, 5, 0, 5);

	  vector<int> indexFirstTight=getBJetIndex(0.898, indexLooseB.at(0), indexLooseB.at(1),jetsP4, csv,lrm , 40.,2.1,false);
	  if(indexFirstTight.size()>0 && indexFirstTight.size()==(indexLooseB.size()-2)) {
	    plot1D("h_TABLE_1l",  2  ,       evtweight, h_1d, 5, 0, 5);
	    //	    plot1D("h_TABLE_1l",  3  ,       evtweight, h_1d, 5, 0, 5);
	    if(jetsP4.size()>=5) plot1D("h_TABLE_1l",  4  ,       evtweight, h_1d, 5, 0, 5);

	    classify3B(evtweight, "_met50_4j_3bLT" , h_1d_qg, h_2d_qg, isData,false);

	  }
	}
      }

    
      ///////                                                                                                                                                  
      /// SIGNAL REGION OS + 3b                                                                                                                                
      /// 
      if(stopt.ngoodlep()==2
         && passDileptonSelection(isData)
         && (abs(stopt.id1()) != abs(stopt.id2()) || fabs( stopt.dilmass() - 91.) > 15. )
         //      && stopt.t1metphicorrmt()>150                                                                                                                 
	 && stopt.t1metphicorr()>50
         && n_ljets>=4) {

        float mt2_standard = MT2( stopt.t1metphicorr(), stopt.t1metphicorrphi(), stopt.lep1() , stopt.lep2() );
       
	vector<int> indexLooseB=getBJetIndex(0.240, -1., -1., jetsP4, csv,lrm, 40.,2.1,false);

        if(indexLooseB.size()<3) plot1D("h_mt2_OS_12",  mt2_standard  ,       evtweight, h_1d, 50, 0, 200);
        if(indexLooseB.size()>=3) plot1D("h_mt2_OS_3b",  mt2_standard  ,       evtweight, h_1d, 50, 0, 200);
 
	plot1D("h_TABLE_2l",  0  ,       evtweight, h_1d, 5, 0, 5);

        if(indexLooseB.size()>2) {
	  plot1D("h_TABLE_2l",  1  ,       evtweight, h_1d, 5, 0, 5);
	  vector<int> indexFirstTight=getBJetIndex(0.898, indexLooseB.at(0), indexLooseB.at(1),jetsP4, csv,lrm, 40.,2.1,false);
	  if(indexFirstTight.size()>0 && indexFirstTight.size()==(indexLooseB.size()-2)) {
	    plot1D("h_TABLE_2l",  2  ,       evtweight, h_1d, 5, 0, 5);
	    if(mt2_standard>80) plot1D("h_TABLE_2l",  3  ,       evtweight, h_1d, 5, 0, 5);
	    if(n_ljets>=5) plot1D("h_TABLE_2l",  4  ,       evtweight, h_1d, 5, 0, 5);
	    
	    classify3B(evtweight, "_OS_met50_4j_3bLT" , h_1d_qg, h_2d_qg, isData,false);

	  }

	}

      }


      ///////                                                                                                                                                  
      /// CONTROL REGION OS + 3b                                                                                                                                
      /// 
      if ( dataset_CR4
           && stopt.t1metphicorr()>50
           && stopt.ngoodlep()==2
           && passDileptonSelection(isData)
           && (abs(stopt.id1()) != abs(stopt.id2()) || fabs( stopt.dilmass() - 91.) > 15. )
           && n_ljets>=4
           )
        //         && btag.size()>0 )                                                                                                                                          
        {
          float trigweight = isData ? 1. : getdltrigweight(stopt.id1(), stopt.id2());

	  vector<int> indexLooseB=getBJetIndex(0.240, -1., -1., jetsP4, csv, lrm , 40.,2.1, false);
	  vector<int> indexMediumB=getBJetIndex(0.679, -1., -1., jetsP4, csv, lrm , 40.,2.1, false);
	  
	  if(indexMediumB.size()>0)    plotCR4(evtweight*trigweight, "_1bM", h_1d_cr4, isData);

	  if(indexLooseB.size()>2) {
	    
	    vector<int> indexFirstTight=getBJetIndex(0.898, indexLooseB.at(0), indexLooseB.at(1),jetsP4, csv,lrm , 40.,2.1,false);
	    if(indexFirstTight.size()>0 && indexFirstTight.size()==(indexLooseB.size()-2)) {
	      plotCR4(evtweight*trigweight, "_3bLT", h_1d_cr4, isData);
	    }

	  }

	}


      ///////                                                                                                                                                  
      /// CONTROL REGION CR5 + 3b                                                                                                                                
      /// 
      if ( dataset_1l                                                                                                                                          
	   && stopt.t1metphicorr()>50                                                                                                                          
	   && passLepPlusIsoTrkSelection(isData)
	   && n_ljets>=4
	   //	   && passLepPlusTauSelection(isData)
	   //	   && n_taujets >=4
	   )
	{
	  
	  float trigweight = isData ? 1. : getsltrigweight(stopt.id1(), stopt.lep1().Pt(), stopt.lep1().Eta());

	  vector<int> indexLooseB=getBJetIndex(0.240, -1., -1., jetsP4, csv, lrm , 40.,2.1, false);
	  vector<int> indexMediumB=getBJetIndex(0.679, -1., -1., jetsP4, csv, lrm , 40.,2.1, false);
	  
	  if(indexMediumB.size()>0)   {
	    studyIsoTrack(evtweight*trigweight, "_1bM", h_1d_cr5, h_2d, isData);
	    studyIsoTaus(evtweight*trigweight, "_1bM", h_1d_cr5, h_2d, isData);
	    //	    plotCR5(evtweight*trigweight, "_1bM", h_1d_cr5, isData);

	  }

	  if(indexLooseB.size()>2) {
	    
	    vector<int> indexFirstTight=getBJetIndex(0.898, indexLooseB.at(0), indexLooseB.at(1),jetsP4, csv,lrm , 40.,2.1,false);
	    if(indexFirstTight.size()>0 && indexFirstTight.size()==(indexLooseB.size()-2)) {
	      studyIsoTrack(evtweight*trigweight, "_3bLT", h_1d_cr5, h_2d, isData);
	      studyIsoTaus(evtweight*trigweight, "_3bLT", h_1d_cr5, h_2d, isData);
	      //	      plotCR5(evtweight*trigweight, "_3bLT", h_1d_cr5, isData);
	    }
	  }
	}


      ///////                                                                                                                                                  
      /// CONTROL REGION CR1 + 3b                                                                                                                                
      /// 
      if ( dataset_1l                                                                                                                                          
	   && stopt.t1metphicorr()>50                                                                                                                          
	   && stopt.t1metphicorrmt()>50                                                                                                                          
	   && stopt.t1metphicorrmt()<100                                                                                                                          
	   && passSingleLeptonSelection(isData)
           && passIsoTrkVeto_v4() && passTauVeto()
	   && jetsP4.size()>=4
	   //	   && passLepPlusTauSelection(isData)
	   //	   && n_taujets >=4
	   )
	{
	  
	  float trigweight = isData ? 1. : getsltrigweight(stopt.id1(), stopt.lep1().Pt(), stopt.lep1().Eta());

	  vector<int> indexLooseB=getBJetIndex(0.240, -1., -1., jetsP4, csv, lrm , 40.,2.1, false);
	  vector<int> indexMediumB=getBJetIndex(0.679, -1., -1., jetsP4, csv, lrm , 40.,2.1, false);
	  
	  if(indexMediumB.size()>0)   {
	    plotCR1(evtweight*trigweight, "_1bM", h_1d_cr1, isData);
	  }

	  if(indexLooseB.size()>2) {
	    
	    vector<int> indexFirstTight=getBJetIndex(0.898, indexLooseB.at(0), indexLooseB.at(1),jetsP4, csv,lrm , 40.,2.1,false);
	    if(indexFirstTight.size()>0 && indexFirstTight.size()==(indexLooseB.size()-2)) {
	      plotCR1(evtweight*trigweight, "_3bLT", h_1d_cr1, isData);
	      //	      plotCR5(evtweight*trigweight, "_3bLT", h_1d_cr5, isData);
	    }
	  }
	}

      ///////                                                                                                                                                  
      /// CONTROL REGION CR3 + 3b                                                                                                                                
      /// 

      if ( dataset_CR4
           && stopt.t1metphicorr()>50
           && stopt.ngoodlep()==2
           && passDileptonSelection(isData)
           && (abs(stopt.id1()) != abs(stopt.id2()) || fabs( stopt.dilmass() - 91.) > 15. )
           && n_ljets>=4
           )
	{

          float trigweight = isData ? 1. : getdltrigweight(stopt.id1(), stopt.id2());

	  vector<int> indexLooseB=getBJetIndex(0.240, -1., -1., jetsP4, csv, lrm , 40.,2.1, false);

	  if(indexLooseB.size()>0 && indexLooseB.size()<3) plotCR3(evtweight*trigweight, "_12bL", h_1d_cr3, isData); 


	}


      ///////                                                                                                                                                  
      /// CONTROL REGION CR2 + 3b                                                                                                                                
      /// 

      if ( dataset_CR4
           && stopt.t1metphicorr()>50
           && stopt.ngoodlep()==2
           && passDileptonSelection(isData)
           && (abs(stopt.id1()) != abs(stopt.id2()) || fabs( stopt.dilmass() - 91.) > 15. )
           && n_ljets>=4
           )
	{

	  float mt2_standard = MT2( stopt.t1metphicorr(), stopt.t1metphicorrphi(), stopt.lep1() , stopt.lep2() );
	  if(mt2_standard>80) continue;

          float trigweight = isData ? 1. : getdltrigweight(stopt.id1(), stopt.id2());

	  vector<int> indexLooseB=getBJetIndex(0.240, -1., -1., jetsP4, csv, lrm , 40.,2.1, false);

	  if(indexLooseB.size()>2) {
	    
	    vector<int> indexFirstTight=getBJetIndex(0.898, indexLooseB.at(0), indexLooseB.at(1),jetsP4, csv,lrm , 40.,2.1,false);
	    if(indexFirstTight.size()>0 && indexFirstTight.size()==(indexLooseB.size()-2)) {
	      plotCR2(evtweight*trigweight, "_3bLT", h_1d_cr2, isData);
	      //	      plotCR5(evtweight*trigweight, "_3bLT", h_1d_cr5, isData);
	    }
	  }



	}




    } // end event loop

    // delete tree;
    
  } // end file loop
  
    //
    // finish
    //
  
  TFile outfile(m_outfilename_.c_str(),"RECREATE") ; 
  printf("[StopTreeLooper::loop] Saving histograms to %s\n", m_outfilename_.c_str());
  
  std::map<std::string, TH1F*>::iterator it1d;
  for(it1d=h_1d.begin(); it1d!=h_1d.end(); it1d++) {
    it1d->second->Write(); 
    delete it1d->second;
  }
  

  outfile.mkdir("QGTag");
  outfile.cd("QGTag");

  for(it1d=h_1d_qg.begin(); it1d!=h_1d_qg.end(); it1d++) {

    it1d->second->Write();
    delete it1d->second;
  }

  std::map<std::string, TH2F*>::iterator it2d;
  for(it2d=h_2d_qg.begin(); it2d!=h_2d_qg.end(); it2d++) {
    it2d->second->Write();
    delete it2d->second;

  }

  /////
  /////
  /////

  outfile.mkdir("CR2");
  outfile.cd("CR2");

  for(it1d=h_1d_cr2.begin(); it1d!=h_1d_cr2.end(); it1d++) {

    it1d->second->Write();
    delete it1d->second;
  }

  outfile.mkdir("CR1");
  outfile.cd("CR1");

  for(it1d=h_1d_cr1.begin(); it1d!=h_1d_cr1.end(); it1d++) {

    it1d->second->Write();
    delete it1d->second;
  }

  outfile.mkdir("CR3");
  outfile.cd("CR3");

  for(it1d=h_1d_cr3.begin(); it1d!=h_1d_cr3.end(); it1d++) {

    it1d->second->Write();
    delete it1d->second;
  }

  outfile.mkdir("CR4");
  outfile.cd("CR4");

  for(it1d=h_1d_cr4.begin(); it1d!=h_1d_cr4.end(); it1d++) {

    it1d->second->Write();
    delete it1d->second;
  }

  outfile.mkdir("CR5");
  outfile.cd("CR5");

  for(it1d=h_1d_cr5.begin(); it1d!=h_1d_cr5.end(); it1d++) {

    it1d->second->Write();
    delete it1d->second;
  }

  /////
  /////
  /////

  outfile.Write();
  outfile.Close();

  already_seen.clear();

  gROOT->cd();

}


void StopTreeLooper::classify3B(float evtweight, string tag_selection, std::map<std::string, TH1F*> &h_1d, std::map<std::string, TH2F*> &h_2d, bool isData, bool doLRM) {

  int nLooseb=0;
  int nb=0;
  int nW=0;
  int nG=0;
  int nL=0;
  int nT=0;
  int nUnk=0;

  for( unsigned int i = 0 ; i < jetsP4.size() ; ++i ){

    if( jetsP4.at(i).pt()<40 )  continue;
    if( fabs(jetsP4.at(i).eta())>2.1 )  continue;

    if( csv.at(i)<0.240 )  continue;
    if( doLRM && (lrm.at(i) < passLRM(jetsP4.at(i).pt()))) continue;  

    float yProd=stopt.pfjets().at(i).Rapidity()*stopt.lep1().Rapidity();

    if(abs(mc3.at(i))==1)             plot1D("h_yProd_b"+tag_selection,       yProd,       evtweight, h_1d, 100, -10, 10.);
    if(abs(mc3.at(i))==2 || abs(mc3.at(i))==5)   plot1D("h_yProd_W"+tag_selection,       yProd,       evtweight, h_1d, 100, -10, 10.);
    if(abs(mc3.at(i))==3 || abs(mc3.at(i))==4)   plot1D("h_yProd_G"+tag_selection,       yProd,       evtweight, h_1d, 100, -10, 10.);
    if(abs(mc3.at(i))==11 || abs(mc3.at(i))==13) plot1D("h_yProd_L"+tag_selection,       yProd,       evtweight, h_1d, 100, -10, 10.);
    if(abs(mc3.at(i))==15)            plot1D("h_yProd_Tau"+tag_selection,     yProd,       evtweight, h_1d, 100, -10, 10.);
    if(abs(mc3.at(i))==9)             plot1D("h_yProd_uncl"+tag_selection,    yProd,       evtweight, h_1d, 100, -10, 10.);

    float dr=deltaR(stopt.pfjets().at(i).eta() , stopt.pfjets().at(i).phi() , stopt.lep1().eta(), stopt.lep1().phi());

    if(abs(mc3.at(i))==1)             plot1D("h_drbl_b"+tag_selection,       dr,       evtweight, h_1d, 100, 0, 6.5);
    if(abs(mc3.at(i))==2 || abs(mc3.at(i))==5)   plot1D("h_drbl_W"+tag_selection,       dr,       evtweight, h_1d, 100, 0, 6.5);
    if(abs(mc3.at(i))==3 || abs(mc3.at(i))==4)   plot1D("h_drbl_G"+tag_selection,       dr,       evtweight, h_1d, 100, 0, 6.5);
    if(abs(mc3.at(i))==11 || abs(mc3.at(i))==13) plot1D("h_drbl_L"+tag_selection,       dr,       evtweight, h_1d, 100, 0, 6.5);
    if(abs(mc3.at(i))==15)            plot1D("h_drbl_Tau"+tag_selection,     dr,       evtweight, h_1d, 100, 0, 6.5);
    if(abs(mc3.at(i))==9)             plot1D("h_brbl_uncl"+tag_selection,    dr,       evtweight, h_1d, 100, 0, 6.5);

    ///    if(yProd<(-0.5)) continue;

    nLooseb++;

    plot1D("h_mc3"+tag_selection,   abs(mc3.at(i)),       evtweight, h_1d, 25, 0, 25);
    plot1D("h_mc1"+tag_selection,   abs(mc1.at(i)),       evtweight, h_1d, 25, 0, 25);

    if(abs(mc3.at(i))==1) nb++;
    if(abs(mc3.at(i))==2 || abs(mc3.at(i))==5) nW++;
    if(abs(mc3.at(i))==3 || abs(mc3.at(i))==4) nG++;
    if(abs(mc3.at(i))==11 || abs(mc3.at(i)) ==13) nL++;
    if(abs(mc3.at(i))==15) nT++;
    if(abs(mc3.at(i))==9) nUnk++;



  }

  if(nLooseb<3) return;


  int myClass=0;

  if(nb==3 && nUnk==0) myClass=1;
  if(nb==2 && nW==1 && nUnk==0) myClass=2;
  if(nb==2 && nW==2 && nUnk==0) myClass=3;
  if(nb==1 && nW==2 && nUnk==0) myClass=4;
  if(nb==2 && nT==1 && nUnk==0) myClass=5; // bbtau                                                                                                            
  if(nb==2 && nL==1 && nUnk==0) myClass=6;

  //------                                                                                                                                                     

  if(nb==1 && nW==1 && nG==1 && nUnk==0) myClass=7; //bWG                                                                                                      
  if(nb==1 && nG==2 && nUnk==0) myClass=8; // bGG                                                                                                              
  if(nb==2 && (nG==1 || nG==2) && nUnk==0) myClass=9; //bbG + bbGG                                                                                             
  if((nb==1 || nb==2) && nT==1 && nG==1 && nUnk==0) myClass=10; //btauG + bbtauG                                                                               
  if((nb==1 || nb==2) && nL==1 && nG==1 && nUnk==0) myClass=11; // bLG                                                                                         

  if(nUnk!=0) myClass=15;

  plot1D("h_myclass"+tag_selection,    myClass,       evtweight, h_1d, 16, 0, 16);

}

void StopTreeLooper::plotCR2(float evtweight, string tag_selection , std::map<std::string, TH1F*> &h_1d, bool isData) {


  float mt2_standard = MT2( stopt.t1metphicorr(), stopt.t1metphicorrphi(), stopt.lep1() , stopt.lep2() );

  plot1D("cr2_mt2"+tag_selection,   mt2_standard ,       evtweight, h_1d, 100,  0, 200);
  plot1D("cr2_met"+tag_selection,   stopt.t1metphicorr() ,       evtweight, h_1d, 100, 0, 300);
  plot1D("cr2_njets"+tag_selection,  jetsP4.size() ,       evtweight, h_1d, 10,  1., 11.);

}


void StopTreeLooper::plotCR3(float evtweight, string tag_selection , std::map<std::string, TH1F*> &h_1d, bool isData) {


  float mt2_standard = MT2( stopt.t1metphicorr(), stopt.t1metphicorrphi(), stopt.lep1() , stopt.lep2() );

  plot1D("cr3_mt2"+tag_selection,   mt2_standard ,       evtweight, h_1d, 100,  0, 200);
  plot1D("cr3_met"+tag_selection,   stopt.t1metphicorr() ,       evtweight, h_1d, 100, 0, 300);
  plot1D("cr3_njets"+tag_selection,  jetsP4.size() ,       evtweight, h_1d, 10,  1., 11.);

}


void StopTreeLooper::plotCR4(float evtweight, string tag_selection , std::map<std::string, TH1F*> &h_1d, bool isData) {

  int btagLoose=0;

  for( unsigned int i = 0 ; i < jetsP4.size() ; ++i ){

    if( jetsP4.at(i).pt()<30 )  continue;
    if( fabs(jetsP4.at(i).eta())>2.4 )  continue;

    if (csv.at(i) > 0.240) btagLoose++; 
    if (csv.at(i) > 0.240) plot1D("cr4_eta_jets"+tag_selection,  jetsP4.at(i).eta() ,       evtweight, h_1d, 100,  -5., 5.);
    if (csv.at(i) > 0.240) plot1D("cr4_pt_jets"+tag_selection,  jetsP4.at(i).pt() ,       evtweight, h_1d, 100,  0., 300.);
    if (csv.at(i) > 0.240) plot1D("cr4_discr_jets"+tag_selection,  csv.at(i),       evtweight, h_1d, 100,  0., 1.);

    if (csv.at(i) < 0.240) plot1D("cr4_eta_nonbjets"+tag_selection,  jetsP4.at(i).eta() ,       evtweight, h_1d, 100,  -5., 5.);

    if (csv.at(i) > 0.240) {

      float yProd= stopt.lep1().Rapidity()*jetsP4.at(i).Rapidity() ;
      plot1D("cr4_yProd"+tag_selection,  yProd ,       evtweight, h_1d, 100, -5., 5.);

    }

  }

  plot1D("cr4_mt"+tag_selection,   stopt.mt() ,       evtweight, h_1d, 100,  0, 300);
  plot1D("cr4_met"+tag_selection,   stopt.t1metphicorr() ,       evtweight, h_1d, 100, 0, 300);

  plot1D("cr4_njets"+tag_selection,  jetsP4.size() ,       evtweight, h_1d, 10,  1., 11.);
  plot1D("cr4_nbjets"+tag_selection,   btagLoose ,       evtweight, h_1d, 10,  1., 11.);

  float pt12=(stopt.lep1()+stopt.lep2()).pt();
  float m12=(stopt.lep1()+stopt.lep2()).M();
  float Rap12=(stopt.lep1()+stopt.lep2()).Rapidity();

  float dr=deltaR(stopt.lep1().eta() , stopt.lep1().phi() , stopt.lep2().eta(), stopt.lep2().phi());
  float dphi=deltaPhi(stopt.lep1().phi(), stopt.lep2().phi());

  plot1D("cr4_dilep_m"+tag_selection,  m12 ,       evtweight, h_1d, 100, 0, 300.);
  plot1D("cr4_dilep_pt"+tag_selection,  pt12 ,       evtweight, h_1d, 100, 0, 500.);
  plot1D("cr4_dilep_rap"+tag_selection,  Rap12 ,       evtweight, h_1d, 100, -3., 3.);
  plot1D("cr4_dilep_dr"+tag_selection,  dr ,       evtweight, h_1d, 100, 0, 6.);
  plot1D("cr4_dilep_dphi"+tag_selection,  dphi ,       evtweight, h_1d, 100, 0., 3.14);


}

void StopTreeLooper::plotCR1(float evtweight, string tag_selection , std::map<std::string, TH1F*> &h_1d, bool isData) {

  int btagLoose=0;

  for( unsigned int i = 0 ; i < jetsP4.size() ; ++i ){

    if( jetsP4.at(i).pt()<30 )  continue;
    if( fabs(jetsP4.at(i).eta())>2.4 )  continue;

    if (csv.at(i) > 0.240) btagLoose++; 
    if (csv.at(i) > 0.240) plot1D("cr1_eta_jets"+tag_selection,  jetsP4.at(i).eta() ,       evtweight, h_1d, 100,  -5., 5.);
    if (csv.at(i) > 0.240) plot1D("cr1_pt_jets"+tag_selection,  jetsP4.at(i).pt() ,       evtweight, h_1d, 100,  0., 300.);
    if (csv.at(i) > 0.240) plot1D("cr1_discr_jets"+tag_selection,  csv.at(i),       evtweight, h_1d, 100,  0., 1.);

    if (csv.at(i) < 0.240) plot1D("cr1_eta_nonbjets"+tag_selection,  jetsP4.at(i).eta() ,       evtweight, h_1d, 100,  -5., 5.);

  }

  plot1D("cr1_mt"+tag_selection,   stopt.mt() ,       evtweight, h_1d, 100,  0, 300);
  plot1D("cr1_met"+tag_selection,   stopt.t1metphicorr() ,       evtweight, h_1d, 100, 0, 300);

  plot1D("cr1_njets"+tag_selection,  jetsP4.size() ,       evtweight, h_1d, 10,  1., 11.);
  plot1D("cr1_nbjets"+tag_selection,   btagLoose ,       evtweight, h_1d, 10,  1., 11.);

}


void StopTreeLooper::studyIsoTaus(float evtweight, string tag_selection, std::map<std::string, TH1F*> &h_1d, std::map<std::string, TH2F*> &h_2d, bool isData) {  


  plot1D("cr5_njets_taus"+tag_selection, jetsP4.size() ,           evtweight, h_1d, 10,  1., 11.);                                                            

  //    n_taujets--;

}


void StopTreeLooper::studyIsoTrack(float evtweight, string tag_selection, std::map<std::string, TH1F*> &h_1d, std::map<std::string, TH2F*> &h_2d, bool isData) {  
                                                                                                                                                                  
  //  if(abs(stopt.pfcandidOS10())!=211) cout << "particle ID " << stopt.pfcandidOS10() << endl;                                                                  
                         
  float charge=(stopt.pfcandidOS10looseZ()*stopt.id1());                                                                                                                
  // charge < 0 is a SS , charge > 0 is a OS for e/mu; need to flip for pions                                                                                     
  if((abs(stopt.pfcandidOS10looseZ())!=11) && (abs(stopt.pfcandidOS10looseZ())!=13)) charge*=(-1);                                                                            
                                                                                                                                                                  
  string tag_charge="";                                                                                                                                              
  if(charge>0) tag_charge="_OS";                                                                                                                                  
  if(charge<0) tag_charge="_SS";                                                                                                                                  
  if(stopt.pfcandptOS10looseZ() >=9998.) tag_charge="_undefined";
  ////pfcands_particleId().at(ipf)0) tag_charge="_SS";                                                                                                            

  string tag_flav="";                                                                                                                                              
  if(abs(stopt.pfcandidOS10looseZ())==11 || abs(stopt.pfcandidOS10looseZ())==13 ) {
    tag_flav="_EMU";                                                                                                                                  
  } else {
    tag_flav="_noEMU";                                                                                                                                  
  }

                                                                                                                                                                  
  //  if((abs(stopt.pfcandidOS10())!=11) && (abs(stopt.pfcandidOS10())!=13)  && (abs(stopt.pfcandidOS10())!=211)) cout << 'undefined ' << stopt.pfcandid10() << endl;
                                                                                                                                                                  
  plot1D("cr5_njets"+tag_selection+tag_flav+tag_charge, jetsP4.size() ,           evtweight, h_1d, 10,  1., 11.);                                                            
                                                                                                                                                                  
  plot1D("cr5_iso"+tag_selection+tag_flav+tag_charge, stopt.pfcandisoOS10looseZ() ,           evtweight, h_1d, 100, 0., 1.);                                                     
  plot1D("cr5_pt"+tag_selection+tag_flav+tag_charge, stopt.pfcandptOS10looseZ() ,           evtweight, h_1d, 10,  1., 100.);                                                     
                                                                                                                                                                  
}     
