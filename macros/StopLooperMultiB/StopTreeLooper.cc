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
//#include "../../Tools/BTagReshaping/BTagReshaping.h"
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
  std::map<std::string, TH1F*> h_1d_cr6;
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

      dataset_1l=false;

      if((isData) && name.Contains("muo") && (abs(stopt.id1()) == 13 ))  dataset_1l=true;
      if((isData) && name.Contains("ele") && (abs(stopt.id1()) == 11 ))  dataset_1l=true;

      if(!isData) dataset_1l=true;

      dataset_CR4=false;

      /*
	if((isData) && name.Contains("dimu") && (abs(stopt.id1()) == 13 ) && (abs(stopt.id2())==13) && (fabs( stopt.dilmass() - 91.) > 15.)) dataset_CR4=true;
	if((isData) && name.Contains("diel") && (abs(stopt.id1()) == 11 ) && (abs(stopt.id2())==11) && (fabs( stopt.dilmass() - 91.) > 15.)) dataset_CR4=true;
	if((isData) && name.Contains("mueg") && abs(stopt.id1()) != abs(stopt.id2())) dataset_CR4=true;
      */

      if((isData) && name.Contains("dimu") && (abs(stopt.id1()) == 13 ) && (abs(stopt.id2())==13)) dataset_CR4=true;
      if((isData) && name.Contains("diel") && (abs(stopt.id1()) == 11 ) && (abs(stopt.id2())==11)) dataset_CR4=true;
      if((isData) && name.Contains("mueg") && abs(stopt.id1()) != abs(stopt.id2())) dataset_CR4=true;

      if(!isData) dataset_CR4=true;

      // remove the runs that are bad for the pixel alignement
      //      if(stopt.run()>=207883 && stopt.run()<=208307) continue;

      if(stopt.indexfirstGoodVertex_()) continue;

      //---------------------------------------------------------------------------- 
      // do WBB
      //---------------------------------------------------------------------------- 

      bool goodWjbb=true;
      //      if (name.Contains("wjets") && stopt.nbs() >=2 ) goodWjbb = false;
      //      if (name.Contains("wbb") && stopt.nbs() <2 ) goodWjbb = false;

      //---------------------------------------------------------------------------- 
      // determine event weight
      // make 2 example histograms of nvtx and corresponding weight
      //---------------------------------------------------------------------------- 

      float puweight = 1;

      if (!name.Contains("T2") && !name.Contains("T6")) {
        puweight = vtxweight_n( stopt.ntruepu(), h_pu_wgt, isData );
      }

      float evtweight = isData ? 1. : ( stopt.weight() * 19.5 * puweight) * goodWjbb;

      float nEts=100000;
      float BR=0.60*0.60;
      float lumi= 19.5;
      /* // for inclusive sample
      if(name.Contains("T6tthh_450")) evtweight = isData ? 1. : ( ( BR * 0.169668 * 1000.0 * lumi ) / nEts);
      if(name.Contains("T6tthh_350")) evtweight = isData ? 1. : ( ( BR * 0.807323 * 1000.0 * lumi ) / nEts);
      */
      if(name.Contains("T6tthh450")) evtweight = isData ? 1. : ( ( 1. * 0.169668 * 1000.0 * lumi ) / 499969);
      if(name.Contains("T6tthh350")) evtweight = isData ? 1. : ( ( 1. * 0.807323 * 1000.0 * lumi ) / 412464);

      if(name.Contains("T6ttzz_450")) evtweight = isData ? 1. : ( ( 0.169668 * 1000.0 * lumi ) / (5*nEts));
      if(name.Contains("T6ttzz_350")) evtweight = isData ? 1. : ( ( 0.807323 * 1000.0 * lumi ) / (5*nEts));

      // to reweight from file - also need to comment stuff before
      //      float vtxweight = vtxweight_n( nvtx, h_vtx_wgt, isData );

      plot1D("h_vtx",       stopt.nvtx(),       evtweight, h_1d, 40, 0, 40);
      plot1D("h_vtxweight", stopt.nvtxweight(), evtweight, h_1d, 41, -4., 4.);

      //----------------------------------------------------------------------------
      // apply preselection:
      // rho 0-40 GeV, MET filters, >=1 good lepton, veto 2 leptons dR < 0.1
      //----------------------------------------------------------------------------

      if ( !passEvtSelection(name,false) ) continue;

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
      indexJetsP4.clear();
      csv.clear();
      sigma.clear();
      mc3.clear();
      mc1.clear();
      lrm.clear();
      chm.clear();
      neu.clear();

      csvNonbMax=(-1);
      csvNonbTauMax=(-1);

      minMassbbCut=100;
      maxMassbbCut=150;

      bool dostudyMass = false;

      n_ljets=0;
      n_taujets=0;
      int n_trackjets=0;

      /*
	int indexF=-1.;
	int indexB=-1.;
	float etaF=0;
	float etaB=0;
	detaFB=-9999.;
      */

      bool rejectTOBTEC = false;

      for( unsigned int i = 0 ; i < stopt.pfjets().size() ; ++i ){

        if( stopt.pfjets().at(i).pt()<30 )  continue;
        if( fabs(stopt.pfjets().at(i).eta())>2.4 )  continue;

	//	if( stopt.pfjets_beta2_0p5().at(i) < 0.2) continue;
	bool pass5xPUid = passMVAJetId(stopt.pfjets().at(i).pt(), stopt.pfjets().at(i).eta(),stopt.pfjets_mva5xPUid().at(i),0);
        if(!pass5xPUid) continue;              
       
	float csv_nominal=isData ? stopt.pfjets_csv().at(i)
	  : nominalShape->reshape( stopt.pfjets().at(i).eta(),
				   stopt.pfjets().at(i).pt(),
				   stopt.pfjets_csv().at(i),
				   stopt.pfjets_mcflavorAlgo().at(i) ); 
       
	//        float csv_nominal = stopt.pfjets_csv().at(i) ;

        // fill again because of the beta and csv                                                                                                                      
	indexJetsP4.push_back(i);
        jetsP4.push_back( stopt.pfjets().at(i) );
        csv.push_back( csv_nominal );
        mc3.push_back( stopt.pfjets_mc3().at(i) );
        mc1.push_back( stopt.pfjets_mcflavorAlgo().at(i) );
        sigma.push_back( stopt.pfjets_sigma().at(i) );
        lrm.push_back( stopt.pfjets_lrm().at(i));
        chm.push_back( stopt.pfjets_chm().at(i));
        neu.push_back( stopt.pfjets_neu().at(i));

	/*	
	  if(csv_nominal>0.240 && stopt.pfjets().at(i).eta()>etaF) {
          etaF=stopt.pfjets().at(i).eta();
          indexF=i;
	  }
	  
	  if(csv_nominal>0.240 && stopt.pfjets().at(i).eta()<etaB) {
          etaB=stopt.pfjets().at(i).eta();
          indexB=i;
	  }
	*/
	
	if(abs(stopt.pfjets().at(i).eta())>1 && (stopt.pfjets_chm().at(i) - stopt.pfjets_neu().at(i)) > 40 ) rejectTOBTEC = true;

        n_ljets++;
        n_taujets++;
        n_trackjets++;

	if(stopt.id2()!=(-1)) {
	  if(deltaR(stopt.pfjets().at(i).eta() , stopt.pfjets().at(i).phi() , stopt.lep2().eta(), stopt.lep2().phi())>0.4) continue; 
	  n_ljets--;
	}

	if(stopt.pfTau_leadPtcandID()!=(-1)) {
	  if(deltaR(stopt.pfjets().at(i).eta() , stopt.pfjets().at(i).phi() , stopt.pfTau().eta(), stopt.pfTau().phi())>0.4) continue; 
	  n_taujets--;
	}

	if(stopt.pfcandptOS10looseZ() <9998.) {
          if(deltaR(stopt.pfjets().at(i).eta() , stopt.pfjets().at(i).phi() , stopt.pfcandOS10looseZ().eta(), stopt.pfcandOS10looseZ().phi()) > 0.4) continue;
          n_trackjets--;
        }
	
      }

      /*
      // vetoing events with second lepton, so this is redundant

      // TAU and second lepton overlap                                                                                                                                                                 
      if(stopt.id2()!=(-1)) {
        if(deltaR(stopt.lep2().eta() , stopt.lep2().phi() , stopt.pfTau().eta(), stopt.pfTau().phi())>0.4) n_taujets--;
      }
      */

      //  cout << " indexP4 " << indexJetsP4.size() << " P4 " << jetsP4.size() << endl;

      if(rejectTOBTEC) continue;

      //      if(indexF!=(-1) && indexB!=(-1)) detaFB=etaF-etaB;

      //----------------------------------------------------------------------------
      // GET bjet collection
      //----------------------------------------------------------------------------

      bool doLRM=false;
      tightDiscr=0.85;
      looseDiscr=0.679;
      vector<int> indexMediumB=getBJetIndex(0.679, -1., -1., jetsP4, csv, lrm , 30.,2.4, doLRM, false); // not used
      //      vector<int> indexRank=getRankIndex(looseDiscr, jetsP4, csv, 40.,2.1); // not used

      /////

      //      indexLooseB30=getBJetIndex(looseDiscr, -1., -1., jetsP4, csv, lrm , 30.,2.1, doLRM, false); // this is for the 1l+4b and also the 2l +4b      
      indexLooseB30=getBJetIndex(looseDiscr, -1., -1., jetsP4, csv, lrm , 40.,2.1, doLRM, false); // this is for the 1l+4b and also the 2l +4b      
      indexLooseB40=getBJetIndex(looseDiscr, -1., -1., jetsP4, csv, lrm , 40.,2.1, doLRM, false); // this is for the 2l+3b
      indexTightB40=getBJetIndex(tightDiscr, -1., -1., jetsP4, csv, lrm , 40., 2.1, doLRM, false); // this is for the 1l + 3b

      //      indexLooseB30Tau=getBJetIndex(looseDiscr, -1., -1., jetsP4, csv, lrm , 30.,2.1, doLRM, true);
      indexLooseB30Tau=getBJetIndex(looseDiscr, -1., -1., jetsP4, csv, lrm , 40.,2.1, doLRM, true);
      indexLooseB40Tau=getBJetIndex(looseDiscr, -1., -1., jetsP4, csv, lrm , 40.,2.1, doLRM, true);
      indexTightB40Tau=getBJetIndex(tightDiscr, -1., -1., jetsP4, csv, lrm , 40., 2.1, doLRM, true); // this is for the 1l + 3b

      csvNonbMax = getCSVNonb(jetsP4, csv, looseDiscr, 40., 2.1, false);
      csvNonbTauMax = getCSVNonb(jetsP4, csv, looseDiscr, 40., 2.1, true);

      //////////                                                                                                                                                     
      /// invariant mass selection                                                                                                                                  
      /////                                                                                                                                                         

      double deltaRHbb=-1.;
      massHbb=-1.;
      massHbb30=-1.;
      double ptHbb=-1.;
      double dilepVec=-1.;
      double deltaRap=-1.;
      double deltaPhibb=-1.;

      if(indexLooseB40.size()>1) {

	pair<int, int> index=getIndexPair(indexLooseB40,jetsP4,true,-1,-1);

        if(index.first!=(-1)) {
	  deltaRHbb= deltaR(jetsP4.at(index.first).eta(),jetsP4.at(index.first).phi() , jetsP4.at(index.second).eta(),jetsP4.at(index.second).phi());
	  massHbb= (jetsP4.at(index.first)+jetsP4.at(index.second)).M();
	  ptHbb= (jetsP4.at(index.first)+jetsP4.at(index.second)).pt();
	  deltaRap= abs(jetsP4.at(index.first).Rapidity()-jetsP4.at(index.second).Rapidity());
	  dilepVec=(jetsP4.at(index.first)+jetsP4.at(index.second)+stopt.lep1()+stopt.lep2()).pt();
	  deltaPhibb=deltaPhi(jetsP4.at(index.first).phi(),jetsP4.at(index.second).phi());
        }
      }


      if(indexTightB40.size()>1) {
	pair<int, int> index=getIndexPair(indexTightB40,jetsP4,true,-1,-1);

        if(index.first!=(-1)) {
	  massHbbTight= (jetsP4.at(index.first)+jetsP4.at(index.second)).M();
        }
      }

      /////

      massHbbTau=-1.;
      massHbbTauTight=-1.;

      if(indexLooseB40Tau.size()>1) {
        pair<int, int> indexTau=getIndexPair(indexLooseB40Tau,jetsP4,true,-1,-1);
        if(indexTau.first!=(-1)) {
	  massHbbTau = (jetsP4.at(indexTau.first)+jetsP4.at(indexTau.second)).M();
        }
      }

      if(indexTightB40Tau.size()>1) {
        pair<int, int> indexTau=getIndexPair(indexTightB40Tau,jetsP4,true,-1,-1);
        if(indexTau.first!=(-1)) {
	  massHbbTauTight = (jetsP4.at(indexTau.first)+jetsP4.at(indexTau.second)).M();
        }
      }

      //----------------------------------------------------------------------------
      // Fill tables Yield 
      //----------------------------------------------------------------------------

      fillSignalYieldTable(evtweight, ""  , h_1d, h_1d_qg, h_2d_qg, isData);

      //----------------------------------------------------------------------------
      // CONTROL REGION study  --- Single Lepton
      //----------------------------------------------------------------------------

      if ( dataset_1l
           && stopt.t1metphicorr()>50 
           && passSingleLeptonSelection(isData)
           && passIsoTrkVeto_v4() && passTauVeto()
           && jetsP4.size()>=4
	   && stopt.ngoodlep()==1
           )
        {

          float trigweight = isData ? 1. : getsltrigweight(stopt.id1(), stopt.lep1().Pt(), stopt.lep1().Eta());

	  if(indexLooseB30.size()>0 && indexLooseB30.size()<3) plot1D("cr1_mt_12bM",   stopt.t1metphicorrmt() ,       evtweight*trigweight , h_1d_cr1, 100,  0, 300);

          if(!isData && indexTightB40.size()>2 && jetsP4.size()>=5)   plot1D("h_cr1_mt_3bT",   stopt.t1metphicorrmt() ,       evtweight*trigweight , h_1d_cr1, 100,  0, 300);
          if(!isData && indexLooseB30.size()>3)   plot1D("h_cr1_mt_4bM",   stopt.t1metphicorrmt() ,       evtweight*trigweight , h_1d_cr1, 100,  0, 300);

          if(stopt.t1metphicorrmt()>50 && stopt.t1metphicorrmt()<100) {
            if(indexLooseB30.size()>0) plotCR1(evtweight*trigweight, "_1bM", h_1d_cr1, isData);
            if(indexLooseB30.size()==2) plotCR1(evtweight*trigweight, "_eq2bM", h_1d_cr1, isData);
            if(indexLooseB40.size()==3) plotCR1(evtweight*trigweight, "_eq3bM", h_1d_cr1, isData);

            if(indexLooseB30.size()>3) plotCR1(evtweight*trigweight, "_4bM", h_1d_cr1, isData);
            if(indexLooseB40.size()==3 && csvNonbMax>0.240) plotCR1(evtweight*trigweight, "_3bM", h_1d_cr1, isData); // this is for the 2l search 
	    //            if(indexTightB40.size()>2 && jetsP4.size()>=5) plotCR1(evtweight*trigweight, "_3bT", h_1d_cr1, isData); // this is for the 1l search 
          }

	  //	  if(indexLooseB40.size()==3 && dostudyMass) studyMassbb(evtweight*trigweight, "_1l_eq3bLoose40" , h_1d, h_2d, isData) ;
	  //	  if(indexLooseB40.size()>3 && dostudyMass) studyMassbb(evtweight*trigweight, "_1l_4bLoose40" , h_1d, h_2d, isData) ;

        }

      //----------------------------------------------------------------------------
      // CONTROL REGION study --- dilepton 
      //----------------------------------------------------------------------------

      // do outside of the loop to gain statistics
      if(indexLooseB40.size()==3 && dostudyMass) studyMassbb(1, "_2l_eq3bLoose40" , h_1d, h_2d, isData) ;
      if(indexLooseB40.size()>3 && dostudyMass) studyMassbb(1, "_2l_4bLoose40" , h_1d, h_2d, isData) ;

      ///////                                                                                                                                                  
      /// CONTROL REGION OS + 3b                                                                                                                                
      /// 
      if ( dataset_CR4
	   && stopt.t1metphicorr()>50
           && passDileptonSelection(isData)
	   //           && (abs(stopt.id1()) != abs(stopt.id2()) || fabs( stopt.dilmass() - 91.) > 15. )
	   && n_ljets>=4
           )
        {
          float trigweight = isData ? 1. : getdltrigweight(stopt.id1(), stopt.id2());
	  
	  
	  if(indexLooseB30.size()>0 && indexLooseB30.size()<3) plot1D("cr4_mt_12bM",   stopt.t1metphicorrmt() ,       evtweight*trigweight , h_1d_cr4, 100,  0, 300);
	  
	  bool doTAU=false;
	  
	  if(indexLooseB30.size()==2)  plotCR4(evtweight*trigweight, "_eq2bM", h_1d_cr4, isData,doLRM,doTAU);
	  
	  if(!isData && indexLooseB40.size()==3 && n_ljets>=5) plot1D("h_cr4_massHbb_3bM",  massHbb ,       evtweight*trigweight, h_1d_cr4, 51, 0, 500);
	  if(!isData && indexLooseB30.size()>3 ) plot1D("h_cr4_massHbb_4bM",  massHbbTau ,       evtweight*trigweight, h_1d_cr4, 51, 0, 500);
	  
	  /*
	  if(!isData && indexLooseB40.size()==3 && n_ljets==4) plot2D("h_njets_vs_nb", 0, 0,       evtweight, h_2d, 2, -0.5, 1.5,2, -0.5, 1.5);
	  if(!isData && indexLooseB40.size()>3 && n_ljets==4) plot2D("h_njets_vs_nb", 1, 0,       evtweight, h_2d, 2, -0.5, 1.5,2, -0.5, 1.5);
	  
	  if(!isData && indexLooseB40.size()==3 && n_ljets>4) plot2D("h_njets_vs_nb", 0, 1,       evtweight, h_2d, 2, -0.5, 1.5,2, -0.5, 1.5);
	  if(!isData && indexLooseB40.size()>3 && n_ljets>4) plot2D("h_njets_vs_nb", 1, 1,       evtweight, h_2d, 2, -0.5, 1.5,2, -0.5, 1.5);
	  */

	  if((massHbb < minMassbbCut || massHbb > maxMassbbCut) && indexLooseB30.size()>3 ) plotCR4(evtweight*trigweight, "_4bM_lowMassbb", h_1d_cr4, isData,doLRM,doTAU); // this is for the 2l search
	  if((massHbb < minMassbbCut) && indexLooseB40.size()==3 && n_ljets>=5) plotCR4(evtweight*trigweight, "_3bM_lowMassbb", h_1d_cr4, isData,doLRM,doTAU); // this is for the 2l search

	  if(massHbb < minMassbbCut || massHbb > maxMassbbCut) {

	    /*
	    if(indexLooseB40.size()==3 && n_ljets==4) plot2D("h_njets_vs_nb_lowMassbb", 0, 0,       evtweight, h_2d, 2, -0.5, 1.5,2, -0.5, 1.5);
	    if(indexLooseB40.size()>3 && n_ljets==4) plot2D("h_njets_vs_nb_lowMassbb", 1, 0,       evtweight, h_2d, 2, -0.5, 1.5,2, -0.5, 1.5);

	    if(indexLooseB40.size()==3 && n_ljets>4) plot2D("h_njets_vs_nb_lowMassb", 0, 1,       evtweight, h_2d, 2, -0.5, 1.5,2, -0.5, 1.5);
	    if(indexLooseB40.size()>3 && n_ljets>4) plot2D("h_njets_vs_nb_lowMassbb", 1, 1,       evtweight, h_2d, 2, -0.5, 1.5,2, -0.5, 1.5);
	    */

	    //	    if(indexTightB40.size()>2 && n_ljets>=5) plotCR4(evtweight*trigweight, "_3bT_lowMassbb", h_1d_cr4, isData,doLRM,doTAU); // this is for the 1l search

	    /*
	    if(n_ljets==4 && indexLooseB30.size()==4) {

	      cout << "========= found bad event in the CR  " << endl;
	      cout << "RUN " << stopt.run() << " lumi " << stopt.lumi() << " event " << stopt.event() << endl;
	      cout << "id1 " << stopt.id1() << " id2 " << stopt.id2() << endl;
	      cout << "pt1 " << stopt.lep1().pt() << " pt2 " << stopt.lep2().pt() << endl;
	      cout << "eta1 " << stopt.lep1().eta() << " eta2 " << stopt.lep2().eta() << endl;
	      cout << "phi1 " << stopt.lep1().phi() << " phi2 " << stopt.lep2().phi() << endl;
	      cout << "met " <<  stopt.t1metphicorr() << endl;

	      if(abs(stopt.id1())==11 && abs(stopt.lep1().eta())>=1.4) cout << " ++ electron1 in the endcap" << endl;
	      if(abs(stopt.id2())==11 && abs(stopt.lep2().eta())>=1.4) cout << " ++ electron2 in the endcap" << endl;
	      cout << "number good leptons " << stopt.ngoodlep() << endl;

	    }
	    */

	  }
	}

    
      ///////                                                                                                                                                  
      /// CONTROL REGION 1l + tau + 3b                                                                                                                                
      /// 
      
      if ( dataset_1l
           && stopt.t1metphicorr()>50
           && passLepPlusTauSelection_v2(isData)
           && n_taujets>=4
           )
        {

	  bool doTAU=true;

          float trigweight = isData ? 1. : getsltrigweight(stopt.id1(), stopt.lep1().Pt(), stopt.lep1().Eta());

	  if(indexLooseB30Tau.size()==2)  plotCR4(evtweight*trigweight, "_tau_eq2bM", h_1d_cr4, isData,doLRM,doTAU);

	  if(!isData && indexLooseB40Tau.size()>2 && n_taujets>=5) plot1D("h_cr4_tau_massHbb_3bM",  massHbbTau ,       evtweight*trigweight, h_1d_cr4, 51, 0, 500);
	  if(!isData && indexLooseB30Tau.size()>3 ) plot1D("h_cr4_tau_massHbb_4bM",  massHbbTau ,       evtweight*trigweight, h_1d_cr4, 51, 0, 500);
	 
          if( massHbbTau < minMassbbCut || massHbbTau > maxMassbbCut ) {

            if(indexLooseB30Tau.size()>3) plotCR4(evtweight*trigweight, "_tau_4bM_lowMassbb", h_1d_cr4, isData,doLRM,doTAU);
            if(indexLooseB40Tau.size()>2 && n_taujets>=5) plotCR4(evtweight*trigweight, "_tau_3bM_lowMassbb", h_1d_cr4, isData,doLRM,doTAU); // this is for the 2l search
            if(indexTightB40Tau.size()>2 && n_taujets>=5) plotCR4(evtweight*trigweight, "_tau_3bT_lowMassbb", h_1d_cr4, isData,doLRM,doTAU); // this is for the 1l search

          }
	
        }
    
      //----------------------------------------------------------------------------
      // END CONTROL REGION study
      //----------------------------------------------------------------------------

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

  outfile.mkdir("CR6");
  outfile.cd("CR6");

  for(it1d=h_1d_cr6.begin(); it1d!=h_1d_cr6.end(); it1d++) {

    it1d->second->Write();
    delete it1d->second;
  }

  outfile.mkdir("MASS");
  outfile.cd("MASS");

  for(it2d=h_2d.begin(); it2d!=h_2d.end(); it2d++) {

    it2d->second->Write();
    delete it2d->second;
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

  int nb_1=0;
  int nc_1=0;
  int nl_1=0;
  int n0=0;

  for( unsigned int i = 0 ; i < jetsP4.size() ; ++i ){

    if( jetsP4.at(i).pt()<40 )  continue;
    if( fabs(jetsP4.at(i).eta())>2.1 )  continue;

    if( csv.at(i)<0.698 )  continue;
    if( doLRM && (lrm.at(i) < passLRM(jetsP4.at(i).pt())) ) continue;  

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

    if(abs(mc1.at(i))==5) nb_1++;
    if(abs(mc1.at(i))==4) nc_1++;
    if((abs(mc1.at(i))>0 && abs(mc1.at(i))<=3) || abs(mc1.at(i))==21) nl_1++;
    if((abs(mc1.at(i)))==0 ) n0++;

  }

  if(nLooseb<3) return;

  int myClass1=0;

  if(nb_1>=3)  myClass1=1; // three b
  if(nb_1==2 && nc_1==1 )  myClass1=2; // bbc
  if(nb_1==2 && nc_1==0 && (nl_1>=1))  myClass1=3; // bbl
  if(nb_1<2)  myClass1=4; // 1b two mistag
  if(nb_1==2 && n0==1)  myClass1=5; // 2b plus one unkonw

  int myClass=0;

  if(nb==3 && nUnk==0) myClass=1;
  if(nb==2 && nW==1 && nUnk==0) myClass=2;
  if(nb==2 && nW==2 && nUnk==0) myClass=3;
  if(nb==1 && nW==2 && nUnk==0) myClass=4;
  if(nb==2 && (nT==1 || nT==2) && nUnk==0) myClass=5; // bbtau + bbtautau                                                                                                     
  if(nb==1 && nT==2 && nUnk==0) myClass=5; // btautau                                                                                                            
  if(nb==2 && nL==1 && nUnk==0) myClass=6;

  //------                                                                                                                                                     

  if(nb==1 && nW==1 && nG==1 && nUnk==0) myClass=7; //bWG                                                                                                      
  if(nb==1 && nG==2 && nUnk==0) myClass=8; // bGG                                                                                                              
  if(nb==2 && (nG==1 || nG==2) && nUnk==0) myClass=9; //bbG + bbGG                                                                                             
  if((nb==1 || nb==2) && nT==1 && nG==1 && nUnk==0) myClass=10; //btauG + bbtauG                                                                               
  if((nb==1 || nb==2) && nL==1 && nG==1 && nUnk==0) myClass=11; // bLG                                                                                         

  if(nUnk!=0) myClass=15;

  plot1D("h_myclass"+tag_selection,    myClass,       evtweight, h_1d, 16, 0, 16);

  plot1D("h_myclass1"+tag_selection,    myClass1,       evtweight, h_1d, 10, 0, 10);

}

void StopTreeLooper::fillSignalYieldTable(float evtweight, string tag_selection , std::map<std::string, TH1F*> &h_1d, std::map<std::string, TH1F*> &h_1d_qg,std::map<std::string, TH2F*> &h_2d_qg, bool isData) {

  //  cout <<"Mass bb "  << massHbb << endl;

  bool doLRM=false;
  
  ///////                                                                                                                                                     
  /// SIGNAL REGION 1l + 3b                                                                                                                                   
  ///                                                                                                                                                         
  
  float mtCut=150;
  
  if(passSingleLeptonSelection(isData) 
     && stopt.t1metphicorr()>50 
     && stopt.t1metphicorrmt()>120 
     && jetsP4.size()>=4 
     && stopt.ngoodlep()==1
     && passIsoTrkVeto_v4() && passTauVeto()) 
    {
    
    float trigweight = isData ? 1. : getsltrigweight(stopt.id1(), stopt.lep1().Pt(), stopt.lep1().Eta());
    
    if(stopt.t1metphicorrmt()>mtCut && jetsP4.size()>=5) plot1D("h_TABLE_1l",  0  ,       evtweight*trigweight, h_1d, 5, 0, 5);
    if(stopt.t1metphicorrmt()>mtCut && jetsP4.size()>=5 && indexLooseB40.size()>2) plot1D("h_TABLE_1l",  1  ,       evtweight*trigweight, h_1d, 5, 0, 5);
    if(stopt.t1metphicorrmt()>mtCut && jetsP4.size()>=5 && indexTightB40.size()>2) plot1D("h_TABLE_1l",  2  ,       evtweight*trigweight, h_1d, 5, 0, 5);
    if(stopt.t1metphicorrmt()>mtCut && indexLooseB40.size()==3 && jetsP4.size()>=5 && csvNonbMax>0.240) plot1D("h_TABLE_1l",  3  ,       evtweight*trigweight, h_1d, 5, 0, 5);
    if(indexLooseB30.size()>3) plot1D("h_TABLE_1l",  4  ,       evtweight*trigweight, h_1d, 5, 0, 5);
    
    if(indexLooseB40.size()>2 && jetsP4.size()>=5) plot1D("h_massHbb_1l"+tag_selection, massHbb ,       evtweight, h_1d, 102, -10, 500);
    if(indexLooseB40.size()>2) classify3B(evtweight, "_1l_met50_4j_3bM" , h_1d_qg, h_2d_qg, isData,doLRM);
    
    }
  
  ///////                                                                                                                                                     
  /// SIGNAL REGION OS + 3b                                                                                                                                   
  ///                                                                                                                                                         
  
  if(
     passDileptonSelection(isData)
     && stopt.t1metphicorr()>50
     && jetsP4.size()>=4)
    {
      
      float trigweight = isData ? 1. : getdltrigweight(stopt.id1(), stopt.id2());

      plot1D("h_TABLE_2l",  0  ,       evtweight*trigweight, h_1d, 5, 0, 5);
      if(indexLooseB40.size()>2) plot1D("h_TABLE_2l",  1  ,       evtweight*trigweight, h_1d, 5, 0, 5);
      if(indexLooseB40.size()>2 && jetsP4.size()>=5) plot1D("h_TABLE_2l",  2  ,       evtweight*trigweight, h_1d, 5, 0, 5);

      if(indexLooseB40.size()==3 && jetsP4.size()>=5 && massHbb>minMassbbCut && massHbb<maxMassbbCut) plot1D("h_TABLE_2l",  3  ,       evtweight*trigweight, h_1d, 5, 0, 5);
      if(indexLooseB30.size()>3 && massHbb>minMassbbCut && massHbb<maxMassbbCut) plot1D("h_TABLE_2l",  4  ,       evtweight*trigweight, h_1d, 5, 0, 5);

      if(indexLooseB40.size()>2 && jetsP4.size()>=5) plot1D("h_massHbb"+tag_selection, massHbb ,       evtweight, h_1d, 102, -10, 500);
      if(indexLooseB40.size()==3 && jetsP4.size()>=5) plot1D("h_massHbb_eq3b"+tag_selection, massHbb ,       evtweight, h_1d, 102, -10, 500);
      if(indexLooseB40.size()>3) plot1D("h_massHbb_ge4b"+tag_selection, massHbb ,       evtweight, h_1d, 102, -10, 500);
      if(indexLooseB40.size()>2) classify3B(evtweight, "_OS_met50_4j_3bM" , h_1d_qg, h_2d_qg, isData,doLRM);

      //      if(indexLooseB40.size()==3 && dostudyMass ) studyMassbb(evtweight*trigweight, "_2l_eq3bLoose40" , h_1d, h_2d, isData) ;
      //      if(indexLooseB40.size()==3 && dostudyMass && massHbb>minMassbbCut && massHbb<maxMassbbCut) studyMassbb(evtweight*trigweight, "_2l_eq3bLoose40_Mass150" , h_1d, h_2d, isData) ;

    }
    
  ///////                                                                                                                                                     
  /// SIGNAL REGION l+tau + 3b                                                                                                                                   
  ///                                                                                                                                                         
  if ( dataset_1l
       && stopt.t1metphicorr()>50
       && passLepPlusTauSelection_v2(isData)
       && n_taujets >=4
       )
    {
      
      float trigweight = isData ? 1. : getsltrigweight(stopt.id1(), stopt.lep1().Pt(), stopt.lep1().Eta());
      
      plot1D("h_TABLE_lt",  0  ,       evtweight*trigweight, h_1d, 5, 0, 5);
      if( indexLooseB40Tau.size()>2 ) plot1D("h_TABLE_lt",  1  ,       evtweight*trigweight, h_1d, 5, 0, 5);
      if( indexLooseB40Tau.size()>2 && n_taujets>=5 ) plot1D("h_TABLE_lt",  2  ,       evtweight*trigweight, h_1d, 5, 0, 5);
      if(indexLooseB40Tau.size()==3 && n_taujets>=5 && csvNonbTauMax>0.240 && massHbbTau>minMassbbCut && massHbbTau<maxMassbbCut) plot1D("h_TABLE_lt",  3  ,       evtweight*trigweight, h_1d, 5, 0, 5);
      if( indexLooseB30Tau.size()>3 && massHbbTau>minMassbbCut && massHbbTau<maxMassbbCut) plot1D("h_TABLE_lt",  4  ,       evtweight*trigweight, h_1d, 5, 0, 5);
  
      if(indexLooseB40Tau.size()>2 && n_taujets>=5) plot1D("h_massHbb"+tag_selection, massHbbTau ,       evtweight, h_1d, 102, -10, 500);

    }

}

void StopTreeLooper::studyMassbb(float evtweight, string tag_selection , std::map<std::string, TH1F*> &h_1d, std::map<std::string, TH2F*> &h_2d, bool isData) {

  //indexJetsP4
  //indexLooseB40

  pair<int, int> index;
  pair<int, int> indexOS_noClean;

  if(indexLooseB40.size()>1) {

    index=getIndexPair(indexLooseB40,jetsP4,true,-1,-1);
    indexOS_noClean=getIndexPair(indexLooseB40,jetsP4,false,-1,-1);

  }

  //  cout << " indexP4 " << indexJetsP4.size() << " P4 " << jetsP4.size() << endl;

  if(index.first!=(-1.)) {

    float massHbb1 = (jetsP4.at(index.first)+jetsP4.at(index.second)).M();

    bool goodMatch1 = ((mc3.at(index.first)==7 && mc3.at(index.second)==-7) || (mc3.at(index.first)==-7 && mc3.at(index.second)==7));
    bool goodMatchZ1 = ((mc3.at(index.first)==6 && mc3.at(index.second)==-6) || (mc3.at(index.first)==-6 && mc3.at(index.second)==6));
    //    bool goodMatch1 = (abs(mc3.at(index.first))==7 && abs(mc3.at(index.second)==7));

    float theta = getThetaStar(jetsP4.at(index.first), jetsP4.at(index.second));

    if(!goodMatch1) {

      plot1D("h_index_Pair1_NotgoodMatch"+tag_selection, index.first ,  evtweight, h_1d, 10, 0., 10.);
      plot1D("h_index_Pair1_NotgoodMatch"+tag_selection, index.second ,  evtweight, h_1d, 10, 0., 10.);

      plot1D("h_theta_Pair1_NotgoodMatch"+tag_selection, theta ,  evtweight, h_1d, 100, -1., 1.);

      plot1D("h_massHbb_Pair1_NotgoodMatch"+tag_selection, massHbb1 ,       evtweight, h_1d, 102, -10, 500);
      plot1D("h_deltaRap_Pair1_NotgoodMatch"+tag_selection, fabs(jetsP4.at(index.first).Rapidity()-jetsP4.at(index.second).Rapidity()) , evtweight, h_1d, 100, 0, 10.);
      plot1D("h_deltaPhi_Pair1_NotgoodMatch"+tag_selection, deltaPhi(jetsP4.at(index.first).phi(),jetsP4.at(index.second).phi()) ,  evtweight, h_1d, 100, 0, 10);

      plot1D("h_RapHbb_Pair1_NotgoodMatch"+tag_selection, (jetsP4.at(index.first)+jetsP4.at(index.second)).Rapidity() , evtweight, h_1d, 100, -5., 5.);

      plot1D("h_B1_mc3_pair_NotgoodMatch"+tag_selection,  abs(mc3.at(index.first))  ,       evtweight, h_1d, 16, 0, 16);
      plot1D("h_B2_mc3_pair_NotgoodMatch"+tag_selection,  abs(mc3.at(index.second))  ,       evtweight, h_1d, 16, 0, 16);

      plot1D("h_B1_mc1_pair_NotgoodMatch"+tag_selection,  abs(mc1.at(index.first))  ,       evtweight, h_1d, 30, 0, 30);
      plot1D("h_B2_mc1_pair_NotgoodMatch"+tag_selection,  abs(mc1.at(index.second))  ,       evtweight, h_1d, 30, 0, 30);
      
      plot1D("h_B1_yProd_pair_NotgoodMatch"+tag_selection,   jetsP4.at(index.first).Rapidity()*stopt.lep1().Rapidity(),       evtweight, h_1d, 100, -10, 10.);
      plot1D("h_B2_yProd_pair_NotgoodMatch"+tag_selection,   jetsP4.at(index.second).Rapidity()*stopt.lep1().Rapidity(),       evtweight, h_1d, 100, -10, 10.);
      
      plot1D("h_B1_pt_pair_NotgoodMatch"+tag_selection,  jetsP4.at(index.first).pt()  ,       evtweight, h_1d, 60, 0, 300);
      plot1D("h_B2_pt_pair_NotgoodMatch"+tag_selection,  jetsP4.at(index.second).pt()  ,       evtweight, h_1d, 60, 0, 300);
      
      plot1D("h_B1_eta_pair_NotgoodMatch"+tag_selection,  jetsP4.at(index.first).eta()  ,       evtweight, h_1d, 100, 5., 5.);
      plot1D("h_B2_eta_pair_NotgoodMatch"+tag_selection,  jetsP4.at(index.second).eta()  ,       evtweight, h_1d, 100, 5., 5.);

      /*      
      plot1D("h_B1_massJet_masspair_NotgoodMatch"+tag_selection,  jetsP4.at(index.first).M()/massHbb1  ,       evtweight, h_1d, 100, 0, 2);
      plot1D("h_B2_massJet_masspair_NotgoodMatch"+tag_selection,  jetsP4.at(index.second).M()/massHbb1  ,       evtweight, h_1d, 100, 0, 2);

      plot1D("h_B1_massJet_ptpair_NotgoodMatch"+tag_selection,  jetsP4.at(index.first).M()/jetsP4.at(index.first).pt()  ,       evtweight, h_1d, 100, 0, 2);
      plot1D("h_B2_massJet_ptpair_NotgoodMatch"+tag_selection,  jetsP4.at(index.second).M()/jetsP4.at(index.second).pt()  ,       evtweight, h_1d, 100, 0, 2);
      */

    }

    if(goodMatch1 || goodMatchZ1) {
      
      plot1D("h_index_Pair1_goodMatch"+tag_selection, index.first ,  evtweight, h_1d, 10, 0., 10.);
      plot1D("h_index_Pair1_goodMatch"+tag_selection, index.second ,  evtweight, h_1d, 10, 0., 10.);

      plot1D("h_theta_Pair1_goodMatch"+tag_selection, theta ,  evtweight, h_1d, 100, -1., 1.);

      plot1D("h_massHbb_Pair1_goodMatch"+tag_selection, massHbb1 ,       evtweight, h_1d, 102, -10, 500);
      plot1D("h_deltaRap_Pair1_goodMatch"+tag_selection, fabs(jetsP4.at(index.first).Rapidity()-jetsP4.at(index.second).Rapidity()) , evtweight, h_1d, 100, 0, 10);
      plot1D("h_deltaPhi_Pair1_goodMatch"+tag_selection, deltaPhi(jetsP4.at(index.first).phi(),jetsP4.at(index.second).phi()) ,  evtweight, h_1d, 100, 0, 10);
      plot1D("h_RapHbb_Pair1_goodMatch"+tag_selection, (jetsP4.at(index.first)+jetsP4.at(index.second)).Rapidity() , evtweight, h_1d, 100, -5., 5.);

      /*
      plot1D("h_B1_massJet_masspair_goodMatch"+tag_selection,  jetsP4.at(index.first).M()/massHbb1  ,       evtweight, h_1d, 100, 0, 2);
      plot1D("h_B2_massJet_masspair_goodMatch"+tag_selection,  jetsP4.at(index.second).M()/massHbb1  ,       evtweight, h_1d, 100, 0, 2);

      plot1D("h_B1_massJet_ptpair_goodMatch"+tag_selection,  jetsP4.at(index.first).M()/jetsP4.at(index.first).pt()  ,       evtweight, h_1d, 100, 0, 2);
      plot1D("h_B2_massJet_ptpair_goodMatch"+tag_selection,  jetsP4.at(index.second).M()/jetsP4.at(index.second).pt()  ,       evtweight, h_1d, 100, 0, 2);
      */

      plot1D("h_B1_pt_pair_goodMatch"+tag_selection,  jetsP4.at(index.first).pt()  ,       evtweight, h_1d, 60, 0, 300);
      plot1D("h_B2_pt_pair_goodMatch"+tag_selection,  jetsP4.at(index.second).pt()  ,       evtweight, h_1d, 60, 0, 300);

      plot1D("h_B1_eta_pair_goodMatch"+tag_selection,  jetsP4.at(index.first).eta()  ,       evtweight, h_1d, 100, 5., 5.);
      plot1D("h_B2_eta_pair_goodMatch"+tag_selection,  jetsP4.at(index.second).eta()  ,       evtweight, h_1d, 100, 5., 5.);

      plot1D("h_B1_yProd_pair_goodMatch"+tag_selection,   jetsP4.at(index.first).Rapidity()*stopt.lep1().Rapidity(),       evtweight, h_1d, 100, -10, 10.);
      plot1D("h_B2_yProd_pair_goodMatch"+tag_selection,   jetsP4.at(index.second).Rapidity()*stopt.lep1().Rapidity(),       evtweight, h_1d, 100, -10, 10.);

    }

    float pt1 = (jetsP4.at(index.first)+jetsP4.at(index.second)).pt();
    if(goodMatch1) plot1D("h_pt1_goodMatch"+tag_selection,  pt1  ,       evtweight, h_1d, 50, 0, 500);
    if(!goodMatch1) plot1D("h_pt1_nongoodMatch"+tag_selection,  pt1 ,       evtweight, h_1d, 50, 0, 500);

    if(goodMatch1) plot1D("h_deltaPhi_H1_lep1_goodMatch"+tag_selection, deltaPhi((jetsP4.at(index.first)+jetsP4.at(index.second)).phi(),stopt.lep1().phi()) ,  evtweight, h_1d, 100, 0, 10);
    if(!goodMatch1) plot1D("h_deltaPhi_H1_lep1_nongoodMatch"+tag_selection, deltaPhi((jetsP4.at(index.first)+jetsP4.at(index.second)).phi(),stopt.lep1().phi()) ,  evtweight, h_1d, 100, 0, 10);

    float ptdilepton=(stopt.lep1()+stopt.lep2()).pt();
    plot1D("h_ptll"+tag_selection,  ptdilepton  ,       evtweight, h_1d, 50, 0, 500);

    if(goodMatch1) plot1D("h_ptbb_DileVec_diff_goodMatch"+tag_selection,  (jetsP4.at(index.first)+jetsP4.at(index.second)+stopt.lep1()+stopt.lep2()).pt()  , evtweight, h_1d, 100, -100, 1000);
    if(!goodMatch1) plot1D("h_ptbb_DileVec_diff_nongoodMatch"+tag_selection,  (jetsP4.at(index.first)+jetsP4.at(index.second)+stopt.lep1()+stopt.lep2()).pt()  , evtweight, h_1d, 100, -100, 1000);


    pair<int, int> index_Pair2=getIndexPair(indexLooseB40,jetsP4,true,index.first,index.second);    

    if(index_Pair2.first!=(-1.)) {

      /*
    if(indexJetsP4.size()==3) {

      cout << "3 b events found a second PAIR ?????"<< endl;

      cout << "indexPair1 are: " << index.first << "" << index.second << endl;
      cout << "indexPair2 are: " << index_Pair2.first << "" << index_Pair2.second << endl;

    }
      */

      float massHbb2 = (jetsP4.at(index_Pair2.first)+jetsP4.at(index_Pair2.second)).M();

      float pt2 = (jetsP4.at(index_Pair2.first)+jetsP4.at(index_Pair2.second)).pt();

      //      float average = fabs(pt2-pt1)/max(pt2,pt1);

      if((mc3.at(index_Pair2.first)==7 && mc3.at(index_Pair2.second)==-7) || (mc3.at(index_Pair2.first)==-7 && mc3.at(index_Pair2.second)==7) ) {
      //      if(abs(mc3.at(index_Pair2.first))==7 && abs(mc3.at(index_Pair2.second))==7) {

	//	plot1D("h_averagePtPair_goodMatch"+tag_selection, average ,       evtweight, h_1d, 100, 0, 2);
	//	if(goodMatch1) plot1D("h_pt2_goodMatch"+tag_selection,  pt2  ,       evtweight, h_1d, 50, 0, 500);

	plot1D("h_index_Pair2_goodMatch"+tag_selection, index_Pair2.first ,  evtweight, h_1d, 10, 0., 10.);
	plot1D("h_index_Pair2_goodMatch"+tag_selection, index_Pair2.second ,  evtweight, h_1d, 10, 0., 10.);

	plot1D("h_massHbb_Pair2_goodMatch"+tag_selection, massHbb2 ,       evtweight, h_1d, 102, -10, 500);
	plot1D("h_deltaRap_Pair2_goodMatch"+tag_selection, fabs(jetsP4.at(index_Pair2.first).Rapidity()-jetsP4.at(index_Pair2.second).Rapidity()) , evtweight, h_1d, 100, 0, 10);
	plot1D("h_deltaPhi_Pair2_goodMatch"+tag_selection, deltaPhi(jetsP4.at(index_Pair2.first).phi(),jetsP4.at(index_Pair2.second).phi()) ,  evtweight, h_1d, 100, 0, 10);

	plot1D("h_RapHbb_Pair2_goodMatch"+tag_selection, (jetsP4.at(index_Pair2.first)+jetsP4.at(index_Pair2.second)).Rapidity() , evtweight, h_1d, 100, -5., 5.);

	if(goodMatch1) plot1D("h_deltaPhi_H2_lep1_goodMatch"+tag_selection, deltaPhi((jetsP4.at(index_Pair2.first)+jetsP4.at(index_Pair2.second)).phi(),stopt.lep1().phi()) ,  evtweight, h_1d, 100, 0, 10);
	if(!goodMatch1) plot1D("h_deltaPhi_H2_lep1_nongoodMatch"+tag_selection, deltaPhi((jetsP4.at(index_Pair2.first)+jetsP4.at(index_Pair2.second)).phi(),stopt.lep1().phi()) ,  evtweight, h_1d, 100, 0, 10);

	plot2D("h_mass1_vs_mass2_goodMatch"+tag_selection, massHbb1 , massHbb2,    evtweight, h_2d, 50, 0, 500, 50, 0, 500);
	plot2D("h_mass1_vs_DMmass12_goodMatch"+tag_selection, massHbb1 , massHbb2-massHbb1,    evtweight, h_2d, 50, 0, 500, 50, 0, 500);

      } else {

	//	plot1D("h_averagePtPair_notgoodMatch"+tag_selection, average ,       evtweight, h_1d, 100, 0, 2);
	//	if(!goodMatch1) plot1D("h_pt2_nongoodMatch"+tag_selection,  pt2  ,       evtweight, h_1d, 50, 0, 500);

	plot1D("h_massHbb_Pair2_NotgoodMatch"+tag_selection, massHbb2 ,       evtweight, h_1d, 102, -10, 500);
	plot1D("h_deltaRap_Pair2_NotgoodMatch"+tag_selection, fabs(jetsP4.at(index_Pair2.first).Rapidity()-jetsP4.at(index_Pair2.second).Rapidity()) , evtweight, h_1d, 100, 0, 10);
	plot1D("h_deltaPhi_Pair2_NotgoodMatch"+tag_selection, deltaPhi(jetsP4.at(index_Pair2.first).phi(),jetsP4.at(index_Pair2.second).phi()) ,  evtweight, h_1d, 100, 0, 10);

	plot2D("h_mass1_vs_mass2_NotgoodMatch"+tag_selection, massHbb1 , massHbb2,    evtweight, h_2d, 50, 0, 500, 50, 0, 500);
	plot2D("h_mass1_vs_DMmass12_NotgoodMatch"+tag_selection, massHbb1 , massHbb2-massHbb1,    evtweight, h_2d, 50, 0, 500, 50, 0, 500);

      }
    }
  }


  ////// not cleaned

  if(indexOS_noClean.first!=(-1.)) {
    
    float deltaRHbb_noclean = deltaR(jetsP4.at(indexOS_noClean.first).eta(),jetsP4.at(indexOS_noClean.first).phi(),jetsP4.at(indexOS_noClean.second).eta(),jetsP4.at(indexOS_noClean.second).phi());
    float massHbb_noclean = (jetsP4.at(indexOS_noClean.first)+jetsP4.at(indexOS_noClean.second)).M();
    float ptHbb_noclean = (jetsP4.at(indexOS_noClean.first)+jetsP4.at(indexOS_noClean.second)).pt();
    float deltaRap_noclean = fabs(jetsP4.at(indexOS_noClean.first).Rapidity()-jetsP4.at(indexOS_noClean.second).Rapidity());
    float deltaPhibb_noclean = deltaPhi(jetsP4.at(indexOS_noClean.first).phi(),jetsP4.at(indexOS_noClean.second).phi());

    float deltaRapHb1_noclean = fabs((jetsP4.at(indexOS_noClean.first)+jetsP4.at(indexOS_noClean.second)).Rapidity()-jetsP4.at(indexOS_noClean.first).Rapidity());
    float deltaRapHb2_noclean = fabs((jetsP4.at(indexOS_noClean.first)+jetsP4.at(indexOS_noClean.second)).Rapidity()-jetsP4.at(indexOS_noClean.second).Rapidity());

    bool goodMatch1 = ((mc3.at(indexOS_noClean.first)==7 && mc3.at(indexOS_noClean.second)==-7) || (mc3.at(indexOS_noClean.first)==-7 && mc3.at(indexOS_noClean.second)==7));

    float dropPt=min(jetsP4.at(indexOS_noClean.first).pt(),jetsP4.at(indexOS_noClean.second).pt())/((jetsP4.at(indexOS_noClean.first)+jetsP4.at(indexOS_noClean.second)).pt());
    float dropAngle=(massHbb_noclean/ptHbb_noclean)/deltaRHbb_noclean;
    
    //    if(massHbb_noclean > 100 && massHbb_noclean<150) {
      
    float minInvPt=min(1/(jetsP4.at(indexOS_noClean.first).pt()*jetsP4.at(indexOS_noClean.first).pt()), 1/(jetsP4.at(indexOS_noClean.second).pt()*jetsP4.at(indexOS_noClean.second).pt()));
    float dropKt=minInvPt*(deltaRHbb_noclean*deltaRHbb_noclean);                                             
    
    if(goodMatch1) {
      
      //              plot2D("h_DR_vs_ptbb", ptHbb_noclean , deltaRHbb,    evtweight, h_2d, 50, 0, 500, 100, 0, 6);                                                       
      
      plot1D("h_minInvPt_noclean_goodMatch"+tag_selection,  minInvPt  ,       evtweight, h_1d, 100, 0, 0.1);                                                                       
      plot1D("h_deltaR_noclean_goodMatch"+tag_selection,  deltaRHbb_noclean  ,       evtweight, h_1d, 100, 0, 5);
      plot1D("h_dropKt_noclean_goodMatch"+tag_selection,  dropKt  ,       evtweight, h_1d, 100, 0, 10);
      plot1D("h_deltaRap_noclean_goodMatch"+tag_selection,  deltaRap_noclean  ,       evtweight, h_1d, 100, 0, 5.);
      
      plot1D("h_deltaRapHb1_noclean_goodMatch"+tag_selection,  deltaRapHb1_noclean  ,       evtweight, h_1d, 100, 0, 5.);
      plot1D("h_deltaRapHb2_noclean_goodMatch"+tag_selection,  deltaRapHb2_noclean  ,       evtweight, h_1d, 100, 0, 5.);
      
      plot1D("h_minJetPt_noclean_goodMatch"+tag_selection,  min(jetsP4.at(indexOS_noClean.first).pt(),jetsP4.at(indexOS_noClean.second).pt())  ,       evtweight, h_1d, 100, 0, 500);
      plot1D("h_ptbb_noclean_goodMatch"+tag_selection,  ptHbb_noclean  ,       evtweight, h_1d, 100, 0, 500);
      plot1D("h_deltaPhi_noclean_goodMatch"+tag_selection,  deltaPhibb_noclean ,       evtweight, h_1d, 100, 0, 4);
      
      plot1D("h_dropPt_noclean_goodMatch"+tag_selection,  dropPt  ,       evtweight, h_1d, 100, 0, 2);
      
      plot1D("h_dropAngle_noclean_goodMatch"+tag_selection,  dropAngle  ,       evtweight, h_1d, 100, 0, 10);
      plot2D("h_dR_vs_Pt_noclean_goodMatch"+tag_selection,  ptHbb_noclean, deltaRHbb_noclean  ,       evtweight, h_2d, 100, 0, 500, 100, 0, 10);
      plot2D("h_dropAngle_vs_Pt_noclean_goodMatch"+tag_selection,  ptHbb_noclean, dropAngle  ,       evtweight, h_2d, 100, 0, 500, 100, 0, 2.5);
      
    } else {
      
	plot1D("h_minInvPt_noclean_NotgoodMatch"+tag_selection,  minInvPt  ,       evtweight, h_1d, 100, 0, 0.1);
	plot1D("h_dropKt_noclean_NotgoodMatch"+tag_selection,  dropKt  ,       evtweight, h_1d, 100, 0, 10);
	
	plot1D("h_deltaRapHb1_noclean_NotgoodMatch"+tag_selection,  deltaRapHb1_noclean  ,       evtweight, h_1d, 100, 0, 5.);
	plot1D("h_deltaRapHb2_noclean_NotgoodMatch"+tag_selection,  deltaRapHb2_noclean  ,       evtweight, h_1d, 100, 0, 5.);

	plot1D("h_deltaR_noclean_NotgoodMatch"+tag_selection,  deltaRHbb_noclean  ,       evtweight, h_1d, 100, 0, 5.);
	plot1D("h_deltaRap_noclean_NotgoodMatch"+tag_selection,  deltaRap_noclean  ,       evtweight, h_1d, 100, 0, 10);
	plot1D("h_minJetPt_noclean_NotgoodMatch"+tag_selection,  min(jetsP4.at(indexOS_noClean.first).pt(),jetsP4.at(indexOS_noClean.second).pt())  ,       evtweight, h_1d, 100, 0, 500);
	plot1D("h_ptMax_noclean_NotgoodMatch"+tag_selection,  ptHbb_noclean  ,       evtweight, h_1d, 100, 0, 500);
	plot1D("h_deltaPhibb_noclean_NotgoodMatch"+tag_selection,  deltaPhibb_noclean ,       evtweight, h_1d, 100, 0, 4);
	
	plot1D("h_dropPt_noclean_NotgoodMatch"+tag_selection,  dropPt  ,       evtweight, h_1d, 100, 0, 2);
	
	plot1D("h_dropAngle_noclean_NotgoodMatch"+tag_selection,  dropAngle  ,       evtweight, h_1d, 100, 0, 10);
	plot2D("h_dR_vs_Pt_noclean_NotgoodMatch"+tag_selection,  ptHbb_noclean, deltaRHbb_noclean  ,       evtweight, h_2d, 100, 0, 500, 100, 0, 10);
	plot2D("h_dropAngle_vs_Pt_noclean_NotgoodMatch"+tag_selection,  ptHbb_noclean, dropAngle  ,       evtweight, h_2d, 100, 0, 500, 100, 0, 2.5);
	
      }
      //    }
    
  }
  
}

void StopTreeLooper::plotCR4(float evtweight, string tag_selection , std::map<std::string, TH1F*> &h_1d, bool isData, bool doLRM, bool doTAU) {

  vector<int> i_LooseB40; 
  vector<int> i_LooseB30; 
  vector<int> i_TightB40; 

  if(doTAU) i_LooseB30 = indexLooseB30Tau;
  if(doTAU) i_LooseB40 = indexLooseB40Tau;
  if(doTAU) i_TightB40 = indexTightB40Tau;

  if(!doTAU) i_LooseB30 = indexLooseB30;
  if(!doTAU) i_LooseB40 = indexLooseB40;
  if(!doTAU) i_TightB40 = indexTightB40;

  TString name(tag_selection);

  if(name.Contains("4bM")) {

    plot1D("cr4_nbjets"+tag_selection,   i_LooseB30.size() ,       evtweight, h_1d, 10,  1., 11.);

    for( unsigned int i = 0 ; i < i_LooseB30.size() ; ++i ){
      
      //      if( doLRM && (lrm.at(i) < passLRM(jetsP4.at(i).pt()))) continue;

      plot1D("cr4_eta_bjets"+tag_selection,  jetsP4.at(i_LooseB30.at(i)).eta() ,       evtweight, h_1d, 100,  -5., 5.);
      plot1D("cr4_phi_bjets"+tag_selection,  jetsP4.at(i_LooseB30.at(i)).phi() ,       evtweight, h_1d, 100,  -5., 5.);
      plot1D("cr4_pt_bjets"+tag_selection,  jetsP4.at(i_LooseB30.at(i)).pt() ,       evtweight, h_1d, 100,  0., 300.);
      plot1D("cr4_discr_bjets"+tag_selection,  csv.at(i_LooseB30.at(i)),       evtweight, h_1d, 100,  0., 1.);
      plot1D("cr4_lrm_bjets"+tag_selection,  lrm.at(i_LooseB30.at(i)),       evtweight, h_1d, 100,  0., 1.);
      plot1D("cr4_neu_bjets"+tag_selection,  neu.at(i_LooseB30.at(i)),       evtweight, h_1d, 100,  0., 100.);
      plot1D("cr4_chm_bjets"+tag_selection,  chm.at(i_LooseB30.at(i)),       evtweight, h_1d, 100,  0., 100.);

      float yProd= stopt.lep1().Rapidity()*jetsP4.at(i_LooseB30.at(i)).Rapidity() ;
      plot1D("cr4_yProd"+tag_selection,  yProd ,       evtweight, h_1d, 100, -5., 5.);

      float dphi_bjetMet=deltaPhi(jetsP4.at(i_LooseB30.at(i)).phi(), stopt.t1metphicorrphi());
      plot1D("cr4_bjetsMetdphi"+tag_selection,  dphi_bjetMet ,       evtweight, h_1d, 100, 0, TMath::Pi());

    }

  }

  if(name.Contains("3bM")) {

    plot1D("cr4_nbjets"+tag_selection,   i_LooseB40.size() ,       evtweight, h_1d, 10,  1., 11.);

    for( unsigned int i = 0 ; i < i_LooseB40.size() ; ++i ){
      
      //      if( doLRM && (lrm.at(i) < passLRM(jetsP4.at(i).pt()))) continue;

      plot1D("cr4_eta_bjets"+tag_selection,  jetsP4.at(i_LooseB40.at(i)).eta() ,       evtweight, h_1d, 100,  -5., 5.);
      plot1D("cr4_phi_bjets"+tag_selection,  jetsP4.at(i_LooseB40.at(i)).phi() ,       evtweight, h_1d, 100,  -5., 5.);
      plot1D("cr4_pt_bjets"+tag_selection,  jetsP4.at(i_LooseB40.at(i)).pt() ,       evtweight, h_1d, 100,  0., 300.);
      plot1D("cr4_discr_bjets"+tag_selection,  csv.at(i_LooseB40.at(i)),       evtweight, h_1d, 100,  0., 1.);
      plot1D("cr4_lrm_bjets"+tag_selection,  lrm.at(i_LooseB40.at(i)),       evtweight, h_1d, 100,  0., 1.);
      plot1D("cr4_neu_bjets"+tag_selection,  neu.at(i_LooseB40.at(i)),       evtweight, h_1d, 100,  0., 100.);
      plot1D("cr4_chm_bjets"+tag_selection,  chm.at(i_LooseB40.at(i)),       evtweight, h_1d, 100,  0., 10.);

      float yProd= stopt.lep1().Rapidity()*jetsP4.at(i_LooseB40.at(i)).Rapidity() ;
      plot1D("cr4_yProd"+tag_selection,  yProd ,       evtweight, h_1d, 100, -5., 5.);

      float dphi_bjetMet=deltaPhi(jetsP4.at(i_LooseB40.at(i)).phi(), stopt.t1metphicorrphi());
      plot1D("cr4_bjetsMetdphi"+tag_selection,  dphi_bjetMet ,       evtweight, h_1d, 100, 0, TMath::Pi());

    }
  }


  if(name.Contains("3bT")) {

    plot1D("cr4_nbjets"+tag_selection,   i_TightB40.size() ,       evtweight, h_1d, 10,  1., 11.);

    for( unsigned int i = 0 ; i < i_TightB40.size() ; ++i ){
      
      plot1D("cr4_eta_bjets"+tag_selection,  jetsP4.at(i_TightB40.at(i)).eta() ,       evtweight, h_1d, 100,  -5., 5.);
      plot1D("cr4_phi_bjets"+tag_selection,  jetsP4.at(i_TightB40.at(i)).phi() ,       evtweight, h_1d, 100,  -5., 5.);
      plot1D("cr4_pt_bjets"+tag_selection,  jetsP4.at(i_TightB40.at(i)).pt() ,       evtweight, h_1d, 100,  0., 300.);
      plot1D("cr4_discr_bjets"+tag_selection,  csv.at(i_TightB40.at(i)),       evtweight, h_1d, 100,  0., 1.);

      plot1D("cr4_lrm_bjets"+tag_selection,  lrm.at(i_TightB40.at(i)),       evtweight, h_1d, 100,  0., 1.);
      plot1D("cr4_neu_bjets"+tag_selection,  neu.at(i_TightB40.at(i)),       evtweight, h_1d, 100,  0., 1.);
      plot1D("cr4_chm_bjets"+tag_selection,  chm.at(i_TightB40.at(i)),       evtweight, h_1d, 100,  0., 1.);

      float yProd= stopt.lep1().Rapidity()*jetsP4.at(i_TightB40.at(i)).Rapidity() ;
      plot1D("cr4_yProd"+tag_selection,  yProd ,       evtweight, h_1d, 100, -5., 5.);

      float dphi_bjetMet=deltaPhi(jetsP4.at(i_TightB40.at(i)).phi(), stopt.t1metphicorrphi());
      plot1D("cr4_bjetsMetdphi"+tag_selection,  dphi_bjetMet ,       evtweight, h_1d, 100, 0, TMath::Pi());
      
    }
  }

  plot1D("cr4_mt"+tag_selection,   stopt.mt() ,       evtweight, h_1d, 100,  0, 300);

  if(doTAU) {
    float mtpftau = getMT( stopt.pfTau().pt() , stopt.pfTau().phi() , stopt.t1metphicorr() , stopt.t1metphicorrphi());
    plot1D("cr4_mtpftau"+tag_selection,   mtpftau ,       evtweight, h_1d, 100,  0, 300);
  }
  

  plot1D("cr4_met"+tag_selection,   stopt.t1metphicorr() ,       evtweight, h_1d, 100, 0, 300);  
  if(!doTAU) plot1D("cr4_njets"+tag_selection,  jetsP4.size() ,       evtweight, h_1d, 10,  1., 11.);
  if(doTAU) plot1D("cr4_njets"+tag_selection,  n_taujets ,       evtweight, h_1d, 10,  1., 11.);

  //// plotCR4

  if(doTAU && (stopt.pfTau_leadPtcandID()!=-1)) {
    
    plot1D("cr4_pfTau_pt"+tag_selection,  stopt.pfTau().pt() ,       evtweight, h_1d, 150, 0, 150.);
    plot1D("cr4_pfTau_eta"+tag_selection,  stopt.pfTau().eta() ,       evtweight, h_1d, 100, -3, 3.);  
    
    float mlt=(stopt.lep1()+stopt.pfTau()).M();
    plot1D("cr4_tauLep_m"+tag_selection,  mlt,       evtweight, h_1d, 100, 0, 300.);
    
    float drlt=deltaR(stopt.lep1().eta() , stopt.lep1().phi() , stopt.pfTau().eta(), stopt.pfTau().phi());
    plot1D("cr4_tauLep_dr"+tag_selection,  drlt ,       evtweight, h_1d, 100, 0, 6.);

    float ptlt=(stopt.lep1()+stopt.pfTau()).pt();
    plot1D("cr4_tauLep_pt"+tag_selection,  ptlt ,       evtweight, h_1d, 100, 0, 500.);

  }

  //// plotCR4

  if(!doTAU) {

    float pt12=(stopt.lep1()+stopt.lep2()).pt();
    float m12=(stopt.lep1()+stopt.lep2()).M();
    float Rap12=(stopt.lep1()+stopt.lep2()).Rapidity();
    
    float dr=deltaR(stopt.lep1().eta() , stopt.lep1().phi() , stopt.lep2().eta(), stopt.lep2().phi());
    float dphi=deltaPhi(stopt.lep1().phi(), stopt.lep2().phi());
    float dphiL1Met=deltaPhi(stopt.lep1().phi(), stopt.t1metphicorrphi());
    float dphiL2Met=deltaPhi(stopt.lep2().phi(), stopt.t1metphicorrphi());

    float mt2_standard = MT2( stopt.t1metphicorr(), stopt.t1metphicorrphi(), stopt.lep1() , stopt.lep2() );
  
    plot1D("cr4_dilep_mt2"+tag_selection,  mt2_standard ,       evtweight, h_1d, 100, 0, 300.);
    plot1D("cr4_dilep_m"+tag_selection,  m12 ,       evtweight, h_1d, 100, 0, 300.);
    plot1D("cr4_dilep_pt"+tag_selection,  pt12 ,       evtweight, h_1d, 100, 0, 500.);
    plot1D("cr4_dilep_rap"+tag_selection,  Rap12 ,       evtweight, h_1d, 100, -3., 3.);
    plot1D("cr4_dilep_dr"+tag_selection,  dr ,       evtweight, h_1d, 100, 0, 6.);
    plot1D("cr4_dilep_dphi"+tag_selection,  dphi ,       evtweight, h_1d, 100, 0., TMath::Pi());

    plot1D("cr4_lep2_pt"+tag_selection,  stopt.lep2().pt() ,       evtweight, h_1d, 100, 0, 100.);
    plot1D("cr4_lep2_eta"+tag_selection,  stopt.lep2().eta() ,       evtweight, h_1d, 100, -3, 3.);

    plot1D("cr4_lep1Met"+tag_selection,  dphiL1Met ,       evtweight, h_1d, 100, 0, TMath::Pi());
    plot1D("cr4_lep2Met"+tag_selection,  dphiL2Met ,       evtweight, h_1d, 100, 0, TMath::Pi());
    
  }  

  plot1D("cr4_lep1_pt"+tag_selection,  stopt.lep1().pt() ,       evtweight, h_1d, 150, 0, 150.);
  plot1D("cr4_lep1_eta"+tag_selection,  stopt.lep1().eta() ,       evtweight, h_1d, 100, -3, 3.);  

  if(!doTAU) plot1D("cr4_massHbb"+tag_selection,  massHbb ,       evtweight, h_1d, 51, 0, 500);
  if(doTAU) plot1D("cr4_massHbb"+tag_selection,  massHbbTau ,       evtweight, h_1d, 51, 0, 500);


  /*

  ///// ??????????  need to add for the non b jets
 

    }

    if (csv.at(i) < 0.240) {
      if(stopt.id2()!=(-1)) {
	if(deltaR(stopt.pfjets().at(i).eta() , stopt.pfjets().at(i).phi() , stopt.lep2().eta(), stopt.lep2().phi())<0.4) continue; 
     
	plot1D("cr4_eta_nonbjets"+tag_selection,  jetsP4.at(i).eta() ,       evtweight, h_1d, 100,  -5., 5.);
	plot1D("cr4_discr_nonbjets"+tag_selection,  csv.at(i),       evtweight, h_1d, 100,  -1.5, 1.);
	plot1D("cr4_lrm_nonbjets"+tag_selection,  lrm.at(i),       evtweight, h_1d, 100,  0., 0.5);
	plot1D("cr4_chm_nonbjets"+tag_selection,  chm.at(i),       evtweight, h_1d, 100,  0., 100.);
	plot1D("cr4_neu_nonbjets"+tag_selection,  neu.at(i),       evtweight, h_1d, 100,  0., 100.);

	float dphinonbjetMet=deltaPhi(jetsP4.at(i).phi(), stopt.t1metphicorrphi());
	plot1D("cr4_nonbjetsMet"+tag_selection,  dphinonbjetMet ,       evtweight, h_1d, 100, 0, TMath::Pi());

      }
    }

    //  if(isData)  cout << "RUN " << stopt.run() << " lumi " << stopt.lumi() << " event " << stopt.event() << endl;

  }
*/

}

void StopTreeLooper::plotCR1(float evtweight, string tag_selection , std::map<std::string, TH1F*> &h_1d, bool isData) {

  //  if(detaFB!=-9999.) plot1D("cr1_detaFB"+tag_selection,  detaFB  ,       evtweight, h_1d, 100,  0, 10.);

  TString name(tag_selection);

  if(name.Contains("4bM")) {

    plot1D("cr1_nbjets"+tag_selection, indexLooseB30.size()   ,       evtweight, h_1d, 10,  1., 11.);

    for( unsigned int i = 0 ; i < indexLooseB30.size() ; ++i ){
      
      plot1D("cr1_eta_bjets"+tag_selection,  jetsP4.at(indexLooseB30.at(i)).eta() ,       evtweight, h_1d, 100,  -5., 5.);
      plot1D("cr1_pt_bjets"+tag_selection,  jetsP4.at(indexLooseB30.at(i)).pt() ,       evtweight, h_1d, 100,  0., 300.);
      plot1D("cr1_discr_bjets"+tag_selection,  csv.at(indexLooseB30.at(i)),       evtweight, h_1d, 100,  0., 1.);
    }
  }


  if(name.Contains("3bM")) {

    plot1D("cr1_nbjets"+tag_selection, indexLooseB40.size()   ,       evtweight, h_1d, 10,  1., 11.);

    for( unsigned int i = 0 ; i < indexLooseB40.size() ; ++i ){
      
      plot1D("cr1_eta_bjets"+tag_selection,  jetsP4.at(indexLooseB40.at(i)).eta() ,       evtweight, h_1d, 100,  -5., 5.);
      plot1D("cr1_pt_bjets"+tag_selection,  jetsP4.at(indexLooseB40.at(i)).pt() ,       evtweight, h_1d, 100,  0., 300.);
      plot1D("cr1_discr_bjets"+tag_selection,  csv.at(indexLooseB40.at(i)),       evtweight, h_1d, 100,  0., 1.);
    }
  }

  if(name.Contains("3bT")) {

    plot1D("cr1_nbjets"+tag_selection, indexTightB40.size()   ,       evtweight, h_1d, 10,  1., 11.);

    for( unsigned int i = 0 ; i < indexTightB40.size() ; ++i ){
      
      plot1D("cr1_eta_bjets"+tag_selection,  jetsP4.at(indexTightB40.at(i)).eta() ,       evtweight, h_1d, 100,  -5., 5.);
      plot1D("cr1_pt_bjets"+tag_selection,  jetsP4.at(indexTightB40.at(i)).pt() ,       evtweight, h_1d, 100,  0., 300.);
      plot1D("cr1_discr_bjets"+tag_selection,  csv.at(indexTightB40.at(i)),       evtweight, h_1d, 100,  0., 1.);
    }
  }
 

  plot1D("cr1_met"+tag_selection,   stopt.t1metphicorr() ,       evtweight, h_1d, 100, 0, 300);
  plot1D("cr1_mt"+tag_selection,   stopt.mt() ,       evtweight, h_1d, 100,  0, 300);
  plot1D("cr1_massHbb"+tag_selection,   massHbb ,       evtweight, h_1d,  102, -10, 500);

  plot1D("cr1_njets"+tag_selection,  jetsP4.size() ,       evtweight, h_1d, 10,  1., 11.);

  plot1D("cr1_lep1_pt"+tag_selection,  stopt.lep1().pt() ,       evtweight, h_1d, 150, 0, 150.);
  plot1D("cr1_lep1_eta"+tag_selection,  stopt.lep1().eta() ,       evtweight, h_1d, 100, -3, 3.);


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
                                                                                                                                                                  
  plot1D("cr5_iso"+tag_selection+tag_flav+tag_charge, stopt.pfcandisoOS10looseZ()*stopt.pfcandOS10looseZ().pt() ,           evtweight, h_1d, 100, 0., 10.);                                                     
  plot1D("cr5_isoRel"+tag_selection+tag_flav+tag_charge, stopt.pfcandisoOS10looseZ() ,           evtweight, h_1d, 100, 0., 1.);                                                     
  plot1D("cr5_dz"+tag_selection+tag_flav+tag_charge, stopt.pfcanddzOS10looseZ() ,           evtweight, h_1d, 100,  0., 0.1);                                                     
  plot1D("cr5_pt"+tag_selection+tag_flav+tag_charge, stopt.pfcandOS10looseZ().pt() ,           evtweight, h_1d, 20,  0., 100.);                                                     
  plot1D("cr5_eta"+tag_selection+tag_flav+tag_charge, stopt.pfcandOS10looseZ().eta() ,           evtweight, h_1d, 100,  -3., 3.);                                                     

                                                                                                                                                                  
}     
