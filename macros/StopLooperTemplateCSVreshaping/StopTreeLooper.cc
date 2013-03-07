#include "StopTreeLooper.h"

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
#include "../../Tools/BTagReshaping/BTagReshaping.h"
#include "../looper/BtagFuncs.h"
#include "TRandom3.h"

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


  TRandom3 *random3_ = new TRandom3(1);
  //------------------------------------------------------------------------------------------------------
  // set csv discriminator reshaping
  //------------------------------------------------------------------------------------------------------
  
  BTagShapeInterface * nominalShape=0;

  //alternative reshaping for systematics
  //BTagShapeInterface * downBCShape=0;
  //BTagShapeInterface * upBCShape=0;
  //BTagShapeInterface * downLShape=0;
  //BTagShapeInterface * upLShape=0;  
  //BTagShapeInterface * nominalShape4p=0;
  //BTagShapeInterface * downBCShape4p=0;
  //BTagShapeInterface * upBCShape4p=0;

  nominalShape = new BTagShapeInterface("../../Tools/BTagReshaping/csvdiscr.root", 0.0, 0.0); 

  //alternative reshaping for systematics
  //upBCShape = new BTagShapeInterface("../../Tools/BTagReshaping/csvdiscr.root", 1.0, 0.0); 
  //downBCShape = new BTagShapeInterface("../../Tools/BTagReshaping/csvdiscr.root", -1.0, 0.0); 
  //upLShape = new BTagShapeInterface("../../Tools/BTagReshaping/csvdiscr.root", 0.0, 1.0); 
  //downLShape = new BTagShapeInterface("../../Tools/BTagReshaping/csvdiscr.root", 0.0, -1.0); 
  //nominalShape4p = new BTagShapeInterface("../../Tools/BTagReshaping/csvdiscr.root", 0.0, 0.0,true,1.003,1.003); 
  //upBCShape4p = new BTagShapeInterface("../../Tools/BTagReshaping/csvdiscr.root", 0.0, 0.0,true,1.001,1.001); 
  //downBCShape4p = new BTagShapeInterface("../../Tools/BTagReshaping/csvdiscr.root", 0.0, 0.0,true,1.005,1.005); 


  //------------------------------
  // set up histograms
  //------------------------------

  gROOT->cd();

  cout << "[StopTreeLooper::loop] setting up histos" << endl;

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

      //---------------------------------------------------------------------------- 
      // determine event weight
      // make 2 example histograms of nvtx and corresponding weight
      //---------------------------------------------------------------------------- 

      float evtweight = isData ? 1. : ( stopt.weight() * 19.5 * stopt.nvtxweight() * stopt.mgcor() );
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
      // Do the csv discriminator reshaping
      //----------------------------------------------------------------------------

      int nbtagscsvm_default = 0; //count default csvm tags to cross-check with results from looper
      int nbtagscsvm_reshaped = 0; //number of tags after reshaping

      int nbtagscsvm_random = 0;  //number of csvm tags after random upgrading/downgrading, using default SFs
      int nbtagscsvm_random_new = 0; //number of csvm tags after random upgrading/downgrading, using the same SFs as for reshaping

      for (unsigned int ijet = 0 ; ijet < stopt.pfjets_csv().size() ; ijet++) {

        if(       stopt.pfjets()[ijet].Pt()    < 30. )           continue;
        if( fabs( stopt.pfjets()[ijet].Eta() ) > 2.4 )           continue;


        float csv_nominal=nominalShape->reshape(stopt.pfjets()[ijet].Eta(),stopt.pfjets()[ijet].Pt(),stopt.pfjets_csv()[ijet],stopt.pfjets_mcflavorAlgo()[ijet]);
        //float csv_nominal=nominalShape->reshape(stopt.pfjets()[ijet].Eta(),stopt.pfjets()[ijet].Pt(),stopt.pfjets_csv()[ijet],(stopt.pfjets_mcflavorAlgo()[ijet]==5?5:0)); //only reshape for b jets

        //alternative reshaping for systematics
        //float csv_upBC=upBCShape->reshape(stopt.pfjets()[ijet].Eta(),stopt.pfjets()[ijet].Pt(),stopt.pfjets_csv()[ijet],stopt.pfjets_mcflavorAlgo()[ijet]);
        //float csv_downBC=downBCShape->reshape(stopt.pfjets()[ijet].Eta(),stopt.pfjets()[ijet].Pt(),stopt.pfjets_csv()[ijet],stopt.pfjets_mcflavorAlgo()[ijet]);
        //float csv_upL=upLShape->reshape(stopt.pfjets()[ijet].Eta(),stopt.pfjets()[ijet].Pt(),stopt.pfjets_csv()[ijet],stopt.pfjets_mcflavorAlgo()[ijet]);
        //float csv_downL=downLShape->reshape(stopt.pfjets()[ijet].Eta(),stopt.pfjets()[ijet].Pt(),stopt.pfjets_csv()[ijet],stopt.pfjets_mcflavorAlgo()[ijet]);
        //float csv_nominal4p=nominalShape4p->reshape(stopt.pfjets()[ijet].Eta(),stopt.pfjets()[ijet].Pt(),stopt.pfjets_csv()[ijet],stopt.pfjets_mcflavorAlgo()[ijet]);
        //float csv_upBC4p=upBCShape4p->reshape(stopt.pfjets()[ijet].Eta(),stopt.pfjets()[ijet].Pt(),stopt.pfjets_csv()[ijet],stopt.pfjets_mcflavorAlgo()[ijet]);
        //float csv_downBC4p=downBCShape4p->reshape(stopt.pfjets()[ijet].Eta(),stopt.pfjets()[ijet].Pt(),stopt.pfjets_csv()[ijet],stopt.pfjets_mcflavorAlgo()[ijet]);


        //fill csv discriminator shape histos for comparison with csvdiscr.root
        if ( stopt.pfjets_mcflavorAlgo()[ijet] == 5 ) plot1D("hb",       stopt.pfjets_csv()[ijet],       1., h_1d, 2000, -1, 1);
        else if ( stopt.pfjets_mcflavorAlgo()[ijet] == 4 ) plot1D("hc",       stopt.pfjets_csv()[ijet],       1., h_1d, 2000, -1, 1);
        else if ( stopt.pfjets_mcflavorAlgo()[ijet] == 3 || stopt.pfjets_mcflavorAlgo()[ijet] == 2 || stopt.pfjets_mcflavorAlgo()[ijet] == 1 ) plot1D("hl",       stopt.pfjets_csv()[ijet],       1., h_1d, 2000, -1, 1);

        //look at how the shapes are changed by the reshaping
        if ( stopt.pfjets_mcflavorAlgo()[ijet] == 5 ) plot1D("hb_afterreshaping",       csv_nominal,       1., h_1d, 2000, -1, 1);
        else if ( stopt.pfjets_mcflavorAlgo()[ijet] == 4 ) plot1D("hc_afterreshaping",       csv_nominal,       1., h_1d, 2000, -1, 1);
        else if ( stopt.pfjets_mcflavorAlgo()[ijet] == 3 || stopt.pfjets_mcflavorAlgo()[ijet] == 2 || stopt.pfjets_mcflavorAlgo()[ijet] == 1 ) plot1D("hl_afterreshaping",       csv_nominal,       1., h_1d, 2000, -1, 1);

        //calculate number of b-tags pre- and post-reshaping
        bool isbtagcsvm = false;
        if( stopt.pfjets_csv()[ijet]  > 0.679 ) {isbtagcsvm = true; nbtagscsvm_default++;}
        if( csv_nominal  > 0.679 ) nbtagscsvm_reshaped++;


        //now use random upgrading/downgrading method for comparison

        //set seed for random number generator based on event number and jet phi
        int randseed = stopt.event()+(int)(stopt.pfjets()[ijet].Phi()*1000.);
        random3_->SetSeed(randseed);
        float rand = random3_->Uniform(1.);  
        
        //start by reproducing result stored in babies

        //use old flavour definition using status 3 matching for consistency with method used in babies
        int pdgid = (stopt.pfjets_flav()[ijet]==5?5:0);
        
        //SFs used for result in babies
        float SFb_csvm  = getBtagSF( stopt.pfjets()[ijet].Pt(), stopt.pfjets()[ijet].Eta(), "CSVM");
        float SFl_csvm  = getMistagSF( stopt.pfjets()[ijet].Pt(), stopt.pfjets()[ijet].Eta(), "CSVM");
        float Effl_csvm = getMistags( stopt.pfjets()[ijet].Pt(), stopt.pfjets()[ijet].Eta(), "CSVM");
        bool iscorrbtagcsvm = getCorrBtag(isbtagcsvm, pdgid, SFb_csvm, SFl_csvm, Effl_csvm, rand);
        if (iscorrbtagcsvm) nbtagscsvm_random++;
        

        //now recalculate using SFs and flav used by reshaping method
        pdgid = (stopt.pfjets_mcflavorAlgo()[ijet]==5?5:0);
        SFb_csvm  = beff::CSVM_SFb( stopt.pfjets()[ijet].Pt() );
        SFl_csvm  = mistag_CSVM( fabs(stopt.pfjets()[ijet].Eta()),stopt.pfjets()[ijet].Pt(),0.);
        bool iscorrbtagcsvm_new = getCorrBtag(isbtagcsvm, pdgid, SFb_csvm, SFl_csvm, Effl_csvm, rand);
        if (iscorrbtagcsvm_new) nbtagscsvm_random_new++;


      }

      //first plot values stored in tree to cross-check results
      plot1D("h_nbtagscsvm",       stopt.nbtagscsvm()>5?5:stopt.nbtagscsvm(),       evtweight, h_1d, 6, 0, 6);
      plot1D("h_nbtagscsvmcorr",       stopt.nbtagscsvmcorr()>5?5:stopt.nbtagscsvmcorr(),       evtweight, h_1d, 6, 0, 6);

      //now calculated default and reshaped values
      plot1D("h_nbtagscsvm_default",       nbtagscsvm_default>5?5:nbtagscsvm_default,       evtweight, h_1d, 6, 0, 6);
      plot1D("h_nbtagscsvm_reshaped",       nbtagscsvm_reshaped>5?5:nbtagscsvm_reshaped,       evtweight, h_1d, 6, 0, 6);

      //now random method, using original and new SFs
      plot1D("h_nbtagscsvm_random",       nbtagscsvm_random>5?5:nbtagscsvm_random,       evtweight, h_1d, 6, 0, 6);
      plot1D("h_nbtagscsvm_random_new",       nbtagscsvm_random_new>5?5:nbtagscsvm_random_new,       evtweight, h_1d, 6, 0, 6);

      //----------------------------------------------------------------------------
      // ADD CODE BELOW THIS LINE
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
  
  outfile.Write();
  outfile.Close();

  already_seen.clear();

  gROOT->cd();

}

