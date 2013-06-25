//#include "../../CORE/jetSmearingTools.h"
#include "../../CORE/Thrust.h"
#include "../../CORE/EventShape.h"

#include "StopTreeLooper.h"

//#include "../Core/STOPT.h"
#include "../Core/stopUtils.h"
#include "../Plotting/PlotUtilities.h"
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
#include "TRandom3.h"
#include "TVector2.h"
#include "TLorentzVector.h"

#include <algorithm>
#include <utility>
#include <map>
#include <set>
#include <list>

using namespace Stop;

std::set<DorkyEventIdentifier> already_seen; 
std::set<DorkyEventIdentifier> events_lasercalib; 
std::set<DorkyEventIdentifier> events_hcallasercalib; 

struct SUSYGenParticle { // To be filled with status-3 genParticles
      int pdgId; // PDG identifier (with sign, please)
      int firstMother; // first mother, set to <0 if no mothers
      double energy; // energy [GeV]
      double pt; // pt [GeV]
      double eta; // eta
      double phi; // phi
};

double Reweight_Stop_to_TopChi0 (std::vector<SUSYGenParticle> genParticles, double referenceTopPolarization, double requestedTopPolarization, TString prefix);
double Reweight_T2bW (double thetaChi_eff, double thetaW_eff, std::vector<SUSYGenParticle> genParticles);

StopTreeLooper::StopTreeLooper()
{
    m_outfilename_ = "histos.root";
    // t1metphicorr = -9999.;
    // t1metphicorrphi = -9999.;
    // t1metphicorrmt = -9999.;
    // min_mtpeak = -9999.;
    // max_mtpeak = -9999.; 

    JET_PT = 30.;
    JET_ETA = 2.4;
    BJET_ETA = 2.4;

    N_JETS_TO_CONSIDER = 6;
    NJETS_CUT = 4;

    MET_CUT = 50.0;

    m_minibabylabel_ = "";

    __apply_mva = true;
}

StopTreeLooper::~StopTreeLooper()
{
}

void StopTreeLooper::setOutFileName(string filename)
{
    m_outfilename_ = filename;
    jets.clear();
    nonbjets.clear();
    jets_up.clear();
    jets_down.clear();
    bjets.clear();
    jets_btag.clear();
    jets_up_btag.clear();
    jets_down_btag.clear();
    jets_bup_btag.clear();
    jets_bdown_btag.clear();
    jets_sigma.clear();
    jets_up_sigma.clear();
    jets_down_sigma.clear();
    metphi  = -99999.;
}

void StopTreeLooper::doWHNtuple()
{
    JET_PT = 30.;
    JET_ETA = 4.7;
    BJET_ETA = 2.4;

    N_JETS_TO_CONSIDER = 6;
    NJETS_CUT = 2;

    MET_CUT = 50.0;

    m_minibabylabel_ = "_whmet";
}

void StopTreeLooper::disableMVA(){
  __apply_mva = false;
}

void StopTreeLooper::loop(TChain *chain, TString name)
{
    TRandom3 aleatorio;

    printf("[StopTreeLooper::loop] %s\n", name.Data());

    load_badlaserevents("../Core/badlaser_events.txt", events_lasercalib);
    load_badlaserevents("../Core/badhcallaser_events.txt", events_hcallasercalib);

    //---------------------------------
    // check for valid chain
    //---------------------------------

    TObjArray *listOfFiles = chain->GetListOfFiles();
    TIter fileIter(listOfFiles);
    if (listOfFiles->GetEntries() == 0) {
        cout << "[StopTreeLooper::loop] no files in chain" << endl;
        return;
    }

    //------------------------------------------------------------------------------------------------------
    // set csv discriminator reshaping
    //------------------------------------------------------------------------------------------------------
  
    BTagShapeInterface * nominalShape = new BTagShapeInterface("../../Tools/BTagReshaping/csvdiscr.root",  0.0 ,  0.0);
    //systematic variations for payloads
    BTagShapeInterface * upBCShape    = new BTagShapeInterface("../../Tools/BTagReshaping/csvdiscr.root",  1.0 ,  0.0);
    BTagShapeInterface * downBCShape  = new BTagShapeInterface("../../Tools/BTagReshaping/csvdiscr.root", -1.0 ,  0.0);
    BTagShapeInterface * upLShape     = new BTagShapeInterface("../../Tools/BTagReshaping/csvdiscr.root",  0.0 ,  1.0);
    BTagShapeInterface * downLShape   = new BTagShapeInterface("../../Tools/BTagReshaping/csvdiscr.root",  0.0 , -1.0);

    //---------------------------------
    // set up histograms
    //---------------------------------

    gROOT->cd();

    makeTree(name.Data(), chain);

    TH2F *h_nsig, *h_nsig25, *h_nsig75 ;
    TH2F *h_nsig_masslessLSP;

    if( name.Contains("T2") ){
        char* h_nsig_filename             = "";
        char* h_nsig_filename_masslessLSP = "";

        if( name.Contains("T2tt") && name.Contains("Filt") ){
	  h_nsig_filename             = "/nfs-7/userdata/stop/cms2V05-03-26_stoplooperV00-02-25/T2tt_mad/myMassDB_lepFilter_rescaled.root";
	  h_nsig_filename_masslessLSP = h_nsig_filename;
	  cout << "[StopTreeLooper::loop] opening mass TH2 file (massive LSP)  " << h_nsig_filename << endl;
	  cout << "[StopTreeLooper::loop] opening mass TH2 file (massless LSP) " << h_nsig_filename_masslessLSP << endl;
	}

        else if( name.Contains("T2tt") ){
	  h_nsig_filename             = "/nfs-7/userdata/stop/cms2V05-03-26_stoplooperV00-02-24/T2tt_mad/myMassDB_T2tt_massiveLSP.root";
	  h_nsig_filename_masslessLSP = "/nfs-7/userdata/stop/cms2V05-03-26_stoplooperV00-02-24/T2tt_mad/myMassDB_T2tt_masslessLSP.root";
	  cout << "[StopTreeLooper::loop] opening mass TH2 file (massive LSP)  " << h_nsig_filename << endl;
	  cout << "[StopTreeLooper::loop] opening mass TH2 file (massless LSP) " << h_nsig_filename_masslessLSP << endl;
	}

        else if( name.Contains("T2bw") ){

	  if( name.Contains("fine") ){
	    h_nsig_filename = "/nfs-7/userdata/stop/cms2V05-03-26_stoplooperV00-02-23/T2bw_pythia_fine/myMassDB_T2bw_fine.root";
	  }

	  else if( name.Contains("coarse") ){
	    h_nsig_filename = "/nfs-7/userdata/stop/cms2V05-03-26_stoplooperV00-02-25/T2bw_pythiaCoarse/myMassDB_T2bw_coarse.root";
	  }

	  else if( name.Contains("x050") ){
	    h_nsig_filename = "/nfs-7/userdata/stop/cms2V05-03-26_stoplooperV00-02-25/T2bw_mad/myMassDB_mStop_combined_x50.root";
	  }

	  else if( name.Contains("x025") ){
	    h_nsig_filename = "/nfs-7/userdata/stop/cms2V05-03-26_stoplooperV00-02-25/T2bw_mad/myMassDB_mStop_combined_x25.root";
	  }

	  else if( name.Contains("x075") ){
	    h_nsig_filename = "/nfs-7/userdata/stop/cms2V05-03-26_stoplooperV00-02-25/T2bw_mad/myMassDB_mStop_combined_x75.root";
	  }

	  else{
	    cout << __FILE__ << " " << __LINE__ << " ERROR! I don't recognize " << name << ", quitting" << endl;
	    exit(0);
	  }

	  cout << "[StopTreeLooper::loop] opening mass TH2 file  " << h_nsig_filename << endl;
	}

        TFile *f_nsig = TFile::Open(h_nsig_filename);

        assert(f_nsig);

        h_nsig = (TH2F*) f_nsig->Get("masses");

        if( name.Contains("T2bw") ){
            h_nsig25 = (TH2F*) f_nsig->Get("masses25");
            h_nsig75 = (TH2F*) f_nsig->Get("masses75");
        }

        if( name.Contains("T2tt") ){
	  TFile *f_nsig_masslessLSP = TFile::Open(h_nsig_filename_masslessLSP);
	  h_nsig_masslessLSP = (TH2F*) f_nsig_masslessLSP->Get("masses");
	}

    }

    const int NREG_T2tt = 6;
    const int NREG_T2bw = 5;
    TMVA::Reader* reader_T2tt[NREG_T2tt];
    TMVA::Reader* reader_T2bw[3][NREG_T2bw];

    TMVA::Reader* reader_T2tt_up[NREG_T2tt];
    TMVA::Reader* reader_T2bw_up[3][NREG_T2bw];

    TMVA::Reader* reader_T2tt_down[NREG_T2tt];
    TMVA::Reader* reader_T2bw_down[3][NREG_T2bw];

    TMVA::Reader* reader_T2tt_bup[NREG_T2tt];
    TMVA::Reader* reader_T2bw_bup[3][NREG_T2bw];

    TMVA::Reader* reader_T2tt_bdown[NREG_T2tt];
    TMVA::Reader* reader_T2bw_bdown[3][NREG_T2bw];

    if ( __apply_mva ) {
        for (int i=1; i < NREG_T2tt ; i++){
            reader_T2tt[i] = new TMVA::Reader( "!Color:!Silent" );
            reader_T2tt[i]->AddVariable("mini_met", &met_);
            reader_T2tt[i]->AddVariable("mini_mt2w", &mt2w_);
            if ( i != 5 )
                reader_T2tt[i]->AddVariable("mini_chi2", &chi2_);
            reader_T2tt[i]->AddVariable("mini_htssm/(mini_htosm+mini_htssm)", &htratiom_);
            reader_T2tt[i]->AddVariable("mini_dphimjmin", &dphimjmin_);
            if ( i == 5 )
                reader_T2tt[i]->AddVariable("mini_pt_b", &pt_b_);

	    // JES up
	    reader_T2tt_up[i] = new TMVA::Reader( "!Color:!Silent" );
	    reader_T2tt_up[i]->AddVariable("mini_met", &metup_);
	    reader_T2tt_up[i]->AddVariable("mini_mt2w", &mt2wup_);
	    if ( i != 5 )
	      reader_T2tt_up[i]->AddVariable("mini_chi2", &chi2up_);
	    reader_T2tt_up[i]->AddVariable("mini_htssm/(mini_htosm+mini_htssm)", &htratiomup_);
	    reader_T2tt_up[i]->AddVariable("mini_dphimjmin", &dphimjmin_);
	    if ( i == 5 )
	      reader_T2tt_up[i]->AddVariable("mini_pt_b", &pt_b_up_);

	    // JES down
            reader_T2tt_down[i] = new TMVA::Reader( "!Color:!Silent" );
            reader_T2tt_down[i]->AddVariable("mini_met", &metdown_);
            reader_T2tt_down[i]->AddVariable("mini_mt2w", &mt2wdown_);
            if ( i != 5 )
                reader_T2tt_down[i]->AddVariable("mini_chi2", &chi2down_);
            reader_T2tt_down[i]->AddVariable("mini_htssm/(mini_htosm+mini_htssm)", &htratiomdown_);
            reader_T2tt_down[i]->AddVariable("mini_dphimjmin", &dphimjmin_);
            if ( i == 5 )
                reader_T2tt_down[i]->AddVariable("mini_pt_b", &pt_b_down_);

	    // btagging up
	    reader_T2tt_bup[i] = new TMVA::Reader( "!Color:!Silent" );
	    reader_T2tt_bup[i]->AddVariable("mini_met", &met_);
	    reader_T2tt_bup[i]->AddVariable("mini_mt2w", &mt2wbup_);
	    if ( i != 5 )
	      reader_T2tt_bup[i]->AddVariable("mini_chi2", &chi2bup_);
	    reader_T2tt_bup[i]->AddVariable("mini_htssm/(mini_htosm+mini_htssm)", &htratiom_);
	    reader_T2tt_bup[i]->AddVariable("mini_dphimjmin", &dphimjmin_);
	    if ( i == 5 )
	      reader_T2tt_bup[i]->AddVariable("mini_pt_b", &pt_b_bup_);

	    // btagging down
            reader_T2tt_bdown[i] = new TMVA::Reader( "!Color:!Silent" );
            reader_T2tt_bdown[i]->AddVariable("mini_met", &met_);
            reader_T2tt_bdown[i]->AddVariable("mini_mt2w", &mt2wbdown_);
            if ( i != 5 )
                reader_T2tt_bdown[i]->AddVariable("mini_chi2", &chi2bdown_);
            reader_T2tt_bdown[i]->AddVariable("mini_htssm/(mini_htosm+mini_htssm)", &htratiom_);
            reader_T2tt_bdown[i]->AddVariable("mini_dphimjmin", &dphimjmin_);
            if ( i == 5 )
                reader_T2tt_bdown[i]->AddVariable("mini_pt_b", &pt_b_bdown_);

            TString dir, prefix;
            dir    = "/home/users/magania/stop/SingleLepton2012/MVA/weights/";
            prefix = "classification_T2tt_";
            prefix += i;
            prefix += "_BDT";

            TString weightfile = dir + prefix + TString(".weights.xml");
            reader_T2tt[i]      ->BookMVA( "BDT" , weightfile );
            reader_T2tt_up[i]   ->BookMVA( "BDT" , weightfile );
            reader_T2tt_down[i] ->BookMVA( "BDT" , weightfile );
            reader_T2tt_bup[i]  ->BookMVA( "BDT" , weightfile );
            reader_T2tt_bdown[i]->BookMVA( "BDT" , weightfile );
        }

        for (int j=0; j < 3; j++){
            float x;
            switch(j){
                case 0: x=0.25;break;
                case 1: x=0.50;break;
                case 2: x=0.75;break;
                default: x=0.;
            }
            for (int i=1; i < NREG_T2bw ; i++){
                if ( j==0 && i==1 ) continue;
                reader_T2bw[j][i] = new TMVA::Reader( "!Color:!Silent" );
                reader_T2bw[j][i]->AddVariable("mini_met", &met_);
                reader_T2bw[j][i]->AddVariable("mini_mt2w", &mt2w_);
                if ( j!=2 && i==4 )
                reader_T2bw[j][i]->AddVariable("mini_lep1pt", &lep1pt_);
                reader_T2bw[j][i]->AddVariable("mini_htssm/(mini_htosm+mini_htssm)", &htratiom_);
                reader_T2bw[j][i]->AddVariable("mini_dphimjmin", &dphimjmin_);
                reader_T2bw[j][i]->AddVariable("mini_pt_b", &pt_b_);
                reader_T2bw[j][i]->AddVariable("mini_dRleptB1", &dRleptB1_);

		// JES up
                reader_T2bw_up[j][i] = new TMVA::Reader( "!Color:!Silent" );
                reader_T2bw_up[j][i]->AddVariable("mini_met", &metup_);
                reader_T2bw_up[j][i]->AddVariable("mini_mt2w", &mt2wup_);
                if ( j!=2 && i==4 )
                reader_T2bw_up[j][i]->AddVariable("mini_lep1pt", &lep1pt_);
                reader_T2bw_up[j][i]->AddVariable("mini_htssm/(mini_htosm+mini_htssm)", &htratiomup_);
                reader_T2bw_up[j][i]->AddVariable("mini_dphimjmin", &dphimjmin_);
                reader_T2bw_up[j][i]->AddVariable("mini_pt_b", &pt_b_up_);
                reader_T2bw_up[j][i]->AddVariable("mini_dRleptB1", &dRleptB1_);

		// JES down
                reader_T2bw_down[j][i] = new TMVA::Reader( "!Color:!Silent" );
                reader_T2bw_down[j][i]->AddVariable("mini_met", &metdown_);
                reader_T2bw_down[j][i]->AddVariable("mini_mt2w", &mt2wdown_);
                if ( j!=2 && i==4 )
                reader_T2bw_down[j][i]->AddVariable("mini_lep1pt", &lep1pt_);
                reader_T2bw_down[j][i]->AddVariable("mini_htssm/(mini_htosm+mini_htssm)", &htratiomdown_);
                reader_T2bw_down[j][i]->AddVariable("mini_dphimjmin", &dphimjmin_);
                reader_T2bw_down[j][i]->AddVariable("mini_pt_b", &pt_b_down_);
                reader_T2bw_down[j][i]->AddVariable("mini_dRleptB1", &dRleptB1_);

		// btagging up
                reader_T2bw_bup[j][i] = new TMVA::Reader( "!Color:!Silent" );
                reader_T2bw_bup[j][i]->AddVariable("mini_met", &met_);
                reader_T2bw_bup[j][i]->AddVariable("mini_mt2w", &mt2wbup_);
                if ( j!=2 && i==4 )
                reader_T2bw_bup[j][i]->AddVariable("mini_lep1pt", &lep1pt_);
                reader_T2bw_bup[j][i]->AddVariable("mini_htssm/(mini_htosm+mini_htssm)", &htratiom_);
                reader_T2bw_bup[j][i]->AddVariable("mini_dphimjmin", &dphimjmin_);
                reader_T2bw_bup[j][i]->AddVariable("mini_pt_b", &pt_b_bup_);
                reader_T2bw_bup[j][i]->AddVariable("mini_dRleptB1", &dRleptB1_bup_);

		// btagging down
                reader_T2bw_bdown[j][i] = new TMVA::Reader( "!Color:!Silent" );
                reader_T2bw_bdown[j][i]->AddVariable("mini_met", &met_);
                reader_T2bw_bdown[j][i]->AddVariable("mini_mt2w", &mt2wbdown_);
                if ( j!=2 && i==4 )
                reader_T2bw_bdown[j][i]->AddVariable("mini_lep1pt", &lep1pt_);
                reader_T2bw_bdown[j][i]->AddVariable("mini_htssm/(mini_htosm+mini_htssm)", &htratiom_);
                reader_T2bw_bdown[j][i]->AddVariable("mini_dphimjmin", &dphimjmin_);
                reader_T2bw_bdown[j][i]->AddVariable("mini_pt_b", &pt_b_bdown_);
                reader_T2bw_bdown[j][i]->AddVariable("mini_dRleptB1", &dRleptB1_bdown_);

                TString dir    = "/home/users/magania/stop/SingleLepton2012/MVA/weights/";
                TString weightfile = Form("classification_T2bw_%d_%.2f_BDT.weights.xml",i,x);

                reader_T2bw[j][i]      ->BookMVA( "BDT" , dir+weightfile );
                reader_T2bw_up[j][i]   ->BookMVA( "BDT" , dir+weightfile );
                reader_T2bw_down[j][i] ->BookMVA( "BDT" , dir+weightfile );
                reader_T2bw_bup[j][i]  ->BookMVA( "BDT" , dir+weightfile );
                reader_T2bw_bdown[j][i]->BookMVA( "BDT" , dir+weightfile );
            }
        }

    }

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

    unsigned int nEventsPass=0;
    unsigned int nEventsChain=0;
    unsigned int nEvents = chain->GetEntries();
    nEventsChain = nEvents;
    ULong64_t nEventsTotal = 0;
    //  int i_permille_old = 0;

    bool isData = name.Contains("data") ? true : false;

    cout << "[StopTreeLooper::loop] running over chain with total entries " << nEvents << endl;

    cout << "[StopTreeLooper::loop] N jets requirement " << NJETS_CUT
        << " with pT>"<<JET_PT<<" GeV and |eta|<" << JET_ETA
        << endl;

    while (TChainElement *currentFile = (TChainElement*)fileIter.Next()) {

        //---------------------------------
        // load the stop baby tree
        //---------------------------------
        TFile *file = new TFile( currentFile->GetTitle() );
        TTree *tree = (TTree*)file->Get("t");
        stopt.Init(tree);

        if ( __add_babies )
            tree->CopyAddresses(outTree_);

        //---------------------------------
        // event loop
        //---------------------------------

        ULong64_t nEvents = tree->GetEntries();
        //nEvents = 1000;

        for(ULong64_t event = 0; event < nEvents; ++event) {
            stopt.GetEntry(event);
            tree->GetEntry(event);

            //---------------------------------
            // increment counters
            //---------------------------------

            ++nEventsTotal;
            //if (nEventsTotal%10000==0) {
            if (nEventsTotal%1000==0) {
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

                if (is_duplicate(id, already_seen) ) continue;
                if (is_badLaserEvent(id,events_lasercalib) ) continue;
                if (is_badLaserEvent(id,events_hcallasercalib) ) continue;

            }

            //------------------------------------------ 
            // init baby
            //------------------------------------------ 

            initBaby(); // set all branches to -1

            //------------------------------------------ 
            // event weight
            //------------------------------------------ 

	    float puweight = vtxweight_n( stopt.ntruepu(), h_pu_wgt, isData );

            weight_ = isData ? 1. : 
                ( stopt.weight() * 19.5 * puweight * stopt.mgcor() );

	    nsigevents_ = -1;

	    x_ = stopt.x();

            if( name.Contains("T2bw") ) {
	      if( stopt.mg() < -1 || stopt.ml() < -1 || stopt.x() < -1 ){
		cout << "ERROR! negative SUSY mass! Skip event!" << endl;
		cout << "run lumi event " << stopt.run() << " " << stopt.lumi() << " " << stopt.event() << endl;
		cout << "mstop " << stopt.mg() << endl;
		cout << "mLSP  " << stopt.ml() << endl;
		cout << "x     " << stopt.x()  << endl;
		continue;
	      }
	    }

            if( name.Contains("T2tt") ) {
	      if( stopt.mg() < -1 || stopt.ml() < -1 ){
		cout << "ERROR! negative SUSY mass! Skip event!" << endl;
		cout << "run lumi event " << stopt.run() << " " << stopt.lumi() << " " << stopt.event() << endl;
		cout << "mstop " << stopt.mg() << endl;
		cout << "mLSP  " << stopt.ml() << endl;
		continue;
	      }
	    }


            if( name.Contains("T2tt") ) {
                int bin = h_nsig->FindBin(stopt.mg(),stopt.ml());
                float nevents = h_nsig->GetBinContent(bin);
		if( stopt.ml() < 10.0 )
		  nevents = h_nsig_masslessLSP->GetBinContent(bin);

		// skip events with LSP mass = 0 because they're buggy
		// the fixed slice has LSP mass = 1 GeV
		if( stopt.ml() < 0.5 ) continue;

                //NOTE::need to add vtx. reweighting for the signal sample
                weight_  = stopt.xsecsusy() * 1000.0 / nevents * 19.5; 

                // cout << "mg " << tree->mg_ << " ml " << tree->ml_ 
                //      << " bin " << bin << " nevents " << nevents 
                //      << " xsec " << tree->xsecsusy_ 
                //      << " weight " << evtweight << endl;

		nsigevents_ = (int) nevents;
            }

            if( name.Contains("T2bw") ) {

	      //---------------------------------------------------------------------
	      // if x is buggy, calculate it by hand using the genparticle masses
	      //---------------------------------------------------------------------

	      float nevents = 0;
	      if ( x_ == 0.25 ){
		int bin = h_nsig25->FindBin(stopt.mg(),stopt.ml());
		nevents = h_nsig25->GetBinContent(bin);
	      } 

	      else if ( x_ == 0.50 ){
		int bin = h_nsig->FindBin(stopt.mg(),stopt.ml());
		nevents = h_nsig->GetBinContent(bin);
	      } 

	      else if ( x_ == 0.75 ){
		int bin = h_nsig75->FindBin(stopt.mg(),stopt.ml());
		nevents = h_nsig75->GetBinContent(bin);
	      }

	      if( nevents == 0 ){
		cout << __FILE__ << " " << __LINE__ << " ERROR! couldn't get nevents" << endl;
		cout << "event number " << stopt.event() << endl;
		cout << stopt.mg() << " " << stopt.ml() << " " << stopt.x() << endl;
	      }

	      assert ( nevents > 0 );

	      //NOTE::need to add vtx. reweighting for the signal sample
	      weight_ =  stopt.xsecsusy() * 1000.0 / nevents * 19.5; 

	      nsigevents_ = (int) nevents;
	      // cout << "mg " << stopt.mg() << " ml " << stopt.ml() 
	      //      << " x " << stopt.x() 
	      //      << " nevents " << nevents 
	      //      << " xsec " << stopt.xsecsusy() 
	      //      << " weight " << evtweight << endl;
            }

	    if( name.Contains("HHWWbb") ){
	      weight_ = stopt.weight() * 19.5;
	    }

	    //	    nvtxweight_ = stopt.nvtxweight();
	    nvtxweight_ = puweight;

            sltrigeff_   = isData ? 1. : 
	      getsltrigweight(stopt.id1(), stopt.lep1().Pt(), stopt.lep1().Eta());
            dltrigeff_ = isData ? 1. : 
                getdltrigweight(stopt.id1(), stopt.id2());

            //------------------------------------------ 
            // polarization reweighting
            //------------------------------------------ 

	    mini_weightleft_  = 1.0;
	    mini_weightright_ = 1.0;

	    mini_t2bwweight_lr_ = 1.0;
	    mini_t2bwweight_ls_ = 1.0;
	    mini_t2bwweight_ll_ = 1.0;

	    mini_t2bwweight_sr_ = 1.0;
	    mini_t2bwweight_ss_ = 1.0;
	    mini_t2bwweight_sl_ = 1.0;

	    mini_t2bwweight_rr_ = 1.0;
	    mini_t2bwweight_rs_ = 1.0;
	    mini_t2bwweight_rl_ = 1.0;

            if( name.Contains("T2") ) {

	      std::vector<SUSYGenParticle> genParticles;
	
	      for (unsigned int ig=0; ig < genps_pdgId().size(); ++ig) {
		
		SUSYGenParticle part;
		part.pdgId       = stopt.genps_pdgId().at(ig);
		part.energy      = stopt.genps_energy().at(ig);
		part.pt          = stopt.genps_pt().at(ig);
		part.eta         = stopt.genps_eta().at(ig);
		part.phi         = stopt.genps_phi().at(ig);
		part.firstMother = stopt.genps_firstMother().at(ig);
		genParticles.push_back(part);
	      }

	      if( name.Contains("T2tt") ) {	      
		mini_weightleft_  = Reweight_Stop_to_TopChi0 (genParticles, 0., -1, name);
		mini_weightright_ = Reweight_Stop_to_TopChi0 (genParticles, 0.,  1, name);
	      }

	      if( name.Contains("T2bw") ) {
		float pi = acos(-1.0);

		mini_t2bwweight_lr_ = Reweight_T2bW( 0 , 0      , genParticles );
		mini_t2bwweight_ls_ = Reweight_T2bW( 0 , pi/4.0 , genParticles );
		mini_t2bwweight_ll_ = Reweight_T2bW( 0 , pi/2.0 , genParticles );

		mini_t2bwweight_sr_ = Reweight_T2bW( pi/4.0 , 0      , genParticles );
		mini_t2bwweight_ss_ = Reweight_T2bW( pi/4.0 , pi/4.0 , genParticles );
		mini_t2bwweight_sl_ = Reweight_T2bW( pi/4.0 , pi/2.0 , genParticles );

		mini_t2bwweight_rr_ = Reweight_T2bW( pi/2.0 , 0      , genParticles );
		mini_t2bwweight_rs_ = Reweight_T2bW( pi/2.0 , pi/4.0 , genParticles );
		mini_t2bwweight_rl_ = Reweight_T2bW( pi/2.0 , pi/2.0 , genParticles );
	      }

	      // cout << endl << endl;
	      // cout << "weightleft  " << stopt.weightleft()  << " " << mini_weightleft_  << endl;
	      // cout << "weightright " << stopt.weightright() << " " << mini_weightright_ << endl;
	    }

            //------------------------------------------ 
            // selection criteria
            //------------------------------------------ 

            // >=1 lepton, rho cut, MET filters, veto 2 nearby leptons
            if ( !passEvtSelection(name) ) continue; 

            // MET > cut value (100 GeV for stop)
            met_       = stopt.t1metphicorr();
            metphi     = stopt.t1metphicorrphi();
	    
	    metup_     = stopt.t1metphicorrup();
	    metupphi   = stopt.t1metphicorrphiup();

	    metdown_   = stopt.t1metphicorrdn();
	    metdownphi = stopt.t1metphicorrphidn();

	    // warning for possibly bad MET up/down variables
	    // if( fabs( stopt.t1metphicorrdn() - met_ ) > 50.0 || fabs( stopt.t1metphicorrup() - met_ ) > 50.0 ){
	    //   cout << __FILE__ << " " << __LINE__ << " WARNING! large difference in MET JES up/down variables" << endl;
	    //   cout << "Nominal MET " << met_     << endl;
	    //   cout << "MET up      " << metup_   << endl;
	    //   cout << "MET down    " << metdown_ << endl;
	    // }

            if ( met_ < MET_CUT && metup_ < MET_CUT && metdown_ < MET_CUT ) continue; 

            // jet selection
            jets.clear();
	    nonbjets.clear();
            jets_up.clear();
            jets_down.clear();
            bjets.clear();
            jets_btag.clear();
            jets_up_btag.clear();
            jets_down_btag.clear();
            jets_bup_btag.clear();
            jets_bdown_btag.clear();
            jets_sigma.clear();
            jets_up_sigma.clear();
            jets_down_sigma.clear();
            njets_ = 0;
            njets_up_ = 0;
            njets_down_ = 0;
            nb_ = 0;
            nb_upBCShape_ = 0;
            nb_downBCShape_ = 0;
            nb_upLShape_ = 0;
            nb_downLShape_ = 0;
	    nnonbjets_ = 0;

            // kinematic variables
            htssl_ = 0.;
            htosl_ = 0.;
            htssm_ = 0.;
            htosm_ = 0.;

	    float htssmup   = 0.0;
	    float htssmdown = 0.0;
	    float htosmup   = 0.0;
	    float htosmdown = 0.0;

	    int nb_bup   = 0;
	    int nb_bdown = 0;

            for (unsigned int i =0; i<stopt.pfjets().size(); i++){

	      bool passTightPUid = passMVAJetId(stopt.pfjets().at(i).pt(), stopt.pfjets().at(i).eta(),stopt.pfjets_mva5xPUid().at(i),0);

	      float dPhiM = getdphi( metphi, stopt.pfjets().at(i).phi() );

	      // b-tagging information
	       
	      //float csv_nominal=(isData || name.Contains("T2")) ? stopt.pfjets_csv().at(i)
	      float csv_nominal=(isData) ? stopt.pfjets_csv().at(i)
		: nominalShape->reshape( stopt.pfjets().at(i).eta(),
					 stopt.pfjets().at(i).pt(),
					 stopt.pfjets_csv().at(i),
					 stopt.pfjets_mcflavorAlgo().at(i) ); 

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

	      // jets with up and down variations
	      float unc     = stopt.pfjets_uncertainty().at(i);
	      float pt_up   = stopt.pfjets().at(i).pt() * ( 1 + unc );
	      float pt_down = stopt.pfjets().at(i).pt() * ( 1 - unc );
	      float eta     = fabs(stopt.pfjets().at(i).eta());

	      if( pt_up > JET_PT && eta < JET_ETA && passTightPUid ){
	      	njets_up_++;
	      	jets_up.push_back      ( (1+unc)*stopt.pfjets().at(i) );
                jets_up_sigma.push_back( stopt.pfjets_sigma().at(i)   );
		jets_up_btag.push_back ( csv_nominal                  );

		if ( dPhiM  < (TMath::Pi()/2) ) htssmup = htssmup + (1+unc)*stopt.pfjets().at(i).pt();
		else                            htosmup = htosmup + (1+unc)*stopt.pfjets().at(i).pt();

	      }

	      if( pt_down > JET_PT && eta < JET_ETA && passTightPUid ){
	      	njets_down_++;
	      	jets_down.push_back      ( (1-unc)*stopt.pfjets().at(i) );
                jets_down_sigma.push_back( stopt.pfjets_sigma().at(i)   );
		jets_down_btag.push_back ( csv_nominal                  );

		if ( dPhiM  < (TMath::Pi()/2) ) htssmdown = htssmdown + (1-unc)*stopt.pfjets().at(i).pt();
		else                            htosmdown = htosmdown + (1-unc)*stopt.pfjets().at(i).pt();
	      }
	      
	      if ( stopt.pfjets().at(i).pt()        < JET_PT  ) continue;
	      if ( fabs(stopt.pfjets().at(i).eta()) > JET_ETA ) continue;
	      // if ( (fabs(stopt.pfjets().at(i).eta()) <= 2.5) 
	      //        && (stopt.pfjets_beta2_0p5().at(i) < 0.2) ) continue;
	      // bool passMediumPUid = passMVAJetId(stopt.pfjets().at(i).pt(), stopt.pfjets().at(i).eta(),stopt.pfjets_mvaPUid().at(i),1);

	      if(!passTightPUid) continue;

	      // jet information
	      njets_++;
	      if (njets_==1) {
		pt_J1_   = stopt.pfjets().at(i).pt();
		dphimj1_ = getdphi(metphi, stopt.pfjets().at(i).phi() );
	      }
	      if (njets_==2) {
		pt_J2_ = stopt.pfjets().at(i).pt();
		dphimj2_    = getdphi(metphi, stopt.pfjets().at(i).phi() );
		dphimjmin_  = TMath::Min( dphimj1_ , dphimj2_ );
	      }
	      jets.push_back(stopt.pfjets().at(i));
	      jets_sigma.push_back(stopt.pfjets_sigma().at(i));

 

	      if ( (fabs(stopt.pfjets().at(i).eta()) <= BJET_ETA) && (csv_nominal > 0.679) ) {
		nb_++;
		bjets.push_back(stopt.pfjets().at(i));
		// leading b-jet variables
		if (nb_==1) {
		  pt_b_      = stopt.pfjets().at(i).pt();
		  pt_b_up_   = (1+unc)*stopt.pfjets().at(i).pt();
		  pt_b_down_ = (1-unc)*stopt.pfjets().at(i).pt();
		  dRleptB1_  = deltaR(stopt.pfjets().at(i).eta() , 
				      stopt.pfjets().at(i).phi() , 
				      stopt.lep1().eta(), stopt.lep1().phi());
		}
	      }

	      else{
		nnonbjets_++;
		nonbjets.push_back(stopt.pfjets().at(i));
	      }

	      // bjet pt and dR(b,lep) with btagging up
	      if ( (fabs(stopt.pfjets().at(i).eta()) <= BJET_ETA) && (csv_upBCShape > 0.679) ) {
		nb_bup++;
		if( nb_bup == 1 ){
		  pt_b_bup_      = stopt.pfjets().at(i).pt();
		  dRleptB1_bup_  = deltaR(stopt.pfjets().at(i).eta() , 
					  stopt.pfjets().at(i).phi() , 
					  stopt.lep1().eta(), stopt.lep1().phi());
		}
	      }

	      // bjet pt and dR(b,lep) with btagging down
	      if ( (fabs(stopt.pfjets().at(i).eta()) <= BJET_ETA) && (csv_downBCShape > 0.679) ) {
		nb_bdown++;
		if( nb_bdown == 1 ){
		  pt_b_bdown_      = stopt.pfjets().at(i).pt();
		  dRleptB1_bdown_  = deltaR(stopt.pfjets().at(i).eta() , 
					    stopt.pfjets().at(i).phi() , 
					    stopt.lep1().eta(), stopt.lep1().phi());
		}
	      }
		
	      // for the cr1 we use the leading jet
	      if (nb_==0) {
		pt_b_      = stopt.pfjets().at(0).pt();
		pt_b_up_   = (1+unc)*stopt.pfjets().at(0).pt();
		pt_b_down_ = (1-unc)*stopt.pfjets().at(0).pt();
		dRleptB1_  = deltaR(stopt.pfjets().at(0).eta() , 
				    stopt.pfjets().at(0).phi() , 
				    stopt.lep1().eta(), stopt.lep1().phi());
	      }

	      if( nb_bup == 0 ){
		pt_b_bup_      = stopt.pfjets().at(0).pt();
		dRleptB1_bup_  = deltaR(stopt.pfjets().at(0).eta() , 
					stopt.pfjets().at(0).phi() , 
					stopt.lep1().eta(), stopt.lep1().phi());
	      }

	      if( nb_bdown == 0 ){
		pt_b_bdown_      = stopt.pfjets().at(0).pt();
		dRleptB1_bdown_  = deltaR(stopt.pfjets().at(0).eta() , 
					  stopt.pfjets().at(0).phi() , 
					  stopt.lep1().eta(), stopt.lep1().phi());
	      }

	      jets_btag.push_back(csv_nominal);
	      jets_bup_btag.push_back(csv_upBCShape);
	      jets_bdown_btag.push_back(csv_downBCShape);

	      if ( (fabs(stopt.pfjets().at(i).eta()) <= BJET_ETA) && (csv_upBCShape > 0.679)   ) nb_upBCShape_++;
	      if ( (fabs(stopt.pfjets().at(i).eta()) <= BJET_ETA) && (csv_downBCShape > 0.679) ) nb_downBCShape_++;
	      if ( (fabs(stopt.pfjets().at(i).eta()) <= BJET_ETA) && (csv_upLShape > 0.679)    ) nb_upLShape_++;
	      if ( (fabs(stopt.pfjets().at(i).eta()) <= BJET_ETA) && (csv_downLShape > 0.679)  ) nb_downLShape_++;

	      float dPhiL = getdphi( stopt.lep1().Phi(), stopt.pfjets().at(i).phi() );


	      if ( dPhiL  < (TMath::Pi()/2) ) htssl_=htssl_+stopt.pfjets().at(i).pt();
	      else htosl_=htosl_+stopt.pfjets().at(i).pt();
	      if ( dPhiM  < (TMath::Pi()/2) ) htssm_=htssm_+stopt.pfjets().at(i).pt();
	      else htosm_=htosm_+stopt.pfjets().at(i).pt();

            }

	    // store events with >=4 jets with JES up
	    if ( njets_up_ < NJETS_CUT ) continue; 

            // event shapes: HT in hemispheres
            htratiol_   = htssl_ / (htosl_ + htssl_);
            htratiom_   = htssm_ / (htosm_ + htssm_);

            htratiomup_   = htssmup   / (htosmup   + htssmup  );
            htratiomdown_ = htssmdown / (htosmdown + htssmdown);

            //------------------------------------------ 
            // kinematic variables
            //------------------------------------------ 

            mt_     = (float)getMT(stopt.lep1().pt(), stopt.lep1().phi(), met_     , metphi     );
            mtup_   = (float)getMT(stopt.lep1().pt(), stopt.lep1().phi(), metup_   , metupphi   );
            mtdown_ = (float)getMT(stopt.lep1().pt(), stopt.lep1().phi(), metdown_ , metdownphi );
            mt2b_   = (float)calculateMT2w(jets, jets_btag, stopt.lep1(), met_, metphi, MT2b);
            mt2bl_  = (float)calculateMT2w(jets, jets_btag, stopt.lep1(), met_, metphi, MT2bl);
            mt2w_   = (float)calculateMT2w(jets, jets_btag, stopt.lep1(), met_, metphi, MT2w);
            chi2_   = (float)calculateChi2SNT(jets, jets_sigma, jets_btag);

	    // chi2 with JES up/down
            chi2up_    = (float)calculateChi2SNT(jets_up   , jets_up_sigma   , jets_up_btag  );
            chi2down_  = (float)calculateChi2SNT(jets_down , jets_down_sigma , jets_down_btag);

	    // chi2 with btagging up/down
            chi2bup_   = (float)calculateChi2SNT(jets, jets_sigma, jets_bup_btag  );
            chi2bdown_ = (float)calculateChi2SNT(jets, jets_sigma, jets_bdown_btag);

	    // MT2W with JES up/down
            mt2wup_   = (float)calculateMT2w(jets_up   , jets_up_btag   , stopt.lep1() , metup_   , metupphi   , MT2w);
            mt2wdown_ = (float)calculateMT2w(jets_down , jets_down_btag , stopt.lep1() , metdown_ , metdownphi , MT2w);

	    // MT2W with btagging up/down
            mt2wbup_  = (float)calculateMT2w(jets, jets_bup_btag   , stopt.lep1(), met_, metphi, MT2w);
            mt2wbdown_= (float)calculateMT2w(jets, jets_bdown_btag , stopt.lep1(), met_, metphi, MT2w);

            // for WH+MET ntuple
            TVector2 lep(stopt.lep1().px(), stopt.lep1().py());
            TVector2 met;
            met.SetMagPhi(met_, metphi);
            TVector2 w = lep+met; 
            wpt_ = w.Mod();
            if (nb_ >= 2) {
                LorentzVector bb = bjets.at(0) + bjets.at(1);
                bbmass_ = bb.M();
                bbpt_ = bb.pt();
                bbwdphi_ = fabs(TVector2::Phi_mpi_pi(bb.phi() - w.Phi()));
            }

	    // for HH->WWbb analysis
	    if( nnonbjets_ >= 2 ){
	      LorentzVector dijet = nonbjets.at(0) + nonbjets.at(1);

	      jjmaxpt_mass_ = dijet.M();
	      jjmaxpt_pt_   = dijet.pt();

	      float mindm = 99999;

	      for( unsigned int i = 0 ; i < nonbjets.size() ; i++ ){
		for( unsigned int j = i+1 ; j < nonbjets.size() ; j++ ){

		  LorentzVector this_dijet = nonbjets.at(i) + nonbjets.at(j);
		  
		  if( fabs( this_dijet.M() - 81.0 ) < mindm ){
		    mindm     = fabs( this_dijet.M() - 81.0 );
		    jjw_mass_ = this_dijet.M();
		    jjw_pt_   = this_dijet.pt();
		  }

		}
	      }

	    }

            //------------------------------------------ 
            // datasets bit
            //------------------------------------------ 

            bool dataset_1l=false;

            if((isData) && name.Contains("muo") 
                    && (abs(stopt.id1()) == 13 ))  dataset_1l=true;
            if((isData) && name.Contains("ele") 
                    && (abs(stopt.id1()) == 11 ))  dataset_1l=true;

            if(!isData) dataset_1l=true;

            bool dataset_CR4=false;

            if((isData) && name.Contains("dimu") 
                    && (abs(stopt.id1()) == 13 ) && (abs(stopt.id2())==13) 
                    && (fabs( stopt.dilmass() - 91.) > 15.)) dataset_CR4=true;
            if((isData) && name.Contains("diel") 
                    && (abs(stopt.id1()) == 11 ) && (abs(stopt.id2())==11) 
                    && (fabs( stopt.dilmass() - 91.) > 15.)) dataset_CR4=true;
            if((isData) && name.Contains("mueg") 
                    && abs(stopt.id1()) != abs(stopt.id2())) dataset_CR4=true;

            if(!isData) dataset_CR4=true;

            //------------------------------------------ 
            // variables to add to baby
            //------------------------------------------ 

            // isolated track veto selection
            passisotrk_   = passIsoTrkVeto_v4() ? 1 : 0; 
            // tau veto selection 
	    passtauveto_  = (passTauVeto()) ? 1 : 0;

            // which selections are passed
            // pass signal region preselection
            sig_ = ( dataset_1l 
		     && passOneLeptonSelection(isData) 
		     && passisotrk_ == 1
		     && passtauveto_ == 1
		     && nb_>=1 ) ? 1 : 0; 
            // pass CR1 (b-veto) control region preselection
            cr1_ = ( dataset_1l 
		     && passOneLeptonSelection(isData) 
		     && passisotrk_ == 1
		     && passtauveto_ == 1
		     && nb_==0 ) ? 1 : 0; 
            // pass CR4 (dilepton) control region preselection
            cr4_ = ( dataset_CR4 
		     && passDileptonSelection(isData) 
		     && nb_>=1) ? 1 : 0; 
            // pass CR1 (lepton+isotrack) control region preselection
            cr5_ = ( dataset_1l 
		     && ( passLepPlusIsoTrkSelection(isData) || passLepPlusTauSelection(isData) )
		     && nb_>=1) ? 1 : 0; 
            // -- WH+MET analysis regions
            // signal region
            whsig_      = ( dataset_1l 
                    && passOneLeptonSelection(isData) 
                    && passisotrk_ == 1
                    && nb_ == 2 
                    && njets_ == 2
                    && (bbmass_ >= 100. && bbmass_ <= 140.) ) ? 1 : 0; 

            // pass CR1: high bbmass
            whcr1_      = ( dataset_1l 
                    && passOneLeptonSelection(isData) 
                    && passisotrk_ == 1
                    && nb_ == 2 
                    && njets_ == 2
                    && bbmass_ >= 150. ) ? 1 : 0; 

            // Cut & Count analysis Selections
	    // tag_T2tt_LM -- dphi and chi2 selection
	    t2ttLM_ = (dphimjmin_>0.8 && chi2_<5) ? true : false;
	    // tag_T2tt_HM -- add mt2w requirement
	    t2ttHM_ = (t2ttLM_ && mt2w_>200) ? true : false;
	    //until the functions are updated
            // t2ttLM_     = pass_T2tt_LM(isData,name);
            // t2ttHM_     = pass_T2tt_HM(isData,name);

	    // susy vars
	    mstop_       = stopt.mg();                   // stop mass
	    mlsp_        = stopt.ml();                   // LSP mass
	    // I moved this line higher up
	    //x_           = stopt.x();                    // chargino mass parameter x
	    xsecsusy_    = stopt.xsecsusy();

            // number of analysis selected leptons
            nlep_       = stopt.ngoodlep();              
            lep1pt_     = stopt.lep1().pt();             // 1st lepton pt
            lep1eta_    = fabs( stopt.lep1().eta() );    // 1st lepton eta

            if( nlep_ > 1 ){
                lep2pt_    = stopt.lep1().pt();            // 2nd lepton pt
                lep2eta_   = stopt.lep1().eta();           // 2nd lepton eta
                dilmass_   = stopt.dilmass();              // dilepton mass
            }

            rand_       = aleatorio.Uniform(0.0,1.0);

	    float isrboost = 0.0;

	    if     ( name.Contains("T2") ) isrboost = (stopt.stop_t() + stopt.stop_tbar()).pt();
	    else if( name.Contains("ttsl") || name.Contains("ttdl") || name.Contains("ttall") || name.Contains("tt_") ) isrboost = (stopt.ttbar()).pt();

	    isrweight_ = 1.0;
	    if( isrboost > 120.0 && isrboost < 150.0 ) isrweight_ = 0.95;
	    if( isrboost > 150.0 && isrboost < 250.0 ) isrweight_ = 0.90;
	    if( isrboost > 250.0                     ) isrweight_ = 0.80;

            for ( int i=0; i < NREG_T2tt; i++){

	      float bdtval     = 0;
	      float bdtvalup   = 0;
	      float bdtvaldown = 0;
	      float bdtvalbup   = 0;
	      float bdtvalbdown = 0;

	      if ( __apply_mva && i>0 ){
		bdtval      = reader_T2tt[i]      ->EvaluateMVA( "BDT" );
		bdtvalup    = reader_T2tt_up[i]   ->EvaluateMVA( "BDT" );
		bdtvaldown  = reader_T2tt_down[i] ->EvaluateMVA( "BDT" );
		bdtvalbup   = reader_T2tt_bup[i]  ->EvaluateMVA( "BDT" );
		bdtvalbdown = reader_T2tt_bdown[i]->EvaluateMVA( "BDT" );
	      }

	      bdt_     .push_back(bdtval);
	      bdtup_   .push_back(bdtvalup);
	      bdtdown_ .push_back(bdtvaldown);
	      bdtbup_  .push_back(bdtvalbup);
	      bdtbdown_.push_back(bdtvalbdown);
            }

            for ( int j=0; j < 3; j++){
	      for ( int i=0; i < NREG_T2bw; i++){

		float bdtval      = 0;
		float bdtvalup    = 0;
		float bdtvaldown  = 0;
		float bdtvalbup   = 0;
		float bdtvalbdown = 0;

		if ( __apply_mva && i>0 && !(j==0&&i==1) ){
		  bdtval      = reader_T2bw[j][i]      ->EvaluateMVA( "BDT" );
		  bdtvalup    = reader_T2bw_up[j][i]   ->EvaluateMVA( "BDT" );
		  bdtvaldown  = reader_T2bw_down[j][i] ->EvaluateMVA( "BDT" );
		  bdtvalbup   = reader_T2bw_bup[j][i]  ->EvaluateMVA( "BDT" );
		  bdtvalbdown = reader_T2bw_bdown[j][i]->EvaluateMVA( "BDT" );
		}

		bdt_     .push_back(bdtval);
		bdtup_   .push_back(bdtvalup);
		bdtdown_ .push_back(bdtvaldown);
		bdtbup_  .push_back(bdtvalbup);
		bdtbdown_.push_back(bdtvalbdown);

	      }
	    }

            // fill me up
            nEventsPass++;
            outTree_->Fill();


            } // end event loop

            delete tree;
        } // end file loop

        //-------------------------
        // finish and clean up
        //-------------------------

        cout << "[StopTreeLooper::loop] saving mini-baby with total entries " << nEventsPass << endl;

        outFile_->cd();
        outTree_->Write();
        outFile_->Close();
        delete outFile_;

        already_seen.clear();

        gROOT->cd();
    }

    //--------------------------------------------
    // create the tree and set branch addresses
    //--------------------------------------------

    void StopTreeLooper::makeTree(const char *prefix, TChain* chain){
        TDirectory *rootdir = gDirectory->GetDirectory("Rint:");
        rootdir->cd();

//        outFile_   = new TFile(Form("output/%s%s.root", prefix, m_minibabylabel_.c_str()), "RECREATE");
        outFile_   = new TFile(Form("output_V00-03-13/%s%s.root", prefix, m_minibabylabel_.c_str()), "RECREATE");
	//        outFile_   = new TFile(Form("/nfs-7/userdata/stop/output_V00-02-21_2012_4jskim/Minibabies/%s%s.root", prefix, m_minibabylabel_.c_str()), "RECREATE");
        outFile_->cd();

        if ( __add_babies )
            outTree_ = chain->CloneTree(0);
        else
            outTree_ = new TTree("t","Tree"); 

        if ( __mini_branches) {
            outTree_->Branch("mini_mt"        , &mt_        ,  "mini_mt/F"	 );
            outTree_->Branch("mini_mtup"      , &mtup_      ,  "mini_mtup/F"	 );
            outTree_->Branch("mini_mtdown"    , &mtdown_    ,  "mini_mtdown/F"	 );

            outTree_->Branch("mini_isrweight" , &isrweight_ ,  "mini_isrweight/F");

            outTree_->Branch("mini_weightleft" , &mini_weightleft_  ,  "mini_weightleft/F" );
            outTree_->Branch("mini_weightright", &mini_weightright_ ,  "mini_weightright/F");

            outTree_->Branch("mini_t2bwweight_lr", &mini_t2bwweight_lr_ ,  "mini_t2bwweight_lr/F");
            outTree_->Branch("mini_t2bwweight_ls", &mini_t2bwweight_ls_ ,  "mini_t2bwweight_ls/F");
            outTree_->Branch("mini_t2bwweight_ll", &mini_t2bwweight_ll_ ,  "mini_t2bwweight_ll/F");

            outTree_->Branch("mini_t2bwweight_sr", &mini_t2bwweight_sr_ ,  "mini_t2bwweight_sr/F");
            outTree_->Branch("mini_t2bwweight_ss", &mini_t2bwweight_ss_ ,  "mini_t2bwweight_ss/F");
            outTree_->Branch("mini_t2bwweight_sl", &mini_t2bwweight_sl_ ,  "mini_t2bwweight_sl/F");

            outTree_->Branch("mini_t2bwweight_rr", &mini_t2bwweight_rr_ ,  "mini_t2bwweight_rr/F");
            outTree_->Branch("mini_t2bwweight_rs", &mini_t2bwweight_rs_ ,  "mini_t2bwweight_rs/F");
            outTree_->Branch("mini_t2bwweight_rl", &mini_t2bwweight_rl_ ,  "mini_t2bwweight_rl/F");

            outTree_->Branch("mini_met"       , &met_       ,  "mini_met/F"	 );
            outTree_->Branch("mini_metup"     , &metup_     ,  "mini_metup/F"	 );
            outTree_->Branch("mini_metdown"   , &metdown_   ,  "mini_metdown/F"	 );

            outTree_->Branch("mini_sig"       , &sig_       ,  "mini_sig/I"	 );
            outTree_->Branch("mini_cr1"       , &cr1_       ,  "mini_cr1/I"	 );
            outTree_->Branch("mini_cr4"       , &cr4_       ,  "mini_cr4/I"	 );
            outTree_->Branch("mini_cr5"       , &cr5_       ,  "mini_cr5/I"	 );

            outTree_->Branch("mini_whsig"     , &whsig_     ,  "mini_whsig/I"	 );
            outTree_->Branch("mini_whcr1"     , &whcr1_     ,  "mini_whcr1/I"	 );

            outTree_->Branch("mini_chi2"      , &chi2_      ,  "mini_chi2/F"      );
            outTree_->Branch("mini_chi2up"    , &chi2up_    ,  "mini_chi2up/F"    );
            outTree_->Branch("mini_chi2down"  , &chi2down_  ,  "mini_chi2down/F"  );
            outTree_->Branch("mini_chi2bup"   , &chi2bup_   ,  "mini_chi2bup/F"   );
            outTree_->Branch("mini_chi2bdown" , &chi2bdown_ ,  "mini_chi2bdown/F" );
            outTree_->Branch("mini_mt2b"      , &mt2b_      ,  "mini_mt2b/F"      );
            outTree_->Branch("mini_mt2bl"     , &mt2bl_     ,  "mini_mt2bl/F"     );
            outTree_->Branch("mini_mt2w"      , &mt2w_      ,  "mini_mt2w/F"      );
            outTree_->Branch("mini_mt2wup"    , &mt2wup_    ,  "mini_mt2wup/F"    );
            outTree_->Branch("mini_mt2wdown"  , &mt2wdown_  ,  "mini_mt2wdown/F"  );
            outTree_->Branch("mini_mt2wbup"   , &mt2wbup_   ,  "mini_mt2wbup/F"   );
            outTree_->Branch("mini_mt2wbdown" , &mt2wbdown_ ,  "mini_mt2wbdown/F" );

            outTree_->Branch("mini_weight"    , &weight_    ,  "mini_weight/F"	 );
            outTree_->Branch("mini_nvtxweight", &nvtxweight_,  "mini_nvtxweight/F"	 );
            outTree_->Branch("mini_sltrigeff" , &sltrigeff_ ,  "mini_sltrigeff/F" );
            outTree_->Branch("mini_dltrigeff" , &dltrigeff_ ,  "mini_dltrigeff/F" );
            outTree_->Branch("mini_nsigevents", &nsigevents_,  "mini_nsigevents/I" );

            outTree_->Branch("mini_nb"        , &nb_        ,   "mini_nb/I"	         );
            outTree_->Branch("mini_nbupBC"    , &nb_upBCShape_  ,  "mini_nbupBC/I"	);
            outTree_->Branch("mini_nbdownBC"  , &nb_downBCShape_,  "mini_nbdownBC/I"	);
            outTree_->Branch("mini_nbupL"     , &nb_upLShape_   ,  "mini_nbupL/I"	);
            outTree_->Branch("mini_nbdownL"   , &nb_downLShape_ ,  "mini_nbdownL/I"	);
            outTree_->Branch("mini_njets"     , &njets_      ,  "mini_njets/I"  	 );
            outTree_->Branch("mini_njetsup"   , &njets_up_   ,  "mini_njetsup/I"	 );
            outTree_->Branch("mini_njetsdown" , &njets_down_ ,  "mini_njetsdown/I"	 );

            outTree_->Branch("mini_passisotrk" , &passisotrk_,  "mini_passisotrk/I");
	    outTree_->Branch("mini_passtauveto", &passtauveto_,  "mini_passtauveto/I");

            outTree_->Branch("mini_nlep"      , &nlep_      ,  "mini_nlep/I"	 );
            outTree_->Branch("mini_lep1pt"    , &lep1pt_    ,  "mini_lep1pt/F"	 );
            outTree_->Branch("mini_lep1eta"   , &lep1eta_   ,  "mini_lep1eta/F"	 );
            outTree_->Branch("mini_lep2pt"    , &lep2pt_    ,  "mini_lep2pt/F"	 );
            outTree_->Branch("mini_lep2eta"   , &lep2eta_   ,  "mini_lep2eta/F"	 );
            outTree_->Branch("mini_dilmass"   , &dilmass_   ,  "mini_dilmass/F"	 );
	    outTree_->Branch("mini_mstop"     , &mstop_     ,  "mini_mstop/F"	 );
	    outTree_->Branch("mini_mlsp"      , &mlsp_      ,  "mini_mlsp/F"	 );
	    outTree_->Branch("mini_x"         , &x_         ,  "mini_x/F"	 );
	    outTree_->Branch("mini_xsecsusy"  , &xsecsusy_  ,  "mini_xsecsusy/F" );
            outTree_->Branch("mini_dRleptB1"  , &dRleptB1_  ,  "mini_dRleptB1/F" );
            outTree_->Branch("mini_dRleptB1_bup"   , &dRleptB1_bup_    ,  "mini_dRleptB1_bup/F"	  );
            outTree_->Branch("mini_dRleptB1_bdown" , &dRleptB1_bdown_  ,  "mini_dRleptB1_bdown/F" );

            outTree_->Branch("mini_htssl"     , &htssl_     ,  "mini_htssl/F"     );
            outTree_->Branch("mini_htosl"     , &htosl_     ,  "mini_htosl/F"     );
            outTree_->Branch("mini_htratiol"  , &htratiol_  ,  "mini_htraiol/F"   );
            outTree_->Branch("mini_htssm"     , &htssm_     ,  "mini_htssm/F"     );
            outTree_->Branch("mini_htosm"     , &htosm_     ,  "mini_htosm/F"     );
            outTree_->Branch("mini_htratiom"  , &htratiom_  ,  "mini_htraiom/F"   );
            outTree_->Branch("mini_htratiomup"    , &htratiomup_    ,  "mini_htraiomup/F"     );
            outTree_->Branch("mini_htratiomdown"  , &htratiomdown_  ,  "mini_htraiomdown/F"   );
            outTree_->Branch("mini_dphimj1"   , &dphimj1_   ,  "mini_dphimj1/F"   );
            outTree_->Branch("mini_dphimj2"   , &dphimj2_   ,  "mini_dphimj2/F"   );
            outTree_->Branch("mini_dphimjmin" , &dphimjmin_ ,  "mini_dphimjmin/F" );

            outTree_->Branch("mini_pt_b"      , &pt_b_      ,  "mini_pt_b/F"      );
            outTree_->Branch("mini_pt_b_up"   , &pt_b_up_   ,  "mini_pt_b_up/F"   );
            outTree_->Branch("mini_pt_b_down" , &pt_b_down_ ,  "mini_pt_b_down/F" );
            outTree_->Branch("mini_pt_b_bup"  , &pt_b_bup_  ,  "mini_pt_b_bup/F"  );
            outTree_->Branch("mini_pt_b_bdown", &pt_b_bdown_,  "mini_pt_b_bdown/F");
            outTree_->Branch("mini_pt_J1"     , &pt_J1_     ,  "mini_pt_J1/F"     );
            outTree_->Branch("mini_pt_J2"     , &pt_J2_     ,  "mini_pt_J2/F"     );

            outTree_->Branch("mini_bbmass"    , &bbmass_    ,  "mini_bbmass/F"	 );
            outTree_->Branch("mini_bbpt"      , &bbpt_      ,  "mini_bbpt/F"	 );
            outTree_->Branch("mini_wpt"       , &wpt_       ,  "mini_wpt/F"	 );
            outTree_->Branch("mini_bbwdphi"   , &bbwdphi_   ,  "mini_bbwdphi/F"	 );

            outTree_->Branch("mini_rand"      , &rand_      ,  "mini_rand/F"      );

            outTree_->Branch("mini_t2ttLM"    , &t2ttLM_    ,  "mini_t2ttLM/F"    );
            outTree_->Branch("mini_t2ttHM"    , &t2ttHM_    ,  "mini_t2ttHM/F"    );

            outTree_->Branch("mini_nnonbjets"    , &nnonbjets_    ,  "mini_t2ttHM/I"        );
            outTree_->Branch("mini_jjmaxpt_mass" , &jjmaxpt_mass_ ,  "mini_jjmaxpt_mass/F"  );
            outTree_->Branch("mini_jjmaxpt_pt"   , &jjmaxpt_pt_   ,  "mini_jjmaxpt_pt/F"    );
            outTree_->Branch("mini_jjw_mass"     , &jjw_mass_     ,  "mini_jjw_mass/F"      );
            outTree_->Branch("mini_jjw_pt"       , &jjw_pt_       ,  "mini_jjw_pt/F"        );
        }

        if (__apply_mva){
	  outTree_->Branch("mini_bdt"       , "std::vector<float>" , &bdt_      );
	  outTree_->Branch("mini_bdtup"     , "std::vector<float>" , &bdtup_    );
	  outTree_->Branch("mini_bdtdown"   , "std::vector<float>" , &bdtdown_  );
	  outTree_->Branch("mini_bdtbup"    , "std::vector<float>" , &bdtbup_   );
	  outTree_->Branch("mini_bdtbdown"  , "std::vector<float>" , &bdtbdown_ );
	}

    }

    //--------------------------------------------
    // set all branches to -1
    //--------------------------------------------

    void StopTreeLooper::initBaby(){

        // which selections are passed
        sig_        = -1;
        cr1_        = -1;
        cr4_        = -1;
        cr5_        = -1; 

        whsig_      = -1;
        whcr1_      = -1;

        // kinematic variables
        mt_         = -1.0;
        mtup_       = -1.0;
        mtdown_     = -1.0;
        met_        = -1.0;
        metup_      = -1.0;
        metdown_    = -1.0;
        chi2_       = -1.0;
        chi2up_     = -1.0;
        chi2down_   = -1.0;
        chi2bup_    = -1.0;
        chi2bdown_  = -1.0;
        mt2b_       = -1.0;
        mt2bl_      = -1.0;
        mt2w_       = -1.0;
        mt2wup_     = -1.0;
        mt2wdown_   = -1.0;
	mt2wbup_    = -1.0;
        mt2wbdown_  = -1.0;

        // event shapes
        htssl_      = -1.0;
        htosl_      = -1.0;
        htratiol_   = -1.0;
        htssm_      = -1.0;
        htosm_      = -1.0;
        htratiom_   = -1.0;
        htratiomup_     = -1.0;
        htratiomdown_   = -1.0;

        // weights
        weight_     = -1.0;
        nvtxweight_ = -1.0;
        sltrigeff_  = -1.0;
        dltrigeff_  = -1.0;

        // jet counting
        nb_         = -1;
        nb_upBCShape_   = -1;
        nb_downBCShape_ = -1;
        nb_upLShape_    = -1;
        nb_downLShape_  = -1;
        njets_      = -1;
        njets_up_   = -1;
        njets_down_ = -1;

        // lepton variables
        passisotrk_ = -1;
        passtauveto_ = -1;

        nlep_       = -1;
        dRleptB1_   = -1;
        dRleptB1_bup_   = -1;
        dRleptB1_bdown_ = -1;
        lep1pt_     = -1.0;
        lep1eta_    = -1.0;
        lep2pt_     = -1.0;
        lep2eta_    = -1.0;
        dilmass_    = -1.0;

	// susy variables
	mstop_      = -1.0;
	mlsp_       = -1.0;
	x_          = -1.0;
	xsecsusy_   = -1.0;

        // jet kinematics
        pt_b_       = -1.0;
        pt_b_up_    = -1.0;
        pt_b_down_  = -1.0;
        pt_b_bup_   = -1.0;
        pt_b_bdown_ = -1.0;
        pt_J1_      = -1.0;
        pt_J2_      = -1.0;

        // wh event kinematics
        bbmass_      = -1.0;
        bbpt_        = -1.0;
        wpt_         = -1.0;
        bbwdphi_     = -1.0;


        rand_  = -1.0;

        t2ttLM_ =-1;
        t2ttLM_ =-1;

        bdt_.clear();
	bdtup_.clear();
	bdtdown_.clear();
	bdtbup_.clear();
	bdtbdown_.clear();
    }



// The following reweighting only makes sense for on-shell stop, top and chi0
// In the off-shell case top and anti-top may get very different polarizations
double Reweight_Stop_to_TopChi0 (std::vector<SUSYGenParticle> genParticles, double referenceTopPolarization, double requestedTopPolarization, TString prefix) {

  if( !prefix.Contains("T2tt") ) return 1.0;

  double weight = 1.;
  int nFoundStops = 0;

  unsigned int ngen = genParticles.size();

  for (unsigned int ig=0; ig<ngen; ++ig) {
    const SUSYGenParticle& gen = genParticles[ig];
    if (gen.firstMother<0) continue;
    if (abs(gen.pdgId)>20) continue; // expect quarks or leptons from W decay

    // Navigate upwards in the stop->top->W->fermion decay chain
    const SUSYGenParticle& genW = genParticles[gen.firstMother];
    if (genW.firstMother<0) continue;
    if (abs(genW.pdgId)!=24) continue;
    const SUSYGenParticle& genTop = genParticles[genW.firstMother];
    if (abs(genTop.pdgId)!=6) continue;

    // We only care about the down-type fermion
    if (genTop.pdgId*gen.pdgId>0) continue;

    // We also need a stop
    if (genTop.firstMother<0) continue;
    const SUSYGenParticle& genStop = genParticles[genTop.firstMother];
    if (abs(genStop.pdgId)!=1000006) continue;

    // Move top and fermion to the stop center-of-mass frame
    TLorentzVector stop4;
    stop4.SetPtEtaPhiE(genStop.pt, genStop.eta, genStop.phi, genStop.energy);
    TVector3 betaV(-stop4.Px()/stop4.Energy(),-stop4.Py()/stop4.Energy(),-stop4.Pz()/stop4.Energy());

    TLorentzVector top4;
    top4.SetPtEtaPhiE(genTop.pt, genTop.eta, genTop.phi, genTop.energy);
    top4.Boost(betaV);

    TLorentzVector ferm4;
    ferm4.SetPtEtaPhiE(gen.pt, gen.eta, gen.phi, gen.energy);
    ferm4.Boost(betaV);

    // Do not reweight if by any reason top/fermion directions are undefined
    // This should be pathological if things are fine
    if (top4.P()<=0 || ferm4.P()<=0) {
      printf("Warning: particles at rest, no weight applied: ptop: %.3e, pf: %.3e\n", top4.P(), ferm4.P());
      continue; 
    }

    double costh = (top4.Px()*ferm4.Px()+top4.Py()*ferm4.Py()+top4.Pz()*ferm4.Pz())/top4.P()/ferm4.P();
      
    double weight_L = (top4.Energy()+top4.P())*(1-costh);
    double weight_R = (top4.Energy()-top4.P())*(1+costh);
    weight *= ((1+requestedTopPolarization)*weight_R+(1-requestedTopPolarization)*weight_L)/((1+referenceTopPolarization)*weight_R+(1-referenceTopPolarization)*weight_L);

    nFoundStops++;
  }

  if( nFoundStops!=2 ) cout << __FILE__ << " " << __LINE__ << " WARNING: found " << nFoundStops << " stops, should be 2." << endl;

  return weight;

};

double Norm(double M1, double M2, double MV, double CL, double CR) {
      double lambda = pow(M1,4) + pow(M2,4) + pow(MV,4) - 2*pow(M1*M2,2) - 2*pow(M1*MV,2) - 2*pow(M2*MV,2);
      double norm = (CL*CL+CR*CR)*(lambda + 3*MV*MV*(M1*M1+M2*M2-MV*MV)) - 12*CL*CR*M1*M2*MV*MV;
      norm /= 3.;
      return norm;
}

void Boost_To_Stop_Rest_Frame(TLorentzVector& stop4, TLorentzVector& chargino4, TLorentzVector& b4, TLorentzVector& neutralino4, TLorentzVector& W4, TLorentzVector& up4, TLorentzVector& down4, TLorentzVector& s4){
      TVector3 betaV(-stop4.Px()/stop4.Energy(),-stop4.Py()/stop4.Energy(),-stop4.Pz()/stop4.Energy());
      stop4.Boost(betaV);
      chargino4.Boost(betaV);
      b4.Boost(betaV);
      neutralino4.Boost(betaV);
      W4.Boost(betaV);
      up4.Boost(betaV);
      down4.Boost(betaV);
      s4.SetE(chargino4.P()/chargino4.M());
      s4.SetVect(chargino4.Vect().Unit()*chargino4.Gamma());
}

double Reweight_T2bW (double thetaChi_eff, double thetaW_eff, std::vector<SUSYGenParticle> genParticles) {
    double weight = 1.;

    unsigned int ngen = genParticles.size();

    for (unsigned int i_stop=0; i_stop<ngen; ++i_stop) {
      // Look for stops
      const SUSYGenParticle& gen = genParticles[i_stop];
      if (abs(gen.pdgId)!=1000006) continue;

      // Look for stop decay products
      int i_b = -1;
      int i_chargino = -1;
      for (unsigned int ig=i_stop+1; ig<ngen; ++ig) {
            const SUSYGenParticle& gen = genParticles[ig];
            if (abs(gen.firstMother)!=i_stop) continue;
            if (abs(gen.pdgId)==5) i_b = ig;
            else if (abs(gen.pdgId)==1000024) i_chargino = ig;
            if (i_b>=0 && i_chargino>=0) break;
      }
      if (i_b<0 || i_chargino<0) continue;

      int i_neutralino = -1;
      int i_W = -1;
      for (unsigned int ig=i_chargino+1; ig<ngen; ++ig) {
            const SUSYGenParticle& gen = genParticles[ig];
            if (abs(gen.firstMother)!=i_chargino) continue;
            if (abs(gen.pdgId)==24) i_W = ig;
            else if (abs(gen.pdgId)==1000022) i_neutralino = ig;
            if (i_W>=0 && i_neutralino>=0) break;
      }
      if (i_W<0 || i_neutralino<0) continue;

      int i_up = -1;
      int i_down = -1;
      for (unsigned int ig=i_W+1; ig<ngen; ++ig) {
            const SUSYGenParticle& gen = genParticles[ig];
            if (abs(gen.firstMother)!=i_W) continue;
            if (abs(gen.pdgId)%2==0) i_up = ig;
            else if (abs(gen.pdgId)%2==1) i_down = ig;
            if (i_up>=0 && i_down>=0) break;
      }
      if (i_up<0 || i_down<0) continue;

      const SUSYGenParticle& gen_stop = genParticles[i_stop];
      const SUSYGenParticle& gen_b = genParticles[i_b];
      const SUSYGenParticle& gen_chargino = genParticles[i_chargino];
      const SUSYGenParticle& gen_W = genParticles[i_W];
      const SUSYGenParticle& gen_neutralino = genParticles[i_neutralino];
      const SUSYGenParticle& gen_up = genParticles[i_up];
      const SUSYGenParticle& gen_down = genParticles[i_down];

      // Fill Lorentz four-vectors
      TLorentzVector stop4, chargino4, b4, neutralino4, W4, up4, down4;

      stop4.SetPtEtaPhiE(gen_stop.pt, gen_stop.eta, gen_stop.phi, gen_stop.energy);
      chargino4.SetPtEtaPhiE(gen_chargino.pt, gen_chargino.eta, gen_chargino.phi, gen_chargino.energy);
      b4.SetPtEtaPhiE(gen_b.pt, gen_b.eta, gen_b.phi, gen_b.energy);
      neutralino4.SetPtEtaPhiE(gen_neutralino.pt, gen_neutralino.eta, gen_neutralino.phi, gen_neutralino.energy);
      W4.SetPtEtaPhiE(gen_W.pt, gen_W.eta, gen_W.phi, gen_W.energy);
      up4.SetPtEtaPhiE(gen_up.pt, gen_up.eta, gen_up.phi, gen_up.energy);
      down4.SetPtEtaPhiE(gen_down.pt, gen_down.eta, gen_down.phi, gen_down.energy);

      // Reference spin four-vector along the chargino direction (filled in after boost)
      TLorentzVector s4;

      // Move everything to the stop center-of-mass frame
      Boost_To_Stop_Rest_Frame(stop4, chargino4, b4, neutralino4, W4, up4, down4, s4);

      double c_L = sin(thetaW_eff);
      double c_R = cos(thetaW_eff);
      double norm_target = Norm(chargino4.M(), neutralino4.M(), W4.M(), c_L, c_R);
      double target = 0;
      for (int hel = -1; hel<2; hel += 2) {
            TLorentzVector t4 = s4*hel;
            TLorentzVector chargino4_plus = chargino4 + t4*chargino4.M();
            TLorentzVector chargino4_minus = chargino4 - t4*chargino4.M();
            target += (1. - chargino4.M()*cos(2*thetaChi_eff)*(b4*t4)/((b4*chargino4) - (b4.M()*chargino4.M())*sin(2*thetaChi_eff)))/2 *
           (8*c_L*c_L*(down4*chargino4_minus)*(up4*neutralino4)
                      + 8*c_R*c_R*(up4*chargino4_plus)*(down4*neutralino4)
                      - 4*c_L*c_R*neutralino4.M()*(pow(W4.M(),2)*chargino4.M()-2*(down4*t4)*(up4*chargino4)+2*(down4*chargino4)*(up4*t4)))/norm_target; 
      }

      weight *= target;

    }

    return weight;

};
