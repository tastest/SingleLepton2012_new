//#include "../../CORE/jetSmearingTools.h"
#include "../../CORE/Thrust.h"
#include "../../CORE/EventShape.h"

#include "StopTreeLooper.h"

//#include "../Core/STOPT.h"
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
#include "Riostream.h"
#include "TFitter.h"
#include "TRandom3.h"

#include <algorithm>
#include <utility>
#include <map>
#include <set>
#include <list>

using namespace Stop;

StopTreeLooper::StopTreeLooper()
{
  m_outfilename_ = "histos.root";
  // t1metphicorr = -9999.;
  // t1metphicorrphi = -9999.;
  // t1metphicorrmt = -9999.;
  // min_mtpeak = -9999.;
  // max_mtpeak = -9999.; 
}

StopTreeLooper::~StopTreeLooper()
{
}

void StopTreeLooper::setOutFileName(string filename)
{
  m_outfilename_ = filename;

}
//--------------------------------------------------------------------

struct DorkyEventIdentifier {
  // this is a workaround for not having unique event id's in MC
  unsigned long int run, event,lumi;
  bool operator < (const DorkyEventIdentifier &) const;
  bool operator == (const DorkyEventIdentifier &) const;
};

//--------------------------------------------------------------------

bool DorkyEventIdentifier::operator < (const DorkyEventIdentifier &other) const
{
  if (run != other.run)
    return run < other.run;
  if (event != other.event)
    return event < other.event;
  if(lumi != other.lumi)
    return lumi < other.lumi;
  return false;
}

//--------------------------------------------------------------------

bool DorkyEventIdentifier::operator == (const DorkyEventIdentifier &other) const
{
  if (run != other.run)
    return false;
  if (event != other.event)
    return false;
  return true;
}

//--------------------------------------------------------------------

std::set<DorkyEventIdentifier> already_seen; 
bool is_duplicate (const DorkyEventIdentifier &id) {
  std::pair<std::set<DorkyEventIdentifier>::const_iterator, bool> ret =
    already_seen.insert(id);
  return !ret.second;
}

//--------------------------------------------------------------------


std::set<DorkyEventIdentifier> events_hcallasercalib; 
int load_badhcallaserevents  () {

  ifstream in;
  in.open("../Core/badhcallaser_events.txt");

  int run, event, lumi;
  int nlines = 0;

  while (1) {
    in >> run >> event >> lumi;
    if (!in.good()) break;
    nlines++;
    DorkyEventIdentifier id = {run, event, lumi };
    events_hcallasercalib.insert(id);
  }
  printf(" found %d bad events \n",nlines);

  in.close();

  return 0;

}

bool is_badHcalLaserEvent (const DorkyEventIdentifier &id) {
  if (events_hcallasercalib.find(id) != events_hcallasercalib.end()) return true;
  return false;
}


//--------------------------------------------------------------------

std::set<DorkyEventIdentifier> events_lasercalib; 
int load_badlaserevents  () {

  ifstream in;
  in.open("../Core/badlaser_events.txt");

  int run, event, lumi;
  int nlines = 0;

  while (1) {
    in >> run >> event >> lumi;
    if (!in.good()) break;
    nlines++;
    DorkyEventIdentifier id = {run, event, lumi };
    events_lasercalib.insert(id);
  }
  printf(" found %d bad events \n",nlines);

  in.close();

  return 0;

}

bool is_badLaserEvent (const DorkyEventIdentifier &id) {
  if (events_lasercalib.find(id) != events_lasercalib.end()) return true;
  return false;
}

void StopTreeLooper::loop(TChain *chain, TString name)
{
  TRandom3 r;

  printf("[StopTreeLooper::loop] %s\n", name.Data());

  load_badlaserevents  ();
  load_badhcallaserevents();


  //---------------------------------
  // check for valid chain
  //---------------------------------

  TObjArray *listOfFiles = chain->GetListOfFiles();
  TIter fileIter(listOfFiles);
  if (listOfFiles->GetEntries() == 0) {
    cout << "[StopTreeLooper::loop] no files in chain" << endl;
    return;
  }

  //---------------------------------
  // set up histograms
  //---------------------------------

  gROOT->cd();

  makeTree(name.Data(), chain);

  TH2F *h_nsig, *h_nsig25, *h_nsig75 ;
  
  if( name.Contains("T2") ){
    char* h_nsig_filename = "";

    if( name.Contains("T2tt") ){
      h_nsig_filename = "/nfs-3/userdata/stop/cms2V05-03-18_stoplooperV00-02-07/crabT2tt_3/myMassDB_T2tt.root";
    }

    if( name.Contains("T2bw") ){
      h_nsig_filename = "/nfs-3/userdata/stop/cms2V05-03-18_stoplooperV00-02-07/crabT2bw_3/myMassDB_T2bw_fine.root";
    }

    cout << "[StopTreeLooper::loop] opening mass TH2 file " << h_nsig_filename << endl;

    TFile *f_nsig = TFile::Open(h_nsig_filename);

    if( f_nsig == 0 ) cout << "ERROR! didn't find file " << h_nsig_filename << endl;

    h_nsig = (TH2F*) f_nsig->Get("masses");

    if( name.Contains("T2bw") ){
      h_nsig25 = (TH2F*) f_nsig->Get("masses25");
      h_nsig75 = (TH2F*) f_nsig->Get("masses75");
    }
    

  }

      float apply_mva = true;
      TMVA::Reader *reader;
      if ( apply_mva ) {
          reader = new TMVA::Reader( "!Color:!Silent" );
          reader->AddVariable("met", &met_); 
          reader->AddVariable("lep1pt", &lep1pt_); 
          reader->AddVariable("chi2minprob", &chi2minprob_); 
          reader->AddVariable("mt2wmin", &mt2wmin_); 
          reader->AddVariable("htssm/(htosm+htssm)", &htratiom_); 
          reader->AddVariable("dphimjmin", &dphimjmin_); 
          reader->AddVariable("pt_b", &pt_b_); 

          TString dir    = "/home/users/magania/stop/SingleLepton2012/MVA/weights/";
          TString prefix = "classification_T2tt_8_BDT";

          TString weightfile = dir + prefix + TString(".weights.xml");
          reader->BookMVA( "BDT1" , weightfile );
      } 

  // TFile* vtx_file = TFile::Open("vtxreweight/vtxreweight_Summer12_DR53X-PU_S10_9p7ifb_Zselection.root");
  // if( vtx_file == 0 ){
  //   cout << "vtxreweight error, couldn't open vtx file. Quitting!"<< endl;
  //   exit(0);
  // }

  // TH1F* h_vtx_wgt = (TH1F*)vtx_file->Get("hratio");
  // h_vtx_wgt->SetName("h_vtx_wgt");

  //
  // file loop
  //

  unsigned int nEventsPass=0;
  unsigned int nEventsChain=0;
  unsigned int nEvents = chain->GetEntries();
  nEventsChain = nEvents;
  ULong64_t nEventsTotal = 0;
  //  int i_permille_old = 0;

  bool isData = name.Contains("data") ? true : false;

  cout << "[StopTreeLooper::loop] running over chain with total entries " << nEvents << endl;

  while (TChainElement *currentFile = (TChainElement*)fileIter.Next()) {

    //---------------------------------
    // load the stop baby tree
    //---------------------------------
    TFile *file = new TFile( currentFile->GetTitle() );
    TTree *tree = (TTree*)file->Get("t");
    stopt.Init(tree);

    //---------------------------------
    // event loop
    //---------------------------------

    ULong64_t nEvents = tree->GetEntries();
    //nEvents = 1000;

    for(ULong64_t event = 0; event < nEvents; ++event) {
      stopt.GetEntry(event);

      //---------------------------------
      // increment counters
      //---------------------------------

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

      //------------------------------------------ 
      // skip duplicates
      //------------------------------------------ 

      if( isData ) {
        DorkyEventIdentifier id = {stopt.run(), stopt.event(), stopt.lumi() };
        if (is_duplicate(id) ){
          continue;
        }
	if (is_badLaserEvent(id) ){
	  //std::cout<<"Removed bad laser calibration event:" << stopt.run() <<"   "<< stopt.event.() <<"\n";
	  continue;
	}
	if (is_badHcalLaserEvent(id) ){
	  std::cout<< "Removed bad hcal laser calibration event:" << stopt.run() << "   " << stopt.event() <<"\n";
	  continue;
	}

      }

      //------------------------------------------ 
      // event weight
      //------------------------------------------ 

      float evtweight = isData ? 1. : ( stopt.weight() * 19.5 * stopt.nvtxweight() * stopt.mgcor() );

      if( name.Contains("T2tt") ) {
	int bin = h_nsig->FindBin(stopt.mg(),stopt.ml());
	float nevents = h_nsig->GetBinContent(bin);
	evtweight =  stopt.xsecsusy() * 1000.0 / nevents; 

	//cout << "mg " << tree->mg_ << " ml " << tree->ml_ << " bin " << bin << " nevents " << nevents << " xsec " << tree->xsecsusy_ << " weight " << evtweight << endl;
      }

      if( name.Contains("T2bw") ) {
        float nevents = 0;
        if ( stopt.x() == 0.25 ){
	   int bin = h_nsig25->FindBin(stopt.mg(),stopt.ml());
   	   nevents = h_nsig25->GetBinContent(bin);
        } else if ( stopt.x() == 0.50 ){
	   int bin = h_nsig->FindBin(stopt.mg(),stopt.ml());
   	   nevents = h_nsig->GetBinContent(bin);
        } else if ( stopt.x() == 0.75 ){
	   int bin = h_nsig75->FindBin(stopt.mg(),stopt.ml());
   	   nevents = h_nsig75->GetBinContent(bin);
        }

        assert ( nevents > 0 );

	evtweight =  stopt.xsecsusy() * 1000.0 / nevents; 

//	cout << "mg " << stopt.mg() << " ml " << stopt.ml() << " x " << stopt.x() << " nevents " << nevents << " xsec " << stopt.xsecsusy() << " weight " << evtweight << endl;
      }

      float trigweight   = isData ? 1. : getsltrigweight(stopt.id1(), stopt.lep1().Pt(), stopt.lep1().Eta());
      float trigweightdl = isData ? 1. : getdltrigweight(stopt.id1(), stopt.id2());

      //------------------------------------------ 
      // selection criteria
      //------------------------------------------ 

      if ( !passEvtSelection(name) ) continue; // >=1 lepton, rho cut, MET filters, veto 2 nearby leptons
      if ( getNJets() < 4          ) continue; // >=4 jets
      if ( stopt.t1metphicorr() < 50.0    ) continue; // MET > 50 GeV

      //------------------------------------------ 
      // get list of candidates
      //------------------------------------------ 
      
      float met = stopt.t1metphicorr();
      float metphi = stopt.t1metphicorrphi();

      // get list of candidates

      PartonCombinatorics pc (stopt.pfjets(), stopt.pfjets_csv(), stopt.pfjets_sigma(), stopt.pfjets_mc3(), stopt.lep1(), met, metphi, isData);
      MT2CHI2 mc = pc.getMt2Chi2();

      //------------------------------------------ 
      // event shapes
      //------------------------------------------ 

      std::vector<LorentzVector>  myPfJets = stopt.pfjets();
      std::vector<LorentzVector>  jetVector;
      std::vector<LorentzVector>  jetLeptonVector;
      std::vector<LorentzVector>  jetLeptonMetVector;

      jetVector.clear();
      jetLeptonVector.clear();
      jetLeptonMetVector.clear();

      float HT_SSL=0;
      float HT_OSL=0;

      float HT_SSM=0;
      float HT_OSM=0;


      LorentzVector lep1_p4( stopt.lep1().px() , stopt.lep1().py() , 0 , stopt.lep1().E() ); 

      LorentzVector met_p4( stopt.t1metphicorr() * cos( stopt.t1metphicorr()) , 
			    stopt.t1metphicorr() * sin( stopt.t1metphicorr()) , 
			    0,
			    stopt.t1metphicorr()
			    );

      jetLeptonVector.push_back(lep1_p4);

      jetLeptonMetVector.push_back(lep1_p4);
      jetLeptonMetVector.push_back(met_p4);

      for ( unsigned int i=0; i<myPfJets.size() ; i++) {

	if( myPfJets.at(i).pt()<30 )          continue;
	if( fabs(myPfJets.at(i).eta())>3.0 )  continue;

	LorentzVector jet_p4( myPfJets.at(i).px() , myPfJets.at(i).py() , 0 , myPfJets.at(i).E() );
    
	// jetVector.push_back(myPfJets->at(i));
	// jetLeptonVector.push_back(myPfJets->at(i));
	// jetLeptonMetVector.push_back(myPfJets->at(i));

	jetVector.push_back(jet_p4);
	jetLeptonVector.push_back(jet_p4);
	jetLeptonMetVector.push_back(jet_p4);

	float dPhiL = getdphi(stopt.lep1().Phi(), myPfJets.at(i).phi() );
	float dPhiM = getdphi(stopt.t1metphicorrphi(), myPfJets.at(i).phi() );
    
	if(dPhiL<(3.14/2))  HT_SSL=HT_SSL+myPfJets.at(i).pt();
	if(dPhiL>=(3.14/2)) HT_OSL=HT_OSL+myPfJets.at(i).pt();

	if(dPhiM<(3.14/2))  HT_SSM=HT_SSM+myPfJets.at(i).pt();
	if(dPhiM>=(3.14/2)) HT_OSM=HT_OSM+myPfJets.at(i).pt();

      }
  
      // from jets only
      Thrust thrust(jetVector);
      EventShape eventshape(jetVector);

      // from jets and lepton
      Thrust thrustl(jetLeptonVector);
      EventShape eventshapel(jetLeptonVector);

      // from jets and lepton and MET
      Thrust thrustlm(jetLeptonMetVector);
      EventShape eventshapelm(jetLeptonMetVector);

      //------------------------------------------ 
      // datasets bit
      //------------------------------------------ 

      bool dataset_1l=false;

      if((isData) && name.Contains("muo") && (abs(stopt.id1()) == 13 ))  dataset_1l=true;
      if((isData) && name.Contains("ele") && (abs(stopt.id1()) == 11 ))  dataset_1l=true;

      if(!isData) dataset_1l=true;

      bool dataset_CR4=false;

      if((isData) && name.Contains("dimu") && (abs(stopt.id1()) == 13 ) && (abs(stopt.id2())==13) && (fabs( stopt.dilmass() - 91.) > 15.)) dataset_CR4=true;
      if((isData) && name.Contains("diele") && (abs(stopt.id1()) == 11 ) && (abs(stopt.id2())==11) && (fabs( stopt.dilmass() - 91.) > 15.)) dataset_CR4=true;
      if((isData) && name.Contains("mueg") && abs(stopt.id1()) != abs(stopt.id2())) dataset_CR4=true;

      if(!isData) dataset_CR4=true;

      //------------------------------------------ 
      // variables to add to baby
      //------------------------------------------ 
      
      initBaby(); // set all branches to -1

      vector<int> indexBJets=getBJetIndex(0.679,-1,-1);
      if(indexBJets.size()>0) pt_b_ = myPfJets.at(indexBJets.at(0)).pt();

      int J1Index=leadingJetIndex( stopt.pfjets(), -1, -1);
      int J2Index=leadingJetIndex( stopt.pfjets(), J1Index, -1);
      pt_J1_ = myPfJets.at(J1Index).pt();
      pt_J2_ = myPfJets.at(J2Index).pt();

      // which selections are passed
      sig_        = ( dataset_1l && passOneLeptonSelection(isData) && indexBJets.size()>=1 ) ? 1 : 0; // pass signal region preselection
      cr1_        = ( dataset_1l && passOneLeptonSelection(isData) && indexBJets.size()==0 ) ? 1 : 0; // pass CR1 (b-veto) control region preselection
      cr4_        = ( dataset_CR4 && passDileptonSelection(isData) && indexBJets.size()==1) ? 1 : 0; // pass CR4 (dilepton) control region preselection
      cr5_        = ( dataset_1l && passLepPlusIsoTrkSelection(isData) && indexBJets.size()==1) ? 1 : 0; // pass CR1 (lepton+isotrack) control region preselection

      // kinematic variables
      met_        = stopt.t1metphicorr();       // MET (type1, MET-phi corrections)
      mt_         = stopt.t1metphicorrmt();     // MT (type1, MET-phi corrections)

      // chi2 and MT2 variables
      chi2min_			= mc.one_chi2;               // minimum chi2 
      chi2minprob_		= TMath::Prob(chi2min_,1);   // probability of minimum chi2
      chi2min_mt2b_		= mc.two_mt2b;               // minimum MT2b consistent with chi2min
      chi2min_mt2bl_		= mc.two_mt2bl;              // minimum MT2bl consistent with chi2min
      chi2min_mt2w_		= mc.two_mt2w;               // minimum MT2w consistent with chi2min
      mt2bmin_			= mc.three_mt2b;             // minimum MT2b
      mt2blmin_			= mc.three_mt2bl;            // minimum MT2bl
      mt2wmin_			= mc.three_mt2w;             // minimum MT2w
      mt2bmin_chi2_		= mc.four_chi2b;             // minimum chi2 consistent with mt2bmin
      mt2blmin_chi2_		= mc.four_chi2bl;            // minimum chi2 consistent with mt2blmin
      mt2wmin_chi2_		= mc.four_chi2w;             // minimum chi2 consistent with mt2wmin
      mt2bmin_chi2prob_		= mc.four_chi2b;             // probability of minimum chi2 consistent with mt2bmin
      mt2blmin_chi2prob_	= mc.four_chi2bl;            // probability of minimum chi2 consistent with mt2blmin
      mt2wmin_chi2prob_		= mc.four_chi2w;             // probability of minimum chi2 consistent with mt2wmin
      ncand_			= pc.b_candidates_.size();   // number of candidates consisting with btagging info

      // weights
      weight_     = evtweight;                    // event weight
      sltrigeff_  = trigweight;                   // trigger weight (single lepton)
      dltrigeff_  = trigweightdl;                 // trigger weight (dilepton)

      // hadronic variables 
      nb_         = indexBJets.size();            // nbjets (pT > 30, CSVM)
      njets_      = getNJets();             // njets (pT > 30, eta < 2.5)

      // lepton variables
      passisotrk_   = passIsoTrkVeto() ? 1 : 0; // is there an isolated track? (pT>10 GeV, reliso<0.1)
      passisotrkv2_ = passIsoTrkVeto_v2() ? 1 : 0; // is there an isolated track? (pT>10 GeV, reliso<0.1)
      nlep_       = stopt.ngoodlep();              // number of analysis selected leptons
 
      lep1pt_     = stopt.lep1().pt();             // 1st lepton pt
      lep1eta_    = fabs( stopt.lep1().eta() );    // 1st lepton eta

      if( nlep_ > 1 ){
	lep2pt_    = stopt.lep1().pt();            // 2nd lepton pt
	lep2eta_   = stopt.lep1().eta();           // 2nd lepton eta
	dilmass_   = stopt.dilmass();              // dilepton mass
      }
      
      // susy vars
      mstop_       = stopt.mg();                   // stop mass
      mlsp_        = stopt.ml();                   // LSP mass
      x_           = stopt.x();                    // chargino mass parameter x

      // event shapes (from jets only)
      thrjet_    = thrust.thrust();
      apljet_    = eventshape.aplanarity();
      sphjet_    = eventshape.sphericity();
      cirjet_    = eventshape.circularity();

      // event shapes (from jets and lepton)
      thrjetl_   = thrustl.thrust();
      apljetl_   = eventshapel.aplanarity();
      sphjetl_   = eventshapel.sphericity();
      cirjetl_   = eventshapel.circularity();

      // event shapes (from jets, lepton, and MET)
      thrjetlm_   = thrustlm.thrust();
      apljetlm_   = eventshapelm.aplanarity();
      sphjetlm_   = eventshapelm.sphericity();
      cirjetlm_   = eventshapelm.circularity();

      // event shapes: HT in hemisphers
      htssl_      = HT_SSL;
      htosl_      = HT_OSL;
      htratiol_   = HT_SSL / (HT_OSL + HT_SSL);

      htssm_      = HT_SSM;
      htosm_      = HT_OSM;
      htratiom_   = HT_SSM / (HT_OSM + HT_SSM);

      dphimj1_    = getdphi(stopt.t1metphicorrphi(), jetVector.at(0).phi() );
      dphimj2_    = getdphi(stopt.t1metphicorrphi(), jetVector.at(1).phi() );
      dphimjmin_  = TMath::Min( dphimj1_ , dphimj2_ );

      rand_       = r.Uniform(0.0,1.0);

      if ( apply_mva ) {
          bdt_ = reader->EvaluateMVA( "BDT1" );
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

  string revision = "$Revision: 1.26 $";
  string revision_no = revision.substr(11, revision.length() - 13);
  outFile_   = new TFile(Form("output/%s_mini_%s.root",prefix,revision_no.c_str()), "RECREATE");
  outFile_->cd();

//  outTree_ = new TTree("t","Tree");

  outTree_ = chain->CloneTree(0);

  outTree_->Branch("lep1pt"           ,        &lep1pt_          ,         "lep1pt/F"		);
  outTree_->Branch("lep1eta"          ,        &lep1eta_         ,         "lep1eta/F"		);
  outTree_->Branch("sig"              ,        &sig_             ,         "sig/I"		);
  outTree_->Branch("cr1"              ,        &cr1_             ,         "cr1/I"		);
  outTree_->Branch("cr4"              ,        &cr4_             ,         "cr4/I"		);
  outTree_->Branch("cr5"              ,        &cr5_             ,         "cr5/I"		);
  outTree_->Branch("met"              ,        &met_             ,         "met/F"		);
  outTree_->Branch("mt"               ,        &mt_              ,         "mt/F"		);
  outTree_->Branch("weight"           ,        &weight_          ,         "weight/F"		);
  outTree_->Branch("sltrigeff"        ,        &sltrigeff_       ,         "sltrigeff/F"	);
  outTree_->Branch("dltrigeff"        ,        &dltrigeff_       ,         "dltrigeff/F"	);
  outTree_->Branch("nb"               ,        &nb_              ,         "nb/I"		);
  outTree_->Branch("njets"            ,        &njets_           ,         "njets/I"		);
  outTree_->Branch("passisotrk"       ,        &passisotrk_      ,         "passisotrk/I"	);
  outTree_->Branch("passisotrkv2"     ,        &passisotrkv2_    ,         "passisotrkv2/I"	);
  outTree_->Branch("nlep"             ,        &nlep_            ,         "nlep/I"		);
  outTree_->Branch("lep1pt"           ,        &lep1pt_          ,         "lep1pt/F"		);
  outTree_->Branch("lep1eta"          ,        &lep1eta_         ,         "lep1eta/F"		);
  outTree_->Branch("lep2pt"           ,        &lep2pt_          ,         "lep2pt/F"		);
  outTree_->Branch("lep2eta"          ,        &lep2eta_         ,         "lep2eta/F"		);
  outTree_->Branch("dilmass"          ,        &dilmass_         ,         "dilmass/F"		);
  outTree_->Branch("mstop"            ,        &mstop_           ,         "mstop/F"		);
  outTree_->Branch("mlsp"             ,        &mlsp_            ,         "mlsp/F"		);
  outTree_->Branch("x"                ,        &x_               ,         "x/F"		);
  outTree_->Branch("chi2min"          ,        &chi2min_         ,         "chi2min/F"          );
  outTree_->Branch("chi2minprob"      ,        &chi2minprob_     ,         "chi2minprob/F"      );
  outTree_->Branch("chi2min_mt2b"     ,        &chi2min_mt2b_    ,         "chi2min_mt2b/F"     );  
  outTree_->Branch("chi2min_mt2bl"    ,        &chi2min_mt2bl_   ,         "chi2min_mt2bl/F"    );  
  outTree_->Branch("chi2min_mt2w"     ,        &chi2min_mt2w_    ,         "chi2min_mt2w/F"     );  
  outTree_->Branch("mt2bmin"          ,        &mt2bmin_         ,         "mt2bmin/F"          );       
  outTree_->Branch("mt2blmin"         ,        &mt2blmin_        ,         "mt2blmin/F"         );       
  outTree_->Branch("mt2wmin"          ,        &mt2wmin_         ,         "mt2wmin/F"          );       
  outTree_->Branch("mt2bmin_chi2"     ,        &mt2bmin_chi2_    ,         "mt2bmin_chi2/F"     );       
  outTree_->Branch("mt2blmin_chi2"    ,        &mt2blmin_chi2_   ,         "mt2blmin_chi2/F"    );       
  outTree_->Branch("mt2wmin_chi2"     ,        &mt2wmin_chi2_    ,         "mt2wmin_chi2/F"     );       
  outTree_->Branch("mt2bmin_chi2prob" ,        &mt2bmin_chi2prob_    ,         "mt2bmin_chi2prob/F"     );       
  outTree_->Branch("mt2blmin_chi2prob",        &mt2blmin_chi2prob_   ,         "mt2blmin_chi2prob/F"    );       
  outTree_->Branch("mt2wmin_chi2prob" ,        &mt2wmin_chi2prob_    ,         "mt2wmin_chi2prob/F"     );       
  outTree_->Branch("ncand"            ,        &ncand_           ,         "ncand/F"            );       
  outTree_->Branch("thrjet"           ,        &thrjet_          ,         "thrjet/F"           );
  outTree_->Branch("sphjet"           ,        &sphjet_          ,         "sphjet/F"           );
  outTree_->Branch("apljet"           ,        &apljet_          ,         "apljet/F"           );
  outTree_->Branch("cirjet"           ,        &cirjet_          ,         "cirjet/F"           );
  outTree_->Branch("thrjetl"          ,        &thrjetl_         ,         "thrjetl/F"          );
  outTree_->Branch("sphjetl"          ,        &sphjetl_         ,         "sphjetl/F"          );
  outTree_->Branch("apljetl"          ,        &apljetl_         ,         "apljetl/F"          );
  outTree_->Branch("cirjetl"          ,        &cirjetl_         ,         "cirjetl/F"          );
  outTree_->Branch("thrjetlm"         ,        &thrjetlm_        ,         "thrjetlm/F"         );
  outTree_->Branch("sphjetlm"         ,        &sphjetlm_        ,         "sphjetlm/F"         );
  outTree_->Branch("apljetlm"         ,        &apljetlm_        ,         "apljetlm/F"         );
  outTree_->Branch("cirjetlm"         ,        &cirjetlm_        ,         "cirjetlm/F"         );
  outTree_->Branch("htssl"            ,        &htssl_           ,         "htssl/F"            );
  outTree_->Branch("htosl"            ,        &htosl_           ,         "htosl/F"            );
  outTree_->Branch("htratiol"         ,        &htratiol_        ,         "htraiol/F"          );
  outTree_->Branch("htssm"            ,        &htssm_           ,         "htssm/F"            );
  outTree_->Branch("htosm"            ,        &htosm_           ,         "htosm/F"            );
  outTree_->Branch("htratiom"         ,        &htratiom_        ,         "htraiom/F"          );
  outTree_->Branch("dphimj1"          ,        &dphimj1_         ,         "dphimj1/F"          );
  outTree_->Branch("dphimj2"          ,        &dphimj2_         ,         "dphimj2/F"          );
  outTree_->Branch("dphimjmin"        ,        &dphimjmin_       ,         "dphimjmin/F"        );
  outTree_->Branch("pt_b"             ,        &pt_b_            ,         "pt_b/F"             );
  outTree_->Branch("pt_J1"            ,        &pt_J1_           ,         "pt_J1/F"            );
  outTree_->Branch("pt_J2"            ,        &pt_J2_           ,         "pt_J2/F"            );
  outTree_->Branch("rand"             ,        &rand_            ,         "rand/F"             );

  outTree_->Branch("bdt"             ,        &bdt_            ,           "bdt/F"             );

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

  // kinematic variables
  met_        = -1.0;
  mt_         = -1.0;
  chi2_       = -1.0;
  mt2w_       = -1.0;
  mt2b_       = -1.0;
  mt2bl_      = -1.0;

  // "best" chi2 and MT2 variables
  chi2min_        = -1.0;
  chi2min_mt2b_   = -1.0;  
  chi2min_mt2bl_  = -1.0; 
  chi2min_mt2w_   = -1.0;  
  mt2bmin_        = -1.0;       
  mt2blmin_       = -1.0;      
  mt2wmin_        = -1.0;       
  mt2bmin_chi2_   = -1.0;  
  mt2blmin_chi2_  = -1.0; 
  mt2wmin_chi2_   = -1.0;  

  // event shapes
  thrjet_     = -1.0;
  apljet_     = -1.0;
  sphjet_     = -1.0;
  cirjet_     = -1.0;

  thrjetl_    = -1.0;
  apljetl_    = -1.0;
  sphjetl_    = -1.0;
  cirjetl_    = -1.0;

  thrjetlm_   = -1.0;
  apljetlm_   = -1.0;
  sphjetlm_   = -1.0;
  cirjetlm_   = -1.0;

  htssl_      = -1.0;
  htosl_      = -1.0;
  htratiol_   = -1.0;
  htssm_      = -1.0;
  htosm_      = -1.0;
  htratiom_   = -1.0;

  // weights
  weight_     = -1.0;
  sltrigeff_  = -1.0;
  dltrigeff_  = -1.0;

  // hadronic variables
  nb_         = -1;
  njets_      = -1;

  // lepton variables
  passisotrk_ = -1;
  nlep_       = -1;
  lep1pt_     = -1.0;
  lep1eta_    = -1.0;
  lep2pt_     = -1.0;
  lep2eta_    = -1.0;
  dilmass_    = -1.0;

  // susy variables
  mstop_      = -1.0;
  mlsp_       = -1.0;
  x_          = -1.0;

  // jet kinematics
  pt_b_ = -1.0;
  pt_J1_ = -1.0;
  pt_J2_ = -1.0;

  rand_       = -1.0;
}

