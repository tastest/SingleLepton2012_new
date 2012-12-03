
#include "StopTreeLooper.h"
#include "../Core/StopTree.h"
#include "../Plotting/PlotUtilities.h"

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


float StopTreeLooper::vtxweight_n( const int nvertices, TH1F *hist, bool isData ) 
{

  if( isData ) return 1;

  int nvtx = nvertices;
  if( nvtx > hist->GetNbinsX() )
        nvtx = hist->GetNbinsX();

  float weight = 0;
  weight = hist->GetBinContent( hist->FindBin(nvtx) );
  if( weight <= 0 ) //we don't want to kill events bc they have no weight
        weight = 1.;
  //  cout << "nvtx " << nvtx << " weight " << weight << endl;
  return weight;

}

float StopTreeLooper::getdphi( float phi1 , float phi2 ) 
{
  float dphi = fabs( phi1 - phi2 );
  if( dphi > TMath::Pi() ) dphi = TMath::TwoPi() - dphi;
  return dphi;
}

float StopTreeLooper::getMT( float pt1 , float phi1 , float pt2 , float phi2 ) 
{

  float dphi = getdphi(phi1, phi2);
  return sqrt( 2 * ( pt1 * pt2 * (1 - cos( dphi ) ) ) );

}

float StopTreeLooper::dRbetweenVectors(LorentzVector vec1, 
				       LorentzVector vec2 )
{ 

  float dphi = std::min(::fabs(vec1.Phi() - vec2.Phi()), 2 * M_PI - fabs(vec1.Phi() - vec2.Phi()));
  float deta = vec1.Eta() - vec2.Eta();

  return sqrt(dphi*dphi + deta*deta);

}


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



//--------------------------------------------------------------------

void StopTreeLooper::loop(TChain *chain, TString name)
{

  //------------------------------
  // check for valid chain
  //------------------------------

  printf("[StopTreeLooper::loop] %s\n", name.Data());

  load_badlaserevents();

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
  unsigned int nEventsTotal = 0;
  int i_permille_old = 0;

  bool isData = name.Contains("data") ? true : false;

  while (TChainElement *currentFile = (TChainElement*)fileIter.Next()) {

    //----------------------------
    // load the stop baby tree
    //----------------------------

    StopTree *tree = new StopTree();
    tree->LoadTree(currentFile->GetTitle());
    tree->InitTree();

    //----------------------------
    // event loop
    //----------------------------

    ULong64_t nEvents = tree->tree_->GetEntries();
    for(ULong64_t event = 0; event < nEvents; ++event) {
      tree->tree_->GetEntry(event);

      //----------------------------
      // increment counters
      //----------------------------

      ++nEventsTotal;
      if (nEventsTotal%10000==0) {
	int i_permille = (int)floor(1000 * nEventsTotal / float(nEventsChain));
	//if (i_permille != i_permille_old) {//this prints too often!
	// xterm magic from L. Vacavant and A. Cerri
	if (isatty(1)) {
	  printf("\015\033[32m ---> \033[1m\033[31m%4.1f%%"
		 "\033[0m\033[32m <---\033[0m\015", i_permille/10.);
	  fflush(stdout);
	}
	i_permille_old = i_permille;
      }

      //---------------------
      // skip duplicates
      //---------------------

      if( isData ) {
        DorkyEventIdentifier id = {tree->run_,tree->event_, tree->lumi_ };
        if (is_duplicate(id) ){
          continue;
        }
	if (is_badLaserEvent(id) ){
	  //std::cout<<"Removed bad laser calibration event:" <<tree->run_<<"   "<<tree->event_<<"\n";
	  continue;
	}
      }

      //---------------------------------------------------------------------------- 
      // determine event weight
      // make 2 example histograms of nvtx and corresponding weight
      //---------------------------------------------------------------------------- 

      float evtweight = isData ? 1. : ( tree->weight_ * 9.708 * tree->nvtxweight_ * tree->mgcor_ );
      // to reweight from file - also need to comment stuff before
      //      float vtxweight = vtxweight_n( tree->nvtx_, h_vtx_wgt, isData );

      plot1D("h_vtx",       tree->nvtx_,       evtweight, h_1d, 40, 0, 40);
      plot1D("h_vtxweight", tree->nvtxweight_, evtweight, h_1d, 41, -4., 4.);

      //----------------------------------------------------------------------------
      // apply preselection:
      // rho 0-40 GeV, MET filters, >=1 good lepton, veto 2 leptons dR < 0.1
      //----------------------------------------------------------------------------

      if ( !passEvtSelection(tree, name) ) continue;

      //----------------------------------------------------------------------------
      // Function to perform MET phi corrections on-the-fly
      // Default branches are: tree->t1metphicorr_ and tree->t1metphicorrmt_
      //----------------------------------------------------------------------------

      // pair<float, float> p_t1metphicorr = 
      // 	getPhiCorrMET( tree->t1met10_, tree->t1met10phi_, tree->nvtx_, !isData);
      // t1metphicorr    = p_t1metphicorr.first;
      // t1metphicorrphi = p_t1metphicorr.second;
      // t1metphicorrmt  = getMT( tree->lep1_.Pt() , tree->lep1_.Phi() , t1metphicorr , t1metphicorrphi );  


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

bool StopTreeLooper::passEvtSelection(const StopTree *sTree, TString name) 
{

  //rho requirement
  if ( sTree->rhovor_<0. || sTree->rhovor_>=40. ) return false;

  if (!name.Contains("T2")) {
      //met filters
      if ( sTree->csc_      != 0 ) return false;
      if ( sTree->hbhe_     != 1 ) return false;
      if ( sTree->hcallaser_!= 1 ) return false;
      if ( sTree->ecaltp_   != 1 ) return false;
      if ( sTree->trkfail_  != 1 ) return false;
      if ( sTree->eebadsc_  != 1 ) return false;
      if ( sTree->hbhenew_  != 1 ) return false;
    }

  //at least 1 lepton
  if ( sTree->ngoodlep_ < 1 ) return false;

  //if have more than 1 lepton, remove cases where have 2 close together
  if ( sTree->ngoodlep_ > 1 && 
       dRbetweenVectors( sTree->lep1_ ,  sTree->lep2_ )<0.1 ) return false;

  return true;

}

bool StopTreeLooper::passOneLeptonSelection(const StopTree *sTree, bool isData) 
{
  //single lepton selection for 8 TeV 53 analysis
  if ( !passSingleLeptonSelection(sTree, isData) ) return false;

  //pass isolated track veto
  //unfortunately changed default value to 9999.
  if ( sTree->pfcandpt10_ <9998. && sTree->pfcandiso10_ < 0.1 ) return false;

  return true;

}

bool StopTreeLooper::passTwoLeptonSelection(const StopTree *sTree, bool isData) 
{
  //single lepton selection for 8 TeV 53 analysis
  if ( !passDileptonSelection(sTree, isData) ) return false;

  //apply isolated track veto in addition to 2 leptons
  //default value for this one is -999
  if ( sTree->trkpt10loose_ >0. && sTree->trkreliso10loose_ < 0.1 ) return false;

  return true;

}

bool StopTreeLooper::passIsoTrkVeto(const StopTree *sTree) 
{

  //pass isolated track veto
  //unfortunately changed default value to 9999.
  if ( sTree->pfcandpt10_ <9998. && sTree->pfcandiso10_ < 0.1 ) return false;

  return true;

}

bool StopTreeLooper::passSingleLeptonSelection(const StopTree *sTree, bool isData) 
{
  //single lepton selection for 8 TeV 53 analysis

  //at least one lepton
  if ( sTree->ngoodlep_ < 1 ) return false;

  //lepton flavor - trigger, pt and eta requirements
  if ( sTree->lep1_.Pt() < 30 )          return false;
  if ( fabs( sTree->pflep1_.Pt() - sTree->lep1_.Pt() ) > 10. )  return false;
  if ( ( sTree->isopf1_ * sTree->lep1_.Pt() ) > 5. )  return false; 
  
  if ( sTree->leptype_ == 0 ) {

    //pass trigger if data - single electron
    if ( isData && sTree->ele27wp80_ != 1 ) return false;
    //    if ( isData && sTree->trgel1_ != 1 )  return false;
    
    //barrel only electrons
    if ( fabs(sTree->lep1_.Eta() ) > 1.4442) return false;
    if ( sTree->eoverpin_ > 4. ) return false;
    

  } else if ( sTree->leptype_ == 1 ) {

    //pass trigger if data - single muon
    if ( isData && sTree->isomu24_ != 1 ) return false;
    //    if ( isData && sTree->trgmu1_ != 1 )  return false;
    
    if ( fabs(sTree->lep1_.Eta() ) > 2.1)  return false;

  }

  return true;

}


bool StopTreeLooper::passDileptonSelection(const StopTree *sTree, bool isData) 
{
  //two lepton selection for 8 TeV 53 analysis

  //exactly 2 leptons
  if ( sTree->ngoodlep_ != 2 ) return false;

  //opposite sign
  if ( sTree->id1_*sTree->id2_>0 ) return false;

  //pass trigger if data - dilepton
  if ( isData && sTree->mm_ != 1 && sTree->me_ != 1 
       && sTree->em_ != 1 && sTree->ee_ != 1 ) return false;

  //passes pt and eta requirements
  if ( sTree->lep1_.Pt() < 20 )          return false;
  if ( sTree->lep2_.Pt() < 20 )          return false;
  if ( fabs(sTree->lep1_.Eta() ) > 2.4)  return false;
  if ( fabs(sTree->lep2_.Eta() ) > 2.4)  return false;

  //consistency with pf leptons
  if ( fabs( sTree->pflep1_.Pt() - sTree->lep1_.Pt() ) > 10. )  return false;
  if ( fabs( sTree->pflep2_.Pt() - sTree->lep2_.Pt() ) > 10. )  return false;

  //information is only stored for leading lepton
  if ( ( sTree->isopf1_ * sTree->lep1_.Pt() ) > 5. )  return false; 
  if ( fabs(sTree->id1_)==11 && sTree->eoverpin_ > 4. ) return false;

  //barrel only electrons
  if (fabs(sTree->id1_)==11 && fabs(sTree->lep1_.Eta() ) > 1.4442) return false;
  if (fabs(sTree->id2_)==11 && fabs(sTree->lep2_.Eta() ) > 1.4442) return false;
  
  return true;

}

bool StopTreeLooper::passLepPlusIsoTrkSelection(const StopTree *sTree, bool isData) 
{
  //single lepton plus iso trk selection for 8 TeV 53 analysis

  //at least one lepton
  if ( !passSingleLeptonSelection(sTree, isData) ) return false;

  //pass isolated track requirement
  //unfortunately changed default value to 9999.
  if ( sTree->pfcandpt10_ > 9990. || sTree->pfcandiso10_ > 0.1 ) return false;

  return true;

}

pair<float,float> StopTreeLooper::getPhiCorrMET( float met, float metphi, int nvtx, bool ismc){

  //using met phi corrections from C. Veelken (emails from Oct. 4th)
  //previous versions are available here:
  //http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/CMSSW/JetMETCorrections/Type1MET/python/pfMETsysShiftCorrections_cfi.py

  // Data
  // ------
  // x :  "+2.87340e-01 + 3.29813e-01*Nvtx" 
  // y : "-2.27938e-01 - 1.71272e-01*Nvtx"
  // MC
  // ------            
  // x : "+8.72683e-02 - 1.66671e-02*Nvtx"
  // y :  "+1.86650e-01 - 1.21946e-01*Nvtx"
  

  float metx = met * cos( metphi );
  float mety = met * sin( metphi );

  float shiftx = 0.;
  float shifty = 0.;

  //use correction for data vs. mc 
  shiftx = ismc ? (+8.72683e-02 - 1.66671e-02*nvtx)
    : (+2.87340e-01 + 3.29813e-01*nvtx);
  shifty = ismc ? (+1.86650e-01 - 1.21946e-01*nvtx)
    : (-2.27938e-01 - 1.71272e-01*nvtx);
  
  metx -= shiftx;
  mety -= shifty;

  pair<float, float> phicorrmet = make_pair( sqrt( metx*metx + mety*mety ), atan2( mety , metx ) );
  return phicorrmet;
}

float StopTreeLooper::getdltrigweight(int id1, int id2)
{ 
  if (abs(id1)==11 && abs(id2)==11) return 0.95;
  if (abs(id1)==13 && abs(id2)==13) return 0.88;
  if (abs(id1)!=abs(id2)) return 0.92;
  return -999.;

}

float StopTreeLooper::getsltrigweight(int id1, float pt, float eta) {

  //electron efficiencies
  if ( abs(id1)==11 ) {
    if ( fabs(eta)<1.5) {
      if ( pt>=20 && pt<22 ) return 0.00;
      if ( pt>=22 && pt<24 ) return 0.00;
      if ( pt>=24 && pt<26 ) return 0.00;
      if ( pt>=26 && pt<28 ) return 0.08;
      if ( pt>=28 && pt<30 ) return 0.61;
      if ( pt>=30 && pt<32 ) return 0.86;
      if ( pt>=32 && pt<34 ) return 0.88;
      if ( pt>=34 && pt<36 ) return 0.90;
      if ( pt>=36 && pt<38 ) return 0.91;
      if ( pt>=38 && pt<40 ) return 0.92;
      if ( pt>=40 && pt<50 ) return 0.94;
      if ( pt>=50 && pt<60 ) return 0.95;
      if ( pt>=60 && pt<80 ) return 0.96;
      if ( pt>=80 && pt<100 ) return 0.96;
      if ( pt>=100 && pt<150 ) return 0.96;
      if ( pt>=150 && pt<200 ) return 0.97;
      if ( pt>=200 ) return 0.97;
    } else if ( fabs(eta)>=1.5 && fabs(eta)<2.1) {
      if ( pt>=20 && pt<22 ) return 0.00;
      if ( pt>=22 && pt<24 ) return 0.00;
      if ( pt>=24 && pt<26 ) return 0.02;
      if ( pt>=26 && pt<28 ) return 0.18;
      if ( pt>=28 && pt<30 ) return 0.50;
      if ( pt>=30 && pt<32 ) return 0.63;
      if ( pt>=32 && pt<34 ) return 0.68;
      if ( pt>=34 && pt<36 ) return 0.70;
      if ( pt>=36 && pt<38 ) return 0.72;
      if ( pt>=38 && pt<40 ) return 0.74;
      if ( pt>=40 && pt<50 ) return 0.76;
      if ( pt>=50 && pt<60 ) return 0.77;
      if ( pt>=60 && pt<80 ) return 0.78;
      if ( pt>=80 && pt<100 ) return 0.80;
      if ( pt>=100 && pt<150 ) return 0.79;
      if ( pt>=150 && pt<200 ) return 0.76;
      if ( pt>=200 ) return 0.81;
    }
  } else if ( abs(id1)==13 ) {//muon efficiencies

    if ( fabs(eta)<0.8 ) {
      if (pt>=20 && pt<22)  return  0.00;	 
      if (pt>=22 && pt<24)  return  0.03; 	 
      if (pt>=24 && pt<26)  return  0.87; 
      if (pt>=26 && pt<28)  return  0.90; 
      if (pt>=28 && pt<30)  return  0.91; 
      if (pt>=30 && pt<32)  return  0.91; 
      if (pt>=32 && pt<34)  return  0.92; 
      if (pt>=34 && pt<36)  return  0.93; 
      if (pt>=36 && pt<38)  return  0.93; 
      if (pt>=38 && pt<40)  return  0.93; 
      if (pt>=40 && pt<50)  return  0.94; 
      if (pt>=50 && pt<60)  return  0.95; 
      if (pt>=60 && pt<80)  return  0.95; 
      if (pt>=80 && pt<100) return 0.94; 
      if (pt>=100 && pt<150) return 0.94; 
      if (pt>=150 && pt<200) return 0.93; 
      if (pt>=200) return 0.92; 
    } else if ( fabs(eta)>=0.8 && fabs(eta)<1.5 ) {
      if (pt>=20 && pt<22)  return  0.00;
      if (pt>=22 && pt<24)  return  0.05;
      if (pt>=24 && pt<26)  return  0.78;
      if (pt>=26 && pt<28)  return  0.81;
      if (pt>=28 && pt<30)  return  0.81;
      if (pt>=30 && pt<32)  return  0.81;
      if (pt>=32 && pt<34)  return  0.82;
      if (pt>=34 && pt<36)  return  0.82;
      if (pt>=36 && pt<38)  return  0.83;
      if (pt>=38 && pt<40)  return  0.83;
      if (pt>=40 && pt<50)  return  0.84;
      if (pt>=50 && pt<60)  return  0.84;
      if (pt>=60 && pt<80)  return  0.84;
      if (pt>=80 && pt<100) return 0.84; 
      if (pt>=100 && pt<150) return 0.84;
      if (pt>=150 && pt<200) return 0.84;
      if (pt>=200) return 0.82;
    } else if ( fabs(eta)>=1.5 && fabs(eta)<2.1 ) {
      if (pt>=20 && pt<22)  return  0.00;
      if (pt>=22 && pt<24)  return  0.11;
      if (pt>=24 && pt<26)  return  0.76;
      if (pt>=26 && pt<28)  return  0.78;
      if (pt>=28 && pt<30)  return  0.79;
      if (pt>=30 && pt<32)  return  0.80;
      if (pt>=32 && pt<34)  return  0.80;
      if (pt>=34 && pt<36)  return  0.81;
      if (pt>=36 && pt<38)  return  0.81;
      if (pt>=38 && pt<40)  return  0.82;
      if (pt>=40 && pt<50)  return  0.82;
      if (pt>=50 && pt<60)  return  0.83;
      if (pt>=60 && pt<80)  return  0.83;
      if (pt>=80 && pt<100) return 0.83;
      if (pt>=100 && pt<150) return 0.83;
      if (pt>=150 && pt<200) return 0.82;
      if (pt>=200) return 0.82;
    }
  }//end check for muons

  return 1.;

}
