
#include "StopTreeLooper.h"
#include "Core/StopTree.h"
#include "Plotting/PlotUtilities.h"

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


float StopTreeLooper::getHadronicTopMass(const StopTree *tree){ 

  LorentzVector jets[6];
  jets[0] = tree->pfjet1_;
  jets[1] = tree->pfjet2_;
  jets[2] = tree->pfjet3_;
  jets[3] = tree->pfjet4_;
  jets[4] = tree->pfjet5_;
  jets[5] = tree->pfjet6_;

  //loop over all jet pairs and find the one with the smallest dR with inv mass > 60 GeV

  float mhadtop = 0.;
  int i_j1 = -9;
  int i_j2 = -9;
  LorentzVector hadW = LorentzVector();
  float mindRjj = 9999.;
  for ( int i_jet=0; i_jet<6; ++i_jet ) {
    if ( jets[i_jet].Pt()<20. ) continue;
    for ( int j_jet=i_jet+1; j_jet<6; ++j_jet ) {
      if ( jets[j_jet].Pt()<20. ) continue;
      float jjmass = (jets[i_jet]+jets[j_jet]).mass();
      if (jjmass<60.) continue;
      float dRjj = dRbetweenVectors(jets[i_jet],jets[j_jet]);
      if (dRjj<mindRjj) {
	hadW = jets[i_jet]+jets[j_jet];
	mindRjj = dRjj;
	i_j1 = i_jet;
	i_j2 = j_jet;
      }
    }
  }

  //look for third jet closest to hadronic W, with m(jjj) < 130 GeV

  if (mindRjj<9998.) {
    LorentzVector hadtop = LorentzVector();
    float mindRWj = 9999.;
    for ( int i_jet=0; i_jet<6; ++i_jet ) {
      if ( jets[i_jet].Pt()<20. ) continue;
      if (i_jet==i_j1 || i_jet==i_j2) continue;
      float Wjmass = (jets[i_jet]+hadW).mass();
      if (Wjmass<130.) continue;
      float dRWj = dRbetweenVectors(hadW, jets[i_jet]);
      if (dRWj<mindRWj) {
	hadtop = jets[i_jet]+hadW;
	mindRWj = dRWj;
	mhadtop = hadtop.mass();
      }
    }
  } 

  return mhadtop;

}


void makeYieldHistos( char* name , float mt , float met , float evtweight , std::map<string, TH1F*> & h_1d_sig ){

  if( mt > 150 && met > 100 )  plot1D(Form("%s_SRA",name),0.5, evtweight, h_1d_sig, 1, 0, 1);
  if( mt > 120 && met > 150 )  plot1D(Form("%s_SRB",name),0.5, evtweight, h_1d_sig, 1, 0, 1);
  if( mt > 120 && met > 200 )  plot1D(Form("%s_SRC",name),0.5, evtweight, h_1d_sig, 1, 0, 1);
  if( mt > 120 && met > 250 )  plot1D(Form("%s_SRD",name),0.5, evtweight, h_1d_sig, 1, 0, 1);
  if( mt > 120 && met > 300 )  plot1D(Form("%s_SRE",name),0.5, evtweight, h_1d_sig, 1, 0, 1);
  if( mt > 120 && met > 350 )  plot1D(Form("%s_SRF",name),0.5, evtweight, h_1d_sig, 1, 0, 1);
  if( mt > 120 && met > 400 )  plot1D(Form("%s_SRG",name),0.5, evtweight, h_1d_sig, 1, 0, 1);
  
}

void makeATLASYieldHistos( char* name , float mt , float met , float y, float evtweight , std::map<string, TH1F*> & h_1d_sig ){

  if( mt > 120 && met > 150 && y >  7.0 )  plot1D(Form("%s_ATLAS_SRA",name),0.5, evtweight, h_1d_sig, 1, 0, 1);
  if( mt > 120 && met > 150 && y >  9.0 )  plot1D(Form("%s_ATLAS_SRB",name),0.5, evtweight, h_1d_sig, 1, 0, 1);
  if( mt > 120 && met > 150 && y > 11.0 )  plot1D(Form("%s_ATLAS_SRC",name),0.5, evtweight, h_1d_sig, 1, 0, 1);
  if( mt > 130 && met > 225 && y > 11.0 )  plot1D(Form("%s_ATLAS_SRD",name),0.5, evtweight, h_1d_sig, 1, 0, 1);
  if( mt > 140 && met > 275 && y > 11.0 )  plot1D(Form("%s_ATLAS_SRE",name),0.5, evtweight, h_1d_sig, 1, 0, 1);
  
}

void makeATLASYieldHistos_NoMetSignificance( char* name , float mt , float met , float evtweight , std::map<string, TH1F*> & h_1d_sig ){

  if( mt > 120 && met > 150 )  plot1D(Form("%s_ATLAS_SRA",name),0.5, evtweight, h_1d_sig, 1, 0, 1);
  if( mt > 120 && met > 150 )  plot1D(Form("%s_ATLAS_SRB",name),0.5, evtweight, h_1d_sig, 1, 0, 1);
  if( mt > 120 && met > 150 )  plot1D(Form("%s_ATLAS_SRC",name),0.5, evtweight, h_1d_sig, 1, 0, 1);
  if( mt > 130 && met > 225 )  plot1D(Form("%s_ATLAS_SRD",name),0.5, evtweight, h_1d_sig, 1, 0, 1);
  if( mt > 140 && met > 275 )  plot1D(Form("%s_ATLAS_SRE",name),0.5, evtweight, h_1d_sig, 1, 0, 1);
  
}

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
  in.open("badlaser_events.txt");

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

  printf("[StopTreeLooper::loop] %s\n", name.Data());

  load_badlaserevents  ();

  //
  // check for valid chain
  //

  TObjArray *listOfFiles = chain->GetListOfFiles();
  TIter fileIter(listOfFiles);
  if (listOfFiles->GetEntries() == 0) {
    cout << "[StopTreeLooper::loop] no files in chain" << endl;
    return;
  }

  //
  // set up histograms
  //

  gROOT->cd();

  cout << "[StopTreeLooper::loop] setting up histos" << endl;

  //plotting map
  std::map<std::string, TH1F*> h_1d;
  //also for control regions
  std::map<std::string, TH1F*> h_1d_cr1, h_1d_cr2, h_1d_cr4, h_1d_cr5;
  //for signal region 
  std::map<std::string, TH1F*> h_1d_sig;
  //for ttbar dilepton njets distribution
  std::map<std::string, TH1F*> h_1d_nj;
  //z sample for yields etc
  std::map<std::string, TH1F*> h_1d_z;

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

  unsigned int nEventsChain=0;
  unsigned int nEvents = chain->GetEntries();
  nEventsChain = nEvents;
  unsigned int nEventsTotal = 0;
  int i_permille_old = 0;

  bool isData = name.Contains("data") ? true : false;

  //Define peak region 
  min_mtpeak = 50.; 
  max_mtpeak = 80.; 
  printf("[StopTreeLooper::loop] MT PEAK definition %.0f - %.0f GeV \n", min_mtpeak, max_mtpeak);

  cout << "[StopTreeLooper::loop] running over chain with total entries " << nEvents << endl;

  while (TChainElement *currentFile = (TChainElement*)fileIter.Next()) {

    //
    // load the stop baby tree
    //

    StopTree *tree = new StopTree();
    tree->LoadTree(currentFile->GetTitle());
    tree->InitTree();

    //
    // event loop
    //

    ULong64_t nEvents = tree->tree_->GetEntries();
    for(ULong64_t event = 0; event < nEvents; ++event) {
      tree->tree_->GetEntry(event);

      //
      // increment counters
      //

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

      // 
      // event weight
      // 

      float evtweight = isData ? 1. : ( tree->weight_ * 9.708 * tree->nvtxweight_ * tree->mgcor_ );

      if (name.Contains("T2")) evtweight = tree->xsecsusy_ * (1000./50000.) * 9.708;

      // to reweight from file - also need to comment stuff before
      //      float vtxweight = vtxweight_n( tree->nvtx_, h_vtx_wgt, isData );

      plot1D("h_vtx",       tree->nvtx_,       evtweight, h_1d_z, 40, 0, 40);
      plot1D("h_vtxweight", tree->nvtxweight_, evtweight, h_1d_z, 41, -4., 4.);

      // 
      // selection criteria
      // 

      //preselection
      if ( !passEvtSelection(tree, name) ) continue;

      // histogram tags
      //jet multiplicity
      string tag_njets = Form("_nj%i", (tree->npfjets30_<4) ? tree->npfjets30_ : 4);
      //b-tagging
      string tag_btag = (tree->nbtagscsvm_<1) ? "_bveto" : "";
      //iso-trk-veto
      string tag_isotrk = passIsoTrkVeto(tree) ? "" : "_wisotrk";
      //z-peak/veto
      string tag_zcut;
      if ( fabs( tree->dilmass_ - 91.) > 15. ) tag_zcut = "_zveto";
      else if  ( fabs( tree->dilmass_ - 91.) < 10. ) tag_zcut = "_zpeak";
      else tag_zcut = "_ignore";
      //Corrected jet counting with simple overlap removal
      //Check if each jet is matched to a lepton
      int njets_corr = tree->npfjets30_;
      if ( tree->mclep2_.Pt() > 30. ) {
	if ( tree->pfjet1_.Pt()>30. && dRbetweenVectors(tree->mclep2_, tree->pfjet1_) < 0.4 ) njets_corr--;
	if ( tree->pfjet2_.Pt()>30. && dRbetweenVectors(tree->mclep2_, tree->pfjet2_) < 0.4 ) njets_corr--;
	if ( tree->pfjet3_.Pt()>30. && dRbetweenVectors(tree->mclep2_, tree->pfjet3_) < 0.4 ) njets_corr--;
	if ( tree->pfjet4_.Pt()>30. && dRbetweenVectors(tree->mclep2_, tree->pfjet4_) < 0.4 ) njets_corr--;
      } 

     //to make plots for the two Kfactor bins
     //note this is to be used for the case where njets>=4
     string tag_kbin = (njets_corr<4) ? "_K3" : "_K4";

     //event with true truth-level track
     bool hastruetrk = false;
     if (tree->nleps_==2 && abs(tree->mclep2_.Eta())<2.5)  {
       //check if second lepton is e/mu pT>10GeV
       if (abs(tree->mcid2_)<14 && tree->mclep2_.Pt()>10.) hastruetrk = true;
       //if second lepton is tau 
       //check if daughter lepton or single track has pT>10GeV
       if (abs(tree->mcid2_)>14 && tree->mctaudpt2_>10.) {  
         if (tree->mcdecay2_==2) hastruetrk = true;
         if (tree->mcdecay2_==1 && tree->mcndec2_==1) hastruetrk = true;
       }
     }
     string tag_truetrk = hastruetrk ? "_wtruetrk" : "_notruetrk";

      //flavor types
      string flav_tag_sl;
      if ( abs(tree->id1_)==13 ) flav_tag_sl = "_muo";
      else if ( abs(tree->id1_)==11 ) flav_tag_sl = "_ele";
      else flav_tag_sl = "_mysterysl";
      string flav_tag_dl;
      if      ( abs(tree->id1_) == abs(tree->id2_) && abs(tree->id1_) == 13 ) flav_tag_dl = "_dimu";
      else if ( abs(tree->id1_) == abs(tree->id2_) && abs(tree->id1_) == 11 ) flav_tag_dl = "_diel";
      else if ( abs(tree->id1_) != abs(tree->id2_) && abs(tree->id1_) == 13 ) flav_tag_dl = "_muel";
      else if ( abs(tree->id1_) != abs(tree->id2_) && abs(tree->id1_) == 11 ) flav_tag_dl = "_elmu";
      else flav_tag_dl = "_mysterydl";
      string basic_flav_tag_dl = flav_tag_dl;
      if ( abs(tree->id1_) != abs(tree->id2_) ) basic_flav_tag_dl = "_mueg";

      //Apply the MET phi corrections to the type1 met ---- this will eventually be updated in the babies
      pair<float, float> p_t1metphicorr = 
	getPhiCorrMET( tree->t1met10_, tree->t1met10phi_, tree->nvtx_, !isData);
      t1metphicorr    = p_t1metphicorr.first;   // MET VARIABLE
      t1metphicorrphi = p_t1metphicorr.second;  // MET PHI VARIABLE
      t1metphicorrmt  = getMT( tree->lep1_.Pt() , tree->lep1_.Phi() , t1metphicorr , t1metphicorrphi ); // MT VARIABLE

      //
      // SIGNAL REGION - single lepton + b-tag
      //

      float met = t1metphicorr;
      float mt  = t1metphicorrmt;
      float ht  = tree->pfjet1_.pt() + tree->pfjet2_.pt() + tree->pfjet3_.pt() + tree->pfjet4_.pt();
      float y   = met/sqrt(ht);

      if( passATLASSelection(tree,isData,true ,true ,true)  ){
	makeATLASYieldHistos                  ("hyield"             ,mt,met,y,evtweight,h_1d_sig);
	makeATLASYieldHistos_NoMetSignificance("hyield_noMetSig"    ,mt,met,  evtweight,h_1d_sig);
      }
      if( passATLASSelection(tree,isData,false,true ,true)  )  makeATLASYieldHistos("hyield_nohadtop"    ,mt,met,y,evtweight,h_1d_sig);
      if( passATLASSelection(tree,isData,true ,false,true)  )  makeATLASYieldHistos("hyield_nodphi"      ,mt,met,y,evtweight,h_1d_sig);
      if( passATLASSelection(tree,isData,true ,true ,false) )  makeATLASYieldHistos("hyield_noisotrk"    ,mt,met,y,evtweight,h_1d_sig);

      // selection - 1 lepton 
      // Add iso track veto
      // Add b-tag
      if ( passSingleLeptonSelection(tree, isData) 
	   && tree->npfjets30_>=3 
	   && passIsoTrkVeto(tree) 
	   && tree->nbtagscsvm_>0 )
	{
	  //this is the event weight including the trigger SFs
	  evtweight *= isData ? 1. : getsltrigweight(tree->id1_, tree->lep1_.Pt(), tree->lep1_.Eta());

	  //met > 50 GeV requirement 
	  if ( t1metphicorr > 50. ) {

	    //this is where the code for selected events goes

	    //horrible code to attempt to reconstruct the hadronic top following the ATLAS procedure

	    LorentzVector jets[6];
	    jets[0] = tree->pfjet1_;
	    jets[1] = tree->pfjet2_;
	    jets[2] = tree->pfjet3_;
	    jets[3] = tree->pfjet4_;
	    jets[4] = tree->pfjet5_;
	    jets[5] = tree->pfjet6_;


	    //loop over all jet pairs and find the one with the smallest dR with inv mass > 60 GeV
	    /*
	    float mhadtop = 0.;
	    int i_j1 = -9;
	    int i_j2 = -9;
	    LorentzVector hadW = LorentzVector();
	    float mindRjj = 9999.;
	    for ( int i_jet=0; i_jet<6; ++i_jet ) {
	      if ( jets[i_jet].Pt()<20. ) continue;
	      for ( int j_jet=i_jet+1; j_jet<6; ++j_jet ) {
		if ( jets[j_jet].Pt()<20. ) continue;
		float jjmass = (jets[i_jet]+jets[j_jet]).mass();
		if (jjmass<60.) continue;
		float dRjj = dRbetweenVectors(jets[i_jet],jets[j_jet]);
		if (dRjj<mindRjj) {
		  hadW = jets[i_jet]+jets[j_jet];
		  mindRjj = dRjj;
		  i_j1 = i_jet;
		  i_j2 = j_jet;
		}
	      }
	    }

	    //look for third jet closest to hadronic W, with m(jjj) < 130 GeV

	    if (mindRjj<9998.) {
	      LorentzVector hadtop = LorentzVector();
	      float mindRWj = 9999.;
	      for ( int i_jet=0; i_jet<6; ++i_jet ) {
	    	if ( jets[i_jet].Pt()<20. ) continue;
	    	if (i_jet==i_j1 || i_jet==i_j2) continue;
	    	float Wjmass = (jets[i_jet]+hadW).mass();
	    	if (Wjmass<130.) continue;
	    	float dRWj = dRbetweenVectors(hadW, jets[i_jet]);
	    	if (dRWj<mindRWj) {
	    	  hadtop = jets[i_jet]+hadW;
	    	  mindRWj = dRWj;
	    	  mhadtop = hadtop.mass();
	    	}
	      }
	    } 
	    */

	    float mhadtop = getHadronicTopMass(tree);

	    // hadronic top mass = mhadtop

	    //store
	    plot1D("h_mhadtop",min(mhadtop, (float)499.99), evtweight, h_1d_sig, 50, 0, 500.);
	    
	    float met = t1metphicorr;
	    float mt  = t1metphicorrmt;
	    float ht  = tree->pfjet1_.pt() + tree->pfjet2_.pt() + tree->pfjet3_.pt() + tree->pfjet4_.pt();
	    float y   = met/sqrt(ht);

	    //cout << "met " << met << " mt " << mt << endl;

	    float dphi1 = getdphi( t1metphicorrphi , jets[0].phi() );
	    float dphi2 = getdphi( t1metphicorrphi , jets[1].phi() );
	    
	    if( tree->npfjets30_>=4 ){

	      makeYieldHistos("hyield",mt,met,evtweight,h_1d_sig);

	      makeATLASYieldHistos("hyieldCMS",mt,met,y,evtweight,h_1d_sig);
	      
	      if( dphi1 > 0.8 && dphi2 > 0.8 ){
		makeYieldHistos("hyield_dphi",mt,met,evtweight,h_1d_sig);
	      }
	      
	      if( mhadtop > 130.0 && mhadtop < 205.0 ){
		makeYieldHistos("hyield_hadtop",mt,met,evtweight,h_1d_sig);
	      }

	    }

	    if( tree->pfjet1_.pt() > 80.0 && tree->pfjet2_.pt() > 60.0 && tree->pfjet3_.pt() > 40.0 && tree->pfjet4_.pt() > 25.0 ){
	      makeYieldHistos("hyield_jets",mt,met,evtweight,h_1d_sig);
	    }

	  } // end met cut
	} // end signal sample selection

    } // end event loop

    // delete tree;

  } // end file loop

    //
    // finish
    //

  // TFile outfile(m_outfilename_.c_str(),"RECREATE") ; 
  // printf("[StopTreeLooper::loop] Saving histograms to %s\n", m_outfilename_.c_str());

  // std::map<std::string, TH1F*>::iterator it1d;
  // for(it1d=h_1d.begin(); it1d!=h_1d.end(); it1d++) {
  //   it1d->second->Write(); 
  //   delete it1d->second;
  // }

  // outfile.Write();
  // outfile.Close();
    
  //signal region
  //h_1d_sig

  TFile outfile_sig(Form("SIG%s",m_outfilename_.c_str()),"RECREATE") ; 
  printf("[StopTreeLooper::loop] Saving SIG histograms to %s\n", m_outfilename_.c_str());

  std::map<std::string, TH1F*>::iterator it1d_sig;
  for(it1d_sig=h_1d_sig.begin(); it1d_sig!=h_1d_sig.end(); it1d_sig++) {
    it1d_sig->second->Write(); 
    delete it1d_sig->second;
  }

  outfile_sig.Write();
  outfile_sig.Close();

  //control regions
  //h_1d_cr1, h_1d_cr2, h_1d_cr4, h_1d_cr5

  TFile outfile_cr1(Form("CR1%s",m_outfilename_.c_str()),"RECREATE") ; 
  printf("[StopTreeLooper::loop] Saving CR1 histograms to %s\n", m_outfilename_.c_str());

  std::map<std::string, TH1F*>::iterator it1d_cr1;
  for(it1d_cr1=h_1d_cr1.begin(); it1d_cr1!=h_1d_cr1.end(); it1d_cr1++) {
    it1d_cr1->second->Write(); 
    delete it1d_cr1->second;
  }

  outfile_cr1.Write();
  outfile_cr1.Close();

  TFile outfile_cr2(Form("CR2%s",m_outfilename_.c_str()),"RECREATE") ; 
  printf("[StopTreeLooper::loop] Saving CR2 histograms to %s\n", m_outfilename_.c_str());

  std::map<std::string, TH1F*>::iterator it1d_cr2;
  for(it1d_cr2=h_1d_cr2.begin(); it1d_cr2!=h_1d_cr2.end(); it1d_cr2++) {
    it1d_cr2->second->Write(); 
    delete it1d_cr2->second;
  }

  outfile_cr2.Write();
  outfile_cr2.Close();
    
  TFile outfile_cr4(Form("CR4%s",m_outfilename_.c_str()),"RECREATE") ; 
  printf("[StopTreeLooper::loop] Saving CR4 histograms to %s\n", m_outfilename_.c_str());

  std::map<std::string, TH1F*>::iterator it1d_cr4;
  for(it1d_cr4=h_1d_cr4.begin(); it1d_cr4!=h_1d_cr4.end(); it1d_cr4++) {
    it1d_cr4->second->Write(); 
    delete it1d_cr4->second;
  }

  outfile_cr4.Write();
  outfile_cr4.Close();

  TFile outfile_cr5(Form("CR5%s",m_outfilename_.c_str()),"RECREATE") ; 
  printf("[StopTreeLooper::loop] Saving CR5 histograms to %s\n", m_outfilename_.c_str());

  std::map<std::string, TH1F*>::iterator it1d_cr5;
  for(it1d_cr5=h_1d_cr5.begin(); it1d_cr5!=h_1d_cr5.end(); it1d_cr5++) {
    it1d_cr5->second->Write(); 
    delete it1d_cr5->second;
  }

  outfile_cr5.Write();
  outfile_cr5.Close();
    
  TFile outfile_nj(Form("NJ%s",m_outfilename_.c_str()),"RECREATE") ; 
  printf("[StopTreeLooper::loop] Saving NJ histograms to %s\n", m_outfilename_.c_str());

  std::map<std::string, TH1F*>::iterator it1d_nj;
  for(it1d_nj=h_1d_nj.begin(); it1d_nj!=h_1d_nj.end(); it1d_nj++) {
    it1d_nj->second->Write(); 
    delete it1d_nj->second;
  }

  outfile_nj.Write();
  outfile_nj.Close();

  TFile outfile_z(Form("Z%s",m_outfilename_.c_str()),"RECREATE") ; 
  printf("[StopTreeLooper::loop] Saving Z histograms to %s\n", m_outfilename_.c_str());

  std::map<std::string, TH1F*>::iterator it1d_z;
  for(it1d_z=h_1d_z.begin(); it1d_z!=h_1d_z.end(); it1d_z++) {
    it1d_z->second->Write(); 
    delete it1d_z->second;
  }

  outfile_z.Write();
  outfile_z.Close();

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


bool StopTreeLooper::passATLASSelection(const StopTree *sTree, bool isData, bool hadtopcut, bool dphicut, bool isotrkcut) 
{
  //single lepton selection for 8 TeV 53 analysis

  //exactly one lepton
  if ( sTree->ngoodlep_ != 1 ) return false;
  
  if( isotrkcut && !passIsoTrkVeto(sTree) ) return false;

  //e-channel
  if ( sTree->leptype_ == 0 ) {

    //pass trigger if data - single electron
    if ( isData && sTree->ele27wp80_ != 1 ) return false;
    
    // pt > 25 GeV, |eta| < 2.47
    if ( sTree->lep1_.Pt() < 25 )          return false;
    if ( fabs(sTree->lep1_.Eta() ) > 2.47) return false;
    

  }

  //mu-channel
  else if ( sTree->leptype_ == 1 ) {

    //pass trigger if data - single muon
    if ( isData && sTree->isomu24_ != 1 ) return false;

    // pt > 20 GeV, |eta| < 2.4
    if ( sTree->lep1_.Pt() < 20 )          return false;
    if ( fabs(sTree->lep1_.Eta() ) > 2.4 ) return false;

  }

  // hadronic top mass reconstruction
  float mhadtop = getHadronicTopMass(sTree);
  if( hadtopcut && (mhadtop < 130.0 || mhadtop > 205.0) ) return false;
  
  // jet pt's
  if( sTree->pfjet1_.pt() < 80.0 ) return false;
  if( sTree->pfjet2_.pt() < 60.0 ) return false;
  if( sTree->pfjet3_.pt() < 40.0 ) return false;
  if( sTree->pfjet4_.pt() < 25.0 ) return false;

  // >=1 b-tag
  if( sTree->nbtagscsvm_ == 0 ) return false;

  // dphi(jet,met)
  float dphi1 = getdphi( t1metphicorrphi , sTree->pfjet1_.phi() );
  float dphi2 = getdphi( t1metphicorrphi , sTree->pfjet2_.phi() );
  if( dphicut && (dphi1 < 0.8 || dphi2 < 0.8 ) ) return false;


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

void StopTreeLooper::makeCR2Plots( const StopTree *sTree, float evtweight, std::map<std::string, TH1F*> &h_1d, 
				   string tag_selection, string tag_njets, string tag_kbin, string flav_tag_dl, float mtcut ) 
{

  //find positive lepton - this is the one that is combined with the pseudomet to form the mT
  bool isfirstp = (sTree->id1_ > 0) ? true : false;
  
  //recalculate met
  float metx = t1metphicorr * cos( t1metphicorrphi );
  float mety = t1metphicorr * sin( t1metphicorrphi );
          
  //recalculate the MET with the positive lepton
  metx += isfirstp ? sTree->lep1_.px() : sTree->lep2_.px();
  mety += isfirstp ? sTree->lep1_.py() : sTree->lep2_.py();
  
  float t1met10_lep    = sqrt(metx*metx + mety*mety);
  float t1met10phi_lep = atan2( mety , metx );
  
  //recalculate the MT with the negative lepton
  float t1met10mt_lep = isfirstp ?
    getMT( sTree->lep2_.Pt() , sTree->lep2_.Phi() , t1met10_lep , t1met10phi_lep ) :
    getMT( sTree->lep1_.Pt() , sTree->lep1_.Phi() , t1met10_lep , t1met10phi_lep );

  //binning for mT plots
  int nbins = 30;
  float h_xmin = 0.;
  float h_xmax = 300.;
  float x_ovflw = h_xmax-0.001;
  
  float pseudomt_count = -1.;
  if ( t1met10mt_lep > min_mtpeak
       && t1met10mt_lep < max_mtpeak )    pseudomt_count = 0.5;
  else if ( t1met10mt_lep > mtcut ) pseudomt_count = 1.5;
  
  //default met
  plot1D("h_cr2_met"+tag_selection          +flav_tag_dl, min(t1metphicorr, x_ovflw), evtweight, h_1d, nbins, h_xmin, h_xmax);
  plot1D("h_cr2_met"+tag_selection+tag_njets+flav_tag_dl, min(t1metphicorr, x_ovflw), evtweight, h_1d, nbins, h_xmin, h_xmax);
  plot1D("h_cr2_met"+tag_selection+tag_njets+tag_kbin+flav_tag_dl, min(t1metphicorr, x_ovflw), evtweight, h_1d, nbins, h_xmin, h_xmax);
  //pseudo-met
  plot1D("h_cr2_pseudomet"+tag_selection	  +flav_tag_dl, min(t1met10_lep, x_ovflw), evtweight, h_1d, nbins, h_xmin, h_xmax);
  plot1D("h_cr2_pseudomet"+tag_selection+tag_njets+flav_tag_dl, min(t1met10_lep, x_ovflw), evtweight, h_1d, nbins, h_xmin, h_xmax);
  plot1D("h_cr2_pseudomet"+tag_selection+tag_njets+tag_kbin+flav_tag_dl, min(t1met10_lep, x_ovflw), evtweight, h_1d, nbins, h_xmin, h_xmax);
  //positive lepton pt - enters mT calculation
  float leppt = isfirstp ? sTree->lep1_.Pt() : sTree->lep2_.Pt();
  plot1D("h_cr2_leppt"+tag_selection          +flav_tag_dl, min(leppt, x_ovflw), evtweight, h_1d, nbins, h_xmin, h_xmax);
  plot1D("h_cr2_leppt"+tag_selection+tag_njets+flav_tag_dl, min(leppt, x_ovflw), evtweight, h_1d, nbins, h_xmin, h_xmax);
  plot1D("h_cr2_leppt"+tag_selection+tag_njets+tag_kbin+flav_tag_dl, min(leppt, x_ovflw), evtweight, h_1d, nbins, h_xmin, h_xmax);
  //angle between pos-lep and pseudopseudomet
  float dphi_pseudometlep = isfirstp ?
    getdphi( sTree->lep2_.Phi() , t1met10phi_lep ) :
    getdphi( sTree->lep1_.Phi() , t1met10phi_lep );
  plot1D("h_cr2_dphi_pseudometlep"+tag_selection          +flav_tag_dl, dphi_pseudometlep, evtweight, h_1d, 15, 0., TMath::Pi());
  plot1D("h_cr2_dphi_pseudometlep"+tag_selection+tag_njets+flav_tag_dl, dphi_pseudometlep, evtweight, h_1d, 15, 0., TMath::Pi());
  plot1D("h_cr2_dphi_pseudometlep"+tag_selection+tag_njets+tag_kbin+flav_tag_dl, dphi_pseudometlep, evtweight, h_1d, 15, 0., TMath::Pi());
  //pseudo-mt
  plot1D("h_cr2_pseudomt"      +tag_selection          +flav_tag_dl, min(t1met10mt_lep, x_ovflw), evtweight, h_1d, nbins, h_xmin, h_xmax);
  plot1D("h_cr2_pseudomt"      +tag_selection+tag_njets+flav_tag_dl, min(t1met10mt_lep, x_ovflw), evtweight, h_1d, nbins, h_xmin, h_xmax);
  plot1D("h_cr2_pseudomt"      +tag_selection+tag_njets+tag_kbin+flav_tag_dl, min(t1met10mt_lep, x_ovflw), evtweight, h_1d, nbins, h_xmin, h_xmax);
  plot1D("h_cr2_pseudomt_count"+tag_selection          +flav_tag_dl, pseudomt_count, evtweight, h_1d, 2, 0, 2);
  plot1D("h_cr2_pseudomt_count"+tag_selection+tag_njets+flav_tag_dl, pseudomt_count, evtweight, h_1d, 2, 0, 2);
  plot1D("h_cr2_pseudomt_count"+tag_selection+tag_njets+tag_kbin+flav_tag_dl, pseudomt_count, evtweight, h_1d, 2, 0, 2);

}

void StopTreeLooper::makeCR4Plots( const StopTree *sTree, float evtweight, std::map<std::string, TH1F*> &h_1d, 
				   string tag_selection, string tag_njets, string tag_kbin, string flav_tag_dl, float mtcut ) 
{
  int nbins = 50;
  float h_xmin = 0.;
  float h_xmax = 500.;
  float x_ovflw = h_xmax-0.001;
  
  float mt_count = -1.;
  if ( t1metphicorrmt > min_mtpeak 
       && t1metphicorrmt < max_mtpeak )    mt_count = 0.5;
  else if ( t1metphicorrmt > mtcut ) mt_count = 1.5;
  
  //default met
  plot1D("h_cr4_met"+tag_selection                   +flav_tag_dl, min(t1metphicorr, x_ovflw), evtweight, h_1d, nbins-5, 50., h_xmax);
  plot1D("h_cr4_met"+tag_selection+tag_njets         +flav_tag_dl, min(t1metphicorr, x_ovflw), evtweight, h_1d, nbins-5, 50., h_xmax);
  plot1D("h_cr4_met"+tag_selection+tag_njets+tag_kbin+flav_tag_dl, min(t1metphicorr, x_ovflw), evtweight, h_1d, nbins-5, 50., h_xmax);
  //leading lepton pt - enters mT calculation
  plot1D("h_cr4_leppt"+tag_selection                   +flav_tag_dl, min(sTree->lep1_.Pt(), x_ovflw), evtweight, h_1d, nbins, h_xmin, h_xmax);
  plot1D("h_cr4_leppt"+tag_selection+tag_njets         +flav_tag_dl, min(sTree->lep1_.Pt(), x_ovflw), evtweight, h_1d, nbins, h_xmin, h_xmax);
  plot1D("h_cr4_leppt"+tag_selection+tag_njets+tag_kbin+flav_tag_dl, min(sTree->lep1_.Pt(), x_ovflw), evtweight, h_1d, nbins, h_xmin, h_xmax);
  //leading lepton eta
  plot1D("h_cr4_lepeta"+tag_selection                   +flav_tag_dl, sTree->lep1_.Eta(), evtweight, h_1d, 21, -2.1, 2.1);
  plot1D("h_cr4_lepeta"+tag_selection+tag_njets         +flav_tag_dl, sTree->lep1_.Eta(), evtweight, h_1d, 21, -2.1, 2.1);
  plot1D("h_cr4_lepeta"+tag_selection+tag_njets+tag_kbin+flav_tag_dl, sTree->lep1_.Eta(), evtweight, h_1d, 21, -2.1, 2.1);
  //subleading lepton pt
  plot1D("h_cr4_subleadleppt"+tag_selection                   +flav_tag_dl, min(sTree->lep2_.Pt(), x_ovflw), evtweight, h_1d, nbins, h_xmin, h_xmax);
  plot1D("h_cr4_subleadleppt"+tag_selection+tag_njets         +flav_tag_dl, min(sTree->lep2_.Pt(), x_ovflw), evtweight, h_1d, nbins, h_xmin, h_xmax);
  plot1D("h_cr4_subleadleppt"+tag_selection+tag_njets+tag_kbin+flav_tag_dl, min(sTree->lep2_.Pt(), x_ovflw), evtweight, h_1d, nbins, h_xmin, h_xmax);
  //angle between lead-lep and met
  float dphi_metlep = getdphi( sTree->lep1_.Phi() , t1metphicorrphi );
  plot1D("h_cr4_dphi_metlep"+tag_selection                   +flav_tag_dl, dphi_metlep, evtweight, h_1d, 15, 0., TMath::Pi());
  plot1D("h_cr4_dphi_metlep"+tag_selection+tag_njets         +flav_tag_dl, dphi_metlep, evtweight, h_1d, 15, 0., TMath::Pi());
  plot1D("h_cr4_dphi_metlep"+tag_selection+tag_njets+tag_kbin+flav_tag_dl, dphi_metlep, evtweight, h_1d, 15, 0., TMath::Pi());
  //MT
  //binning for mT plots
  nbins = 30;
  h_xmin = 0.;
  h_xmax = 300.;
  x_ovflw = h_xmax-0.001;
  plot1D("h_cr4_mt"+tag_selection                   +flav_tag_dl, min(t1metphicorrmt, x_ovflw), evtweight, h_1d, nbins, h_xmin, h_xmax);
  plot1D("h_cr4_mt"+tag_selection+tag_njets         +flav_tag_dl, min(t1metphicorrmt, x_ovflw), evtweight, h_1d, nbins, h_xmin, h_xmax);
  plot1D("h_cr4_mt"+tag_selection+tag_njets+tag_kbin+flav_tag_dl, min(t1metphicorrmt, x_ovflw), evtweight, h_1d, nbins, h_xmin, h_xmax);
  plot1D("h_cr4_mt_count"+tag_selection                   +flav_tag_dl, mt_count, evtweight, h_1d, 2, 0, 2);
  plot1D("h_cr4_mt_count"+tag_selection+tag_njets         +flav_tag_dl, mt_count, evtweight, h_1d, 2, 0, 2);
  plot1D("h_cr4_mt_count"+tag_selection+tag_njets+tag_kbin+flav_tag_dl, mt_count, evtweight, h_1d, 2, 0, 2);

  //Plot more angles between various objects
  //angle between 2 leptons
  float dphi_dilep = getdphi( sTree->lep1_.Phi() ,  sTree->lep2_.Phi() );
  plot1D("h_cr4_dphi_dilep"+tag_selection                   +flav_tag_dl, dphi_dilep, evtweight, h_1d, 15, 0., TMath::Pi());
  plot1D("h_cr4_dphi_dilep"+tag_selection+tag_njets         +flav_tag_dl, dphi_dilep, evtweight, h_1d, 15, 0., TMath::Pi());
  plot1D("h_cr4_dphi_dilep"+tag_selection+tag_njets+tag_kbin+flav_tag_dl, dphi_dilep, evtweight, h_1d, 15, 0., TMath::Pi());
  //dR between 2 leptons
  float dR_dilep = dRbetweenVectors( sTree->lep1_ ,  sTree->lep2_ );
  plot1D("h_cr4_dR_dilep"+tag_selection                   +flav_tag_dl, min(dR_dilep, (float)4.999), evtweight, h_1d, 15, 0., 5.);
  plot1D("h_cr4_dR_dilep"+tag_selection+tag_njets         +flav_tag_dl, min(dR_dilep, (float)4.999), evtweight, h_1d, 15, 0., 5.);
  plot1D("h_cr4_dR_dilep"+tag_selection+tag_njets+tag_kbin+flav_tag_dl, min(dR_dilep, (float)4.999), evtweight, h_1d, 15, 0., 5.);

}

void StopTreeLooper::makeCR5Plots( const StopTree *sTree, float evtweight, std::map<std::string, TH1F*> &h_1d, 
				   string tag_selection, string tag_njets, string tag_kbin, string flav_tag, float mtcut ) 
{

  int nbins = 50;
  float h_xmin = 0.;
  float h_xmax = 500.;
  float x_ovflw = h_xmax-0.001;
  
  float mt_count = -1.;
  if ( t1metphicorrmt > min_mtpeak 
       && t1metphicorrmt < max_mtpeak )    mt_count = 0.5;
  else if ( t1metphicorrmt > mtcut ) mt_count = 1.5;
  
  //default met
  plot1D("h_cr5_met"+tag_selection                   +flav_tag, min(t1metphicorr, x_ovflw), evtweight, h_1d, nbins-5, 50., h_xmax);
  plot1D("h_cr5_met"+tag_selection+tag_njets         +flav_tag, min(t1metphicorr, x_ovflw), evtweight, h_1d, nbins-5, 50., h_xmax);
  plot1D("h_cr5_met"+tag_selection+tag_njets+tag_kbin+flav_tag, min(t1metphicorr, x_ovflw), evtweight, h_1d, nbins-5, 50., h_xmax);
  //leading lepton pt - enters mT calculation
  plot1D("h_cr5_leppt"+tag_selection                   +flav_tag, min(sTree->lep1_.Pt(), x_ovflw), evtweight, h_1d, nbins, h_xmin, h_xmax);
  plot1D("h_cr5_leppt"+tag_selection+tag_njets         +flav_tag, min(sTree->lep1_.Pt(), x_ovflw), evtweight, h_1d, nbins, h_xmin, h_xmax);
  plot1D("h_cr5_leppt"+tag_selection+tag_njets+tag_kbin+flav_tag, min(sTree->lep1_.Pt(), x_ovflw), evtweight, h_1d, nbins, h_xmin, h_xmax);
  //isolated track pt
  plot1D("h_cr5_isotrkpt"+tag_selection                   +flav_tag, min(sTree->pfcand10_.Pt(), x_ovflw), evtweight, h_1d, nbins, h_xmin, h_xmax);
  plot1D("h_cr5_isotrkpt"+tag_selection+tag_njets         +flav_tag, min(sTree->pfcand10_.Pt(), x_ovflw), evtweight, h_1d, nbins, h_xmin, h_xmax);
  plot1D("h_cr5_isotrkpt"+tag_selection+tag_njets+tag_kbin+flav_tag, min(sTree->pfcand10_.Pt(), x_ovflw), evtweight, h_1d, nbins, h_xmin, h_xmax);
  //angle between lead-lep and met
  float dphi_metlep = getdphi( sTree->lep1_.Phi() , t1metphicorrphi );
  plot1D("h_cr5_dphi_metlep"+tag_selection                   +flav_tag, dphi_metlep, evtweight, h_1d, 15, 0., TMath::Pi());
  plot1D("h_cr5_dphi_metlep"+tag_selection+tag_njets         +flav_tag, dphi_metlep, evtweight, h_1d, 15, 0., TMath::Pi());
  plot1D("h_cr5_dphi_metlep"+tag_selection+tag_njets+tag_kbin+flav_tag, dphi_metlep, evtweight, h_1d, 15, 0., TMath::Pi());
  //MT
  //binning for mT plots
  nbins = 30;
  h_xmin = 0.;
  h_xmax = 300.;
  x_ovflw = h_xmax-0.001;
  plot1D("h_cr5_mt"+tag_selection                   +flav_tag, min(t1metphicorrmt, x_ovflw), evtweight, h_1d, nbins, h_xmin, h_xmax);
  plot1D("h_cr5_mt"+tag_selection+tag_njets         +flav_tag, min(t1metphicorrmt, x_ovflw), evtweight, h_1d, nbins, h_xmin, h_xmax);
  plot1D("h_cr5_mt"+tag_selection+tag_njets+tag_kbin+flav_tag, min(t1metphicorrmt, x_ovflw), evtweight, h_1d, nbins, h_xmin, h_xmax);
  plot1D("h_cr5_mt_count"+tag_selection                   +flav_tag, mt_count, evtweight, h_1d, 2, 0, 2);
  plot1D("h_cr5_mt_count"+tag_selection+tag_njets         +flav_tag, mt_count, evtweight, h_1d, 2, 0, 2);
  plot1D("h_cr5_mt_count"+tag_selection+tag_njets+tag_kbin+flav_tag, mt_count, evtweight, h_1d, 2, 0, 2);

  //Plot more angles between various objects
  //angle between lepton and isolated track
  float dphi_leptrk = getdphi( sTree->lep1_.Phi() ,  sTree->pfcand10_.Phi() );
  plot1D("h_cr5_dphi_leptrk"+tag_selection                   +flav_tag, dphi_leptrk, evtweight, h_1d, 15, 0., TMath::Pi());
  plot1D("h_cr5_dphi_leptrk"+tag_selection+tag_njets         +flav_tag, dphi_leptrk, evtweight, h_1d, 15, 0., TMath::Pi());
  plot1D("h_cr5_dphi_leptrk"+tag_selection+tag_njets+tag_kbin+flav_tag, dphi_leptrk, evtweight, h_1d, 15, 0., TMath::Pi());
  //dR between lepton and isolated track
  float dR_leptrk = dRbetweenVectors( sTree->lep1_ ,  sTree->pfcand10_ );
  plot1D("h_cr5_dR_leptrk"+tag_selection                   +flav_tag, min(dR_leptrk, (float)4.999), evtweight, h_1d, 15, 0., 5.);
  plot1D("h_cr5_dR_leptrk"+tag_selection+tag_njets         +flav_tag, min(dR_leptrk, (float)4.999), evtweight, h_1d, 15, 0., 5.);
  plot1D("h_cr5_dR_leptrk"+tag_selection+tag_njets+tag_kbin+flav_tag, min(dR_leptrk, (float)4.999), evtweight, h_1d, 15, 0., 5.);

}

void StopTreeLooper::makeSIGPlots( const StopTree *sTree, float evtweight, std::map<std::string, TH1F*> &h_1d, 
				   string tag_selection, string tag_kbin, string flav_tag, float mtcut ) 
{

  int nbins = 50;
  float h_xmin = 0.;
  float h_xmax = 500.;
  float x_ovflw = h_xmax-0.001;
  
  float mt_count = -1.;
  if ( t1metphicorrmt > min_mtpeak 
       && t1metphicorrmt < max_mtpeak )    mt_count = 0.5;
  else if ( t1metphicorrmt > mtcut ) mt_count = 1.5;
  
  //default met
  plot1D("h_sig_met"+tag_selection         +flav_tag, min(t1metphicorr, x_ovflw), evtweight, h_1d, nbins-5, 50, h_xmax);
  plot1D("h_sig_met"+tag_selection+tag_kbin+flav_tag, min(t1metphicorr, x_ovflw), evtweight, h_1d, nbins-5, 50, h_xmax);
  //lepton pt - enters mT calculation
  plot1D("h_sig_leppt"+tag_selection         +flav_tag, min(sTree->lep1_.Pt(), x_ovflw), evtweight, h_1d, nbins, h_xmin, h_xmax);
  plot1D("h_sig_leppt"+tag_selection+tag_kbin+flav_tag, min(sTree->lep1_.Pt(), x_ovflw), evtweight, h_1d, nbins, h_xmin, h_xmax);
  //angle between lepton and met
  float dphi_metlep = getdphi( sTree->lep1_.Phi() , t1metphicorrphi );
  plot1D("h_sig_dphi_metlep"+tag_selection         +flav_tag, dphi_metlep, evtweight, h_1d, 15, 0., TMath::Pi());
  plot1D("h_sig_dphi_metlep"+tag_selection+tag_kbin+flav_tag, dphi_metlep, evtweight, h_1d, 15, 0., TMath::Pi());
  //MT
  //binning for mT plots
  nbins = 30;
  h_xmin = 0.;
  h_xmax = 300.;
  x_ovflw = h_xmax-0.001;
  plot1D("h_sig_mt"+tag_selection         +flav_tag, min(t1metphicorrmt, x_ovflw), evtweight, h_1d, nbins, h_xmin, h_xmax);
  plot1D("h_sig_mt"+tag_selection+tag_kbin+flav_tag, min(t1metphicorrmt, x_ovflw), evtweight, h_1d, nbins, h_xmin, h_xmax);
  plot1D("h_sig_mt_count"+tag_selection         +flav_tag, mt_count, evtweight, h_1d, 2, 0, 2);
  plot1D("h_sig_mt_count"+tag_selection+tag_kbin+flav_tag, mt_count, evtweight, h_1d, 2, 0, 2);

}

void StopTreeLooper::makeCR1Plots( const StopTree *sTree, float evtweight, std::map<std::string, TH1F*> &h_1d, 
				   string tag_selection, string tag_njets, string tag_kbin, string flav_tag, float mtcut ) 
{

  int nbins = 50;
  float h_xmin = 0.;
  float h_xmax = 500.;
  float x_ovflw = h_xmax-0.001;

  float mt_count = -1.;
  if ( t1metphicorrmt > min_mtpeak 
       && t1metphicorrmt < max_mtpeak )    mt_count = 0.5;
  else if ( t1metphicorrmt > mtcut ) mt_count = 1.5;
  
  plot1D("h_cr1_njets"    +tag_selection+flav_tag, min(sTree->npfjets30_,4),  evtweight, h_1d, 5,0,5);
  plot1D("h_cr1_njets_all"+tag_selection+flav_tag, min(sTree->npfjets30_,9),  evtweight, h_1d, 10, 0, 10);
  //default met
  plot1D("h_cr1_met"+tag_selection                   +flav_tag, min(t1metphicorr, x_ovflw), evtweight, h_1d, nbins-5, 50., h_xmax);
  plot1D("h_cr1_met"+tag_selection+tag_njets         +flav_tag, min(t1metphicorr, x_ovflw), evtweight, h_1d, nbins-5, 50., h_xmax);
  plot1D("h_cr1_met"+tag_selection+tag_njets+tag_kbin+flav_tag, min(t1metphicorr, x_ovflw), evtweight, h_1d, nbins-5, 50., h_xmax);
  //lepton pt - enters mT calculation
  plot1D("h_cr1_leppt"+tag_selection                   +flav_tag, min(sTree->lep1_.Pt(), x_ovflw), evtweight, h_1d, nbins, h_xmin, h_xmax);
  plot1D("h_cr1_leppt"+tag_selection+tag_njets         +flav_tag, min(sTree->lep1_.Pt(), x_ovflw), evtweight, h_1d, nbins, h_xmin, h_xmax);
  plot1D("h_cr1_leppt"+tag_selection+tag_njets+tag_kbin+flav_tag, min(sTree->lep1_.Pt(), x_ovflw), evtweight, h_1d, nbins, h_xmin, h_xmax);
  //lepton phi
  plot1D("h_cr1_lepphi"+tag_selection                   +flav_tag, sTree->lep1_.Phi(), evtweight, h_1d, 30, -1.*TMath::Pi(), TMath::Pi());
  plot1D("h_cr1_lepphi"+tag_selection+tag_njets         +flav_tag, sTree->lep1_.Phi(), evtweight, h_1d, 30, -1.*TMath::Pi(), TMath::Pi());
  plot1D("h_cr1_lepphi"+tag_selection+tag_njets+tag_kbin+flav_tag, sTree->lep1_.Phi(), evtweight, h_1d, 30, -1.*TMath::Pi(), TMath::Pi());
  //met phi
  plot1D("h_cr1_metphi"+tag_selection                   +flav_tag, t1metphicorrphi, evtweight, h_1d, 30, -1.*TMath::Pi(), TMath::Pi());
  plot1D("h_cr1_metphi"+tag_selection+tag_njets         +flav_tag, t1metphicorrphi, evtweight, h_1d, 30, -1.*TMath::Pi(), TMath::Pi());
  plot1D("h_cr1_metphi"+tag_selection+tag_njets+tag_kbin+flav_tag, t1metphicorrphi, evtweight, h_1d, 30, -1.*TMath::Pi(), TMath::Pi());
  //angle between lepton and met
  float dphi_metlep = getdphi( sTree->lep1_.Phi() , t1metphicorrphi );
  plot1D("h_cr1_dphi_metlep"+tag_selection                   +flav_tag, dphi_metlep, evtweight, h_1d, 15, 0., TMath::Pi());
  plot1D("h_cr1_dphi_metlep"+tag_selection+tag_njets         +flav_tag, dphi_metlep, evtweight, h_1d, 15, 0., TMath::Pi());
  plot1D("h_cr1_dphi_metlep"+tag_selection+tag_njets+tag_kbin+flav_tag, dphi_metlep, evtweight, h_1d, 15, 0., TMath::Pi());
  //MT
  //binning for mT plots
  nbins = 30;
  h_xmin = 0.;
  h_xmax = 300.;
  x_ovflw = h_xmax-0.001;
  plot1D("h_cr1_mt"+tag_selection                   +flav_tag, min(t1metphicorrmt, x_ovflw), evtweight, h_1d, nbins, h_xmin, h_xmax);
  plot1D("h_cr1_mt"+tag_selection+tag_njets         +flav_tag, min(t1metphicorrmt, x_ovflw), evtweight, h_1d, nbins, h_xmin, h_xmax);
  plot1D("h_cr1_mt"+tag_selection+tag_njets+tag_kbin+flav_tag, min(t1metphicorrmt, x_ovflw), evtweight, h_1d, nbins, h_xmin, h_xmax);
  plot1D("h_cr1_mt_count"+tag_selection                   +flav_tag, mt_count, evtweight, h_1d, 2, 0, 2);
  plot1D("h_cr1_mt_count"+tag_selection+tag_njets         +flav_tag, mt_count, evtweight, h_1d, 2, 0, 2);
  plot1D("h_cr1_mt_count"+tag_selection+tag_njets+tag_kbin+flav_tag, mt_count, evtweight, h_1d, 2, 0, 2);

}

void StopTreeLooper::makeNJPlots( const StopTree *sTree, float evtweight, std::map<std::string, TH1F*> &h_1d, 
				  string tag_selection, string flav_tag ) 
{

  plot1D("h_njets"    +tag_selection,          min(sTree->npfjets30_,4), evtweight, h_1d, 4,1,5);
  plot1D("h_njets"    +tag_selection+flav_tag, min(sTree->npfjets30_,4), evtweight, h_1d, 4,1,5);
  plot1D("h_njets_all"+tag_selection,          min(sTree->npfjets30_,8), evtweight, h_1d, 7,1,8);
  plot1D("h_njets_all"+tag_selection+flav_tag, min(sTree->npfjets30_,8), evtweight, h_1d, 7,1,8);

}

void StopTreeLooper::makeZPlots( const StopTree *sTree, float evtweight, std::map<std::string, TH1F*> &h_1d, 
				  string tag_selection, string tag_njets, string flav_tag ) 
{

  int nbins = 30;
  float h_xmin = 0.;
  float h_xmax = 300.;
  float x_ovflw = h_xmax-0.001;

  plot1D("h_z_met"+tag_selection+flav_tag, min(t1metphicorr,x_ovflw), evtweight, h_1d, nbins, h_xmin, h_xmax);
  plot1D("h_z_met"+tag_selection+tag_njets+flav_tag, min(t1metphicorr,x_ovflw), evtweight, h_1d, nbins, h_xmin, h_xmax);

  string lep1type =  abs(sTree->id1_)==13 ? "h_muo" : "h_ele";
  string lep2type =  abs(sTree->id2_)==13 ? "h_muo" : "h_ele";
  plot1D(lep1type+"pt"+tag_selection+flav_tag, min(sTree->lep1_.Pt(),(float)199.99), evtweight, h_1d, 40, 20, 200);
  plot1D(lep2type+"pt"+tag_selection+flav_tag, min(sTree->lep2_.Pt(),(float)199.99), evtweight, h_1d, 40, 20, 200);
  plot1D(lep1type+"pt"+tag_selection+tag_njets+flav_tag, min(sTree->lep1_.Pt(),(float)199.99), evtweight, h_1d, 40, 20, 200);
  plot1D(lep2type+"pt"+tag_selection+tag_njets+flav_tag, min(sTree->lep2_.Pt(),(float)199.99), evtweight, h_1d, 40, 20, 200);
  
  plot1D("h_z_leppt"  +tag_selection+flav_tag, min(sTree->lep1_.Pt(),(float)299.99), evtweight, h_1d, 50, 20., 300.);
  plot1D("h_z_lepeta" +tag_selection+flav_tag, sTree->lep1_.Eta(), evtweight, h_1d, 24, -2.4, 2.4);
  plot1D("h_z_lep2pt" +tag_selection+flav_tag, min(sTree->lep2_.Pt(),(float)199.99), evtweight, h_1d, 50, 20., 200.);
  plot1D("h_z_lep2eta"+tag_selection+flav_tag, sTree->lep2_.Eta(), evtweight, h_1d, 24, -2.4, 2.4);
  plot1D("h_z_leppt"  +tag_selection+tag_njets+flav_tag, min(sTree->lep1_.Pt(),(float)299.99), evtweight, h_1d, 50, 20., 300.);
  plot1D("h_z_lepeta" +tag_selection+tag_njets+flav_tag, sTree->lep1_.Eta(), evtweight, h_1d, 24, -2.4, 2.4);
  plot1D("h_z_lep2pt" +tag_selection+tag_njets+flav_tag, min(sTree->lep2_.Pt(),(float)199.99), evtweight, h_1d, 50, 20., 200.);
  plot1D("h_z_lep2eta"+tag_selection+tag_njets+flav_tag, sTree->lep2_.Eta(), evtweight, h_1d, 24, -2.4, 2.4);

  float dphi_metlep = getdphi(sTree->lep1_.Phi(), t1metphicorrphi);
  plot1D("h_z_dphi_metl"+tag_selection+flav_tag, dphi_metlep, evtweight, h_1d, 15, 0., 3.14159);
  plot1D("h_z_dphi_metl"+tag_selection+tag_njets+flav_tag, dphi_metlep, evtweight, h_1d, 15, 0., 3.14159);
  plot1D("h_z_mt"+tag_selection+flav_tag, min(t1metphicorrmt, x_ovflw), evtweight, h_1d, nbins, h_xmin, h_xmax);
  plot1D("h_z_mt"+tag_selection+tag_njets+flav_tag, min(t1metphicorrmt, x_ovflw), evtweight, h_1d, nbins, h_xmin, h_xmax);

  if ( sTree->npfjets30_<1 ) return;
  plot1D("h_z_j1pt" +tag_selection+flav_tag, min(sTree->pfjet1_.Pt(), (float)399.99), evtweight, h_1d, 20, 30., 400.);
  plot1D("h_z_j1eta"+tag_selection+flav_tag, sTree->pfjet1_.Eta(), evtweight, h_1d, 24, -2.4, 2.4);
  if ( sTree->npfjets30_<2 ) return;
  plot1D("h_z_j2pt" +tag_selection+flav_tag, min(sTree->pfjet2_.Pt(), (float)299.99), evtweight, h_1d, 20, 30., 300.);
  plot1D("h_z_j2eta"+tag_selection+flav_tag, sTree->pfjet2_.Eta(), evtweight, h_1d, 24, -2.4, 2.4);
  if ( sTree->npfjets30_<3 ) return;
  plot1D("h_z_j3pt" +tag_selection+flav_tag, min(sTree->pfjet3_.Pt(), (float)199.99), evtweight, h_1d, 20, 30., 200.);
  plot1D("h_z_j3eta"+tag_selection+flav_tag, sTree->pfjet3_.Eta(), evtweight, h_1d, 24, -2.4, 2.4);
  if ( sTree->npfjets30_<4 ) return;
  plot1D("h_z_j4pt" +tag_selection+flav_tag, min(sTree->pfjet4_.Pt(), (float)119.99), evtweight, h_1d, 20, 30., 120.);
  plot1D("h_z_j4eta"+tag_selection+flav_tag, sTree->pfjet4_.Eta(), evtweight, h_1d, 24, -2.4, 2.4);
  
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
