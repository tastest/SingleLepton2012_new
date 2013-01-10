//#include "../../CORE/jetSmearingTools.h"
#include "../../CORE/Thrust.h"
#include "../../CORE/EventShape.h"

#include "StopTreeLooper.h"

#include "../Plotting/PlotUtilities.h"
#include "../Core/MT2Utility.h"
#include "../Core/mt2bl_bisect.h"
#include "../Core/mt2w_bisect.h"
#include "../Core/PartonCombinatorics.h"
#include "../Core/stopUtils.h"

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


/* Reconstruct the hadronic top candidates, select the best candidate and
 * store the chi2 =  (m_jj - m_W)^2/sigma_m_jj + (m_jjj - m_t)^2/sigma_m_jjj
 * return the number of candidates found.
 *
 * n_jets - number of jets.
 * jets - jets
 * btag - b-tagging information of the jets
 * mc - qgjet montecarlo match number for the jets
 * 
 * returns a list of candidates sorted by chi2 ( if __sort = true in .h );
 */ 

list<Candidate> StopTreeLooper::recoHadronicTop(StopTree* tree, bool isData, bool wbtag){

  assert( tree->pfjets_->size() == tree->pfjets_csv_.size() );
  assert( tree->pfjets_->size() == tree->pfjets_sigma_.size() );

  LorentzVector* lep = &tree->lep1_;
  double met = tree->t1metphicorr_;
  double metphi = tree->t1metphicorrphi_;

  vector<LorentzVector> jets;
  vector<float>         btag;
  vector<int>           mc;

  //cout << endl << "baby branches:" << endl;
  for( unsigned int i = 0 ; i < tree->pfjets_->size() ; ++i ){
    //cout << i << " pt csv " << tree->pfjets_->at(i).pt() << " " << tree->pfjets_csv_.at(i) << endl;
    jets.push_back( tree->pfjets_->at(i)    );
    btag.push_back( tree->pfjets_csv_.at(i) );
    if ( !isData ) mc.push_back  ( tree->pfjets_mc3_.at(i) );
    else mc.push_back  ( 0 );
  } 

  // cout << endl << "stored:" << endl;
  // for( unsigned int i = 0 ; i < jets.size() ; ++i ){
  //   cout << i << " pt csv " << jets.at(i).pt() << " " << btag.at(i) << endl;
  // } 

  assert( jets.size() == btag.size() );

  float metx = met * cos( metphi );
  float mety = met * sin( metphi );

  int n_jets = jets.size();

  double sigma_jets[n_jets];
  for (int i=0; i<n_jets; ++i)
    sigma_jets[i] = tree->pfjets_sigma_.at(i);
  
  if ( isData )
    for (int i=0; i<n_jets; ++i)
      sigma_jets[i] *= getDataMCRatio(jets[i].eta());

  int ibl[5];
  int iw1[5];
  int iw2[5];
  int ib[5];

  if ( !isData ){
     
    // Matching MC algoritm search over all conbinations  until the 
    // right combination is found. More than one candidate is suported 
    //  but later only the first is used.
    // 
    int match = 0;
    for (int jbl=0; jbl<n_jets; ++jbl )
      for (int jb=0; jb<n_jets; ++jb )
        for (int jw1=0; jw1<n_jets; ++jw1 )
          for (int jw2=jw1+1; jw2<n_jets; ++jw2 )
            if ( (mc.at(jw2)==2 && mc.at(jw1)==2 && mc.at(jb)==1 && mc.at(jbl)==-1) ||
                 (mc.at(jw2)==-2 && mc.at(jw1)==-2 && mc.at(jb)==-1 && mc.at(jbl)==1) ) {
	      if ( match == 5 ) break; // this avoid the crash in events like ttW
	      ibl[match] = jbl;
	      iw1[match] = jw1;
	      iw2[match] = jw2;
	      ib[match] = jb;
	      match++;
            }
  }

  
  ////////    * Combinatorics. j_1 Pt must be > PTMIN_W1 and so on.
  
  vector<int> v_i, v_j;
  vector<double> v_k1, v_k2;
  for ( int i=0; i<n_jets; ++i )
    for ( int j=i+1; j<n_jets; ++j ){
      double pt_w1 = jets[i].Pt();
      double pt_w2 = jets[j].Pt();
      if (pt_w1 < PTMIN_J1  || pt_w2 < PTMIN_J2)
        continue;

      //
      //  W
      //
      LorentzVector hadW = jets[i] + jets[j];

      //
      //  W Mass Constraint.
      //
      TFitter *minimizer = new TFitter();
      double p1 = -1;

      minimizer->ExecuteCommand("SET PRINTOUT", &p1, 1);
      minimizer->SetFCN(PartonCombinatorics::minuitFunction);
      minimizer->SetParameter(0 , "c1"     , 1.1             , 1 , 0 , 0);
      minimizer->SetParameter(1 , "pt1"    , 1.0             , 1 , 0 , 0);
      minimizer->SetParameter(2 , "sigma1" , sigma_jets[i]   , 1 , 0 , 0);
      minimizer->SetParameter(3 , "pt2"    , 1.0             , 1 , 0 , 0);
      minimizer->SetParameter(4 , "sigma2" , sigma_jets[j]   , 1 , 0 , 0);
      minimizer->SetParameter(5 , "m12"    , jets[i].mass2() , 1 , 0 , 0);
      minimizer->SetParameter(6 , "m22"    , jets[j].mass2() , 1 , 0 , 0);
      minimizer->SetParameter(7 , "m02"    , hadW.mass2()    , 1 , 0 , 0);

      for (unsigned int k = 1; k < 8; k++)
        minimizer->FixParameter(k);

      minimizer->ExecuteCommand("SIMPLEX", 0, 0);
      minimizer->ExecuteCommand("MIGRAD", 0, 0);

      double c1 = minimizer->GetParameter(0);
      if (c1!=c1) {
        cout<<"[StopTreeLooper::recoHadronicTop] ERROR: c1 parameter is NAN! Skipping this parton combination"
	    <<endl;
        continue;
      }
      double c2 = PartonCombinatorics::fc2(c1, jets[i].mass2(), jets[j].mass2(), hadW.mass2(),false);
                
      delete minimizer;

     
      //     * W Mass check :)
      //     *  Never trust a computer you can't throw out a window. 
      //      *  - Steve Wozniak 

      //      cout << "c1 = " <<  c1 << "  c1 = " << c2 << "   M_jj = " 
      //           << ((jets[i] * c1) + (jets[j] * c2)).mass() << endl;
      
      v_i.push_back(i);
      v_j.push_back(j);
      v_k1.push_back(c1);
      v_k2.push_back(c2);
    }


  list<Candidate> chi2candidates;
        
  mt2_bisect::mt2 mt2_event;
  mt2bl_bisect::mt2bl mt2bl_event;
  mt2w_bisect::mt2w mt2w_event;
  
  for ( int b=0; b<n_jets; ++b )
    for (int o=0; o<n_jets; ++o){
      if ( b == o )
        continue;

      if ( wbtag && btag[b] < BTAG_MIN && btag[o] < BTAG_MIN )
        continue;

      double pt_b = jets[b].Pt();
      if ( wbtag && btag[b] >= BTAG_MIN && pt_b < PTMIN_BTAG )
        continue;

      if ( btag[b] < BTAG_MIN && pt_b < PTMIN_B )
        continue;

      double pt_o = jets[o].Pt();
      if ( wbtag && btag[o] >= BTAG_MIN && pt_o < PTMIN_OTAG )
        continue;

      if ( wbtag && btag[o] < BTAG_MIN && pt_o < PTMIN_O)
        continue;

      ///
      //  MT2 Variables
      ///

      double pl[4];     // Visible lepton
      double pb1[4];    // bottom on the same side as the visible lepton
      double pb2[4];    // other bottom, paired with the invisible W
      double pmiss[3];  // <unused>, pmx, pmy   missing pT
      pl[0]= lep->E(); pl[1]= lep->Px(); pl[2]= lep->Py(); pl[3]= lep->Pz();
      pb1[1] = jets[o].Px();  pb1[2] = jets[o].Py();   pb1[3] = jets[o].Pz();
      pb2[1] = jets[b].Px();  pb2[2] = jets[b].Py();   pb2[3] = jets[b].Pz();
      pmiss[0] = 0.; pmiss[1] = metx; pmiss[2] = mety;

      double pmiss_lep[3];
      pmiss_lep[0] = 0.;
      pmiss_lep[1] = pmiss[1]+pl[1]; pmiss_lep[2] = pmiss[2]+pl[2];

      pb1[0] = jets[o].mass();
      pb2[0] = jets[b].mass();
      mt2_event.set_momenta( pb1, pb2, pmiss_lep );
      mt2_event.set_mn( 80.385 );   // Invisible particle mass
      double c_mt2b = mt2_event.get_mt2();

      pb1[0] = jets[o].E();
      pb2[0] = jets[b].E();
      mt2bl_event.set_momenta(pl, pb1, pb2, pmiss);
      double c_mt2bl = mt2bl_event.get_mt2bl();

      mt2w_event.set_momenta(pl, pb1, pb2, pmiss);
      double c_mt2w = mt2w_event.get_mt2w();

      //      cout << b << ":"<< btag[b] << " - " << o << ":" << btag[o] << " = " << c_mt2w << endl;

      for (unsigned int w = 0; w < v_i.size() ; ++w ){
        int i = v_i[w];
        int j = v_j[w];
        if ( i==o || i==b || j==o || j==b )
	  continue;

        double pt_w1 = jets[i].Pt();
        double pt_w2 = jets[j].Pt();

	///
	//  W Mass.
	///
	LorentzVector hadW = jets[i] + jets[j];
	double massW = hadW.mass();

	double c1 = v_k1[w];
	double c2 = v_k2[w];

	///
	// Top Mass.
	///
        LorentzVector hadT = (jets[i] * c1) + (jets[j] * c2) + jets[b];
        double massT = hadT.mass();

        double pt_w = hadW.Pt();
        double sigma_w2 = pt_w1*sigma_jets[i] * pt_w1*sigma_jets[i]
	  + pt_w2*sigma_jets[j] * pt_w2*sigma_jets[j];
        double smw2 = (1.+2.*pt_w*pt_w/massW/massW)*sigma_w2;
        double pt_t = hadT.Pt();
        double sigma_t2 = c1*pt_w1*sigma_jets[i] * c1*pt_w1*sigma_jets[i]
	  + c2*pt_w2*sigma_jets[j] * c2*pt_w2*sigma_jets[j]
	  + pt_b*sigma_jets[b] * pt_b*sigma_jets[b];
        double smtop2 = (1.+2.*pt_t*pt_t/massT/massT)*sigma_t2;

        double c_chi2 = (massT-PDG_TOP_MASS)*(massT-PDG_TOP_MASS)/smtop2
	  + (massW-PDG_W_MASS)*(massW-PDG_W_MASS)/smw2;

        bool c_match = ( !isData &&  iw1[0]==i && iw2[0]==j && ib[0]==b && ibl[0]==o );
  
        Candidate c;
        c.chi2  = c_chi2;
        c.mt2b  = c_mt2b;
        c.mt2w  = c_mt2w;
        c.mt2bl = c_mt2bl;
        c.j1 = i;
        c.j2 = j;
        c.bi = b;
        c.oi = o;
        c.k1 = c1;
        c.k2 = c2;
        c.match = c_match;

        chi2candidates.push_back(c);
      }
    }

  if (__SORT) 
    chi2candidates.sort(PartonCombinatorics::compare_in_chi2);

  return chi2candidates;
}

//--------------------------------------------------------------------

// removes candidates without b-tagging requirement
list<Candidate> StopTreeLooper::getBTaggedCands(list<Candidate> &candidates, StopTree* tree) 
{

  list<Candidate> bcands; 

  list<Candidate>::iterator candIter;

  for(candIter = candidates.begin() ; candIter != candidates.end() ; candIter++ ){

    if( tree->pfjets_csv_.at(candIter->bi) < BTAG_MIN && tree->pfjets_csv_.at(candIter->oi) < BTAG_MIN ) continue;
    bcands.push_back(*candIter);

  }

  return bcands;

}

//--------------------------------------------------------------------

MT2struct StopTreeLooper::Best_MT2Calculator_Ricardo(list<Candidate> candidates, StopTree* tree, bool isData){

  if (candidates.size() == 0){
    MT2struct mfail;
    mfail.mt2w  = -0.999;
    mfail.mt2b  = -0.999;
    mfail.mt2bl = -0.999;
    mfail.chi2  = -0.999;
    return mfail;
  }

  double chi2_min  = 9999;
  double mt2b_min  = 9999;
  double mt2bl_min = 9999;
  double mt2w_min  = 9999;

  // count btags among leading 4 jets
  int n_btag = 0;

  for( int i = 0 ; i < 4 ; i++ ){
    if( tree->pfjets_csv_.at(i) < 0.679 ) n_btag++;
  }

  list<Candidate>::iterator candIter;

  for(candIter = candidates.begin() ; candIter != candidates.end() ; candIter++ ){

    // loop over all candidiates
    //for (unsigned int i=0; i<candidates.size(); ++i) {

    //if ( ! candidates_->at(i).match ) continue;

    // get indices of 4 jets used for hadronic top reconstruction
    // int b  = candidates.at(i).bi;
    // int o  = candidates.at(i).oi;
    // int j1 = candidates.at(i).j1;
    // int j2 = candidates.at(i).j2;

    int b  = (*candIter).bi;
    int o  = (*candIter).oi;
    int j1 = (*candIter).j1;
    int j2 = (*candIter).j2;

    // require the 4 jets used for the hadronic top mass are the 4 leading jets
    if (b>3 || o > 3 || j1 > 3 || j2 > 3) continue;

    // store whether the 2 jets used for MT2 calculation are btagged
    bool b_btag  = (tree->pfjets_csv_.at(b)  > 0.679);
    bool o_btag  = (tree->pfjets_csv_.at(o)  > 0.679);
    //bool j1_btag = (btag_->at(j1) > 0.679);
    //bool j2_btag = (btag_->at(j2) > 0.679);

    // 2 btags: require jets used for MT2 are the 2 b-jets
    if ( n_btag == 2 ){ 
      if ( !b_btag || !o_btag ) continue;
    } 

    // 1 btag: require the bjet and one of the 2 leading non bjets is used for MT2
    else if ( n_btag == 1) {

      if( b_btag ){
	if( b == 3  && o > 1 ) continue;
	if( b <   3 && o > 2 ) continue;
      }

      if( o_btag ){
	if( o == 3  && b > 1 ) continue;
	if( o <   3 && b > 2 ) continue;
      }

      //if (b>1 || o>1) continue;
    } 

    // 0 or >=3 btags: require jets used for MT2 are among 3 leading jets
    else {
      if (b>2 || o>2) continue;
    }

    // double chi2  = candidates.at(i).chi2;
    // double mt2b  = candidates.at(i).mt2b;
    // double mt2bl = candidates.at(i).mt2bl;
    // double mt2w  = candidates.at(i).mt2w;

    double chi2  = (*candIter).chi2;
    double mt2b  = (*candIter).mt2b;
    double mt2bl = (*candIter).mt2bl;
    double mt2w  = (*candIter).mt2w;
 
    //    cout << " " << b << ":" << b_btag  << " " << o << ":" << o_btag << " " << j1 << " " << j2 << " = " << mt2w  << "  " << chi2 <<endl;
    if (chi2  < chi2_min  ) chi2_min  = chi2;
    if (mt2b  < mt2b_min  ) mt2b_min  = mt2b;
    if (mt2bl < mt2bl_min ) mt2bl_min = mt2bl;
    if (mt2w  < mt2w_min  ) mt2w_min  = mt2w;
  }

  if( mt2w_min  > 9000 ) mt2w_min  = -0.999;
  if( mt2b_min  > 9000 ) mt2b_min  = -0.999;
  if( mt2bl_min > 9000 ) mt2bl_min = -0.999;
  if( chi2_min  > 9000 ) chi2_min  = -0.999;

  MT2struct m;
  m.mt2w  = mt2w_min;
  m.mt2b  = mt2b_min;
  m.mt2bl = mt2bl_min;
  m.chi2  = chi2_min;

  return m;

}

void StopTreeLooper::loop(TChain *chain, TString name)
{

  TRandom r;

  printf("[StopTreeLooper::loop] %s\n", name.Data());

  load_badlaserevents  ();

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

  makeTree(name.Data());

  TH2F* h_nsig = new TH2F();
  
  if( name.Contains("T2") ){
    char* h_nsig_filename = "";

    if( name.Contains("T2tt") ){
      h_nsig_filename = "/tas/dalfonso/cms2V05-03-18_stoplooperV00-02-07/crabT2tt_3/myMassDB_T2tt.root";
    }

    cout << "[StopTreeLooper::loop] opening mass TH2 file " << h_nsig_filename << endl;

    TFile *f_nsig = TFile::Open(h_nsig_filename);

    if( f_nsig == 0 ) cout << "ERROR! didn't find file " << h_nsig_filename << endl;

    h_nsig = (TH2F*) f_nsig->Get("masses");
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

    StopTree *tree = new StopTree();
    tree->LoadTree(currentFile->GetTitle());
    tree->InitTree();

    //---------------------------------
    // event loop
    //---------------------------------

    ULong64_t nEvents = tree->tree_->GetEntries();
    //nEvents = 1000;

    for(ULong64_t event = 0; event < nEvents; ++event) {
      tree->tree_->GetEntry(event);

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
        DorkyEventIdentifier id = {tree->run_,tree->event_, tree->lumi_ };
        if (is_duplicate(id) ){
          continue;
        }
	if (is_badLaserEvent(id) ){
	  //std::cout<<"Removed bad laser calibration event:" <<tree->run_<<"   "<<tree->event_<<"\n";
	  continue;
	}
      }

      //------------------------------------------ 
      // event weight
      //------------------------------------------ 

      float evtweight    = isData ? 1. : ( tree->weight_ * tree->nvtxweight_ * tree->mgcor_ );

      if( name.Contains("T2") ) {
	int bin = h_nsig->FindBin(tree->mg_,tree->ml_);
	float nevents = h_nsig->GetBinContent(bin);
	evtweight = tree->xsecsusy_ * 1000.0 / nevents; 

	//cout << "mg " << tree->mg_ << " ml " << tree->ml_ << " bin " << bin << " nevents " << nevents << " xsec " << tree->xsecsusy_ << " weight " << evtweight << endl;
      }


      float trigweight   = isData ? 1. : getsltrigweight(tree->id1_, tree->lep1_.Pt(), tree->lep1_.Eta());
      float trigweightdl = isData ? 1. : getdltrigweight(tree->id1_, tree->id2_);

      //------------------------------------------ 
      // selection criteria
      //------------------------------------------ 
      
      if ( !passEvtSelection(tree, name) ) continue; // >=1 lepton, rho cut, MET filters, veto 2 nearby leptons
      if ( tree->npfjets30_ < 4          ) continue; // >=4 jets
      if ( tree->t1metphicorr_ < 50.0    ) continue; // MET > 50 GeV

      //------------------------------------------ 
      // get list of jets
      //------------------------------------------ 

      n_jets = 0;
      jets.clear();
      btag.clear();
      mc.clear();

      vector<float> sigma_jets;

      for( unsigned int i = 0 ; i < tree->pfjets_->size() ; ++i ){

	if ( i > (N_JETS_TO_CONSIDER-1) ) break;
	if ( tree->pfjets_->at(i).Pt() < JET_PT ) continue;
	if ( fabs(tree->pfjets_->at(i).Eta()) > JET_ETA ) continue;

	n_jets++;
	jets.push_back( tree->pfjets_->at(i)    );
	btag.push_back( tree->pfjets_csv_.at(i) );

	if ( !isData ) mc.push_back  ( tree->pfjets_mc3_.at(i) );
	else mc.push_back  ( 0 );

	float sigma = tree->pfjets_sigma_.at(i) ;
	if ( isData) sigma *= getDataMCRatio(jets[i].eta());
	sigma_jets.push_back(sigma);

      }

      for (int i=0; i<n_jets; ++i){
      }



      //------------------------------------------ 
      // get list of candidates
      //------------------------------------------ 
      
      LorentzVector* lep = &tree->lep1_;
      float met = tree->t1metphicorr_;
      float metphi = tree->t1metphicorrphi_;

      // get list of candidates
      PartonCombinatorics pc (jets, btag, sigma_jets, mc, *lep, met, metphi, isData);
      MT2CHI2 mc = pc.getMt2Chi2();

      //------------------------------------------ 
      // event shapes
      //------------------------------------------ 

      std::vector<LorentzVector>* myPfJets = tree->pfjets_;
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


      LorentzVector lep1_p4( tree->lep1_.px() , tree->lep1_.py() , 0 , tree->lep1_.E() ); 

      LorentzVector met_p4( tree->t1metphicorr_ * cos( tree->t1metphicorr_ ) , 
			    tree->t1metphicorr_ * sin( tree->t1metphicorr_ ) , 
			    0,
			    tree->t1metphicorr_
			    );

      jetLeptonVector.push_back(lep1_p4);

      jetLeptonMetVector.push_back(lep1_p4);
      jetLeptonMetVector.push_back(met_p4);

      for ( unsigned int i=0; i<myPfJets->size() ; i++) {

	if( myPfJets->at(i).pt()<30 )          continue;
	if( fabs(myPfJets->at(i).eta())>3.0 )  continue;

	LorentzVector jet_p4( myPfJets->at(i).px() , myPfJets->at(i).py() , 0 , myPfJets->at(i).E() );
    
	// jetVector.push_back(myPfJets->at(i));
	// jetLeptonVector.push_back(myPfJets->at(i));
	// jetLeptonMetVector.push_back(myPfJets->at(i));

	jetVector.push_back(jet_p4);
	jetLeptonVector.push_back(jet_p4);
	jetLeptonMetVector.push_back(jet_p4);

	float dPhiL = getdphi(tree->lep1_.Phi(), myPfJets->at(i).phi() );
	float dPhiM = getdphi(tree->t1metphicorrphi_, myPfJets->at(i).phi() );
    
	if(dPhiL<(3.14/2))  HT_SSL=HT_SSL+myPfJets->at(i).pt();
	if(dPhiL>=(3.14/2)) HT_OSL=HT_OSL+myPfJets->at(i).pt();

	if(dPhiM<(3.14/2))  HT_SSM=HT_SSM+myPfJets->at(i).pt();
	if(dPhiM>=(3.14/2)) HT_OSM=HT_OSM+myPfJets->at(i).pt();

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
      // variables to add to baby
      //------------------------------------------ 
      
      initBaby(); // set all branches to -1

      list<Candidate> allcandidates = recoHadronicTop(tree, isData , false);
      MT2struct mr                  = Best_MT2Calculator_Ricardo(allcandidates, tree, isData);

      // which selections are passed
      sig_        = ( passOneLeptonSelection(tree, isData) && tree->nbtagscsvm_>=1 ) ? 1 : 0; // pass signal region preselection
      cr1_        = ( passOneLeptonSelection(tree, isData) && tree->nbtagscsvm_==0 ) ? 1 : 0; // pass CR1 (b-veto) control region preselection
      cr4_        = ( passDileptonSelection(tree, isData) && tree->nbtagscsvm_>=1) ? 1 : 0; // pass CR4 (dilepton) control region preselection
      cr5_        = ( passLepPlusIsoTrkSelection(tree, isData) && tree->nbtagscsvm_>=1) ? 1 : 0; // pass CR1 (lepton+isotrack) control region preselection

      // kinematic variables
      met_        = tree->t1metphicorr_;       // MET (type1, MET-phi corrections)
      mt_         = tree->t1metphicorrmt_;     // MT (type1, MET-phi corrections)

      // the following variables are obsolete, to be updated
      chi2_       = mr.chi2;                   
      mt2w_       = mr.mt2w;                   
      mt2b_       = mr.mt2b;
      mt2bl_      = mr.mt2bl;

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
      nb_         = tree->nbtagscsvm_;            // nbjets (pT > 30, CSVM)
      njets_      = tree->npfjets30_;             // njets (pT > 30, eta < 2.5)

      // lepton variables
      passisotrk_   = passIsoTrkVeto(tree) ? 1 : 0; // is there an isolated track? (pT>10 GeV, reliso<0.1)
      passisotrkv2_ = passIsoTrkVeto_v2(tree) ? 1 : 0; // is there an isolated track? (pT>10 GeV, reliso<0.1)
      nlep_       = tree->ngoodlep_;              // number of analysis selected leptons
 
      lep1pt_     = tree->lep1_.pt();             // 1st lepton pt
      lep1eta_    = fabs( tree->lep1_.eta() );    // 1st lepton eta

      if( nlep_ > 1 ){
	lep2pt_    = tree->lep1_.pt();            // 2nd lepton pt
	lep2eta_   = tree->lep1_.eta();           // 2nd lepton eta
	dilmass_   = tree->dilmass_;              // dilepton mass
      }
      
      // susy vars
      mstop_       = tree->mg_;                   // stop mass
      mlsp_        = tree->ml_;                   // LSP mass
      x_           = tree->x_;                    // chargino mass parameter x

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

      dphimj1_    = getdphi(tree->t1metphicorrphi_, jetVector.at(0).phi() );
      dphimj2_    = getdphi(tree->t1metphicorrphi_, jetVector.at(1).phi() );
      dphimjmin_  = TMath::Min( dphimj1_ , dphimj2_ );

      rand_       = r.Uniform(0.0,1.0);

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

void StopTreeLooper::makeTree(const char *prefix){

  TDirectory *rootdir = gDirectory->GetDirectory("Rint:");
  rootdir->cd();

  outFile_   = new TFile(Form("output/%s_mini.root",prefix), "RECREATE");
  outFile_->cd();

  outTree_ = new TTree("t","Tree");

  outTree_->Branch("lep1pt"           ,        &lep1pt_          ,         "lep1pt/F"		);
  outTree_->Branch("lep1eta"          ,        &lep1eta_         ,         "lep1eta/F"		);
  outTree_->Branch("sig"              ,        &sig_             ,         "sig/I"		);
  outTree_->Branch("cr1"              ,        &cr1_             ,         "cr1/I"		);
  outTree_->Branch("cr4"              ,        &cr4_             ,         "cr4/I"		);
  outTree_->Branch("cr5"              ,        &cr5_             ,         "cr5/I"		);
  outTree_->Branch("met"              ,        &met_             ,         "met/F"		);
  outTree_->Branch("mt"               ,        &mt_              ,         "mt/F"		);
  outTree_->Branch("chi2"             ,        &chi2_            ,         "chi2/F"		);
  outTree_->Branch("mt2w"             ,        &mt2w_            ,         "mt2w/F"		);
  outTree_->Branch("mt2b"             ,        &mt2b_            ,         "mt2b/F"		);
  outTree_->Branch("mt2bl"            ,        &mt2bl_           ,         "mt2bl/F"		);
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
  outTree_->Branch("rand"             ,        &rand_            ,         "rand/F"             );

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

  rand_       = -1.0;
}

