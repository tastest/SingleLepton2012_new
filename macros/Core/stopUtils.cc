#include "stopUtils.h"

#include "Math/LorentzVector.h"
#include "Math/Point3D.h"
#include "TMath.h"
#include "TBranch.h"
#include "TTree.h"
#include "TH1F.h"
#include "TFile.h"
#include "TBits.h"

#include <vector> 
#include <list> 

#include "MT2Utility.h"
#include "mt2bl_bisect.h"
#include "mt2w_bisect.h"

//------------------------------------------------------------------------------------------------

list<Candidate> MT2Calculator(StopTree* tree, bool isData){

  //-----------------------------------------------------------------------------------------
  // This function calculate the MT2 variables without doing the hadronic top reconstruction
  // Function returns a list of candidates with mt2b, mt2bl, mt2w and the 2 b-jet indices
  // All values relevant only for the hadronic top reconstruction are set to -1
  //
  // EXAMPLE USAGE:
  // list<Candidate> mt2candidates = MT2Calculator(tree, isData );
  // list<Candidate>::iterator candIter;
  // for(candIter = mt2candidates.begin() ; candIter != mt2candidates.end() ; candIter++ ){
  // 	cout << "mt2w        " << (*candIter).mt2w   << endl;
  // 	cout << "mt2bl       " << (*candIter).mt2bl  << endl;
  // 	cout << "mt2b        " << (*candIter).mt2b   << endl;
  // 	cout << endl;
  // }

  //-----------------------------------------------------------------------------------------

  //-----------------------------------------
  // set jet pt and btagging thresholds
  //-----------------------------------------

  const double PTMIN_BTAG = 30;
  const double PTMIN_OTAG = 30; 
  const double PTMIN_B    = 30;  //  These two should be tigther than the  
  const double PTMIN_O    = 30;  //  b-tagged versions.
  const double BTAG_MIN = 0.679;

  //-----------------------------------------
  // get the lepton pt and MET
  //-----------------------------------------

  StopTree::LorentzVector* lep = &tree->lep1_;
  double met         = tree->t1metphicorr_;
  double metphi      = tree->t1metphicorrphi_;

  float metx = met * cos( metphi );
  float mety = met * sin( metphi );

  //-----------------------------------------
  // get the jets and b-tagging info
  //-----------------------------------------

  assert( tree->pfjets_->size() == tree->pfjets_csv_.size() );

  vector<StopTree::LorentzVector> jets;
  vector<float> btag;

  for( unsigned int i = 0 ; i < tree->pfjets_->size() ; ++i ){
    jets.push_back( tree->pfjets_->at(i) );
    btag.push_back( tree->pfjets_csv_.at(i) );
  } 

  int n_jets = jets.size();
  
  assert( jets.size() == btag.size() );

  //-----------------------------------------
  // initialize stuff for
  //-----------------------------------------

  list<Candidate> mt2candidates;
        
  mt2_bisect::mt2     mt2_event;
  mt2bl_bisect::mt2bl mt2bl_event;
  mt2w_bisect::mt2w   mt2w_event;

  //-----------------------------------------
  // loop over all pairs of jets
  //-----------------------------------------
  
  for ( int b=0; b<n_jets; ++b ){
    for (int o=0; o<n_jets; ++o){

      //---------------------------------------------------------------------------
      // select jet pairs with thresholds on jet pt's and b-tagging discriminants
      //---------------------------------------------------------------------------
      
      if ( b == o )
        continue;

      if ( btag[b] < BTAG_MIN && btag[o] < BTAG_MIN )
        continue;

      double pt_b = jets[b].Pt();

      if ( btag[b] >= BTAG_MIN && pt_b < PTMIN_BTAG )
        continue;

      if ( btag[b] < BTAG_MIN && pt_b < PTMIN_B )
        continue;

      double pt_o = jets[o].Pt();

      if ( btag[o] >= BTAG_MIN && pt_o < PTMIN_OTAG )
        continue;

      if ( btag[o] < BTAG_MIN && pt_o < PTMIN_O)
        continue;

      //-------------------------------------------------------------------------------
      // for selected jet pairs, construct arrays for lepton pt, bjet pt's, met, etc
      //-------------------------------------------------------------------------------
         
      double pl[4];        // Visible lepton
      double pb1[4];       // bottom on the same side as the visible lepton
      double pb2[4];       // other bottom, paired with the invisible W
      double pmiss[3];     // <unused>, pmx, pmy   missing pT
      double pmiss_lep[3]; // missing ET + lepton pT

      pl[0]= lep->E(); 
      pl[1]= lep->Px(); 
      pl[2]= lep->Py(); 
      pl[3]= lep->Pz();

      pb1[1] = jets[o].Px(); 
      pb1[2] = jets[o].Py(); 
      pb1[3] = jets[o].Pz();

      pb2[1] = jets[b].Px(); 
      pb2[2] = jets[b].Py(); 
      pb2[3] = jets[b].Pz();

      pmiss[0] = 0.; 
      pmiss[1] = metx; 
      pmiss[2] = mety;

      pmiss_lep[0] = 0.;
      pmiss_lep[1] = pmiss[1]+pl[1]; 
      pmiss_lep[2] = pmiss[2]+pl[2];

      //-------------------------------------------------------------------------------
      // calculate MT2 variables in 3 flavors: MT2b, MT2bl, MT2W
      //-------------------------------------------------------------------------------

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

      //-------------------------------------------------------------------------------
      // store a candidate with the MT2 values
      //-------------------------------------------------------------------------------

      Candidate c;
      c.chi2  = -1;
      c.mt2b  = c_mt2b;
      c.mt2w  = c_mt2w;
      c.mt2bl = c_mt2bl;
      c.j1    = -1;
      c.j2    = -1;
      c.bi    = b;
      c.oi    = o;
      c.k1    = -1;
      c.k2    = -1;
      c.match = false;
      
      mt2candidates.push_back(c);
      
    }
  }
  
  //--------------------------------------
  // return the list of candidates
  //--------------------------------------

  return mt2candidates;
}

int leadingJetIndex(vector<StopTree::LorentzVector> jets, int iskip1 = -1, int iskip2 = -1){

  int imaxpt  = -1;
  float maxpt = -1.0;

  // loop over jets
  for ( int i = 0 ; i < (int) jets.size() ; ++i ){      
      
    // skip these indices
    if( i == iskip1 ) continue;
    if( i == iskip2 ) continue;

    // store index of max pt jet
    if( jets.at(i).pt() > maxpt ){
      maxpt  = jets.at(i).pt();
      imaxpt = i;
    }
  }

  return imaxpt;
}

//------------------------------------------------------------------------------------------------

MT2struct Best_MT2Calculator(StopTree* tree, bool isData){

  //--------------------------------------------------------------------------------------------
  // This function returns the "best" MT2 variables over the list of possible b-jet pairs
  //
  // The algorithm depends on whether nbtags = 1, 2, >=3
  // nbtag = 1: use pair with tagged jet, 1 of the 2 leading non-b jets, choose minimum MT2
  // nbtag = 2: use pair with 2 b-jets
  // nbtag > 2: ignore btagging, do all combinations with 3 leading jets, choose minimum MT2
  //
  // For all cases there is a 2-fold ambiguity depending on which b is grouped with the lepton,
  // so we always take the minimum of the 2 jet pair assignments
  //
  //--------------------------------------------------------------------------------------------

  //----------------------------------------------------------------
  // start by counting jets with minimum pt and CSV thresholds
  //----------------------------------------------------------------

  const double PTMIN_BTAG = 30;
  const double BTAG_MIN   = 0.679;

  int nbtags = 0;

  assert( tree->pfjets_->size() == tree->pfjets_csv_.size() );

  vector<StopTree::LorentzVector> jets;
  vector<float> btag;

  for( unsigned int i = 0 ; i < tree->pfjets_->size() ; ++i ){
    jets.push_back( tree->pfjets_->at(i) );
    btag.push_back( tree->pfjets_csv_.at(i) );
  } 

  int n_jets = jets.size();

  vector<int> btagJets;

  for ( int i = 0 ; i < n_jets ; ++i ){      
    if( jets.at(i).pt() < PTMIN_BTAG ) continue;
    if( btag.at(i)      < BTAG_MIN   ) continue;
    nbtags++;
    btagJets.push_back(i);
  }

  //----------------------------------------------------------------
  // get the indices of the leading jets
  //----------------------------------------------------------------
  
  int imaxpt1 = -1;
  int imaxpt2 = -1;
  int imaxpt3 = -1;

  // if 1 btag, we want the indices of the 2 leading non-b jets

  if( nbtags == 1 ){

    imaxpt1 = leadingJetIndex( jets , btagJets.at(0) );           // leading non-b jet
    imaxpt2 = leadingJetIndex( jets , btagJets.at(0) , imaxpt1 ); // 2nd leading non-b jet

    // cout << endl << endl;
    // cout << "btagged jet " << btagJets.at(0) << endl;
    // cout << "jet pt's" << endl;
    
    // for ( int i = 0 ; i < n_jets ; ++i ){      
    //   cout << i << " " << jets.at(i).pt() << endl;
    // }
    
    // cout << "Leading jet     " << imaxpt  << " " << jets.at(imaxpt).pt()  << endl;
    // cout << "2nd leading jet " << imaxpt2 << " " << jets.at(imaxpt2).pt() << endl;
  }

  // if >=3 btags, we want the indices of the 3 leading jets

  if( nbtags >= 3 ){
    imaxpt1 = leadingJetIndex( jets                     ); // leading jet
    imaxpt2 = leadingJetIndex( jets  , imaxpt1          ); // 2nd leading jet
    imaxpt3 = leadingJetIndex( jets  , imaxpt1 , imaxpt2); // 3rd leading jet

    // cout << endl << endl;
    // cout << "jet pt's" << endl;
    
    // for ( int i = 0 ; i < n_jets ; ++i ){      
    //   cout << i << " " << jets.at(i).pt() << endl;
    // }
    
    // cout << "Leading jet     " << imaxpt1 << " " << jets.at(imaxpt1).pt() << endl;
    // cout << "2nd leading jet " << imaxpt2 << " " << jets.at(imaxpt2).pt() << endl;
    // cout << "3rd leading jet " << imaxpt3 << " " << jets.at(imaxpt3).pt() << endl;

  }


  //----------------------------------------------------------------
  // get all the MT2 values for all jet pairings
  //----------------------------------------------------------------

  list<Candidate> mt2candidates = MT2Calculator(tree, isData );

  list<Candidate>::iterator candIter;

  float mt2w_min  = 1e10;
  float mt2b_min  = 1e10;
  float mt2bl_min = 1e10;

  //------------------------------------------------------------------------------
  // loop over jet pairings, choose jet pairing best on algorithm described above
  // bi is the index of the jet grouped with the lepton
  // oi is the index of the jet not grouped with the lepton
  //------------------------------------------------------------------------------
      
  for(candIter = mt2candidates.begin() ; candIter != mt2candidates.end() ; candIter++ ){

    //--------------
    // 1 btags
    //--------------
      
    if( nbtags == 1 ){

      // require that either bi or oi is the btagged jet
      if( (*candIter).bi != btagJets.at(0) && (*candIter).oi != btagJets.at(0) ) continue;

      bool goodCombo = false;

      // bi is btagged, oi is one of the 2 leading non-b jets
      if( (*candIter).bi == btagJets.at(0) && ( (*candIter).oi == imaxpt1 || (*candIter).oi == imaxpt2) ) goodCombo = true;

      // oi is btagged, bi is one of the 2 leading non-b jets
      if( (*candIter).oi == btagJets.at(0) && ( (*candIter).bi == imaxpt1 || (*candIter).bi == imaxpt2) ) goodCombo = true;

      if( !goodCombo) continue;
    }

    //--------------
    // 2 btags
    //--------------

    else if( nbtags == 2 ){
	
      // require bi is one of the 2 btagged jets
      if( (*candIter).bi != btagJets.at(0) && (*candIter).bi != btagJets.at(1) ) continue; 
	
      // require oi is one of the 2 btagged jets
      if( (*candIter).oi != btagJets.at(0) && (*candIter).oi != btagJets.at(1) ) continue; 

    }

    //--------------
    // >=3 btags
    //--------------
      
    else if( nbtags >= 3 ){
	
      // bi must be one of 3 leading jets
      if( (*candIter).bi != imaxpt1 && (*candIter).bi != imaxpt2 && (*candIter).bi != imaxpt3 ) continue;

      // oi must be one of 3 leading jets
      if( (*candIter).oi != imaxpt1 && (*candIter).oi != imaxpt2 && (*candIter).oi != imaxpt3 ) continue;
    }
      
    else{
      cout << __FILE__ << " " << __LINE__ << " should never get here!!!" << endl;
    }

    float mt2w  = (*candIter).mt2w;
    float mt2b  = (*candIter).mt2b;
    float mt2bl = (*candIter).mt2bl;

    if( mt2w  < mt2w_min  ) mt2w_min  = mt2w;
    if( mt2b  < mt2b_min  ) mt2b_min  = mt2b;
    if( mt2bl < mt2bl_min ) mt2bl_min = mt2bl;

  }

  if( mt2w_min  > 9000 ) mt2w_min  = 0.0;
  if( mt2b_min  > 9000 ) mt2b_min  = 0.0;
  if( mt2bl_min > 9000 ) mt2bl_min = 0.0;

  MT2struct m;
  m.mt2w  = mt2w_min;
  m.mt2b  = mt2b_min;
  m.mt2bl = mt2bl_min;
  m.chi2  = -1;

  return m;

}

//-------------------------------------------
// efficiency for dilepton trigger
//-------------------------------------------

float getdltrigweight(int id1, int id2){ 
  if (abs(id1)==11 && abs(id2)==11) return 0.95;
  if (abs(id1)==13 && abs(id2)==13) return 0.88;
  if (abs(id1)!=abs(id2)) return 0.92;
  return -999.;

}

//-------------------------------------------
// efficiency for single lepton trigger
//-------------------------------------------

float getsltrigweight(int id1, float pt, float eta){

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

//-----------------------------------------------------------------------
// basic event selection:
// >=1 good lepton, rho cut, MET filters, remove 2 nearby lepton events
//-----------------------------------------------------------------------

bool passEvtSelection(const StopTree *sTree, TString name) 
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

//-------------------------------------------
// good lepton + veto isolated track
//-------------------------------------------

bool passOneLeptonSelection(const StopTree *sTree, bool isData) 
{
  //single lepton selection for 8 TeV 53 analysis
  if ( !passSingleLeptonSelection(sTree, isData) ) return false;

  //pass isolated track veto
  //unfortunately changed default value to 9999.
  if ( sTree->pfcandpt10_ <9998. && sTree->pfcandiso10_ < 0.1 ) return false;

  return true;

}

//-------------------------------------------
// dilepton selection
//-------------------------------------------

bool passTwoLeptonSelection(const StopTree *sTree, bool isData) 
{
  //single lepton selection for 8 TeV 53 analysis
  if ( !passDileptonSelection(sTree, isData) ) return false;

  //apply isolated track veto in addition to 2 leptons
  //default value for this one is -999
  if ( sTree->trkpt10loose_ >0. && sTree->trkreliso10loose_ < 0.1 ) return false;

  return true;

}

//-------------------------------------------
// the isolated track veto
//-------------------------------------------

bool passIsoTrkVeto(const StopTree *sTree) 
{

  //pass isolated track veto
  //unfortunately changed default value to 9999.
  if ( sTree->pfcandpt10_ <9998. && sTree->pfcandiso10_ < 0.1 ) return false;

  return true;

}

//-------------------------------------------
// the isolated track veto
//-------------------------------------------

bool passIsoTrkVeto_v2(const StopTree *sTree) 
{

  //pass isolated track veto
  //unfortunately changed default value to 9999.
  if ( sTree->pfcandpt10_ <9998. && sTree->pfcandiso10_ < 0.1 ) return false;
  if ( sTree->pfcandpt5_  <9998. && abs(sTree->pfcandid5_)==13 && sTree->pfcandiso5_ < 0.2) return false;
  if ( sTree->pfcandpt5_  <9998. && abs(sTree->pfcandid5_)==11 && sTree->pfcandiso5_ < 0.2) return false;

  return true;

}

//-------------------------------------------
// >=1 selected lepton and trigger
//-------------------------------------------

bool passSingleLeptonSelection(const StopTree *sTree, bool isData) 
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

//-------------------------------------------
// the dilepton selection
//-------------------------------------------

bool passDileptonSelection(const StopTree *sTree, bool isData) 
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

//-------------------------------------------
// lepton + isolated track selection CR5
//-------------------------------------------

bool passLepPlusIsoTrkSelection(const StopTree *sTree, bool isData) 
{
  //single lepton plus iso trk selection for 8 TeV 53 analysis

  //at least one lepton
  if ( !passSingleLeptonSelection(sTree, isData) ) return false;

  //pass isolated track requirement
  //unfortunately changed default value to 9999.
  if ( sTree->pfcandpt10_ > 9990. || sTree->pfcandiso10_ > 0.1 ) return false;

  return true;

}

//-------------------------------------------
// on-the-fly MET phi corrections
//-------------------------------------------

pair<float,float> getPhiCorrMET( float met, float metphi, int nvtx, bool ismc){

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

//---------------------------------------
// utility functions
//---------------------------------------

float vtxweight_n( const int nvertices, TH1F *hist, bool isData ) {

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

float getdphi( float phi1 , float phi2 ){
  float dphi = fabs( phi1 - phi2 );
  if( dphi > TMath::Pi() ) dphi = TMath::TwoPi() - dphi;
  return dphi;
}

float getMT( float pt1 , float phi1 , float pt2 , float phi2 ){

  float dphi = getdphi(phi1, phi2);
  return sqrt( 2 * ( pt1 * pt2 * (1 - cos( dphi ) ) ) );

}

float dRbetweenVectors(StopTree::LorentzVector vec1,StopTree::LorentzVector vec2 ){ 
  float dphi = std::min(::fabs(vec1.Phi() - vec2.Phi()), 2 * M_PI - fabs(vec1.Phi() - vec2.Phi()));
  float deta = vec1.Eta() - vec2.Eta();

  return sqrt(dphi*dphi + deta*deta);
}
