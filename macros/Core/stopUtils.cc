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

      pb1[0] = jets[o].E();  
      pb1[1] = jets[o].Px(); 
      pb1[2] = jets[o].Py(); 
      pb1[3] = jets[o].Pz();

      pb2[0] = jets[b].E();
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

      mt2_event.set_momenta( pb1, pb2, pmiss_lep );
      mt2_event.set_mn( 0.0 );   // Invisible particle mass
      double c_mt2b = mt2_event.get_mt2();
         
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

  MT2struct m;
  m.mt2w  = mt2w_min;
  m.mt2b  = mt2b_min;
  m.mt2bl = mt2bl_min;

  return m;

}

