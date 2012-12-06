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

  //LorentzVector* lep = &tree->lep1_;
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
  //int n_jets = tree->pfjets_->size();

  //assert( jets.size() == btag.size() );

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
      //double pt_b = tree->pfjets_->at(b).pt();

      if ( btag[b] >= BTAG_MIN && pt_b < PTMIN_BTAG )
        continue;

      if ( btag[b] < BTAG_MIN && pt_b < PTMIN_B )
        continue;

      double pt_o = jets[o].Pt();
      //double pt_o = tree->pfjets_->at(o).pt();

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

      // pl[0]= (&tree->lep1_)->E(); 
      // pl[1]= (&tree->lep1_)->Px(); 
      // pl[2]= (&tree->lep1_)->Py(); 
      // pl[3]= (&tree->lep1_)->Pz();

      pb1[0] = jets[o].E();  
      pb1[1] = jets[o].Px(); 
      pb1[2] = jets[o].Py(); 
      pb1[3] = jets[o].Pz();

      pb2[0] = jets[b].E();
      pb2[1] = jets[b].Px(); 
      pb2[2] = jets[b].Py(); 
      pb2[3] = jets[b].Pz();

      // double pt_o = tree->pfjets_->at(o).pt();
      // p

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

//------------------------------------------------------------------------------------------------
