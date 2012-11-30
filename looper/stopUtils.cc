#include "stopUtils.h"

/*
//--------------------------------------------------------------------

double weight3D( int pv1, int pv2, int pv3 ) {

  using std::min;

  int npm1 = min(pv1,49);
  int np0 = min(pv2,49);
  int npp1 = min(pv3,49);

  return Weight3D[npm1][np0][npp1];

}

//--------------------------------------------------------------------                

void weight3D_init( std::string WeightFileName ) { 

  TFile *infile = new TFile(WeightFileName.c_str());
  TH1F *WHist = (TH1F*)infile->Get("WHist");

  
  // Check if the histogram exists           
  if (!WHist) {
    cout << "Error, could not find the histogram WHist in the file "
	 << WeightFileName << ", quitting" << endl;
    exit(0);
  }

  for (int i=0; i<50; i++) 
    for(int j=0; j<50; j++)
      for(int k=0; k<50; k++) {
	Weight3D[i][j][k] = WHist->GetBinContent(i+1,j+1,k+1);
      }

  cout << " 3D Weight Matrix initialized! " << endl;

  delete infile;

  return;


}
*/
//--------------------------------------------------------------------

void checkElectron( int elidx ){

  cout << "Check electron" << endl;
  cout << "Pass all    " << pass_electronSelection( elidx , electronSelection_ssV5			) << endl;
  cout << "Pass ID     " << pass_electronSelection( elidx , electronSelection_ssV5_noIso		) << endl;
  cout << "Pass iso    " << pass_electronSelection( elidx , electronSelection_ssV5_iso	        	) << endl;
  cout << "VBTF90      " << pass_electronSelection( elidx , 1ll<<ELEID_VBTF_90_HLT_CALOIDT_TRKIDVL	) << endl;
  cout << "PV          " << pass_electronSelection( elidx , 1ll<<ELEIP_PV_OSV2				) << endl;
  cout << "nomuon      " << pass_electronSelection( elidx , 1ll<<ELENOMUON_010				) << endl;
  cout << "hitpattern  " << pass_electronSelection( elidx , 1ll<<ELENOTCONV_HITPATTERN			) << endl;
  cout << "convrej     " << pass_electronSelection( elidx , 1ll<<ELENOTCONV_DISTDCOT002			) << endl;
  cout << "pt10        " << pass_electronSelection( elidx , 1ll<<ELEPT_010				) << endl;
  cout << "eta25       " << pass_electronSelection( elidx , 1ll<<ELEETA_250				) << endl;
  cout << "transition  " << pass_electronSelection( elidx , 1ll<<ELE_NOT_TRANSITION			) << endl;
  cout << "HLT iso     " << pass_electronSelection( elidx , 1ll<<ELEISO_ECAL_RELNT020_NPS		) << endl;
  cout << "offline iso " << pass_electronSelection( elidx , 1ll<<ELEISO_RELNT015			) << endl;

}

void checkMuon( int muidx ){

  cout << "Check muon" << endl;
  cout << "Pass all  " <<  muonId(muidx , OSGeneric_v3)                                            << endl;
  cout << "Pass ID   " <<  muonIdNotIsolated(muidx , OSGeneric_v3 )                                << endl;
  cout << "Pass iso  " <<  ( muonIsoValue(muidx,false) < 0.15 )                                    << endl;
  cout << "eta24     " <<  ( TMath::Abs(mus_p4()[muidx].eta()) < 2.4)                         << endl;
  cout << "chi2/ndf  " <<  ( mus_gfit_chi2().at(muidx)/mus_gfit_ndof().at(muidx) < 10)   << endl;
  cout << "global    " <<  ( ((mus_type().at(muidx)) & (1<<1)) != 0)                          << endl;
  cout << "tracker   " <<  ( ((mus_type().at(muidx)) & (1<<2)) != 0)                          << endl;
  cout << "nhits     " <<  ( mus_validHits().at(muidx) > 10)                                  << endl;
  cout << "stahits   " <<  ( mus_gfit_validSTAHits().at(muidx) != 0)                          << endl;
  cout << "d0PV      " <<  ( TMath::Abs(mud0PV_smurfV3(muidx)) < 0.02)                             << endl;
  cout << "dzPV      " <<  ( TMath::Abs(mudzPV_smurfV3(muidx)) < 1  )                              << endl;
  cout << "dpt/pt    " <<  ( mus_ptErr().at(muidx)/mus_p4().at(muidx).pt()<0.1)          << endl;

}

//--------------------------------------------------------------------                                                                                                                                               

pair<float, float> ScaleMET( pair<float, float> p_met, LorentzVector p4_dilep, double rescale){
  float met = p_met.first;
  float metPhi = p_met.second;
  float metx = met*cos(metPhi);
  float mety = met*sin(metPhi);

  float lepx = p4_dilep.Px();
  float lepy = p4_dilep.Py();

  //hadronic component of MET (well, mostly), scaled                                                                                                                                                                 
  float metHx = (metx + lepx)*rescale;
  float metHy = (mety + lepy)*rescale;
  float metNewx = metHx - lepx;
  float metNewy = metHy - lepy;
  float metNewPhi = atan2(metNewy, metNewx);

  pair<float, float> p_met2 = make_pair(sqrt(metNewx*metNewx + metNewy*metNewy), metNewPhi);
  return p_met2;
}

//--------------------------------------------------------------------                                                                                 

pair<float,float> getPhiCorrMET( float met, float metphi, int nvtx, bool ismc ) {                                                 

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

//--------------------------------------------------------------------                                                                                                                                               

pair<float,float> getTrackerMET( P4 *lep, double deltaZCut, bool dolepcorr )
{

  if ( cms2.vtxs_sumpt().empty() ) return make_pair(-999.,-999.);

  float pX = 0.;
  float pY = 0.;

  if (dolepcorr ) {
    pX -= lep->px();
    pY -= lep->py();
  }

  for (unsigned int i=0; i<cms2.pfcands_particleId().size(); ++i){
    if ( cms2.pfcands_charge().at(i)==0 ) continue;

    if ( dolepcorr && dRbetweenVectors( cms2.pfcands_p4().at(i) , *lep ) < 0.1 ) continue;

    int trkIndex = cms2.pfcands_trkidx().at(i);
    if (trkIndex<0) continue;

    //    double dzpv = dzPV(cms2.trks_vertex_p4()[trkIndex], cms2.trks_trk_p4()[trkIndex], cms2.vtxs_position().front());
    // this is neede for the slim
    double dzpv = trks_dz_pv(trkIndex,0).first;

    if ( fabs(dzpv) > deltaZCut) continue;

    pX -= cms2.pfcands_p4().at(i).px();
    pY -= cms2.pfcands_p4().at(i).py();
  }

  pair<float, float> trkmet = make_pair( sqrt( pX*pX + pY*pY ), atan2( pY , pX ) );
  return trkmet;

}


//--------------------------------------------------------------------                                                                                                                                               

pair<float,float> Type1PFMET( VofP4 jets_p4 , vector<float> cors , vector<float> l1cors , float minpt ){

  float metx = evt_pfmet() * cos( evt_pfmetPhi() );
  float mety = evt_pfmet() * sin( evt_pfmetPhi() );

  assert( jets_p4.size() == cors.size() );

  for( unsigned int i = 0 ; i < jets_p4.size() ; ++i ){
    float corrpt = jets_p4.at(i).pt() * cors.at(i);
    if( corrpt < minpt ) continue;
    float l1corr = (l1cors.size()==0) ? 1. : l1cors.at(i);
    metx += jets_p4.at(i).px() * l1corr - jets_p4.at(i).px() * cors.at(i);
    mety += jets_p4.at(i).py() * l1corr - jets_p4.at(i).py() * cors.at(i);
  }

  pair<float, float> type1met = make_pair( sqrt( metx*metx + mety*mety ), atan2( mety , metx ) );
  return type1met;
}

/////--------------------------------------------------------------------

float getDataMCRatio(float eta){
  if (eta >=0.0 && eta < 0.5) return 1.052;
  if (eta >=0.5 && eta < 1.1) return 1.057;
  if (eta >=1.1 && eta < 1.7) return 1.096;
  if (eta >=1.7 && eta < 2.3) return 1.134;
  if (eta >=2.3 && eta < 5.0) return 1.288;
  return 1.0;
}

/////--------------------------------------------------------------------

pair<float,float> Type1PFMETSmear(JetSmearer* jetSmearer, bool isData,
				  VofP4 jets_p4 , float met, float metphi){
  float metx = met * cos( metphi );
  float mety = met * sin( metphi );

  if (!isData){
    UInt_t seed = 0;

    for (unsigned int i=0; i<jets_p4.size(); ++i)
      seed += jets_p4.at(i).Pt()*1000;

    TRandom3 kicker(seed);

    for (unsigned int i=0; i<jets_p4.size(); ++i){
      LorentzVector recJet = jets_p4.at(i);

      float c = getDataMCRatio(recJet.eta());
      double sigma = getJetResolution(recJet, jetSmearer);
      double alpha = kicker.Gaus(0.0, TMath::Sqrt(c*c-1.0)*sigma);

      metx -= alpha*recJet.px();
      mety -= alpha*recJet.py();
    }
  }

  return make_pair( sqrt( metx*metx + mety*mety ), atan2( mety , metx ) );
}

/////--------------------------------------------------------------------  

float getMT( float leppt , float lepphi , float met , float metphi ) {
  float dphi = fabs( lepphi - metphi );
  if( dphi > TMath::Pi() ) dphi = TMath::TwoPi() - dphi;
  return sqrt( 2 * ( leppt * met * (1 - cos( dphi ) ) ) );
}

//--------------------------------------------------------------------                                                                                                                                               
float dz_trk_vtx( const unsigned int trkidx, const unsigned int vtxidx ){

  //  return ((trks_vertex_p4()[trkidx].z()-vtxs_position()[vtxidx].z()) - ((trks_vertex_p4()[trkidx].x()-vtxs_position()[vtxidx].x()) * trks_trk_p4()[trkidx].px() + (trks_vertex_p4()[trkidx].y() - vtxs_position()[vtxidx].y()) * trks_trk_p4()[trkidx].py())/trks_trk_p4()[trkidx].pt() * trks_trk_p4()[trkidx].pz()/trks_trk_p4()[trkidx].pt());

  // this is for the slim
  return trks_dz_pv(trkidx,vtxidx).first;

}

float trackIso( int thisPf , float coneR , float dz_thresh , bool dovtxcut , float pt_thresh ){

  float iso = 0.0;

  for (int ipf = 0; ipf < (int)cms2.pfcands_p4().size(); ipf++) {

    if( ipf == thisPf                 ) continue; // skip this PFCandidate                                                                                                                       
    if( cms2.pfcands_charge().at(ipf) == 0 ) continue; // skip neutrals                                                                                                                          
    if( cms2.pfcands_p4().at(ipf).pt() < pt_thresh ) continue; // skip pfcands below pt threshold                                                                                                

    if( dRbetweenVectors( pfcands_p4().at(ipf) , pfcands_p4().at(thisPf) ) > coneR ) continue;

    int itrk = cms2.pfcands_trkidx().at(ipf);

    if( itrk >= (int)trks_trk_p4().size() || itrk < 0 ){
      //note: this should only happen for electrons which do not have a matched track                                                                                                            
      //currently we are just ignoring these guys                                                                                                                                                
      continue;
    }

    //----------------------------------------                                                                                                                                                   
    // find closest PV and dz w.r.t. that PV                                                                                                                                                     
    //----------------------------------------                                                                                                                                                   

    float mindz = 999.;
    int vtxi    = -1;

    if (dovtxcut) {
      for (unsigned int ivtx = 0; ivtx < cms2.vtxs_position().size(); ivtx++) {

        if(!isGoodVertex(ivtx)) continue;

        float mydz = dz_trk_vtx(itrk,ivtx);
	//        fillOverFlow( h_dz_vtx_trk , mydz );

        if (fabs(mydz) < fabs(mindz)) {
          mindz = mydz;
          vtxi = ivtx;
        }

      }

      //----------------------------------------------------------------------------                                                                                                               
      // require closest PV is signal PV, dz cut, exclude tracks near hyp leptons                                                                                                                  
      //----------------------------------------------------------------------------                                                                                                               

      if ( vtxi != 0 )     continue;
    } else {
      mindz = dz_trk_vtx(itrk,0);
    }
    if ( fabs(mindz) > dz_thresh )     continue;

    //---------------------------------------                                                                                                                                                    
    // passes cuts, add up isolation value                                                                                                                                                       
    //---------------------------------------                                                                                                                                                    

    iso += cms2.pfcands_p4().at(ipf).pt();

  }

  return iso;
}






//--------------------------------------------------------------------                                                                                                                                               

double dRbetweenVectors(const LorentzVector &vec1,
                        const LorentzVector &vec2 ){
  double dphi = std::min(::fabs(vec1.Phi() - vec2.Phi()), 2 * M_PI - fabs(vec1.Phi() - vec2.Phi()));
  double deta = vec1.Eta() - vec2.Eta();
  return sqrt(dphi*dphi + deta*deta);
}

//--------------------------------------------------------------------                                                                                                                                               

bool isGenBMatched ( LorentzVector p4, float dR ) {

  //For now only checking status 3 in dR                                                                                                                                                                             
  for (unsigned int igen = 0; igen < genps_p4().size(); igen++) {

    int id = genps_id().at(igen);
    if( abs(id)!=5 ) continue;

    if( dRbetweenVectors( p4 , genps_p4().at(igen) ) < dR ) return true;

  }
  return false;
}

//--------------------------------------------------------------------                                                                                                                                               

int isGenQGMatched ( LorentzVector p4, float dR ) {
  //Start from the end that seems to have the decay products of the W first                                                                                                                                          
  for (int igen = (genps_p4().size()-1); igen >-1; igen--) {
    float deltaR = dRbetweenVectors( p4 , genps_p4().at(igen) );
    if ( deltaR > dR ) continue;
    int id = genps_id().at(igen);
    int mothid = genps_id_mother().at(igen);
    // cout<<"status 3 particle ID "<<id<<" mother "<<mothid                                                                                                                                                         
    //  <<" dR to jet "<<deltaR<<endl;                                                                                                                                                                               
    if (abs(id)<6 && abs(mothid)==24)
      return (mothid>0) ? 2 : -2;
    if (abs(id)==5 && abs(mothid)==6)
      return (mothid>0) ? 1 : -1;
    if (abs(id)==21) return 3;
    if (abs(id)<6) return 4;
  }
  return -9;
}

//--------------------------------------------------------------------                                                                                                                                               

float dRGenJet ( LorentzVector p4, float genminpt) {

  //return dR to closest gen-jet with pT > genminpt                                                                                                                                                                  
  float mindeltaR = 9999.;
  for (unsigned int igen = 0; igen < genjets_p4().size(); igen++) {
    LorentzVector vgenj = genjets_p4().at(igen);
    if ( vgenj.Pt() < genminpt ) continue;
    float deltaR = dRbetweenVectors( p4 , vgenj );
    if ( deltaR< mindeltaR ) mindeltaR = deltaR;
  }
  return mindeltaR;

}

//--------------------------------------------------------------------                                                                                                                                               

float getminjdr( VofiP4 jets, LorentzVector *particle ) {
  float mindr = 9999.;
  if (jets.size()==0 || particle==0) return mindr;
  for ( unsigned int ijet = 0; ijet<jets.size(); ++ijet ) {
    float partjdr = dRbetweenVectors(jets.at(ijet).p4obj,*particle);
    if ( partjdr<mindr ) mindr = partjdr;
  }
  return mindr;
}


//--------------------------------------------------------------------                                                                                                                                               

int isSSVMTagged ( int ijet ) {

  if ( pfjets_simpleSecondaryVertexHighEffBJetTag().at(ijet) > 1.74 )
    return 1;

  return 0;

}

//--------------------------------------------------------------------                                                                                                                                               

int isCSVTagged ( int ijet ) {

  float discrim = pfjets_combinedSecondaryVertexBJetTag().at(ijet);
  if ( discrim > 0.679 ) return 1;//medium                                                                                                                                                                           
  if ( discrim > 0.244 ) return 2;//loose                                                                                                                                                                            

  return 0;

}

//--------------------------------------------------------------------                                                                                                                                               

bool isBTagged ( LorentzVector p4, VofP4 bJets ) {

  for( int ijet = 0 ; ijet < (int)bJets.size() ; ijet++ ){
    if( dRbetweenVectors( p4 , bJets.at(ijet) ) < 0.4 ) return true;
  }

  return false;

}

//--------------------------------------------------------------------                                                                                                                                               

int getLeptonMatchIndex ( LorentzVector *jet, LorentzVector *lep1, LorentzVector *lep2, float dR ) {
  if (lep1 && dRbetweenVectors( *lep1 , *jet ) < dR) return 1;
  if (lep2 && dRbetweenVectors( *lep2 , *jet ) < dR) return 2;
  return -1;

}

//--------------------------------------------------------------------                                                                                                                                               

int findTriggerIndex(TString trigName)
{
  vector<TString>::const_iterator begin_it = hlt_trigNames().begin();
  vector<TString>::const_iterator end_it = hlt_trigNames().end();
  vector<TString>::const_iterator found_it = find(begin_it, end_it, trigName);
  if(found_it != end_it) return found_it - begin_it;
  return -1;
}

//--------------------------------------------------------------------                                                                                                                                               

bool objectPassTrigger(const LorentzVector &obj, const std::vector<LorentzVector> &trigObjs, float pt)
{

  float drMin = 999.99;
  for (size_t i = 0; i < trigObjs.size(); ++i)
    {
      if (trigObjs[i].Pt() < pt) continue;
      float dr = dRbetweenVectors(trigObjs[i], obj);
      if (dr < drMin) drMin = dr;
    }

  if (drMin < 0.1) return true;
  return false;

}

TString triggerName(TString triggerPattern){

  //-------------------------------------------------------                                                                                                                                                          
  // get exact trigger name corresponding to given pattern                                                                                                                                                           
  //-------------------------------------------------------                                                                                                                                                          

  bool    foundTrigger  = false;
  TString exact_hltname = "";

  for( unsigned int itrig = 0 ; itrig < hlt_trigNames().size() ; ++itrig ){
    if( TString( hlt_trigNames().at(itrig) ).Contains( triggerPattern ) ){
      foundTrigger  = true;
      exact_hltname = hlt_trigNames().at(itrig);
      break;
    }
  }

  if( !foundTrigger) return "TRIGGER_NOT_FOUND";

  return exact_hltname;

}


bool objectPassTrigger(const LorentzVector &obj, char* trigname, float drmax = 0.1 ){

  TString exact_trigname = triggerName( trigname );

  if( exact_trigname.Contains("TRIGGER_NOT_FOUND") ){
    cout << __FILE__ << " " << __LINE__ << " Error! couldn't find trigger name " << trigname << endl;
    return false;
  }

  std::vector<LorentzVector> trigp4 = hlt_trigObjs_p4()[findTriggerIndex(exact_trigname)];

  if( trigp4.size() == 0 ) return false;

  for (unsigned int i = 0; i < trigp4.size(); ++i){
    float dr = dRbetweenVectors(trigp4[i], obj);
    if( dr < drmax ) return true;
  }

  return false;
}


//--------------------------------------------------------------------

void fillUnderOverFlow(TH1F *h1, float value, float weight)
{
  float min = h1->GetXaxis()->GetXmin();
  float max = h1->GetXaxis()->GetXmax();

  if (value > max) value = h1->GetBinCenter(h1->GetNbinsX());
  if (value < min) value = h1->GetBinCenter(1);

  h1->Fill(value, weight);
}

//--------------------------------------------------------------------

void fillUnderOverFlow(TH2F *h2, float xvalue, float yvalue, float weight)
{
  float maxx = h2->GetXaxis()->GetXmax();
  float minx = h2->GetXaxis()->GetXmin();
  float maxy = h2->GetYaxis()->GetXmax();
  float miny = h2->GetYaxis()->GetXmin();

  if (xvalue > maxx) xvalue = h2->GetXaxis()->GetBinCenter(h2->GetNbinsX());
  if (xvalue < minx) xvalue = h2->GetXaxis()->GetBinCenter(1);
  if (yvalue > maxy) yvalue = h2->GetYaxis()->GetBinCenter(h2->GetNbinsY());
  if (yvalue < miny) yvalue = h2->GetYaxis()->GetBinCenter(1);

  h2->Fill(xvalue, yvalue, weight);
}

//--------------------------------------------------------------------

// void fillUnderOverFlow(TProfile *h2, float xvalue, float yvalue)
// {
//   float maxx = h2->GetXaxis()->GetXmax();
//   float minx = h2->GetXaxis()->GetXmin();
//   float maxy = h2->GetYaxis()->GetXmax();
//   float miny = h2->GetYaxis()->GetXmin();

//   if (xvalue > maxx) xvalue = h2->GetXaxis()->GetBinCenter(h2->GetNbinsX());
//   if (xvalue < minx) xvalue = h2->GetXaxis()->GetBinCenter(1);
//   if (yvalue > maxy) yvalue = h2->GetYaxis()->GetBinCenter(h2->GetNbinsY());
//   if (yvalue < miny) yvalue = h2->GetYaxis()->GetBinCenter(1);

//   h2->Fill(xvalue, yvalue);
// }

//--------------------------------------------------------------------

void fillOverFlow(TH1F *h1, float value, float weight)
{
  float max = h1->GetXaxis()->GetXmax();
  if (value > max) value = h1->GetBinCenter(h1->GetNbinsX());
  h1->Fill(value, weight);
}

//--------------------------------------------------------------------

void fillOverFlow(TH2F *h2, float xvalue, float yvalue, float weight)
{
  float maxx = h2->GetXaxis()->GetXmax();
  float maxy = h2->GetYaxis()->GetXmax();

  if (xvalue > maxx) xvalue = h2->GetXaxis()->GetBinCenter(h2->GetNbinsX());
  if (yvalue > maxy) yvalue = h2->GetYaxis()->GetBinCenter(h2->GetNbinsY());

  h2->Fill(xvalue, yvalue, weight);
}

//--------------------------------------------------------------------

void fillHistos(TH1F *h1[4][4],float value, float weight, int myType, int nJetsIdx)
{
  fillUnderOverFlow(h1[myType][nJetsIdx], value, weight);      
  fillUnderOverFlow(h1[myType][3],        value, weight);      
  fillUnderOverFlow(h1[3][nJetsIdx],      value, weight);      
  fillUnderOverFlow(h1[3][3],             value, weight);      
}

//--------------------------------------------------------------------

void fillHistos(TH2F *h2[4][4],float xvalue, float yvalue, float weight, int myType, int nJetsIdx)
{
  fillUnderOverFlow(h2[myType][nJetsIdx], xvalue, yvalue, weight);      
  fillUnderOverFlow(h2[myType][3],        xvalue, yvalue, weight);      
  fillUnderOverFlow(h2[3][nJetsIdx],      xvalue, yvalue, weight);      
  fillUnderOverFlow(h2[3][3],             xvalue, yvalue, weight);      
}

//--------------------------------------------------------------------

void fillHistos(TProfile *h2[4][4],float xvalue, float yvalue, int myType, int nJetsIdx)
{
  h2[myType][nJetsIdx] -> Fill(xvalue, yvalue);      
  h2[myType][3]        -> Fill(xvalue, yvalue);      
  h2[3][nJetsIdx]      -> Fill(xvalue, yvalue);      
  h2[3][3]             -> Fill(xvalue, yvalue);      
}

//--------------------------------------------------------------------
