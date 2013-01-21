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
#include <iostream>
#include <fstream>

#include "MT2Utility.h"
#include "mt2bl_bisect.h"
#include "mt2w_bisect.h"

using namespace Stop;
using namespace std;

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > LorentzVector;


unsigned int getNJets(){

  unsigned int njets=0;

  for ( unsigned int i=0; i<stopt.pfjets().size() ; i++) {
    
    if( stopt.pfjets().at(i).pt()<30 )  continue;
    if( fabs(stopt.pfjets().at(i).eta())>2.4 )  continue;
    //    if(stopt.pfjets_beta2().at(i)<=0.01) continue;

    njets++;

  }

  return njets;

}

vector<int> getBJetIndex(double discr, int iskip1, int iskip2)
{

  vector<int> btagJets;

  for ( unsigned int i=0; i<stopt.pfjets().size() ; i++) {

    if( stopt.pfjets().at(i).pt()<30 )  continue;
    if( fabs(stopt.pfjets().at(i).eta())>2.4 )  continue;
    if( stopt.pfjets_csv().at(i)      < discr   ) continue;

    // skip these indices                                                                                                                                                                                  
    if( i == iskip1 ) continue;
    if( i == iskip2 ) continue;

    //    nbtags++;                                                                                                                                                        
    btagJets.push_back(i);

  }

  return btagJets;

}



int leadingJetIndex(vector<LorentzVector> jets, int iskip1 = -1, int iskip2 = -1){

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
//this is for the jetResolutions
//------------------------------------------------------------------------------------------------

float getDataMCRatio(float eta){
  if (eta >=0.0 && eta < 0.5) return 1.052;
  if (eta >=0.5 && eta < 1.1) return 1.057;
  if (eta >=1.1 && eta < 1.7) return 1.096;
  if (eta >=1.7 && eta < 2.3) return 1.134;
  if (eta >=2.3 && eta < 5.0) return 1.288;
  return 1.0;
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

bool passEvtSelection(TString name) 
{

  //rho requirement
  if ( stopt.rhovor()<0. || stopt.rhovor()>=40. ) return false;

  if (!name.Contains("T2")) {
    //met filters
    if ( stopt.csc()      != 0 ) return false;
    if ( stopt.hbhe()     != 1 ) return false;
    if ( stopt.hcallaser()!= 1 ) return false;
    if ( stopt.ecaltp()   != 1 ) return false;
    if ( stopt.trkfail()  != 1 ) return false;
    if ( stopt.eebadsc()  != 1 ) return false;
    if ( stopt.hbhenew()  != 1 ) return false;
  }

  //at least 1 lepton
  if ( stopt.ngoodlep() < 1 ) return false;

  //if have more than 1 lepton, remove cases where have 2 close together
  if ( stopt.ngoodlep() > 1 && 
       dRbetweenVectors( stopt.lep1() ,  stopt.lep2() )<0.1 ) return false;

  return true;

}

//-------------------------------------------
// good lepton + veto isolated track
//-------------------------------------------

bool passOneLeptonSelection( bool isData) 
{
  //single lepton selection for 8 TeV 53 analysis
  if ( !passSingleLeptonSelection(isData) ) return false;

  //pass isolated track veto
  //unfortunately changed default value to 9999.
  if ( stopt.pfcandpt10() <9998. && stopt.pfcandiso10() < 0.1 ) return false;

  return true;

}

//-------------------------------------------
// dilepton selection
//-------------------------------------------

bool passTwoLeptonSelection(bool isData) 
{
  //single lepton selection for 8 TeV 53 analysis
  if ( !passDileptonSelection(isData) ) return false;

  //apply isolated track veto in addition to 2 leptons
  //default value for this one is -999
  if ( stopt.trkpt10loose() >0. && stopt.trkreliso10loose() < 0.1 ) return false;

  return true;

}

//-------------------------------------------
// the isolated track veto
//-------------------------------------------

bool passIsoTrkVeto() 
{

  //pass isolated track veto
  //unfortunately changed default value to 9999.
  if ( stopt.pfcandpt10() <9998. && stopt.pfcandiso10() < 0.1 ) return false;

  return true;

}

//-------------------------------------------
// the isolated track veto
//-------------------------------------------

bool passIsoTrkVeto_v2() 
{

  //pass isolated track veto
  //unfortunately changed default value to 9999.
  if ( stopt.pfcandpt10() <9998. && stopt.pfcandiso10() < 0.1 ) return false;
  if ( stopt.pfcandpt5()  <9998. && abs(stopt.pfcandid5())==13 && stopt.pfcandiso5() < 0.2) return false;
  if ( stopt.pfcandpt5()  <9998. && abs(stopt.pfcandid5())==11 && stopt.pfcandiso5() < 0.2) return false;

  return true;

}

//-------------------------------------------
// >=1 selected lepton and trigger
//-------------------------------------------

bool passSingleLeptonSelection(bool isData) 
{
  //single lepton selection for 8 TeV 53 analysis

  //at least one lepton
  if ( stopt.ngoodlep() < 1 ) return false;

  //lepton flavor - trigger, pt and eta requirements
  if ( stopt.lep1().Pt() < 30 )          return false;
  if ( fabs( stopt.pflep1().Pt() - stopt.lep1().Pt() ) > 10. )  return false;
  if ( ( stopt.isopf1() * stopt.lep1().Pt() ) > 5. )  return false; 

  if ( stopt.leptype() == 0 ) {

    //pass trigger if data - single electron
    if ( isData && stopt.ele27wp80() != 1 ) return false;
    //    if ( isData && stopt.trgel1() != 1 )  return false;
 
    //barrel only electrons
    if ( fabs(stopt.lep1().Eta() ) > 1.4442) return false;
    if ( stopt.eoverpin() > 4. ) return false;
 

  } else if ( stopt.leptype() == 1 ) {

    //pass trigger if data - single muon
    if ( isData && stopt.isomu24() != 1 ) return false;
    //    if ( isData && stopt.trgmu1() != 1 )  return false;
 
    if ( fabs(stopt.lep1().Eta() ) > 2.1)  return false;

  }

  return true;

}

//-------------------------------------------
// the dilepton selection
//-------------------------------------------

bool passDileptonSelection(bool isData) 
{
  //two lepton selection for 8 TeV 53 analysis

  //exactly 2 leptons
  if ( stopt.ngoodlep() != 2 ) return false;

  //opposite sign
  if ( stopt.id1()*stopt.id2()>0 ) return false;

  //pass trigger if data - dilepton
  if ( isData && stopt.mm() != 1 && stopt.me() != 1 
       && stopt.em() != 1 && stopt.ee() != 1 ) return false;

  //passes pt and eta requirements
  if ( stopt.lep1().Pt() < 20 )          return false;
  if ( stopt.lep2().Pt() < 20 )          return false;
  if ( fabs(stopt.lep1().Eta() ) > 2.4)  return false;
  if ( fabs(stopt.lep2().Eta() ) > 2.4)  return false;

  //consistency with pf leptons
  if ( fabs( stopt.pflep1().Pt() - stopt.lep1().Pt() ) > 10. )  return false;
  if ( fabs( stopt.pflep2().Pt() - stopt.lep2().Pt() ) > 10. )  return false;

  //information is only stored for leading lepton
  if ( ( stopt.isopf1() * stopt.lep1().Pt() ) > 5. )  return false; 
  if ( fabs(stopt.id1())==11 && stopt.eoverpin() > 4. ) return false;

  //barrel only electrons
  if (fabs(stopt.id1())==11 && fabs(stopt.lep1().Eta() ) > 1.4442) return false;
  if (fabs(stopt.id2())==11 && fabs(stopt.lep2().Eta() ) > 1.4442) return false;

  return true;

}

//-------------------------------------------
// lepton + isolated track selection CR5
//-------------------------------------------

bool passLepPlusIsoTrkSelection(bool isData) 
{
  //single lepton plus iso trk selection for 8 TeV 53 analysis

  //at least one lepton
  if ( !passSingleLeptonSelection(isData) ) return false;

  //pass isolated track requirement
  //unfortunately changed default value to 9999.
  if ( pfcandpt10() > 9990. || pfcandiso10() > 0.1 ) return false;

  return true;

}

//-------------------------------------------
// on-the-fly MET phi corrections
//-------------------------------------------

bool passLepPlusIsoTrkSelection_noEMu(bool isData, bool pickMuE) 
{
  //single lepton plus iso trk selection for 8 TeV 53 analysis

  //at least one lepton
  if ( !passSingleLeptonSelection(isData) ) return false;

  //pass isolated track requirement
  //unfortunately changed default value to 9999.
  //////  if ( pfcandpt10() > 9990. || pfcandiso10() > 0.1 ) return false;

  bool foundIsoTrack_noEMU=(stopt.pfcandpt10() <9998. && (abs(stopt.pfcandid10())!=13 && abs(stopt.pfcandid10())!=11 ) && stopt.pfcandiso10() < 0.1); 
  bool foundIsoTrack_EMU=(stopt.pfcandpt5() <9998. && (abs(stopt.pfcandid5())==13 || abs(stopt.pfcandid5())==11 ) && stopt.pfcandiso5() < 0.2); 

  ///  if(foundIsoTrack_EMU) cout << "EMU " << foundIsoTrack_EMU << endl;
  ///  if(foundIsoTrack_noEMU) cout << "noEMU " << foundIsoTrack_noEMU << endl;

  bool muE=(pickMuE && foundIsoTrack_EMU );
  bool noMuE=(!pickMuE && foundIsoTrack_noEMU);

  if(muE) return true;
  if(noMuE) return true;

  return false;

}


//-------------------------------------------
// on-the-fly MET phi corrections
//-------------------------------------------

pair<float,float> getPhiCorrMET( float met, float metphi, int nvtx, bool ismc){

  //using met phi corrections from C. Veelken (revision 1.6)
  //functions are available here:                                                                                                            
  //http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/CMSSW/JetMETCorrections/Type1MET/python/pfMETsysShiftCorrections_cfi.py                               

  float metx = met * cos( metphi );
  float mety = met * sin( metphi );

  float shiftx = 0.;
  float shifty = 0.;

  //use correction for data vs. mc                                                                                                                  
  shiftx = ismc ? (0.1166 + 0.0200*nvtx)
    : (0.2661 + 0.3217*nvtx);
  shifty = ismc ? (0.2764 - 0.1280*nvtx)
    : (-0.2251 - 0.1747*nvtx);

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

float dRbetweenVectors(LorentzVector vec1,LorentzVector vec2 ){ 

  float dphi = std::min(::fabs(vec1.Phi() - vec2.Phi()), 2 * M_PI - fabs(vec1.Phi() - vec2.Phi()));
  float deta = vec1.Eta() - vec2.Eta();

  return sqrt(dphi*dphi + deta*deta);
}

float getMinDphi(float metPhi, LorentzVector vec1, LorentzVector vec2 ) {

  float dphimj1_    = getdphi(metPhi, vec1.phi() );
  float dphimj2_    = getdphi(metPhi, vec2.phi() );
  float dphimjmin_  = TMath::Min( dphimj1_ , dphimj2_ );

  return dphimjmin_;

}


//---------------------------------------
// duplicates and filters
//---------------------------------------

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

bool is_duplicate (const DorkyEventIdentifier &id, std::set<DorkyEventIdentifier> &already_seen) {
  std::pair<std::set<DorkyEventIdentifier>::const_iterator, bool> ret =
    already_seen.insert(id);
  return !ret.second;
}

//--------------------------------------------------------------------

int load_badlaserevents  (char* filename, std::set<DorkyEventIdentifier> &events_lasercalib) {

  ifstream in;
  in.open(filename);

   unsigned long int run, event, lumi;
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

//--------------------------------------------------------------------

bool is_badLaserEvent (const DorkyEventIdentifier &id, std::set<DorkyEventIdentifier> &events_lasercalib) {
  if (events_lasercalib.find(id) != events_lasercalib.end()) return true;
  return false;
}
