#include "stopUtils.h"

#include "Math/LorentzVector.h"
#include "Math/VectorUtil.h"
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

#include "TFitter.h"
#include "MT2Utility.h"
#include "mt2bl_bisect.h"
#include "mt2w_bisect.h"

using namespace Stop;
using namespace std;

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > LorentzVector;

unsigned int getNJets(const float etacut){

  unsigned int njets=0;

  for ( unsigned int i=0; i<stopt.pfjets().size() ; i++) {
    
    if( stopt.pfjets().at(i).pt()<30 )  continue;
    if( fabs(stopt.pfjets().at(i).eta())>etacut )  continue;
    if ( (fabs(stopt.pfjets().at(i).eta()) < 2.5) 
	&& (stopt.pfjets_beta2_0p5().at(i)<0.2) ) continue;

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
    if ( stopt.pfjets_beta2_0p5().at(i) < 0.2 ) continue;
    if( stopt.pfjets_csv().at(i)      < discr   ) continue;

    // skip these indices                                                                                                                                                                                  
    if( int(i) == iskip1 ) continue;
    if( int(i) == iskip2 ) continue;

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
  if (fabs(eta) >=0.0 && fabs(eta) < 0.5) return 1.052;
  if (fabs(eta) >=0.5 && fabs(eta) < 1.1) return 1.057;
  if (fabs(eta) >=1.1 && fabs(eta) < 1.7) return 1.096;
  if (fabs(eta) >=1.7 && fabs(eta) < 2.3) return 1.134;
  if (fabs(eta) >=2.3 && fabs(eta) < 5.0) return 1.288;
  return 1.0;
}

//------------------------------------------------------------------------------------------------
// Fix function to partial correction used in looper before 1.16
//------------------------------------------------------------------------------------------------

float getDataMCRatioFix(float eta){
  if ( eta > 0 ) 
	return 1.0;
  else
     	return getDataMCRatio(eta);
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

  if (!name.Contains("T2")  && !name.Contains("TChiwh")) {
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
// the isolated track veto used at HCP
//-------------------------------------------

bool passTauVeto() {

  if(stopt.pfTau_leadPtcandID()!=(-1)) return false;

  return true;

}

//-------------------------------------------
// the isolated track veto used at HCP
//-------------------------------------------

bool passIsoTrkVeto() 
{

  //pass isolated track veto
  //unfortunately changed default value to 9999.
  if ( stopt.pfcandpt10() <9998. && stopt.pfcandiso10() < 0.1 ) return false;

  return true;

}

//-------------------------------------------
// the isolated track veto looser E+mu
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
// the isolated track veto looser E+mu and OSTrack
//-------------------------------------------

bool passIsoTrkVeto_v3() 
{

  //pass isolated track veto
  //unfortunately changed default value to 9999.
  if ( stopt.pfcandptOS10() <9998. && abs(stopt.pfcandidOS10())!=13 && abs(stopt.pfcandidOS10())!=11 && stopt.pfcandisoOS10() < 0.1 ) return false;
  if ( stopt.pfcandpt5()  <9998. && abs(stopt.pfcandid5())==13 && stopt.pfcandiso5() < 0.2) return false;
  if ( stopt.pfcandpt5()  <9998. && abs(stopt.pfcandid5())==11 && stopt.pfcandiso5() < 0.2) return false;

  return true;

}


//-------------------------------------------
// the isolated track veto with looser dz
//-------------------------------------------

bool passIsoTrkVeto_v4()
{

  //pass isolated track veto
  //unfortunately changed default value to 9999.
  // We want to check for the generic track only there is now good e/mu candidate
  if ( stopt.pfcandptOS10looseZ() <9998. && abs(stopt.pfcandid5looseZ())!=13 && abs(stopt.pfcandid5looseZ())!=11 && stopt.pfcandisoOS10looseZ() < 0.1 ) return false;

  if ( stopt.pfcandpt5looseZ()  <9998. && abs(stopt.pfcandid5looseZ())==13 && stopt.pfcandiso5looseZ() < 0.2) return false;
  if ( stopt.pfcandpt5looseZ()  <9998. && abs(stopt.pfcandid5looseZ())==11 && stopt.pfcandiso5looseZ() < 0.2) return false;

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
// good lepton + veto isolated track
//-------------------------------------------

bool passOneLeptonSelection( bool isData) 
{
  //single lepton selection for 8 TeV 53 analysis
  if ( !passSingleLeptonSelection(isData) ) return false;

  //pass isolated track veto
  //unfortunately changed default value to 9999.
  //  if ( stopt.pfcandpt10() <9998. && stopt.pfcandiso10() < 0.1 ) return false;
  if ( !passIsoTrkVeto_v4() ) return false;

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
  //  if ( passIsoTrkVeto_v4() ) return false;

  return true;

}

bool pass_T2tt_LM(bool isData){

  if ( !passSingleLeptonSelection(isData) ) return false;

  if ( !passIsoTrkVeto_v4() ) return false;

  // this is just a preselection
  if(stopt.t1metphicorr()<100) return false;

  if(stopt.t1metphicorrmt()<120) return false;

  vector<LorentzVector> myJets;
  vector<float> myJetsTag;
  vector<int> myJetsMC;
  vector<float> myJetsSigma;
  int btag=0;

  for (int ijet =0; ijet<(int)stopt.pfjets().size(); ijet++){

    if ( stopt.pfjets().at(ijet).Pt() < 30 ) continue;
    if ( fabs(stopt.pfjets().at(ijet).eta()) > 2.4 ) continue;
    // later add the MVA
    if ( stopt.pfjets_beta2_0p5().at(ijet) < 0.2 ) continue;

    myJets.push_back(stopt.pfjets().at(ijet));
    myJetsTag.push_back(stopt.pfjets_csv().at(ijet));
    if(stopt.pfjets_csv().at(ijet) > 0.679) btag++;

    myJetsSigma.push_back(stopt.pfjets_sigma().at(ijet));

  }

  if(myJets.size()<4) return false;
  if(btag==0) return false;

  // mindPhi
  float dphimjmin=getMinDphi(stopt.t1metphicorrphi(), myJets.at(0),myJets.at(1));
  if(dphimjmin<0.8) return false;

  // chi2
  double chi2 = calculateChi2SNT(myJets, myJetsSigma, myJetsTag);
  if(chi2>5) return false;

  // mt2w
  //  double x_mt2w = calculateMT2w(myJets, myJetsTag, stopt.lep1(), stopt.t1metphicorr(), stopt.t1metphicorrphi());
  return true;

}


bool pass_T2tt_HM(bool isData){

  if ( !passSingleLeptonSelection(isData) ) return false;

  // this is temporary waiting for the new babies                                                                                                                                    
  if ( !passIsoTrkVeto_v4() ) return false;

  // this is just a preselection                                                                                                                                                     
  if(stopt.t1metphicorr()<100) return false;

  if(stopt.t1metphicorrmt()<120) return false;

  vector<LorentzVector> myJets;
  vector<float> myJetsTag;
  vector<int> myJetsMC;
  vector<float> myJetsSigma;
  int btag;

  for (int ijet =0; ijet<(int)stopt.pfjets().size(); ijet++){

      if ( stopt.pfjets().at(ijet).Pt() < 30 ) continue;
      if ( fabs(stopt.pfjets().at(ijet).eta()) > 2.4 ) continue;
      // later add the MVA
      if ( stopt.pfjets_beta2_0p5().at(ijet) < 0.2 ) continue;

      myJets.push_back(stopt.pfjets().at(ijet));
      myJetsTag.push_back(stopt.pfjets_csv().at(ijet));
      if(stopt.pfjets_csv().at(ijet) > 0.679) btag++;
      myJetsSigma.push_back(stopt.pfjets_sigma().at(ijet));
  }

  if(myJets.size()<4) return false;
  if(btag==0) return false;

  // mindPhi                                                                                                                                                                         
  float dphimjmin=getMinDphi(stopt.t1metphicorrphi(), myJets.at(0),myJets.at(1));
  if(dphimjmin<0.8) return false;

  // chi2                                                                                                                                                                            
  double chi2 = calculateChi2SNT(myJets, myJetsSigma, myJetsTag);

  if(chi2>5) return false;

  // mt2w                                                                                                                                                                            
  double x_mt2w = calculateMT2w(myJets, myJetsTag, stopt.lep1(), stopt.t1metphicorr(), stopt.t1metphicorrphi());

  if(x_mt2w<175) return false;

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
  bool foundIsoTrack_noEMU=(stopt.pfcandptOS10looseZ() <9998. && (abs(stopt.pfcandidOS10looseZ())!=13 && abs(stopt.pfcandidOS10looseZ())!=11 ) && stopt.pfcandisoOS10looseZ() < 0.1); 
  bool foundIsoTrack_EMU= (stopt.pfcandpt5looseZ() <9998. && (abs(stopt.pfcandid5looseZ())==13 || abs(stopt.pfcandid5looseZ())==11 ) && stopt.pfcandiso5looseZ() < 0.2);

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

float dRbetweenVectors(LorentzVector& vec1,LorentzVector& vec2 ){ 

  float dphi = std::min(::fabs(vec1.Phi() - vec2.Phi()), 2 * M_PI - fabs(vec1.Phi() - vec2.Phi()));
  float deta = vec1.Eta() - vec2.Eta();

  return sqrt(dphi*dphi + deta*deta);
}

float getMinDphi(float metPhi, LorentzVector& vec1, LorentzVector& vec2 ) {

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


//--------------------------------------------------------------------
double calculateMT2w(vector<LorentzVector>& jets, vector<float>& btag, LorentzVector& lep, float met, float metphi, MT2Type mt2type){

	// I am asumming that jets is sorted by Pt
	assert ( jets.size() == btag.size() );
	// require at least 2 jets
	if ( jets.size()<2 ) return 99999.; 

	// First we count the number of b-tagged jets, and separate those non b-tagged
	std::vector<int> bjets;
	std::vector<int> non_bjets;
	for( unsigned int i = 0 ; i < jets.size() ; i++ ){
	  if( btag.at(i) > BTAG_MED ) {
	    bjets.push_back(i);
	  } else {
	    non_bjets.push_back(i);
	  }
	}	

	int n_btag = (int) bjets.size();
	//	cout << "n_btag = " << n_btag << endl;

	// We do different things depending on the number of b-tagged jets
	// arXiv:1203.4813 recipe

	int nMax=-1;
	if(jets.size()<=3) nMax=non_bjets.size();
	else nMax=3;

	if (n_btag == 0){
	  // If no b-jets select the minimum of the mt2w from all combinations with 
	  // the three leading jets
	  float min_mt2w = 9999;

	  for (int i=0; i<nMax; i++)
	    for (int j=0; j<nMax; j++){
	      if (i == j) continue;
	      float c_mt2w = mt2wWrapper(lep, 
					 jets[non_bjets[i]],
					 jets[non_bjets[j]], met, metphi, mt2type);
	      if (c_mt2w < min_mt2w)
		min_mt2w = c_mt2w;
	    }
	  return min_mt2w;
	} else if (n_btag == 1 ){
	  // if only one b-jet choose the three non-b leading jets and choose the smaller
	  float min_mt2w = 9999;

	  for (int i=0; i<nMax; i++){
	    float c_mt2w = mt2wWrapper(lep, jets[bjets[0]], jets[non_bjets[i]], met, metphi, mt2type);
	    if (c_mt2w < min_mt2w)
	      min_mt2w = c_mt2w;
	  }
	  for (int i=0; i<nMax; i++){
	    float c_mt2w = mt2wWrapper(lep, jets[non_bjets[i]], jets[bjets[0]], met, metphi, mt2type);
	    if (c_mt2w < min_mt2w)
	      min_mt2w = c_mt2w;
	  }
	  return min_mt2w;
	} else if (n_btag >= 2) {
	  // if 3 or more b-jets the paper says ignore b-tag and do like 0-bjets 
	  // but we are going to make the combinations with the b-jets
	  float min_mt2w = 9999;
	  for (int i=0; i<n_btag; i++)
	    for (int j=0; j<n_btag; j++){
	      if (i == j) continue;
	      float c_mt2w = mt2wWrapper(lep, 
					 jets[bjets[i]],
					 jets[bjets[j]], met, metphi, mt2type);
	      if (c_mt2w < min_mt2w)
		min_mt2w = c_mt2w;
	    }
	  return min_mt2w;
	}

	return -1.;
}


//---------------------------------------------------------------------


// This funcion is a wrapper for mt2w_bisect etc that takes LorentzVectors instead of doubles
double mt2wWrapper(LorentzVector& lep, LorentzVector& jet_o, LorentzVector& jet_b, float met, float metphi, MT2Type mt2type){

	// same for all MT2x variables
	float metx = met * cos( metphi );
	float mety = met * sin( metphi );

	double pl[4];     // Visible lepton
	double pb1[4];    // bottom on the same side as the visible lepton
	double pb2[4];    // other bottom, paired with the invisible W
	double pmiss[3];  // <unused>, pmx, pmy   missing pT
	pl[0]= lep.E(); pl[1]= lep.Px(); pl[2]= lep.Py(); pl[3]= lep.Pz();
	pb1[1] = jet_o.Px();  pb1[2] = jet_o.Py();   pb1[3] = jet_o.Pz();
	pb2[1] = jet_b.Px();  pb2[2] = jet_b.Py();   pb2[3] = jet_b.Pz();
	pmiss[0] = 0.; pmiss[1] = metx; pmiss[2] = mety;

	// specifics for each variable
	if (mt2type == MT2b) {
	  double pmiss_lep[3];
	  pmiss_lep[0] = 0.;
	  pmiss_lep[1] = pmiss[1]+pl[1]; pmiss_lep[2] = pmiss[2]+pl[2];

	  pb1[0] = jet_o.mass();
	  pb2[0] = jet_b.mass();

	  mt2_bisect::mt2 mt2_event;
	  mt2_event.set_momenta( pb1, pb2, pmiss_lep );
	  mt2_event.set_mn( 80.385 );   // Invisible particle mass == W mass
	  return mt2_event.get_mt2();
	}

	else {
	  pb1[0] = jet_o.E();
	  pb2[0] = jet_b.E();

	  if (mt2type == MT2bl) {
	    mt2bl_bisect::mt2bl mt2bl_event;
	    mt2bl_event.set_momenta(pl, pb1, pb2, pmiss);
	    return mt2bl_event.get_mt2bl();
	  }
	  else if (mt2type == MT2w) {
	    mt2w_bisect::mt2w mt2w_event;
	    mt2w_event.set_momenta(pl, pb1, pb2, pmiss);
	    return mt2w_event.get_mt2w();
	  }
	}

	// SHOULDN'T GET HERE
	std::cout << "ERROR: " << __FILE__ << " " << __LINE__ << ": mt2wWrapper: shouldn't get here! mt2type == " 
		  << mt2type << std::endl;
	return -1.;
}


//--------------------------------------------------------------------
double fc2 (double c1, double m12, double m22, double m02, bool verbose = false)
{
  if (verbose) {
    printf("c1: %4.2f\n", c1);
    printf("m12: %4.2f\n", m12);
    printf("m22: %4.2f\n", m22);
    printf("m02: %4.2f\n", m02);
  }

  double a = m22;
  double b = (m02 - m12 - m22) * c1;
  double c = m12 * c1 * c1 - PDG_W_MASS * PDG_W_MASS;

  if (verbose) {
    printf("a: %4.2f\n", a);
    printf("b: %4.2f\n", b);
    printf("c: %4.2f\n", c);
  }

  double num = -1. * b + sqrt(b * b - 4 * a * c);
  double den = 2 * a;

  if (verbose) {
    printf("num: %4.2f\n", num);
    printf("den: %4.2f\n", den);
    printf("num/den: %4.2f\n", num/den);
  }

  return (num/den);
}

//--------------------------------------------------------------------
double fchi2 (double c1, double pt1, double sigma1, double pt2, double sigma2,
              double m12, double m22, double m02){
  double rat1 = pt1 * (1 - c1) / sigma1;
  double rat2 = pt2 * (1 - fc2(c1, m12, m22, m02)) / sigma2;

  return ( rat1 * rat1 + rat2 * rat2);
}

//--------------------------------------------------------------------
void minuitFunction(int&, double* , double &result, double par[], int){
  result=fchi2(par[0], par[1], par[2], par[3], par[4], par[5], par[6], par[7]);
}

// This function calculates the hadronic chi2 - atlas version
double calculateChi2(vector<LorentzVector>& jets, vector<float>& sigma_jets){

	assert(jets.size() == sigma_jets.size());

	int n_jets = jets.size();

	int j1=-1;
	int j2=-1;
	int bi=-1;
	double min_dR = 9999;
	for (int i=0; i < n_jets; i++)
		for (int j=0; j < n_jets; j++){
			double dR =  ROOT::Math::VectorUtil::DeltaR(jets.at(i),jets.at(j));
			LorentzVector hadW = jets.at(i) + jets.at(j);
			if ( dR < min_dR && hadW.mass() > 60){
				min_dR = dR;
				j1 = i;
				j2 = j;
			}
		}
	if ( j1 < 0 || j2 < 0 ) return -1;
	
	LorentzVector hadW = jets.at(j1) + jets.at(j2);

	min_dR = 9999;	
	for (int j=0; j < n_jets; j++){
		double dR =  ROOT::Math::VectorUtil::DeltaR(hadW,jets.at(j));
		LorentzVector hadT = hadW + jets.at(j);
		if ( dR < min_dR && hadT.mass() > 130 && hadT.mass()< 205 ){
			min_dR = dR;
			bi = j;
		}
	}

	if ( bi <  0 ) return -1;

/*
	LorentzVector c1_direction = jets.at(0)/jets.at(0).P();
	LorentzVector c2_direction = jets.at(1)/jets.at(1).P();

	int cluster1[15];
	int cluster2[15];
	int n_cluster1, n_cluster2;

	float diff2;
	do {
		n_cluster1 = 0;
		n_cluster2 = 0;
		for(int i=0; i < n_jets; i++){
			double dRc1 = (c1_direction - jets.at(i)).P();
			double dRc2 = (c2_direction - jets.at(i)).P();
//			cout << "RR: " << dRc1 << "  " << dRc2 <<  "   " << jets.at(i) << endl;
			if ( dRc1 < dRc2 )
				cluster1[n_cluster1++] = i;
			else
				cluster2[n_cluster2++] = i;
		}

		LorentzVector c1_new;
		for ( int i=0; i<n_cluster1; i++)
			c1_new += jets.at(cluster1[i]);
		c1_new = c1_new/c1_new.P();

		LorentzVector c2_new;
		for ( int i=0; i<n_cluster2; i++)
			c2_new += jets.at(cluster2[i]);
		c2_new = c2_new/c2_new.P();

//		cout << c1_new << "    " << c2_new << endl;

		diff2 = ( c1_new - c1_direction ).P2() + ( c2_new - c2_direction ).P2();
//		cout << n_cluster1 << "   " << n_cluster2 << "   " << diff2 << endl;

		c1_direction = c1_new;
		c2_direction = c2_new;

		assert ( n_cluster1 > 0 ); 
		assert ( n_cluster2 > 0 ); 
	} while ( diff2 > 0.001 );

	cout << n_cluster1 << "   " << n_cluster2 << "   " << diff2 << endl;

	

	double min_chi2 = 9999;
	for ( int b=0; b<n_jets; ++b )
		for (int j1=0; j1<n_jets; ++j1)
			for (int j2=0; j2<n_jets; ++j2){
				if (b == j1 || b == j2 || j1 == j2) continue;
				
				double c_chi2 = getChi2(jets.at(b), jets.at(j1), jets.at(j2),
						sigma_jets.at(b), sigma_jets.at(j1), sigma_jets.at(j2));
				
				if (c_chi2 < min_chi2) min_chi2 = c_chi2;

			}
*/

	double c_chi2 =  getChi2(jets.at(bi), jets.at(j1), jets.at(j2), 
		sigma_jets.at(bi), sigma_jets.at(j1), sigma_jets.at(j2) ); 
	return c_chi2;
}

// This function calculates the hadronic chi2 - SNT version
double calculateChi2SNT(vector<LorentzVector>& jets, vector<float>& sigma_jets, vector<float>& btag){

  assert(jets.size() == sigma_jets.size());

  //check at most first 6 jets
  int n_jets = jets.size();
  if (n_jets>6) n_jets = 6;
  //consider at least 3 jets
  if (n_jets<3) return 99999.;
  
  vector<int> v_i, v_j;
  vector<double> v_k1, v_k2;
  for ( int i=0; i<n_jets; ++i )
    for ( int j=i+1; j<n_jets; ++j ){

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
      minimizer->SetFCN(minuitFunction);
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
        cout<<"[PartonCombinatorics::recoHadronicTop] ERROR: c1 parameter is NAN! Skipping this parton combination"
	    <<endl;
        continue;
      }
      double c2 = fc2(c1, jets[i].mass2(), jets[j].mass2(), hadW.mass2());

      delete minimizer;


      //     * W Mass check :)
      //     *  Never trust a computer you can't throw out a window.
      //      *  - Steve Wozniak

      // cout << "c1 = " <<  c1 << "  c1 = " << c2 << "   M_jj = "
      // 	   << ((jets[i] * c1) + (jets[j] * c2)).mass() << endl;

      v_i.push_back(i);
      v_j.push_back(j);
      v_k1.push_back(c1);
      v_k2.push_back(c2);
    }
  
  //Apply b-consistency requirement
  int n_btag = 0;
  for( int i = 0 ; i < n_jets ; i++ )
    if( btag.at(i) > BTAG_MED ) n_btag++;

  double chi2min = 99999.;

  //consider b-jet in leading 3 jets
  for ( int b=0; b<n_jets; ++b ) {    

    //if not tagged, consider only 3 leading jets
    if( btag.at(b)<BTAG_MED && b>2 ) continue;

    //require b-tagging if have more than 1 b-tag
    if( n_btag>1 && btag.at(b) < BTAG_MED ) continue;
    double pt_b = jets[b].Pt();
      
    for (unsigned int w = 0; w < v_i.size() ; ++w ) {
      int i = v_i[w];
      int j = v_j[w];
      if ( i==b || j==b ) continue;
      //count number of b-tagged Ws
      int nwb = 0;
      if (btag.at(i) > BTAG_MED) nwb++;
      if (btag.at(j) > BTAG_MED) nwb++;
      //no btagged jets in W if have few btags
      if ( n_btag<3  && nwb>0 ) continue;
      //In 3 b-tag case, allow for 1 W jet to be tagged
      // If have more b-tags then btagging information not useful
      if ( n_btag==3 && nwb>1 ) continue;
 
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
      double sigma_w2 = pow(pt_w1*sigma_jets[i], 2)
	+ pow(pt_w2*sigma_jets[j], 2);
      double smw2 = (1. + 2.*pow(pt_w,2)/pow(massW,2))*sigma_w2;
      double pt_t = hadT.Pt();
      double sigma_t2 = pow(c1*pt_w1*sigma_jets[i],2)
	+ pow(c2*pt_w2*sigma_jets[j],2)
	+ pow(pt_b*sigma_jets[b],2);
      double smtop2 = (1. + 2.*pow(pt_t,2)/pow(massT,2))*sigma_t2;
      
      double c_chi2 = pow(massT-PDG_TOP_MASS, 2)/smtop2
	+ pow(massW-PDG_W_MASS, 2)/smw2;
      if (c_chi2<chi2min) chi2min = c_chi2;

    }
  }

  return chi2min;
}


double getChi2(LorentzVector& jets_b, LorentzVector& jets_j1, LorentzVector& jets_j2,
		float sigma_b, float sigma_j1, float sigma_j2){

	LorentzVector hadW = jets_j1 + jets_j2;
	LorentzVector hadT = jets_b + hadW;

	double massT = hadT.mass();
	double massW = hadW.mass();

	double pt_w1 = jets_j1.Pt();
	double pt_w2 = jets_j2.Pt();
	double pt_b = jets_b.Pt();

	double pt_w = hadW.Pt();

	double sigma_w2 = pt_w1*sigma_j1 * pt_w1*sigma_j1
		+ pt_w2*sigma_j2 * pt_w2*sigma_j2;

	double smw2 = (1.+2.*pt_w*pt_w/massW/massW)*sigma_w2;

	double pt_t = hadT.Pt();

	double sigma_t2 = pt_w1*sigma_j1 * pt_w1*sigma_j1
		+ pt_w2*sigma_j2 * pt_w2*sigma_j2
		+ pt_b*sigma_b * pt_b*sigma_b;

	double smtop2 = (1.+2.*pt_t*pt_t/massT/massT)*sigma_t2;

	double c_chi2 = (massT-PDG_TOP_MASS)*(massT-PDG_TOP_MASS)/smtop2
		+ (massW-PDG_W_MASS)*(massW-PDG_W_MASS)/smw2;

	return c_chi2;
}


