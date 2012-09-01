#include <algorithm>
#include <iostream>
#include <map>
#include <vector>
#include <set>
#include <sstream>
#include "TChain.h"
#include "TChainElement.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TProfile.h"
#include "TTree.h"
#include "TVector2.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"
#include "TRandom3.h"
#include "Math/LorentzVector.h"
#include "histtools.h"
#include "singleLeptonLooper.h"
#include "TTreeCache.h"
#include "TDatabasePDG.h"

#include "../CORE/CMS2.h"
#include "../CORE/utilities.h"
#include "../CORE/ssSelections.h"
#include "../CORE/electronSelections.h"
#include "../CORE/electronSelectionsParameters.h"
#include "../CORE/MITConversionUtilities.h"
#include "../CORE/muonSelections.h"
#include "../CORE/eventSelections.h"
#include "../CORE/trackSelections.h"
#include "../CORE/metSelections.h"
#include "../CORE/jetcorr/FactorizedJetCorrector.h"
#include "../CORE/jetcorr/JetCorrectionUncertainty.h"
#include "../CORE/jetSelections.h"
#include "../CORE/photonSelections.h"
#include "../CORE/triggerUtils.h"
#include "../CORE/triggerSuperModel.h"
#include "../CORE/mcSelections.h"
#include "../CORE/susySelections.h"
#include "../CORE/mcSUSYkfactor.h"
#include "../CORE/SimpleFakeRate.h"
#include "../Tools/goodrun.h"
#include "../Tools/vtxreweight.h"
#include "../Tools/msugraCrossSection.h"
#include "BtagFuncs.h"
//#include "../Tools/bTagEff_BTV.h"

bool verbose              = false;
bool doTenPercent         = false;
bool vetoTransition       = true;
bool useOldIsolation      = false;

//#include "../CORE/topmass/getTopMassEstimate.icc" // REPLACETOPMASS
//#include "../CORE/triggerUtils.cc"

using namespace std;
using namespace tas;

typedef vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > VofP4;

//mSUGRA scan parameters-----------------------------

const int   nm0points    = 100;
const float m0min        = 20.;
const float m0max        = 2020.;
const int   nm12points   = 38;
const float m12min       = 20.;
const float m12max       = 780.;

//---------------------------------------------------

void fillUnderOverFlow(TH1F *h1, float value, float weight = 1.);
void fillUnderOverFlow(TH2F *h2, float xvalue, float yvalue, float weight = 1.);
//void fillUnderOverFlow(TProfile *h2, float xvalue, float yvalue);
void fillOverFlow(TH1F *h1, float value, float weight = 1.);
void fillOverFlow(TH2F *h2, float xvalue, float yvalue, float weight = 1.);
void fillHistos(TH1F *h1[4][4],float value, float weight, int myType, int nJetsIdx);
void fillHistos(TH2F *h2[4][4],float xvalue, float yvalue, float weight, int myType, int nJetsIdx);
void fillHistos(TProfile *h2[4][4],float xvalue, float yvalue,  int myType, int nJetsIdx);

//--------------------------------------------------------------------

void singleLeptonLooper::weight3D_init( std::string WeightFileName ) { 

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

//--------------------------------------------------------------------

double singleLeptonLooper::weight3D( int pv1, int pv2, int pv3 ) {

  using std::min;

  int npm1 = min(pv1,49);
  int np0 = min(pv2,49);
  int npp1 = min(pv3,49);

  return Weight3D[npm1][np0][npp1];

}

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
  cout << "eta24     " <<  ( TMath::Abs(cms2.mus_p4()[muidx].eta()) < 2.4)                         << endl;
  cout << "chi2/ndf  " <<  ( cms2.mus_gfit_chi2().at(muidx)/cms2.mus_gfit_ndof().at(muidx) < 10)   << endl;
  cout << "global    " <<  ( ((cms2.mus_type().at(muidx)) & (1<<1)) != 0)                          << endl;
  cout << "tracker   " <<  ( ((cms2.mus_type().at(muidx)) & (1<<2)) != 0)                          << endl;
  cout << "nhits     " <<  ( cms2.mus_validHits().at(muidx) > 10)                                  << endl;
  cout << "stahits   " <<  ( cms2.mus_gfit_validSTAHits().at(muidx) != 0)                          << endl;
  cout << "d0PV      " <<  ( TMath::Abs(mud0PV_smurfV3(muidx)) < 0.02)                             << endl;
  cout << "dzPV      " <<  ( TMath::Abs(mudzPV_smurfV3(muidx)) < 1  )                              << endl;
  cout << "dpt/pt    " <<  ( cms2.mus_ptErr().at(muidx)/cms2.mus_p4().at(muidx).pt()<0.1)          << endl;

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
    for (unsigned int igen = 0; igen < cms2.genps_p4().size(); igen++) {
      
      int id = cms2.genps_id().at(igen);
      if( abs(id)!=5 ) continue;
      
      if( dRbetweenVectors( p4 , cms2.genps_p4().at(igen) ) < dR ) return true;

    } 
    return false;
}

//--------------------------------------------------------------------

int isGenQGMatched ( LorentzVector p4, float dR ) {
  //Start from the end that seems to have the decay products of the W first
  for (int igen = (cms2.genps_p4().size()-1); igen >-1; igen--) {
    float deltaR = dRbetweenVectors( p4 , cms2.genps_p4().at(igen) );
    if ( deltaR > dR ) continue;
    int id = abs(cms2.genps_id().at(igen));
    int mothid = abs(genps_id_mother().at(igen));
    // cout<<"status 3 particle ID "<<id<<" mother "<<mothid
    // 	<<" dR to jet "<<deltaR<<endl;
    if (id<6 && mothid==24) return 2;
    if (id==5 && mothid==6) return 1;
    if (id==21) return 3;
    if (id<6) return 4;
  }
  return -1;
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

//--------------------------------------------------------------------

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

//--------------------------------------------------------------------

bool objectPassTrigger(const LorentzVector &obj, char* trigname, float drmax = 0.1 ){

  TString exact_trigname = triggerName( trigname );

  if( exact_trigname.Contains("TRIGGER_NOT_FOUND") ){
    cout << __FILE__ << " " << __LINE__ << " Error! couldn't find trigger name " << trigname << endl;
    return false;
  }

  std::vector<LorentzVector> trigp4 = cms2.hlt_trigObjs_p4()[findTriggerIndex(exact_trigname)];

  if( trigp4.size() == 0 ) return false;

  for (unsigned int i = 0; i < trigp4.size(); ++i){
    float dr = dRbetweenVectors(trigp4[i], obj);
    if( dr < drmax ) return true;
  }

  return false;
}

//--------------------------------------------------------------------

singleLeptonLooper::singleLeptonLooper()
{
  g_susybaseline = false;
  g_createTree   = false;
  g_useBitMask   = false;
  random3_ = new TRandom3(1);
  initialized = false;
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

int getIndexFromM0(float m0){
  
  float binsize = (m0max - m0min) / (float) nm0points;
  int index     = (int)((m0 - m0min) / binsize);
  return index;
    
}

//--------------------------------------------------------------------

int getIndexFromM12(float m12){
  
  float binsize = (m12max - m12min) / (float) nm12points;
  int index     = (int)((m12 - m12min) / binsize);
  return index;
    
}

//--------------------------------------------------------------------

void singleLeptonLooper::InitBaby(){

  //pdf variables
  pdfid1_ = -999;
  pdfid2_ = -999;
  pdfQ_   = -99999.;
  pdfx1_  = -99999.;
  pdfx2_  = -99999.;

  mutrigweight_ = 1.;

  jetid_	= 1;
  jetid30_	= 1;

  // btag variables
  nbtagsssv_	= 0;
  nbtagstcl_	= 0;
  nbtagstcm_	= 0;
  nbtagscsvl_   = 0;
  nbtagscsvm_   = 0;
  nbtagscsvt_   = 0;
  nbtagsssvcorr_    = 0;
  nbtagstclcorr_    = 0;
  nbtagstcmcorr_    = 0;
  nbtagscsvlcorr_   = 0;
  nbtagscsvmcorr_   = 0;
  nbtagscsvtcorr_   = 0;

  // njets with JEC variation
  njetsUp_	= 0;
  njetsDown_	= 0;
  ht_           = -999.;
  htUp_		= -999.;
  htDown_	= -999.;

  // Type1 pfmet
  t1met10_	=-999.;
  t1met20_	=-999.;
  t1met30_	=-999.;
  t1met10phi_	=-999.;
  t1met20phi_	=-999.;
  t1met30phi_	=-999.;
  t1met10mt_	=-999.;
  t1met20mt_	=-999.;
  t1met30mt_	=-999.;
  lepmetpt_	=-999.;
  lept1met10pt_	=-999.;
  
  //phi corrected type1 mets
  t1metphicorr_	      =-999.;
  t1metphicorrphi_    =-999.;
  t1metphicorrlep_    =-999.;
  t1metphicorrlepphi_ =-999.;
  t1metphicorrmt_     =-999.;
  t1metphicorrlepmt_  =-999.;
  
  // pfjet vars
  npfjets30_	= 0;
  npfjets35_	= 0;
  npfjets40_	= 0;
  npfjets45_	= 0;
  npfjets30lepcorr_ = 0;
  knjets_       = 1.;

  htpf30_	= 0.;
  htpf35_	= 0.;
  htpf40_	= 0.;
  htpf45_	= 0.;
  htpfres30_	= 0.;
  htpfres35_	= 0.;
  htpfres40_	= 0.;
  htpfres45_	= 0.;

  //iso trk vars
  trkpt5_ 	    = -999.;
  mleptrk5_ 	    = -999.;
  trkreliso5_ 	    = -999.;
  trkpt10_ 	    = -999.;
  mleptrk10_ 	    = -999.;
  trkreliso10_ 	    = -999.;
  trkpt5loose_ 	    = -999.;
  trkreliso5loose_  = -999.;
  trkpt10loose_     = -999.;
  trkreliso10loose_ = -999.;

  // MC truth info
  mcid1_	= -1;
  mcid2_	= -1;
  mclep1_	=  0;
  mclep2_	=  0;
  mctaud1_      =  0;
  mctaud2_      =  0;
  mctaudvis1_   =  0;
  mctaudvis2_   =  0;
  mcdecay1_	= -1;
  mcdecay2_	= -1;
  mcdr1_	= -1;
  mcdr2_	= -1;
  
  mlepid_       = -1;
  mlep_         =  0;
  mleppassid_   = -1;
  mleppassiso_  = -1;
  mlepiso_      = -1.0;

  mllgen_	= -1;
  pthat_	= -1;
  qscale_	= -1;
  genmet_       = -9999;
  gensumet_     = -9999;
  genmetphi_    = -9999;   
  m0_		= -9999;
  m12_		= -9999;
  mG_		= -9999;
  mL_		= -9999;
  x_		= -9999;
  ksusy_	= -999;
  ksusyup_	= -999;
  ksusydn_	= -999;
  xsecsusy_	= -999;
  xsecsusy2_	= -999;
  w1_		= -999;
  w2_		= -999;
  acc_2010_	= -999;
  acc_highmet_  = -999;
  acc_highht_	= -999;
  nels_		= -1;
  nmus_		= -1;
  ntaus_	= -1;
  nleps_	= -1;
  ptjetraw_	= -9999.;
  ptjet23_	= -9999.;
  ptjetF23_	= -9999.;
  ptjetO23_	= -9999.;
  cosphijz_	= -9999.;
  dilep_	= 0;
  jet_		= 0;

  // jet p4's
  pfjet1_	= 0;
  pfjet2_	= 0;
  pfjet3_	= 0;
  pfjet4_	= 0;
  pfjet5_	= 0;
  pfjet6_	= 0;

  bjet1_      = -1; 
  bjet2_      = -1; 
  bjet3_      = -1; 
  bjet4_      = -1; 
  bjet5_      = -1; 
  bjet6_      = -1; 
  lepjet1_    = -1; 
  lepjet2_    = -1; 
  lepjet3_    = -1; 
  lepjet4_    = -1; 
  lepjet5_    = -1; 
  lepjet6_    = -1; 
  qgjet1_     = -1; 
  qgjet2_     = -1; 
  qgjet3_     = -1; 
  qgjet4_     = -1; 
  qgjet5_     = -1; 
  qgjet6_     = -1; 

  lep1chi2ndf_	= -9999.;
  lep2chi2ndf_	= -9999.;
  lep1dpt_	= -9999.;
  lep2dpt_	= -9999.;
  leptype1_	= -999;
  leptype2_	= -999;
  lep1_		= 0;
  lep2_		= 0;
  trklep1_	= 0;
  trklep2_	= 0;
  gfitlep1_	= 0;
  gfitlep2_	= 0;
  lepp_		= 0;
  lepm_		= 0;
  pflep1_	= 0;
  pflep2_	= 0;
  leppfjet1_	= 0;
  leppfjet2_	= 0;
  pflep_        = 0;
  pftaud_       = 0;
  mbb_		= -9999.;
  mcmln_	= -9999.;
  mcmtln_	= -9999.;
  pflepmindrj_   = 9999.;
  pftaudmindrj_  = 9999.;
  lep1pfjetdr_   = 9999.;
  lep2pfjetdr_   = 9999.;

  pfcand5_        = 0;
  pfcand10_       = 0;
  pfcandiso5_     = 9999.;     
  pfcandiso10_    = 9999.;     
  pfcandpt5_      = 9999.;
  pfcandpt10_     = 9999.;
  pfcandmindrj5_  = 9999.;
  pfcandmindrj10_ = 9999.;

  trkpt10pt0p1_	    = 9999.;
  trkreliso10pt0p1_ = 9999.;
  trkpt10pt0p2_	    = 9999.;
  trkreliso10pt0p2_ = 9999.;
  trkpt10pt0p3_	    = 9999.;
  trkreliso10pt0p3_ = 9999.;
  trkpt10pt0p4_	    = 9999.;
  trkreliso10pt0p4_ = 9999.;
  trkpt10pt0p5_	    = 9999.;
  trkreliso10pt0p5_ = 9999.;
  trkpt10pt0p6_	    = 9999.;
  trkreliso10pt0p6_ = 9999.;
  trkpt10pt0p7_	    = 9999.;
  trkreliso10pt0p7_ = 9999.;
  trkpt10pt0p8_	    = 9999.;
  trkreliso10pt0p8_ = 9999.;
  trkpt10pt0p9_	    = 9999.;
  trkreliso10pt0p9_ = 9999.;
  trkpt10pt1p0_	    = 9999.;
  trkreliso10pt1p0_ = 9999.;

  pfcandpt10pt0p1_  = 9999.;
  pfcandiso10pt0p1_ = 9999.;
  pfcandpt10pt0p2_  = 9999.;
  pfcandiso10pt0p2_ = 9999.;
  pfcandpt10pt0p3_  = 9999.;
  pfcandiso10pt0p3_ = 9999.;
  pfcandpt10pt0p4_  = 9999.;
  pfcandiso10pt0p4_ = 9999.;
  pfcandpt10pt0p5_  = 9999.;
  pfcandiso10pt0p5_ = 9999.;
  pfcandpt10pt0p6_  = 9999.;
  pfcandiso10pt0p6_ = 9999.;
  pfcandpt10pt0p7_  = 9999.;
  pfcandiso10pt0p7_ = 9999.;
  pfcandpt10pt0p8_  = 9999.;
  pfcandiso10pt0p8_ = 9999.;
  pfcandpt10pt0p9_  = 9999.;
  pfcandiso10pt0p9_ = 9999.;
  pfcandpt10pt1p0_  = 9999.;
  pfcandiso10pt1p0_ = 9999.;

}

//--------------------------------------------------------------------

int getProcessType(char *prefix)
{
  int proc = -1;

  if(strcmp(prefix,"data")   == 0) proc = 0;
  if(strcmp(prefix,"Zjets")  == 0) proc = 1;
  if(strcmp(prefix,"ttdil")  == 0) proc = 2;
  if(strcmp(prefix,"ttotr")  == 0) proc = 3;
  if(strcmp(prefix,"ww")     == 0) proc = 4;
  if(strcmp(prefix,"wz")     == 0) proc = 5;
  if(strcmp(prefix,"zz")     == 0) proc = 6;
  if(strcmp(prefix,"wjets")  == 0) proc = 7;
  if(strcmp(prefix,"tW")     == 0) proc = 8;
  if(strcmp(prefix,"LM0")    == 0) proc = 10;
  if(strcmp(prefix,"LM1")    == 0) proc = 11;
  if(strcmp(prefix,"LM2")    == 0) proc = 12;
  if(strcmp(prefix,"LM3")    == 0) proc = 13;
  if(strcmp(prefix,"LM4")    == 0) proc = 14;
  if(strcmp(prefix,"LM5")    == 0) proc = 15;
  if(strcmp(prefix,"LM6")    == 0) proc = 16;
  if(strcmp(prefix,"LM7")    == 0) proc = 17;
  if(strcmp(prefix,"LM8")    == 0) proc = 18;
  if(strcmp(prefix,"LM9")    == 0) proc = 19;
  if(strcmp(prefix,"LM10")   == 0) proc = 20;
  if(strcmp(prefix,"LM11")   == 0) proc = 21;
  if(strcmp(prefix,"LM12")   == 0) proc = 22;

  return proc;
}

//--------------------------------------------------------------------

pair<float, float> ScaleMET( pair<float, float> p_met, LorentzVector p4_dilep, double rescale = 1.0 ){
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

void singleLeptonLooper::closeTree()
{
  outFile->cd();
  outTree->Write();
  outFile->Close();
  delete outFile;
}

//--------------------------------------------------------------------

float singleLeptonLooper::stopPairCrossSection( float stopmass ){

  int   bin  = stop_xsec_hist->FindBin(stopmass);
  float xsec = stop_xsec_hist->GetBinContent(bin);
  return xsec;

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

//--------------------------------------------------------------------

float getMT( float leppt , float lepphi , float met , float metphi ) {
  float dphi = fabs( lepphi - metphi );
      if( dphi > TMath::Pi() ) dphi = TMath::TwoPi() - dphi;
      return sqrt( 2 * ( leppt * met * (1 - cos( dphi ) ) ) );
}

//--------------------------------------------------------------------

bool passSingleMuTrigger2011_pt30( bool isData , int lepType ) {
  
  //----------------------------
  // single muon triggers
  //----------------------------

  // no triggers required for MC
  //if( !isData ) return true;

  // false for electron channel
  if( lepType == 0 ){
    return false;
  }

  // muon channel
  else if( lepType == 1 ){    
    if( passUnprescaledHLTTriggerPattern("HLT_IsoMu30_v") )          return true; //  < 173212
    if( passUnprescaledHLTTriggerPattern("HLT_IsoMu30_eta2p1_v") )   return true; // >= 173212
  }

  else{
    cout << __FILE__ << " " << __LINE__ << " ERROR unrecognized lepType " << lepType << ", quitting" << endl;
    exit(0);
  }

  return false;
}

//--------------------------------------------------------------------

float getMuTriggerWeight( float pt, float eta ) {
  //Trigger efficiency for single muon triggers averaged over full 2011 dataset
  //From AN2011-456 Table 28
  float trigweights[2][3] = {{0.9002, 0.8352, 0.8266},
			     {0.9440, 0.8821, 0.8611}};
  int i_pt = -1;
  if ( pt > 30 && pt < 40 ) i_pt = 0;
  else if ( pt > 40 ) i_pt = 1;
  if ( i_pt < 0 ) return 1.;

  int i_eta = -1;
  if ( fabs(eta) < 0.8 ) i_eta = 0;
  else if ( fabs(eta) >= 0.8 && abs(eta) < 1.5 ) i_eta = 1;
  else if ( fabs(eta) >= 1.5 && abs(eta) < 2.1 ) i_eta = 2;
  if ( i_eta < 0 ) return 1.;

  return trigweights[i_pt][i_eta];

}


//--------------------------------------------------------------------

float getMuTriggerWeightNew( float pt, float eta ) {

  float aeta = fabs(eta);

  if( aeta < 0.8 ){
    if( pt >  30.0 && pt <   40.0 ) return 0.89;
    if( pt >  40.0 && pt <   60.0 ) return 0.89;
    if( pt >  60.0 && pt <  100.0 ) return 0.88;
    if( pt > 100.0                ) return 0.87;
  }
  
  else if( aeta > 0.8 && aeta < 1.3 ){
    if( pt >  30.0 && pt <   40.0 ) return 0.81;
    if( pt >  40.0 && pt <   60.0 ) return 0.82;
    if( pt >  60.0 && pt <  100.0 ) return 0.81;
    if( pt > 100.0                ) return 0.80;
  }

  else if( aeta > 1.3 && aeta < 1.8 ){
    if( pt >  30.0 && pt <   40.0 ) return 0.82;
    if( pt >  40.0 && pt <   60.0 ) return 0.84;
    if( pt >  60.0 && pt <  100.0 ) return 0.82;
    if( pt > 100.0                ) return 0.80;
  }

  else if( aeta > 1.8 && aeta < 2.0 ){
    if( pt >  30.0 && pt <   40.0 ) return 0.79;
    if( pt >  40.0 && pt <   60.0 ) return 0.81;
    if( pt >  60.0 && pt <  100.0 ) return 0.80;
    if( pt > 100.0                ) return 0.79;
  }

  else if( aeta > 2.0 && aeta < 2.1 ){
    if( pt >  30.0 && pt <   40.0 ) return 0.69;
    if( pt >  40.0 && pt <   60.0 ) return 0.71;
    if( pt >  60.0 && pt <  100.0 ) return 0.70;
    if( pt > 100.0                ) return 0.70;
  }


  return 1;
}

//--------------------------------------------------------------------

float getminjdr( VofP4 jets, LorentzVector *particle ) {
  float mindr = 9999.;
  if (jets.size()==0 || particle==0) return mindr;
  for ( unsigned int ijet = 0; ijet<jets.size(); ++ijet ) {
    float partjdr = dRbetweenVectors(jets.at(ijet),*particle);
    if ( partjdr<mindr ) mindr = partjdr;
  }
  return mindr;
}

//--------------------------------------------------------------------

int singleLeptonLooper::ScanChain(TChain* chain, char *prefix, float kFactor, int prescale, float lumi,
				  FREnum frmode, bool doFakeApp)

{

  bool isLM = TString(prefix).Contains("LM");
  bool isData = false;
  if( TString(prefix).Contains("data") || TString(prefix).Contains("2012") 
      || TString(prefix).Contains("dimu") || TString(prefix).Contains("diel")
      || TString(prefix).Contains("mueg") ){
    cout << "DATA!!!" << endl;
    isData       = true;
    doTenPercent = false;
  }
  if( doTenPercent ) cout << "Processing 10% of MC" << endl;

  //------------------------------------------------------------------------------------------------------
  // set json, vertex reweighting function and msugra cross section files
  //------------------------------------------------------------------------------------------------------
  
  if( !initialized ){

    //set json
    cout << "setting json " << g_json << endl;
    set_goodrun_file( g_json );

    if( TString(prefix).Contains("ttall_massivebin") ) 
      set_vtxreweight_rootfile("vtxreweight/vtxreweight_Summer12MC_PUS6.root",true);
    else
      set_vtxreweight_rootfile("vtxreweight/vtxreweight_Summer12MC_PUS7_5p1fb_Zselection.root",true);

    /*
    //set vtx reweighting hist - depends on pileup scenario
    if( TString(prefix).Contains("ttall") || TString(prefix).Contains("mcatnlo") 
     || TString(prefix).Contains("scaleup") || TString(prefix).Contains("scaledw") 
     || TString(prefix).Contains("matchup") )
      set_vtxreweight_rootfile("vtxreweight/vtxreweight_Fall11MC_PUS6_4p7fb_Zselection.root",true);
    else if( TString(prefix).Contains("pythia") )
      set_vtxreweight_rootfile("vtxreweight/vtxreweight_Summer11-PU_S3_START42_V11_4p7fb_ttbar.root",true);
    else
      set_vtxreweight_rootfile("vtxreweight/vtxreweight_Summer11MC_PUS4_4p7fb_Zselection.root",true);
    */

    //    weight3D_init( "vtxreweight/Weight3D.root" );

    //set msugra cross section file
    set_msugra_file("goodModelNames_tanbeta10.txt");


    initialized = true;
  }

  //------------------------------------------------------------------------------------------------------
  // latest-and-greatest JEC
  //------------------------------------------------------------------------------------------------------

  std::vector<std::string> jetcorr_filenames_pfL1FastJetL2L3;
  FactorizedJetCorrector *jet_corrector_pfL1FastJetL2L3;

  jetcorr_filenames_pfL1FastJetL2L3.clear();
  
  string pfUncertaintyFile;
  //string caloUncertaintyFile;

  if ( isData ) {
    jetcorr_filenames_pfL1FastJetL2L3.push_back  ("jetCorrections/GR_R_52_V9_L1FastJet_AK5PF.txt");
    jetcorr_filenames_pfL1FastJetL2L3.push_back  ("jetCorrections/GR_R_52_V9_L2Relative_AK5PF.txt");
    jetcorr_filenames_pfL1FastJetL2L3.push_back  ("jetCorrections/GR_R_52_V9_L3Absolute_AK5PF.txt");
    jetcorr_filenames_pfL1FastJetL2L3.push_back  ("jetCorrections/GR_R_52_V9_L2L3Residual_AK5PF.txt");

    pfUncertaintyFile = "jetCorrections/GR_R_52_V9_Uncertainty_AK5PF.txt";
  } 
  else {
    jetcorr_filenames_pfL1FastJetL2L3.push_back  ("jetCorrections/START52_V9B_L1FastJet_AK5PF.txt");
    jetcorr_filenames_pfL1FastJetL2L3.push_back  ("jetCorrections/START52_V9B_L2Relative_AK5PF.txt");
    jetcorr_filenames_pfL1FastJetL2L3.push_back  ("jetCorrections/START52_V9B_L3Absolute_AK5PF.txt");
    
    pfUncertaintyFile = "jetCorrections/START52_V9B_Uncertainty_AK5PF.txt";
  }

  jet_corrector_pfL1FastJetL2L3  = makeJetCorrector(jetcorr_filenames_pfL1FastJetL2L3);

  JetCorrectionUncertainty *pfUncertainty   = new JetCorrectionUncertainty( pfUncertaintyFile   );

  //------------------------------------------------
  // set stop cross section file
  //------------------------------------------------

  stop_xsec_file = TFile::Open("stop_xsec.root");
  
  if( !stop_xsec_file->IsOpen() ){
    cout << "Error, could not open stop cross section TFile, quitting" << endl;
    exit(0);
  }
  
  stop_xsec_hist        = (TH1D*) stop_xsec_file->Get("h_stop_xsec");
  
  if( stop_xsec_hist == 0 ){
    cout << "Error, could not retrieve stop cross section hist, quitting" << endl;
    exit(0);
  }

  // instanciate topmass solver REPLACETOPMASS
  //ttdilepsolve * d_llsol = new ttdilepsolve;

  //instantiate SimpleFakeRate class for electrons and muons
  //this is the default, can change it below if needed
  SimpleFakeRate* mufr = 0;
  SimpleFakeRate* elfr = 0;

  if(doFakeApp) {

    cout << "NOT CURRENTLY SET UP FOR FAKES!!!!! QUITTING!" << endl;
    exit(0);

    std::cout<<"**************************"<<std::endl;
    std::cout<<"Running FR application job"<<std::endl;
    std::cout<<"**************************"<<std::endl;

    if(isData) {
      std::cout<<"Using data derived FR files"<<std::endl;
      mufr = new SimpleFakeRate("fr_os7June2011.root", "fr_mu_OSGV3" );
      elfr = new SimpleFakeRate("fr_os7June2011.root", "fr_el_OSGV3" );
    }
    else {
      std::cout<<"Using data derived FR files"<<std::endl;
      std::cout<<"CURRENTLY USING DATA FR FOR MC FIXME!!!!!" <<std::endl;
      mufr = new SimpleFakeRate("fr_os7June2011.root", "fr_mu_OSGV3" );
      elfr = new SimpleFakeRate("fr_os7June2011.root", "fr_el_OSGV3" );
    }
  }

  TDirectory *rootdir = gDirectory->GetDirectory("Rint:");
  rootdir->cd();

  BookHistos(prefix);
  
  unsigned int nEventsChain = chain->GetEntries();
  unsigned int nEventsTotal = 0;
  // map isn't needed for this purpose, vector is sufficient
  // better would be to use a struct with run, lb, event
  map<int,int> m_events;

  // loop over files
  TObjArray *listOfFiles = chain->GetListOfFiles();
  TIter fileIter(listOfFiles);
  TChainElement* currentFile = 0;

  int nSkip_els_conv_dist = 0;

  float netot  = 0.;
  float nmtot  = 0.;
  float nepass = 0.;
  float nmpass = 0.;

  if(g_createTree) makeTree(prefix, doFakeApp, frmode);

  while((currentFile = (TChainElement*)fileIter.Next())) {
    TFile* f = new TFile(currentFile->GetTitle());

    cout << currentFile->GetTitle() << endl;

    if( !f || f->IsZombie() ) {
      cout << "Skipping bad input file: " << currentFile->GetTitle() << endl;
      continue; //exit(1);                                                                                             
    }

    TTree *tree = (TTree*)f->Get("Events");

    //Matevz
    TTreeCache::SetLearnEntries(100);
    tree->SetCacheSize(128*1024*1024);

    cms2.Init(tree);
      
    unsigned int nEntries = tree->GetEntries();

    for(unsigned int z = 0; z < nEntries; ++z) {
      ++nEventsTotal;

      if( doTenPercent ){
	if( !(nEventsTotal%10==0) ) continue;
      }

      // progress feedback to user
      if (nEventsTotal % 1000 == 0){
        
        // xterm magic from L. Vacavant and A. Cerri
        if (isatty(1)){
                
          printf("\015\033[32m ---> \033[1m\033[31m%4.1f%%"
                 "\033[0m\033[32m <---\033[0m\015", (float)nEventsTotal/(nEventsChain*0.01));
          fflush(stdout);
        }
      }

      //Matevz
      tree->LoadTree(z);

      cms2.GetEntry(z);

      // PrintTriggers();
      // exit(0);

      if( cms2.evt_ww_rho_vor() != cms2.evt_ww_rho_vor() ){
	cout << "Skipping event with rho = nan!!!" << endl;
	continue;
      }

      InitBaby();

      isdata_ = isData ? 1 : 0;

      if( verbose ){
	cout << "-------------------------------------------------------"   << endl;
	cout << "Event " << z                                               << endl;
	cout << "File  " << currentFile->GetTitle()                         << endl;
	cout << evt_dataset().at(0) << " " << evt_run() << " " << evt_lumiBlock() << " " << evt_event() << endl;
	cout << "-------------------------------------------------------"   << endl;
      }

      TString datasetname(evt_dataset().at(0));
      bool isperiodA = datasetname.Contains("2011A") ? true : false;
      //      cout<<"dataset: "<<datasetname.Data()<<" isperiodA: "<<isperiodA<<endl;

      // skip stop-pair events with m(stop) > 850 GeV
      // if( TString(prefix).Contains("T2") ){
      // 	if( sparm_mG() > 600.0 ) continue;
      // }

      //---------------------------------------------
      // event cleaning and good run list
      //---------------------------------------------

      if( !cleaning_goodVertexApril2011() )                          continue;
      if( isData && !goodrun(cms2.evt_run(), cms2.evt_lumiBlock()) ) continue;

      //---------------------
      // skip duplicates
      //---------------------

      if( isData ) {
        DorkyEventIdentifier id = { evt_run(),evt_event(), evt_lumiBlock() };
        if (is_duplicate(id) ){
          continue;
        }
      }

      //-------------------------------------
      // skip events with bad els_conv_dist
      //-------------------------------------

      bool skipEvent = false;
      for( unsigned int iEl = 0 ; iEl < els_conv_dist().size() ; ++iEl ){
        if( els_conv_dist().at(iEl) != els_conv_dist().at(iEl) ){
          skipEvent = true;
        }
        if( els_sigmaIEtaIEta().at(iEl) != els_sigmaIEtaIEta().at(iEl) ){
          skipEvent = true;
        }
        if( els_sigmaIEtaIEtaSC().at(iEl) != els_sigmaIEtaIEtaSC().at(iEl) ){
          skipEvent = true;
        }
      }
             
      if( skipEvent ){
        nSkip_els_conv_dist++;
        continue;
      }
   
      //---------------------------------------------
      // find leptons passing analysis selection
      //---------------------------------------------

      VofP4 goodLeptons;
      vector<float> lepchi2ndf;
      vector<float> lepdpt;
      vector<int> lepId;
      vector<int> lepIndex;

      ngoodlep_ = 0;
      ngoodel_  = 0;
      ngoodmu_  = 0;
            
      for( unsigned int iel = 0 ; iel < els_p4().size(); ++iel ){
	if( els_p4().at(iel).pt() < 10 )                                                                continue;
        if( !passElectronSelection_Stop2012_v2( iel , vetoTransition,vetoTransition,useOldIsolation) )  continue;

	goodLeptons.push_back( els_p4().at(iel) );
	lepchi2ndf.push_back( -9999. );
	lepdpt.push_back( -9999. );
	lepId.push_back( els_charge().at(iel) * 11 );
	lepIndex.push_back(iel);
	ngoodel_++;
	ngoodlep_++;

	//cout << "Found electron " << ngoodlep_ << " pt " << els_p4().at(iel).pt() << endl;
      }
          
      for( unsigned int imu = 0 ; imu < mus_p4().size(); ++imu ){
	if( mus_p4().at(imu).pt() < 10 )           continue;
	if( !muonId( imu , ZMet2012_v1 ))         continue;

	goodLeptons.push_back( mus_p4().at(imu) );
	//in original OSGeneric_v3 version, cut on chi2ndf is at 10
	lepchi2ndf.push_back( cms2.mus_gfit_chi2().at(imu)/cms2.mus_gfit_ndof().at(imu) );
	//in original OSGeneric_v3, cut on dpt/pt is at 0.1
	lepdpt.push_back( mus_ptErr().at(imu) / mus_p4().at(imu).pt() );
	lepId.push_back( mus_charge().at(imu) * 13 );
	lepIndex.push_back(imu);
	ngoodmu_++;
	ngoodlep_++;

	//cout << "Found muon " << ngoodlep_ << " pt " << mus_p4().at(imu).pt() << endl;
      }  

      // REQUIRE AT LEAST 1 GOOD LEPTON!!!
      if( goodLeptons.size() < 1 ) continue;

      //---------------------------------------------
      // find leading lepton
      //---------------------------------------------
      float maxpt   = -1;
      int   imaxpt  = -1;
      int   imaxpt2 = -1;

      for( unsigned int ilep = 0 ; ilep < goodLeptons.size() ; ilep++ ){
	if( goodLeptons.at(ilep).pt() > maxpt ){
	  maxpt  = goodLeptons.at(ilep).pt();
	  imaxpt = ilep;
	}
      }

      if( imaxpt < 0 ){
	cout << "ERROR! QUITTING imaxpt " << imaxpt << endl;
	exit(2);
      }

      // REQUIRE LEADING LEPTON PT > 20 GEV
      if( maxpt < 20 ) continue;

      id1_         = lepId.at(imaxpt);
      lep1_        = &goodLeptons.at(imaxpt);
      lep1chi2ndf_ = lepchi2ndf.at(imaxpt);
      lep1dpt_     = lepdpt.at(imaxpt);

      int index1 = lepIndex.at(imaxpt);

      // Matching lepton with pfjet
      LorentzVector minvjet1;
      for (unsigned int ijet = 0 ; ijet < pfjets_p4().size() ; ijet++) {
	
	LorentzVector vjet = pfjets_p4().at(ijet);

	if( vjet.pt()  < 10  )             continue;
	if( fabs(vjet.eta()) > 3.0 )       continue;
        
	float dr = dRbetweenVectors( vjet, *lep1_ );

	if( dr < lep1pfjetdr_ ){
	  lep1pfjetdr_    = dr;
	  minvjet1         = vjet;
	}
      }
      if (lep1pfjetdr_<9998.) 
	leppfjet1_ = &minvjet1;

      // Matching lepton with pflepton
      // and with trk and gfit where available
      if( abs(id1_) == 13 ) {

        int ipf1 = mus_pfmusidx().at(index1);
        if( ipf1 >= 0 ) pflep1_ = &(pfmus_p4().at(ipf1));

	trklep1_ = &(mus_trk_p4().at(index1));
	gfitlep1_ = &(mus_gfit_p4().at(index1));

	float dtrkpt1  = fabs( trklep1_->Pt()  - lep1_->Pt() );
	float dgfitpt1 = fabs( gfitlep1_->Pt() - lep1_->Pt() );

	if ( dtrkpt1<=0.1 && dgfitpt1<=0.1 ) leptype1_ =  2;
	else if ( dtrkpt1 <=0.1 )            leptype1_ =  0;
	else if ( dgfitpt1<=0.1 )            leptype1_ =  1;
	else                                 leptype1_ = -1;

      } else if( abs(id1_) == 11 ) {
        int ipf1 = els_pfelsidx().at(index1);
        if( ipf1 >= 0 ) pflep1_ = &(pfels_p4().at(ipf1));
	trklep1_ = &(els_trk_p4().at(index1));
      }

      //cout << "Leading lepton: pt " << lep1_->pt() << " id " << id1_ << endl;

      //---------------------------------------------
      // find 2nd leading lepton (if >=2 leptons)
      //---------------------------------------------

      id2_       = -999;
      lep2_      = 0;
      int index2 = -1;
      dilmass_       = -999;
      dilpt_         = -999.;
      dilrecoil_     = -999.;
      dilrecoilparl_ = -999.;
      dilrecoilperp_ = -999.;

      hyptype_ = -1;

      if( ngoodlep_ > 1 ){

	maxpt = -1;

	for( unsigned int ilep = 0 ; ilep < goodLeptons.size() ; ilep++ ){

	  if( (int)ilep == imaxpt ) continue;

	  if( goodLeptons.at(ilep).pt() > maxpt ){
	    maxpt   = goodLeptons.at(ilep).pt();
	    imaxpt2 = ilep;
	  }
	}

	if( imaxpt2 < 0 ){
	  cout << "ERROR! QUITTING imaxpt2 " << imaxpt2 << endl;
	  exit(3);
	}

	id2_         = lepId.at(imaxpt2);
	lep2_        = &goodLeptons.at(imaxpt2);
	lep2chi2ndf_ = lepchi2ndf.at(imaxpt2);
	lep2dpt_     = lepdpt.at(imaxpt2);
	index2       = lepIndex.at(imaxpt2);
	dilmass_     = sqrt((goodLeptons.at(imaxpt) + goodLeptons.at(imaxpt2)).mass2());

	if     ( abs(id1_) == 11 && abs(id2_) == 11 ) hyptype_ = 3;
	else if( abs(id1_) == 13 && abs(id2_) == 13 ) hyptype_ = 0;
	else if( abs(id1_) == 11 && abs(id2_) == 13 ) hyptype_ = 2;
	else if( abs(id1_) == 13 && abs(id2_) == 11 ) hyptype_ = 1;
	else{
	  cout << "ERROR! " << id1_ << " " << id2_ << endl;
	}

	if( id1_ > 0 && id2_ < 0 ){
	  lepp_ = &goodLeptons.at(imaxpt);
	  lepm_ = &goodLeptons.at(imaxpt2);
	}
	else{
	  lepp_ = &goodLeptons.at(imaxpt2);
	  lepm_ = &goodLeptons.at(imaxpt);
	}

	// Matching lepton with pfjet
	LorentzVector minvjet2;
	for (unsigned int ijet = 0 ; ijet < pfjets_p4().size() ; ijet++) {
	  
	  LorentzVector vjet = pfjets_p4().at(ijet);
	  
	  if( vjet.pt()  < 10  )             continue;
	  if( fabs(vjet.eta()) > 3.0 )       continue;
	  
	  float dr = dRbetweenVectors( vjet, *lep2_ );
	  
	  if( dr < lep2pfjetdr_ ){
	    lep2pfjetdr_    = dr;
	    minvjet2        = vjet;
	  }
	}
	if (lep2pfjetdr_<9998.)
	  leppfjet2_      = &minvjet2;      

	//---------------------------------------------
	// Matching sub-leading lepton with pflepton
	// and to trk and gfit where available
	//---------------------------------------------
	if( abs(id2_) == 13 ) {

	  int ipf2 = mus_pfmusidx().at(index2);
	  if( ipf2 >= 0 ) pflep2_ = &(pfmus_p4().at(ipf2));
 	  
	  trklep2_ = &(mus_trk_p4().at(index2));
	  gfitlep2_ = &(mus_gfit_p4().at(index2));

	  float dtrkpt2  = fabs( trklep2_->Pt()  - lep2_->Pt() );
	  float dgfitpt2 = fabs( gfitlep2_->Pt() - lep2_->Pt() );

	  if ( dtrkpt2<=0.1 && dgfitpt2<=0.1 ) leptype2_ =  2;
	  else if ( dtrkpt2 <=0.1 )            leptype2_ =  0;
	  else if ( dgfitpt2<=0.1 )            leptype2_ =  1;
	  else                                 leptype2_ = -1;

	} else if( abs(id2_) == 11 ) {
	  int ipf2 = els_pfelsidx().at(index2);
	  if( ipf2 >= 0 ) pflep2_ = &(pfels_p4().at(ipf2));
	  trklep2_ = &(els_trk_p4().at(index2));
	}
	
        float metx = evt_pfmet() * cos( evt_pfmetPhi() );
        float mety = evt_pfmet() * sin( evt_pfmetPhi() );
        TVector2 met( metx, mety );

        LorentzVector dilep = goodLeptons.at(imaxpt) + goodLeptons.at(imaxpt2);
        TVector2 dilepv( dilep.X(), dilep.Y() );
        dilpt_ = dilepv.Mod();

        //axes
        TVector2 ptaxis = dilepv/dilpt_;
        TVector2 perpaxis = ptaxis.Rotate(TMath::Pi()/2.);
        //hadronic recoil information
        TVector2 hadrecoil = -1.*(met+dilepv);
        dilrecoil_ = hadrecoil.Mod();
        dilrecoilparl_ = hadrecoil*ptaxis;
        dilrecoilperp_ = hadrecoil*perpaxis;

	//cout << "2nd deading lepton: pt " << lep2_->pt() << " id " << id2_ << " mass " << dilmass_ << endl;
      }

      //--------------------------------
      // store dilepton type in myType
      //--------------------------------

      leptype_ = 99;
      if      ( abs(id1_) == 11 )                   leptype_ = 0; // e
      else if ( abs(id1_) == 13 )                   leptype_ = 1; // m
      else{
	cout << "Skipping unknown lepton type = " << id1_ << endl;
	continue;
      }

      //--------------------------------
      // require trigger
      //--------------------------------

      //int hypType = (leptype_==1) ? 0 : 3;//lepton type for dilepton triggers
      //if( !passSingleLepSUSYTrigger2011_v1( isData , leptype_ ) && !passSUSYTrigger2011_v1( isData , hypType , true ) ) continue;

      //-----------------------------------------------------------------
      // number of OS generic electrons *in addition to* primary lepton
      //-----------------------------------------------------------------

      nosel_ = 0;

      for( unsigned int iel = 0 ; iel < els_p4().size(); ++iel ){
	if( els_p4().at(iel).pt() < 10 )                                                 continue;
	if( !pass_electronSelection( iel , electronSelection_el_OSV3 , false , false ) ) continue;
	if( dRbetweenVectors( *lep1_ , els_p4().at(iel) ) < 0.1 )                        continue;
	nosel_++;
      }
      
      //--------------------------------
      // get MC quantities
      //--------------------------------
      
      int nels       =  0;
      int nmus       =  0;
      int ntaus      =  0;
      int nleps      =  0;
      float dilptgen = -1;

      ptwgen_   = -1;
      ptzgen_   = -1;
      ptttbar_  = -1;
      ptt_      = -1;
      pttbar_   = -1;
      mttbar_   = -1;
      etattbar_ = -999;
      t_        = 0;
      tbar_     = 0;
      ttbar_    = 0;

      npartons_    =  0;
      maxpartonpt_ = -1;

      mgcor_ = 1.0;
      wflav_ = -1;

      if( !isData ){

	w1_     = leptonOrTauIsFromW( index1 , id1_ , isLM );
	pthat_  = genps_pthat();
	qscale_ = genps_qScale();
	
	//store W flavor history
	if (!TString(prefix).Contains("mcatnlo"))
	  wflav_ = (int)genps_flavorHistoryFilterResult();

	//splitting ttbar into ttdil/ttotr
	//nleps = leptonGenpCount_lepTauDecays(nels, nmus, ntaus);
	nleps = leptonGenpCount(nels, nmus, ntaus);
	
	nels_  = nels;
	nmus_  = nmus;
	ntaus_ = ntaus;
	nleps_ = nleps;

	// this is a weight which corrects for the wrong MG W->lnu BF
	if( TString(prefix).Contains("ttall") ){
	  if( nleps == 0 ) mgcor_ = 1.028;
	  if( nleps == 1 ) mgcor_ = 0.986;
	  if( nleps == 2 ) mgcor_ = 0.945;
	}

	if( strcmp(prefix,"ttem")  == 0 && ( nels + nmus ) != 2 ) continue;
	if( strcmp(prefix,"ttdil") == 0 && nleps != 2           ) continue;
	if( strcmp(prefix,"ttotr") == 0 && nleps == 2           ) continue;
	
	LorentzVector vdilepton(0,0,0,0);
	LorentzVector vttbar(0,0,0,0);
	int ntops = 0;

	for ( int igen = 0 ; igen < (int)genps_id().size() ; igen++ ) { 
	  if ( abs( cms2.genps_id().at(igen) ) == 11) vdilepton += genps_p4().at(igen); 
	  if ( abs( cms2.genps_id().at(igen) ) == 13) vdilepton += genps_p4().at(igen); 

	  int id = cms2.genps_id().at(igen);

	  if( id == 6 ){
	    t_         = &(genps_p4().at(igen));
	    ptt_       = genps_p4().at(igen).pt();
	    vttbar    += genps_p4().at(igen);
	    ntops++;
	  }
	  if( id == -6 ){
	    tbar_      = &(genps_p4().at(igen));
	    pttbar_    = genps_p4().at(igen).pt();
	    vttbar    += genps_p4().at(igen);
	    ntops++;
	  }

	  // store W or Z pT 
	  // ignoring cases where have more than 1 boson for now
	  if ( abs(id) == 24 )
	    ptwgen_ = genps_p4().at(igen).pt();
	  if ( abs(id) == 23 )
	    ptzgen_ = genps_p4().at(igen).pt();

	  // skip lines up to t and tbar
	  if( igen < 8 ) continue;

	  // require particle is a quark or a gluon
	  int pid = abs( genps_id().at(igen) );
	  if( !( pid==1 || pid==2 || pid==3 || pid==4 || pid==5 || pid==6 || pid == 21 ) ) continue;

	  // require mother is not a top or W
	  int mothid = abs(genps_id_mother().at(igen));
	  if( mothid == 6 || mothid == 24) continue;

	  // found additional parton
	  npartons_ ++;
	  if( genps_p4().at(igen).pt() > maxpartonpt_ ) maxpartonpt_ = genps_p4().at(igen).pt();
	  //	  cout << "found parton, igen " << igen << " id " << pid << " motherid " << mothid << " pt " << genps_p4().at(igen).pt() << endl;

	}

	// if( npartons_ > 0 ){
	//   cout << endl << endl;
	//   dumpDocLines();
	//   cout << endl << endl;
	//   cout << "number of partons " << npartons_    << endl;
	//   cout << "max parton pt     " << maxpartonpt_ << endl;
	// }

	//count tops and only get two
	//ttbar_    = &(vttbar);
	if (ntops==2) {
	  LorentzVector ttpair = *t_ + *tbar_;
	  ttbar_    = &ttpair;
	  ptttbar_  = ttbar_->pt();
	  mttbar_   = ttbar_->mass();
	  etattbar_ = ttbar_->eta();
	}
	
	if( nels + nmus == 2) dilptgen = vdilepton.pt();
        
	if ( strcmp(prefix , "DYee"     ) == 0 &&  nels  != 2  ) continue;
	if ( strcmp(prefix , "DYmm"     ) == 0 &&  nmus  != 2  ) continue;
	if ( strcmp(prefix , "DYtautau" ) == 0 &&  ntaus != 2  ) continue;
	
	//splice together the DY samples - if its madgraph, then we do nothing
	if(TString(prefix).Contains("DY") && TString(evt_dataset().at(0)).Contains("madgraph") == false) {	
	  bool doNotContinue = false;
	  for(unsigned int i = 0; i < genps_p4().size(); i++){
	    if(abs(genps_id()[i]) == 23 && genps_p4()[i].M() > 50.)
	      doNotContinue = true;
	  }
	  if(doNotContinue)
	    continue;	
	}
	
	//extract pthat
	if(TString(prefix).Contains("DY")){
	  int nz = 0;
	  for(unsigned int i = 0; i < genps_p4().size(); i++){
	    if(abs(genps_id()[i]) == 23){
	      mllgen_ = genps_p4()[i].M();
	      nz++;
	    }
	  }
	  if(nz != 1 ) cout << "ERROR NZ " << nz << endl;
	}
            
	//-----------------------------------------------------
	// track gen lepton information here
	//-----------------------------------------------------
	
	mcid1_      = -1;
	mcid2_      = -1;
	mclep1_     =  0;
	mclep2_     =  0;
	mcdecay1_   = -1;
	mcdecay2_   = -1;
	mcdr1_      = -1;
	mcdr2_      = -1;
	mcndec1_    =  0;
	mcndec2_    =  0;
	mcndeckls1_ =  0;
	mcndeckls2_ =  0;
	mcndecem1_  =  0;
	mcndecem2_  =  0;
	mctaudpt1_  = -1;
	mctaudpt2_  = -1;
	mctaud1_    =  0;
	mctaud2_    =  0;
	mctaudvis1_   =  0;
	mctaudvis2_   =  0;
	mctaudid1_  = -1;
	mctaudid2_  = -1;
	//variables for single lepton studies
	mcnu_       =  0;
	mclep_      =  0;

	//-----------------------------------------------------
	// store single gen lepton info
	//-----------------------------------------------------

	if( nleps_ == 1 ){


	  int nfoundleps = 0;

	  for ( int igen = 0 ; igen < (int)cms2.genps_id().size() ; igen++ ) { 

	    int id = genps_id().at(igen);

	    if( !( abs(id)==11 || abs(id)==13 || abs(id)==15 ) ) continue;

	    nfoundleps++;
	    mcid1_   = id;
	    mclep1_  = &genps_p4().at(igen);
	    mcdr1_   = dRbetweenVectors( *lep1_ , *mclep1_ );

	    if( abs(id)==15 ){
	      mcdecay1_ = 1;
	      //variables to calculate the visible tau energy
	      LorentzVector taudvis; 
	      for(unsigned int kk = 0; kk < cms2.genps_lepdaughter_id().at(igen).size(); kk++) {
		int daughter = abs(cms2.genps_lepdaughter_id()[igen][kk]);
		//4-vector of visible part of tau
		if ( daughter != 12 && daughter != 14 && daughter != 16 ) {
		  if (taudvis.pt()==0)
		    taudvis = genps_lepdaughter_p4()[igen][kk];
		  else 
		    taudvis = genps_lepdaughter_p4()[igen][kk]+taudvis;
		}
		// get tau daughter pt
		if( daughter == 11 || daughter == 13 || daughter == 211 || daughter == 321 ){
		  if( genps_lepdaughter_p4()[igen][kk].pt() > mctaudpt1_ ) {
		    mctaudpt1_ = genps_lepdaughter_p4()[igen][kk].pt();
		    mctaudid1_ = cms2.genps_lepdaughter_id()[igen][kk];
		    mctaud1_   = &genps_lepdaughter_p4()[igen][kk];
		  }
		}

		if( daughter == 211 || daughter == 321 ) mcndec1_  ++;  // count charged hadrons
		if( daughter == 130 || daughter == 310 ) mcndeckls1_  ++;  // count K_L and K_S
		if( daughter ==  22 || daughter == 111 ) mcndecem1_  ++;  // count photons and pi0
		if( daughter ==  12 || daughter ==  14 ) mcdecay1_ = 2; // check for nu_e or nu_mu 
	      }// end loop over daughters
	      mctaudvis1_ = &taudvis;
	      
	    }
	  } 


	  if( nfoundleps != 1 ) cout << "ERROR! expected 1 lepton, found " << nfoundleps << endl;
	}

	//-----------------------------------------------------	
	// store both gen lepton info
	//-----------------------------------------------------
	
	else if( nleps_ == 2 ){
	  
	  float drmin      = 9999;
	  int   igenmin    =   -1;
	  int   nfoundleps =    0;
	  
	  //-----------------------------------------------------
	  // find gen lepton closest to reco lepton
	  //-----------------------------------------------------
	  
	  for ( int igen = 0 ; igen < (int)cms2.genps_id().size() ; igen++ ) { 

	    int id = genps_id().at(igen);

	    if( !( abs(id)==11 || abs(id)==13 || abs(id)==15 ) ) continue;

	    nfoundleps++;

	    if( dRbetweenVectors( *lep1_ , genps_p4().at(igen) ) < drmin ){
	      drmin   = dRbetweenVectors( *lep1_ , genps_p4().at(igen) );
	      igenmin = igen;
	    }
	  }

	  if( nfoundleps != 2 ) cout << "ERROR! expected 2 leptons, found " << nfoundleps << endl;

	  //-----------------------------------------------------
	  // store info for closest gen lepton
	  //-----------------------------------------------------

	  mcid1_   = genps_id().at(igenmin);
	  mclep1_  = &genps_p4().at(igenmin);
	  mcdr1_   = dRbetweenVectors( *lep1_ , *mclep1_ );

	  if( abs(mcid1_)==15 ){
	    mcdecay1_ = 1;
	    
	    //variables to calculate the visible tau energy
	    LorentzVector taudvis; 
	    
	    for(unsigned int kk = 0; kk < cms2.genps_lepdaughter_id().at(igenmin).size(); kk++) {
	      int daughter = abs(cms2.genps_lepdaughter_id()[igenmin][kk]);
	      //4-vector of visible part of tau
	      if ( daughter != 12 && daughter != 14 && daughter != 16 ) {
		if (taudvis.pt()==0)
		  taudvis = genps_lepdaughter_p4()[igenmin][kk];
		else 
		  taudvis = genps_lepdaughter_p4()[igenmin][kk]+taudvis;
	      }
	      // get tau daughter pt
	      if( daughter == 11 || daughter == 13 || daughter == 211 || daughter == 321 ){
		if( genps_lepdaughter_p4()[igenmin][kk].pt() > mctaudpt1_ ) {
		  mctaudpt1_ = genps_lepdaughter_p4()[igenmin][kk].pt();
		  mctaudid1_ = cms2.genps_lepdaughter_id()[igenmin][kk];
		  mctaud1_   = &genps_lepdaughter_p4()[igenmin][kk];
		}
	      }
	      
	      if( daughter == 211 || daughter == 321 ) mcndec1_  ++;  // count charged hadrons
	      if( daughter == 130 || daughter == 310 ) mcndeckls1_  ++;  // count K_L and K_S
	      if( daughter ==  22 || daughter == 111 ) mcndecem1_  ++;  // count photons and pi0
	      if( daughter ==  12 || daughter ==  14 ) mcdecay1_ = 2; // check for nu_e or nu_mu 
	    }
	    mctaudvis1_ = &taudvis;
	  }

	  //-----------------------------------------------------
	  // find 2nd lepton
	  //-----------------------------------------------------

	  int igenmin2 = -1;

	  for ( int igen = 0 ; igen < (int)cms2.genps_id().size() ; igen++ ) { 

	    if( igen == igenmin ) continue; //skip closest lepton

	    int id = genps_id().at(igen);

	    if( !( abs(id)==11 || abs(id)==13 || abs(id)==15 ) ) continue;

	    igenmin2 = igen;

	    mcid2_   = id;
	    mclep2_  = &genps_p4().at(igen);
	    if( ngoodlep_ > 1 ) mcdr2_   = dRbetweenVectors( *lep2_ , *mclep2_ );

	    if( abs(id)==15 ){
	      mcdecay2_ = 1;
	      //variables to calculate the visible tau energy
	      LorentzVector taudvis; 
	      
	      for(unsigned int kk = 0; kk < cms2.genps_lepdaughter_id().at(igen).size(); kk++) {
		int daughter = abs(cms2.genps_lepdaughter_id()[igen][kk]);


		//4-vector of visible part of tau
		if ( daughter != 12 && daughter != 14 && daughter != 16 ) {
		  if (taudvis.pt()==0)
		    taudvis = genps_lepdaughter_p4()[igen][kk];
		  else 
		    taudvis = genps_lepdaughter_p4()[igen][kk]+taudvis;
		}
		// get tau daughter pt
		if( daughter == 11 || daughter == 13 || daughter == 211 || daughter == 321 ){
		  if( genps_lepdaughter_p4()[igen][kk].pt() > mctaudpt1_ ) {
		    mctaudpt2_ = genps_lepdaughter_p4()[igen][kk].pt();
		    mctaudid2_ = cms2.genps_lepdaughter_id()[igen][kk];
		    mctaud2_   = &genps_lepdaughter_p4()[igen][kk];
		  }
		}

		if( daughter == 211 || daughter == 321 ) mcndec2_  ++;  // count charged hadrons
		if( daughter == 130 || daughter == 310 ) mcndeckls2_  ++;  // count K_L and K_S
		if( daughter ==  22 || daughter == 111 ) mcndecem2_  ++;  // count photons and pi0
		if( daughter == 12  || daughter == 14  ) mcdecay2_ = 2; // check for nu_e or nu_mu
	      }//end loop over daughters
	      mctaudvis2_ = &taudvis; 
	    } 
	  } 

	  if( igenmin2 < 0 ) cout << __FILE__ << " " << __LINE__ << " Error! unable to find 2nd gen lepton" << endl;

	  //-----------------------------------------------------
	  // find reco lepton corresponding to 2nd gen lepton
	  // check if found, if pass ID and/or iso
	  //-----------------------------------------------------

	  int  nMatchLeptons =  0;
	  int  imatch        = -1;
	  int  ID            = -1;

	  float drminlep = 999;

	  for( unsigned int iel = 0 ; iel < els_p4().size(); ++iel ){
	    int mc3idx = els_mc3idx().at(iel);

	    if( mc3idx != igenmin2 ) continue;
	    nMatchLeptons++;

	    float dr = els_mc3dr().at(iel);

	    if( dr < drminlep ){
	      drminlep   = dr;
	      imatch     = iel;
	      ID         = 1;
	    }
	  }
          
	  for( unsigned int imu = 0 ; imu < mus_p4().size(); ++imu ){
	    int mc3idx = mus_mc3idx().at(imu);

	    if( mc3idx != igenmin2 ) continue;
	    nMatchLeptons++;

	    float dr = mus_mc3dr().at(imu);

	    if( dr < drminlep ){
	      drminlep   = dr;
	      imatch     = imu;
	      ID         = 2;
	    }

	  }

	  mlepid_       = -1;
	  mlep_         =  0;
	  mleppassid_   = -1;
	  mleppassiso_  = -1;
	  mlepiso_      = -1.0;
	  mlepdr_       = -1.0;

	  if( nMatchLeptons > 0 ){

	    // found matched electron
	    if( ID == 1 ){
	      mlepid_       = 11 * els_charge().at(imatch);
	      mlep_         = &els_p4().at(imatch);
	      mleppassid_   = passElectronSelection_Stop2012_v2_NoIso( imatch , vetoTransition,vetoTransition,useOldIsolation) ? 1 : 0;
	      mleppassiso_  = passElectronSelection_Stop2012_v2_Iso  ( imatch , vetoTransition,vetoTransition,useOldIsolation) ? 1 : 0;
	      mlepiso_      = electronIsoValuePF2012_FastJetEffArea( imatch , 0.3 , 0 );
	    }

	    // found matched muon
	    else if( ID == 2 ){
	      mlepid_       = 13 * mus_charge().at(imatch);
	      mlep_         = &mus_p4().at(imatch);
	      mleppassid_   = muonIdNotIsolated( imatch , ZMet2012_v1 )   ? 1 : 0;
	      mleppassiso_  = muonIsoValuePF2012_deltaBeta(imatch) < 0.15 ? 1 : 0;
	      mlepiso_      = muonIsoValuePF2012_deltaBeta(imatch);
	    }

	    mlepdr_ = dRbetweenVectors( *mlep_ , *mclep2_ );
	  }
	  
	}

	else if( nleps_ < 0 || nleps_ > 2 ){
	  cout << "ERROR nleptons = " << nleps_ << endl;
	}

      }

      for (unsigned int ipf = 0; ipf < cms2.pfcands_p4().size(); ipf++) {

	if( pfcands_charge().at(ipf) == 0   ) continue;

 	int itrk = cms2.pfcands_trkidx().at(ipf);
	
 	if( itrk < (int)trks_trk_p4().size() && itrk >= 0 ){
 	  if( fabs( dz_trk_vtx(itrk,0) ) > 0.2 ){
 	    fillOverFlow( h_PU_trkpt , pfcands_p4().at(ipf).pt() );
 	  }
 	}
      }

      //------------------------------------------------------
      // store closest pf cand information for 2nd lepton
      //------------------------------------------------------

      pflepdr_       =  999.;
      pflepiso_      = -999.;
      pfleppt_       = -999.;
      pftauddr_      =  999.;
      pftaudiso_     = -999.;
      pftaudpt_      = -999.;
      pflepmindrj_   = 9999.;
      pftaudmindrj_  = 9999.;
      pflep_         = 0;
      pftaud_        = 0;
      if ( nleps_ == 2 )  {
	float pflepdr_out  =  999.;
	float pflepiso_out = -999.;
	float pfleppt_out  = -999.;
	float pftauddr_out  =  999.;
	float pftaudiso_out = -999.;
	float pftaudpt_out  = -999.;
	LorentzVector *pflep_out = 0;
	LorentzVector *pftaud_out = 0;
	for (unsigned int ipf = 0; ipf < cms2.pfcands_p4().size(); ipf++) {
	  
	  //	  if( pfcands_p4().at(ipf).pt() < 5   ) continue;
	  if( pfcands_charge().at(ipf) == 0   ) continue;
	  
	  int itrk = cms2.pfcands_trkidx().at(ipf);
	  
	  if( itrk < (int)trks_trk_p4().size() && itrk >= 0 ){
	    if( fabs( dz_trk_vtx(itrk,0) ) > 0.2 ) continue;
	  }
	  //Only remove leading lepton to see what happens to the sub-leading lepton
	  // bool isGoodLepton = false;
	  // for( int ilep = 0 ; ilep < goodLeptons.size() ; ilep++ ){
	  //   if( dRbetweenVectors( pfcands_p4().at(ipf) , goodLeptons.at(ilep) ) < 0.1 ) 
	  //     isGoodLepton = true;  
	  // }
	  // if( isGoodLepton ) continue;
	  if( dRbetweenVectors( pfcands_p4().at(ipf) , goodLeptons.at(imaxpt) ) < 0.1 ) continue;

	  //Store highest pT track within match radius or if none is found closest pfcand
	  float matchR = 0.15;
	  float drpf = dRbetweenVectors( pfcands_p4().at(ipf) , *mclep2_ );
	  float iso = trackIso(ipf) / pfcands_p4().at(ipf).pt();
	  if ( drpf < matchR && pfcands_p4().at(ipf).pt() > pfleppt_ ) {
	    pflepdr_ = drpf;
	    pfleppt_ = pfcands_p4().at(ipf).pt();
	    pflepiso_ = iso; 
	    pflep_ = &pfcands_p4().at(ipf);
	  } else if ( drpf > matchR && drpf < pflepdr_out ) {
	    pflepdr_out = drpf;
	    pfleppt_out = pfcands_p4().at(ipf).pt();
	    pflepiso_out = iso;
	    pflep_out = &pfcands_p4().at(ipf);
	  }
	  //check for tau decay and store information for daughter
	  if (mctaudid2_==-1) continue;
	  float taudrpf = dRbetweenVectors( pfcands_p4().at(ipf) , *mctaud2_ );
	  float tauiso = trackIso(ipf) / pfcands_p4().at(ipf).pt();
	  if ( taudrpf < matchR && pfcands_p4().at(ipf).pt() > pftaudpt_ ) {
	    pftauddr_ = taudrpf;
	    pftaudpt_ = pfcands_p4().at(ipf).pt();
	    pftaudiso_ = tauiso;
	    pftaud_ = &pfcands_p4().at(ipf);
	  } else if ( taudrpf > matchR && taudrpf < pftauddr_out ) {
	    pftauddr_out = taudrpf;
	    pftaudpt_out = pfcands_p4().at(ipf).pt();
	    pftaudiso_out = tauiso;
	    pftaud_out = &pfcands_p4().at(ipf);
	  }
	}// end loop over pf cands

	//If no match in cone is found, store closest track
	if (pfleppt_<0) {
	  pflepdr_ = pflepdr_out;
	  pfleppt_ = pfleppt_out;
	  pflepiso_ = pflepiso_out;
	  pflep_ = pflep_out;
	}
	if (pftaudpt_<0) {
	  pftauddr_ = pftauddr_out;
	  pftaudpt_ = pftaudpt_out;
	  pftaudiso_ = pftaudiso_out;
	  pftaud_ = pftaud_out;
	}
      }// end check for second lepton

      //------------------------------------------------------
      // track isolation variable definition
      //------------------------------------------------------
      float dz_cut = 0.05;
      float dz_cut_loose = 0.2;

      //------------------------------------------------------
      // store pt and iso for most isolated track (pt>10 GeV)
      //------------------------------------------------------

      trkpt10_           = -1.0;
      trkreliso10_       = 1000.;
      mleptrk10_         = -1.0;
      float miniso10     = 999;
      trkpt10loose_      = -1.0;
      trkreliso10loose_  = 1000.;

      for (unsigned int ipf = 0; ipf < cms2.pfcands_p4().size(); ipf++) {

	if( pfcands_p4().at(ipf).pt() < 10  ) continue;
	if( pfcands_charge().at(ipf) == 0   ) continue;

 	int itrk = cms2.pfcands_trkidx().at(ipf);
	
	bool isGoodLepton = false;
	for( int ilep = 0 ; ilep < (int)goodLeptons.size() ; ilep++ ){
	  if( dRbetweenVectors( pfcands_p4().at(ipf) , goodLeptons.at(ilep) ) < 0.1 ) 
	    isGoodLepton = true;  
	}
	bool isLeadLepton = ( dRbetweenVectors( pfcands_p4().at(ipf) , 
						goodLeptons.at(imaxpt) ) < 0.1 ) ? true : false;

	//store loose definition to compare with previous results
	float iso = trackIso(ipf, 0.3, dz_cut_loose, true) / pfcands_p4().at(ipf).pt();

 	if( itrk < (int)trks_trk_p4().size() && itrk >= 0 ){
 	  if( fabs( dz_trk_vtx(itrk,0) ) > dz_cut_loose ) continue;
 	}

	if( iso < trkreliso10loose_ && !isGoodLepton ){
	  trkpt10loose_       = pfcands_p4().at(ipf).pt();
	  trkreliso10loose_   = iso;
	}

	//tighten dz cut
 	if( itrk < (int)trks_trk_p4().size() && itrk >= 0 ){
 	  if( fabs( dz_trk_vtx(itrk,0) ) > dz_cut ) continue;
 	}

	//recalculated definition of the isolation
	iso = trackIso(ipf) / pfcands_p4().at(ipf).pt();

	if( iso < miniso10 && !isGoodLepton ){
	  miniso10       = iso;
	  trkpt10_       = pfcands_p4().at(ipf).pt();
	  mleptrk10_     = (*lep1_+pfcands_p4().at(ipf)).pt();
	  trkreliso10_   = iso;
	}

	if( iso < pfcandiso10_ && !isLeadLepton ){
	  pfcandiso10_ = iso;
	  pfcandpt10_ = pfcands_p4().at(ipf).pt();
	  pfcand10_ = &pfcands_p4().at(ipf);
	}
	
	//add all the variables with various pt thresholds

	if (isLeadLepton) continue;

	float iso0p1 = trackIso(ipf, 0.3, dz_cut, false, 0.1) / pfcands_p4().at(ipf).pt();
	if( iso0p1 < trkreliso10pt0p1_ && !isGoodLepton ){
	  trkpt10pt0p1_       = pfcands_p4().at(ipf).pt();
	  trkreliso10pt0p1_   = iso0p1;
	}
	if( iso0p1 < pfcandiso10pt0p1_ ){
	  pfcandpt10pt0p1_ = pfcands_p4().at(ipf).pt();
	  pfcandiso10pt0p1_ = iso0p1;
	}

	float iso0p2 = trackIso(ipf, 0.3, dz_cut, false, 0.2) / pfcands_p4().at(ipf).pt();
	if( iso0p2 < trkreliso10pt0p2_ && !isGoodLepton ){
	  trkpt10pt0p2_       = pfcands_p4().at(ipf).pt();
	  trkreliso10pt0p2_   = iso0p2;
	}
	if( iso0p2 < pfcandiso10pt0p2_ ){
	  pfcandpt10pt0p2_ = pfcands_p4().at(ipf).pt();
	  pfcandiso10pt0p2_ = iso0p2;
	}

	float iso0p3 = trackIso(ipf, 0.3, dz_cut, false, 0.3) / pfcands_p4().at(ipf).pt();
	if( iso0p3 < trkreliso10pt0p3_ && !isGoodLepton ){
	  trkpt10pt0p3_       = pfcands_p4().at(ipf).pt();
	  trkreliso10pt0p3_   = iso0p3;
	}
	if( iso0p3 < pfcandiso10pt0p3_ ){
	  pfcandpt10pt0p3_ = pfcands_p4().at(ipf).pt();
	  pfcandiso10pt0p3_ = iso0p3;
	}

	float iso0p4 = trackIso(ipf, 0.3, dz_cut, false, 0.4) / pfcands_p4().at(ipf).pt();
	if( iso0p4 < trkreliso10pt0p4_ && !isGoodLepton ){
	  trkpt10pt0p4_       = pfcands_p4().at(ipf).pt();
	  trkreliso10pt0p4_   = iso0p4;
	}
	if( iso0p4 < pfcandiso10pt0p4_ ){
	  pfcandpt10pt0p4_ = pfcands_p4().at(ipf).pt();
	  pfcandiso10pt0p4_ = iso0p4;
	}

	float iso0p5 = trackIso(ipf, 0.3, dz_cut, false, 0.5) / pfcands_p4().at(ipf).pt();
	if( iso0p5 < trkreliso10pt0p5_ && !isGoodLepton ){
	  trkpt10pt0p5_       = pfcands_p4().at(ipf).pt();
	  trkreliso10pt0p5_   = iso0p5;
	}
	if( iso0p5 < pfcandiso10pt0p5_ ){
	  pfcandpt10pt0p5_ = pfcands_p4().at(ipf).pt();
	  pfcandiso10pt0p5_ = iso0p5;
	}

	float iso0p6 = trackIso(ipf, 0.3, dz_cut, false, 0.6) / pfcands_p4().at(ipf).pt();
	if( iso0p6 < trkreliso10pt0p6_ && !isGoodLepton ){
	  trkpt10pt0p6_       = pfcands_p4().at(ipf).pt();
	  trkreliso10pt0p6_   = iso0p6;
	}
	if( iso0p6 < pfcandiso10pt0p6_ ){
	  pfcandpt10pt0p6_ = pfcands_p4().at(ipf).pt();
	  pfcandiso10pt0p6_ = iso0p6;
	}

	float iso0p7 = trackIso(ipf, 0.3, dz_cut, false, 0.7) / pfcands_p4().at(ipf).pt();
	if( iso0p7 < trkreliso10pt0p7_ && !isGoodLepton ){
	  trkpt10pt0p7_       = pfcands_p4().at(ipf).pt();
	  trkreliso10pt0p7_   = iso0p7;
	}
	if( iso0p7 < pfcandiso10pt0p7_ ){
	  pfcandpt10pt0p7_ = pfcands_p4().at(ipf).pt();
	  pfcandiso10pt0p7_ = iso0p7;
	}

	float iso0p8 = trackIso(ipf, 0.3, dz_cut, false, 0.8) / pfcands_p4().at(ipf).pt();
	if( iso0p8 < trkreliso10pt0p8_ && !isGoodLepton ){
	  trkpt10pt0p8_       = pfcands_p4().at(ipf).pt();
	  trkreliso10pt0p8_   = iso0p8;
	}
	if( iso0p8 < pfcandiso10pt0p8_ ){
	  pfcandpt10pt0p8_ = pfcands_p4().at(ipf).pt();
	  pfcandiso10pt0p8_ = iso0p8;
	}

	float iso0p9 = trackIso(ipf, 0.3, dz_cut, false, 0.9) / pfcands_p4().at(ipf).pt();
	if( iso0p9 < trkreliso10pt0p9_ && !isGoodLepton ){
	  trkpt10pt0p9_       = pfcands_p4().at(ipf).pt();
	  trkreliso10pt0p9_   = iso0p9;
	}
	if( iso0p9 < pfcandiso10pt0p9_ ){
	  pfcandpt10pt0p9_ = pfcands_p4().at(ipf).pt();
	  pfcandiso10pt0p9_ = iso0p9;
	}

	float iso1p0 = trackIso(ipf, 0.3, dz_cut, false, 1.0) / pfcands_p4().at(ipf).pt();
	if( iso1p0 < trkreliso10pt1p0_ && !isGoodLepton ){
	  trkpt10pt1p0_       = pfcands_p4().at(ipf).pt();
	  trkreliso10pt1p0_   = iso1p0;
	}
	if( iso1p0 < pfcandiso10pt1p0_ ){
	  pfcandpt10pt1p0_ = pfcands_p4().at(ipf).pt();
	  pfcandiso10pt1p0_ = iso1p0;
	}


      }

      //------------------------------------------------------
      // store pt and iso for most isolated track (pt>5 GeV)
      //------------------------------------------------------

      trkpt5_          = -1.0;
      trkreliso5_      = 1000.;
      mleptrk5_        = -1.0;
      float miniso5    = 999;
      trkpt5loose_     = -1.0;
      trkreliso5loose_ = 1000.;

      for (unsigned int ipf = 0; ipf < cms2.pfcands_p4().size(); ipf++) {

	if( pfcands_p4().at(ipf).pt() < 5   ) continue;
	if( pfcands_charge().at(ipf) == 0   ) continue;

 	int itrk = cms2.pfcands_trkidx().at(ipf);
	
	bool isGoodLepton = false;
	for( int ilep = 0 ; ilep < (int)goodLeptons.size() ; ilep++ ){
	  if( dRbetweenVectors( pfcands_p4().at(ipf) , goodLeptons.at(ilep) ) < 0.1 ) isGoodLepton = true;  
	}
	bool isLeadLepton = ( dRbetweenVectors( pfcands_p4().at(ipf) , goodLeptons.at(imaxpt) ) < 0.1 ) ? true : false;

 	if( itrk < (int)trks_trk_p4().size() && itrk >= 0 ){
 	  if( fabs( dz_trk_vtx(itrk,0) ) > dz_cut_loose ) continue;
 	}

	float iso = trackIso(ipf, 0.3, dz_cut_loose, true) / pfcands_p4().at(ipf).pt();

	if( iso < trkreliso5loose_ && !isGoodLepton ){
	  trkpt5loose_     = pfcands_p4().at(ipf).pt();
	  trkreliso5loose_ = iso;
	}

	//tighten dz cut
 	if( itrk < (int)trks_trk_p4().size() && itrk >= 0 ){
 	  if( fabs( dz_trk_vtx(itrk,0) ) > dz_cut ) continue;
 	}
 
	iso = trackIso(ipf) / pfcands_p4().at(ipf).pt();

	if( iso < miniso5 && !isGoodLepton ){
	  miniso5     = iso;
	  trkpt5_     = pfcands_p4().at(ipf).pt();
	  mleptrk5_   = (*lep1_+pfcands_p4().at(ipf)).pt();
	  trkreliso5_ = iso;
	  //itrk       = ipf;
	}

	if( iso < pfcandiso5_ && !isLeadLepton ){
	  pfcandiso5_ = iso;
	  pfcandpt5_ = pfcands_p4().at(ipf).pt();
	  pfcand5_ = &pfcands_p4().at(ipf);
	}

      }

      //----------------------------------------
      // nvertex variables
      //----------------------------------------

      nvtx_ = 0;
    
      for (size_t v = 0; v < cms2.vtxs_position().size(); ++v){
	if(isGoodVertex(v)) ++nvtx_;
      }

      npu_ = 0;
      npuMinusOne_ = 0;
      npuPlusOne_ = 0;
      if ( !isData ) {
        //Information for out-of-time PU
        for (unsigned int nbc=0;nbc<cms2.puInfo_nPUvertices().size();++nbc) {
	  if (cms2.puInfo_bunchCrossing().at(nbc)==0) npu_ = cms2.puInfo_nPUvertices().at(nbc);
	  else if (cms2.puInfo_bunchCrossing().at(nbc)==-1) npuMinusOne_ = cms2.puInfo_nPUvertices().at(nbc);
	  else if (cms2.puInfo_bunchCrossing().at(nbc)==+1) npuPlusOne_ = cms2.puInfo_nPUvertices().at(nbc);
	}
	//remove 3d-vtx reweighting for the moment, can do on the fly
	//	n3dvtxweight_ = weight3D( npuMinusOne_, npu_, npuPlusOne_ );
	n3dvtxweight_ = 0.;
      } else 
	n3dvtxweight_ = 1.;
      
      //----------------------------------------
      // PDF Information
      //----------------------------------------
      if ( !isData ) {
	pdfid1_ = int(cms2.pdfinfo_id1());
	pdfid2_ = int(cms2.pdfinfo_id2()); 
	pdfQ_   = cms2.pdfinfo_scale();
	pdfx1_  = cms2.pdfinfo_x1();
	pdfx2_  = cms2.pdfinfo_x2();
      }
      //-------------------------------------
      // jet counting
      //-------------------------------------

      VofP4 mediumBJets;

      int   imaxjet   = -1;
      float maxjetpt  = -1.;
      float ht_ = 0.;

      VofP4 vpfjets_p4;
      VofP4 vpfrawjets_p4;
      vpfjets_p4.clear();
      vpfrawjets_p4.clear();

      vector<float> fullcors;
      vector<float> l2l3cors;
      vector<float> rescors;
      vector<float> l1cors;
      fullcors.clear();
      l2l3cors.clear();
      rescors.clear();
      l1cors.clear();

      rhovor_ = cms2.evt_ww_rho_vor();

      float dmetx  = 0.0;
      float dmety  = 0.0;
      float jetptx = 0.0;
      float jetpty = 0.0;

      for (unsigned int ijet = 0 ; ijet < pfjets_p4().size() ; ijet++) {

	// skip jets with |eta| > 5.0
	if( fabs( pfjets_p4().at(ijet).eta() ) > 5.0 ) continue;

	// get L1FastL2L3Residual total correction
	jet_corrector_pfL1FastJetL2L3->setRho   ( cms2.evt_ww_rho_vor()           );
	jet_corrector_pfL1FastJetL2L3->setJetA  ( cms2.pfjets_area().at(ijet)     );
	jet_corrector_pfL1FastJetL2L3->setJetPt ( cms2.pfjets_p4().at(ijet).pt()  );
	jet_corrector_pfL1FastJetL2L3->setJetEta( cms2.pfjets_p4().at(ijet).eta() );
	double corr = jet_corrector_pfL1FastJetL2L3->getCorrection();

	// get L1Fast, L2, L3, Residual individual corrections
	jet_corrector_pfL1FastJetL2L3->setRho   ( cms2.evt_ww_rho_vor()           );
	jet_corrector_pfL1FastJetL2L3->setJetA  ( cms2.pfjets_area().at(ijet)     );
	jet_corrector_pfL1FastJetL2L3->setJetPt ( cms2.pfjets_p4().at(ijet).pt()  );
	jet_corrector_pfL1FastJetL2L3->setJetEta( cms2.pfjets_p4().at(ijet).eta() );
	vector<float> factors = jet_corrector_pfL1FastJetL2L3->getSubCorrections();

	// get residual correction only
	float rescorr = 1;
	if( isData ){
	  if( factors.size() == 4 ) rescorr = factors.at(3) / factors.at(2);
	  else                      cout << "ERROR! " << factors.size() << " jetSubCorrections" << endl;
	}

	LorentzVector vjet      = corr    * pfjets_p4().at(ijet);

	//---------------------------------------------------------------------------
	// get JES uncertainty
	//---------------------------------------------------------------------------
	
	pfUncertainty->setJetEta(vjet.eta());
	pfUncertainty->setJetPt(vjet.pt());   // here you must use the CORRECTED jet pt
	double unc = pfUncertainty->getUncertainty(true);

	LorentzVector vjetUp   = corr * pfjets_p4().at(ijet) * ( 1 + unc );
	LorentzVector vjetDown = corr * pfjets_p4().at(ijet) * ( 1 - unc );

	//LorentzVector vjetUp    = corr    * pfjets_p4().at(ijet) * 1.075; // over-estimate...
	//LorentzVector vjetDown  = corr    * pfjets_p4().at(ijet) * 0.925; // over-estimate...

	// lepton-jet overlap removal
	bool rejectJet = false;
	for( int ilep = 0 ; ilep < (int)goodLeptons.size() ; ilep++ ){
	  if( dRbetweenVectors( vjet , goodLeptons.at(ilep) ) < 0.4 ) rejectJet = true;  
	}
	if( rejectJet ) continue;
	          
	// PFJetID
	if( !passesPFJetID(ijet) ){
	  jetid_ = 0;
	  if( vjet.pt() > 30 && fabs( vjet.eta() ) < 2.5 ) jetid30_ = 0;
	  continue;
	}
 
        // store raw pfjet p4's and corrections for type1 pfmet
        // using corr jet pT > 10 GeV and adding |eta|<4.7 as in AN2011/459 to avoid problems with
        // erroneously large corrections
        if( vjet.pt() > 10 && fabs(vjet.eta()) < 4.7 ){
	  float l1cor = factors.at(0);
          vpfrawjets_p4.push_back( pfjets_p4().at(ijet) );
          fullcors.push_back( corr );
	  l2l3cors.push_back( corr / l1cor );
          rescors.push_back( rescorr );
          l1cors.push_back( l1cor );
	  // cout<<"l1corr: "<<l1cor<<" (rho: "<<rhovor_
	  //     <<", nvts: "<<ndavtx_
	  //     <<") l2l3resi: "<<corr/l1cor
	  //     <<" all: "<<corr
	  //     <<" pt: "<<vjet.pt()<<" eta: "<<vjet.eta()
	  //     <<endl;
        }

	//------------------------------------------------------------------------------------------------------------
	// MET correction quantities
	// here we store 2 quantities:
	// the delta(METx/y) you get by varying the jet with pT > 10 GeV by their uncertainties (dmetx,dmety
	// the vector sum of pT > 10 GeV selected jets (jetptx,jetpty) --> use this to calculate unclustered energy
	//------------------------------------------------------------------------------------------------------------

	if( vjet.pt() > 10 ){
	  dmetx  += vjetUp.px() - vjet.px();
	  dmety  += vjetUp.py() - vjet.py();
	  jetptx += vjet.px();
	  jetpty += vjet.py();
	}

	// store L1FastL2L3Residual jet p4's pt > 15 GeV
	if( vjet.pt() > 15 && fabs( vjet.eta() ) < 2.5 ){
	  vpfjets_p4.push_back( vjet );
	}

	// njets JEC up
	if( vjetUp.pt() > 30. && fabs( vjetUp.eta() ) < 2.5 ){
	  njetsUp_++;
	  htUp_ += vjetUp.pt();
	}

	// njets JEC down
	if( vjetDown.pt() > 30. && fabs( vjetDown.eta() ) < 2.5 ){
	  njetsDown_++;
	  htDown_ += vjetDown.pt();
	}

	// njets: L1FastL2L3Residual, pt > 30 GeV
	if( vjet.pt() > 30 && fabs( vjet.eta() ) < 2.5 ){
	  npfjets30_ ++;
	  //count jets that are not overlapping with second lepton
	  npfjets30lepcorr_ ++;
	  if (nleps_==2 &&  mclep2_->Pt() > 30.
	      && dRbetweenVectors(*mclep2_, vjet) < 0.4 ) 
	    npfjets30lepcorr_ --;
	  htpf30_ += vjet.pt();
	}

	// njets: L1FastL2L3Residual, pt > 35 GeV
	if( vjet.pt() > 35 && fabs( vjet.eta() ) < 2.5 ){
	  npfjets35_ ++;
	  htpf35_ += vjet.pt();
	}

	// njets: L1FastL2L3Residual, pt > 40 GeV
	if( vjet.pt() > 40 && fabs( vjet.eta() ) < 2.5 ){
	  npfjets40_ ++;
	  htpf40_ += vjet.pt();
	}

	// njets: L1FastL2L3Residual, pt > 45 GeV
	if( vjet.pt() > 45 && fabs( vjet.eta() ) < 2.5 ){
	  npfjets45_ ++;
	  htpf45_ += vjet.pt();
	}

	ht_ += vjet.pt();

	if(       vjet.pt()    < 30. )           continue;
	if( fabs( vjet.eta() ) > 2.5 )           continue;

	
	//-------------------------------------
	// b-tag counting
	//-------------------------------------
	
	// btag variables: SSV
	float discrimssv = pfjets_simpleSecondaryVertexHighEffBJetTag().at(ijet);
	bool isbtagssv = ( discrimssv > 1.74 ) ? true : false;
	if (isbtagssv)     {
	  nbtagsssv_++;
	  mediumBJets.push_back(vjet);
	}

	// btag variables: TCHEL
	float discrimtche = pfjets_trackCountingHighEffBJetTag().at(ijet);
	bool isbtagtcl = ( discrimtche > 1.7 ) ? true: false;
	if (isbtagtcl)     nbtagstcl_++;

	// btag variables: TCHEM
	bool isbtagtcm = ( discrimtche > 3.3 ) ? true: false;
	if (isbtagtcm)     nbtagstcm_++;

	// btag variables: CSVL
	float discrimcsv = pfjets_combinedSecondaryVertexBJetTag().at(ijet);
	bool isbtagcsvl = ( discrimcsv > 0.244 ) ? true: false;
	if (isbtagcsvl)     nbtagscsvl_++;

	// btag variables: CSVM
	bool isbtagcsvm = ( discrimcsv > 0.679 ) ? true: false;
	if (isbtagcsvm)     nbtagscsvm_++;

	// btag variables: CSVT
	bool isbtagcsvt = ( discrimcsv > 0.898 ) ? true: false;
	if (isbtagcsvt)     nbtagscsvt_++;

	// in MC apply b-tagging corrections 
	if ( !isData ) {

	  //set seed for random number generator based on event number and jet phi
	  int randseed = evt_event()+(int)(vjet.phi()*1000.);
	  random3_->SetSeed(randseed);
	  float rand = random3_->Uniform(1.);  
	  
	  bool isbmatched = isGenBMatched(vjet, 0.5);
	  int pdgid = isbmatched ? 5 : 0;

	  // btag variables: SSV
	  float SFb_ssv  = getBtagSF(  vjet.pt(), vjet.eta(), "SSVHEM");
	  float SFl_ssv  = getMistagSF(vjet.pt(), vjet.eta(), "SSVHEM");
	  float Effl_ssv = getMistags( vjet.pt(), vjet.eta(), "SSVHEM");
	  bool iscorrbtagssv = getCorrBtag(isbtagssv, pdgid, SFb_ssv, SFl_ssv, Effl_ssv, rand);
	  if (iscorrbtagssv) nbtagsssvcorr_++;
	  
	  // btag variables: TCHEL
	  float SFb_tcl  = getBtagSF(  vjet.pt(), vjet.eta(), "TCHEL");
	  float SFl_tcl  = getMistagSF(vjet.pt(), vjet.eta(), "TCHEL");
	  float Effl_tcl = getMistags( vjet.pt(), vjet.eta(), "TCHEL");
	  bool iscorrbtagtcl = getCorrBtag(isbtagtcl, pdgid, SFb_tcl, SFl_tcl, Effl_tcl, rand);
	  if (iscorrbtagtcl) nbtagstclcorr_++;

	  // btag variables: TCHEM
	  float SFb_tcm  = getBtagSF(  vjet.pt(), vjet.eta(), "TCHEM");
	  float SFl_tcm  = getMistagSF(vjet.pt(), vjet.eta(), "TCHEM");
	  float Effl_tcm = getMistags( vjet.pt(), vjet.eta(), "TCHEM");
	  bool iscorrbtagtcm = getCorrBtag(isbtagtcm, pdgid, SFb_tcm, SFl_tcm, Effl_tcm, rand);
	  if (iscorrbtagtcm) nbtagstcmcorr_++;

	  // btag variables: CSVL
	  float SFb_csvl  = getBtagSF(  vjet.pt(), vjet.eta(), "CSVL");
	  float SFl_csvl  = getMistagSF(vjet.pt(), vjet.eta(), "CSVL");
	  float Effl_csvl = getMistags( vjet.pt(), vjet.eta(), "CSVL");
	  bool iscorrbtagcsvl = getCorrBtag(isbtagcsvl, pdgid, SFb_csvl, SFl_csvl, Effl_csvl, rand);
	  if (iscorrbtagcsvl) nbtagscsvlcorr_++;

	  // btag variables: CSVM
	  float SFb_csvm  = getBtagSF(  vjet.pt(), vjet.eta(), "CSVM");
	  float SFl_csvm  = getMistagSF(vjet.pt(), vjet.eta(), "CSVM");
	  float Effl_csvm = getMistags( vjet.pt(), vjet.eta(), "CSVM");
	  bool iscorrbtagcsvm = getCorrBtag(isbtagcsvm, pdgid, SFb_csvm, SFl_csvm, Effl_csvm, rand);
	  if (iscorrbtagcsvm) nbtagscsvmcorr_++;
	  
	  // btag variables: CSVT
	  float SFb_csvt  = getBtagSF(  vjet.pt(), vjet.eta(), "CSVT");
	  float SFl_csvt  = getMistagSF(vjet.pt(), vjet.eta(), "CSVT");
	  float Effl_csvt = getMistags( vjet.pt(), vjet.eta(), "CSVT");
	  bool iscorrbtagcsvt = getCorrBtag(isbtagcsvt, pdgid, SFb_csvt, SFl_csvt, Effl_csvt, rand);
	  if (iscorrbtagcsvt) nbtagscsvtcorr_++;

	} 
	
	// store max jet pt
	if( vjet.pt() > maxjetpt ){
	  maxjetpt = vjet.pt();
	  imaxjet  = ijet;
	}
      }

      // no b-tagging corrections in data, so make counts same
      if (isData) {
	nbtagsssvcorr_  = nbtagsssv_;
	nbtagstclcorr_  = nbtagstcl_;
	nbtagstcmcorr_  = nbtagstcm_;
	nbtagscsvlcorr_ = nbtagscsvl_;
	nbtagscsvmcorr_ = nbtagscsvm_;
	nbtagscsvtcorr_ = nbtagscsvt_;
      }

      //store njets kscaling - for ttbar dilepton
      //K3 = 0.92 pm 0.03 and K4 = 0.83 pm 0.03
     if( TString(prefix).Contains("ttall") && nleps_==2 && npfjets30_>3){
       if (npfjets30lepcorr_==3)      knjets_=0.92;
       else if (npfjets30lepcorr_>=4) knjets_=0.83;
     }

      // type1 met's
      pair<float, float> p_t1met10     = Type1PFMET( vpfrawjets_p4 , fullcors , l1cors , 10.0 );
      pair<float, float> p_t1met20     = Type1PFMET( vpfrawjets_p4 , fullcors , l1cors , 20.0 );
      pair<float, float> p_t1met30     = Type1PFMET( vpfrawjets_p4 , fullcors , l1cors , 30.0 );
      t1met10_        = p_t1met10.first;
      t1met20_        = p_t1met20.first;
      t1met30_        = p_t1met30.first;
      t1met10phi_     = p_t1met10.second;
      t1met20phi_     = p_t1met20.second;	  
      t1met30phi_     = p_t1met30.second;	  

      //phi-corrected type1 met
      pair<float, float> p_t1metphicorr = 
	getPhiCorrMET( t1met10_, t1met10phi_, evt_pfsumet(), !isData, true);
      t1metphicorr_    = p_t1metphicorr.first;
      t1metphicorrphi_ = p_t1metphicorr.second;

      //---------------------------------------
      // now calculate METup and METdown
      //---------------------------------------

      float pfmetx = t1metphicorr_ * cos( t1metphicorrphi_ );
      float pfmety = t1metphicorr_ * sin( t1metphicorrphi_ );

      //--------------------------------------------------------
      // calculate unclustered energy x and y components
      // unclustered energy = -1 X ( MET + jets + leptons )
      //--------------------------------------------------------

      float unclustered_x = -1 * ( pfmetx + jetptx );
      float unclustered_y = -1 * ( pfmety + jetpty );

      for( unsigned int ilep = 0 ; ilep < goodLeptons.size() ; ilep++ ){
	unclustered_x -= goodLeptons.at(ilep).px();
	unclustered_y -= goodLeptons.at(ilep).py();
      }
      
      //------------------------------------------------------------------------------
      // now vary jets according to JEC uncertainty, vary unclustered energy by 10%
      //------------------------------------------------------------------------------

      float pfmetx_up = pfmetx - dmetx - 0.1 * unclustered_x; 
      float pfmety_up = pfmety - dmety - 0.1 * unclustered_y; 

      // pfmet DOWN
      t1metphicorrup_    = sqrt( pfmetx_up * pfmetx_up + pfmety_up * pfmety_up );
      t1metphicorrphiup_ = atan2( pfmety_up , pfmetx_up );

      float pfmetx_dn = pfmetx + dmetx + 0.1 * unclustered_x; 
      float pfmety_dn = pfmety + dmety + 0.1 * unclustered_y; 

      // pfmet UP
      t1metphicorrdn_    = sqrt( pfmetx_dn * pfmetx_dn + pfmety_dn * pfmety_dn );
      t1metphicorrphidn_ = atan2( pfmety_dn , pfmetx_dn );

      //-------------------------------------------------------------------------------
      // for dilepton events, calculate MET where negative lepton is added back to MET
      //-------------------------------------------------------------------------------

      t1metphicorrlep_    = 0;
      t1metphicorrlepphi_ = 0;

      if( ngoodlep_ > 1 && lepm_->pt() > 0 ){
	float metx = t1metphicorr_ * cos( t1metphicorrphi_ );
	float mety = t1metphicorr_ * sin( t1metphicorrphi_ );

	metx += lepm_->px();
	mety += lepm_->py();

	t1metphicorrlep_    = sqrt(metx*metx + mety*mety);
	t1metphicorrlepphi_ = atan2( mety , metx );
      }

      // store L1FastL2L3Residual pfjets
      // check if jet is b-tagged
      sort(vpfjets_p4.begin(), vpfjets_p4.end(), sortByPt);
      if( vpfjets_p4.size() > 0 ) {
	pfjet1_  = &vpfjets_p4.at(0);
	bjet1_ = isBTagged( vpfjets_p4.at(0), mediumBJets ) ? 1 : 0;
	if (!isData) {
	  lepjet1_ = getLeptonMatchIndex ( pfjet1_, mclep1_, mclep2_, 0.4 );
	  qgjet1_ = isGenQGMatched( vpfjets_p4.at(0), 0.4 );
	}
      }
      if( vpfjets_p4.size() > 1 ) {
	pfjet2_  = &vpfjets_p4.at(1);
	bjet2_ = isBTagged( vpfjets_p4.at(1), mediumBJets ) ? 1 : 0;
	if (!isData) {
	  lepjet2_ = getLeptonMatchIndex ( pfjet2_, mclep1_, mclep2_, 0.4 );
	  qgjet2_ = isGenQGMatched( vpfjets_p4.at(1), 0.4 );
	}
      }
      if( vpfjets_p4.size() > 2 ) {
	pfjet3_  = &vpfjets_p4.at(2);
	bjet3_ = isBTagged( vpfjets_p4.at(2), mediumBJets ) ? 1 : 0;
	if (!isData) {
	  lepjet3_ = getLeptonMatchIndex ( pfjet3_, mclep1_, mclep2_, 0.4 );
	  qgjet3_ = isGenQGMatched( vpfjets_p4.at(2), 0.4 );
	}
      }
      if( vpfjets_p4.size() > 3 ) {
	pfjet4_  = &vpfjets_p4.at(3);
	bjet4_ = isBTagged( vpfjets_p4.at(3), mediumBJets ) ? 1 : 0;
	if (!isData) {
	  lepjet4_ = getLeptonMatchIndex ( pfjet4_, mclep1_, mclep2_, 0.4 );
	  qgjet4_ = isGenQGMatched( vpfjets_p4.at(3), 0.4 );
	}
      }
      if( vpfjets_p4.size() > 4 ) {
	pfjet5_  = &vpfjets_p4.at(4);
	bjet5_ = isBTagged( vpfjets_p4.at(4), mediumBJets ) ? 1 : 0;
	if (!isData) {
	  lepjet5_ = getLeptonMatchIndex ( pfjet5_, mclep1_, mclep2_, 0.4 );
	  qgjet5_ = isGenQGMatched( vpfjets_p4.at(4), 0.4 );
	}
      }
      if( vpfjets_p4.size() > 5 ) {
	pfjet6_  = &vpfjets_p4.at(5);
	bjet6_ = isBTagged( vpfjets_p4.at(5), mediumBJets ) ? 1 : 0;
	if (!isData) {
	  lepjet6_ = getLeptonMatchIndex ( pfjet6_, mclep1_, mclep2_, 0.4 );
	  qgjet6_ = isGenQGMatched( vpfjets_p4.at(5), 0.4 );
	}
      }

      //store distance to closest jet for pfcand
      if ( nleps_ == 2 ) {
	pflepmindrj_  = getminjdr( vpfjets_p4, pflep_ );
	pftaudmindrj_ = getminjdr( vpfjets_p4, pftaud_ );
      }
      pfcandmindrj5_  = getminjdr( vpfjets_p4, pfcand5_ );
      pfcandmindrj10_ = getminjdr( vpfjets_p4, pfcand10_ );

      // max jet variables
      if( imaxjet > -1 ){ 
	LorentzVector vjetcorr = pfjets_corL1FastL2L3().at(imaxjet) * pfjets_p4().at(imaxjet);
	jet_ = &vjetcorr;

	LorentzVector vjetraw = pfjets_p4().at(imaxjet);

	ptjetraw_     = vjetraw.pt();
	ptjet23_      = pfjets_cor().at(imaxjet)           * vjetraw.pt();
	ptjetF23_     = pfjets_corL1FastL2L3().at(imaxjet) * vjetraw.pt();
	ptjetO23_     = pfjets_corL1L2L3().at(imaxjet)     * vjetraw.pt();
	//cosphijz_     = -1 * cos( vjetraw.phi() - hyp_p4()[hypIdx].phi() );
	
	LorentzVector vjet = pfjets_corL1FastL2L3().at(imaxjet) * pfjets_p4().at(imaxjet);
	dphijm_ = acos(cos(vjet.phi()-evt_pfmetPhi()));
      }

      emjet10_     = -1;
      emjet20_     = -1;

      //-------------------------------------------------------------
      // find jet with max EM fraction outside tracker acceptance
      //-------------------------------------------------------------

      for (unsigned int ijet = 0; ijet < jets_p4().size(); ijet++) {
	
	LorentzVector vjet = jets_p4().at(ijet) * jets_corL1FastL2L3().at(ijet);

	if( fabs( vjet.eta() ) < 2.5 )         continue;

	bool rejectJet = false;
	for( int ilep = 0 ; ilep < (int)goodLeptons.size() ; ilep++ ){
	  if( dRbetweenVectors( vjet , goodLeptons.at(ilep) ) < 0.4 ) rejectJet = true;  
	}
	if( rejectJet ) continue;

	float emfrac = jets_emFrac().at(ijet);

	if( vjet.pt() < 10. ) continue;

	if( emfrac > emjet10_ ) emjet10_ = emfrac;

	if( vjet.pt() < 20. ) continue;

	if( emfrac > emjet20_ ) emjet20_ = emfrac;

      }

      //--------------------------------
      // get non-isolated leptons
      //--------------------------------

      nonisoel_ = 0;
      nonisomu_ = 0;

      float maxelpt  = -1.;
      int   imaxelpt = -1;

      for( unsigned int iel = 0 ; iel < els_p4().size(); ++iel ){

	if( els_p4().at(iel).pt() < 10 )                                                                         continue;
	if( !passElectronSelection_Stop2012_v2_NoIso( iel , vetoTransition,vetoTransition,useOldIsolation))      continue;

	// don't count the leptons that we already counted as good
	bool isGoodLepton = false;
	for( int ilep = 0 ; ilep < (int)goodLeptons.size() ; ilep++ ){
	  if( dRbetweenVectors( els_p4().at(iel) , goodLeptons.at(ilep) ) < 0.1 ) isGoodLepton = true;  
	}
	if( isGoodLepton ) continue;

	// don't count leptons near b-jets (SSVM)
	bool nearBJet = false;
	for( int ijet = 0 ; ijet < (int)mediumBJets.size() ; ijet++ ){
	  if( dRbetweenVectors( els_p4().at(iel) , mediumBJets.at(ijet) ) < 0.4 ) nearBJet = true;
	}
	if( nearBJet ) continue;
	
	if( els_p4().at(iel).pt() > maxelpt ){
	  maxelpt  = els_p4().at(iel).pt();
	  imaxelpt = iel;
	}
      }

      if( imaxelpt >= 0 ) nonisoel_ = &(els_p4().at(imaxelpt));

      float maxmupt  = -1.;
      int   imaxmupt = -1;

      for( unsigned int imu = 0 ; imu < mus_p4().size(); ++imu ){

	if( mus_p4().at(imu).pt() < 10 )                   continue;
	if( !muonIdNotIsolated(imu , ZMet2012_v1 ) )       continue;

	// don't count the leptons that we already counted as good
	bool isGoodLepton = false;
	for( int ilep = 0 ; ilep < (int)goodLeptons.size() ; ilep++ ){
	  if( dRbetweenVectors( mus_p4().at(imu) , goodLeptons.at(ilep) ) < 0.1 ) isGoodLepton = true;  
	}
	if( isGoodLepton ) continue;

	// don't count leptons near b-jets SSVM)
	bool nearBJet = false;
	for( int ijet = 0 ; ijet < (int)mediumBJets.size() ; ijet++ ){
	  if( dRbetweenVectors( mus_p4().at(imu) , mediumBJets.at(ijet) ) < 0.4 ) nearBJet = true;
	}
	if( nearBJet ) continue;
	
	if( mus_p4().at(imu).pt() > maxmupt ){
	  maxmupt  = mus_p4().at(imu).pt();
	  imaxmupt = imu;
	}
      }

      if( imaxmupt >= 0 ) nonisomu_ = &(mus_p4().at(imaxmupt));

      //---------------------------------
      // jet mass variables
      //---------------------------------

      // mjj_ = -1;

      // if( nbctcm_ >= 2 && goodCaloJets.size() >= 2 ){
      // 	sort( goodCaloJets.begin(), goodCaloJets.end(), sortByPt);
      // 	mjj_ = ( goodCaloJets.at(0) + goodCaloJets.at(1) ).mass();
      // }

      //---------------------------------
      // L1offset jets
      //---------------------------------
      
      htoffset_    = 0.;
      njetsoffset_ = 0;

      for (unsigned int ijet = 0 ; ijet < pfjets_p4().size() ; ijet++) {
          
	LorentzVector vjet = pfjets_corL1L2L3().at(ijet) * pfjets_p4().at(ijet);

	bool rejectJet = false;
	for( int ilep = 0 ; ilep < (int)goodLeptons.size() ; ilep++ ){
	  if( dRbetweenVectors( vjet , goodLeptons.at(ilep) ) < 0.4 ) rejectJet = true;  
	}
	if( rejectJet ) continue;
          
	if( fabs( vjet.eta() ) > 2.5 )           continue;
	if( !passesPFJetID(ijet) )               continue;
	if( vjet.pt() < 30. )                    continue;

	njetsoffset_++;
	htoffset_ += vjet.pt();

      }

      //---------------------------------
      // uncorrected jets
      //---------------------------------

      htuncor_    = 0.;
      njetsuncor_ = 0;

      for (unsigned int ijet = 0 ; ijet < pfjets_p4().size() ; ijet++) {
          
	LorentzVector vjet = pfjets_cor().at(ijet) * pfjets_p4().at(ijet);

	bool rejectJet = false;
	for( int ilep = 0 ; ilep < (int)goodLeptons.size() ; ilep++ ){
	  if( dRbetweenVectors( vjet , goodLeptons.at(ilep) ) < 0.4 ) rejectJet = true;  
	}
	if( rejectJet ) continue;
          
	if( fabs( vjet.eta() ) > 2.5 )           continue;
	if( !passesPFJetID(ijet) )               continue;
	if( vjet.pt() < 15. )                    continue;

	if( vjet.pt() < 30. )                    continue;

	njetsuncor_++;
	htuncor_ += vjet.pt();

      }

      ngenjets_ = 0;
      htgen_    = 0;

      if( !isData ){
	for (unsigned int igjet = 0 ; igjet < genjets_p4().size() ; igjet++) {
	    
	  LorentzVector vgjet = genjets_p4().at(igjet);

	  bool rejectJet = false;
	  for( int ilep = 0 ; ilep < (int)goodLeptons.size() ; ilep++ ){
	    if( dRbetweenVectors( vgjet , goodLeptons.at(ilep) ) < 0.4 ) rejectJet = true;  
	  }
	  if( rejectJet ) continue;
	    
	  if( vgjet.pt() < 30.                   )  continue;
	  if( fabs( vgjet.eta() ) > 2.5          )  continue;
	    
	  ngenjets_++;
	  htgen_ += vgjet.pt();
	}
      }
     
      if( !isData ){
	genmet_     = gen_met();
	gensumet_   = gen_sumEt();
	genmetphi_  = gen_metPhi();

      }

      //----------------------------
      // MET flavors
      //----------------------------

      pair<float, float> p_met; //met and met phi
      //p_met = getMet( "tcMET"    , hypIdx);
      p_met = make_pair( evt_tcmet() , evt_tcmetPhi() );

      tcmet_    = p_met.first;
      tcmetphi_ = p_met.second;
      tcsumet_  = evt_tcsumet();

      pfmet_    = evt_pfmet();
      pfmetphi_ = evt_pfmetPhi();
      pfsumet_  = evt_pfsumet();

      //p_met = getMet( "pfMET"    , hypIdx);
      p_met = make_pair( evt_pfmet() , evt_pfmetPhi() );      

      // pair<float, float> pfmetUp   = ScaleMET( p_met , hyp_p4().at(hypIdx) , 1.075 );
      // pair<float, float> pfmetDown = ScaleMET( p_met , hyp_p4().at(hypIdx) , 0.925 );
      // pair<float, float> pfmetTest = ScaleMET( p_met , hyp_p4().at(hypIdx) , 1.000 );

      // pfmetUp_      = pfmetUp.first;
      // pfmetDown_    = pfmetDown.first;
      // pfmetTest_    = pfmetTest.first;

      meff_ = ht_ + pfmet_ + lep1_->pt();

      m_events.insert(pair<int,int>(evt_event(), 1));

      //---------------------------
      // set event weight
      //---------------------------

      weight_ = -1.;

      if( TString(prefix).Contains("T2") ){
	mG_ = -999; //sparm_mG();
	mL_ = -999; //sparm_mL();
	x_  = -999; //sparm_mf();

	weight_ = -999; //lumi * stopPairCrossSection(mG_) * (1000./50000.);

	if( doTenPercent )	  weight_ *= 10;
      }

      else if(strcmp(prefix,"LMscan") == 0){ 

	m0_  = -999; //sparm_m0();
	m12_ = -999; //sparm_m12();

	ksusy_     = kfactorSUSY(m0_,m12_,"tanbeta10");
	ksusyup_   = kfactorSUSY(m0_,m12_,"tanbeta10Scale20");
	ksusydn_   = kfactorSUSY(m0_,m12_,"tanbeta10Scale05");
	xsecsusy_  = cmssm_loxsec(m0_,m12_);
	xsecsusy2_ = getMsugraCrossSection(m0_,m12_,10);

	weight_ = lumi * ksusy_ * xsecsusy_ * (1000. / 10000.); // k * xsec / nevents

	if( doTenPercent )	  weight_ *= 10;
      }

      else if( isData ){
	weight_ = 1;
      }

      else{

	weight_ = kFactor * evt_scale1fb() * lumi;
        //do a signed weight for mcatnlo
        if ( TString(prefix).Contains("mcatnlo") && genps_weight()<0) weight_ *= -1.;

	if( doTenPercent )	  weight_ *= 10;

	if( TString(prefix).Contains("LM") ){
	  if( strcmp( prefix , "LM0" )  == 0 ) weight_ *= kfactorSUSY( "lm0" );
	  if( strcmp( prefix , "LM1" )  == 0 ) weight_ *= kfactorSUSY( "lm1" );
	  if( strcmp( prefix , "LM2" )  == 0 ) weight_ *= kfactorSUSY( "lm2" );
	  if( strcmp( prefix , "LM3" )  == 0 ) weight_ *= kfactorSUSY( "lm3" );
	  if( strcmp( prefix , "LM4" )  == 0 ) weight_ *= kfactorSUSY( "lm4" );
	  if( strcmp( prefix , "LM5" )  == 0 ) weight_ *= kfactorSUSY( "lm5" );
	  if( strcmp( prefix , "LM6" )  == 0 ) weight_ *= kfactorSUSY( "lm6" );
	  if( strcmp( prefix , "LM7" )  == 0 ) weight_ *= kfactorSUSY( "lm7" );
	  if( strcmp( prefix , "LM8" )  == 0 ) weight_ *= kfactorSUSY( "lm8" );
	  if( strcmp( prefix , "LM9" )  == 0 ) weight_ *= kfactorSUSY( "lm9" );
	  if( strcmp( prefix , "LM10" ) == 0 ) weight_ *= kfactorSUSY( "lm10");
	  if( strcmp( prefix , "LM11" ) == 0 ) weight_ *= kfactorSUSY( "lm11");
	  if( strcmp( prefix , "LM12" ) == 0 ) weight_ *= kfactorSUSY( "lm12");
	  if( strcmp( prefix , "LM13" ) == 0 ) weight_ *= kfactorSUSY( "lm13");
	}
      }


      //tranverse mass leading lepton & met
      dphilm_ = fabs( lep1_->phi() - pfmetphi_ );
      if( dphilm_ > TMath::Pi() ) dphilm_ = TMath::TwoPi() - dphilm_;

      mt_ = sqrt( 2 * ( lep1_->pt() * pfmet_ * (1 - cos( dphilm_ ) ) ) );

      //transverse mass for leading lepton & type1 mets
      t1met10mt_     = getMT( lep1_->pt() , lep1_->phi() , t1met10_     , t1met10phi_ );
      t1met20mt_     = getMT( lep1_->pt() , lep1_->phi() , t1met20_     , t1met20phi_ );
      t1met30mt_     = getMT( lep1_->pt() , lep1_->phi() , t1met30_     , t1met30phi_ );
      //phi-corrected met
      t1metphicorrmt_   = getMT( lep1_->pt() , lep1_->phi() , t1metphicorr_ , t1metphicorrphi_ );

      t1metphicorrmtup_ = getMT( lep1_->pt() , lep1_->phi() , t1metphicorrup_ , t1metphicorrphiup_ );
      t1metphicorrmtdn_ = getMT( lep1_->pt() , lep1_->phi() , t1metphicorrdn_ , t1metphicorrphidn_ );

      if( ngoodlep_ >1 && lepp_->pt() > 0 && lepm_->pt() > 0 ){
	t1metphicorrlepmt_ = getMT( lepp_->pt() , lepp_->phi() , t1metphicorrlep_ , t1metphicorrlepphi_ );
      }

      //pt of the leading lepton and met system
      lepmetpt_ = sqrt( pow((lep1_->px() + pfmet_*cos(pfmetphi_)),2) 
			+ pow((lep1_->py() + pfmet_*sin(pfmetphi_)),2) );
      lept1met10pt_ = sqrt( pow((lep1_->px() + t1met10_*cos(t1met10phi_)),2) 
			    + pow((lep1_->py() + t1met10_*sin(t1met10phi_)),2) );
      //dijet mass two bs highest pT b-tagged jets
      if (mediumBJets.size()>1) {
	mbb_ = (mediumBJets.at(0)+mediumBJets.at(1)).M();
      }

      //ugly bit of code to store information about the lepton, neutrino and W
      //for single lepton events
      if( !isData ){
	bool foundlep = false;
	bool foundnu  = false;
	for (int igen = (cms2.genps_p4().size()-1); igen >-1; igen--) {
	  int id = genps_id().at(igen);
	  if (abs(id)==11 || abs(id)==13) {
	    foundlep = true;
	    mclep_ = &cms2.genps_p4()[igen];
	  }
	  if (abs(id)==12 || abs(id)==14) {
	    foundnu = true;
	    mcnu_  = &cms2.genps_p4()[igen];
	  }
	}
	if (foundlep && foundnu) {
	  mcmln_ = (*mclep_+*mcnu_).mass();
	  mcmtln_ = getMT( mclep_->Pt() , mclep_->Phi() , mcnu_->Pt() , mcnu_->Phi() );
	}
      }

      /*
      //REPLACETOPMASS
      // calculate the top mass
      float topMass = getTopMassEstimate(d_llsol, hypIdx, vpfjets_p4, evt_pfmet(), evt_pfmetPhi());
      if(topMass != -999 && 42 != 42) std::cout<<"And top mass from exteral: "<<topMass<<std::endl;

      vector<float> topMassAllComb;
      VofP4 vjpts_p4_Comb;
      int topMassCounter = 0;
      for(int jet1 =  0; jet1 < vjpts_p4.size(); ++jet1) {
      for(int jet2 =  jet1+1; jet2 < vjpts_p4.size(); ++jet2) {
      vjpts_p4_Comb.clear();
      vjpts_p4_Comb.push_back(vjpts_p4.at(jet1));
      vjpts_p4_Comb.push_back(vjpts_p4.at(jet2));
      topMassAllComb.push_back( getTopMassEstimate(d_llsol, hypIdx, vjpts_p4_Comb, tcmet, tcmetphi) );
      //             std::cout<<
      //               "We have "<<vjpts_p4.size()<<
      //               " jets. This is combination: "<<jet1<<
      //               " with "<<jet2<<
      //               "  and topmass "<<topMassAllComb.at(topMassCounter)<<
      //               std::endl;
      ++topMassCounter;
      }
      }
      */

      //other vars for baby

      pass_         = ( npfjets30_ >= 3 && pfmet_ > 25. ) ? 1 : 0;
      proc_         = getProcessType(prefix);       //integer specifying sample
      topmass_      = -999;//topMass;                      //topepton mass //REPLACE TOPMASS
      y_	    = pfmet_ / sqrt( ht_ ); //y=MET/sqrt(HT)

      strcpy(dataset_, cms2.evt_dataset().at(0).Data());  //dataset name
      run_          = evt_run();                          //run
      lumi_         = evt_lumiBlock();                    //lumi
      event_        = evt_event();                        //event
      nvtxweight_   = vtxweight(isData);

      csc_       = cms2.evt_cscTightHaloId();
      hbhe_      = cms2.evt_hbheFilter();
      hcallaser_ = cms2.filt_hcalLaser();
      ecaltp_    = cms2.filt_ecalTP();
      trkfail_   = cms2.filt_trackingFailure();
      eebadsc_   = 1;
      if( isData ) eebadsc_ = cms2.filt_eeBadSc();
      hbhenew_   = passHBHEFilter();

      k_ = 1;
      if     ( strcmp( prefix , "LM0"  )    == 0 ) k_ = kfactorSUSY( "lm0"  );
      else if( strcmp( prefix , "LM1"  )    == 0 ) k_ = kfactorSUSY( "lm1"  );
      else if( strcmp( prefix , "LM2"  )    == 0 ) k_ = kfactorSUSY( "lm2"  );
      else if( strcmp( prefix , "LM3"  )    == 0 ) k_ = kfactorSUSY( "lm3"  );
      else if( strcmp( prefix , "LM4"  )    == 0 ) k_ = kfactorSUSY( "lm4"  );
      else if( strcmp( prefix , "LM5"  )    == 0 ) k_ = kfactorSUSY( "lm5"  );
      else if( strcmp( prefix , "LM6"  )    == 0 ) k_ = kfactorSUSY( "lm6"  );
      else if( strcmp( prefix , "LM7"  )    == 0 ) k_ = kfactorSUSY( "lm7"  );
      else if( strcmp( prefix , "LM8"  )    == 0 ) k_ = kfactorSUSY( "lm8"  );
      else if( strcmp( prefix , "LM9"  )    == 0 ) k_ = kfactorSUSY( "lm9"  );
      else if( strcmp( prefix , "LM10" )    == 0 ) k_ = kfactorSUSY( "lm10" );
      else if( strcmp( prefix , "LM11" )    == 0 ) k_ = kfactorSUSY( "lm11" );
      else if( strcmp( prefix , "LM12" )    == 0 ) k_ = kfactorSUSY( "lm12" );
      else if( strcmp( prefix , "LMscan" )  == 0 ) k_ = kfactorSUSY(m0_,m12_,"tanbeta10");
      

      ecalveto1_ = -1.;
      ecalveto2_ = -1.;
      hcalveto1_ = -1.;
      hcalveto2_ = -1.;
      
      if( leptype_ == 0 ){
	iso1_   = electronIsolation_rel   ( index1 , true ); //truncated
	isont1_ = electronIsolation_rel_v1( index1 , true ); //non-truncated
	etasc1_ = els_etaSC()[index1];
      }
      else if( leptype_ == 1 ){
	iso1_   = muonIsoValue( index1 , true  ); //truncated 
	isont1_ = muonIsoValue( index1 , false ); //non-truncated
	etasc1_ = -999;
	
	ecalveto1_ = mus_iso_ecalvetoDep().at(index1);
	hcalveto1_ = mus_iso_hcalvetoDep().at(index1);
      }
            
      if     ( leptype_ == 0 ) netot += weight_;
      else if( leptype_ == 1 ) nmtot += weight_;

      if( pass_ == 1 ){
	if     ( leptype_ == 0 ) nepass += weight_;
	else if( leptype_ == 1 ) nmpass += weight_;
      }

      //-------------------------------------
      // triggers
      //-------------------------------------
      if (evt_run()<193806 || !isData)
	isomu24_   = passUnprescaledHLTTriggerPattern("HLT_IsoMu24_eta2p1_v")  ? 1 : 0;
      else   
	isomu24_   = passUnprescaledHLTTriggerPattern("HLT_IsoMu24_v")  ? 1 : 0;
      //      isomu24_   = passUnprescaledHLTTriggerPattern("HLT_IsoMu24_v"   )  ? 1 : 0;
      ele27wp80_ = passUnprescaledHLTTriggerPattern("HLT_Ele27_WP80_v")  ? 1 : 0;
      mm_        = passUnprescaledHLTTriggerPattern("HLT_Mu17_Mu8_v")                                     ? 1 : 0;
      mmtk_      = passUnprescaledHLTTriggerPattern("HLT_Mu17_TkMu8_v")                                   ? 1 : 0;
      me_        = passUnprescaledHLTTriggerPattern("HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v") ? 1 : 0;
      em_        = passUnprescaledHLTTriggerPattern("HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v") ? 1 : 0;
      ee_        = passUnprescaledHLTTriggerPattern("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v") ? 1 : 0;
   
      char* isomutrigname = (evt_run()<193806 || !isData) ? 
	(char*) "HLT_IsoMu24_eta2p1_v" : (char*) "HLT_IsoMu24_v";
      // char* isomutrigname = (char*) "HLT_IsoMu24_v";
      // if( !isData ) isomutrigname = (char*) "HLT_IsoMu24_eta2p1_v";

      trgmu1_    = objectPassTrigger( *lep1_ , isomutrigname        , 0.1 ) ? 1 : 0;
      if( ngoodlep_ > 1 )
	trgmu2_    = objectPassTrigger( *lep2_ , isomutrigname      , 0.1 ) ? 1 : 0;

      trgel1_    = objectPassTrigger( *lep1_ , (char*) "HLT_Ele27_WP80_v"     , 0.1 ) ? 1 : 0;
      if( ngoodlep_ > 1 )
	trgel2_    = objectPassTrigger( *lep2_ , (char*) "HLT_Ele27_WP80_v"   , 0.1 ) ? 1 : 0;

      // trgmu30_  = -1;
      // trg2mu30_ = -1;

      // if( isData ){
      // 	if( cms2.evt_run() >= 173212 ) trgmu30_ = objectPassTrigger( *lep1_ , (char*) "HLT_IsoMu30_eta2p1_v" , 0.1 ) ? 1 : 0;
      // 	else                           
	
      // 	if( ngoodlep_ > 1 ){
      // 	  if( cms2.evt_run() >= 173212 ) trg2mu30_ = objectPassTrigger( *lep2_ , (char*) "HLT_IsoMu30_eta2p1_v" , 0.1 ) ? 1 : 0;
      // 	  else                           trg2mu30_ = objectPassTrigger( *lep2_ , (char*) "HLT_IsoMu30_v"        , 0.1 ) ? 1 : 0;
      // 	}

      // }

      // else{

      // 	// for MC, some samples have IsoMu30 and some have IsoMu30_eta2p1, so we need to check

      // 	TString isoMuName    = triggerName("HLT_IsoMu30_v");
      // 	TString isoMu2p1Name = triggerName("HLT_IsoMu30_eta2p1_v");

      // 	if( !isoMuName.Contains("TRIGGER_NOT_FOUND") ){
      // 	  trgmu30_ = objectPassTrigger( *lep1_ , (char*) "HLT_IsoMu30_v" , 0.1 ) ? 1 : 0;
      // 	  if( ngoodlep_ > 1 ) 	  trg2mu30_ = objectPassTrigger( *lep2_ , (char*) "HLT_IsoMu30_v" , 0.1 ) ? 1 : 0;
      // 	}
	
      // 	else if( !isoMu2p1Name.Contains("TRIGGER_NOT_FOUND") ){
      // 	  trgmu30_ = objectPassTrigger( *lep1_ , (char*) "HLT_IsoMu30_eta2p1_v" , 0.1 ) ? 1 : 0;
      // 	  if( ngoodlep_ > 1 ) 	  trg2mu30_ = objectPassTrigger( *lep2_ , (char*) "HLT_IsoMu30_eta2p1_v" , 0.1 ) ? 1 : 0;
      // 	}

      // 	else{
      // 	  trgmu30_ = -2;
      // 	}

      // }

      //set trigger weight
      mutrigweight_  = getMuTriggerWeight   ( lep1_->pt() , lep1_->eta() );
      mutrigweight2_ = getMuTriggerWeightNew( lep1_->pt() , lep1_->eta() );
      
      outTree->Fill();
    
    } // entries

    delete f;
  } // currentFile

  if( nSkip_els_conv_dist > 0 )
    cout << "Skipped " << nSkip_els_conv_dist << " events due to nan in els_conv_dist" << endl;

  cout << endl;
  cout << "Sample: " << prefix << endl;
  cout << endl;
  cout << "-----------------------" << endl;
  cout << "| Lepton yields       |" << endl;
  cout << "-----------------------" << endl;
  cout << "ne  " << netot       << endl;
  cout << "nm  " << nmtot       << endl;
  cout << "tot " << netot+nmtot << endl;

  cout << endl;
  cout << "-----------------------" << endl;
  cout << "| Preselection yields |" << endl;
  cout << "-----------------------" << endl;
  cout << "ne  " << nepass        << endl;
  cout << "nm  " << nmpass        << endl;
  cout << "tot " << nepass+nmpass << endl;
  cout << endl;

  if(g_createTree) closeTree();
  
  already_seen.clear();

  if (nEventsChain != nEventsTotal)
    std::cout << "ERROR: number of events from files is not equal to total number of events" << std::endl;

  stop_xsec_file->Close();

  //delete d_llsol; //REPLACETOPMASS

  return 0;

}


//--------------------------------------------------------------------
 
void singleLeptonLooper::BookHistos(char *prefix)
{
  // Prefix comes from the sample and it is passed to the scanning function
  // Suffix is "ee" "em" "em" "all" which depends on the final state
  // For example: histogram named tt_hnJet_ee would be the Njet distribution
  // for the ee final state in the ttbar sample.
  // MAKE SURE TO CAL SUMW2 FOR EACH 1D HISTOGRAM BEFORE FILLING!!!!!!
  cout << "Begin book histos..." << endl;

  TDirectory *rootdir = gDirectory->GetDirectory("Rint:");
  rootdir->cd();

  h_PU_trkpt = new TH1F(Form("%s_PU_trkpt",prefix),"track pt from PU interactions",100,0,100);
  h_dz_vtx_trk = new TH1F(Form("%s_dz_vtx_trk",prefix),"dZ between vtx and tracks",200,-0.1,0.1);
 
  cout << "End book histos..." << endl;
}// CMS2::BookHistos()

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

void singleLeptonLooper::makeTree(char *prefix, bool doFakeApp, FREnum frmode ){
  TDirectory *rootdir = gDirectory->GetDirectory("Rint:");
  rootdir->cd();


  //char* dir = "";
  //if     ( g_trig == e_lowpt  ) dir = "lowpt";
  //else if( g_trig == e_highpt ) dir = "highpt";

  //Super compressed ntuple here
  char* frsuffix = (char*) "";
  if( doFakeApp ){
    if ( frmode == e_qcd   ) frsuffix = (char*) "_doubleFake";
    if ( frmode == e_wjets ) frsuffix = (char*) "_singleFake";
  }

  char* tpsuffix = (char*) "";
  if( doTenPercent ) tpsuffix = (char*) "_tenPercent";

  outFile   = new TFile(Form("output/%s_smallTree%s%s.root",prefix,frsuffix,tpsuffix), "RECREATE");
  //  outFile   = new TFile(Form("output/%s/%s_smallTree%s%s.root",g_version,prefix,frsuffix,tpsuffix), "RECREATE");
  //outFile   = new TFile("temp.root","RECREATE");
  outFile->cd();
  outTree = new TTree("t","Tree");

  //Set branch addresses
  //variables must be declared in singleLeptonLooper.h
  outTree->Branch("acc_2010",        &acc_2010_,         "acc_2010/I");
  outTree->Branch("acc_highmet",     &acc_highmet_,      "acc_highmet/I");
  outTree->Branch("acc_highht",      &acc_highht_,       "acc_highht/I");

  outTree->Branch("csc"       ,  &csc_       ,  "csc/I");  
  outTree->Branch("hbhe"      ,  &hbhe_      ,  "hbhe/I");  
  outTree->Branch("hbhenew"   ,  &hbhenew_   ,  "hbhenew/I");  
  outTree->Branch("hcallaser" ,  &hcallaser_ ,  "hcallaser/I");  
  outTree->Branch("ecaltp"    ,  &ecaltp_    ,  "ecaltp/I");  
  outTree->Branch("trkfail"   ,  &trkfail_   ,  "trkfail/I");  
  outTree->Branch("eebadsc"   ,  &eebadsc_   ,  "eebadsc/I");  

  outTree->Branch("isdata",          &isdata_,           "isdata/I");
  outTree->Branch("jetid",           &jetid_,            "jetid/I");
  outTree->Branch("jetid30",         &jetid30_,          "jetid30/I");
  outTree->Branch("json",            &json_,             "json/I");
  outTree->Branch("htoffset",        &htoffset_,         "htoffset/F");
  outTree->Branch("htuncor",         &htuncor_,          "htuncor/F");
  outTree->Branch("ptt",             &ptt_,              "ptt/F");
  outTree->Branch("pttbar",          &pttbar_,           "pttbar/F");
  outTree->Branch("ptttbar",         &ptttbar_,          "ptttbar/F");
  outTree->Branch("mttbar",          &mttbar_,           "mttbar/F");
  outTree->Branch("npartons",        &npartons_,         "npartons/I");
  outTree->Branch("hyptype",         &hyptype_,          "hyptype/I");
  outTree->Branch("maxpartonpt",     &maxpartonpt_,      "maxpartonpt/F");
  outTree->Branch("etattbar",        &etattbar_,         "etatbar/F");
  outTree->Branch("njetsoffset",     &njetsoffset_,      "njetsoffset/I");
  outTree->Branch("njetsuncor",      &njetsuncor_,       "njetsuncor/I");
  outTree->Branch("costhetaweight",  &costhetaweight_,   "costhetaweight/F");
  outTree->Branch("weight",          &weight_,           "weight/F");
  outTree->Branch("mutrigweight",    &mutrigweight_,     "mutrigweight/F");
  outTree->Branch("mutrigweight2",   &mutrigweight2_,    "mutrigweight2/F");
  outTree->Branch("trgeff",          &trgeff_,           "trgeff/F");
  outTree->Branch("pthat",           &pthat_,            "pthat/F");
  outTree->Branch("qscale",          &qscale_,           "qscale/F");
  outTree->Branch("mgcor",           &mgcor_,            "mgcor/F");
  outTree->Branch("wflav",           &wflav_,            "wflav/I");
  outTree->Branch("ksusy",           &ksusy_,            "ksusy/F");
  outTree->Branch("ksusyup",         &ksusyup_,          "ksusyup/F");
  outTree->Branch("ksusydn",         &ksusydn_,          "ksusydn/F");
  outTree->Branch("xsecsusy",        &xsecsusy_,         "xsecsusy/F");
  outTree->Branch("xsecsusy2",       &xsecsusy2_,        "xsecsusy2/F");
  outTree->Branch("smeff",           &smeff_,            "smeff/F");
  outTree->Branch("k",               &k_,                "k/F");
  outTree->Branch("mllgen",          &mllgen_,           "mllgen/F");
  outTree->Branch("ptwgen",          &ptwgen_,           "ptwgen/F");
  outTree->Branch("ptzgen",          &ptzgen_,           "ptzgen/F");
  outTree->Branch("nlep",            &nlep_,             "nlep/I");
  outTree->Branch("nosel",           &nosel_,            "nosel/I");
  outTree->Branch("ngoodlep",        &ngoodlep_,         "ngoodlep/I");
  outTree->Branch("ngoodel",         &ngoodel_,          "ngoodel/I");
  outTree->Branch("ngoodmu",         &ngoodmu_,          "ngoodmu/I");
  outTree->Branch("mull",            &mull_,             "mull/I");
  outTree->Branch("mult",            &mult_,             "mult/I");
  //outTree->Branch("eltrijet",        &eltrijet_,         "eltrijet/I");
  //outTree->Branch("mutrijet",        &mutrijet_,         "mutrijet/I");
  //outTree->Branch("ldi",             &ldi_,              "ldi/I");
  //outTree->Branch("ltri",            &ltri_,             "ltri/I");
  //outTree->Branch("smu",             &smu_,              "smu/I");
  //outTree->Branch("smu30",           &smu30_,            "smu30/I");
  //outTree->Branch("trgmu30",         &trgmu30_,          "trgmu30/I");
  //outTree->Branch("trg2mu30",        &trg2mu30_,         "trg2mu30/I");
  //outTree->Branch("dil",             &dil_,              "dil/I");
  outTree->Branch("mullgen",         &mullgen_,          "mullgen/I");
  outTree->Branch("multgen",         &multgen_,          "multgen/I");
  outTree->Branch("proc",            &proc_,             "proc/I");
  outTree->Branch("leptype",         &leptype_,          "leptype/I");
  outTree->Branch("topmass",         &topmass_,          "topmass/F");
  outTree->Branch("dilmass",         &dilmass_,          "dilmass/F");
  outTree->Branch("dilrecoil",       &dilrecoil_,        "dilrecoil/F");
  outTree->Branch("dilrecoilparl",   &dilrecoilparl_,    "dilrecoilparl/F");
  outTree->Branch("dilrecoilperp",   &dilrecoilperp_,    "dilrecoilperp/F");
  outTree->Branch("tcmet",           &tcmet_,            "tcmet/F");
  outTree->Branch("genmet",          &genmet_,           "genmet/F");
  outTree->Branch("gensumet",        &gensumet_,         "gensumet/F");
  outTree->Branch("genmetphi",       &genmetphi_,        "genmetphi/F");
  outTree->Branch("pfmet",           &pfmet_,            "pfmet/F");
  outTree->Branch("pfmetveto",       &pfmetveto_,        "pfmetveto/F");
  outTree->Branch("pfmetsig",        &pfmetsig_,         "pfmetsig/F");
  outTree->Branch("pfmetphi",        &pfmetphi_,         "pfmetphi/F");
  outTree->Branch("pfsumet",         &pfsumet_,          "pfsumet/F");
  outTree->Branch("mucormet",        &mucormet_,         "mucormet/F");
  outTree->Branch("mucorjesmet",     &mucorjesmet_,      "mucorjesmet/F");
  outTree->Branch("tcmet35X",        &tcmet_35X_,        "tcmet35X/F");
  outTree->Branch("tcmetevent",      &tcmet_event_,      "tcmetevent/F");
  outTree->Branch("tcmetlooper",     &tcmet_looper_,     "tcmetlooper/F");
  outTree->Branch("tcmetphi",        &tcmetphi_,         "tcmetphi/F");
  outTree->Branch("tcsumet",         &tcsumet_,          "tcsumet/F");
  outTree->Branch("tcmetUp",         &tcmetUp_,          "tcmetUp/F");
  outTree->Branch("tcmetDown",       &tcmetDown_,        "tcmetDown/F");
  outTree->Branch("tcmetTest",       &tcmetTest_,        "tcmetTest/F");
  outTree->Branch("pfmetUp",         &pfmetUp_,          "pfmetUp/F");
  outTree->Branch("pfmetDown",       &pfmetDown_,        "pfmetDown/F");
  outTree->Branch("pfmetTest",       &pfmetTest_,        "pfmetTest/F");
  outTree->Branch("sumjetpt",        &sumjetpt_,         "sumjetpt/F");
  outTree->Branch("dileta",          &dileta_,           "dileta/F");
  outTree->Branch("dilpt",           &dilpt_,            "dilpt/F");
  outTree->Branch("dildphi",         &dildphi_,          "dildphi/F");
  outTree->Branch("ngenjets",        &ngenjets_,         "ngenjets/I");
  outTree->Branch("njpt",            &njpt_,             "njpt/I");

  outTree->Branch("trgmu1"         ,  &trgmu1_          ,    "trgmu1/I"      );
  outTree->Branch("trgmu2"         ,  &trgmu2_          ,    "trgmu2/I"      );
  outTree->Branch("trgel1"         ,  &trgel1_          ,    "trgel1/I"      );
  outTree->Branch("trgel2"         ,  &trgel2_          ,    "trgel2/I"      );

  outTree->Branch("isomu24"        ,  &isomu24_         ,    "isomu24/I"     );
  outTree->Branch("ele27wp80"      ,  &ele27wp80_       ,    "ele27wp80/I"   );
  outTree->Branch("mm"             ,  &mm_              ,    "mm/I"          );
  outTree->Branch("mmtk"           ,  &mmtk_            ,    "mmtk/I"        );
  outTree->Branch("me"             ,  &me_              ,    "me/I"          );
  outTree->Branch("em"             ,  &em_              ,    "em/I"          );
  outTree->Branch("mu"             ,  &mu_              ,    "mu/I"          );
  outTree->Branch("ee"             ,  &ee_              ,    "ee/I"          );

  // pfjets L1FastL2L3Res
  outTree->Branch("npfjets30",        &npfjets30_,        "npfjets30/I");
  outTree->Branch("npfjets35",        &npfjets35_,        "npfjets35/I");
  outTree->Branch("npfjets40",        &npfjets40_,        "npfjets40/I");
  outTree->Branch("npfjets45",        &npfjets45_,        "npfjets45/I");
  outTree->Branch("npfjets30lepcorr", &npfjets30lepcorr_, "npfjets30lepcorr/I");
  outTree->Branch("knjets",           &knjets_,           "knjets/F");

  //rho correction
  outTree->Branch("rhovor",          &rhovor_,           "rhovor/F");

  outTree->Branch("htpf30",          &htpf30_,           "htpf30/F");
  outTree->Branch("htpf35",          &htpf35_,           "htpf35/F");
  outTree->Branch("htpf40",          &htpf40_,           "htpf40/F");
  outTree->Branch("htpf45",          &htpf45_,           "htpf45/F");

  // type1 met flavors
  outTree->Branch("t1met10",         &t1met10_,          "t1met10/F");
  outTree->Branch("t1met20",         &t1met20_,          "t1met20/F");
  outTree->Branch("t1met30",         &t1met30_,          "t1met30/F");
  outTree->Branch("t1met10phi",      &t1met10phi_,       "t1met10phi/F");
  outTree->Branch("t1met20phi",      &t1met20phi_,       "t1met20phi/F");
  outTree->Branch("t1met30phi",      &t1met30phi_,       "t1met30phi/F");
  outTree->Branch("t1met10mt",       &t1met10mt_,        "t1met10mt/F");
  outTree->Branch("t1met20mt",       &t1met20mt_,        "t1met20mt/F");
  outTree->Branch("t1met30mt",       &t1met30mt_,        "t1met30mt/F");
  outTree->Branch("lepmetpt",        &lepmetpt_,         "lepmetpt/F");
  outTree->Branch("lept1met10pt",    &lept1met10pt_,     "lept1met10pt/F");

  //met variables with phi correction
  outTree->Branch("t1metphicorr"       , &t1metphicorr_       , "t1metphicorr/F");
  outTree->Branch("t1metphicorrup"     , &t1metphicorrup_     , "t1metphicorrup/F");
  outTree->Branch("t1metphicorrdn"     , &t1metphicorrdn_     , "t1metphicorrdn/F");
  outTree->Branch("t1metphicorrphi"    , &t1metphicorrphi_    , "t1metphicorrphi/F");
  outTree->Branch("t1metphicorrphiup"  , &t1metphicorrphiup_  , "t1metphicorrphiup/F");
  outTree->Branch("t1metphicorrphidn"  , &t1metphicorrphidn_  , "t1metphicorrphidn/F");
  outTree->Branch("t1metphicorrlep"    , &t1metphicorrlep_    , "t1metphicorrlep/F");
  outTree->Branch("t1metphicorrlepphi" , &t1metphicorrlepphi_ , "t1metphicorrlepphi/F");
  outTree->Branch("t1metphicorrmt"     , &t1metphicorrmt_     , "t1metphicorrmt/F");
  outTree->Branch("t1metphicorrmtup"   , &t1metphicorrmtup_   , "t1metphicorrmtup/F");
  outTree->Branch("t1metphicorrmtdn"   , &t1metphicorrmtdn_   , "t1metphicorrmtdn/F");
  outTree->Branch("t1metphicorrlepmt"  , &t1metphicorrlepmt_  , "t1metphicorrlepmt/F");

  outTree->Branch("htpfres30",        &htpfres30_,        "htpfres30/F");
  outTree->Branch("htpfres35",        &htpfres35_,        "htpfres35/F");
  outTree->Branch("htpfres40",        &htpfres40_,        "htpfres40/F");
  outTree->Branch("htpfres45",        &htpfres45_,        "htpfres45/F");
				      
  // btag variables		      
  outTree->Branch("nbtagsssv",        &nbtagsssv_,        "nbtagsssv/I");
  outTree->Branch("nbtagstcl",        &nbtagstcl_,        "nbtagstcl/I");
  outTree->Branch("nbtagstcm",        &nbtagstcm_,        "nbtagstcm/I");
  outTree->Branch("nbtagscsvl",       &nbtagscsvl_,       "nbtagscsvl/I");
  outTree->Branch("nbtagscsvm",       &nbtagscsvm_,       "nbtagscsvm/I");
  outTree->Branch("nbtagscsvt",       &nbtagscsvt_,       "nbtagscsvt/I");
  outTree->Branch("nbtagsssvcorr",    &nbtagsssvcorr_,    "nbtagsssvcorr/I");
  outTree->Branch("nbtagstclcorr",    &nbtagstclcorr_,    "nbtagstclcorr/I");
  outTree->Branch("nbtagstcmcorr",    &nbtagstcmcorr_,    "nbtagstcmcorr/I");
  outTree->Branch("nbtagscsvlcorr",   &nbtagscsvlcorr_,   "nbtagscsvlcorr/I");
  outTree->Branch("nbtagscsvmcorr",   &nbtagscsvmcorr_,   "nbtagscsvmcorr/I");
  outTree->Branch("nbtagscsvtcott",   &nbtagscsvtcorr_,   "nbtagscsvtcorr/I");
  outTree->Branch("bjet1",            &bjet1_,            "bjet1/I");
  outTree->Branch("bjet2",            &bjet2_,            "bjet2/I");
  outTree->Branch("bjet3",            &bjet3_,            "bjet3/I");
  outTree->Branch("bjet4",            &bjet4_,            "bjet4/I");
  outTree->Branch("bjet5",            &bjet5_,            "bjet5/I");
  outTree->Branch("bjet6",            &bjet6_,            "bjet6/I");
  outTree->Branch("lepjet1",          &lepjet1_,          "lepjet1/I");
  outTree->Branch("lepjet2",          &lepjet2_,          "lepjet2/I");
  outTree->Branch("lepjet3",          &lepjet3_,          "lepjet3/I");
  outTree->Branch("lepjet4",          &lepjet4_,          "lepjet4/I");
  outTree->Branch("lepjet5",          &lepjet5_,          "lepjet5/I");
  outTree->Branch("lepjet6",          &lepjet6_,          "lepjet6/I");
  outTree->Branch("qgjet1",           &qgjet1_,           "qgjet1/I");
  outTree->Branch("qgjet2",           &qgjet2_,           "qgjet2/I");
  outTree->Branch("qgjet3",           &qgjet3_,           "qgjet3/I");
  outTree->Branch("qgjet4",           &qgjet4_,           "qgjet4/I");
  outTree->Branch("qgjet5",           &qgjet5_,           "qgjet5/I");
  outTree->Branch("qgjet6",           &qgjet6_,           "qgjet6/I");
  outTree->Branch("njetsUp",          &njetsUp_,          "njetsUp/I");
  outTree->Branch("njetsDown",        &njetsDown_,        "njetsDown/I");
  outTree->Branch("htUp",             &htUp_,             "htUp/F");
  outTree->Branch("htDown",           &htDown_,           "htDown/F");
  outTree->Branch("npu",              &npu_,              "npu/I");
  outTree->Branch("npuMinusOne",      &npuMinusOne_,      "npuMinusOne/I");
  outTree->Branch("npuPlusOne",       &npuPlusOne_,       "npuPlusOne/I");
  outTree->Branch("nvtx",             &nvtx_,             "nvtx/I");
  outTree->Branch("nvtxweight",       &nvtxweight_,       "nvtxweight/F");
  outTree->Branch("n3dvtxweight",     &n3dvtxweight_,     "n3dvtxweight/F");
  outTree->Branch("pdfid1",           &pdfid1_,           "pdfid1/I");
  outTree->Branch("pdfid2",           &pdfid2_,           "pdfid2/I");
  outTree->Branch("pdfx1",            &pdfx1_,            "pdfx1/F");
  outTree->Branch("pdfx2",            &pdfx2_,            "pdfx2/F");
  outTree->Branch("pdfQ",             &pdfQ_,             "pdfQ/F");
  outTree->Branch("vecjetpt",         &vecjetpt_,         "vecjetpt/F");
  outTree->Branch("pass",             &pass_,             "pass/I");
  outTree->Branch("passz",            &passz_,            "passz/I");
  outTree->Branch("m0",               &m0_,               "m0/F");
  outTree->Branch("mg",               &mG_,               "mg/F");
  outTree->Branch("ml",               &mL_,               "ml/F");
  outTree->Branch("x",                &x_,                "x/F");
  outTree->Branch("m12",              &m12_,              "m12/F");
  outTree->Branch("lep1chi2ndf",      &lep1chi2ndf_,      "lep1chi2ndf/F");
  outTree->Branch("lep2chi2ndf",      &lep2chi2ndf_,      "lep2chi2ndf/F");
  outTree->Branch("lep1dpt",          &lep1dpt_,          "lep1dpt/F");
  outTree->Branch("lep2dpt",          &lep2dpt_,          "lep2dpt/F");
  outTree->Branch("id1",              &id1_,              "id1/I");
  outTree->Branch("id2",              &id2_,              "id2/I");
  outTree->Branch("leptype1",         &leptype1_,         "leptype1/I");
  outTree->Branch("leptype2",         &leptype2_,         "leptype2/I");
  outTree->Branch("w1",               &w1_,               "w1/I");
  outTree->Branch("w2",               &w2_,               "w2/I");
  outTree->Branch("iso1",             &iso1_,             "iso1/F");
  outTree->Branch("isont1",           &isont1_,           "isont1/F");
  outTree->Branch("etasc1",           &etasc1_,           "etasc1/F");
  outTree->Branch("etasc2",           &etasc2_,           "etasc2/F");
  outTree->Branch("iso2",             &iso2_,             "iso2/F");
  outTree->Branch("ecalveto1",        &ecalveto1_,        "ecalveto1/F");
  outTree->Branch("ecalveto2",        &ecalveto2_,        "ecalveto2/F");
  outTree->Branch("hcalveto1",        &hcalveto1_,        "hcalveto1/F");
  outTree->Branch("hcalveto2",        &hcalveto2_,        "hcalveto2/F");
  outTree->Branch("isont2",           &isont2_,           "isont2/F");
  outTree->Branch("ptl1",             &ptl1_,             "ptl1/F");
  outTree->Branch("ptl2",             &ptl2_,             "ptl2/F");
  outTree->Branch("etal1",            &etal1_,            "etal1/F");
  outTree->Branch("etal2",            &etal2_,            "etal2/F");
  outTree->Branch("phil1",            &phil1_,            "phil1/F");
  outTree->Branch("phil2",            &phil2_,            "phil2/F");
  outTree->Branch("meff",             &meff_,             "meff/F");
  outTree->Branch("mt",               &mt_,               "mt/F");
  outTree->Branch("dataset",          &dataset_,          "dataset[200]/C");
  outTree->Branch("run",              &run_,              "run/I");
  outTree->Branch("lumi",             &lumi_,             "lumi/I");
  outTree->Branch("event",            &event_,            "event/I");
  outTree->Branch("y",                &y_,                "y/F");  
  outTree->Branch("ht",               &ht_,               "ht/F");  
  outTree->Branch("htgen",            &htgen_,            "htgen/F");  
  outTree->Branch("htjpt",            &htjpt_,            "htjpt/F");  
  outTree->Branch("nels",             &nels_,             "nels/I");  
  outTree->Branch("nmus",             &nmus_,             "nmus/I");  
  outTree->Branch("ntaus",            &ntaus_,            "ntaus/I");  
  outTree->Branch("nleps",            &nleps_,            "nleps/I");  
  outTree->Branch("dphijm",           &dphijm_,           "dphijm/F");  
  outTree->Branch("ptjetraw",         &ptjetraw_,         "ptjetraw/F");  
  outTree->Branch("ptjet23",          &ptjet23_,          "ptjet23/F");  
  outTree->Branch("ptjetF23",         &ptjetF23_,         "ptjetF23/F");  
  outTree->Branch("ptjetO23",         &ptjetO23_,         "ptjetO23/F");  
  //outTree->Branch("cosphijz",         &cosphijz_,         "cosphijz/F");  
  outTree->Branch("mcid1",            &mcid1_,            "mcid1/I");  
  outTree->Branch("mcdr1",            &mcdr1_,            "mcdr1/F");  
  outTree->Branch("mcdecay1",         &mcdecay1_,         "mcdecay1/I");  
  outTree->Branch("mcndec1",          &mcndec1_,          "mcndec1/I");  
  outTree->Branch("mcndec2",          &mcndec2_,          "mcndec2/I");  
  outTree->Branch("mcndeckls1",       &mcndeckls1_,       "mcndeckls1/I");  
  outTree->Branch("mcndeckls2",       &mcndeckls2_,       "mcndeckls2/I");  
  outTree->Branch("mcndecem1",        &mcndecem1_,        "mcndecem1/I");  
  outTree->Branch("mcndecem2",        &mcndecem2_,        "mcndecem2/I");  
  outTree->Branch("mcid2",            &mcid2_,            "mcid2/I");  
  outTree->Branch("mcdr2",            &mcdr2_,            "mcdr2/F");  
  outTree->Branch("mcdecay2",         &mcdecay2_,         "mcdecay2/I");  
  outTree->Branch("mctaudpt1",        &mctaudpt1_,        "mctaudpt1/F");  
  outTree->Branch("mctaudpt2",        &mctaudpt2_,        "mctaudpt2/F");  
  outTree->Branch("mctaudid1",        &mctaudid1_,        "mctaudid1/I");  
  outTree->Branch("mctaudid2",        &mctaudid2_,        "mctaudid2/I");  
  outTree->Branch("mlepid",           &mlepid_,           "mlepid/I");  
  outTree->Branch("mleppassid",       &mleppassid_,       "mleppassid/I");  
  outTree->Branch("mleppassiso",      &mleppassiso_,      "mleppassiso/I");  
  outTree->Branch("mlepiso",          &mlepiso_,          "mlepiso/F");  
  outTree->Branch("mlepdr",           &mlepdr_,           "mlepdr/F");  
  outTree->Branch("pflepiso",         &pflepiso_,         "pflepiso/F");  
  outTree->Branch("pflepdr",          &pflepdr_,          "pflepdr/F");  
  outTree->Branch("pfleppt",          &pfleppt_,          "pfleppt/F");  
  outTree->Branch("pflepmindrj",      &pflepmindrj_,      "pflepmindrj/F");  
  outTree->Branch("pftaudiso",        &pftaudiso_,        "pftaudiso/F");  
  outTree->Branch("pftauddr",         &pftauddr_,         "pftauddr/F");  
  outTree->Branch("pftaudpt",         &pftaudpt_,         "pftaudpt/F");  
  outTree->Branch("pftaudmindrj",     &pftaudmindrj_,     "pftaudmindrj/F");  
  outTree->Branch("pfcandiso5",       &pfcandiso5_,       "pfcandiso5/F");  
  outTree->Branch("pfcandpt5",        &pfcandpt5_,        "pfcandpt5/F");  
  outTree->Branch("pfcandmindrj5",    &pfcandmindrj5_,    "pfcandmindrj5/F");  
  outTree->Branch("pfcandiso10",      &pfcandiso10_,      "pfcandiso10/F");  
  outTree->Branch("pfcandpt10",       &pfcandpt10_,       "pfcandpt10/F");  
  outTree->Branch("pfcandmindrj10",   &pfcandmindrj10_,   "pfcandmindrj10/F");  
  outTree->Branch("emjet10",          &emjet10_,          "emjet10/F");  
  outTree->Branch("mjj",              &mjj_,              "mjj/F");  
  outTree->Branch("emjet20",          &emjet20_,          "emjet20/F");  
  outTree->Branch("trkpt5",           &trkpt5_,           "trkpt5/F");  
  outTree->Branch("trkpt10",          &trkpt10_,          "trkpt10/F");  
  outTree->Branch("mleptrk5",         &mleptrk5_,         "mleptrk5/F");  
  outTree->Branch("mleptrk10",        &mleptrk10_,        "mleptrk10/F");  
  outTree->Branch("trkreliso5",       &trkreliso5_,       "trkreliso5/F");  
  outTree->Branch("trkreliso10",      &trkreliso10_,      "trkreliso10/F");  
  outTree->Branch("trkpt5loose",      &trkpt5loose_,      "trkpt5loose/F");  
  outTree->Branch("trkpt10loose",     &trkpt10loose_,     "trkpt10loose/F");  
  outTree->Branch("trkreliso5loose",  &trkreliso5loose_,  "trkreliso5loose/F");  
  outTree->Branch("trkreliso10loose", &trkreliso10loose_, "trkreliso10loose/F");  

  outTree->Branch("trkpt10pt0p1",     &trkpt10pt0p1_,      "trkpt10pt0p1/F");  
  outTree->Branch("trkpt10pt0p2",     &trkpt10pt0p2_,      "trkpt10pt0p2/F");  
  outTree->Branch("trkpt10pt0p3",     &trkpt10pt0p3_,      "trkpt10pt0p3/F");  
  outTree->Branch("trkpt10pt0p4",     &trkpt10pt0p4_,      "trkpt10pt0p4/F");  
  outTree->Branch("trkpt10pt0p5",     &trkpt10pt0p5_,      "trkpt10pt0p5/F");  
  outTree->Branch("trkpt10pt0p6",     &trkpt10pt0p6_,      "trkpt10pt0p6/F");  
  outTree->Branch("trkpt10pt0p7",     &trkpt10pt0p7_,      "trkpt10pt0p7/F");  
  outTree->Branch("trkpt10pt0p8",     &trkpt10pt0p8_,      "trkpt10pt0p8/F");  
  outTree->Branch("trkpt10pt0p9",     &trkpt10pt0p9_,      "trkpt10pt0p9/F");  
  outTree->Branch("trkpt10pt1p0",     &trkpt10pt1p0_,      "trkpt10pt1p0/F");  
  outTree->Branch("trkreliso10pt0p1", &trkreliso10pt0p1_,  "trkreliso10pt0p1/F");  
  outTree->Branch("trkreliso10pt0p2", &trkreliso10pt0p2_,  "trkreliso10pt0p2/F");  
  outTree->Branch("trkreliso10pt0p3", &trkreliso10pt0p3_,  "trkreliso10pt0p3/F");  
  outTree->Branch("trkreliso10pt0p4", &trkreliso10pt0p4_,  "trkreliso10pt0p4/F");  
  outTree->Branch("trkreliso10pt0p5", &trkreliso10pt0p5_,  "trkreliso10pt0p5/F");  
  outTree->Branch("trkreliso10pt0p6", &trkreliso10pt0p6_,  "trkreliso10pt0p6/F");  
  outTree->Branch("trkreliso10pt0p7", &trkreliso10pt0p7_,  "trkreliso10pt0p7/F");  
  outTree->Branch("trkreliso10pt0p8", &trkreliso10pt0p8_,  "trkreliso10pt0p8/F");  
  outTree->Branch("trkreliso10pt0p9", &trkreliso10pt0p9_,  "trkreliso10pt0p9/F");  
  outTree->Branch("trkreliso10pt1p0", &trkreliso10pt1p0_,  "trkreliso10pt1p0/F");  

  outTree->Branch("pfcandpt10pt0p1",  &pfcandpt10pt0p1_,   "pfcandpt10pt0p1/F");  
  outTree->Branch("pfcandpt10pt0p2",  &pfcandpt10pt0p2_,   "pfcandpt10pt0p2/F");  
  outTree->Branch("pfcandpt10pt0p3",  &pfcandpt10pt0p3_,   "pfcandpt10pt0p3/F");  
  outTree->Branch("pfcandpt10pt0p4",  &pfcandpt10pt0p4_,   "pfcandpt10pt0p4/F");  
  outTree->Branch("pfcandpt10pt0p5",  &pfcandpt10pt0p5_,   "pfcandpt10pt0p5/F");  
  outTree->Branch("pfcandpt10pt0p6",  &pfcandpt10pt0p6_,   "pfcandpt10pt0p6/F");  
  outTree->Branch("pfcandpt10pt0p7",  &pfcandpt10pt0p7_,   "pfcandpt10pt0p7/F");  
  outTree->Branch("pfcandpt10pt0p8",  &pfcandpt10pt0p8_,   "pfcandpt10pt0p8/F");  
  outTree->Branch("pfcandpt10pt0p9",  &pfcandpt10pt0p9_,   "pfcandpt10pt0p9/F");  
  outTree->Branch("pfcandpt10pt1p0",  &pfcandpt10pt1p0_,   "pfcandpt10pt1p0/F");  
  outTree->Branch("pfcandiso10pt0p1", &pfcandiso10pt0p1_,  "pfcandiso10pt0p1/F");  
  outTree->Branch("pfcandiso10pt0p2", &pfcandiso10pt0p2_,  "pfcandiso10pt0p2/F");  
  outTree->Branch("pfcandiso10pt0p3", &pfcandiso10pt0p3_,  "pfcandiso10pt0p3/F");  
  outTree->Branch("pfcandiso10pt0p4", &pfcandiso10pt0p4_,  "pfcandiso10pt0p4/F");  
  outTree->Branch("pfcandiso10pt0p5", &pfcandiso10pt0p5_,  "pfcandiso10pt0p5/F");  
  outTree->Branch("pfcandiso10pt0p6", &pfcandiso10pt0p6_,  "pfcandiso10pt0p6/F");  
  outTree->Branch("pfcandiso10pt0p7", &pfcandiso10pt0p7_,  "pfcandiso10pt0p7/F");  
  outTree->Branch("pfcandiso10pt0p8", &pfcandiso10pt0p8_,  "pfcandiso10pt0p8/F");  
  outTree->Branch("pfcandiso10pt0p9", &pfcandiso10pt0p9_,  "pfcandiso10pt0p9/F");  
  outTree->Branch("pfcandiso10pt1p0", &pfcandiso10pt1p0_,  "pfcandiso10pt1p0/F");  

  outTree->Branch("mbb",             &mbb_,              "mbb/F");
  outTree->Branch("lep1pfjetdr",     &lep1pfjetdr_,      "lep1pfjetdr/F");  
  outTree->Branch("lep2pfjetdr",     &lep2pfjetdr_,      "lep2pfjetdr/F");  

  outTree->Branch("mclep"    , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &mclep_      );
  outTree->Branch("mcnu"     , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &mcnu_       );
  outTree->Branch("mcmln",           &mcmln_,              "mcmln/F");
  outTree->Branch("mcmtln",          &mcmtln_,             "mcmtln/F");

  outTree->Branch("mlep"      , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &mlep_	);
  outTree->Branch("lep1"      , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &lep1_	);
  outTree->Branch("lep2"      , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &lep2_	);
  outTree->Branch("trklep1"   , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &trklep1_	);
  outTree->Branch("trklep2"   , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &trklep2_	);
  outTree->Branch("gfitlep1"  , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &gfitlep1_	);
  outTree->Branch("gfitlep2"  , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &gfitlep2_	);
  outTree->Branch("lepp"      , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &lepp_	);
  outTree->Branch("lepm"      , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &lepm_	);
  outTree->Branch("pflep1"    , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &pflep1_	);
  outTree->Branch("pflep2"    , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &pflep2_	);
  outTree->Branch("leppfjet1" , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &leppfjet1_	);
  outTree->Branch("leppfjet2" , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &leppfjet2_	);
  outTree->Branch("mclep1"    , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &mclep1_	);
  outTree->Branch("mclep2"    , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &mclep2_	);
  outTree->Branch("mctaud1"   , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &mctaud1_	);
  outTree->Branch("mctaud2"   , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &mctaud2_	);
  outTree->Branch("mctaudvis1"   , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &mctaudvis1_	);
  outTree->Branch("mctaudvis2"   , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &mctaudvis2_	);
  outTree->Branch("pflep"     , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &pflep_	);
  outTree->Branch("pftaud"    , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &pftaud_	);
  outTree->Branch("pfcand5"   , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &pfcand5_	);
  outTree->Branch("pfcand10"  , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &pfcand10_	);
  outTree->Branch("jet"	      , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &jet_	);

  outTree->Branch("pfjet1"    , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &pfjet1_	);
  outTree->Branch("pfjet2"    , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &pfjet2_	);
  outTree->Branch("pfjet3"    , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &pfjet3_	);
  outTree->Branch("pfjet4"    , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &pfjet4_	);
  outTree->Branch("pfjet5"    , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &pfjet5_	);
  outTree->Branch("pfjet6"    , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &pfjet6_	);

  outTree->Branch("nonisoel"  , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &nonisoel_	);
  outTree->Branch("nonisomu"  , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &nonisomu_	);
  outTree->Branch("t"         , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &t_   	);
  outTree->Branch("tbar"      , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &tbar_   	);
  outTree->Branch("ttbar"     , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &ttbar_   	);


}

//--------------------------------------------------------------------

float singleLeptonLooper::dz_trk_vtx( const unsigned int trkidx, const unsigned int vtxidx ){
  
  return ((cms2.trks_vertex_p4()[trkidx].z()-cms2.vtxs_position()[vtxidx].z()) - ((cms2.trks_vertex_p4()[trkidx].x()-cms2.vtxs_position()[vtxidx].x()) * cms2.trks_trk_p4()[trkidx].px() + (cms2.trks_vertex_p4()[trkidx].y() - cms2.vtxs_position()[vtxidx].y()) * cms2.trks_trk_p4()[trkidx].py())/cms2.trks_trk_p4()[trkidx].pt() * cms2.trks_trk_p4()[trkidx].pz()/cms2.trks_trk_p4()[trkidx].pt());
  
}

float singleLeptonLooper::trackIso( int thisPf , float coneR , float dz_thresh , bool dovtxcut , float pt_thresh ){

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
	fillOverFlow( h_dz_vtx_trk , mydz );
	
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

// std::vector<float> trackIsoPtRanges( int thisPf , float coneR , float dz_thresh ){

//   float iso[10];
//   for (int i=0; i<10; ++i) iso[i] = 0.0;
  
//   for (int ipf = 0; ipf < (int)cms2.pfcands_p4().size(); ipf++) {

//     if( ipf == thisPf                 ) continue; // skip this PFCandidate
//     if( cms2.pfcands_charge().at(ipf) == 0 ) continue; // skip neutrals

//     if( dRbetweenVectors( pfcands_p4().at(ipf) , pfcands_p4().at(thisPf) ) > coneR ) continue;

//     int itrk = cms2.pfcands_trkidx().at(ipf);
    
//     if( itrk >= (int)trks_trk_p4().size() || itrk < 0 ){
//       //note: this should only happen for electrons which do not have a matched track
//       //currently we are just ignoring these guys
//       continue;
//     }
    
//     float mindz = dz_trk_vtx(itrk,0);
//     if ( fabs(mindz) > dz_thresh )     continue;

//     //---------------------------------------
//     // passes cuts, add up isolation value
//     //---------------------------------------
    
//     //figure out which pT threshold this track passes
//     float pfcandpt = cms2.pfcands_p4().at(ipf).pt();
//     for int (i=0; i<10; ++i) {
// 	if ( pfcandpt > 0.1*i+0.1 ) 
// 	  iso[i] += pfcandpt;	  
//       }
//   } 

//   std::vector<float> isovec;
//   for (int i=0; i<10; ++i) isovec.push_back(iso[i]);
  
//   return isovec;

// }

std::vector<float> singleLeptonLooper::totalIso( int thisPf , float coneR , float dz_thresh ){

  float iso = 0.; 

  float emiso = 0.;
  float nhiso = 0.;
  float chiso = 0.;

  for (int ipf = 0; ipf < (int)cms2.pfcands_p4().size(); ipf++) {

    if( ipf == thisPf                 ) continue; // skip this PFCandidate
    float dR = dRbetweenVectors( pfcands_p4().at(ipf) , pfcands_p4().at(thisPf) );
    if( dR > coneR ) continue;
    float pfpt = cms2.pfcands_p4().at(ipf).pt();

    //----------------------------------------
    // neutrals
    //----------------------------------------

    if( cms2.pfcands_charge().at(ipf) == 0 ) {
      // skip neutrals with pT < 1 GeV to reduce pileup dependence
      if ( pfpt < 1. ) continue; 
      int pfid = abs(cms2.pfcands_particleId().at(ipf));
        // get isolation parameters
        
        float pfeta = cms2.pfcands_p4().at(ipf).eta();
        float deta = fabs(pfeta - cms2.pfcands_p4().at(thisPf).eta());
	//photons 
	if (pfid == 22) {
	  // to remove pi0s and possible radiation from their photons
	  if (deta <= 0.1) continue;
	  iso += pfpt;
	  emiso += pfpt;
	}
	else {
	  iso += pfpt;
	  nhiso += pfpt;
	}
    }

    //----------------------------------------
    // charged candidates
    //----------------------------------------

    else {
      
      int itrk = cms2.pfcands_trkidx().at(ipf);
    
      if( itrk >= (int)trks_trk_p4().size() || itrk < 0 ){
	//note: this should only happen for electrons which do not have a matched track
	//currently we are just ignoring these guys
	continue; 
      }
    
      //----------------------------------------
      // find closest PV and dz w.r.t. that PV
      //----------------------------------------
    
      // float mindz = 999.;
      // int vtxi    = -1;
      
      // for (unsigned int ivtx = 0; ivtx < cms2.davtxs_position().size(); ivtx++) {
	
      // 	if(!isGoodDAVertex(ivtx)) continue;

      // 	float mydz = dz_trk_vtx(itrk,ivtx);
      
      // 	if (fabs(mydz) < fabs(mindz)) {
      // 	  mindz = mydz;
      // 	  vtxi = ivtx;
      // 	}
         
      // }
    
      //----------------------------------------------------------------------------
      // require closest PV is signal PV, dz cut, exclude tracks near hyp leptons
      //----------------------------------------------------------------------------
    
      // if ( vtxi != 0               )     continue;
      float mindz = dz_trk_vtx(itrk,0);
      if ( fabs(mindz) > dz_thresh )     continue;

      //---------------------------------------
      // passes cuts, add up isolation value
      //---------------------------------------

      iso += pfpt;
      chiso += pfpt;
    }

  }//end loop over pfcands

  // cout<<"total iso: "<<iso<<" sum of emiso + nhiso + chiso: "
  //     <<emiso<<" + "<<nhiso<<" + "<<chiso<<" = "<<emiso+nhiso+chiso
  //     <<endl;
  //  return iso;
  std::vector<float> isos;
  isos.push_back(chiso);
  isos.push_back(emiso);
  isos.push_back(nhiso);
  return isos;
}

/*
pair<float,float> singleLeptonLooper::getPhiCorrMET( float met, float metphi, float sumet, bool ismc, bool isA)
{
  float metx = met * cos( metphi );
  float mety = met * sin( metphi );
  
  //Values for 2011A
  float cx0A = ismc ? -0.09389 : -0.3365;
  float cx1A = ismc ? 0.0001815 : 0.004801;
  float cy0A = ismc ? 0.1571 :  0.2578; 
  float cy1A = ismc ? -0.00371 : -0.006124;
  //Values for 2011B
  float cx0B = ismc ? -0.1070 : -0.3265;
  float cx1B = ismc ? 0.00009587 : 0.005162;
  float cy0B = ismc ? 0.01517 : -0.1956; 
  float cy1B = ismc ? -0.003357 : -0.006299;

  float cx0 = 0.;
  float cx1 = 0.;
  float cy0 = 0.;
  float cy1 = 0.;
  if (ismc) {
    //For MC Use lumi weighted av
    float lumiA = 2.1;
    float lumiB = 2.6;
    cx0 = (cx0A*lumiA + cx0B*lumiB)/(lumiA+lumiB);
    cx1 = (cx1A*lumiA + cx1B*lumiB)/(lumiA+lumiB);
    cy0 = (cy0A*lumiA + cy0B*lumiB)/(lumiA+lumiB);
    cy1 = (cy1A*lumiA + cy1B*lumiB)/(lumiA+lumiB);
  } else {
    cx0 = isA ? cx0A : cx0B;    
    cx1 = isA ? cx1A : cx1B;
    cy0 = isA ? cy0A : cy0B;
    cy1 = isA ? cy1A : cy1B;
  }

  metx -= cx0 + cx1 * sumet;
  mety -= cy0 + cy1 * sumet;

  pair<float, float> phicorrmet = make_pair( sqrt( metx*metx + mety*mety ), atan2( mety , metx ) );
  return phicorrmet;
}
*/

//VERSION FOR RUNNING ON 8 TEV AND COMBINED 7 TEV
pair<float,float> singleLeptonLooper::getPhiCorrMET( float met, float metphi, float sumet, bool ismc, bool is8TeV ){

  //using met phi correction values from location
  //http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/CMSSW/JetMETCorrections/Type1MET/python/pfMETsysShiftCorrections_cfi.py?revision=1.3&view=markup

  float metx = met * cos( metphi );
  float mety = met * sin( metphi );

  float shiftx = 0.;
  float shifty = 0.;

  //use correction for data vs. mc and 7 TeV vs. 8 TeV
  if (is8TeV) {
    // cout<<"RUNNING WITH 8 TEV SETTINGS! This is a check that you are paying attention! Comment out line "
    // 	<<__LINE__<<" and rerun"<<endl;
    shiftx = ismc ? (+1.77344e-01 - 1.34333e-03*sumet)
      : (-7.67892e-01 + 5.76983e-03*sumet);
    shifty = ismc ? (+8.08402e-01 - 2.84264e-03*sumet)
      : (+5.54005e-01 - 2.94046e-03*sumet);
  } else {
    shiftx = ismc ? (-4.53909e-02 - 2.55863e-05*sumet)
      : (-5.65217e-01 + 5.42436e-03*sumet);
    shifty = ismc ? (+1.27947e-01 - 3.62604e-03*sumet)
      : (+4.54054e-01 - 6.73607e-03*sumet);
  }

  metx -= shiftx;
  mety -= shifty;

  pair<float, float> phicorrmet = make_pair( sqrt( metx*metx + mety*mety ), atan2( mety , metx ) );
  return phicorrmet;
}
