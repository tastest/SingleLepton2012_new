#include "histtools.h"
#include "singleLeptonLooper.h"
#include "TTreeCache.h"
#include "TDatabasePDG.h"
#include "TLorentzVector.h"

#include "../Tools/goodrun.h"
#include "../Tools/vtxreweight.h"
#include "../Tools/msugraCrossSection.h"
#include "BtagFuncs.h"

//#include "stopUtils.h"

bool verbose              = false;
bool doTenPercent         = false;
bool vetoTransition       = true;
bool useOldIsolation      = false;

using namespace std;
using namespace tas;

//--------------------------------------------------------------------

struct SUSYGenParticle { // To be filled with status-3 genParticles
      int pdgId; // PDG identifier (with sign, please)
      int firstMother; // first mother, set to <0 if no mothers
      double energy; // energy [GeV]
      double pt; // pt [GeV]
      double eta; // eta
      double phi; // phi
};

//--------------------------------------------------------------------

int getMotherIndex(int motherid){
  for(int i = 0; i < genps_id().size() ; i++){
    if( motherid == genps_id().at(i) ) return i;
  }

  return -1;
}

// The following reweighting only makes sense for on-shell stop, top and chi0
// In the off-shell case top and anti-top may get very different polarizations
double Reweight_Stop_to_TopChi0 (std::vector<SUSYGenParticle> genParticles, double referenceTopPolarization, double requestedTopPolarization, char* prefix) {

  if( !TString(prefix).Contains("T2") ) return 1.0;

  double weight = 1.;
  int nFoundStops = 0;

  unsigned int ngen = genParticles.size();

  for (unsigned int ig=0; ig<ngen; ++ig) {
    const SUSYGenParticle& gen = genParticles[ig];
    if (gen.firstMother<0) continue;
    if (abs(gen.pdgId)>20) continue; // expect quarks or leptons from W decay

    // Navigate upwards in the stop->top->W->fermion decay chain
    const SUSYGenParticle& genW = genParticles[gen.firstMother];
    if (genW.firstMother<0) continue;
    if (abs(genW.pdgId)!=24) continue;
    const SUSYGenParticle& genTop = genParticles[genW.firstMother];
    if (abs(genTop.pdgId)!=6) continue;

    // We only care about the down-type fermion
    if (genTop.pdgId*gen.pdgId>0) continue;

    // We also need a stop
    if (genTop.firstMother<0) continue;
    const SUSYGenParticle& genStop = genParticles[genTop.firstMother];
    if (abs(genStop.pdgId)!=1000006) continue;

    // Move top and fermion to the stop center-of-mass frame
    TLorentzVector stop4;
    stop4.SetPtEtaPhiE(genStop.pt, genStop.eta, genStop.phi, genStop.energy);
    TVector3 betaV(-stop4.Px()/stop4.Energy(),-stop4.Py()/stop4.Energy(),-stop4.Pz()/stop4.Energy());

    TLorentzVector top4;
    top4.SetPtEtaPhiE(genTop.pt, genTop.eta, genTop.phi, genTop.energy);
    top4.Boost(betaV);

    TLorentzVector ferm4;
    ferm4.SetPtEtaPhiE(gen.pt, gen.eta, gen.phi, gen.energy);
    ferm4.Boost(betaV);

    // Do not reweight if by any reason top/fermion directions are undefined
    // This should be pathological if things are fine
    if (top4.P()<=0 || ferm4.P()<=0) {
      printf("Warning: particles at rest, no weight applied: ptop: %.3e, pf: %.3e\n", top4.P(), ferm4.P());
      continue; 
    }

    double costh = (top4.Px()*ferm4.Px()+top4.Py()*ferm4.Py()+top4.Pz()*ferm4.Pz())/top4.P()/ferm4.P();
      
    double weight_L = (top4.Energy()+top4.P())*(1-costh);
    double weight_R = (top4.Energy()-top4.P())*(1+costh);
    weight *= ((1+requestedTopPolarization)*weight_R+(1-requestedTopPolarization)*weight_L)/((1+referenceTopPolarization)*weight_R+(1-referenceTopPolarization)*weight_L);

    nFoundStops++;
  }

  if( nFoundStops!=2 ) cout << __FILE__ << " " << __LINE__ << " WARNING: found " << nFoundStops << " stops, should be 2." << endl;

  return weight;

};

//--------------------------------------------------------------------

singleLeptonLooper::singleLeptonLooper()
{

  std::cout << " construct " << std::endl;
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

void singleLeptonLooper::InitBaby(){

  weightleft_  = -1.0;
  weightright_ = -1.0;

  //pdf variables
  pdfid1_ = -999;
  pdfid2_ = -999;
  pdfQ_   = -99999.;
  pdfx1_  = -99999.;
  pdfx2_  = -99999.;

  mutrigweight_ = 1.;
  sltrigweight_ = 1.;
  dltrigweight_ = 1.;

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

  // chi2 and MT2
  chi2min_=-1.;
  chi2minprob_=-1.;

  mt2bmin_=-1.;
  mt2blmin_=-1.;
  mt2wmin_=-1.;

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
  
  t1met10s_	=-999.;
  t1met10sphi_	=-999.;
  t1met10smt_	=-999.;

  t1met_off_		= -999.;
  t1metphi_off_		= -999.;
  t1metmt_off_		= -999.;
  t1metphicorr_off_	= -999.;
  t1metphicorrphi_off_	= -999.;
  t1metphicorrmt_off_	= -999.;

  //calomet
  calomet_              =-999.;
  calometphi_           =-999.;

  //trkmet
  trkmet_              =-999.;
  trkmetphi_           =-999.;
  trkmet_nolepcorr_    =-999.;
  trkmetphi_nolepcorr_ =-999.;

  //phi corrected type1 mets
  t1metphicorr_	      =-999.;
  t1metphicorrphi_    =-999.;
  t1metphicorrlep_    =-999.;
  t1metphicorrlepphi_ =-999.;
  t1metphicorrmt_     =-999.;
  t1metphicorrlepmt_  =-999.;
  
  // pfjet vars
  npfjets30_	= 0;
  npfjets30lepcorr_ = 0;
  knjets_       = 1.;

  htpf30_	= 0.;

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
  lep_t_id_	= -1;
  lep_tbar_id_  = -1;
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
  nbs_	        = -1;
  ptjetraw_	= -9999.;
  ptjet23_	= -9999.;
  ptjetF23_	= -9999.;
  ptjetO23_	= -9999.;
  cosphijz_	= -9999.;
  dilep_	= 0;
  jet_		= 0;

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
  mcnu_       =  0;
  mclep_      =  0;
  mcmln_	= -9999.;
  mcmtln_	= -9999.;
  pflepmindrj_   = 9999.;
  pftaudmindrj_  = 9999.;
  lep1pfjetdr_   = 9999.;
  lep2pfjetdr_   = 9999.;

  pfcand5_        = 0;
  pfcand10_       = 0;
  pfcanddir10_       = 0;
  pfcandveto10_       = 0;
  pfcandvetoL10_       = 0;
  pfcandOS10_       = 0;
  pfcandOS10looseZ_       = 0;
  pfcand5looseZ_        = 0;

  pfTau_       = 0;
  pfTau_leadPtcand_  = 0;
  pfTau_leadPtcandID_ = -1;

  pfTau15_       = 0;
  pfTau15_leadPtcand_  = 0;
  pfTau15_leadPtcandID_ = -1;

  pfTauLoose_       = 0;
  pfTauLoose_leadPtcand_  = 0;
  pfTauLoose_leadPtcandID_ = -1;

  pfcandid5_       =-1; 
  pfcandid10_      =-1;
  pfcanddirid10_      =-1;
  pfcandvetoid10_      =-1;
  pfcandvetoLid10_      =-1;
  pfcandidOS10_      =-1;
  pfcandidOS10looseZ_      =-1;
  pfcandid5looseZ_       =-1; 

  pfcandiso5_     = 9999.;     
  pfcandiso10_    = 9999.;     
  pfcanddiriso10_    = 9999.;     
  pfcandvetoiso10_    = 9999.;     
  pfcandvetoLiso10_    = 9999.;     
  pfcandisoOS10_    = 9999.;     
  pfcandisoOS10looseZ_    = 9999.;     
  pfcandiso5looseZ_     = 9999.;     

  pfcandpt5_      = 9999.;
  pfcandpt10_     = 9999.;
  pfcanddirpt10_     = 9999.;
  pfcandvetopt10_     = 9999.;
  pfcandvetoLpt10_     = 9999.;
  pfcandptOS10_     = 9999.;
  pfcandptOS10looseZ_     = 9999.;
  pfcandpt5looseZ_      = 9999.;

  pfcanddz5_      = 9999.;
  pfcanddz10_     = 9999.;
  pfcanddirdz10_     = 9999.;
  pfcandvetodz10_     = 9999.;
  pfcandvetoLdz10_     = 9999.;
  pfcanddzOS10_     = 9999.;
  pfcanddzOS10looseZ_     = 9999.;
  pfcanddz5looseZ_      = 9999.;

  pfcandmindrj5_  = 9999.;
  pfcandmindrj10_ = 9999.;
  pfcanddirmindrj10_ = 9999.;
  pfcandvetomindrj10_ = 9999.;
  pfcandvetoLmindrj10_ = 9999.;

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

  //lepton variables
  iso1_   = -9999; 
  isont1_ = -9999;
  isopfold1_ = -9999;
  isopf1_ = -9999;
  etasc1_ = -9999;
  eoverpin_  = -9999;
  eoverpout_ = -9999;
  dEtaIn_ = -9999;
  dPhiIn_ = -9999;
  sigmaIEtaIEta_ = -9999;
  hOverE_ = -9999;
  ooemoop_ = -9999;
  d0vtx_ = -9999;
  dzvtx_ = -9999;
  expinnerlayers_ = -9999;
  fbrem_ = -9999;
  pfisoch_ = -9999;
  pfisoem_ = -9999;
  pfisonh_ = -9999;
  eSC_ = -9999;
  phiSC_ = -9999;
  eSCRaw_ = -9999;
  eSCPresh_ = -9999;  
  ecalveto1_ = -9999;
  hcalveto1_ = -9999;

  iso2_   = -9999;
  isont2_ = -9999;
  isopf2_ = -9999;
  etasc2_ = -9999;
  eoverpin2_  = -9999;
  eoverpout2_ = -9999;
  dEtaIn2_ = -9999;
  dPhiIn2_ = -9999;
  sigmaIEtaIEta2_ = -9999;
  hOverE2_ = -9999;
  ooemoop2_ = -9999;
  d0vtx2_ = -9999;
  dzvtx2_ = -9999;
  expinnerlayers2_ = -9999;
  fbrem2_ = -9999;
  pfisoch2_ = -9999;
  pfisoem2_ = -9999;
  pfisonh2_ = -9999;
  eSC2_ = -9999;
  phiSC2_ = -9999;
  eSCRaw2_ = -9999;
  eSCPresh2_ = -9999;  
  ecalveto2_ = -9999;
  hcalveto2_ = -9999;

  lep1_scslasercormean_ = -9999.;
  lep1_scslasercormax_  = -9999.;
  lep2_scslasercormean_ = -9999.;
  lep2_scslasercormax_  = -9999.;
  scslasercormax_ = -9999.;
  scslasercormax_pt_ = -9999.;
  scslasercormax_eta_ = -9999.;
  
  //clear vectors
  pfjets_.clear();
  pfjets_genJet_.clear();
  pfjets_csv_.clear();
  pfjets_chm_.clear();
  pfjets_neu_.clear();
  pfjets_l1corr_.clear();
  pfjets_corr_.clear();
  pfjets_mc3_.clear();
  pfjets_mcflavorAlgo_.clear();
  pfjets_mcflavorPhys_.clear();
  pfjets_uncertainty_.clear();

  pfjets_flav_.clear();
  pfjets_lrm_.clear();
  pfjets_lrm2_.clear();
  pfjets_qgtag_.clear();
  pfjets_genJetDr_.clear();
  pfjets_sigma_.clear();
  pfjets_lepjet_.clear();
  pfjets_beta_.clear();
  pfjets_beta2_.clear();
  pfjets_beta_0p1_.clear();
  pfjets_beta_0p2_.clear();
  pfjets_beta2_0p1_.clear();
  pfjets_beta2_0p5_.clear();
  pfjets_mvaPUid_.clear();
  pfjets_mva5xPUid_.clear();
  // pfjets_beta_0p15_.clear();
  // pfjets_beta2_0p15_.clear();
  // pfjets_beta_0p2_.clear();
  // pfjets_beta2_0p2_.clear();
  pfjets_mvaBeta_.clear();

  genps_pdgId_.clear();
  genps_firstMother_.clear();
  genps_energy_.clear();
  genps_pt_.clear();
  genps_eta_.clear();
  genps_phi_.clear();

  genbs_.clear();

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

/* ------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------- */

int singleLeptonLooper::ScanChain(TChain* chain, char *prefix, float kFactor, int prescale, float lumi,
				  FREnum frmode, bool doFakeApp)

{

  //  cout << "ciao " << isData << endl;

  bool isLM = TString(prefix).Contains("LM");
  bool isData = false;
  if( TString(prefix).Contains("data") || TString(prefix).Contains("2012") 
      || TString(prefix).Contains("dimu") || TString(prefix).Contains("diel")
      || TString(prefix).Contains("mueg") ){
    cout << "DATA!!!" << endl;
    isData       = true;
    doTenPercent = false;
  }

  cout << "IS DATA: " << isData << endl;

  if( doTenPercent ) cout << "Processing 10% of MC" << endl;


  //------------------------------------------------------------------------------------------------------
  // set json, vertex reweighting function and msugra cross section files
  //------------------------------------------------------------------------------------------------------
  
  if( !initialized ){

    //set json
    cout << "setting json " << g_json << endl;
    set_goodrun_file( g_json );

    //    if( TString(prefix).Contains("ttall_massivebin") ) 
    set_vtxreweight_rootfile("vtxreweight/vtxreweight_Summer12MC_PUS10_19fb_Zselection.root",true);

    //   weight3D_init( "vtxreweight/Weight3D.root" );

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
    jetcorr_filenames_pfL1FastJetL2L3.push_back  ("jetCorrections/GR_P_V39_AN3_L1FastJet_AK5PF.txt");
    jetcorr_filenames_pfL1FastJetL2L3.push_back  ("jetCorrections/GR_P_V39_AN3_L2Relative_AK5PF.txt");
    jetcorr_filenames_pfL1FastJetL2L3.push_back  ("jetCorrections/GR_P_V39_AN3_L3Absolute_AK5PF.txt");
    jetcorr_filenames_pfL1FastJetL2L3.push_back  ("jetCorrections/GR_P_V39_AN3_L2L3Residual_AK5PF.txt");

    pfUncertaintyFile = "jetCorrections/GR_P_V39_AN3_Uncertainty_AK5PF.txt";
  } 
  else {
    jetcorr_filenames_pfL1FastJetL2L3.push_back  ("jetCorrections/DESIGN53_V15_L1FastJet_AK5PF.txt");
    jetcorr_filenames_pfL1FastJetL2L3.push_back  ("jetCorrections/DESIGN53_V15_L2Relative_AK5PF.txt");
    jetcorr_filenames_pfL1FastJetL2L3.push_back  ("jetCorrections/DESIGN53_V15_L3Absolute_AK5PF.txt");
    
    pfUncertaintyFile = "jetCorrections/DESIGN53_V15_Uncertainty_AK5PF.txt";
  }

  jet_corrector_pfL1FastJetL2L3  = makeJetCorrector(jetcorr_filenames_pfL1FastJetL2L3);

  JetCorrectionUncertainty *pfUncertainty   = new JetCorrectionUncertainty( pfUncertaintyFile   );

  MetCorrector *met_corrector_pfL1FastJetL2L3 = new MetCorrector(jetcorr_filenames_pfL1FastJetL2L3);

  /*
   *  Jet Smearer Object to obtain the jet pt uncertainty.
   */

  std::vector<std::string> list_of_file_names;
  list_of_file_names.push_back("jetSmearData/Spring10_PtResolution_AK5PF.txt");
  list_of_file_names.push_back("jetSmearData/Spring10_PhiResolution_AK5PF.txt");
  list_of_file_names.push_back("jetSmearData/jet_resolutions.txt");
  JetSmearer *jetSmearer = makeJetSmearer(list_of_file_names);
 
  QGLikelihoodCalculator *qglikeli_ = new QGLikelihoodCalculator("QGTaggerConfig/QGTaggerConfig_nCharged_AK5PF.txt",
								 "QGTaggerConfig/QGTaggerConfig_nNeutral_AK5PF.txt",
								 "QGTaggerConfig/QGTaggerConfig_ptD_AK5PF.txt");

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


  cout << " done with initialization "  << endl;
  
  unsigned int nEventsChain = chain->GetEntries();
  unsigned int nEventsTotal = 0;
  // map isn't needed for this purpose, vector is sufficient
  // better would be to use a struct with run, lb, event
  map<int,int> m_events;

  // loop over files
  TObjArray *listOfFiles = chain->GetListOfFiles();
  TIter fileIter(listOfFiles);
  TChainElement* currentFile = 0;

  // test chain
  if (!chain)
    {
      throw std::invalid_argument("at::ScanChain: chain is NULL!");
    }
  if (chain->GetListOfFiles()->GetEntries()<1)
    {
      throw std::invalid_argument("at::ScanChain: chain has no files!");
    }
  if (not chain->GetFile())
    {
      throw std::invalid_argument("at::ScanChain: chain has no files or file path is invalid!");
    }
  int nSkip_els_conv_dist = 0;

  float netot  = 0.;
  float nmtot  = 0.;
  float nepass = 0.;
  float nmpass = 0.;

  if(g_createTree) makeTree(prefix, doFakeApp, frmode);

  while((currentFile = (TChainElement*)fileIter.Next())) {
    TFile* f = new TFile(currentFile->GetTitle());

    cout << currentFile->GetTitle() << endl;

    if (!f || f->IsZombie()) {
      throw std::runtime_error(Form("ERROR::File from TChain is invalid or corrupt: %s", currentFile->GetTitle()));
    }
    
    // get the trees in each file
    // TTree *tree = (TTree*)f->Get("Events");
    TTree *tree = dynamic_cast<TTree*>(f->Get("Events"));
    if (!tree || tree->IsZombie()) {
      throw std::runtime_error(Form("ERROR::File from TChain has an invalid TTree or is corrupt: %s", currentFile->GetTitle()));
    }

    //Matevz
    TTreeCache::SetLearnEntries(100);
    tree->SetCacheSize(128*1024*1024);

    cms2.Init(tree);
      
    unsigned int nEntries = tree->GetEntries();

    for(unsigned int z = 0; z < nEntries; ++z) {
      ++nEventsTotal;

      /////////      cout << nEventsTotal << endl;

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

      if( evt_ww_rho_vor() != evt_ww_rho_vor() ){
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
      //      cout<<"dataset: "<<datasetname.Data()<<" isperiodA: "<<isperiodA<<endl;

      // skip stop-pair events with m(stop) > 850 GeV
      // if( TString(prefix).Contains("T2") ){
      // 	if( sparm_mG() > 600.0 ) continue;
      // }

      //---------------------------------------------
      // event cleaning and good run list
      //---------------------------------------------

      if( !cleaning_goodVertexApril2011() )                          continue;
      if( isData && !goodrun(evt_run(), evt_lumiBlock()) ) continue;

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
            
      //if( TString(prefix).Contains("T2") ) useOldIsolation = true;
      for( unsigned int iel = 0 ; iel < els_p4().size(); ++iel ){

	//-------------------------------------------
	// FASTSIM BUG: check for duplicate electrons
	//-------------------------------------------

	bool foundDuplicate = false;

	for( unsigned int i2 = 0 ; i2 < goodLeptons.size() ; ++i2 ){
	  if( fabs( els_p4().at(iel).pt()  - goodLeptons.at(i2).pt()  ) < 0.001 &&
	      fabs( els_p4().at(iel).eta() - goodLeptons.at(i2).eta() ) < 0.001 &&
	      fabs( els_p4().at(iel).phi() - goodLeptons.at(i2).phi() ) < 0.001 ){
	    eldup_ = 1;
	    //cout << "WARNING! FOUND DUPLICATE ELECTRON " << evt_run() << " " << evt_lumiBlock() << " " << evt_event() << endl;
	    foundDuplicate = true;
	  }
	}
	if( foundDuplicate ) continue;

	if( els_p4().at(iel).pt() < 10 )                                                                continue;
        if( !passElectronSelection_Stop2012_v3( iel , vetoTransition,vetoTransition,useOldIsolation) )  continue;

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
	lepchi2ndf.push_back( mus_gfit_chi2().at(imu)/mus_gfit_ndof().at(imu) );
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
      LorentzVector minvjet1;

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
      LorentzVector minvjet2;

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
	if( ROOT::Math::VectorUtil::DeltaR( *lep1_ , els_p4().at(iel) ) < 0.1 )                        continue;
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
      lep_t_	= 0;   
      lep_tbar_	= 0;   
      stop_t_	= 0;   
      stop_tbar_ = 0;   
      neutralino_t_ = 0;   
      neutralino_tbar_ = 0;  

      npartons_    =  0;
      nwzpartons_  = -9;
      maxpartonpt_ = -1;

      mgcor_ = 1.0;
      wflav_ = -1;

      if( !isData ){

	bool foundwz = false;
	bool foundlep = false;
	bool foundnu  = false;

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
	if( TString(prefix).Contains("ttall") ||
	    TString(prefix).Contains("tt_") ){
	  if( nleps == 0 ) mgcor_ = 1.028;
	  if( nleps == 1 ) mgcor_ = 0.986;
	  if( nleps == 2 ) mgcor_ = 0.945;
	}
	if( TString(prefix).Contains("powheg") ||
	    TString(prefix).Contains("sherpa") ) 
	  mgcor_ = 1.0;

	if( strcmp(prefix,"ttem")  == 0 && ( nels + nmus ) != 2 ) continue;
	if( strcmp(prefix,"ttdil") == 0 && nleps != 2           ) continue;
	if( strcmp(prefix,"ttotr") == 0 && nleps == 2           ) continue;
	
	LorentzVector vdilepton(0,0,0,0);
	LorentzVector vttbar(0,0,0,0);
	int ntops = 0;
	nbs_ = 0;

	for ( int igen = 0 ; igen < (int)genps_id().size() ; igen++ ) { 
	  if ( abs( genps_id().at(igen) ) == 11) vdilepton += genps_p4().at(igen); 
	  if ( abs( genps_id().at(igen) ) == 13) vdilepton += genps_p4().at(igen); 

	  int id = genps_id().at(igen);
	  int pid = abs( genps_id().at(igen) );
	  int mothid = abs(genps_id_mother().at(igen));

	  // only count/store massive b quarks
	  // - in massless b samples, final state b's will have the mass set to 4.8 GeV
	  //    while initial state b's have mass 0
	  if ( abs(id) == 5 ) {
	    // mass (from dumpDocLines)
	    float m2 = cms2.genps_p4().at(igen).M2();
	    float m  = m2 >= 0 ? sqrt(m2) : 0.0;
	    if (m > 0.) {
	      ++nbs_;
	      genbs_.push_back(genps_p4().at(igen));
	    }
	  }

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

	  //store stop
	  if ( id == 1000006)
	    stop_t_ = &(genps_p4().at(igen));   
	  else if ( id == -1000006 )
	    stop_tbar_ = &(genps_p4().at(igen));   

    //store neutralino
    if ( genps_id_mother().at(igen) == 1000006  && ( abs(id) == 1000022 ) ) {
      neutralino_t_ = &(genps_p4().at(igen));
    }
    if ( genps_id_mother().at(igen) == -1000006 && ( abs(id) == 1000022 ) ) {
      neutralino_tbar_ = &(genps_p4().at(igen));
    }

	  //store daughter lepton
	  if ( abs(mothid) == 24 && (abs(id) == 11 || abs(id) == 13 || abs(id) ==15)) {

	    if (genps_id_mother().at(igen)>0) {
	      // lept 1 is the particle 
	      lep_t_id_ = genps_id().at(igen);
	      lep_t_ = &(genps_p4().at(igen));
	    } else {
	      // lept 2 is the anti-particle
	      lep_tbar_id_ = genps_id().at(igen);
	      lep_tbar_ = &(genps_p4().at(igen));
	    }
	  }

	  // store lepton, neutrino and W for single lepton events                     
	  if (pid==11 || pid==13) {
	    foundlep = true;
	    mclep_ = &genps_p4()[igen];
	  }
	  if (pid==12 || pid==14) {
	    foundnu = true;
	    mcnu_  = &genps_p4()[igen];
	  }

	  // store W or Z pT 
	  // ignoring cases where have more than 1 boson for now
	  if ( pid == 24 ) {
	    ptwgen_ = genps_p4().at(igen).pt();
	    foundwz = true;
	    nwzpartons_  = 0;
	  }
	  if ( pid == 23 ) {
	    ptzgen_ = genps_p4().at(igen).pt();
	    foundwz = true;
	    nwzpartons_  = 0;
	  }

	  if (foundwz && ( pid == 1 || pid == 2 || pid == 3 || pid == 4 || pid == 5 || pid == 6 || pid == 21 ) )   
	    nwzpartons_++;

	  // skip lines up to t and tbar
	  if( igen < 8 ) continue;

	  // require particle is a quark or a gluon
	  if( !( pid==1 || pid==2 || pid==3 || pid==4 || pid==5 || pid==6 || pid == 21 ) ) continue;

	  // require mother is not a top or W
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

	if (foundlep && foundnu) {
	  mcmln_ = (*mclep_+*mcnu_).mass();
	  mcmtln_ = getMT( mclep_->Pt() , mclep_->Phi() , mcnu_->Pt() , mcnu_->Phi() );
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

	//-----------------------------------------------------
	// store single gen lepton info
	//-----------------------------------------------------

	if( nleps_ == 1 ){

	  int nfoundleps = 0;

	  for ( int igen = 0 ; igen < (int)genps_id().size() ; igen++ ) { 

	    int id = genps_id().at(igen);

	    if( !( abs(id)==11 || abs(id)==13 || abs(id)==15 ) ) continue;

	    nfoundleps++;
	    mcid1_   = id;
	    mclep1_  = &genps_p4().at(igen);
	    mcdr1_   = ROOT::Math::VectorUtil::DeltaR( *lep1_ , *mclep1_ );

	    if( abs(id)==15 ){
	      mcdecay1_ = 1;
	      //variables to calculate the visible tau energy
	      LorentzVector taudvis; 
	      for(unsigned int kk = 0; kk < genps_lepdaughter_id().at(igen).size(); kk++) {
		int daughter = abs(genps_lepdaughter_id()[igen][kk]);
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
		    mctaudid1_ = genps_lepdaughter_id()[igen][kk];
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
	  
	  for ( int igen = 0 ; igen < (int)genps_id().size() ; igen++ ) { 

	    int id = genps_id().at(igen);

	    if( !( abs(id)==11 || abs(id)==13 || abs(id)==15 ) ) continue;

	    nfoundleps++;

	    if( ROOT::Math::VectorUtil::DeltaR( *lep1_ , genps_p4().at(igen) ) < drmin ){
	      drmin   = ROOT::Math::VectorUtil::DeltaR( *lep1_ , genps_p4().at(igen) );
	      igenmin = igen;
	    }
	  }

	  if( nfoundleps != 2 ) cout << "ERROR! expected 2 leptons, found " << nfoundleps << endl;

	  //-----------------------------------------------------
	  // store info for closest gen lepton
	  //-----------------------------------------------------

	  mcid1_   = genps_id().at(igenmin);
	  mclep1_  = &genps_p4().at(igenmin);
	  mcdr1_   = ROOT::Math::VectorUtil::DeltaR( *lep1_ , *mclep1_ );

	  if( abs(mcid1_)==15 ){
	    mcdecay1_ = 1;
	    
	    //variables to calculate the visible tau energy
	    LorentzVector taudvis; 
	    
	    for(unsigned int kk = 0; kk < genps_lepdaughter_id().at(igenmin).size(); kk++) {
	      int daughter = abs(genps_lepdaughter_id()[igenmin][kk]);
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
		  mctaudid1_ = genps_lepdaughter_id()[igenmin][kk];
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

	  for ( int igen = 0 ; igen < (int)genps_id().size() ; igen++ ) { 

	    if( igen == igenmin ) continue; //skip closest lepton

	    int id = genps_id().at(igen);

	    if( !( abs(id)==11 || abs(id)==13 || abs(id)==15 ) ) continue;

	    igenmin2 = igen;

	    mcid2_   = id;
	    mclep2_  = &genps_p4().at(igen);
	    if( ngoodlep_ > 1 ) mcdr2_   = ROOT::Math::VectorUtil::DeltaR( *lep2_ , *mclep2_ );

	    if( abs(id)==15 ){
	      mcdecay2_ = 1;
	      //variables to calculate the visible tau energy
	      LorentzVector taudvis; 
	      
	      for(unsigned int kk = 0; kk < genps_lepdaughter_id().at(igen).size(); kk++) {
		int daughter = abs(genps_lepdaughter_id()[igen][kk]);


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
		    mctaudid2_ = genps_lepdaughter_id()[igen][kk];
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
	      mleppassid_   = passElectronSelection_Stop2012_v3_NoIso( imatch , vetoTransition,vetoTransition,useOldIsolation) ? 1 : 0;
	      mleppassiso_  = passElectronSelection_Stop2012_v3_Iso  ( imatch , vetoTransition,vetoTransition,useOldIsolation) ? 1 : 0;
	      mlepiso_      = electronIsoValuePF2012_FastJetEffArea_v3( imatch , 0.3 , 0 );
	    }

	    // found matched muon
	    else if( ID == 2 ){
	      mlepid_       = 13 * mus_charge().at(imatch);
	      mlep_         = &mus_p4().at(imatch);
	      mleppassid_   = muonIdNotIsolated( imatch , ZMet2012_v1 )   ? 1 : 0;
	      mleppassiso_  = muonIsoValuePF2012_deltaBeta(imatch) < 0.15 ? 1 : 0;
	      mlepiso_      = muonIsoValuePF2012_deltaBeta(imatch);
	    }

	    mlepdr_ = ROOT::Math::VectorUtil::DeltaR( *mlep_ , *mclep2_ );
	  }
	  
	}

	/*
	else if( nleps_ < 0 || nleps_ > 2 ){
	  cout << "ERROR nleptons = " << nleps_ << endl;
	}
	*/


	//------------------------------------------
	// stop reweighting code
	//------------------------------------------

	std::vector<SUSYGenParticle> genParticles;

	
	for (unsigned int ig=0; ig<genps_id().size(); ++ig) {

	  if( genps_status().at(ig) != 3 ) continue;

	  SUSYGenParticle part;
	  part.pdgId       = genps_id().at(ig);
	  part.energy      = genps_p4().at(ig).E();
	  part.pt          = genps_p4().at(ig).pt();
	  part.eta         = genps_p4().at(ig).eta();
	  part.phi         = genps_p4().at(ig).phi();
	  part.firstMother = getMotherIndex( genps_id_mother().at(ig) );
	  genParticles.push_back(part);

	  genps_pdgId_      .push_back(part.pdgId);
	  genps_firstMother_.push_back(part.firstMother);
	  genps_energy_     .push_back(part.energy);
	  genps_pt_         .push_back(part.pt);
	  genps_eta_        .push_back(part.eta);
	  genps_phi_        .push_back(part.phi);
	}

	weightleft_  = Reweight_Stop_to_TopChi0 (genParticles, 0., -1, prefix);
	weightright_ = Reweight_Stop_to_TopChi0 (genParticles, 0.,  1, prefix);
      }

      /*
      for (unsigned int ipf = 0; ipf < pfcands_p4().size(); ipf++) {

	if( pfcands_charge().at(ipf) == 0   ) continue;

 	int itrk = pfcands_trkidx().at(ipf);
	
 	if( itrk < (int)trks_trk_p4().size() && itrk >= 0 ){
 	  if( fabs( trks_dz_pv(itrk,0).first ) > 0.2 ){
 	    fillOverFlow( h_PU_trkpt , pfcands_p4().at(ipf).pt() );
 	  }
 	}
      }
      */

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
	for (unsigned int ipf = 0; ipf < pfcands_p4().size(); ipf++) {
	  
	  //	  if( pfcands_p4().at(ipf).pt() < 5   ) continue;
	  if( pfcands_charge().at(ipf) == 0   ) continue;
	  
	  int itrk = pfcands_trkidx().at(ipf);
	  
	  if( itrk < (int)trks_trk_p4().size() && itrk >= 0 ){
	    if( fabs( trks_dz_pv(itrk,0).first ) > 0.2 ) continue;
	  }
	  //Only remove leading lepton to see what happens to the sub-leading lepton
	  // bool isGoodLepton = false;
	  // for( int ilep = 0 ; ilep < goodLeptons.size() ; ilep++ ){
	  //   if( ROOT::Math::VectorUtil::DeltaR( pfcands_p4().at(ipf) , goodLeptons.at(ilep) ) < 0.1 ) 
	  //     isGoodLepton = true;  
	  // }
	  // if( isGoodLepton ) continue;
	  if( ROOT::Math::VectorUtil::DeltaR( pfcands_p4().at(ipf) , goodLeptons.at(imaxpt) ) < 0.1 ) continue;

	  //Store highest pT track within match radius or if none is found closest pfcand
	  float matchR = 0.15;
	  float drpf = ROOT::Math::VectorUtil::DeltaR( pfcands_p4().at(ipf) , *mclep2_ );

	  struct myTrackIso myTrackIso=trackIso(ipf);
	  float iso = myTrackIso.iso_dr03_dz005_pt00 / pfcands_p4().at(ipf).pt();

	  //	  float iso = trackIso(ipf) / pfcands_p4().at(ipf).pt();
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
	  float taudrpf = ROOT::Math::VectorUtil::DeltaR( pfcands_p4().at(ipf) , *mctaud2_ );
	  float tauiso = myTrackIso.iso_dr03_dz005_pt00 / pfcands_p4().at(ipf).pt();
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
      float dz_cut_loose = 0.1;
      //------------------------------------------------------
      // store pt and iso for most isolated track (pt>10 GeV) and (pt>5 GeV)
      //------------------------------------------------------

      trkpt10_           = -1.0;
      trkreliso10_       = 1000.;
      mleptrk10_         = -1.0;
      float miniso10     = 999;
      trkpt10loose_      = -1.0;
      trkreliso10loose_  = 1000.;

      trkpt5_          = -1.0;
      trkreliso5_      = 1000.;
      mleptrk5_        = -1.0;
      float miniso5    = 999;
      trkpt5loose_     = -1.0;
      trkreliso5loose_ = 1000.;


      for (unsigned int ipf = 0; ipf < pfcands_p4().size(); ipf++) {

	if( pfcands_p4().at(ipf).pt() < 5  ) continue;
	if( pfcands_charge().at(ipf) == 0   ) continue;

	struct myTrackIso myTrackIso=trackIso(ipf);

	// 	int itrk = pfcands_trkidx().at(ipf);
	
	bool isGoodLepton = false;
	for( int ilep = 0 ; ilep < (int)goodLeptons.size() ; ilep++ ){
	  if( ROOT::Math::VectorUtil::DeltaR( pfcands_p4().at(ipf) , goodLeptons.at(ilep) ) < 0.1 ) 
	    isGoodLepton = true;  
	}
	bool isLeadLepton = ( ROOT::Math::VectorUtil::DeltaR( pfcands_p4().at(ipf) , 
						goodLeptons.at(imaxpt) ) < 0.1 ) ? true : false;

	/// this is with the OS requirement
	float charge=(pfcands_particleId().at(ipf)*id1_);
	// charge < 0 is a SS , charge > 0 is a OS for e/mu; need to flip for pions
	if((abs(pfcands_particleId().at(ipf))!=11) && (abs(pfcands_particleId().at(ipf))!=13)) charge*=(-1); 

	//////////

	int itrk = -1;
	float mindz=999;

	if (abs(pfcands_particleId().at(ipf))!=11) {
	  itrk = pfcands_trkidx().at(ipf);
	  if( itrk >= (int)trks_trk_p4().size() || itrk < 0 ) continue;
	  mindz=trks_dz_pv(itrk,0).first;
	}

	if (abs(pfcands_particleId().at(ipf))==11 && pfcands_pfelsidx().at(ipf)>=0) {
	  itrk = els_gsftrkidx().at(pfcands_pfelsidx().at(ipf));
	  if( itrk >= (int)gsftrks_p4().size() || itrk < 0 ) continue;
	  mindz=gsftrks_dz_pv(itrk,0).first;
	}

	// start with loose dz cut
	if( fabs(mindz) > dz_cut_loose ) continue;

	//store loose definition to compare with previous results
	float iso = myTrackIso.iso_dr03_dz010_pt00 / pfcands_p4().at(ipf).pt();

	if(pfcands_p4().at(ipf).pt()>=10 && iso < trkreliso10loose_ && !isGoodLepton ){
	  trkpt10loose_       = pfcands_p4().at(ipf).pt();
	  trkreliso10loose_   = iso;
	}


	if( pfcands_p4().at(ipf).pt()>=5 && iso < trkreliso5loose_ && !isGoodLepton ){
	  trkpt5loose_     = pfcands_p4().at(ipf).pt();
	  trkreliso5loose_ = iso;
	}

	if( pfcands_p4().at(ipf).pt()>=5 && iso < pfcandiso5looseZ_ && !isLeadLepton){
	  pfcandiso5looseZ_ = iso;
	  pfcanddz5looseZ_ = mindz;
	  pfcandpt5looseZ_ = pfcands_p4().at(ipf).pt();
	  pfcand5looseZ_ = &pfcands_p4().at(ipf);
          pfcandid5looseZ_ =  pfcands_particleId().at(ipf);
	}

	if( pfcands_p4().at(ipf).pt()>=10 && iso < pfcandisoOS10looseZ_ && !isLeadLepton && charge>0){
	  pfcandisoOS10looseZ_ = iso;
	  pfcanddzOS10looseZ_ = mindz;
	  pfcandptOS10looseZ_ = pfcands_p4().at(ipf).pt();
	  pfcandOS10looseZ_ = &pfcands_p4().at(ipf);
	  pfcandidOS10looseZ_ =  pfcands_particleId().at(ipf);
	}

	//tighten dz cut
	if( fabs(mindz) > dz_cut ) continue;

	iso=myTrackIso.isoDir_dr03_dz005_pt00/pfcands_p4().at(ipf).pt();

	if( pfcands_p4().at(ipf).pt()>=10 && iso < pfcanddiriso10_ && !isLeadLepton ){

	  pfcanddiriso10_ = iso;
          pfcanddirpt10_ = pfcands_p4().at(ipf).pt();
          pfcanddir10_ = &pfcands_p4().at(ipf);
          pfcanddirid10_ =  pfcands_particleId().at(ipf);

	}

	// with veto cone
	iso = myTrackIso.iso_dr0503_dz005_pt00 / pfcands_p4().at(ipf).pt();

	if( pfcands_p4().at(ipf).pt()>=10 && iso < pfcandvetoiso10_ && !isLeadLepton ){

	  pfcandvetoiso10_ = iso;
          pfcandvetopt10_ = pfcands_p4().at(ipf).pt();
          pfcandveto10_ = &pfcands_p4().at(ipf);
          pfcandvetoid10_ =  pfcands_particleId().at(ipf);

	}

	// with veto cone
	iso = myTrackIso.iso_dr01503_dz005_pt00 / pfcands_p4().at(ipf).pt();

	if( pfcands_p4().at(ipf).pt()>=10 && iso < pfcandvetoLiso10_ && !isLeadLepton ){

	  pfcandvetoLiso10_ = iso;
          pfcandvetoLpt10_ = pfcands_p4().at(ipf).pt();
          pfcandvetoL10_ = &pfcands_p4().at(ipf);
          pfcandvetoLid10_ =  pfcands_particleId().at(ipf);

	}

	//recalculated definition of the isolation with the default values
	iso = myTrackIso.iso_dr03_dz005_pt00 / pfcands_p4().at(ipf).pt();

	if( pfcands_p4().at(ipf).pt()>=10 && iso < miniso10 && !isGoodLepton ){
	  miniso10       = iso;
	  trkpt10_       = pfcands_p4().at(ipf).pt();
	  mleptrk10_     = (*lep1_+pfcands_p4().at(ipf)).pt();
	  trkreliso10_   = iso;
	}

	if( pfcands_p4().at(ipf).pt()>=5 && iso < miniso5 && !isGoodLepton ){
	  miniso5     = iso;
	  trkpt5_     = pfcands_p4().at(ipf).pt();
	  mleptrk5_   = (*lep1_+pfcands_p4().at(ipf)).pt();
	  trkreliso5_ = iso;
	  //itrk       = ipf;
	}

	if (isLeadLepton) continue;

	if( pfcands_p4().at(ipf).pt()>=5 && iso < pfcandiso5_){
	  pfcandiso5_ = iso;
	  pfcanddz5_ = mindz;
	  pfcandpt5_ = pfcands_p4().at(ipf).pt();
	  pfcand5_ = &pfcands_p4().at(ipf);
          pfcandid5_ =  pfcands_particleId().at(ipf);
	}

	if( pfcands_p4().at(ipf).pt() < 10  ) continue;

	if( iso < pfcandisoOS10_ && charge>0){
	  pfcandisoOS10_ = iso;
	  pfcanddzOS10_ = mindz;
	  pfcandptOS10_ = pfcands_p4().at(ipf).pt();
	  pfcandOS10_ = &pfcands_p4().at(ipf);
	  pfcandidOS10_ =  pfcands_particleId().at(ipf);
	}

	/// this is the default case Track with dz<0.05 ; pt>10 ; notLeadingLepton
	if( iso < pfcandiso10_ ){
	  pfcandiso10_ = iso;
	  pfcanddz10_ = mindz;
	  pfcandpt10_ = pfcands_p4().at(ipf).pt();
	  pfcand10_ = &pfcands_p4().at(ipf);
	  pfcandid10_ =  pfcands_particleId().at(ipf);
	}

	/// Pt scan

	float iso0p1 = myTrackIso.iso_dr03_dz005_pt01 / pfcands_p4().at(ipf).pt();
	if( iso0p1 < trkreliso10pt0p1_ && !isGoodLepton ){
	  trkpt10pt0p1_       = pfcands_p4().at(ipf).pt();
	  trkreliso10pt0p1_   = iso0p1;
	}
	if( iso0p1 < pfcandiso10pt0p1_ ){
	  pfcandpt10pt0p1_ = pfcands_p4().at(ipf).pt();
	  pfcandiso10pt0p1_ = iso0p1;
	}

	float iso0p2 = myTrackIso.iso_dr03_dz005_pt02  / pfcands_p4().at(ipf).pt();
	if( iso0p2 < trkreliso10pt0p2_ && !isGoodLepton ){
	  trkpt10pt0p2_       = pfcands_p4().at(ipf).pt();
	  trkreliso10pt0p2_   = iso0p2;
	}
	if( iso0p2 < pfcandiso10pt0p2_ ){
	  pfcandpt10pt0p2_ = pfcands_p4().at(ipf).pt();
	  pfcandiso10pt0p2_ = iso0p2;
	}

	float iso0p3 =myTrackIso.iso_dr03_dz005_pt03  / pfcands_p4().at(ipf).pt();
	if( iso0p3 < trkreliso10pt0p3_ && !isGoodLepton ){
	  trkpt10pt0p3_       = pfcands_p4().at(ipf).pt();
	  trkreliso10pt0p3_   = iso0p3;
	}
	if( iso0p3 < pfcandiso10pt0p3_ ){
	  pfcandpt10pt0p3_ = pfcands_p4().at(ipf).pt();
	  pfcandiso10pt0p3_ = iso0p3;
	}

	float iso0p4 = myTrackIso.iso_dr03_dz005_pt04  / pfcands_p4().at(ipf).pt();
	if( iso0p4 < trkreliso10pt0p4_ && !isGoodLepton ){
	  trkpt10pt0p4_       = pfcands_p4().at(ipf).pt();
	  trkreliso10pt0p4_   = iso0p4;
	}
	if( iso0p4 < pfcandiso10pt0p4_ ){
	  pfcandpt10pt0p4_ = pfcands_p4().at(ipf).pt();
	  pfcandiso10pt0p4_ = iso0p4;
	}

	float iso0p5 = myTrackIso.iso_dr03_dz005_pt05  / pfcands_p4().at(ipf).pt();
	if( iso0p5 < trkreliso10pt0p5_ && !isGoodLepton ){
	  trkpt10pt0p5_       = pfcands_p4().at(ipf).pt();
	  trkreliso10pt0p5_   = iso0p5;
	}
	if( iso0p5 < pfcandiso10pt0p5_ ){
	  pfcandpt10pt0p5_ = pfcands_p4().at(ipf).pt();
	  pfcandiso10pt0p5_ = iso0p5;
	}

	float iso0p6 = myTrackIso.iso_dr03_dz005_pt06  / pfcands_p4().at(ipf).pt();
	if( iso0p6 < trkreliso10pt0p6_ && !isGoodLepton ){
	  trkpt10pt0p6_       = pfcands_p4().at(ipf).pt();
	  trkreliso10pt0p6_   = iso0p6;
	}
	if( iso0p6 < pfcandiso10pt0p6_ ){
	  pfcandpt10pt0p6_ = pfcands_p4().at(ipf).pt();
	  pfcandiso10pt0p6_ = iso0p6;
	}

	float iso0p7 = myTrackIso.iso_dr03_dz005_pt07  / pfcands_p4().at(ipf).pt();
	if( iso0p7 < trkreliso10pt0p7_ && !isGoodLepton ){
	  trkpt10pt0p7_       = pfcands_p4().at(ipf).pt();
	  trkreliso10pt0p7_   = iso0p7;
	}
	if( iso0p7 < pfcandiso10pt0p7_ ){
	  pfcandpt10pt0p7_ = pfcands_p4().at(ipf).pt();
	  pfcandiso10pt0p7_ = iso0p7;
	}

	float iso0p8 = myTrackIso.iso_dr03_dz005_pt08  / pfcands_p4().at(ipf).pt();
	if( iso0p8 < trkreliso10pt0p8_ && !isGoodLepton ){
	  trkpt10pt0p8_       = pfcands_p4().at(ipf).pt();
	  trkreliso10pt0p8_   = iso0p8;
	}
	if( iso0p8 < pfcandiso10pt0p8_ ){
	  pfcandpt10pt0p8_ = pfcands_p4().at(ipf).pt();
	  pfcandiso10pt0p8_ = iso0p8;
	}

	float iso0p9 = myTrackIso.iso_dr03_dz005_pt09  / pfcands_p4().at(ipf).pt();
	if( iso0p9 < trkreliso10pt0p9_ && !isGoodLepton ){
	  trkpt10pt0p9_       = pfcands_p4().at(ipf).pt();
	  trkreliso10pt0p9_   = iso0p9;
	}
	if( iso0p9 < pfcandiso10pt0p9_ ){
	  pfcandpt10pt0p9_ = pfcands_p4().at(ipf).pt();
	  pfcandiso10pt0p9_ = iso0p9;
	}

	float iso1p0 = myTrackIso.iso_dr03_dz005_pt10  / pfcands_p4().at(ipf).pt();
	if( iso1p0 < trkreliso10pt1p0_ && !isGoodLepton ){
	  trkpt10pt1p0_       = pfcands_p4().at(ipf).pt();
	  trkreliso10pt1p0_   = iso1p0;
	}
	if( iso1p0 < pfcandiso10pt1p0_ ){
	  pfcandpt10pt1p0_ = pfcands_p4().at(ipf).pt();
	  pfcandiso10pt1p0_ = iso1p0;
	}


      }

      //----------------------------------------
      // TAU
      //----------------------------------------

      int indexTauMax=-1;
      double ptTauMax=0.;

      int indexTauLooseMax=-1;
      double ptTauLooseMax=0.;

      int indexTauMax15=-1;
      double ptTauMax15=0.;


      for (unsigned int itau=0; itau < taus_pf_p4().size(); itau++) {

	if(taus_pf_p4().at(itau).pt()<15) continue;
	/// Use only OS charge
	if((taus_pf_charge().at(itau)*id1_)>0) continue;

	///
	bool isLeadLepton = ( ROOT::Math::VectorUtil::DeltaR( taus_pf_p4().at(itau) ,
							      goodLeptons.at(imaxpt) ) < 0.4 ) ? true : false;

	if(isLeadLepton) continue;
	if(!taus_pf_byDecayModeFinding().at(itau)) continue;


	// isolation Medium ; pt > 15
	if(taus_pf_byMediumIsolationMVA2().at(itau)) {
	  if(taus_pf_p4().at(itau).pt()>ptTauMax15) {
	    ptTauMax15 = taus_pf_p4().at(itau).pt();
	    indexTauMax15 = itau;
	  }	
	}


	if(taus_pf_p4().at(itau).pt()<20) continue;

	// isolation Medium ; pt > 20    
	if(taus_pf_byMediumIsolationMVA2().at(itau)) {
	  if(taus_pf_p4().at(itau).pt()>ptTauMax) {
	    ptTauMax = taus_pf_p4().at(itau).pt();
	    indexTauMax = itau;
	  }	
	}

	// isolation Loose ; pt > 20
	if(taus_pf_byLooseIsolationMVA2().at(itau)) {
	  if(taus_pf_p4().at(itau).pt()>ptTauLooseMax) {
	    ptTauLooseMax = taus_pf_p4().at(itau).pt();
	    indexTauLooseMax = itau;
	  }	
	}

	/*
	std::cout << "pt "<< taus_pf_p4().at(itau).pt() << " decayModeFinding " << taus_pf_byDecayModeFinding().at(itau)
		  << " isoMVA2 " << taus_pf_byMediumIsolationMVA2().at(itau) 
		  << " pfcandSize " << taus_pf_pfcandIndicies().at(itau).size() << std::endl;
	*/

      }


      if(indexTauMax15!=-1) {

	pfTau15_=&taus_pf_p4().at(indexTauMax15);
	for(int ipf=0; ipf<taus_pf_pfcandIndicies().at(indexTauMax15).size(); ipf ++) {	  

	  int index=(taus_pf_pfcandIndicies().at(indexTauMax15)).at(ipf);
	  //	  cout << "part  " << ipf << " with pt " << pfcands_p4().at(index).pt() << " eta " << pfcands_p4().at(index).eta()  << endl;

	}	

	if(taus_pf_pfcandIndicies().at(indexTauMax15).size()>0) {
	  int leadingPtCand_index=(taus_pf_pfcandIndicies().at(indexTauMax15)).at(0);
	  pfTau15_leadPtcand_= &(pfcands_p4().at(leadingPtCand_index));
	  pfTau15_leadPtcandID_= pfcands_particleId().at(leadingPtCand_index);

	}
      }


      if(indexTauMax!=-1) {

	pfTau_=&taus_pf_p4().at(indexTauMax);
	for(int ipf=0; ipf<taus_pf_pfcandIndicies().at(indexTauMax).size(); ipf ++) {	  

	  int index=(taus_pf_pfcandIndicies().at(indexTauMax)).at(ipf);
	  //	  cout << "part  " << ipf << " with pt " << pfcands_p4().at(index).pt() << " eta " << pfcands_p4().at(index).eta()  << endl;

	}	

	if(taus_pf_pfcandIndicies().at(indexTauMax).size()>0) {
	  int leadingPtCand_index=(taus_pf_pfcandIndicies().at(indexTauMax)).at(0);
	  pfTau_leadPtcand_= &(pfcands_p4().at(leadingPtCand_index));
	  pfTau_leadPtcandID_= pfcands_particleId().at(leadingPtCand_index);

	}
      }


      if(indexTauLooseMax!=-1) {

	pfTauLoose_=&taus_pf_p4().at(indexTauLooseMax);
	for(int ipf=0; ipf<taus_pf_pfcandIndicies().at(indexTauLooseMax).size(); ipf ++) {	  

	  int index=(taus_pf_pfcandIndicies().at(indexTauLooseMax)).at(ipf);
	  //	  cout << "part  " << ipf << " with pt " << pfcands_p4().at(index).pt() << " eta " << pfcands_p4().at(index).eta()  << endl;

	}	

	if(taus_pf_pfcandIndicies().at(indexTauLooseMax).size()>0) {
	  int leadingPtCand_index=(taus_pf_pfcandIndicies().at(indexTauLooseMax)).at(0);
	  pfTauLoose_leadPtcand_= &(pfcands_p4().at(leadingPtCand_index));
	  pfTauLoose_leadPtcandID_= pfcands_particleId().at(leadingPtCand_index);

	}
      }

      //----------------------------------------
      // nvertex variables
      //----------------------------------------

      nvtx_ = 0;
    
      for (size_t v = 0; v < vtxs_position().size(); ++v){
	if(isGoodVertex(v)) ++nvtx_;
      }

      indexfirstGoodVertex_=firstGoodVertex();

      ntruepu_ = 0;
      npu_ = 0;
      npuMinusOne_ = 0;
      npuPlusOne_ = 0;
      if ( !isData ) {
        //Information for out-of-time PU
        for (unsigned int nbc=0;nbc<puInfo_nPUvertices().size();++nbc) {
	  if (puInfo_bunchCrossing().at(nbc)==0) {
	    npu_ = puInfo_nPUvertices().at(nbc);
	    ntruepu_ = (int)puInfo_trueNumInteractions().at(nbc);
	  } else if (puInfo_bunchCrossing().at(nbc)==-1) npuMinusOne_ = puInfo_nPUvertices().at(nbc);
	  else if (puInfo_bunchCrossing().at(nbc)==+1) npuPlusOne_ = puInfo_nPUvertices().at(nbc);
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
	pdfid1_ = int(pdfinfo_id1());
	pdfid2_ = int(pdfinfo_id2()); 
	pdfQ_   = pdfinfo_scale();
	pdfx1_  = pdfinfo_x1();
	pdfx2_  = pdfinfo_x2();
      }
      //-------------------------------------
      // jet counting
      //-------------------------------------

      VofP4 mediumBJets;

      int   imaxjet   = -1;
      float maxjetpt  = -1.;
      float ht_ = 0.;

      VofP4 vpfrawjets_p4;
      vpfrawjets_p4.clear();
      VofiP4 vipfjets_p4;
      vipfjets_p4.clear();
      //corrections for all 
      vector<float> l1cors_all;
      l1cors_all.clear();

      vector<float> fullcors;
      vector<float> l2l3cors;
      vector<float> rescors;
      vector<float> l1cors;
      fullcors.clear();
      l2l3cors.clear();
      rescors.clear();
      l1cors.clear();

      rhovor_ = evt_ww_rho_vor();

      float dmetx  = 0.0;
      float dmety  = 0.0;
      float jetptx = 0.0;
      float jetpty = 0.0;

      for (unsigned int ijet = 0 ; ijet < pfjets_p4().size() ; ijet++) {
	
	// skip jets with |eta| > 5.0
	if( fabs( pfjets_p4().at(ijet).eta() ) > 5.0 ) {
	  l1cors_all.push_back( -99999. );
	  continue;
	}

	// get L1FastL2L3Residual total correction
	jet_corrector_pfL1FastJetL2L3->setRho   ( evt_ww_rho_vor()           );
	jet_corrector_pfL1FastJetL2L3->setJetA  ( pfjets_area().at(ijet)     );
	jet_corrector_pfL1FastJetL2L3->setJetPt ( pfjets_p4().at(ijet).pt()  );
	jet_corrector_pfL1FastJetL2L3->setJetEta( pfjets_p4().at(ijet).eta() );
	double corr = jet_corrector_pfL1FastJetL2L3->getCorrection();

	// get L1Fast, L2, L3, Residual individual corrections
	jet_corrector_pfL1FastJetL2L3->setRho   ( evt_ww_rho_vor()           );
	jet_corrector_pfL1FastJetL2L3->setJetA  ( pfjets_area().at(ijet)     );
	jet_corrector_pfL1FastJetL2L3->setJetPt ( pfjets_p4().at(ijet).pt()  );
	jet_corrector_pfL1FastJetL2L3->setJetEta( pfjets_p4().at(ijet).eta() );
	vector<float> factors = jet_corrector_pfL1FastJetL2L3->getSubCorrections();

	// get residual correction only
	float rescorr = 1;
	if( isData ){
	  if( factors.size() == 4 ) rescorr = factors.at(3) / factors.at(2);
	  else                      cout << "ERROR! " << factors.size() << " jetSubCorrections" << endl;
	}

	l1cors_all.push_back( factors.at(0) );

	LorentzVector vjet      = corr    * pfjets_p4().at(ijet);
	indP4 ivjet = { vjet, ijet };

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

	//---------------------------------------------------------------------------                                                                                 
        // Matching leptons with pfjet 
        //---------------------------------------------------------------------------                                                                                 
	float dr1 = ROOT::Math::VectorUtil::DeltaR( vjet, *lep1_ );

	if( dr1 < lep1pfjetdr_ ){
	  lep1pfjetdr_    = dr1;
	  minvjet1         = vjet;
	}

	if(lep2_) {
	  float dr2 = ROOT::Math::VectorUtil::DeltaR( vjet, *lep2_ );
	  
	  if( dr2 < lep2pfjetdr_ ){
	    lep2pfjetdr_    = dr2;
	    minvjet2        = vjet;
	  }
	}

	// lepton-jet overlap removal
	bool rejectJet = false;
	for( int ilep = 0 ; ilep < (int)goodLeptons.size() ; ilep++ ){
	  if( ROOT::Math::VectorUtil::DeltaR( vjet , goodLeptons.at(ilep) ) < 0.4 ) rejectJet = true;  
	}
	if( rejectJet ) continue;
	          
	// PFJetID
	if( !passesPFJetID(ijet) ){
	  jetid_ = 0;
	  if( vjet.pt() > 30 && fabs( vjet.eta() ) < 2.4 ) jetid30_ = 0;
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
	if( vjet.pt() > 20 && fabs( vjet.eta() ) < 4.7 ){
	  vipfjets_p4.push_back( ivjet );
	}

	// njets JEC up
	if( vjetUp.pt() > 30. && fabs( vjetUp.eta() ) < 2.4 ){
	  njetsUp_++;
	  htUp_ += vjetUp.pt();
	}

	// njets JEC down
	if( vjetDown.pt() > 30. && fabs( vjetDown.eta() ) < 2.4 ){
	  njetsDown_++;
	  htDown_ += vjetDown.pt();
	}

	// njets: L1FastL2L3Residual, pt > 30 GeV
	if(       vjet.pt()    < 30. )           continue;
	if( fabs( vjet.eta() ) > 2.4 )           continue;

	htpf30_ += vjet.pt();
        npfjets30_ ++;
       
        //count jets that are not overlapping with second lepton                                                                                             
        npfjets30lepcorr_ ++;
	if (nleps_==2 &&  mclep2_->Pt() > 30.
	    && ROOT::Math::VectorUtil::DeltaR(*mclep2_, vjet) < 0.4 )
	npfjets30lepcorr_ --;
     

	//-------------------------------------
	// b-tag counting
	//-------------------------------------
	
	// btag variables: CSVM
	float discrimcsv = pfjets_combinedSecondaryVertexBJetTag().at(ijet);
	bool isbtagcsvm = ( discrimcsv > 0.679 ) ? true: false;
	if (isbtagcsvm) {  
	  nbtagscsvm_++;
	  mediumBJets.push_back(vjet);
	}

	// btag variables: CSVL
	bool isbtagcsvl = ( discrimcsv > 0.244 ) ? true: false;
	if (isbtagcsvl)     nbtagscsvl_++;

	// btag variables: CSVT
	bool isbtagcsvt = ( discrimcsv > 0.898 ) ? true: false;
	if (isbtagcsvt)     nbtagscsvt_++;

	// btag variables: SSV -- not supported anymore
	float discrimssv = pfjets_simpleSecondaryVertexHighEffBJetTag().at(ijet);
	bool isbtagssv = ( discrimssv > 1.74 ) ? true : false;
	if (isbtagssv) nbtagsssv_++;

	// btag variables: TCHEL -- not supported anymore
	float discrimtche = pfjets_trackCountingHighEffBJetTag().at(ijet);
	bool isbtagtcl = ( discrimtche > 1.7 ) ? true: false;
	if (isbtagtcl)     nbtagstcl_++;

	// btag variables: TCHEM -- not supported anymore
	bool isbtagtcm = ( discrimtche > 3.3 ) ? true: false;
	if (isbtagtcm)     nbtagstcm_++;

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


      // store the lepton jet matching
      
      if (lep2pfjetdr_<9998.) leppfjet2_      = &minvjet2;      
      if (lep1pfjetdr_<9998.) leppfjet1_ = &minvjet1;

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
	getPhiCorrMET( t1met10_, t1met10phi_, nvtx_, !isData);
      t1metphicorr_    = p_t1metphicorr.first;
      t1metphicorrphi_ = p_t1metphicorr.second;

      //official prescription
      std::pair<float, float> p_t1met_off = 
	met_corrector_pfL1FastJetL2L3->getCorrectedMET();
      t1met_off_	= p_t1met_off.first;
      t1metphi_off_	= p_t1met_off.second;
      t1metmt_off_ =
	getMT( lep1_->pt() , lep1_->phi() , t1met_off_ , t1metphi_off_ );

      std::pair<float, float> p_t1metphicorr_off = 
	getPhiCorrMET( t1met_off_, t1metphi_off_, nvtx_, !isData);	
      
      t1metphicorr_off_		= p_t1metphicorr_off.first;
      t1metphicorrphi_off_     	= p_t1metphicorr_off.second;
      t1metphicorrmt_off_ = 
	getMT( lep1_->pt() , lep1_->phi() , t1metphicorr_off_ , t1metphicorrphi_off_ );

      // MET after Jet PT smearing.
      pair<float, float> p_t1met10Smear = 
	Type1PFMETSmearRec(jetSmearer, isData, vpfrawjets_p4 , fullcors, t1metphicorr_, t1metphicorrphi_);
      t1met10s_ = p_t1met10Smear.first;
      t1met10sphi_ = p_t1met10Smear.second;

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

      //--------------------------------
      // Store vector of jets
      //--------------------------------

      // store L1FastL2L3Residual pfjets
      // check if jet is b-tagged
      //      sort(vpfjets_p4.begin(), vpfjets_p4.end(), sortByPt);
      sort(vipfjets_p4.begin(), vipfjets_p4.end(), sortIP4ByPt);

      for( int i = 0 ; i < (int)vipfjets_p4.size() ; ++i ){
	pfjets_.push_back(vipfjets_p4.at(i).p4obj);
	pfjets_csv_.push_back(pfjets_combinedSecondaryVertexBJetTag().at(vipfjets_p4.at(i).p4ind));
	pfjets_chm_.push_back(pfjets_chargedMultiplicity().at(vipfjets_p4.at(i).p4ind));
	pfjets_neu_.push_back(pfjets_neutralMultiplicity().at(vipfjets_p4.at(i).p4ind));
	pfjets_qgtag_.push_back(QGtagger(vipfjets_p4.at(i).p4obj,vipfjets_p4.at(i).p4ind,qglikeli_));
	pfjets_lrm_.push_back(getLRM(vipfjets_p4.at(i).p4ind,1));
	pfjets_lrm2_.push_back(getLRM(vipfjets_p4.at(i).p4ind,2));

	pfUncertainty->setJetEta(vipfjets_p4.at(i).p4obj.eta());
	pfUncertainty->setJetPt(vipfjets_p4.at(i).p4obj.pt());   // here you must use the CORRECTED jet pt
	double unc = pfUncertainty->getUncertainty(true);

	pfjets_uncertainty_.push_back(unc);

	//	cout << "lrm " << getLRM(vipfjets_p4.at(i).p4ind) << endl;

	pfjets_corr_.push_back(vipfjets_p4.at(i).p4obj.pt()/pfjets_p4().at(vipfjets_p4.at(i).p4ind).pt());
	pfjets_l1corr_.push_back(l1cors_all.at(vipfjets_p4.at(i).p4ind));

        // Save the jet resolution. if it is MC scale by a factor given by
        // https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetResolution
        float sigma =  getJetResolution(vipfjets_p4.at(i).p4obj, jetSmearer );
        if ( !isData ) sigma *= getDataMCRatio(vipfjets_p4.at(i).p4obj.eta());
	pfjets_sigma_.push_back( sigma );

	if( isData ){

	  pfjets_mc3_.push_back( -1 );
	  pfjets_flav_.push_back( -1 );
	  pfjets_mcflavorAlgo_.push_back( -1 );
	  pfjets_mcflavorPhys_.push_back( -1 );
	  pfjets_genJetDr_.push_back( -1.0 );
	  pfjets_lepjet_.push_back( -1 );

	}
	else{

	  pfjets_mc3_.push_back(isGenQGLMatched( vipfjets_p4.at(i).p4obj, 0.4 ));
	  pfjets_mcflavorAlgo_.push_back( pfjets_mcflavorAlgo().at(vipfjets_p4.at(i).p4ind) );
	  pfjets_mcflavorPhys_.push_back( pfjets_mcflavorPhys().at(vipfjets_p4.at(i).p4ind) );
	  pfjets_flav_.push_back((isGenBMatched(vipfjets_p4.at(i).p4obj, 0.5) ? 5 : (isGenCMatched(vipfjets_p4.at(i).p4obj, 0.5) ? 4 : 0)));
	  pfjets_genJetDr_.push_back(dRGenJet ( vipfjets_p4.at(i).p4obj ));
	  pfjets_lepjet_.push_back(getLeptonMatchIndex( &vipfjets_p4.at(i).p4obj,mclep1_, mclep2_, 0.4 ));
	  if(genjets_p4().size()>indexGenJet(vipfjets_p4.at(i).p4obj)) pfjets_genJet_.push_back(genjets_p4().at(indexGenJet(vipfjets_p4.at(i).p4obj)));

	} 

	//beta variables
	pfjets_beta_.push_back(pfjet_beta(vipfjets_p4.at(i).p4ind,1));
	pfjets_beta2_.push_back(pfjet_beta(vipfjets_p4.at(i).p4ind,2));
        pfjets_beta_0p1_.push_back(  pfjet_beta( vipfjets_p4.at(i).p4ind, 1, 0.1 ) );
        pfjets_beta_0p2_.push_back(  pfjet_beta( vipfjets_p4.at(i).p4ind, 1, 0.2 ) );
        pfjets_beta2_0p1_.push_back(  pfjet_beta( vipfjets_p4.at(i).p4ind, 2, 0.1 ) );
	pfjets_beta2_0p5_.push_back(  pfjet_beta( vipfjets_p4.at(i).p4ind, 2, 0.5 ) );
        // pfjets_beta_0p15_.push_back( pfjet_beta( vipfjets_p4.at(i).p4ind, 1, 0.15) );
        // pfjets_beta2_0p15_.push_back( pfjet_beta( vipfjets_p4.at(i).p4ind, 2, 0.15) );
        // pfjets_beta_0p2_.push_back(  pfjet_beta( vipfjets_p4.at(i).p4ind, 1, 0.2 ) );
        // pfjets_beta2_0p2_.push_back(  pfjet_beta( vipfjets_p4.at(i).p4ind, 2, 0.2 ) );	
	pfjets_mvaPUid_.push_back(pfjets_full53xmvavalue().at(vipfjets_p4.at(i).p4ind));
	pfjets_mva5xPUid_.push_back(pfjets_full5xmvavalue().at(vipfjets_p4.at(i).p4ind));
	pfjets_mvaBeta_.push_back(pfjets_full53xmva_beta().at(vipfjets_p4.at(i).p4ind));

      }
      l1cors_all.clear();

      //store distance to closest jet for pfcand
      if ( nleps_ == 2 ) {
	pflepmindrj_  = getminjdr( vipfjets_p4, pflep_ );
	pftaudmindrj_ = getminjdr( vipfjets_p4, pftaud_ );
      }
      pfcandmindrj5_  = getminjdr( vipfjets_p4, pfcand5_ );
      pfcandmindrj10_ = getminjdr( vipfjets_p4, pfcand10_ );

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

      //--------------------------------
      // Hadronic Top and MT2
      //--------------------------------

      vector<LorentzVector> jets;
      vector<float> sigma_jets;

      int n_jets = pfjets_.size();
      double sigma[n_jets];

      for( unsigned int i = 0 ; i < n_jets ; ++i ){

	jets.push_back( pfjets_.at(i));

	sigma[i] = pfjets_sigma_.at(i);
        if ( !isData ) sigma[i] *= getDataMCRatio(jets[i].eta());
	sigma_jets.push_back(sigma[i]);

      } 

      // get list of candidates
      PartonCombinatorics pc (jets , pfjets_csv_, sigma_jets, pfjets_mc3_, *lep1_, t1metphicorr_, t1metphicorrphi_, isData);
      MT2CHI2 mc = pc.getMt2Chi2();

      // chi2 and MT2 variables
      chi2min_= mc.one_chi2;               // minimum chi2 
      chi2minprob_= TMath::Prob(chi2min_,1);   // probability of minimum chi2

      mt2bmin_= mc.three_mt2b;             // minimum MT2b
      mt2blmin_= mc.three_mt2bl;            // minimum MT2bl
      mt2wmin_= mc.three_mt2w;             // minimum MT2w

      //--------------------------------
      // get non-isolated leptons
      //--------------------------------

      nonisoel_ = 0;
      nonisomu_ = 0;

      float maxelpt  = -1.;
      int   imaxelpt = -1;

      for( unsigned int iel = 0 ; iel < els_p4().size(); ++iel ){

	if( els_p4().at(iel).pt() < 10 )                                                                         continue;
	if( !passElectronSelection_Stop2012_v3_NoIso( iel , vetoTransition,vetoTransition,useOldIsolation))      continue;

	// don't count the leptons that we already counted as good
	bool isGoodLepton = false;
	for( int ilep = 0 ; ilep < (int)goodLeptons.size() ; ilep++ ){
	  if( ROOT::Math::VectorUtil::DeltaR( els_p4().at(iel) , goodLeptons.at(ilep) ) < 0.1 ) isGoodLepton = true;  
	}
	if( isGoodLepton ) continue;

	// don't count leptons near b-jets (SSVM)
	bool nearBJet = false;
	for( int ijet = 0 ; ijet < (int)mediumBJets.size() ; ijet++ ){
	  if( ROOT::Math::VectorUtil::DeltaR( els_p4().at(iel) , mediumBJets.at(ijet) ) < 0.4 ) nearBJet = true;
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
	  if( ROOT::Math::VectorUtil::DeltaR( mus_p4().at(imu) , goodLeptons.at(ilep) ) < 0.1 ) isGoodLepton = true;  
	}
	if( isGoodLepton ) continue;

	// don't count leptons near b-jets SSVM)
	bool nearBJet = false;
	for( int ijet = 0 ; ijet < (int)mediumBJets.size() ; ijet++ ){
	  if( ROOT::Math::VectorUtil::DeltaR( mus_p4().at(imu) , mediumBJets.at(ijet) ) < 0.4 ) nearBJet = true;
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
      // gen jets
      //---------------------------------

      ngenjets_ = 0;
      htgen_    = 0;

      if( !isData ){
	for (unsigned int igjet = 0 ; igjet < genjets_p4().size() ; igjet++) {
	    
	  LorentzVector vgjet = genjets_p4().at(igjet);

	  bool rejectJet = false;
	  for( int ilep = 0 ; ilep < (int)goodLeptons.size() ; ilep++ ){
	    if( ROOT::Math::VectorUtil::DeltaR( vgjet , goodLeptons.at(ilep) ) < 0.4 ) rejectJet = true;  
	  }
	  if( rejectJet ) continue;
	    
	  if( vgjet.pt() < 30.                   )  continue;
	  if( fabs( vgjet.eta() ) > 2.4          )  continue;
	    
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

      // calo met
      calomet_ = evt_metMuonCorr();
      calometphi_ = evt_metMuonCorrPhi() ;

      // track met
      pair<float, float> trkMET = getTrackerMET(lep1_); 
      trkmet_=trkMET.first;
      trkmetphi_=trkMET.second;

      pair<float, float> trkMET_nolepcorr = getTrackerMET(lep1_, 0.1, false); 
      trkmet_nolepcorr_=trkMET_nolepcorr.first;
      trkmetphi_nolepcorr_=trkMET_nolepcorr.second;

      //---------------------------
      // set event weight
      //---------------------------

      weight_ = -1.;

      if( TString(prefix).Contains("T2") ){
	mG_ = -999; //sparm_mG();
        mL_ = -999; //sparm_mL();
        x_  = -999; //sparm_mf();

	if( TString(prefix).Contains("T2tt") ){
	  for (int i=0; i<(int)sparm_values().size(); ++i) {
	    if (sparm_names().at(i).Contains("mstop")) mG_ = sparm_values().at(i);
	    if (sparm_names().at(i).Contains("mlsp")) mL_ = sparm_values().at(i);
	  }
	}
	
	if( TString(prefix).Contains("T2bw") ){
	  for (int i=0; i<(int)sparm_values().size(); ++i) {
	    if (sparm_names().at(i).Contains("x")) x_ = sparm_values().at(i);
	    if (sparm_names().at(i).Contains("mstop")) mG_ = sparm_values().at(i);
	    if (sparm_names().at(i).Contains("mlsp")) mL_ = sparm_values().at(i);
	  }
        }
	
        xsecsusy_  = mG_ > 0. ? stopPairCrossSection(mG_) : -999;
        weight_ = xsecsusy_ > 0. ? lumi * xsecsusy_ * (1000./50000.) : -999.;

	if( doTenPercent )	  weight_ *= 10;
      }

      else if(TString(prefix).Contains("TChiwh")) {
        // need to put some reasonable numbers here..
	// currently applying xsecs downstream
	weight_ = 1.;
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

      t1met10smt_   = getMT( lep1_->pt() , lep1_->phi() , t1met10s_     , t1met10sphi_ );

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

      strcpy(dataset_, evt_dataset().at(0).Data());  //dataset name
      run_          = evt_run();                          //run
      lumi_         = evt_lumiBlock();                    //lumi
      event_        = evt_event();                        //event
      nvtxweight_   = vtxweight(isData);

      csc_       = evt_cscTightHaloId();
      hbhe_      = evt_hbheFilter();
      hcallaser_ = filt_hcalLaser();
      ecaltp_    = filt_ecalTP();
      trkfail_   = filt_trackingFailure();
      eebadsc_   = 1;
      if( isData ) eebadsc_ = filt_eeBadSc();
      hbhenew_   = passHBHEFilter();

      //for SCs above 20 GeV, store the max. ECAL laser correction and the eta/pT of this SC
      if( !TString(prefix).Contains("T2") ) {
	scslasercormax_ = -99999.;
	for (int isc= 0; isc<scs_p4().size(); ++isc) {
	  if (scs_p4()[isc].pt()<20.) continue;
	  if ( scs_laserCorMax()[isc] > scslasercormax_ ) {
	    scslasercormax_ = scs_laserCorMax()[isc];
	    scslasercormax_pt_ = scs_p4()[isc].pt();
	    scslasercormax_eta_ = scs_eta()[isc];
	  }
	}
      }

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
      lep1_badecallaser_ = 0;
      lep2_badecallaser_ = 0;

      if( leptype_ == 0 ){
	iso1_   = electronIsolation_rel   ( index1 , true ); //truncated
	isont1_ = electronIsolation_rel_v1( index1 , true ); //non-truncated
	isopfold1_ = electronIsoValuePF2012_FastJetEffArea( index1 , 0.3 , 0 );
	isopf1_ = electronIsoValuePF2012_FastJetEffArea_v3( index1 , 0.3 , 0 , false);
	etasc1_ = els_etaSC()[index1];
	eoverpin_  = els_eOverPIn ()[index1];
	eoverpout_ = els_eOverPOut()[index1];
	dEtaIn_ = els_dEtaIn()[index1];
	dPhiIn_ = els_dPhiIn()[index1];
	sigmaIEtaIEta_ = els_sigmaIEtaIEta()[index1];
	hOverE_ = els_hOverE()[index1];
	ooemoop_ = fabs( (1.0/els_ecalEnergy()[index1]) - (els_eOverPIn()[index1]/els_ecalEnergy()[index1]) );
	d0vtx_ = electron_d0PV_smurfV3(index1);
	dzvtx_ = electron_dzPV_smurfV3(index1);
	expinnerlayers_ = els_exp_innerlayers()[index1];
	fbrem_ = els_fbrem()[index1];
	pfisoch_ = els_iso03_pf2012ext_ch().at(index1);
	pfisoem_ = els_iso03_pf2012ext_em().at(index1);
	pfisonh_ = els_iso03_pf2012ext_nh().at(index1);
	eSC_ = els_eSC()[index1];
	phiSC_ = els_phiSC()[index1];
	eSCRaw_ = els_eSCRaw()[index1];
	eSCPresh_ = els_eSCPresh()[index1];
	if( !TString(prefix).Contains("T2") ) {
	  //      lep1_scslasercormean_ = scs_laserCorMean()[index1];
	  //      lep1_scslasercormax_  = scs_laserCorMax()[index1];
	  lep1_scslasercormean_ = scs_laserCorMean()[els_scindex()[index1]];
	  lep1_scslasercormax_  = scs_laserCorMax()[els_scindex()[index1]];
	  //      cout<<scs_eta()[els_scindex()[index1]]<<" "<<lep1_scslasercormax_<<endl;
	  if ((fabs(scs_eta()[els_scindex()[index1]])<1.4442)&&(lep1_scslasercormax_>2.0) ||
	      (fabs(scs_eta()[els_scindex()[index1]])>1.4442)&&(lep1_scslasercormax_>3.0))
	    lep1_badecallaser_ = 1;
	}
      }
      else if( leptype_ == 1 ){
	iso1_   = muonIsoValue( index1 , true  ); //truncated 
	isont1_ = muonIsoValue( index1 , false ); //non-truncated
	isopf1_ = muonIsoValuePF2012_deltaBeta( index1 );
	ecalveto1_ = mus_iso_ecalvetoDep().at(index1);
	hcalveto1_ = mus_iso_hcalvetoDep().at(index1);
      }
      
      //fill variables for second lepton if it exists
      if(abs(id2_) == 11) {
	iso2_   = electronIsolation_rel   ( index2 , true ); //truncated
	isont2_ = electronIsolation_rel_v1( index2 , true ); //non-truncated
	isopf2_ = electronIsoValuePF2012_FastJetEffArea_v2( index2 , 0.3 , 0 , false);
	etasc2_ = els_etaSC()[index2];
	eoverpin2_  = els_eOverPIn ()[index2];
	eoverpout2_ = els_eOverPOut()[index2];
	dEtaIn2_ = els_dEtaIn()[index2];
	dPhiIn2_ = els_dPhiIn()[index2];
	sigmaIEtaIEta2_ = els_sigmaIEtaIEta()[index2];
	hOverE2_ = els_hOverE()[index2];
	ooemoop2_ = fabs( (1.0/els_ecalEnergy()[index2]) - (els_eOverPIn()[index2]/els_ecalEnergy()[index2]) );
	d0vtx2_ = electron_d0PV_smurfV3(index2);
	dzvtx2_ = electron_dzPV_smurfV3(index2);
	expinnerlayers2_ = els_exp_innerlayers()[index2];
	fbrem2_ = els_fbrem()[index2];
	pfisoch2_ = els_iso03_pf2012ext_ch().at(index2);
	pfisoem2_ = els_iso03_pf2012ext_em().at(index2);
	pfisonh2_ = els_iso03_pf2012ext_nh().at(index2);
	eSC2_ = els_eSC()[index2];
	phiSC2_ = els_phiSC()[index2];
	eSCRaw2_ = els_eSCRaw()[index2];
	eSCPresh2_ = els_eSCPresh()[index2];
	if( !TString(prefix).Contains("T2") ) {
	  lep2_scslasercormean_ = scs_laserCorMean()[els_scindex()[index2]];
	  lep2_scslasercormax_  = scs_laserCorMax()[els_scindex()[index2]];
	  if ((fabs(scs_eta()[els_scindex()[index2]])<1.4442)&&(lep2_scslasercormax_>2.0) ||
	      (fabs(scs_eta()[els_scindex()[index2]])>1.4442)&&(lep2_scslasercormax_>3.0))
	    lep2_badecallaser_ = 1;
	}
      }
      else if(abs(id2_) == 13) {
	iso2_   = muonIsoValue( index2 , true  ); //truncated 
	isont2_ = muonIsoValue( index2 , false ); //non-truncated
	isopf2_ = muonIsoValuePF2012_deltaBeta( index2 );
	
	ecalveto2_ = mus_iso_ecalvetoDep().at(index2);
	hcalveto2_ = mus_iso_hcalvetoDep().at(index2);
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
      sltrigweight_ = isData ? 1. : getsltrigweight( id1_, lep1_->pt() , lep1_->eta() );
      dltrigweight_ = (!isData && ngoodlep_>1) ? getdltrigweight( id1_, id2_ ) : 1.;
    
      /// hadronic stop reconstruction
      //      candidates_.clear(); 

      //      for (list<Candidate>::iterator it = candidates.begin(); it != candidates.end(); ++it)
      //          candidates_.push_back(*it);

      jets_.clear();
      btag_.clear();
      for (unsigned int i = 0; i < vpfrawjets_p4.size(); i++){
         jets_.push_back(vpfrawjets_p4.at(i));
         btag_.push_back( pfjets_combinedSecondaryVertexBJetTag().at(i) );
      }


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
  //  outFile   = new TFile("baby.root","RECREATE");
  outFile->cd();
  outTree = new TTree("t","Tree");

  //Set branch addresses
  //variables must be declared in singleLeptonLooper.h
  outTree->Branch("acc_2010",        &acc_2010_,         "acc_2010/I");
  outTree->Branch("acc_highmet",     &acc_highmet_,      "acc_highmet/I");
  outTree->Branch("acc_highht",      &acc_highht_,       "acc_highht/I");

  outTree->Branch("eldup"     ,  &eldup_     ,  "eldup/I");  
  outTree->Branch("csc"       ,  &csc_       ,  "csc/I");  
  outTree->Branch("hbhe"      ,  &hbhe_      ,  "hbhe/I");  
  outTree->Branch("hbhenew"   ,  &hbhenew_   ,  "hbhenew/I");  
  outTree->Branch("hcallaser" ,  &hcallaser_ ,  "hcallaser/I");  
  outTree->Branch("ecaltp"    ,  &ecaltp_    ,  "ecaltp/I");  
  outTree->Branch("trkfail"   ,  &trkfail_   ,  "trkfail/I");  
  outTree->Branch("eebadsc"   ,  &eebadsc_   ,  "eebadsc/I");  
  outTree->Branch("lep1_badecallaser" ,  &lep1_badecallaser_ ,  "lep1_badecallaser/I");  
  outTree->Branch("lep2_badecallaser" ,  &lep2_badecallaser_ ,  "lep2_badecallaser/I");  

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
  outTree->Branch("nwzpartons",      &nwzpartons_,       "nwzpartons/I");
  outTree->Branch("hyptype",         &hyptype_,          "hyptype/I");
  outTree->Branch("maxpartonpt",     &maxpartonpt_,      "maxpartonpt/F");
  outTree->Branch("etattbar",        &etattbar_,         "etatbar/F");
  outTree->Branch("njetsoffset",     &njetsoffset_,      "njetsoffset/I");
  outTree->Branch("njetsuncor",      &njetsuncor_,       "njetsuncor/I");
  outTree->Branch("costhetaweight",  &costhetaweight_,   "costhetaweight/F");
  outTree->Branch("weight",          &weight_,           "weight/F");
  outTree->Branch("weightleft",      &weightleft_,       "weightleft/F");
  outTree->Branch("weightright",     &weightright_,      "weightright/F");
  outTree->Branch("mutrigweight",    &mutrigweight_,     "mutrigweight/F");
  outTree->Branch("mutrigweight2",   &mutrigweight2_,    "mutrigweight2/F");
  outTree->Branch("sltrigweight",    &sltrigweight_,     "sltrigweight/F");
  outTree->Branch("dltrigweight",    &dltrigweight_,     "dltrigweight/F");
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
  outTree->Branch("calomet",         &calomet_,          "calomet/F");
  outTree->Branch("calometphi",      &calometphi_,       "calometphi/F");
  outTree->Branch("trkmet",          &trkmet_,           "trkmet/F");
  outTree->Branch("trkmetphi",       &trkmetphi_,        "trkmetphi/F");
  outTree->Branch("trkmet_nolepcorr",    &trkmet_nolepcorr_,    "trkmet_nolepcorr/F");
  outTree->Branch("trkmetphi_nolepcorr", &trkmetphi_nolepcorr_, "trkmetphi_nolepcorr/F");
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
  outTree->Branch("npfjets30lepcorr", &npfjets30lepcorr_, "npfjets30lepcorr/I");
  outTree->Branch("knjets",           &knjets_,           "knjets/F");

  //rho correction
  outTree->Branch("rhovor",          &rhovor_,           "rhovor/F");
  outTree->Branch("htpf30",          &htpf30_,           "htpf30/F");

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

  outTree->Branch("t1met10s",         &t1met10s_,          "t1met10s/F");
  outTree->Branch("t1met10sphi",      &t1met10sphi_,       "t1met10sphi/F");
  outTree->Branch("t1met10smt",       &t1met10smt_,       "t1met10smt/F");

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

  //official prescription
  outTree->Branch("t1met_off"           , &t1met_off_           , "t1met_off/F");
  outTree->Branch("t1metphi_off"        , &t1metphi_off_        , "t1metphi_off/F");
  outTree->Branch("t1metmt_off"         , &t1metmt_off_         , "t1metmt_off/F");
  outTree->Branch("t1metphicorr_off"    , &t1metphicorr_off_    , "t1metphicorr_off/F");
  outTree->Branch("t1metphicorrphi_off" , &t1metphicorrphi_off_ , "t1metphicorrphi_off/F");
  outTree->Branch("t1metphicorrmt_off"  , &t1metphicorrmt_off_  , "t1metphicorrmt_off/F");

  // MT2 and CHI2
  outTree->Branch("mt2bmin"           , &mt2bmin_           , "mt2bmin/F");
  outTree->Branch("mt2blmin"          , &mt2blmin_          , "mt2blmin/F");
  outTree->Branch("mt2wmin"           , &mt2wmin_           , "mt2wmin/F");
  outTree->Branch("chi2min"           , &chi2min_           , "chi2min/F");
  outTree->Branch("chi2minprob"       , &chi2minprob_       , "chi2minprob/F");

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
  outTree->Branch("njetsUp",          &njetsUp_,          "njetsUp/I");
  outTree->Branch("njetsDown",        &njetsDown_,        "njetsDown/I");
  outTree->Branch("htUp",             &htUp_,             "htUp/F");
  outTree->Branch("htDown",           &htDown_,           "htDown/F");
  outTree->Branch("ntruepu",          &ntruepu_,          "ntruepu/I");
  outTree->Branch("npu",              &npu_,              "npu/I");
  outTree->Branch("npuMinusOne",      &npuMinusOne_,      "npuMinusOne/I");
  outTree->Branch("npuPlusOne",       &npuPlusOne_,       "npuPlusOne/I");
  outTree->Branch("nvtx",             &nvtx_,             "nvtx/I");
  outTree->Branch("indexfirstGoodVertex_",             &indexfirstGoodVertex_,             "indexfirstGoodVertex/I");
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
  outTree->Branch("isopfold1",    &isopfold1_,     "isopfold1/F");
  outTree->Branch("isopf1",           &isopf1_,           "isopf1/F");
  outTree->Branch("etasc1",           &etasc1_,           "etasc1/F");
  outTree->Branch("etasc2",           &etasc2_,           "etasc2/F");
  outTree->Branch("eoverpin",         &eoverpin_,         "eoverpin/F");
  outTree->Branch("eoverpout",        &eoverpout_,        "eoverpout/F");
  outTree->Branch("dEtaIn", &dEtaIn_, "dEtaIn/F");
  outTree->Branch("dPhiIn", &dPhiIn_, "dPhiIn/F");
  outTree->Branch("sigmaIEtaIEta", &sigmaIEtaIEta_, "sigmaIEtaIEta/F");
  outTree->Branch("hOverE", &hOverE_, "hOverE/F");
  outTree->Branch("ooemoop", &ooemoop_, "ooemoop/F");
  outTree->Branch("d0vtx", &d0vtx_, "d0vtx/F");
  outTree->Branch("dzvtx", &dzvtx_, "dzvtx/F");
  outTree->Branch("expinnerlayers", &expinnerlayers_, "expinnerlayers/F");
  outTree->Branch("fbrem", &fbrem_, "fbrem/F");
  outTree->Branch("pfisoch", &pfisoch_, "pfisoch/F");
  outTree->Branch("pfisoem", &pfisoem_, "pfisoem/F");
  outTree->Branch("pfisonh", &pfisonh_, "pfisonh/F");
  outTree->Branch("eSC", & eSC_, "eSC/F");
  outTree->Branch("phiSC", & phiSC_, "phiSC/F");
  outTree->Branch("eSCRaw", & eSCRaw_, "eSCRaw/F");
  outTree->Branch("eSCPresh", & eSCPresh_, "eSCPresh/F");
  outTree->Branch("lep1_scslasercormean", &lep1_scslasercormean_, "lep1_scslasercormean/F");
  outTree->Branch("lep1_scslasercormax", &lep1_scslasercormax_, "lep1_scslasercormax/F");

  outTree->Branch("eoverpin2",         &eoverpin2_,         "eoverpin2/F");
  outTree->Branch("eoverpout2",        &eoverpout2_,        "eoverpout2/F");
  outTree->Branch("dEtaIn2", &dEtaIn2_, "dEtaIn2/F");
  outTree->Branch("dPhiIn2", &dPhiIn2_, "dPhiIn2/F");
  outTree->Branch("sigmaIEtaIEta2", &sigmaIEtaIEta2_, "sigmaIEtaIEta2/F");
  outTree->Branch("hOverE2", &hOverE2_, "hOverE2/F");
  outTree->Branch("ooemoop2", &ooemoop2_, "ooemoop2/F");
  outTree->Branch("d0vtx2", &d0vtx2_, "d0vtx2/F");
  outTree->Branch("dzvtx2", &dzvtx2_, "dzvtx2/F");
  outTree->Branch("expinnerlayers2", &expinnerlayers2_, "expinnerlayers2/F");
  outTree->Branch("fbrem2", &fbrem2_, "fbrem2/F");
  outTree->Branch("pfisoch2", &pfisoch2_, "pfisoch2/F");
  outTree->Branch("pfisoem2", &pfisoem2_, "pfisoem2/F");
  outTree->Branch("pfisonh2", &pfisonh2_, "pfisonh2/F");
  outTree->Branch("eSC2", & eSC2_, "eSC2/F");
  outTree->Branch("phiSC2", & phiSC2_, "phiSC2/F");
  outTree->Branch("eSCRaw2", & eSCRaw2_, "eSCRaw2/F");
  outTree->Branch("eSCPresh2", & eSCPresh2_, "eSCPresh2/F");
  outTree->Branch("lep2_scslasercormean", &lep2_scslasercormean_, "lep2_scslasercormean/F");
  outTree->Branch("lep2_scslasercormax", &lep2_scslasercormax_, "lep2_scslasercormax/F");
  outTree->Branch("scslasercormax", &scslasercormax_, "scslasercormax/F");
  outTree->Branch("scslasercormax_pt", &scslasercormax_pt_, "scslasercormax_pt/F");
  outTree->Branch("scslasercormax_eta", &scslasercormax_eta_, "scslasercormax_eta/F");

  outTree->Branch("iso2",             &iso2_,             "iso2/F");
  outTree->Branch("ecalveto1",        &ecalveto1_,        "ecalveto1/F");
  outTree->Branch("ecalveto2",        &ecalveto2_,        "ecalveto2/F");
  outTree->Branch("hcalveto1",        &hcalveto1_,        "hcalveto1/F");
  outTree->Branch("hcalveto2",        &hcalveto2_,        "hcalveto2/F");
  outTree->Branch("isont2",           &isont2_,           "isont2/F");
  outTree->Branch("isopf2",           &isopf2_,           "isopf2/F");
  outTree->Branch("ptl1",             &ptl1_,             "ptl1/F");
  outTree->Branch("ptl2",             &ptl2_,             "ptl2/F");
  outTree->Branch("etal1",            &etal1_,            "etal1/F");
  outTree->Branch("etal2",            &etal2_,            "etal2/F");
  outTree->Branch("phil1",            &phil1_,            "phil1/F");
  outTree->Branch("phil2",            &phil2_,            "phil2/F");
  outTree->Branch("meff",             &meff_,             "meff/F");
  outTree->Branch("mt",               &mt_,               "mt/F");
  outTree->Branch("dataset",          &dataset_,          "dataset[200]/C");
  outTree->Branch("run",              &run_,              "run/i");
  outTree->Branch("lumi",             &lumi_,             "lumi/i");
  outTree->Branch("event",            &event_,            "event/i");
  outTree->Branch("y",                &y_,                "y/F");  
  outTree->Branch("ht",               &ht_,               "ht/F");  
  outTree->Branch("htgen",            &htgen_,            "htgen/F");  
  outTree->Branch("htjpt",            &htjpt_,            "htjpt/F");  
  outTree->Branch("nels",             &nels_,             "nels/I");  
  outTree->Branch("nmus",             &nmus_,             "nmus/I");  
  outTree->Branch("ntaus",            &ntaus_,            "ntaus/I");  
  outTree->Branch("nleps",            &nleps_,            "nleps/I");  
  outTree->Branch("nbs",              &nbs_,              "nbs/I");  
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

  outTree->Branch("pfcandid5",        &pfcandid5_,        "pfcandid5/I");
  outTree->Branch("pfcandiso5",       &pfcandiso5_,       "pfcandiso5/F");  
  outTree->Branch("pfcandpt5",        &pfcandpt5_,        "pfcandpt5/F");  
  outTree->Branch("pfcanddz5",        &pfcanddz5_,        "pfcanddz5/F");  
  outTree->Branch("pfcandmindrj5",    &pfcandmindrj5_,    "pfcandmindrj5/F");  

  outTree->Branch("pfcandid10",       &pfcandid10_,       "pfcandid10/I");
  outTree->Branch("pfcandiso10",      &pfcandiso10_,      "pfcandiso10/F");  
  outTree->Branch("pfcandpt10",       &pfcandpt10_,       "pfcandpt10/F");  
  outTree->Branch("pfcanddz10",       &pfcanddz10_,       "pfcanddz10/F");  
  outTree->Branch("pfcandmindrj10",   &pfcandmindrj10_,   "pfcandmindrj10/F");  

  outTree->Branch("pfcandidOS10",       &pfcandidOS10_,       "pfcandidOS10/I");
  outTree->Branch("pfcandisoOS10",      &pfcandisoOS10_,      "pfcandisoOS10/F");  
  outTree->Branch("pfcandptOS10",       &pfcandptOS10_,       "pfcandptOS10/F");  
  outTree->Branch("pfcanddzOS10",       &pfcanddzOS10_,       "pfcanddzOS10/F");  

  outTree->Branch("pfcandid5looseZ",        &pfcandid5looseZ_,        "pfcandid5looseZ/I");
  outTree->Branch("pfcandiso5looseZ",       &pfcandiso5looseZ_,       "pfcandiso5looseZ/F");  
  outTree->Branch("pfcandpt5looseZ",        &pfcandpt5looseZ_,        "pfcandpt5looseZ/F");  
  outTree->Branch("pfcanddz5looseZ",        &pfcanddz5looseZ_,        "pfcanddz5looseZ/F");  

  outTree->Branch("pfcandidOS10looseZ",       &pfcandidOS10looseZ_,       "pfcandidOS10looseZ/I");
  outTree->Branch("pfcandisoOS10looseZ",      &pfcandisoOS10looseZ_,      "pfcandisoOS10looseZ/F");  
  outTree->Branch("pfcandptOS10looseZ",       &pfcandptOS10looseZ_,       "pfcandptOS10looseZ/F");  
  outTree->Branch("pfcanddzOS10looseZ",       &pfcanddzOS10looseZ_,       "pfcanddzOS10looseZ/F");  

  outTree->Branch("pfcanddirid10",       &pfcanddirid10_,       "pfcanddirid10/I");
  outTree->Branch("pfcanddiriso10",      &pfcanddiriso10_,      "pfcanddiriso10/F");  
  outTree->Branch("pfcanddirpt10",       &pfcanddirpt10_,       "pfcanddirpt10/F");  
  outTree->Branch("pfcanddirmindrj10",   &pfcanddirmindrj10_,   "pfcanddirmindrj10/F");  

  outTree->Branch("pfcandvetoid10",       &pfcandvetoid10_,       "pfcandvetoid10/I");
  outTree->Branch("pfcandvetoiso10",      &pfcandvetoiso10_,      "pfcandvetoiso10/F");  
  outTree->Branch("pfcandvetopt10",       &pfcandvetopt10_,       "pfcandvetopt10/F");  
  outTree->Branch("pfcandvetomindrj10",   &pfcandvetomindrj10_,   "pfcandvetomindrj10/F");  

  outTree->Branch("pfcandvetoLid10",       &pfcandvetoLid10_,       "pfcandvetoLid10/I");
  outTree->Branch("pfcandvetoLiso10",      &pfcandvetoLiso10_,      "pfcandvetoLiso10/F");  
  outTree->Branch("pfcandvetoLpt10",       &pfcandvetoLpt10_,       "pfcandvetoLpt10/F");  
  outTree->Branch("pfcandvetoLmindrj10",   &pfcandvetoLmindrj10_,   "pfcandvetoLmindrj10/F");  

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

  outTree->Branch("pfTau15_leadPtcandID",        &pfTau15_leadPtcandID_,        "pfTau15_leadPtcandID/I");
  outTree->Branch("pfTau15"  , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &pfTau15_	);
  outTree->Branch("pfTau15_leadPtcand"  , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &pfTau15_leadPtcand_	);

  outTree->Branch("pfTau_leadPtcandID",        &pfTau_leadPtcandID_,        "pfTau_leadPtcandID/I");
  outTree->Branch("pfTau"  , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &pfTau_	);
  outTree->Branch("pfTau_leadPtcand"  , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &pfTau_leadPtcand_	);

  outTree->Branch("pfTauLoose_leadPtcandID",        &pfTauLoose_leadPtcandID_,        "pfTau_leadPtcandID/I");
  outTree->Branch("pfTauLoose"  , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &pfTauLoose_	);
  outTree->Branch("pfTauLoose_leadPtcand"  , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &pfTauLoose_leadPtcand_	);

  outTree->Branch("pfcandOS10"  , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &pfcandOS10_	);
  outTree->Branch("pfcandOS10looseZ"  , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &pfcandOS10looseZ_	);
  outTree->Branch("pfcand5looseZ"   , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &pfcand5looseZ_	);
  outTree->Branch("pfcanddir10"  , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &pfcanddir10_	);
  outTree->Branch("pfcandveto10"  , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &pfcandveto10_	);
  outTree->Branch("pfcandvetoL10"  , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &pfcandvetoL10_	);
  outTree->Branch("jet"	      , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &jet_	);

  outTree->Branch("nonisoel"  , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &nonisoel_	);
  outTree->Branch("nonisomu"  , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &nonisomu_	);
  outTree->Branch("t"         , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &t_   	);
  outTree->Branch("tbar"      , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &tbar_   	);
  outTree->Branch("ttbar"     , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &ttbar_   	);

  outTree->Branch("lep_t"     , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &lep_t_   	);
  outTree->Branch("lep_tbar"  , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &lep_tbar_  );
  outTree->Branch("stop_t"    , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &stop_t_   	);
  outTree->Branch("stop_tbar" , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &stop_tbar_ );
  outTree->Branch("neutralino_t"    , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &neutralino_t_    );
  outTree->Branch("neutralino_tbar" , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &neutralino_tbar_ );
  outTree->Branch("lep_t_id",            &lep_t_id_,            "lep_t_id/I");  
  outTree->Branch("lep_tbar_id",         &lep_tbar_id_,         "lep_tbar_id/I");  

  //  outTree->Branch("candidates", "std::vector<Candidate>", &candidates_);
  //  outTree->Branch("jets", "std::vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > >", &jets_ );
  //  outTree->Branch("btag", "std::vector<float>", &btag_ );

  outTree->Branch("pfjets"    , "std::vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > >", &pfjets_ );
  outTree->Branch("pfjets_genJet_"    , "std::vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > >", &pfjets_genJet_ );
  outTree->Branch("pfjets_csv", "std::vector<float>", &pfjets_csv_ );
  outTree->Branch("pfjets_chm", "std::vector<float>", &pfjets_chm_ );
  outTree->Branch("pfjets_neu", "std::vector<float>", &pfjets_neu_ );
  outTree->Branch("pfjets_l1corr",  "std::vector<float>", &pfjets_l1corr_   );
  outTree->Branch("pfjets_corr",    "std::vector<float>", &pfjets_corr_     );
  outTree->Branch("pfjets_mc3",     "std::vector<int>"  , &pfjets_mc3_      ); 
  outTree->Branch("pfjets_mcflavorAlgo", "std::vector<int>", &pfjets_mcflavorAlgo_ ); 
  outTree->Branch("pfjets_mcflavorPhys", "std::vector<int>", &pfjets_mcflavorPhys_ ); 
  outTree->Branch("pfjets_uncertainty" , "std::vector<float>", &pfjets_uncertainty_ );

  outTree->Branch("pfjets_flav",    "std::vector<int>"  , &pfjets_flav_     ); 
  outTree->Branch("pfjets_lrm", "std::vector<float>", &pfjets_lrm_ );
  outTree->Branch("pfjets_lrm2", "std::vector<float>", &pfjets_lrm2_ );
  outTree->Branch("pfjets_qgtag",   "std::vector<float>", &pfjets_qgtag_    );
  outTree->Branch("pfjets_genJetDr","std::vector<float>", &pfjets_genJetDr_ );
  outTree->Branch("pfjets_sigma",   "std::vector<float>", &pfjets_sigma_    );
  outTree->Branch("pfjets_lepjet",  "std::vector<int>"  , &pfjets_lepjet_   );
 
  outTree->Branch("pfjets_beta",      "std::vector<float>", &pfjets_beta_      );
  outTree->Branch("pfjets_beta2",     "std::vector<float>", &pfjets_beta2_     );
  outTree->Branch("pfjets_beta_0p1",  "std::vector<float>", &pfjets_beta_0p1_  );
  outTree->Branch("pfjets_beta_0p2",  "std::vector<float>", &pfjets_beta_0p2_  );
  outTree->Branch("pfjets_beta2_0p1", "std::vector<float>", &pfjets_beta2_0p1_ );
  outTree->Branch("pfjets_beta2_0p5", "std::vector<float>", &pfjets_beta2_0p5_ );
  // outTree->Branch("pfjets_beta_0p15", "std::vector<float>", &pfjets_beta_0p15_ );
  // outTree->Branch("pfjets_beta2_0p15","std::vector<float>", &pfjets_beta2_0p15_);
  // outTree->Branch("pfjets_beta_0p2",  "std::vector<float>", &pfjets_beta_0p2_  );
  // outTree->Branch("pfjets_beta2_0p2", "std::vector<float>", &pfjets_beta2_0p2_ );
  outTree->Branch("pfjets_mvaPUid",      "std::vector<float>", &pfjets_mvaPUid_ );
  outTree->Branch("pfjets_mva5xPUid",      "std::vector<float>", &pfjets_mva5xPUid_ );
  outTree->Branch("pfjets_mvaBeta",      "std::vector<float>", &pfjets_mvaBeta_ );

  outTree->Branch("genbs"    , "std::vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > >", &genbs_ );

  outTree->Branch("genps_pdgId"        ,   "std::vector<int>"    , &genps_pdgId_         );
  outTree->Branch("genps_firstMother"  ,   "std::vector<int>"    , &genps_firstMother_   );
  outTree->Branch("genps_energy"       ,   "std::vector<float>"  , &genps_energy_        );
  outTree->Branch("genps_pt"           ,   "std::vector<float>"  , &genps_pt_            );
  outTree->Branch("genps_eta"          ,   "std::vector<float>"  , &genps_eta_           );
  outTree->Branch("genps_phi"          ,   "std::vector<float>"  , &genps_phi_           );

}

//--------------------------------------------------------------------

/*

std::vector<float> singleLeptonLooper::totalIso( int thisPf , float coneR , float dz_thresh ){

  float iso = 0.; 

  float emiso = 0.;
  float nhiso = 0.;
  float chiso = 0.;

  for (int ipf = 0; ipf < (int)cms2.pfcands_p4().size(); ipf++) {

    if( ipf == thisPf                 ) continue; // skip this PFCandidate
    float dR = ROOT::Math::VectorUtil::DeltaR( pfcands_p4().at(ipf) , pfcands_p4().at(thisPf) );
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
      
      // if(!isGoodDAVertex(ivtx)) continue;

      // float mydz = dz_trk_vtx(itrk,ivtx);
      
      // if (fabs(mydz) < fabs(mindz)) {
      //   mindz = mydz;
      //   vtxi = ivtx;
      // }
         
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

*/


float getdltrigweight(int id1, int id2)
{ 
  if (abs(id1)==11 && abs(id2)==11) return 0.95;
  if (abs(id1)==13 && abs(id2)==13) return 0.88;
  if (abs(id1)!=abs(id2)) return 0.92;
  return -999.;

}

float getsltrigweight(int id1, float pt, float eta) 
{

  //--------------------------------------------------------------------------------------
  // This function returns the trigger efficiencies for HLT_IsoMu24 and HLT_Ele27_WP80
  // Efficiencies are calculated with the full 2012 data sample using Z tag-and-probe
  // 3 arguments:
  // id1 : lepton PDG ID, abs(id1)=11 (electron),  abs(id1)=13 (muon)  
  // pt  : lepton pt
  // eta : lepton eta
  //--------------------------------------------------------------------------------------

  //electron efficiencies
  if ( abs(id1)==11 ) {
    
    if ( fabs(eta)<1.5) {
      if ( pt>=20  && pt<22  ) return 0.00;
      if ( pt>=22  && pt<24  ) return 0.00;
      if ( pt>=24  && pt<26  ) return 0.00;
      if ( pt>=26  && pt<28  ) return 0.07;
      if ( pt>=28  && pt<30  ) return 0.57;
      if ( pt>=30  && pt<32  ) return 0.85;
      if ( pt>=32  && pt<34  ) return 0.88;
      if ( pt>=34  && pt<36  ) return 0.89;
      if ( pt>=36  && pt<38  ) return 0.91;
      if ( pt>=38  && pt<40  ) return 0.92;
      if ( pt>=40  && pt<50  ) return 0.94;
      if ( pt>=50  && pt<60  ) return 0.95;
      if ( pt>=60  && pt<80  ) return 0.96;
      if ( pt>=80  && pt<100 ) return 0.96;
      if ( pt>=100 && pt<150 ) return 0.97;
      if ( pt>=150 && pt<200 ) return 0.97;
      if ( pt>=200           ) return 0.97;
    } 

    else if ( fabs(eta)>=1.5 && fabs(eta)<2.1) {
      if ( pt>=20  && pt<22  ) return 0.00; 
      if ( pt>=22  && pt<24  ) return 0.00; 
      if ( pt>=24  && pt<26  ) return 0.03; 
      if ( pt>=26  && pt<28  ) return 0.22; 
      if ( pt>=28  && pt<30  ) return 0.52; 
      if ( pt>=30  && pt<32  ) return 0.65; 
      if ( pt>=32  && pt<34  ) return 0.70; 
      if ( pt>=34  && pt<36  ) return 0.72; 
      if ( pt>=36  && pt<38  ) return 0.74; 
      if ( pt>=38  && pt<40  ) return 0.75; 
      if ( pt>=40  && pt<50  ) return 0.77; 
      if ( pt>=50  && pt<60  ) return 0.79; 
      if ( pt>=60  && pt<80  ) return 0.79; 
      if ( pt>=80  && pt<100 ) return 0.80; 
      if ( pt>=100 && pt<150 ) return 0.82; 
      if ( pt>=150 && pt<200 ) return 0.83; 
      if ( pt>=200           ) return 0.85; 
    }

    // else{
    //   std::cout << "WARNING: electron eta " << eta << " is out of range. Return trigger efficiency = 1" << std::endl;
    // }

  } 

  //muon efficiencies
  else if ( abs(id1)==13 ) {

    if ( fabs(eta)<0.8 ) {
      if (pt>=20  && pt<22  ) return 0.00;
      if (pt>=22  && pt<24  ) return 0.02;
      if (pt>=24  && pt<26  ) return 0.87;
      if (pt>=26  && pt<28  ) return 0.90;
      if (pt>=28  && pt<30  ) return 0.91;
      if (pt>=30  && pt<32  ) return 0.91;
      if (pt>=32  && pt<34  ) return 0.92;
      if (pt>=34  && pt<36  ) return 0.93;
      if (pt>=36  && pt<38  ) return 0.93;
      if (pt>=38  && pt<40  ) return 0.93;
      if (pt>=40  && pt<50  ) return 0.94;
      if (pt>=50  && pt<60  ) return 0.95;
      if (pt>=60  && pt<80  ) return 0.95;
      if (pt>=80  && pt<100 ) return 0.94;
      if (pt>=100 && pt<150 ) return 0.94;
      if (pt>=150 && pt<200 ) return 0.93;
      if (pt>=200           ) return 0.92;
    } 

    else if ( fabs(eta)>=0.8 && fabs(eta)<1.5 ) {
      if (pt>=20  && pt<22  ) return 0.00;
      if (pt>=22  && pt<24  ) return 0.05;
      if (pt>=24  && pt<26  ) return 0.78;
      if (pt>=26  && pt<28  ) return 0.80;
      if (pt>=28  && pt<30  ) return 0.81;
      if (pt>=30  && pt<32  ) return 0.82;
      if (pt>=32  && pt<34  ) return 0.82;
      if (pt>=34  && pt<36  ) return 0.82;
      if (pt>=36  && pt<38  ) return 0.83;
      if (pt>=38  && pt<40  ) return 0.83;
      if (pt>=40  && pt<50  ) return 0.84;
      if (pt>=50  && pt<60  ) return 0.84;
      if (pt>=60  && pt<80  ) return 0.84;
      if (pt>=80  && pt<100 ) return 0.84;
      if (pt>=100 && pt<150 ) return 0.84;
      if (pt>=150 && pt<200 ) return 0.83;
      if (pt>=200           ) return 0.83;
    } 

    else if ( fabs(eta)>=1.5 && fabs(eta)<2.1 ) {
      if (pt>=20 && pt<22   ) return 0.00;
      if (pt>=22 && pt<24   ) return 0.10;
      if (pt>=24 && pt<26   ) return 0.76;
      if (pt>=26 && pt<28   ) return 0.78;
      if (pt>=28 && pt<30   ) return 0.79;
      if (pt>=30 && pt<32   ) return 0.80;
      if (pt>=32 && pt<34   ) return 0.81;
      if (pt>=34 && pt<36   ) return 0.81;
      if (pt>=36 && pt<38   ) return 0.82;
      if (pt>=38 && pt<40   ) return 0.82;
      if (pt>=40 && pt<50   ) return 0.83;
      if (pt>=50 && pt<60   ) return 0.83;
      if (pt>=60 && pt<80   ) return 0.84;
      if (pt>=80 && pt<100  ) return 0.84;
      if (pt>=100 && pt<150 ) return 0.84;
      if (pt>=150 && pt<200 ) return 0.82;
      if (pt>=200           ) return 0.83;
    }

    // else{
    //   std::cout << "WARNING: muon eta " << eta << " is out of range. Return trigger efficiency = 1" << std::endl;
    // }
  }

  else{
    std::cout << "WARNING: unrecognized lepton id " << id1 << ". Return trigger efficiency = 1" << std::endl;
  }

  return 1.;

}
