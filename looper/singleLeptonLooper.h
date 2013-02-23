#ifndef singleLeptonLooper_h
#define singleLeptonLooper_h

#include "TFitter.h"

#include <vector>
#include <list>
#include <string>
#include <map>
#include "Math/LorentzVector.h"
#include "Math/PxPyPzE4D.h"

#include "stopUtils.h"

#include "../macros/Core/MT2Utility.h"
#include "../macros/Core/mt2bl_bisect.h"
#include "../macros/Core/mt2w_bisect.h"
//#include "Candidate.h"                                                                                                                                                                 
class Candidate : public TObject {
 public:
  float chi2, mt2w, mt2bl, mt2b;
  int j1, j2, bi, oi;
  float k1, k2;
  bool match;

  ClassDef(Candidate, 2)
    };


#include "../macros/Core/PartonCombinatorics.h"

#include "../CORE/SimpleFakeRate.h" // will .h be ok? lets see.. 101007
#include "../CORE/QuarkGluonTagger/QuarkGluonTagger.h"
#include "../CORE/metSelections.h"

//#include "../CORE/topmass/ttdilepsolve.h" REPLACETOPMASS

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > LorentzVector;
//typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > P4;
typedef vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > VofP4;
//typedef map<unsigned int, unsigned int> m_uiui;

class  TChain;
class  TH1F;
class  TH2F;
class  TRandom3;
class  TTree;
struct metStruct;

/* ------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------- */
const double PTMIN_J1   = 25;
const double PTMIN_J2   = 25;
const double PTMIN_BTAG = 30;
const double PTMIN_OTAG = 30;
const double PTMIN_B    = 30;  //  This Two should be tigther than the  
const double PTMIN_O    = 30;  //  b-tagged versions.
const double PDG_TOP_MASS = 173.5;
const double PDG_W_MASS = 80.385;
const double BTAG_MIN = 0.679;

const bool __SORT = true;

/* ------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------- */


class singleLeptonLooper
{
    public: 
        singleLeptonLooper();
        ~singleLeptonLooper() {}

        /* enum TrigEnum { e_highpt = 0, e_lowpt };   */
        /* // e_highpt  :   high pt dilepton triggers, 20,10   */
        /* // e_lowpt   :   dilepton-HT cross triggers, 10,5             */
        /* enum JetTypeEnum { e_JPT = 0, e_calo , e_pfjet }; */
        /* // e_JPT     :   jpt jets */
        /* // e_calo    :   l1 and l2 corrected calo jets */
        /* // e_pfjet   :   corrected pfjets */
        /* enum MetTypeEnum { e_tcmet = 0, e_muon, e_muonjes , e_pfmet }; */
        /* // e_tcmet   :   track corrected met */
        /* // e_muon    :   calo met with muon corrections */
        /* // e_muonjes :   calo met with muon and jet energy scale corrections */
        /* // e_pfmet   :   particle-flow met */
        /* enum ZVetoEnum   { e_standard = 0, e_allzveto, e_nozveto, e_selectz }; */
        /* // e_standard:   apply Z-veto to same-flavor pairs */
        /* // e_allzveto:   apply Z-veto regardless of lepton flavor */
        /* // e_nozveto :   no Z-veto */
        /* // e_selectz :   select Z by requiring SF OS pair in Z mass window */
        enum FREnum   { e_qcd = 0, e_wjets };
        // e_qcd     :   derive prediction for 2 fake leptons
        // e_wjets   :   derive prediction for 1 real and one fake lepton
	
        int  ScanChain(TChain *chain, char *prefix = "", float kFactor = 1., 
		       int prescale = 1., float lumi = 1.,
                       FREnum frmode  = e_wjets,
                       bool doFakeApp = false
                       );
        void BookHistos (char *prefix);
	void InitBaby();
	//	float dz_trk_vtx( const unsigned int trkidx, const unsigned int vtxidx = 0 );

        // Set globals
        void set_susybaseline (bool  b)    { g_susybaseline = b; }
        void set_createTree   (bool  b)    { g_createTree   = b; }
        void set_useBitMask   (bool  b)    { g_useBitMask   = b; }
        void set_version      (const char* v)    { g_version      = v; }
	void set_json         (const char* v)    { g_json         = v; }        
        //void set_trigger      (TrigEnum t) { g_trig         = t; } 

        // Baby ntuple methods
        void makeTree (char *prefix,bool doFakeApp, FREnum frmode );
	float stopPairCrossSection( float stopmass );
        void closeTree ();
	//	float trackIso( int thisPf , float coneR = 0.3 , float dz_thresh = 0.05 , bool dovtxcut = false , float pt_thresh = 0.0);
	std::vector<float> trackIsoPtRanges( int thisPf , float coneR = 0.3 , float dz_thresh = 0.05 );
	std::vector<float> totalIso( int thisPf , float coneR = 0.3 , float dz_thresh = 0.05 );
	//pair<float,float> getPhiCorrMET( float met, float metphi, float sumet, bool ismc, bool is8TeV = false);
	//	pair<float,float> getPhiCorrMET( float met, float metphi, int nvtx, bool ismc);
	//	pair<float,float> getTrackerMET( P4 *lep, double deltaZCut = 0.1, bool dolepcorr = true );
	bool initialized;
	TH1D*   stop_xsec_hist;
	TFile*  stop_xsec_file;
	//3D Vertex weight
	//	double Weight3D[50][50][50];

    private:

        // Globals
        bool  g_susybaseline;
        bool  g_createTree;
        bool  g_useBitMask;
        const char* g_version;
	const char* g_json;      
	//TrigEnum g_trig;
        TRandom3 *random3_;

	//PDF information
        Float_t pdfQ_;
        Float_t pdfx1_; 
	Float_t pdfx2_;
        Int_t pdfid1_; 
	Int_t pdfid2_;

	Int_t   eldup_;    

	// MC truth lepton info
	Int_t   mcid1_;    
	Int_t   mcid2_;    
	Int_t   lep_t_id_;
	Int_t   lep_tbar_id_;
	Int_t	mcdecay1_; 
	Int_t   mcdecay2_; 
	Int_t	mcndec1_; 
	Int_t   mcndec2_;
	Int_t	mcndeckls1_; 
	Int_t   mcndeckls2_; 
	Int_t	mcndecem1_; 
	Int_t   mcndecem2_;  
	Float_t mctaudpt1_;
	Float_t mctaudpt2_;
	Int_t   mctaudid1_;    
	Int_t   mctaudid2_;    
	Float_t mcdr1_;    
	Float_t mcdr2_;    

	// isolated track vars
	Float_t trkpt5_;
	Float_t mleptrk5_;
	Float_t trkreliso5_;
	Float_t trkpt10_;
	Float_t mleptrk10_;
	Float_t trkreliso10_;

	Float_t trkpt5loose_;
	Float_t trkreliso5loose_;
	Float_t trkpt10loose_;
	Float_t trkreliso10loose_;

	//min pt cut on particles 
	//included in the isolation
	Float_t trkpt10pt0p1_;
	Float_t trkreliso10pt0p1_;
	Float_t trkpt10pt0p2_;
	Float_t trkreliso10pt0p2_;
	Float_t trkpt10pt0p3_;
	Float_t trkreliso10pt0p3_;
	Float_t trkpt10pt0p4_;
	Float_t trkreliso10pt0p4_;
	Float_t trkpt10pt0p5_;
	Float_t trkreliso10pt0p5_;
	Float_t trkpt10pt0p6_;
	Float_t trkreliso10pt0p6_;
	Float_t trkpt10pt0p7_;
	Float_t trkreliso10pt0p7_;
 	Float_t trkpt10pt0p8_;
	Float_t trkreliso10pt0p8_;
	Float_t trkpt10pt0p9_;
	Float_t trkreliso10pt0p9_;
	Float_t trkpt10pt1p0_;
	Float_t trkreliso10pt1p0_;

	// extra pfcand vars 
        Int_t pfcandid5_;     
        Int_t pfcandid10_;     
        Int_t pfcanddirid10_;     
        Int_t pfcandvetoid10_;     
        Int_t pfcandvetoLid10_;     
        Int_t pfcandidOS10_;  
        Int_t pfcandidOS10looseZ_;  
        Int_t pfcandid5looseZ_;     
   
        Float_t pfcandiso5_;     
        Float_t pfcandiso10_;     
        Float_t pfcanddiriso10_;     
        Float_t pfcandvetoiso10_;     
        Float_t pfcandvetoLiso10_;     
        Float_t pfcandisoOS10_;     
        Float_t pfcandisoOS10looseZ_;     
        Float_t pfcandiso5looseZ_;     

        Float_t pfcandpt5_;
        Float_t pfcandpt10_;
        Float_t pfcanddirpt10_;
        Float_t pfcandvetopt10_;
        Float_t pfcandvetoLpt10_;
        Float_t pfcandptOS10_;
        Float_t pfcandptOS10looseZ_;
        Float_t pfcandpt5looseZ_;

        Float_t pfcanddz5_;
        Float_t pfcanddz10_;
        Float_t pfcanddirdz10_;
        Float_t pfcandvetodz10_;
        Float_t pfcandvetoLdz10_;
        Float_t pfcanddzOS10_;
        Float_t pfcanddzOS10looseZ_;
        Float_t pfcanddz5looseZ_;

        Float_t pfcandmindrj5_;
        Float_t pfcandmindrj10_;
        Float_t pfcanddirmindrj10_;
        Float_t pfcandvetomindrj10_;
        Float_t pfcandvetoLmindrj10_;

	Float_t pfcandpt10pt0p1_;
	Float_t pfcandiso10pt0p1_;
	Float_t pfcandpt10pt0p2_;
	Float_t pfcandiso10pt0p2_;
	Float_t pfcandpt10pt0p3_;
	Float_t pfcandiso10pt0p3_;
	Float_t pfcandpt10pt0p4_;
	Float_t pfcandiso10pt0p4_;
	Float_t pfcandpt10pt0p5_;
	Float_t pfcandiso10pt0p5_;
	Float_t pfcandpt10pt0p6_;
	Float_t pfcandiso10pt0p6_;
	Float_t pfcandpt10pt0p7_;
	Float_t pfcandiso10pt0p7_;
 	Float_t pfcandpt10pt0p8_;
	Float_t pfcandiso10pt0p8_;
	Float_t pfcandpt10pt0p9_;
	Float_t pfcandiso10pt0p9_;
	Float_t pfcandpt10pt1p0_;
	Float_t pfcandiso10pt1p0_;

	// btag variables
	Int_t   nbtagsssv_;     
	Int_t   nbtagstcl_;     
	Int_t   nbtagstcm_;     
	Int_t   nbtagscsvl_;    
	Int_t   nbtagscsvm_;    
	Int_t   nbtagscsvt_;    
	Int_t   nbtagsssvcorr_;     
	Int_t   nbtagstclcorr_;     
	Int_t   nbtagstcmcorr_;     
	Int_t   nbtagscsvlcorr_;    
	Int_t   nbtagscsvmcorr_;    
	Int_t   nbtagscsvtcorr_;    

	// pfjet counters
	Int_t   npfjets30_;
	Int_t   npfjets30lepcorr_;
	Float_t knjets_;

	//rho correction
	Float_t rhovor_;

	// pfht vars
	Float_t htpf30_;

	// matched lepton vars
	Int_t   mlepid_;
	Int_t   mleppassid_;
	Int_t   mleppassiso_;
	Float_t mlepiso_;
	Float_t mlepdr_;

	Float_t pflepdr_;
	Float_t pflepiso_;
	Float_t pfleppt_;
	Float_t pflepmindrj_;
	Float_t pftauddr_;
	Float_t pftaudiso_;
	Float_t pftaudpt_;
	Float_t pftaudmindrj_;
	
	//lepton jet matching variable
	Float_t lep1pfjetdr_;
	Float_t lep2pfjetdr_;

	// HLT variables
	/* Int_t   ldi_; */
	/* Int_t   ltri_; */
	/* Int_t   smu_; */
	/* Int_t   smu30_; */
	Int_t   trgmu1_;
	Int_t   trgmu2_;
	Int_t   trgel1_;
	Int_t   trgel2_;
	//Int_t   dil_;

	// MC truth vars
	Int_t   npartons_;
	Int_t   nwzpartons_;
	Float_t maxpartonpt_;
	Float_t ptt_;
	Float_t pttbar_;
	Float_t ptttbar_;
	Float_t mttbar_;
	Float_t etattbar_;
	Float_t mgcor_;
	Int_t   wflav_;

	// Type1 pfmet
	Float_t t1met10_;
	Float_t t1met20_;
	Float_t t1met30_;
	Float_t t1met10phi_;
	Float_t t1met20phi_;
	Float_t t1met30phi_;
	Float_t t1met10mt_;
	Float_t t1met20mt_;
	Float_t t1met30mt_;
	Float_t lepmetpt_;
	Float_t lept1met10pt_;

	Float_t t1met10s_;
	Float_t t1met10sphi_;
	Float_t t1met10smt_;

	//phi corrected type1 mets
	Float_t t1metphicorr_;
	Float_t t1metphicorrphi_;
	Float_t t1metphicorrphiup_;
	Float_t t1metphicorrphidn_;
	Float_t t1metphicorrmt_;
	Float_t t1metphicorrlepmt_;
	Float_t t1metphicorrlep_;
	Float_t t1metphicorrlepphi_;
	Float_t t1metphicorrup_;
	Float_t t1metphicorrdn_;
	Float_t t1metphicorrmtup_;
	Float_t t1metphicorrmtdn_;

	//smeared and phicorred type1 mets
	Float_t t1metphicorr_smear_gen_;
	Float_t t1metphicorrphi_smear_gen_;
	Float_t t1metphicorrmt_smear_gen_;

	//official prescription
	Float_t t1met_off_;
	Float_t t1metphi_off_;
	Float_t t1metmt_off_;
	Float_t t1metphicorr_off_;
	Float_t t1metphicorrphi_off_;
	Float_t t1metphicorrmt_off_;

	// assorted p4's
	LorentzVector*  t_;   
	LorentzVector*  tbar_;   
	LorentzVector*  lep_t_;   
	LorentzVector*  lep_tbar_;   
	LorentzVector*  stop_t_;   
	LorentzVector*  stop_tbar_;
	LorentzVector*  neutralino_t_;   
	LorentzVector*  neutralino_tbar_;    
	LorentzVector*  ttbar_;   
	LorentzVector*  mlep_;   
	LorentzVector*  mclep1_;   
	LorentzVector*  mclep2_;   
	LorentzVector*  mctaud1_;   
	LorentzVector*  mctaud2_;  
	LorentzVector*  mctaudvis1_;   
	LorentzVector*  mctaudvis2_;  
        LorentzVector*  pflep_;
        LorentzVector*  pftaud_;
        LorentzVector*  pfcand5_;
        LorentzVector*  pfcand10_;
        LorentzVector*  pfcanddir10_;
        LorentzVector*  pfcandveto10_;
        LorentzVector*  pfcandvetoL10_;
        LorentzVector*  pfcandOS10_;
        LorentzVector*  pfcandOS10looseZ_;
        LorentzVector*  pfcand5looseZ_;
        LorentzVector*  lep1_;
        LorentzVector*  lep2_;
        LorentzVector*  trklep1_;
        LorentzVector*  trklep2_;
        LorentzVector*  gfitlep1_;
        LorentzVector*  gfitlep2_;
        LorentzVector*  lepp_;
        LorentzVector*  lepm_;
        LorentzVector*  pflep1_;
        LorentzVector*  pflep2_;
        LorentzVector*  leppfjet1_;
        LorentzVector*  leppfjet2_;
        LorentzVector*  dilep_;
        LorentzVector*  jet_; 
	LorentzVector*  mcnu_;
	LorentzVector*  mclep_;

	//        std::vector<Candidate>  candidates_;
        VofP4 jets_;
	std::vector<float> btag_;

        VofP4 pfjets_;
        VofP4 pfjets_genJet_;
        std::vector<float> pfjets_csv_;
	// CHM and NEU
        std::vector<float> pfjets_chm_;
        std::vector<float> pfjets_neu_;
	//jet corrections and ID
	std::vector<float> pfjets_l1corr_;
	std::vector<float> pfjets_corr_;
	//status 3 parton matching
	std::vector<int>   pfjets_mc3_;
	std::vector<int>   pfjets_flav_;
	std::vector<int>   pfjets_mcflavorAlgo_;
	std::vector<int>   pfjets_mcflavorPhys_;

	//qg tagging
	std::vector<float> pfjets_qgtag_;
	//jet shape
	std::vector<float> pfjets_lrm_;
	std::vector<float> pfjets_lrm2_;
	//gen jet matching
	std::vector<float> pfjets_genJetDr_;
	//truth lepton matching
        std::vector<int>   pfjets_lepjet_;
	//jet resolution
	std::vector<float> pfjets_sigma_;
	//beta variables
	std::vector<float> pfjets_beta_;
	std::vector<float> pfjets_beta2_;
        std::vector<float> pfjets_beta_0p1_;
        std::vector<float> pfjets_beta_0p2_;
        std::vector<float> pfjets_beta2_0p5_;
        std::vector<float> pfjets_beta2_0p1_;
        /* std::vector<float> pfjets_beta_0p15_;  */
        /* std::vector<float> pfjets_beta2_0p15_; */
        /* std::vector<float> pfjets_beta_0p2_;   */
        /* std::vector<float> pfjets_beta2_0p2_;  */
	std::vector<float> pfjets_mvaPUid_;
	std::vector<float> pfjets_mvaBeta_;

	float chi2min_;
	float chi2minprob_;

	float mt2bmin_;
	float mt2blmin_;
	float mt2wmin_;


	Int_t hyptype_;

	Int_t mm_;
	Int_t mmtk_;
	Int_t me_;
	Int_t em_;
	Int_t mu_;
	Int_t ee_;
	Int_t isomu24_;
	Int_t ele27wp80_;
	
 	LorentzVector*  nonisoel_;   
 	LorentzVector*  nonisomu_;   

        // Baby ntuple variables
 	Float_t mcmln_;
 	Float_t mcmtln_;
	Float_t mjj_;
	Float_t dphilm_;
	Float_t mG_;
	Float_t x_;
	Float_t mL_;
	Float_t ecalveto1_;
	Float_t ecalveto2_;
	Float_t hcalveto1_;
	Float_t hcalveto2_;
        TFile  *outFile;
        TTree  *outTree;
	Int_t   acc_2010_;
	Int_t   acc_highmet_;
	Int_t   acc_highht_;
	Float_t pthat_;
	Float_t qscale_;
        Float_t weight_;
        Float_t trgeff_;
        Float_t mutrigweight_;
        Float_t mutrigweight2_;
        Float_t sltrigweight_;
        Float_t dltrigweight_;
        Float_t pfmetsig_;
        Float_t smeff_;
        Float_t k_;
        Float_t mllgen_;
        Float_t ptzgen_;
        Float_t ptwgen_;
        Float_t dphijm_;
        Float_t costhetaweight_;
	Int_t   njpt_;
	Float_t htjpt_;

	Int_t   csc_;
	Int_t   hcallaser_;
	Int_t   ecaltp_;
	Int_t   trkfail_;
	Int_t   eebadsc_;
	Int_t   hbhenew_;
        Int_t   hbhe_;

	Int_t   jetid_;
	Int_t   jetid30_;
        Int_t   mull_;
        Int_t   json_;
        Int_t   mult_;
        Int_t   mullgen_;
        Int_t   multgen_;
        Int_t   nlep_;
        Int_t   ngoodlep_;
        Int_t   ngoodel_;
        Int_t   nosel_;
        Int_t   ngoodmu_;
        Int_t   proc_;
        Int_t   leptype_;
        Int_t   ngenjets_;
        Int_t   njetsUp_;
        Int_t   npfjets25_;
        Int_t   njetsDown_;

	Float_t calomet_;
	Float_t calometphi_;
	Float_t trkmet_nolepcorr_;
	Float_t trkmetphi_nolepcorr_;
	Float_t trkmet_;
	Float_t trkmetphi_;
	Float_t trkmetproj_;
	Float_t trkmet4_;
	Float_t trkmet4phi_;
	Float_t trkmet4proj_;
	Float_t trkmet8_;
	Float_t trkmet8phi_;
	Float_t trkmet8proj_;
	Float_t trkjetmet_;
	Float_t trkjetmetphi_;
	Float_t trkjetmetproj_;
        Float_t htUp_;
        Float_t htDown_;
        Int_t   ntruepu_;
        Int_t   npu_;
        Int_t   npuMinusOne_;
        Int_t   npuPlusOne_;
        Int_t   nvtx_;
        Int_t	indexfirstGoodVertex_;
        Float_t dilmass_;
        Float_t topmass_;
        Float_t tcmet_;
        Float_t tcmet00_;
        Float_t tcmet10_;
        Float_t tcmet20_;
        Float_t tcmet30_;
        Float_t tcmet40_;
        Float_t tcmet50_;
        Float_t genmet_;
        Float_t gensumet_;
        Float_t genmetphi_;
        Float_t mucormet_;
        Float_t mucorjesmet_;
        Float_t pfmet_;
        Float_t pfsumet_;
        Float_t pfmetveto_;
        Float_t pfmetphi_;
        Float_t tcmet_35X_;
        Float_t tcmet_event_;
        Float_t tcmet_looper_;
        Float_t tcmetUp_;
        Float_t tcmetDown_;
        Float_t tcmetTest_;
        Float_t pfmetUp_;
        Float_t pfmetDown_;
        Float_t pfmetTest_;
        Float_t tcsumet_;
        Float_t tcmetphi_;
        Float_t mt2_;
        Float_t mt2j_;
        Float_t mt2jcore_;
        Float_t sumjetpt_;
        Float_t dileta_;
        Float_t dilpt_;
        Float_t dildphi_;
        Float_t vecjetpt_;
        Int_t   pass_;
        Int_t   passz_;
        Float_t m0_;
        Float_t m12_;
        Float_t ptl1_;
	Float_t lep1chi2ndf_;
	Float_t lep2chi2ndf_;
	Float_t lep1dpt_;
	Float_t lep2dpt_;
	Int_t   id1_;
	Int_t   id2_;
	Int_t   leptype1_;
	Int_t   leptype2_;
	Int_t   w1_;
	Int_t   w2_;
	Float_t iso1_;
	Float_t isont1_;
	Float_t isopf1_;
	Float_t isopfold1_;
	Float_t iso2_;
	Float_t isont2_;
	Float_t isopf2_;
        Float_t ptl2_;
        Float_t etal1_;
        Float_t etal2_;
        Float_t phil1_;
        Float_t phil2_;
        Float_t meff_;
        Float_t mt_;
        char    dataset_[200];
        UInt_t  run_;
        UInt_t  lumi_;
        UInt_t  event_;
	Float_t y_;
	Float_t ht_;
	Float_t htoffset_;
	Float_t htuncor_;
	Int_t   njetsuncor_;
	Int_t   njetsoffset_;
	Float_t htgen_;
	Float_t htpf_;
	Float_t ptjetraw_;
	Float_t ptjet23_;
	Float_t ptjetF23_;
	Float_t ptjetO23_;
	Float_t cosphijz_;
	Int_t   nels_;
	Int_t   nmus_;
	Int_t   ntaus_;
	Int_t   nleps_;
	Float_t nvtxweight_;
	Float_t n3dvtxweight_;
	Float_t etasc1_;
	Float_t etasc2_;
	Float_t eoverpin_;
	Float_t eoverpout_;
	Float_t dEtaIn_;
	Float_t dPhiIn_;
	Float_t sigmaIEtaIEta_;
	Float_t hOverE_;
	Float_t ooemoop_;
	Float_t d0vtx_;
	Float_t dzvtx_;
	Float_t expinnerlayers_;
	Float_t fbrem_;
	Float_t pfisoch_;
	Float_t pfisoem_;
	Float_t pfisonh_;
	Float_t eSC_;
	Float_t phiSC_;
	Float_t eSCRaw_;
	Float_t eSCPresh_;
	Float_t emjet10_;
	Float_t emjet20_;
	Float_t eoverpin2_;
	Float_t eoverpout2_;
	Float_t dEtaIn2_;
	Float_t dPhiIn2_;
	Float_t sigmaIEtaIEta2_;
	Float_t hOverE2_;
	Float_t ooemoop2_;
	Float_t d0vtx2_;
	Float_t dzvtx2_;
	Float_t expinnerlayers2_;
	Float_t fbrem2_;
	Float_t pfisoch2_;
	Float_t pfisoem2_;
	Float_t pfisonh2_;
	Float_t eSC2_;
	Float_t phiSC2_;
	Float_t eSCRaw2_;
	Float_t eSCPresh2_;



	//recoil
	Float_t dilrecoil_;
	Float_t dilrecoilparl_;
	Float_t dilrecoilperp_;

	Float_t ksusy_;
	Float_t ksusyup_;
	Float_t ksusydn_;
	Float_t xsecsusy_;
	Float_t xsecsusy2_;

	Float_t mbb_;
	Int_t   isdata_;

        // Lots and lots of histograms
        TH1F* h_PU_trkpt;
        TH1F* h_dz_vtx_trk;
};

#endif
