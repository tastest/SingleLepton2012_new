#ifndef StopTree_H
#define StopTree_H

#include "TFile.h"
#include "TTree.h"
#include "TError.h"
#include "TObject.h"

#include "Math/LorentzVector.h"

#include <cmath>
#include <vector>
#include "assert.h"

//#include "../../looper/Candidate.h"

using namespace std;
using namespace ROOT::Math;

//
// Ntuple structure:
//
// Ntuple content:

struct MT2struct {
  float mt2w;
  float mt2b;
  float mt2bl;
  float chi2;
};

class Candidate : public TObject {
public:
 float chi2, mt2w, mt2bl, mt2b;
 int j1, j2, bi, oi;
 float k1, k2;
 bool match;

 ClassDef(Candidate, 2)
};

typedef vector<Candidate> CANDIDATES;

class StopTree {

    public:
  typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > LorentzVector;

        /// variables
	char          dataset_[200];
	unsigned int  run_;
	unsigned int  lumi_;
	unsigned int  event_;
	int           nvtx_;
	float         nvtxweight_;
	float         weight_;
	float         xsecsusy_;
	float         mutrigweight_;
	float         rhovor_;
	float         mgcor_;
	int	      csc_;
	int	      hbhe_;
	int	      hcallaser_;
	int	      ecaltp_;
	int	      trkfail_;
	int	      eebadsc_;
	int	      hbhenew_;
	int           isomu24_;
	int           ele27wp80_;
	int           trgmu1_;
	int           trgel1_;
	int           mm_;
	int           me_;
	int           em_;
	int           ee_;
	int           ngoodlep_;
	int           leptype_;
	int           id1_;
	int           id2_;
	int           lep_t_id_;
	int           lep_tbar_id_;
	int           mcid1_;
	int           mcid2_;
	int           mcdecay2_;
	int           mcndec2_;
	int           npfjets30_;
	int           npfjets40_;
	float 	      mctaudpt2_;	
	float 	      lep1chi2ndf_;	
	float 	      lep1dpt_;	
	float         dilmass_;
	float         pfmet_;
	float         mt_;
	float         pfmetphi_;
	float         t1met10_;
	float         t1met10mt_;
	float         t1met10phi_;
	float         t1metphicorr_;
	float         t1metphicorrmt_;
	float         t1metphicorrphi_;
	float         t1metphicorrlep_;
	float         t1metphicorrlepmt_;
	float         t1metphicorrlepphi_;
	int           nbtagscsvm_;
	int           nbtagscsvmcorr_;
	float         pfcandpt10_;
	float         pfcandiso10_;
	float         trkpt10loose_;
	float         trkreliso10loose_;
	int           nleps_;
	int           nmus_;
	int           nels_;
	int           ntaus_;

	Float_t         mcmtln_;
	Int_t         nwzpartons_;
	Int_t         nbtagscsvl_;
	Float_t         trkmet_nolepcorr_;
	Float_t         trkmetphi_nolepcorr_;
	Float_t         trkmet_;
	Float_t         trkmetphi_;
	Float_t         isopf1_;
	Float_t         isopfold1_;
	Float_t         isopf2_;
	Float_t         eoverpin_;
	Float_t         eoverpout_;
	Float_t         dEtaIn_;
	Float_t         dPhiIn_;
	Float_t         sigmaIEtaIEta_;
	Float_t         hOverE_;
	Float_t         ooemoop_;
	Float_t         d0vtx_;
	Float_t         dzvtx_;
	Float_t         expinnerlayers_;
	Float_t         fbrem_;
	Float_t         pfisoch_;
	Float_t         pfisoem_;
	Float_t         pfisonh_;
	Float_t         eSC_;
	Float_t         phiSC_;
	Float_t         eSCRaw_;
	Float_t         eSCPresh_;
	Float_t         etasc1_;
	Float_t         iso1_;
	Float_t         isont1_;
	Float_t         iso2_;
	Float_t         isont2_;
	         
	LorentzVector lep1_;
	LorentzVector lep2_;
	LorentzVector t_;
	LorentzVector tbar_;
	LorentzVector lep_t_;
	LorentzVector lep_tbar_;
	LorentzVector mclep1_;
	LorentzVector mclep2_;
	LorentzVector stop_t_;
	LorentzVector stop_tbar_;
	LorentzVector pfcand10_;

	LorentzVector pflep1_;
	LorentzVector pflep2_;

        CANDIDATES* candidates_;
	//        vector<LorentzVector>* jets_;
        vector<LorentzVector>* pfjets_;
	vector<float> pfjets_csv_;
	vector<float> pfjets_qgtag_;
	vector<float> pfjets_mc3_;
	vector<float> pfjets_beta2_;
	vector<float> pfjets_beta_;
	vector<int>   pfjets_lepjet_;

    public:
        /// this is the main element
        TTree *tree_;
        TFile *f_;

        /// hold the names of variables to facilitate things (filled during Init)
        vector<string> variables_;
	
        /// default constructor  
	StopTree() :  lep1Ptr_(&lep1_), lep2Ptr_(&lep2_), tPtr_(&t_), tbarPtr_(&tbar_), stop_tPtr_(&stop_t_), stop_tbarPtr_(&stop_tbar_), lep_tPtr_(&lep_t_), lep_tbarPtr_(&lep_tbar_), mclep1Ptr_(&mclep1_), mclep2Ptr_(&mclep2_), pfcand10Ptr_(&pfcand10_), pflep1Ptr_(&pflep1_), pflep2Ptr_(&pflep2_), pfjets_csv_Ptr_(&pfjets_csv_), pfjets_qgtag_Ptr_(&pfjets_qgtag_), pfjets_mc3_Ptr_(&pfjets_mc3_), pfjets_beta_Ptr_(&pfjets_beta_), pfjets_beta2_Ptr_(&pfjets_beta2_), pfjets_lepjet_Ptr_(&pfjets_lepjet_)  {}
 //StopTree() :  lep1Ptr_(&lep1_), lep2Ptr_(&lep2_), pfcand10Ptr_(&pfcand10_), jet1Ptr_(&pfjet1_), jet2Ptr_(&pfjet2_), jet3Ptr_(&pfjet3_), jet4Ptr_(&pfjet4_), jet5Ptr_(&pfjet5_), jet6Ptr_(&pfjet6_) {}
        /// default destructor
        ~StopTree(){ 
	  cout << "~StopTree()" << endl;
	  if (f_) f_->Close();  
	  cout << "~StopTree() done" << endl;
	  
        };

        /// initialize varibles and fill list of available variables
        void InitVariables();

        /// load a StopTree
        void LoadTree(const char* file){
            f_ = TFile::Open(file);
            assert(f_);
            tree_ = dynamic_cast<TTree*>(f_->Get("t"));
            assert(tree_);
        }

        /// create a StopTree
        void CreateTree(){
            tree_ = new TTree("t","stop babytuple");
            f_ = 0;
            InitVariables();
            //book the branches
	    tree_->Branch("dataset", 	        &dataset_, 		"dataset/I");
	    tree_->Branch("run", 	        &run_, 		        "run/I");
	    tree_->Branch("lumi", 	        &lumi_, 		"lumi/I");
	    tree_->Branch("event", 	        &event_, 		"event/I");
	    tree_->Branch("nvtx", 		&nvtx_, 		"nvtx/I");
	    tree_->Branch("nvtxweight", 	&nvtxweight_, 		"nvtxweight/F");    
	    tree_->Branch("weight", 		&weight_, 		"weight/F");	      	      
	    tree_->Branch("xsecsusy", 		&xsecsusy_, 		"xsecsusy/F");	      	      
	    tree_->Branch("mutrigweight", 	&mutrigweight_, 	"mutrigweight/F");    
	    tree_->Branch("rhovor", 		&rhovor_, 		"rhovor/F");	      	      
	    tree_->Branch("mgcor", 		&mgcor_, 		"mgcor/F");	 
	    tree_->Branch("csc", 		&csc_, 			"csc/I");	      	      
	    tree_->Branch("hbhe", 		&hbhe_, 		"hbhe/I");	      	      
	    tree_->Branch("hcallaser", 		&hcallaser_, 		"hcallaser/I");	      	      
	    tree_->Branch("ecaltp", 		&ecaltp_, 		"ecaltp/I");	      	      
	    tree_->Branch("trkfail", 		&trkfail_, 		"trkfail/I");	      	      
	    tree_->Branch("eebadsc", 		&eebadsc_, 		"eebadsc/I");	      	      
	    tree_->Branch("hbhenew", 		&hbhenew_, 		"hbhenew/I");	      	      
	    tree_->Branch("isomu24", 		&isomu24_, 		"isomu24/I");	      	      
	    tree_->Branch("ele27wp80", 		&ele27wp80_, 		"ele27wp80/I");	      	      
	    tree_->Branch("trgmu1", 		&trgmu1_, 		"trgmu1/I");	      	      
	    tree_->Branch("trgel1", 		&trgel1_, 		"trgel1/I");	      	      
	    tree_->Branch("mm", 		&mm_, 			"mm/I");	      	      
	    tree_->Branch("me", 		&me_, 			"me/I");	      	      
	    tree_->Branch("em", 		&em_, 			"em/I");	      	      
	    tree_->Branch("ee", 		&ee_, 			"ee/I");	      	      
	    tree_->Branch("ngoodlep", 		&ngoodlep_, 		"ngoodlep/I");            
	    tree_->Branch("leptype", 		&leptype_, 		"leptype/I");	      	      
	    tree_->Branch("id1", 		&id1_, 			"id1/I");	      	      
	    tree_->Branch("id2", 		&id2_, 			"id2/I");	      	      
	    tree_->Branch("lep_t_id", 		&lep_t_id_, 		"lep_t_id/I");	      	      
	    tree_->Branch("lep_tbar_id", 		&lep_tbar_id_, 		"lep_tbar_id/I");	      	      
	    tree_->Branch("mcid1", 		&mcid1_, 		"mcid1/I");	      	      
	    tree_->Branch("mcid2", 		&mcid2_, 		"mcid2/I");	      	      
	    tree_->Branch("mcdecay2", 		&mcdecay2_, 		"mcdecay2/I");	      	      
	    tree_->Branch("mcndec2", 		&mcndec2_, 		"mcndec2/I");	      	      
	    tree_->Branch("npfjets30", 		&npfjets30_, 		"npfjets30/I");        
	    tree_->Branch("npfjets40", 		&npfjets40_, 		"npfjets40/I");            
	    tree_->Branch("mctaudpt2", 		&mctaudpt2_, 		"mctaudpt2/F");	      	      
	    tree_->Branch("lep1chi2ndf", 	&lep1chi2ndf_, 		"lep1chi2ndf/F");	      	      
	    tree_->Branch("lep1dpt", 		&lep1dpt_, 		"lep1dpt/F");
	    tree_->Branch("dilmmas", 		&dilmass_, 		"dilmass/F");	      	      
	    tree_->Branch("pfmet", 		&pfmet_, 		"pfmet/F");	      	      
	    tree_->Branch("mt", 		&mt_, 			"mt/F");          
	    tree_->Branch("pfmetphi", 		&pfmetphi_, 		"pfmetphi/F");        	      	      
	    tree_->Branch("t1met10", 		&t1met10_, 		"t1met10/F");	      	      
	    tree_->Branch("t1met10mt", 		&t1met10mt_, 		"t1met10mt/F");          
	    tree_->Branch("t1met10phi", 	&t1met10phi_, 		"t1met10phi/F");        	      	      
	    tree_->Branch("nbtagscsvm", 	&nbtagscsvm_, 		"nbtagscsvm/I");          
	    tree_->Branch("nbtagscsvmcorr", 	&nbtagscsvmcorr_, 	"nbtagscsvmcorr/I");  
	    tree_->Branch("t1metphicorr", 	&t1metphicorr_, 	"t1metphicorr/F");	      	      
	    tree_->Branch("t1metphicorrmt", 	&t1metphicorrmt_, 	"t1metphicorrmt/F");          
	    tree_->Branch("t1metphicorrphi", 	&t1metphicorrphi_, 	"t1metphicorrphi/F");        
	    tree_->Branch("t1metphicorrlep", 	&t1metphicorrlep_, 	"t1metphicorrlep/F");	      	      
	    tree_->Branch("t1metphicorrlepmt", 	&t1metphicorrlepmt_, 	"t1metphicorrlepmt/F");          
	    tree_->Branch("t1metphicorrlepphi", &t1metphicorrlepphi_, 	"t1metphicorrlepphi/F");        
	    tree_->Branch("pfcandpt10", 	&pfcandpt10_, 		"pfcandpt10/F");        
	    tree_->Branch("pfcandiso10", 	&pfcandiso10_, 		"pfcandiso10/F");      
	    tree_->Branch("trkpt10loose", 	&trkpt10loose_, 	"trkpt10loose/F");        
	    tree_->Branch("trkreliso10loose", 	&trkreliso10loose_, 	"trkreliso10loose/F");      
	    tree_->Branch("nleps", 		&nleps_, 		"nleps/I");
	    tree_->Branch("nmus",               &nmus_,                 "nmuus/I");
	    tree_->Branch("nels",               &nels_,                 "nels/I");
	    tree_->Branch("ntaus",              &ntaus_,                "ntaus/I");
	    

	    tree_->Branch("mcmtlns", &mcmtln_, "mcmtln/F");
	    tree_->Branch("nwzpartons", &nwzpartons_, "nwzpartons/F");
	    tree_->Branch("nbtagscsvl", &nbtagscsvl_, "nbtagscsvl/F");
	    tree_->Branch("trkmet_nolepcorr", &trkmet_nolepcorr_, "trkmet_nolepcorr/F");
	    tree_->Branch("trkmetphi_nolepcorr", &trkmetphi_nolepcorr_, "trkmetphi_nolepcorr/F");
	    tree_->Branch("trkmet", &trkmet_, "trkmet/F");
	    tree_->Branch("trkmetphi", &trkmetphi_, "trkmetphi/F");
	    tree_->Branch("isopf1", &isopf1_, "isopf1/F");
	    tree_->Branch("isopfold1", &isopfold1_, "isopfold1/F");
	    tree_->Branch("isopf2", &isopf2_, "isopf2/F");
	    tree_->Branch("eoverpin", &eoverpin_, "eoverpin/F");
	    tree_->Branch("eoverpout", &eoverpout_, "eoverpout/F");
	    tree_->Branch("dEtaIn", &dEtaIn_, "dEtaIn/F");
	    tree_->Branch("dPhiIn", &dPhiIn_, "dPhiIn/F");
	    tree_->Branch("sigmaIEtaIEta", &sigmaIEtaIEta_, "sigmaIEtaIEta/F");
	    tree_->Branch("hOverE", &hOverE_, "hOverE/F");
	    tree_->Branch("ooemoop", &ooemoop_, "ooemoop/F");
	    tree_->Branch("d0vtx", &d0vtx_, "d0vtx/F");
	    tree_->Branch("dzvtx", &dzvtx_, "dzvtx/F");
	    tree_->Branch("expinnerlayers", &expinnerlayers_, "expinnerlayers/F");
	    tree_->Branch("fbrem", &fbrem_, "fbrem/F");
	    tree_->Branch("pfisoch", &pfisoch_, "pfisoch/F");
	    tree_->Branch("pfisoem", &pfisoem_, "pfisoem/F");
	    tree_->Branch("pfisonh", &pfisonh_, "pfisonh/F");
	    tree_->Branch("eSC", &eSC_, "eSC/F");
	    tree_->Branch("phiSC", &phiSC_, "phiSC/F");
	    tree_->Branch("eSCRaw", &eSCRaw_, "eSCRaw/F");
	    tree_->Branch("eSCPresh", &eSCPresh_, "eSCPresh/F");
	    tree_->Branch("etasc1", &etasc1_, "etasc1/F");
	    tree_->Branch("iso1", &iso1_, "iso1/F");
	    tree_->Branch("isont1", &isont1_, "isont1/F");
	    tree_->Branch("iso2", &iso2_, "iso2/F");
	    tree_->Branch("isont2", &isont2_, "isont2/F");
	    
            tree_->Branch("lep1",    "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", &lep1Ptr_);
            tree_->Branch("lep2",    "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", &lep2Ptr_);
            tree_->Branch("t",       "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", &tPtr_);
            tree_->Branch("tbar",    "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", &tbarPtr_);
            tree_->Branch("stop_t",  "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", &stop_tPtr_);
            tree_->Branch("stop_tbar",  "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", &stop_tbarPtr_);
            tree_->Branch("lep_t",  "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", &lep_tPtr_);
            tree_->Branch("lep_tbar",  "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", &lep_tbarPtr_);
            tree_->Branch("mclep1",  "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", &mclep1Ptr_);
            tree_->Branch("mclep2",  "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", &mclep2Ptr_);
            tree_->Branch("pfcand10","ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", &pfcand10Ptr_);
	    /*
            tree_->Branch("pfjet1",  "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", &jet1Ptr_);
            tree_->Branch("pfjet2",  "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", &jet2Ptr_);
            tree_->Branch("pfjet3",  "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", &jet3Ptr_);
            tree_->Branch("pfjet4",  "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", &jet4Ptr_);
            tree_->Branch("pfjet5",  "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", &jet5Ptr_);
            tree_->Branch("pfjet6",  "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", &jet6Ptr_);
	    */
            tree_->Branch("pflep1",  "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", &pflep1Ptr_);
            tree_->Branch("pflep2",  "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", &pflep2Ptr_);

            tree_->Branch("candidates", "std::vector<Candidate>", &candidates_);
	    //            tree_->Branch("jets", "vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >", &jets_);
            tree_->Branch("pfjets"       , "vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >", &pfjets_);
            tree_->Branch("pfjets_csv"   , "std::vector<float>" , &pfjets_csv_Ptr_    );
            tree_->Branch("pfjets_qgtag" , "std::vector<float>" , &pfjets_qgtag_Ptr_  );
            tree_->Branch("pfjets_mc3"   , "std::vector<float>" , &pfjets_mc3_Ptr_    );
            tree_->Branch("pfjets_beta"  , "std::vector<float>" , &pfjets_beta_Ptr_   );
            tree_->Branch("pfjets_beta2" , "std::vector<float>" , &pfjets_beta2_Ptr_  );
            tree_->Branch("pfjets_lepjet", "std::vector<int>"   , &pfjets_lepjet_Ptr_ );

        }

        // initialze a StopTree
        void InitTree(){
            assert(tree_);
            // don't forget to set pointers to zero before you set address
            // or you will fully appreciate that "ROOT sucks" :)
            InitVariables();
            //Set branch address
            Int_t currentState = gErrorIgnoreLevel;
            // gErrorIgnoreLevel = kError;
            gErrorIgnoreLevel = kBreak;
	    tree_->SetBranchAddress("dataset", 		  &dataset_);	      	      
	    tree_->SetBranchAddress("run", 		  &run_);	      	      
	    tree_->SetBranchAddress("lumi", 		  &lumi_);	      	      
	    tree_->SetBranchAddress("event", 		  &event_);	      	      
	    tree_->SetBranchAddress("nvtx", 		  &nvtx_);	      	      
	    tree_->SetBranchAddress("nvtxweight", 	  &nvtxweight_);    
	    tree_->SetBranchAddress("weight", 		  &weight_);	      	      
	    tree_->SetBranchAddress("xsecsusy",		  &xsecsusy_);	      	      
	    tree_->SetBranchAddress("mutrigweight", 	  &mutrigweight_);    
	    tree_->SetBranchAddress("rhovor", 		  &rhovor_);	      	      
	    tree_->SetBranchAddress("mgcor", 		  &mgcor_);	      
	    tree_->SetBranchAddress("csc", 		  &csc_);	      	      
	    tree_->SetBranchAddress("hbhe", 		  &hbhe_);	      	      
	    tree_->SetBranchAddress("hcallaser", 	  &hcallaser_);	      	      
	    tree_->SetBranchAddress("ecaltp", 		  &ecaltp_);	      	      
	    tree_->SetBranchAddress("trkfail", 		  &trkfail_);	      	      
	    tree_->SetBranchAddress("eebadsc", 		  &eebadsc_);	      	      
	    tree_->SetBranchAddress("hbhenew", 		  &hbhenew_);	      	      	      
	    tree_->SetBranchAddress("isomu24", 		  &isomu24_);	      	      
	    tree_->SetBranchAddress("ele27wp80", 	  &ele27wp80_);	      	      
	    tree_->SetBranchAddress("trgmu1", 		  &trgmu1_);	      	      
	    tree_->SetBranchAddress("trgel1", 		  &trgel1_);	      	      
	    tree_->SetBranchAddress("mm", 		  &mm_);	      	      
	    tree_->SetBranchAddress("me", 		  &me_);	      	      
	    tree_->SetBranchAddress("em", 		  &em_);	      	      
	    tree_->SetBranchAddress("ee", 		  &ee_);	      	      
	    tree_->SetBranchAddress("ngoodlep", 	  &ngoodlep_);            
	    tree_->SetBranchAddress("leptype", 		  &leptype_);	      	      
	    tree_->SetBranchAddress("id1", 		  &id1_);	      	      
	    tree_->SetBranchAddress("id2", 		  &id2_);	      	      
	    tree_->SetBranchAddress("lep_t_id", 		  &lep_t_id_);	      	      
	    tree_->SetBranchAddress("lep_tbar_id", 		  &lep_tbar_id_);	      	      
	    tree_->SetBranchAddress("mcid1", 		  &mcid1_);	      	      
	    tree_->SetBranchAddress("mcid2", 		  &mcid2_);	      	      
	    tree_->SetBranchAddress("mcdecay2", 	  &mcdecay2_);	      	      
	    tree_->SetBranchAddress("mcndec2", 		  &mcndec2_);	      	      
	    tree_->SetBranchAddress("npfjets30", 	  &npfjets30_);          
	    tree_->SetBranchAddress("npfjets40", 	  &npfjets40_);          
	    tree_->SetBranchAddress("mctaudpt2", 	  &mctaudpt2_);	      	      
	    tree_->SetBranchAddress("lep1chi2ndf", 	  &lep1chi2ndf_);	      	      
	    tree_->SetBranchAddress("lep1dpt", 		  &lep1dpt_);	      	      
	    tree_->SetBranchAddress("dilmass", 		  &dilmass_);	      	      
	    tree_->SetBranchAddress("pfmet", 		  &pfmet_);	      	      
	    tree_->SetBranchAddress("mt", 	          &mt_);          
	    tree_->SetBranchAddress("pfmetphi", 	  &pfmetphi_);        
	    tree_->SetBranchAddress("t1met10", 		  &t1met10_);	      	      
	    tree_->SetBranchAddress("t1met10mt", 	  &t1met10mt_);          
	    tree_->SetBranchAddress("t1met10phi", 	  &t1met10phi_);        
	    tree_->SetBranchAddress("t1metphicorr", 	  &t1metphicorr_);	      	      
	    tree_->SetBranchAddress("t1metphicorrmt", 	  &t1metphicorrmt_);          
	    tree_->SetBranchAddress("t1metphicorrphi", 	  &t1metphicorrphi_);        
	    tree_->SetBranchAddress("t1metphicorrlep", 	  &t1metphicorrlep_);	      	      
	    tree_->SetBranchAddress("t1metphicorrlepmt",  &t1metphicorrlepmt_);          
	    tree_->SetBranchAddress("t1metphicorrlepphi", &t1metphicorrlepphi_);        
	    tree_->SetBranchAddress("nbtagscsvm", 	  &nbtagscsvm_);          
	    tree_->SetBranchAddress("nbtagscsvmcorr", 	  &nbtagscsvmcorr_);  
	    tree_->SetBranchAddress("pfcandpt10", 	  &pfcandpt10_);        
	    tree_->SetBranchAddress("pfcandiso10", 	  &pfcandiso10_);      
	    tree_->SetBranchAddress("trkpt10loose", 	  &trkpt10loose_);        
	    tree_->SetBranchAddress("trkreliso10loose",   &trkreliso10loose_);      
	    tree_->SetBranchAddress("nleps", 		  &nleps_);
	    tree_->SetBranchAddress("nmus",               &nmus_);       
	    tree_->SetBranchAddress("nels",               &nels_);       
	    tree_->SetBranchAddress("ntaus",              &ntaus_);       

	    
	    tree_->SetBranchAddress("mcmtlns",   &mcmtln_);
	    tree_->SetBranchAddress("nwzpartons",   &nwzpartons_);
	    tree_->SetBranchAddress("nbtagscsvl",   &nbtagscsvl_);
	    tree_->SetBranchAddress("trkmet_nolepcorr",   &trkmet_nolepcorr_);
	    tree_->SetBranchAddress("trkmetphi_nolepcorr",   &trkmetphi_nolepcorr_);
	    tree_->SetBranchAddress("trkmet",   &trkmet_);
	    tree_->SetBranchAddress("trkmetphi",   &trkmetphi_);
	    tree_->SetBranchAddress("isopf1",   &isopf1_);
	    tree_->SetBranchAddress("isopfold1",   &isopfold1_);
	    tree_->SetBranchAddress("isopf2",   &isopf2_);
	    tree_->SetBranchAddress("eoverpin",   &eoverpin_);
	    tree_->SetBranchAddress("eoverpout",   &eoverpout_);
	    tree_->SetBranchAddress("dEtaIn",   &dEtaIn_);
	    tree_->SetBranchAddress("dPhiIn",   &dPhiIn_);
	    tree_->SetBranchAddress("sigmaIEtaIEta",   &sigmaIEtaIEta_);
	    tree_->SetBranchAddress("hOverE",   &hOverE_);
	    tree_->SetBranchAddress("ooemoop",   &ooemoop_);
	    tree_->SetBranchAddress("d0vtx",   &d0vtx_);
	    tree_->SetBranchAddress("dzvtx",   &dzvtx_);
	    tree_->SetBranchAddress("expinnerlayers",   &expinnerlayers_);
	    tree_->SetBranchAddress("fbrem",   &fbrem_);
	    tree_->SetBranchAddress("pfisoch",   &pfisoch_);
	    tree_->SetBranchAddress("pfisoem",   &pfisoem_);
	    tree_->SetBranchAddress("pfisonh",   &pfisonh_);
	    tree_->SetBranchAddress("eSC",   &eSC_);
	    tree_->SetBranchAddress("phiSC",   &phiSC_);
	    tree_->SetBranchAddress("eSCRaw",   &eSCRaw_);
	    tree_->SetBranchAddress("eSCPresh",   &eSCPresh_);
	    tree_->SetBranchAddress("etasc1",   &etasc1_);
	    tree_->SetBranchAddress("iso1",   &iso1_);
	    tree_->SetBranchAddress("isont1",   &isont1_);
	    tree_->SetBranchAddress("iso2",   &iso2_);
	    tree_->SetBranchAddress("isont2",   &isont2_);

            tree_->SetBranchAddress("lep1",   		  &lep1Ptr_);
            tree_->SetBranchAddress("lep2",   		  &lep2Ptr_);
            tree_->SetBranchAddress("t",   		  &tPtr_);
            tree_->SetBranchAddress("tbar",   		  &tbarPtr_);
            tree_->SetBranchAddress("stop_t",   	  &stop_tPtr_);
            tree_->SetBranchAddress("stop_tbar",   	  &stop_tbarPtr_);
            tree_->SetBranchAddress("lep_t",   	  &lep_tPtr_);
            tree_->SetBranchAddress("lep_tbar",   	  &lep_tbarPtr_);
            tree_->SetBranchAddress("mclep1",   	  &mclep1Ptr_);
            tree_->SetBranchAddress("mclep2",   	  &mclep2Ptr_);
            tree_->SetBranchAddress("pfcand10",   	  &pfcand10Ptr_);
	    /*
            tree_->SetBranchAddress("pfjet1", 		  &jet1Ptr_);
            tree_->SetBranchAddress("pfjet2", 		  &jet2Ptr_);
            tree_->SetBranchAddress("pfjet3", 		  &jet3Ptr_);
            tree_->SetBranchAddress("pfjet4", 		  &jet4Ptr_);
            tree_->SetBranchAddress("pfjet5", 		  &jet5Ptr_);
            tree_->SetBranchAddress("pfjet6", 		  &jet6Ptr_);
	    */
            tree_->SetBranchAddress("pflep1",   	  &pflep1Ptr_);
            tree_->SetBranchAddress("pflep2",   	  &pflep2Ptr_);

            tree_->SetBranchAddress("candidates"    ,  &candidates_);
	    //            tree_->SetBranchAddress("jets",  &jets_);
            tree_->SetBranchAddress("pfjets"        ,  &pfjets_            );
            tree_->SetBranchAddress("pfjets_csv"    ,  &pfjets_csv_Ptr_    );
            tree_->SetBranchAddress("pfjets_qgtag"  ,  &pfjets_qgtag_Ptr_  );
            tree_->SetBranchAddress("pfjets_mc3"    ,  &pfjets_mc3_Ptr_    );
            tree_->SetBranchAddress("pfjets_beta"   ,  &pfjets_beta_Ptr_   );
            tree_->SetBranchAddress("pfjets_beta2"  ,  &pfjets_beta2_Ptr_  );
            tree_->SetBranchAddress("pfjets_lepjet" ,  &pfjets_lepjet_Ptr_ );

            gErrorIgnoreLevel = currentState;
        }

        /// get a built in type variable by name
        double Get(string value);
        /// compare two StopTrees for a given event on a given level of precision; 
        /// returns the variables that failed the comparison 
        vector<string> Compare(StopTree* value, double prec=0.005);

    private:

        LorentzVector *lep1Ptr_;
        LorentzVector *lep2Ptr_;
        LorentzVector *tPtr_;
        LorentzVector *tbarPtr_;
        LorentzVector *stop_tPtr_;
        LorentzVector *stop_tbarPtr_;
        LorentzVector *lep_tPtr_;
        LorentzVector *lep_tbarPtr_;
        LorentzVector *mclep1Ptr_;
        LorentzVector *mclep2Ptr_;
        LorentzVector *pfcand10Ptr_;
	/*
        LorentzVector *jet1Ptr_;
        LorentzVector *jet2Ptr_;
        LorentzVector *jet3Ptr_;
        LorentzVector *jet4Ptr_;
        LorentzVector *jet5Ptr_;
        LorentzVector *jet6Ptr_;
	*/
        LorentzVector *pflep1Ptr_;
        LorentzVector *pflep2Ptr_;

	vector<float>* pfjets_csv_Ptr_;
	vector<float>* pfjets_qgtag_Ptr_;
	vector<float>* pfjets_mc3_Ptr_;
	vector<float>* pfjets_beta_Ptr_;
	vector<float>* pfjets_beta2_Ptr_;
	vector<int>*   pfjets_lepjet_Ptr_;

}; 

inline void 
StopTree::InitVariables(){
    // create list of available variables
    if(variables_.empty()){
        //make sure that this is only done once
	variables_.push_back(string("run"		));
	variables_.push_back(string("dataset"		));
	variables_.push_back(string("lumi"		));
	variables_.push_back(string("event"		));
	variables_.push_back(string("nvtx"		));
	variables_.push_back(string("nvtxweight"	));
	variables_.push_back(string("weight"		));
	variables_.push_back(string("xsecsusy"		));
	variables_.push_back(string("mutrigweight"	));
	variables_.push_back(string("rhovor"		));
	variables_.push_back(string("mgcor"		));
	variables_.push_back(string("csc" 		));
	variables_.push_back(string("hbhe" 		));
	variables_.push_back(string("hcallaser" 	));
	variables_.push_back(string("ecaltp" 		));
	variables_.push_back(string("trkfail" 		));
	variables_.push_back(string("eebadsc" 		));
	variables_.push_back(string("hbhenew" 		));
	variables_.push_back(string("isomu24"		));
	variables_.push_back(string("ele27wp80"		));
	variables_.push_back(string("trgmu1"		));
	variables_.push_back(string("trgel1"		));
	variables_.push_back(string("mm"		));
	variables_.push_back(string("me"		));
	variables_.push_back(string("em"		));
	variables_.push_back(string("ee"		));
	variables_.push_back(string("ngoodlep"		));
	variables_.push_back(string("leptype"		));
	variables_.push_back(string("id1"		));
	variables_.push_back(string("id2"		));
	variables_.push_back(string("lep_t_id"		));
	variables_.push_back(string("lep_tbar_id"		));
	variables_.push_back(string("mcid1"		));
	variables_.push_back(string("mcid2"		));
	variables_.push_back(string("mcdecay2"		));
	variables_.push_back(string("mcndec2"		));
	variables_.push_back(string("npfjets30"		));
	variables_.push_back(string("npfjets40"		));
	variables_.push_back(string("mctaudpt2"		));
	variables_.push_back(string("lep1chi2ndf"	));
	variables_.push_back(string("lep1dpt"		));
	variables_.push_back(string("dilmass"		));
	variables_.push_back(string("pfmet"		));
	variables_.push_back(string("pfmt"		));
	variables_.push_back(string("pfmetphi"		));
	variables_.push_back(string("t1met10"		));
	variables_.push_back(string("t1met10mt"		));
	variables_.push_back(string("t1met10phi"	));
	variables_.push_back(string("t1metphicorr"	));
	variables_.push_back(string("t1metphicorrmt"	));
	variables_.push_back(string("t1metphicorrphi"	));
	variables_.push_back(string("t1metphicorrlep"	));
	variables_.push_back(string("t1metphicorrlepmt"	));
	variables_.push_back(string("t1metphicorrlepphi"));
	variables_.push_back(string("nbtagscsvm"	));
	variables_.push_back(string("nbtagscsvmcorr"	));
	variables_.push_back(string("pfcandpt10"	));
	variables_.push_back(string("pfcandiso10"	));
	variables_.push_back(string("trkpt10loose"	));
	variables_.push_back(string("trkreliso10loose"	));
	variables_.push_back(string("nleps"		));         
	variables_.push_back(string("nmus"		));         
	variables_.push_back(string("nels"		));         
	variables_.push_back(string("ntaus"		));         
	
	variables_.push_back(string("mcmtln"));
	variables_.push_back(string("nwzpartons"));
	variables_.push_back(string("nbtagscsvl"));
	variables_.push_back(string("trkmet_nolepcorr"));
	variables_.push_back(string("trkmetphi_nolepcorr"));
	variables_.push_back(string("trkmet"));
	variables_.push_back(string("trkmetphi"));
	variables_.push_back(string("isopf1"));
	variables_.push_back(string("isopfold1"));
	variables_.push_back(string("isopf2"));
	variables_.push_back(string("eoverpin"));
	variables_.push_back(string("eoverpout"));
	variables_.push_back(string("dEtaIn"));
	variables_.push_back(string("dPhiIn"));
	variables_.push_back(string("sigmaIEtaIEta"));
	variables_.push_back(string("hOverE"));
	variables_.push_back(string("ooemoop"));
	variables_.push_back(string("d0vtx"));
	variables_.push_back(string("dzvtx"));
	variables_.push_back(string("expinnerlayers"));
	variables_.push_back(string("fbrem"));
	variables_.push_back(string("pfisoch"));
	variables_.push_back(string("pfisoem"));
	variables_.push_back(string("pfisonh"));
	variables_.push_back(string("eSC"));
	variables_.push_back(string("phiSC"));
	variables_.push_back(string("eSCRaw"));
	variables_.push_back(string("eSCPresh"));
	variables_.push_back(string("etasc1"));
	variables_.push_back(string("iso1"));
	variables_.push_back(string("isont1"));
	variables_.push_back(string("iso2"));
	variables_.push_back(string("isont2"));

	variables_.push_back(string("lep1"		));
	variables_.push_back(string("lep2"		));
	variables_.push_back(string("t"			));
	variables_.push_back(string("tbar"		));
	variables_.push_back(string("stop_t"		));
	variables_.push_back(string("stop_tbar"		));
	variables_.push_back(string("lep_t"		));
	variables_.push_back(string("lep_tbar"		));
	variables_.push_back(string("mclep1"		));
	variables_.push_back(string("mclep2"		));
	variables_.push_back(string("pfcand10"		));
	variables_.push_back(string("pfjet1"		));
	variables_.push_back(string("pfjet2"		));
	variables_.push_back(string("pfjet3"		));
	variables_.push_back(string("pfjet4"		));
	variables_.push_back(string("pfjet5"		));
	variables_.push_back(string("pfjet6"		));
	variables_.push_back(string("pflep1"		));
	variables_.push_back(string("pflep2"		));
    }

    // inizialize variables
    //    dataset_		= "";
    run_		= 9;
    event_		= 9999;
    lumi_		= 9999;
    nvtx_		= 9999;
    nvtxweight_		= -999.;
    weight_		= -999.;
    xsecsusy_		= -999.;
    mutrigweight_	= -999.;
    rhovor_		= -999.;
    mgcor_		= -999.;
    csc_		= 999;
    hbhe_		= 999;
    hcallaser_		= 999;
    ecaltp_		= 999;
    trkfail_		= 999;
    eebadsc_		= 999;
    hbhenew_		= 999;
    isomu24_		= 999;
    ele27wp80_		= 999;
    trgmu1_		= 999;
    trgel1_		= 999;
    mm_			= 999;
    me_			= 999;
    em_			= 999;
    ee_			= 999;
    ngoodlep_		= 999;
    leptype_		= 999;
    id1_		= 999;
    id2_		= 999;
    lep_t_id_		= 999;
    lep_tbar_id_		= 999;
    mcid1_		= 999;
    mcid2_		= 999;
    mcdecay2_		= 999;
    mcndec2_		= 999;
    npfjets30_		= 999;
    npfjets40_		= 999;
    mctaudpt2_		= -999.;
    lep1chi2ndf_	= -999.;
    lep1dpt_		= -999.;
    dilmass_		= -999.;
    pfmet_		= -999.;
    mt_			= -999.;
    pfmetphi_		= -999.;
    t1met10_		= -999.;
    t1met10mt_		= -999.;
    t1met10phi_		= -999.;
    t1metphicorr_	= -999.;
    t1metphicorrmt_	= -999.;
    t1metphicorrphi_	= -999.;
    t1metphicorrlep_	= -999.;
    t1metphicorrlepmt_	= -999.;
    t1metphicorrlepphi_	= -999.;
    nbtagscsvm_		= 999;
    nbtagscsvmcorr_	= 999;
    pfcandpt10_		= -999.;
    pfcandiso10_	= -999.;
    trkpt10loose_	= -999.;
    trkreliso10loose_	= -999.;
    nleps_		= 999;
    nmus_		= 999;
    nels_		= 999;
    ntaus_		= 999;
    
    mcmtln_= -999;
    nwzpartons_= -999;
    nbtagscsvl_= -999;
    trkmet_nolepcorr_= -999.;
    trkmetphi_nolepcorr_= -999.;
    trkmet_= -999.;
    trkmetphi_= -999.;
    isopf1_= -999.;
    isopfold1_= -999.;
    isopf2_= -999.;
    eoverpin_= -999.;
    eoverpout_= -999.;
    dEtaIn_= -999.;
    dPhiIn_= -999.;
    sigmaIEtaIEta_= -999.;
    hOverE_= -999.;
    ooemoop_= -999.;
    d0vtx_= -999.;
    dzvtx_= -999.;
    expinnerlayers_= -999.;
    fbrem_= -999.;
    pfisoch_= -999.;
    pfisoem_= -999.;
    pfisonh_= -999.;
    eSC_= -999.;
    phiSC_= -999.;
    eSCRaw_= -999.;
    eSCPresh_= -999.;
    etasc1_= -999.;
    iso1_= -999.;
    isont1_= -999.;
    iso2_= -999.;
    isont2_= -999.;

    lep1_		= LorentzVector();
    lep2_		= LorentzVector();
    t_			= LorentzVector();
    tbar_		= LorentzVector();
    stop_t_		= LorentzVector();
    stop_tbar_		= LorentzVector();
    lep_t_		= LorentzVector();
    lep_tbar_		= LorentzVector();
    mclep1_		= LorentzVector();
    mclep2_		= LorentzVector();
    pfcand10_		= LorentzVector();

    pflep1_		= LorentzVector();
    pflep2_		= LorentzVector();

    candidates_ = 0;
    //    jets_ = 0;
    pfjets_ = 0;
    pfjets_csv_.clear();
    pfjets_qgtag_.clear();
    pfjets_mc3_.clear();
    pfjets_beta_.clear();
    pfjets_beta2_.clear();
    pfjets_lepjet_.clear();

}

inline double
StopTree::Get(string value)
{

  //  if(value=="dataset" 		) { return this->dataset_;		}	      	      
  if(value=="run" 		) { return this->run_;			}	      	      
  if(value=="lumi" 		) { return this->lumi_;			}	      	      
  if(value=="event" 		) { return this->event_;		}	      	      
  if(value=="nvtx" 		) { return this->nvtx_;			}	      	      
  if(value=="nvtxweight" 	) { return this->nvtxweight_;		}    
  if(value=="weight" 		) { return this->weight_;		}	      	      
  if(value=="xsecsusy" 		) { return this->xsecsusy_;		}	      	      
  if(value=="mutrigweight" 	) { return this->mutrigweight_;		}    
  if(value=="rhovor" 		) { return this->rhovor_;		}	      	      
  if(value=="mgcor" 		) { return this->mgcor_;		}
  if(value=="csc" 		) { return this->csc_;			}	      	      
  if(value=="hbhe" 		) { return this->hbhe_;			}	      	      
  if(value=="hcallaser" 	) { return this->hcallaser_;		}	      	      
  if(value=="ecaltp" 		) { return this->ecaltp_;		}	      	      
  if(value=="trkfail" 		) { return this->trkfail_;		}	      	      
  if(value=="eebadsc" 		) { return this->eebadsc_;		}	      	      
  if(value=="hbhenew" 		) { return this->hbhenew_;		}	      	      	      
  if(value=="isomu24" 		) { return this->isomu24_;		}	      	      
  if(value=="ele27wp80" 	) { return this->ele27wp80_;		}	      	      
  if(value=="trgmu1" 		) { return this->trgmu1_;		}	      	      
  if(value=="trgel1" 		) { return this->trgel1_;		}	      	      
  if(value=="mm" 		) { return this->mm_;			}	      	      
  if(value=="me" 		) { return this->me_;			}	      	      
  if(value=="em" 		) { return this->em_;			}	      	      
  if(value=="ee" 		) { return this->ee_;			}	      	      
  if(value=="ngoodlep" 		) { return this->ngoodlep_;		}            
  if(value=="leptype" 		) { return this->leptype_;		}	      	      
  if(value=="id1" 		) { return this->id1_;			}	      	      
  if(value=="id2" 		) { return this->id2_;			}	      	      
  if(value=="lep_t_id" 		) { return this->lep_t_id_;		}	      	      
  if(value=="lep_tbar_id" 	) { return this->lep_tbar_id_;		}	      	      
  if(value=="mcid1" 		) { return this->mcid1_;		}	      	      
  if(value=="mcid2" 		) { return this->mcid2_;		}	      	      
  if(value=="mcdecay2" 		) { return this->mcdecay2_;		}	      	      
  if(value=="mcndec2" 		) { return this->mcndec2_;		}	      	      
  if(value=="npfjets30" 	) { return this->npfjets30_;		}          
  if(value=="npfjets40" 	) { return this->npfjets40_;		}          
  if(value=="mctaudpt2" 	) { return this->mctaudpt2_;		}	      	      
  if(value=="lep1chi2ndf" 	) { return this->lep1chi2ndf_;		}	      	      
  if(value=="lep1dpt" 		) { return this->lep1dpt_;		}	      	      
  if(value=="dilmass" 		) { return this->dilmass_;		}	      	      
  if(value=="pfmet" 		) { return this->pfmet_;		}	      	      
  if(value=="mt" 		) { return this->mt_;			}          
  if(value=="pfmetphi" 		) { return this->pfmetphi_;		}        
  if(value=="t1met10" 		) { return this->t1met10_;		}	      	      
  if(value=="t1met10mt" 	) { return this->t1met10mt_;		}          
  if(value=="t1met10phi" 	) { return this->t1met10phi_;		}        
  if(value=="t1metphicorr" 	) { return this->t1metphicorr_;		}	      	      
  if(value=="t1metphicorrmt" 	) { return this->t1metphicorrmt_;	}          
  if(value=="t1metphicorrphi" 	) { return this->t1metphicorrphi_;	}        
  if(value=="t1metphicorrlep" 	) { return this->t1metphicorrlep_;	}	      	      
  if(value=="t1metphicorrlepmt" ) { return this->t1metphicorrlepmt_;	}          
  if(value=="t1metphicorrlepphi") { return this->t1metphicorrlepphi_;	}        
  if(value=="nbtagscsvm" 	) { return this->nbtagscsvm_;		}          
  if(value=="nbtagscsvmcorr" 	) { return this->nbtagscsvmcorr_;	}  
  if(value=="pfcandpt10" 	) { return this->pfcandpt10_;		}        
  if(value=="pfcandiso10" 	) { return this->pfcandiso10_;		}      
  if(value=="trkpt10loose" 	) { return this->trkpt10loose_;		}        
  if(value=="trkreliso10loose" 	) { return this->trkreliso10loose_;	}      
  if(value=="nleps" 		) { return this->nleps_;		}     
  if(value=="nmus"              ) { return this->nmus_;} 
  if(value=="nels"              ) { return this->nels_;} 
  if(value=="ntaus"             ) { return this->ntaus_;} 

  if(value=="mcmtln" ) { return this->mcmtln_; }
  if(value=="nwzpartons" ) { return this->nwzpartons_; }
  if(value=="nbtagscsvl" ) { return this->nbtagscsvl_; }
  if(value=="trkmet_nolepcorr" ) { return this->trkmet_nolepcorr_; }
  if(value=="trkmetphi_nolepcorr" ) { return this->trkmetphi_nolepcorr_; }
  if(value=="trkmet" ) { return this->trkmet_; }
  if(value=="trkmetphi" ) { return this->trkmetphi_; }
  if(value=="isopf1" ) { return this->isopf1_; }
  if(value=="isopfold1" ) { return this->isopfold1_; }
  if(value=="isopf2" ) { return this->isopf2_; }
  if(value=="eoverpin" ) { return this->eoverpin_; }
  if(value=="eoverpout" ) { return this->eoverpout_; }
  if(value=="dEtaIn" ) { return this->dEtaIn_; }
  if(value=="dPhiIn" ) { return this->dPhiIn_; }
  if(value=="sigmaIEtaIEta" ) { return this->sigmaIEtaIEta_; }
  if(value=="hOverE" ) { return this->hOverE_; }
  if(value=="ooemoop" ) { return this->ooemoop_; }
  if(value=="d0vtx" ) { return this->d0vtx_; }
  if(value=="dzvtx" ) { return this->dzvtx_; }
  if(value=="expinnerlayers" ) { return this->expinnerlayers_; }
  if(value=="fbrem" ) { return this->fbrem_; }
  if(value=="pfisoch" ) { return this->pfisoch_; }
  if(value=="pfisoem" ) { return this->pfisoem_; }
  if(value=="pfisonh" ) { return this->pfisonh_; }
  if(value=="eSC" ) { return this->eSC_; }
  if(value=="phiSC" ) { return this->phiSC_; }
  if(value=="eSCRaw" ) { return this->eSCRaw_; }
  if(value=="eSCPresh" ) { return this->eSCPresh_; }
  if(value=="etasc1" ) { return this->etasc1_; }
  if(value=="iso1" ) { return this->iso1_; }
  if(value=="isont1" ) { return this->isont1_; }
  if(value=="iso2" ) { return this->iso2_; }
  if(value=="isont2" ) { return this->isont2_; }
    
  /* if(value=="lep1"   		) { return this->lep1_;		} */
  /* if(value=="lep2"   		) { return this->lep2_;		} */
  /* if(value=="lep_tbar"   		) { return this->mclep2_;	} */
  /* if(value=="pfcand10"   		) { return this->pfcand10_;	} */
  /* if(value=="pfjet1" 		) { return this->pfjet1_;	} */
  /* if(value=="pfjet2" 		) { return this->pfjet2_;	} */
  /* if(value=="pfjet3" 		) { return this->pfjet3_;	} */
  /* if(value=="pfjet4" 		) { return this->pfjet4_;	} */
  /* if(value=="pfjet5" 		) { return this->pfjet5_;	} */
  /* if(value=="pfjet6" 		) { return this->pfjet6_;	} */

  return -9999.; 

}

inline vector<string> 
StopTree::Compare(StopTree* value, double prec){
    vector<string> fails;
    // this should always fit with ultimate precision
    /* if( this->event_ != value->event_ ){ fails.push_back( "event" ); } */
    /* if( this->run_   != value->run_   ){ fails.push_back( "run"   ); } */
    /* if( this->lumi_  != value->lumi_  ){ fails.push_back( "lumi"  ); } */

    // check within  (relative) precision
    for(vector<string>::const_iterator var=variables_.begin(); var!=variables_.end(); ++var){
        if( fabs(Get(*var)-value->Get(*var))/Get(*var)>prec ) fails.push_back(*var);
    }
    return fails;
}


#endif
