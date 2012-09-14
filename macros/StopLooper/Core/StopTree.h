#ifndef StopTree_H
#define StopTree_H

#include "TFile.h"
#include "TTree.h"
#include "TError.h"

#include "Math/LorentzVector.h"

#include <cmath>
#include "assert.h"

using namespace std;
using namespace ROOT::Math;

//
// Ntuple structure:
//
// Ntuple content:

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
	int           npfjets30_;
	int           npfjets40_;
	float 	      lep1chi2ndf_;	
	float 	      lep1dpt_;	
	float         dilmass_;
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
	int           nleps_;         
	LorentzVector lep1_;
	LorentzVector lep2_;
	LorentzVector pfcand10_;
	LorentzVector pfjet1_;
	LorentzVector pfjet2_;
	LorentzVector pfjet3_;
	LorentzVector pfjet4_;
	LorentzVector pfjet5_;
	LorentzVector pfjet6_;

    public:
        /// this is the main element
        TTree *tree_;
        TFile *f_;

        /// hold the names of variables to facilitate things (filled during Init)
        vector<string> variables_;
	
        /// default constructor  
 StopTree() :  lep1Ptr_(&lep1_), lep2Ptr_(&lep2_), pfcand10Ptr_(&pfcand10_), jet1Ptr_(&pfjet1_), jet2Ptr_(&pfjet2_), jet3Ptr_(&pfjet3_), jet4Ptr_(&pfjet4_), jet5Ptr_(&pfjet5_), jet6Ptr_(&pfjet6_) {}
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
	    tree_->Branch("npfjets30", 		&npfjets30_, 		"npfjets30/I");        
	    tree_->Branch("npfjets40", 		&npfjets40_, 		"npfjets40/I");            
	    tree_->Branch("lep1chi2ndf", 	&lep1chi2ndf_, 		"lep1chi2ndf/F");	      	      
	    tree_->Branch("lep1dpt", 		&lep1dpt_, 		"lep1dpt/F");
	    tree_->Branch("dilmmas", 		&dilmass_, 		"dilmass/F");	      	      
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
	    tree_->Branch("nleps", 		&nleps_, 		"nleps/I");     
            tree_->Branch("lep1",    "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", &lep1Ptr_);
            tree_->Branch("lep2",    "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", &lep2Ptr_);
            tree_->Branch("pfcand10","ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", &pfcand10Ptr_);
            tree_->Branch("pfjet1",  "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", &jet1Ptr_);
            tree_->Branch("pfjet2",  "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", &jet2Ptr_);
            tree_->Branch("pfjet3",  "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", &jet3Ptr_);
            tree_->Branch("pfjet4",  "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", &jet4Ptr_);
            tree_->Branch("pfjet5",  "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", &jet5Ptr_);
            tree_->Branch("pfjet6",  "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", &jet6Ptr_);
             
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
	    tree_->SetBranchAddress("npfjets30", 	  &npfjets30_);          
	    tree_->SetBranchAddress("npfjets40", 	  &npfjets40_);          
	    tree_->SetBranchAddress("lep1chi2ndf", 	  &lep1chi2ndf_);	      	      
	    tree_->SetBranchAddress("lep1dpt", 		  &lep1dpt_);	      	      
	    tree_->SetBranchAddress("dilmass", 		  &dilmass_);	      	      
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
	    tree_->SetBranchAddress("nleps", 		  &nleps_);     
            tree_->SetBranchAddress("lep1",   		  &lep1Ptr_);
            tree_->SetBranchAddress("lep2",   		  &lep2Ptr_);
            tree_->SetBranchAddress("pfcand10",   	  &pfcand10Ptr_);
            tree_->SetBranchAddress("pfjet1", 		  &jet1Ptr_);
            tree_->SetBranchAddress("pfjet2", 		  &jet2Ptr_);
            tree_->SetBranchAddress("pfjet3", 		  &jet3Ptr_);
            tree_->SetBranchAddress("pfjet4", 		  &jet4Ptr_);
            tree_->SetBranchAddress("pfjet5", 		  &jet5Ptr_);
            tree_->SetBranchAddress("pfjet6", 		  &jet6Ptr_);

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
        LorentzVector *pfcand10Ptr_;
        LorentzVector *jet1Ptr_;
        LorentzVector *jet2Ptr_;
        LorentzVector *jet3Ptr_;
        LorentzVector *jet4Ptr_;
        LorentzVector *jet5Ptr_;
        LorentzVector *jet6Ptr_;

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
	variables_.push_back(string("npfjets30"		));
	variables_.push_back(string("npfjets40"		));
	variables_.push_back(string("lep1chi2ndf"	));
	variables_.push_back(string("lep1dpt"		));
	variables_.push_back(string("dilmass"		));
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
	variables_.push_back(string("nleps"		));         
	variables_.push_back(string("lep1"		));
	variables_.push_back(string("lep2"		));
	variables_.push_back(string("pfcand10"		));
	variables_.push_back(string("pfjet1"		));
	variables_.push_back(string("pfjet2"		));
	variables_.push_back(string("pfjet3"		));
	variables_.push_back(string("pfjet4"		));
	variables_.push_back(string("pfjet5"		));
	variables_.push_back(string("pfjet6"		));
    }

    // inizialize variables
    //    dataset_		= "";
    run_		= 9;
    event_		= 9999;
    lumi_		= 9999;
    nvtx_		= 9999;
    nvtxweight_		= -999.;
    weight_		= -999.;
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
    npfjets30_		= 999;
    npfjets40_		= 999;
    lep1chi2ndf_	= -999.;
    lep1dpt_		= -999.;
    dilmass_		= -999.;
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
    nleps_		= 999;         
    lep1_		= LorentzVector();
    lep2_		= LorentzVector();
    pfcand10_		= LorentzVector();
    pfjet1_		= LorentzVector();
    pfjet2_		= LorentzVector();
    pfjet3_		= LorentzVector();
    pfjet4_		= LorentzVector();
    pfjet5_		= LorentzVector();
    pfjet6_		= LorentzVector();
    
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
  if(value=="npfjets30" 	) { return this->npfjets30_;		}          
  if(value=="npfjets40" 	) { return this->npfjets40_;		}          
  if(value=="lep1chi2ndf" 	) { return this->lep1chi2ndf_;		}	      	      
  if(value=="lep1dpt" 		) { return this->lep1dpt_;		}	      	      
  if(value=="dilmass" 		) { return this->dilmass_;		}	      	      
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
  if(value=="nleps" 		) { return this->nleps_;		}     
  /* if(value=="lep1"   		) { return this->lep1_;		} */
  /* if(value=="lep2"   		) { return this->lep2_;		} */
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
