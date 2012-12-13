#include "../../CORE/jetSmearingTools.h"

#include "StopTreeLooper.h"
#include "../Plotting/PlotUtilities.h"
#include "../Core/MT2Utility.h"
#include "../Core/mt2bl_bisect.h"
#include "../Core/mt2w_bisect.h"
#include "PartonCombinatorics.h"

#include "../Core/stopUtils.h"

#include "TROOT.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TChainElement.h"
#include "TFile.h"
#include "TMath.h"
#include "TChain.h"
#include "Riostream.h"
#include "TFitter.h"

#include <algorithm>
#include <utility>
#include <map>
#include <set>
#include <list>

using namespace TMath;

//--------------------------------------------------------------------
float StopTreeLooper::vtxweight_n( const int nvertices, TH1F *hist, bool isData ) 
{

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

//--------------------------------------------------------------------
float StopTreeLooper::getdphi( float phi1 , float phi2 ) 
{
  float dphi = fabs( phi1 - phi2 );
  if( dphi > TMath::Pi() ) dphi = TMath::TwoPi() - dphi;
  return dphi;
}

//--------------------------------------------------------------------
float StopTreeLooper::getMT( float pt1 , float phi1 , float pt2 , float phi2 ) 
{

  float dphi = getdphi(phi1, phi2);
  return sqrt( 2 * ( pt1 * pt2 * (1 - cos( dphi ) ) ) );

}

//--------------------------------------------------------------------
float StopTreeLooper::dRbetweenVectors(LorentzVector vec1, 
				       LorentzVector vec2 )
{ 

  float dphi = std::min(::fabs(vec1.Phi() - vec2.Phi()), 2 * M_PI - fabs(vec1.Phi() - vec2.Phi()));
  float deta = vec1.Eta() - vec2.Eta();

  return sqrt(dphi*dphi + deta*deta);

}

//--------------------------------------------------------------------
StopTreeLooper::StopTreeLooper()
{
  m_outfilename_ = "histos.root";
  t1metphicorr = -9999.;
  t1metphicorrphi = -9999.;
  t1metphicorrmt = -9999.;
  min_mtpeak = -9999.;
  max_mtpeak = -9999.;
  n_jets = -999;
  jets.clear();
  btag.clear();
  mc.clear();
   
}

//--------------------------------------------------------------------
StopTreeLooper::~StopTreeLooper()
{
}

//--------------------------------------------------------------------
void StopTreeLooper::setOutFileName(string filename)
{
  m_outfilename_ = filename;

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
std::set<DorkyEventIdentifier> events_lasercalib; 
int load_badlaserevents  () {

  ifstream in;
  in.open("../Core/badlaser_events.txt");

  int run, event, lumi;
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
bool is_badLaserEvent (const DorkyEventIdentifier &id) {
  if (events_lasercalib.find(id) != events_lasercalib.end()) return true;
  return false;
}

//--------------------------------------------------------------------
float getDataMCRatio(float eta){
  if (eta >=0.0 && eta < 0.5) return 1.052;
  if (eta >=0.5 && eta < 1.1) return 1.057;
  if (eta >=1.1 && eta < 1.7) return 1.096;
  if (eta >=1.7 && eta < 2.3) return 1.134;
  if (eta >=2.3 && eta < 5.0) return 1.288;
  return 1.0;
}


void plotCandidate(StopTree* tree, MT2CHI2 mc, string tag, string sel, map<string,TH1F*> &h_1d , float evtweight){
	float chi_min = -1;
	float chi_max = 15;
	float mt_min = -10;
	float mt_max = 600;

	plot1D("h_"+tag+"_one_chi2"+sel, Min(mc.one_chi2, (float)(chi_max-0.01)) , evtweight , h_1d , 96  , chi_min ,  chi_max );

	plot1D("h_"+tag+"_two_mt2b"   +sel, Max(Min(mc.two_mt2b, (float)(mt_max-0.01)),(float)(mt_min+0.01)) , evtweight , h_1d ,110 , mt_min , mt_max );
	plot1D("h_"+tag+"_two_mt2bl"  +sel, Max(Min(mc.two_mt2bl, (float)(mt_max-0.01)),(float)(mt_min+0.01)) , evtweight , h_1d ,110 , mt_min , mt_max );
	plot1D("h_"+tag+"_two_mt2w"   +sel, Max(Min(mc.two_mt2w, (float)(mt_max-0.01)),(float)(mt_min+0.01)) , evtweight , h_1d ,110 , mt_min , mt_max );

	plot1D("h_"+tag+"_three_mt2b"   +sel, Max(Min(mc.three_mt2b, (float)(mt_max-0.01)),(float)(mt_min+0.01)) , evtweight , h_1d ,110 , mt_min , mt_max );
	plot1D("h_"+tag+"_three_mt2bl"   +sel, Max(Min(mc.three_mt2bl, (float)(mt_max-0.01)),(float)(mt_min+0.01)) , evtweight , h_1d ,110 , mt_min , mt_max );
	plot1D("h_"+tag+"_three_mt2w"   +sel, Max(Min(mc.three_mt2w, (float)(mt_max-0.01)),(float)(mt_min+0.01)) , evtweight , h_1d ,110 , mt_min , mt_max );

	plot1D("h_"+tag+"_four_chi2b"+sel, Min(mc.four_chi2b, (float)(chi_max-0.01)) , evtweight , h_1d , 96  , chi_min ,  chi_max );
	plot1D("h_"+tag+"_four_chi2bl"+sel, Min(mc.four_chi2bl, (float)(chi_max-0.01)) , evtweight , h_1d , 96  , chi_min ,  chi_max );
	plot1D("h_"+tag+"_four_chi2w"+sel, Min(mc.four_chi2w, (float)(chi_max-0.01)) , evtweight , h_1d , 96  , chi_min ,  chi_max );
}

void StopTreeLooper::loop(TChain *chain, TString name)
{

	printf("[StopTreeLooper::loop] %s\n", name.Data());

	load_badlaserevents  ();

	//
	// check for valid chain
	//

	TObjArray *listOfFiles = chain->GetListOfFiles();
	TIter fileIter(listOfFiles);
	if (listOfFiles->GetEntries() == 0) {
		cout << "[StopTreeLooper::loop] no files in chain" << endl;
		return;
	}

	//
	// set up histograms
	//

	gROOT->cd();

	cout << "[StopTreeLooper::loop] setting up histos" << endl;

	//plotting map
	std::map<std::string, TH1F*> h_1d;//h_cr1, h_cr4, h_cr5;
	std::map<std::string, TH2F*> h_2d;//h_cr1, h_cr4, h_cr5;

	// TFile* vtx_file = TFile::Open("vtxreweight/vtxreweight_Summer12_DR53X-PU_S10_9p7ifb_Zselection.root");
	// if( vtx_file == 0 ){
	//   cout << "vtxreweight error, couldn't open vtx file. Quitting!"<< endl;
	//   exit(0);
	// }

	// TH1F* h_vtx_wgt = (TH1F*)vtx_file->Get("hratio");
	// h_vtx_wgt->SetName("h_vtx_wgt");

	//
	// file loop
	//

	unsigned int nEventsChain=0;
	unsigned int nEvents = chain->GetEntries();
	nEventsChain = nEvents;
	unsigned int nEventsTotal = 0;
	int i_permille_old = 0;

	bool isData = name.Contains("data") ? true : false;

	//Define peak region
	min_mtpeak = 50.;
	max_mtpeak = 80.;
	printf("[StopTreeLooper::loop] MT PEAK definition %.0f - %.0f GeV \n", min_mtpeak, max_mtpeak);

	cout << "[StopTreeLooper::loop] running over chain with total entries " << nEvents << endl;

	while (TChainElement *currentFile = (TChainElement*)fileIter.Next()) {

		//
		// load the stop baby tree
		//

		StopTree *tree = new StopTree();
		tree->LoadTree(currentFile->GetTitle());
		tree->InitTree();

		//
		// event loop
		//

		ULong64_t nEvents = tree->tree_->GetEntries();

		//nEvents = 100;

		for(ULong64_t event = 0; event < nEvents; ++event) {
			tree->tree_->GetEntry(event);

			//
			// increment counters
			//

			++nEventsTotal;
			if (nEventsTotal%10000==0) {
				int i_permille = (int)floor(1000 * nEventsTotal / float(nEventsChain));
				//if (i_permille != i_permille_old) {//this prints too often!
				// xterm magic from L. Vacavant and A. Cerri
				if (isatty(1)) {
					printf("\015\033[32m ---> \033[1m\033[31m%4.1f%%"
							"\033[0m\033[32m <---\033[0m\015", i_permille/10.);
					fflush(stdout);
				}
				i_permille_old = i_permille;
			}

			//---------------------
			// skip duplicates
			//---------------------

			if( isData ) {
				DorkyEventIdentifier id = {tree->run_,tree->event_, tree->lumi_ };
				if (is_duplicate(id) ){
					continue;
				}
				if (is_badLaserEvent(id) ){
					//std::cout<<"Removed bad laser calibration event:" <<tree->run_<<"   "<<tree->event_<<"\n";
					continue;
				}
			}

			//
			// event weight
			//

			float evtweight    = isData ? 1. : ( tree->weight_ * 9.708 * tree->nvtxweight_ * tree->mgcor_ );
			float trigweight   = isData ? 1. : getsltrigweight(tree->id1_, tree->lep1_.Pt(), tree->lep1_.Eta());
			float trigweightdl = isData ? 1. : getdltrigweight(tree->id1_, tree->id2_);

			//
			// selection criteria
			//

			//preselection
			if ( !passEvtSelection(tree, name) ) continue;

			//ask for at least 4 jets
			if ( tree->npfjets30_ < 4 ) continue;

			//baseline met and mt requirements
			if (tree->t1metphicorr_  <100.) continue;
			string mtcut = tree->t1metphicorrmt_<120. ? "_lomt" : "_himt";

			// histogram tags
			//flavor types
			string flav_tag_sl;
			if ( abs(tree->id1_)==13 ) flav_tag_sl = "_muo";
			else if ( abs(tree->id1_)==11 ) flav_tag_sl = "_ele";
			else flav_tag_sl = "_mysterysl";
			string flav_tag_dl;
			if      ( abs(tree->id1_) == abs(tree->id2_) && abs(tree->id1_) == 13 ) flav_tag_dl = "_dimu";
			else if ( abs(tree->id1_) == abs(tree->id2_) && abs(tree->id1_) == 11 ) flav_tag_dl = "_diel";
			else if ( abs(tree->id1_) != abs(tree->id2_) && abs(tree->id1_) == 13 ) flav_tag_dl = "_muel";
			else if ( abs(tree->id1_) != abs(tree->id2_) && abs(tree->id1_) == 11 ) flav_tag_dl = "_elmu";
			else flav_tag_dl = "_mysterydl";
			string basic_flav_tag_dl = flav_tag_dl;
			if ( abs(tree->id1_) != abs(tree->id2_) ) basic_flav_tag_dl = "_mueg";

			//apply the cuts on the jets and check the b-tagging for the hadronic top reconstruction
			assert( tree->pfjets_->size() == tree->pfjets_csv_.size() );

			//
			// SIGNAL REGION - single lepton + b-tag
			//

			// selection - 1 lepton
			// Add iso track veto
			// Add b-tag
			if ( !passSingleLeptonSelection(tree, isData) ) continue;
			if ( !passIsoTrkVeto(tree) ) continue;

			n_jets = 0;
			jets.clear();
			btag.clear();
			mc.clear();

			int n_bm = 0;
			int n_bl = 0;

			for( unsigned int i = 0 ; i < tree->pfjets_->size() ; ++i ){

				if ( i > (N_JETS_TO_CONSIDER-1) ) break;
				if ( tree->pfjets_->at(i).Pt() < JET_PT ) continue;
				if ( fabs(tree->pfjets_->at(i).Eta()) > JET_ETA ) continue;

				//cout << i << " pt csv " << tree->pfjets_->at(i).pt() << " " << tree->pfjets_csv_.at(i) << endl;
				n_jets++;
				jets.push_back( tree->pfjets_->at(i)    );
				btag.push_back( tree->pfjets_csv_.at(i) );
				if ( !isData ) mc.push_back  ( tree->pfjets_mc3_.at(i) );

				if( tree->pfjets_csv_.at(i) > BTAG_MED ) n_bm++;
				if( tree->pfjets_csv_.at(i) > BTAG_LOW ) n_bl++;
			}

			static JetSmearer* jetSmearer = 0;
			if (jetSmearer == 0 ){
				std::vector<std::string> list_of_file_names;
				list_of_file_names.push_back("../../CORE/jetsmear/data/Spring10_PtResolution_AK5PF.txt");
				list_of_file_names.push_back("../../CORE/jetsmear/data/Spring10_PhiResolution_AK5PF.txt");
				list_of_file_names.push_back("../../CORE/jetsmear/data/jet_resolutions.txt");
				jetSmearer = makeJetSmearer(list_of_file_names);
			}

			LorentzVector* lep = &tree->lep1_;
			float met = tree->t1metphicorr_;
			float metphi = tree->t1metphicorrphi_;

                        vector<float> sigma_jets;
			for (int i=0; i<n_jets; ++i){
                                float sigma = getJetResolution(jets[i], jetSmearer);
                                if ( isData) sigma *= getDataMCRatio(jets[i].eta());
				sigma_jets.push_back(sigma);
                        }

			// get list of candidates
			PartonCombinatorics pc (jets, btag, sigma_jets, mc, *lep, met, metphi, isData);
			list<Candidate> allcandidates = pc.candidates_;


			assert( jets.size() == btag.size() );


			// jet info
			plot1D("h_njets" , Min(n_jets,9) , evtweight , h_1d ,10 , 0 , 10 );
			plot1D("h_bmjets", Min(n_bm,4)   , evtweight , h_1d , 5 , 0 ,  5 );
			plot1D("h_njets_minus_bjets" , Min(n_jets-n_bm,9) , evtweight , h_1d ,10 , 0 , 10 );
			plot1D("h_bljets", Min(n_bl,4)   , evtweight , h_1d , 5 , 0 ,  5 );

			plot2D("h_njets_vs_bmjets" , Min(n_jets,9), Min(n_bm,4) , evtweight , h_2d ,10 , 0 , 10 , 5 , 0 ,  5 );
			plot2D("h_njets_vs_bljets" , Min(n_jets,9), Min(n_bl,4) , evtweight , h_2d ,10 , 0 , 10 , 5 , 0 ,  5 );

			//-----------------------------------------------------------------
			// calculate the BEST MT2 variables only using Ricardo's algorithm
			//-----------------------------------------------------------------

			//histogram ranges
			int xnbins = 35;
			float xmin = -100.;
			float xmax = 600;
			int ynbins = 32;
			float ymin = -1.;
			float ymax = 15;



			plot1D("h_ncand", Min((int)allcandidates.size(),49) , evtweight , h_1d , 50, 0, 50);

			list<Candidate>::iterator candIter;
			for(candIter = allcandidates.begin() ; candIter != allcandidates.end() ; candIter++ ){

				float cand_mt2w = candIter->mt2w;
				if (cand_mt2w > xmax-0.0001) cand_mt2w = xmax-0.0001;
				if (cand_mt2w < xmin+0.0001) cand_mt2w = xmin+0.0001;
				float cand_chi2 = candIter->chi2;
				if (cand_chi2 > ymax-0.0001) cand_chi2 = ymax-0.0001;
				if (cand_chi2 < ymin+0.0001) cand_chi2 = ymin+0.0001;

				plot2D("h_xm",                 cand_mt2w, cand_chi2 ,     evtweight , h_2d , xnbins , xmin , xmax , ynbins , ymin , ymax);
				plot2D("h_xm" +mtcut,          cand_mt2w, cand_chi2 ,     evtweight , h_2d , xnbins , xmin , xmax , ynbins , ymin , ymax);
				if (candIter->match){
					plot2D("h_xm_match",                 cand_mt2w, cand_chi2 , evtweight ,     h_2d , xnbins , xmin , xmax , ynbins , ymin , ymax);
					plot2D("h_xm_match" +mtcut,          cand_mt2w, cand_chi2 , evtweight ,     h_2d , xnbins , xmin , xmax , ynbins , ymin , ymax);
				}
			}

			//require at least 1 btag
			if (tree->nbtagscsvm_<1) continue;
			list<Candidate> candidates = pc.b_candidates_;
			MT2CHI2 mc = pc.getMt2Chi2();

			plot1D("h_nbcand", Min((int)candidates.size(),49) , evtweight , h_1d , 50, 0, 50);

			list<Candidate>::iterator candIterB;
			for(candIterB = candidates.begin() ; candIterB != candidates.end() ; candIterB++ ){

				float cand_mt2w = candIterB->mt2w;
				if (cand_mt2w > xmax-0.0001) cand_mt2w = xmax-0.0001;
				if (cand_mt2w < xmin+0.0001) cand_mt2w = xmin+0.0001;
				float cand_chi2 = candIterB->chi2;
				if (cand_chi2 > ymax-0.0001) cand_chi2 = ymax-0.0001;
				if (cand_chi2 < ymin+0.0001) cand_chi2 = ymin+0.0001;

				plot2D("h_bcands_xm",                 cand_mt2w, cand_chi2 , evtweight ,     h_2d , xnbins , xmin , xmax , ynbins , ymin , ymax);
				plot2D("h_bcands_xm" +mtcut,          cand_mt2w, cand_chi2 , evtweight ,     h_2d , xnbins , xmin , xmax , ynbins , ymin , ymax);
				if (candIterB->match){
					plot2D("h_bcands_xm_match",                 cand_mt2w, cand_chi2 , evtweight ,     h_2d , xnbins , xmin , xmax , ynbins , ymin , ymax);
					plot2D("h_bcands_xm_match" +mtcut,          cand_mt2w, cand_chi2 , evtweight ,     h_2d , xnbins , xmin , xmax , ynbins , ymin , ymax);
				}
			}

			plot2D("h_one_three",                 mc.three_mt2w, mc.one_chi2 , evtweight ,     h_2d , xnbins , xmin , xmax , ynbins , ymin , ymax);
			plot2D("h_one_three" +mtcut,          mc.three_mt2w, mc.one_chi2 , evtweight ,     h_2d , xnbins , xmin , xmax , ynbins , ymin , ymax);

			plot2D("h_one_four",                  mc.one_chi2, mc.four_chi2w , evtweight ,     h_2d , ynbins , ymin , ymax , ynbins , ymin , ymax);
			plot2D("h_one_four" +mtcut,           mc.one_chi2, mc.four_chi2w , evtweight ,     h_2d , ynbins , ymin , ymax , ynbins , ymin , ymax);

			plot2D("h_two_three",                 mc.two_mt2w, mc.three_mt2w , evtweight ,     h_2d , xnbins , xmin , xmax , xnbins , xmin , xmax);
			plot2D("h_two_three" +mtcut,          mc.two_mt2w, mc.three_mt2w , evtweight ,     h_2d , xnbins , xmin , xmax , xnbins , xmin , xmax);


			plotCandidate(tree, mc , "sig" ,          "" , h_1d , evtweight*trigweight);
			plotCandidate(tree, mc , "sig" , flav_tag_sl , h_1d , evtweight*trigweight);
			plotCandidate(tree, mc , "sig" ,             mtcut , h_1d , evtweight*trigweight);
			plotCandidate(tree, mc , "sig" , flav_tag_sl+mtcut , h_1d , evtweight*trigweight);

		} // end event loop

		// delete tree;

	} // end file loop

	//
	// finish
	//

	char* outfilename = Form("output/%s.root",name.Data());
	savePlots(h_1d, outfilename);
	char* outfilename2 = Form("output/x_%s.root",name.Data());
	savePlots2(h_2d, outfilename2);

	already_seen.clear();

	gROOT->cd();

}

bool StopTreeLooper::passEvtSelection(const StopTree *sTree, TString name) 
{

  //rho requirement
  if ( sTree->rhovor_<0. || sTree->rhovor_>=40. ) return false;

  if (!name.Contains("T2")) {
    //met filters
    if ( sTree->csc_      != 0 ) return false;
    if ( sTree->hbhe_     != 1 ) return false;
    if ( sTree->hcallaser_!= 1 ) return false;
    if ( sTree->ecaltp_   != 1 ) return false;
    if ( sTree->trkfail_  != 1 ) return false;
    if ( sTree->eebadsc_  != 1 ) return false;
    if ( sTree->hbhenew_  != 1 ) return false;
  }

  //at least 1 lepton
  if ( sTree->ngoodlep_ < 1 ) return false;

  //if have more than 1 lepton, remove cases where have 2 close together
  if ( sTree->ngoodlep_ > 1 && 
       dRbetweenVectors( sTree->lep1_ ,  sTree->lep2_ )<0.1 ) return false;

  return true;

}

bool StopTreeLooper::passOneLeptonSelection(const StopTree *sTree, bool isData) 
{
  //single lepton selection for 8 TeV 53 analysis
  if ( !passSingleLeptonSelection(sTree, isData) ) return false;

  //pass isolated track veto
  //unfortunately changed default value to 9999.
  if ( sTree->pfcandpt10_ <9998. && sTree->pfcandiso10_ < 0.1 ) return false;

  return true;

}

bool StopTreeLooper::passTwoLeptonSelection(const StopTree *sTree, bool isData) 
{
  //single lepton selection for 8 TeV 53 analysis
  if ( !passDileptonSelection(sTree, isData) ) return false;

  //apply isolated track veto in addition to 2 leptons
  //default value for this one is -999
  if ( sTree->trkpt10loose_ >0. && sTree->trkreliso10loose_ < 0.1 ) return false;

  return true;

}

bool StopTreeLooper::passIsoTrkVeto(const StopTree *sTree) 
{

  //pass isolated track veto
  //unfortunately changed default value to 9999.
  if ( sTree->pfcandpt10_ <9998. && sTree->pfcandiso10_ < 0.1 ) return false;

  return true;

}

bool StopTreeLooper::passSingleLeptonSelection(const StopTree *sTree, bool isData) 
{
  //single lepton selection for 8 TeV 53 analysis

  //at least one lepton
  if ( sTree->ngoodlep_ < 1 ) return false;

  //lepton flavor - trigger, pt and eta requirements
  if ( sTree->lep1_.Pt() < 30 )          return false;
  if ( fabs( sTree->pflep1_.Pt() - sTree->lep1_.Pt() ) > 10. )  return false;
  if ( ( sTree->isopf1_ * sTree->lep1_.Pt() ) > 5. )  return false; 
  
  if ( sTree->leptype_ == 0 ) {

    //pass trigger if data - single electron
    if ( isData && sTree->ele27wp80_ != 1 ) return false;
    //    if ( isData && sTree->trgel1_ != 1 )  return false;
    
    //barrel only electrons
    if ( fabs(sTree->lep1_.Eta() ) > 1.4442) return false;
    if ( sTree->eoverpin_ > 4. ) return false;
    

  } else if ( sTree->leptype_ == 1 ) {

    //pass trigger if data - single muon
    if ( isData && sTree->isomu24_ != 1 ) return false;
    //    if ( isData && sTree->trgmu1_ != 1 )  return false;
    
    if ( fabs(sTree->lep1_.Eta() ) > 2.1)  return false;

  }

  return true;

}

bool StopTreeLooper::passDileptonSelection(const StopTree *sTree, bool isData) 
{
  //two lepton selection for 8 TeV 53 analysis

  //exactly 2 leptons
  if ( sTree->ngoodlep_ != 2 ) return false;

  //opposite sign
  if ( sTree->id1_*sTree->id2_>0 ) return false;

  //pass trigger if data - dilepton
  if ( isData && sTree->mm_ != 1 && sTree->me_ != 1 
       && sTree->em_ != 1 && sTree->ee_ != 1 ) return false;

  //passes pt and eta requirements
  if ( sTree->lep1_.Pt() < 20 )          return false;
  if ( sTree->lep2_.Pt() < 20 )          return false;
  if ( fabs(sTree->lep1_.Eta() ) > 2.4)  return false;
  if ( fabs(sTree->lep2_.Eta() ) > 2.4)  return false;

  //consistency with pf leptons
  if ( fabs( sTree->pflep1_.Pt() - sTree->lep1_.Pt() ) > 10. )  return false;
  if ( fabs( sTree->pflep2_.Pt() - sTree->lep2_.Pt() ) > 10. )  return false;

  //information is only stored for leading lepton
  if ( ( sTree->isopf1_ * sTree->lep1_.Pt() ) > 5. )  return false; 
  if ( fabs(sTree->id1_)==11 && sTree->eoverpin_ > 4. ) return false;

  //barrel only electrons
  if (fabs(sTree->id1_)==11 && fabs(sTree->lep1_.Eta() ) > 1.4442) return false;
  if (fabs(sTree->id2_)==11 && fabs(sTree->lep2_.Eta() ) > 1.4442) return false;
  
  return true;

}

bool StopTreeLooper::passLepPlusIsoTrkSelection(const StopTree *sTree, bool isData) 
{
  //single lepton plus iso trk selection for 8 TeV 53 analysis

  //at least one lepton
  if ( !passSingleLeptonSelection(sTree, isData) ) return false;

  //pass isolated track requirement
  //unfortunately changed default value to 9999.
  if ( sTree->pfcandpt10_ > 9990. || sTree->pfcandiso10_ > 0.1 ) return false;

  return true;

}

pair<float,float> StopTreeLooper::getPhiCorrMET( float met, float metphi, int nvtx, bool ismc){

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

float StopTreeLooper::getdltrigweight(int id1, int id2)
{ 
  if (abs(id1)==11 && abs(id2)==11) return 0.95;
  if (abs(id1)==13 && abs(id2)==13) return 0.88;
  if (abs(id1)!=abs(id2)) return 0.92;
  return -999.;

}

float StopTreeLooper::getsltrigweight(int id1, float pt, float eta) {

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

