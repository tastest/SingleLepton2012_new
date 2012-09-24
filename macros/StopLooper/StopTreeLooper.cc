
#include "StopTreeLooper.h"
#include "Core/StopTree.h"
#include "Plotting/PlotUtilities.h"

#include "TROOT.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TChainElement.h"
#include "TFile.h"
#include "TMath.h"
#include "TChain.h"

#include <algorithm>
#include <utility>
#include <map>

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

float StopTreeLooper::getdphi( float phi1 , float phi2 ) 
{
  float dphi = fabs( phi1 - phi2 );
  if( dphi > TMath::Pi() ) dphi = TMath::TwoPi() - dphi;
  return dphi;
}

float StopTreeLooper::getMT( float pt1 , float phi1 , float pt2 , float phi2 ) 
{

  float dphi = getdphi(phi1, phi2);
  return sqrt( 2 * ( pt1 * pt2 * (1 - cos( dphi ) ) ) );

}

float StopTreeLooper::dRbetweenVectors(LorentzVector vec1, 
				       LorentzVector vec2 )
{ 

  float dphi = std::min(::fabs(vec1.Phi() - vec2.Phi()), 2 * M_PI - fabs(vec1.Phi() - vec2.Phi()));
  float deta = vec1.Eta() - vec2.Eta();

  return sqrt(dphi*dphi + deta*deta);

}


StopTreeLooper::StopTreeLooper()
{
  m_outfilename_ = "histos.root";
}

StopTreeLooper::~StopTreeLooper()
{
}

void StopTreeLooper::setOutFileName(string filename)
{
  m_outfilename_ = filename;

}

void StopTreeLooper::loop(TChain *chain, TString name)
{

  printf("[StopTreeLooper::loop] %s\n", name.Data());

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
  std::map<std::string, TH1F*> h_1d;
  //also for control regions
  std::map<std::string, TH1F*> h_1d_cr1, h_1d_cr2, h_1d_cr4, h_1d_cr5;
  //for ttbar dilepton njets distribution
  std::map<std::string, TH1F*> h_1d_nj;
  //z sample for yields etc
  std::map<std::string, TH1F*> h_1d_z;

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

      // 
      // event weight
      // 

      float evtweight = isData ? 1. : ( tree->weight_ * 9.708 * tree->nvtxweight_ * tree->mgcor_ );
      // to reweight from file - also need to comment stuff before
      //      float vtxweight = vtxweight_n( tree->nvtx_, h_vtx_wgt, isData );

      plot1D("h_vtx",       tree->nvtx_,       evtweight, h_1d_z, 40, 0, 40);
      plot1D("h_vtxweight", tree->nvtxweight_, evtweight, h_1d_z, 41, -4., 4.);

      // 
      // selection criteria
      // 

	//cout<<"got to line "<<__LINE__<<endl;
      //preselection
      if ( !passEvtSelection(tree) ) continue;
	//cout<<"got to line "<<__LINE__<<endl;

      // histogram tags
      //jet multiplicity
      string tag_njets = Form("_nj%i", (tree->npfjets30_<4) ? tree->npfjets30_ : 4);
      //b-tagging
      string tag_btag = (tree->nbtagscsvm_<1) ? "_bveto" : "_btag";
      //z-peak/veto
      string tag_zcut;
      if ( fabs( tree->dilmass_ - 91.) > 15. ) tag_zcut = "_zveto";
      else if  ( fabs( tree->dilmass_ - 91.) < 10. ) tag_zcut = "_zpeak";
      else tag_zcut = "_ignore";

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

      //
      // CR1 - single lepton + b-veto
      //

      // selection - 1 lepton + iso track veto
      // Add b-tag veto
      if ( passOneLeptonSelection(tree) )
	{
	//cout<<"got to line "<<__LINE__<<endl;

	  //default 
	  makeCR1Plots( tree, evtweight, h_1d_cr1, "_prebveto", tag_njets, flav_tag_sl, 150. );
	  //met > 50 GeV requirement 
	  if ( tree->t1met10_ > 50. ) 
	    makeCR1Plots( tree, evtweight, h_1d_cr1, "_prebveto_met50", tag_njets, flav_tag_sl, 150. );
	  //met > 100 GeV requirement 
	  if ( tree->t1met10_ > 100. ) 
	    makeCR1Plots( tree, evtweight, h_1d_cr1, "_prebveto_met100", tag_njets, flav_tag_sl, 150. );
	  //met > 150 GeV requirement 
	  if ( tree->t1met10_ > 150. ) 
	    makeCR1Plots( tree, evtweight, h_1d_cr1, "_prebveto_met150", tag_njets, flav_tag_sl, 120. );
	  //met > 200 GeV requirement 
	  if ( tree->t1met10_ > 200. ) 
	    makeCR1Plots( tree, evtweight, h_1d_cr1, "_prebveto_met200", tag_njets, flav_tag_sl, 120. );
	  //met > 250 GeV requirement 
	  if ( tree->t1met10_ > 250. ) 
	    makeCR1Plots( tree, evtweight, h_1d_cr1, "_prebveto_met250", tag_njets, flav_tag_sl, 120. );
	  
	  if ( tree->nbtagscsvm_==0 ) {
	//cout<<"got to line "<<__LINE__<<endl;

	     // if ( tree->t1met10mt_>250 ) 
		//cout<<"MT: "<<tree->t1met10mt_<<" dPhi(lep,MET): "<<getdphi( tree->lep1_.Phi() , tree->t1met10phi_ )<<" pT(lep): "<<tree->lep1_.Pt() <<" MET: "<<tree->t1met10_ <<" * dataset: "<<tree->dataset_
		  //  <<" run: "<<tree->run_<<" lumi: "<<tree->lumi_<<" event: "<<tree->event_<<endl;


	    //default 
	    makeCR1Plots( tree, evtweight, h_1d_cr1, "", tag_njets, flav_tag_sl, 150. );
	    //met > 50 GeV requirement 
	    if ( tree->t1met10_ > 50. ) 
	      makeCR1Plots( tree, evtweight, h_1d_cr1, "_met50", tag_njets, flav_tag_sl, 150. );
	    //met > 100 GeV requirement 
	    if ( tree->t1met10_ > 100. ) 
	      makeCR1Plots( tree, evtweight, h_1d_cr1, "_met100", tag_njets, flav_tag_sl, 150. );
	    //met > 150 GeV requirement 
	    if ( tree->t1met10_ > 150. ) 
	      makeCR1Plots( tree, evtweight, h_1d_cr1, "_met150", tag_njets, flav_tag_sl, 120. );
	    //met > 200 GeV requirement 
	    if ( tree->t1met10_ > 200. ) 
	      makeCR1Plots( tree, evtweight, h_1d_cr1, "_met200", tag_njets, flav_tag_sl, 120. );
	    //met > 250 GeV requirement 
	    if ( tree->t1met10_ > 250. ) 
	      makeCR1Plots( tree, evtweight, h_1d_cr1, "_met250", tag_njets, flav_tag_sl, 120. );
	    
	  }
	}

      //
      // CR2 - Z-peak for yields and mT resolution studies
      // 

      // selection - SF dilepton, in z-peak
      if ( passDileptonSelection(tree) && abs(tree->id1_) == abs(tree->id2_) )
	{
	  //invariant mass
	  plot1D("h_z_dilmass"          +flav_tag_dl, tree->dilmass_, evtweight, h_1d_z,  30 , 76 , 106);
	  plot1D("h_z_dilmass"+tag_njets+flav_tag_dl, tree->dilmass_, evtweight, h_1d_z,  30 , 76 , 106);
	  
	  if ( fabs( tree->dilmass_ - 91.) < 10. ) 
	    {

	      //if (tree->npfjets30_>8) 
		//cout<<"NJETS: "<<tree->npfjets30_<<" * dataset: "<<tree->dataset_
		 //   <<" run: "<<tree->run_<<" lumi: "<<tree->lumi_<<" event: "<<tree->event_<<endl;
	      
	      //z peak plots
	      plot1D("h_z_njets"    +flav_tag_dl, min(tree->npfjets30_,4),  evtweight, h_1d_z, 5,0,5);
	      plot1D("h_z_njets_all"+flav_tag_dl, min(tree->npfjets30_,9),  evtweight, h_1d_z, 10, 0, 10);
	      plot1D("h_z_nbjets"   +flav_tag_dl, min(tree->nbtagscsvm_,3), evtweight, h_1d_z, 4, 0, 4);
	      makeZPlots( tree, evtweight, h_1d_z, "", tag_njets, flav_tag_dl );

	      // Add b-tag veto 
	      if ( tree->nbtagscsvm_==0 && tree->npfjets30_ >= 2 ) {

		//find positive lepton
		bool isfirstp = (tree->id1_ < 0) ? true : false;
		
		//recalculate met
		float metx = tree->t1metphicorr_ * cos( tree->t1metphicorrphi_ );
		float mety = tree->t1metphicorr_ * sin( tree->t1metphicorrphi_ );
		
		//recalculate the MET with the positive lepton
		metx += isfirstp ? tree->lep1_.px() : tree->lep2_.px();
		mety += isfirstp ? tree->lep1_.py() : tree->lep2_.py();
		
		float t1metphicorr_lep    = sqrt(metx*metx + mety*mety);
		
		//default 
		makeCR2Plots( tree, evtweight, h_1d_cr2, "", tag_njets, flav_tag_dl, 150. );
		//pseudomet > 50 GeV requirement 
		if ( t1metphicorr_lep > 50. ) 
		  makeCR2Plots( tree, evtweight, h_1d_cr2, "_met50", tag_njets, flav_tag_dl, 150. );
		//pseudomet > 100 GeV requirement 
		if ( t1metphicorr_lep > 100. ) 
		  makeCR2Plots( tree, evtweight, h_1d_cr2, "_met100", tag_njets, flav_tag_dl, 150. );
		//pseudomet > 150 GeV requirement 
		if ( t1metphicorr_lep > 150. ) 
		  makeCR2Plots( tree, evtweight, h_1d_cr2, "_met150", tag_njets, flav_tag_dl, 120. );
		//pseudomet > 200 GeV requirement 
		if ( t1metphicorr_lep > 200. ) 
		  makeCR2Plots( tree, evtweight, h_1d_cr2, "_met200", tag_njets, flav_tag_dl, 120. );
		//pseudomet > 250 GeV requirement 
		if ( t1metphicorr_lep > 250. ) 
		  makeCR2Plots( tree, evtweight, h_1d_cr2, "_met250", tag_njets, flav_tag_dl, 120. );
	      }
	    }
	}

      //
      // CR4 - ttbar dilepton sample with 2 good leptons
      //
      
      // selection - all dilepton, z-veto for SF dilepton
      // Add b-tag requirement
      if ( passDileptonSelection(tree) 
	   && (abs(tree->id1_) != abs(tree->id2_) || fabs( tree->dilmass_ - 91.) > 15. ) 
	   && tree->nbtagscsvm_>0 ) 
	{

	  //jet multiplicity distributions 
	  //store in separate file since this is used for njet reweighting
	  makeNJPlots( tree, evtweight, h_1d_nj, "", basic_flav_tag_dl);
	  //met > 50 GeV requirement 
	  if ( tree->t1met10_ > 50. ) 
	    makeNJPlots( tree, evtweight, h_1d_nj, "_met50", basic_flav_tag_dl);
	  //met > 100 GeV requirement 
	  if ( tree->t1met10_ > 100. ) 
	    makeNJPlots( tree, evtweight, h_1d_nj, "_met100", basic_flav_tag_dl);
	  //met > 150 GeV requirement 
	  if ( tree->t1met10_ > 150. ) 
	    makeNJPlots( tree, evtweight, h_1d_nj, "_met150", basic_flav_tag_dl);	    
	  //met > 200 GeV requirement 
	  if ( tree->t1met10_ > 200. ) 
	    makeNJPlots( tree, evtweight, h_1d_nj, "_met200", basic_flav_tag_dl);	    
	  //met > 250 GeV requirement 
	  if ( tree->t1met10_ > 250. ) 
	    makeNJPlots( tree, evtweight, h_1d_nj, "_met250", basic_flav_tag_dl);	    

	  if ( tree->npfjets30_ < 2 ) continue; 
	  
	  //default 
	  makeCR4Plots( tree, evtweight, h_1d_cr4, "", tag_njets, flav_tag_dl, 150. );
	  //met > 50 GeV requirement 
	  if ( tree->t1met10_ > 50. ) 
	    makeCR4Plots( tree, evtweight, h_1d_cr4, "_met50", tag_njets, flav_tag_dl, 150. );
	  //met > 100 GeV requirement 
	  if ( tree->t1met10_ > 100. ) 
	    makeCR4Plots( tree, evtweight, h_1d_cr4, "_met100", tag_njets, flav_tag_dl, 150. );
	  //met > 150 GeV requirement 
	  if ( tree->t1met10_ > 150. ) 
	    makeCR4Plots( tree, evtweight, h_1d_cr4, "_met150", tag_njets, flav_tag_dl, 120. );
	  //met > 200 GeV requirement 
	  if ( tree->t1met10_ > 200. ) 
	    makeCR4Plots( tree, evtweight, h_1d_cr4, "_met200", tag_njets, flav_tag_dl, 120. );
	  //met > 250 GeV requirement 
	  if ( tree->t1met10_ > 250. ) 
	    makeCR4Plots( tree, evtweight, h_1d_cr4, "_met250", tag_njets, flav_tag_dl, 120. );
	  
	}

      ////////////////////////////////////////////////////////////////////////////////////////////////////
      // Ask for at least 2 jets from now on
      if ( tree->npfjets30_ < 2 ) continue;


      //
      // Sample before isolated track requirement - for fake rate of requirement
      //

      // selection - at least 1 lepton
      // Add b-tag requirement
      if ( passSingleLeptonSelection(tree) 
	   && tree->nbtagscsvm_>0 ) 
	{
	  
	  //inclusive sample
	  //default 
	  makeCR1Plots( tree, evtweight, h_1d_cr5, "_preveto", tag_njets, flav_tag_sl, 150. );
	  //met > 50 GeV requirement 
	  if ( tree->t1met10_ > 50. ) 
	    makeCR1Plots( tree, evtweight, h_1d_cr5, "_preveto_met50", tag_njets, flav_tag_sl, 150. );
	  //met > 100 GeV requirement 
	  if ( tree->t1met10_ > 100. ) 
	    makeCR1Plots( tree, evtweight, h_1d_cr5, "_preveto_met100", tag_njets, flav_tag_sl, 150. );
	  //met > 150 GeV requirement 
	  if ( tree->t1met10_ > 150. ) 
	    makeCR1Plots( tree, evtweight, h_1d_cr5, "_preveto_met150", tag_njets, flav_tag_sl, 120. );
	  //met > 200 GeV requirement 
	  if ( tree->t1met10_ > 200. ) 
	    makeCR1Plots( tree, evtweight, h_1d_cr5, "_preveto_met200", tag_njets, flav_tag_sl, 120. );
	  //met > 250 GeV requirement 
	  if ( tree->t1met10_ > 250. ) 
	    makeCR1Plots( tree, evtweight, h_1d_cr5, "_preveto_met250", tag_njets, flav_tag_sl, 120. );
	}
      
      //
      // CR5 - lepton + isolated track
      //

      // selection - lepton + isolated track
      // Add b-tag requirement
      if ( passLepPlusIsoTrkSelection(tree) 
	   && tree->nbtagscsvm_>0 ) 
	{

	  //inclusive sample
	  //default 
	  makeCR5Plots( tree, evtweight, h_1d_cr5, "_all", tag_njets, flav_tag_sl, 150. );
	  //met > 50 GeV requirement 
	  if ( tree->t1met10_ > 50. ) 
	    makeCR5Plots( tree, evtweight, h_1d_cr5, "_all_met50", tag_njets, flav_tag_sl, 150. );
	  //met > 100 GeV requirement 
	  if ( tree->t1met10_ > 100. ) 
	    makeCR5Plots( tree, evtweight, h_1d_cr5, "_all_met100", tag_njets, flav_tag_sl, 150. );
	  //met > 150 GeV requirement 
	  if ( tree->t1met10_ > 150. ) 
	    makeCR5Plots( tree, evtweight, h_1d_cr5, "_all_met150", tag_njets, flav_tag_sl, 120. );
	  //met > 200 GeV requirement 
	  if ( tree->t1met10_ > 200. ) 
	    makeCR5Plots( tree, evtweight, h_1d_cr5, "_all_met200", tag_njets, flav_tag_sl, 120. );
	  //met > 250 GeV requirement 
	  if ( tree->t1met10_ > 250. ) 
	    makeCR5Plots( tree, evtweight, h_1d_cr5, "_all_met250", tag_njets, flav_tag_sl, 120. );

	  // sample with only 1 lepton - this is the true CR5
	  if ( tree->ngoodlep_ == 1 ) {

	    //default 
	    makeCR5Plots( tree, evtweight, h_1d_cr5, "", tag_njets, flav_tag_sl, 150. );
	    //met > 50 GeV requirement 
	    if ( tree->t1met10_ > 50. ) 
	      makeCR5Plots( tree, evtweight, h_1d_cr5, "_met50", tag_njets, flav_tag_sl, 150. );
	    //met > 100 GeV requirement 
	    if ( tree->t1met10_ > 100. ) 
	      makeCR5Plots( tree, evtweight, h_1d_cr5, "_met100", tag_njets, flav_tag_sl, 150. );
	    //met > 150 GeV requirement 
	    if ( tree->t1met10_ > 150. ) 
	      makeCR5Plots( tree, evtweight, h_1d_cr5, "_met150", tag_njets, flav_tag_sl, 120. );
	    //met > 200 GeV requirement 
	    if ( tree->t1met10_ > 200. ) 
	      makeCR5Plots( tree, evtweight, h_1d_cr5, "_met200", tag_njets, flav_tag_sl, 120. );
	    //met > 250 GeV requirement 
	    if ( tree->t1met10_ > 250. ) 
	      makeCR5Plots( tree, evtweight, h_1d_cr5, "_met250", tag_njets, flav_tag_sl, 120. );

	  }
	}

    } // end event loop

    // delete tree;

  } // end file loop

    //
    // finish
    //

  TFile outfile(m_outfilename_.c_str(),"RECREATE") ; 
  printf("[StopTreeLooper::loop] Saving histograms to %s\n", m_outfilename_.c_str());

  std::map<std::string, TH1F*>::iterator it1d;
  for(it1d=h_1d.begin(); it1d!=h_1d.end(); it1d++) {
    it1d->second->Write(); 
    delete it1d->second;
  }

  outfile.Write();
  outfile.Close();
    
  //control regions
  //h_1d_cr1, h_1d_cr2, h_1d_cr4, h_1d_cr5

  TFile outfile_cr1(Form("CR1%s",m_outfilename_.c_str()),"RECREATE") ; 
  printf("[StopTreeLooper::loop] Saving CR1 histograms to %s\n", m_outfilename_.c_str());

  std::map<std::string, TH1F*>::iterator it1d_cr1;
  for(it1d_cr1=h_1d_cr1.begin(); it1d_cr1!=h_1d_cr1.end(); it1d_cr1++) {
    it1d_cr1->second->Write(); 
    delete it1d_cr1->second;
  }

  outfile_cr1.Write();
  outfile_cr1.Close();

  TFile outfile_cr2(Form("CR2%s",m_outfilename_.c_str()),"RECREATE") ; 
  printf("[StopTreeLooper::loop] Saving CR2 histograms to %s\n", m_outfilename_.c_str());

  std::map<std::string, TH1F*>::iterator it1d_cr2;
  for(it1d_cr2=h_1d_cr2.begin(); it1d_cr2!=h_1d_cr2.end(); it1d_cr2++) {
    it1d_cr2->second->Write(); 
    delete it1d_cr2->second;
  }

  outfile_cr2.Write();
  outfile_cr2.Close();
    
  TFile outfile_cr4(Form("CR4%s",m_outfilename_.c_str()),"RECREATE") ; 
  printf("[StopTreeLooper::loop] Saving CR4 histograms to %s\n", m_outfilename_.c_str());

  std::map<std::string, TH1F*>::iterator it1d_cr4;
  for(it1d_cr4=h_1d_cr4.begin(); it1d_cr4!=h_1d_cr4.end(); it1d_cr4++) {
    it1d_cr4->second->Write(); 
    delete it1d_cr4->second;
  }

  outfile_cr4.Write();
  outfile_cr4.Close();

  TFile outfile_cr5(Form("CR5%s",m_outfilename_.c_str()),"RECREATE") ; 
  printf("[StopTreeLooper::loop] Saving CR5 histograms to %s\n", m_outfilename_.c_str());

  std::map<std::string, TH1F*>::iterator it1d_cr5;
  for(it1d_cr5=h_1d_cr5.begin(); it1d_cr5!=h_1d_cr5.end(); it1d_cr5++) {
    it1d_cr5->second->Write(); 
    delete it1d_cr5->second;
  }

  outfile_cr5.Write();
  outfile_cr5.Close();
    
  TFile outfile_nj(Form("NJ%s",m_outfilename_.c_str()),"RECREATE") ; 
  printf("[StopTreeLooper::loop] Saving NJ histograms to %s\n", m_outfilename_.c_str());

  std::map<std::string, TH1F*>::iterator it1d_nj;
  for(it1d_nj=h_1d_nj.begin(); it1d_nj!=h_1d_nj.end(); it1d_nj++) {
    it1d_nj->second->Write(); 
    delete it1d_nj->second;
  }

  outfile_nj.Write();
  outfile_nj.Close();

  TFile outfile_z(Form("Z%s",m_outfilename_.c_str()),"RECREATE") ; 
  printf("[StopTreeLooper::loop] Saving Z histograms to %s\n", m_outfilename_.c_str());

  std::map<std::string, TH1F*>::iterator it1d_z;
  for(it1d_z=h_1d_z.begin(); it1d_z!=h_1d_z.end(); it1d_z++) {
    it1d_z->second->Write(); 
    delete it1d_z->second;
  }

  outfile_z.Write();
  outfile_z.Close();

  gROOT->cd();

}

bool StopTreeLooper::passEvtSelection(const StopTree *sTree) 
{
  //lepton + 2j selection

  //rho requirement
  if ( sTree->rhovor_<0. || sTree->rhovor_>=40. ) return false;

  //met filters
  if ( sTree->csc_      != 0 ) return false;
  if ( sTree->hbhe_     != 1 ) return false;
  if ( sTree->hcallaser_!= 1 ) return false;
  if ( sTree->ecaltp_   != 1 ) return false;
  if ( sTree->trkfail_  != 1 ) return false;
  if ( sTree->eebadsc_  != 1 ) return false;
  if ( sTree->hbhenew_  != 1 ) return false;

  //2-jet requirement
  //  if ( sTree->npfjets30_<2) return false;

  //at least 1 lepton
  if ( sTree->ngoodlep_ < 1 ) return false;

  //if have more than 1 lepton, remove cases where have 2 close together
  if ( sTree->ngoodlep_ > 1 && 
       dRbetweenVectors( sTree->lep1_ ,  sTree->lep2_ )<0.1 ) return false;

  return true;

}

bool StopTreeLooper::passSingleMuonSelection(const StopTree *sTree) 
{
  //single lepton muon selection for 8 TeV 53 analysis
  
  //at least one lepton
  if ( sTree->ngoodlep_ < 1 ) return false;

  //lepton flavor - ask for muon
  if ( sTree->leptype_ != 1 ) return false;

  //pass trigger if data - single muon
  if ( sTree->isomu24_ != 1 ) return false;
  if ( sTree->trgmu1_ != 1 )  return false;

  //passes pt and eta requirements
  if ( sTree->lep1_.Pt() < 30 )          return false;
  if ( fabs(sTree->lep1_.Eta() ) > 2.1)  return false;
  
  //pass isolated track veto
  //unfortunately changed default value to 9999.
  if ( sTree->pfcandpt10_ <9998. && sTree->pfcandiso10_ < 0.1 ) return false;

  return true;

}

bool StopTreeLooper::passSingleElecSelection(const StopTree *sTree) 
{
  //single lepton electron selection for 8 TeV 53 analysis

  //at least one lepton
  if ( sTree->ngoodlep_ < 1 ) return false;

  //lepton flavor - ask for electron
  if ( sTree->leptype_ != 0 ) return false;

  //pass trigger if data - single electron
  if ( sTree->ele27wp80_ != 1 ) return false;
  if ( sTree->trgel1_ != 1 )  return false;

  //passes pt and eta requirements
  if ( sTree->lep1_.Pt() < 30 )          return false;
  if ( fabs(sTree->lep1_.Eta() ) > 2.1)  return false;
  
  //pass isolated track veto
  //unfortunately changed default value to 9999.
  if ( sTree->pfcandpt10_ <9998. && sTree->pfcandiso10_ < 0.1 ) return false;

  return true;

}

bool StopTreeLooper::passOneLeptonSelection(const StopTree *sTree) 
{
  //single lepton selection for 8 TeV 53 analysis
  if ( !passSingleLeptonSelection(sTree) ) return false;

  //pass isolated track veto
  //unfortunately changed default value to 9999.
  if ( sTree->pfcandpt10_ <9998. && sTree->pfcandiso10_ < 0.1 ) return false;

  return true;

}

bool StopTreeLooper::passSingleLeptonSelection(const StopTree *sTree) 
{
  //single lepton selection for 8 TeV 53 analysis

  //at least one lepton
  if ( sTree->ngoodlep_ < 1 ) return false;

  //lepton flavor - trigger, pt and eta requirements
  if ( sTree->lep1_.Pt() < 30 )          return false;
  if ( sTree->leptype_ == 0 ) {

    //pass trigger if data - single electron
    if ( sTree->ele27wp80_ != 1 ) return false;
    //    if ( sTree->trgel1_ != 1 )  return false;
    
    if ( fabs(sTree->lep1_.Eta() ) > 2.1)  return false;

  } else if ( sTree->leptype_ == 1 ) {

    //pass trigger if data - single muon
    if ( sTree->isomu24_ != 1 ) return false;
    //    if ( sTree->trgmu1_ != 1 )  return false;
    
    if ( fabs(sTree->lep1_.Eta() ) > 2.1)  return false;

  }

  return true;

}


bool StopTreeLooper::passDileptonSelection(const StopTree *sTree) 
{
  //two lepton selection for 8 TeV 53 analysis

  //exactly 2 leptons
  if ( sTree->ngoodlep_ != 2 ) return false;

  //opposite sign
  if ( sTree->id1_*sTree->id2_>0 ) return false;

  //pass trigger if data - dilepton
  if ( sTree->mm_ != 1 && sTree->me_ != 1 
       && sTree->em_ != 1 && sTree->ee_ != 1 ) return false;

  //passes pt and eta requirements
  if ( sTree->lep1_.Pt() < 20 )          return false;
  if ( sTree->lep2_.Pt() < 20 )          return false;
  if ( fabs(sTree->lep1_.Eta() ) > 2.4)  return false;
  if ( fabs(sTree->lep2_.Eta() ) > 2.4)  return false;
  
  return true;

}

bool StopTreeLooper::passLepPlusIsoTrkSelection(const StopTree *sTree) 
{
  //single lepton plus iso trk selection for 8 TeV 53 analysis

  //at least one lepton
  if ( !passSingleLeptonSelection(sTree) ) return false;

  //pass isolated track requirement
  //unfortunately changed default value to 9999.
  if ( sTree->pfcandpt10_ > 9990. || sTree->pfcandiso10_ > 0.1 ) return false;

  return true;

}

void StopTreeLooper::makeCR2Plots( const StopTree *sTree, float evtweight, std::map<std::string, TH1F*> &h_1d, 
				   string tag_selection, string tag_njets, string flav_tag_dl, float mtcut ) 
{

  //find positive lepton - this is the one that is combined with the pseudomet to form the mT
  bool isfirstp = (sTree->id1_ > 0) ? true : false;
  
  //recalculate met
  float metx = sTree->t1met10_ * cos( sTree->t1met10phi_ );
  float mety = sTree->t1met10_ * sin( sTree->t1met10phi_ );
          
  //recalculate the MET with the positive lepton
  metx += isfirstp ? sTree->lep1_.px() : sTree->lep2_.px();
  mety += isfirstp ? sTree->lep1_.py() : sTree->lep2_.py();
  
  float t1met10_lep    = sqrt(metx*metx + mety*mety);
  float t1met10phi_lep = atan2( mety , metx );
  
  //recalculate the MT with the negative lepton
  float t1met10mt_lep = isfirstp ?
    getMT( sTree->lep2_.Pt() , sTree->lep2_.Phi() , t1met10_lep , t1met10phi_lep ) :
    getMT( sTree->lep1_.Pt() , sTree->lep1_.Phi() , t1met10_lep , t1met10phi_lep );

  //binning for mT plots
  int nbins = 30;
  float h_xmin = 0.;
  float h_xmax = 300.;
  float x_ovflw = h_xmax-0.001;
  
  float pseudomt_count = -1.;
  if ( t1met10mt_lep > 60. 
       && t1met10mt_lep < 100. )    pseudomt_count = 0.5;
  else if ( t1met10mt_lep > mtcut ) pseudomt_count = 1.5;
  
  //default met
  plot1D("h_cr2_met"+tag_selection          +flav_tag_dl, min(sTree->t1met10_, x_ovflw), evtweight, h_1d, nbins, h_xmin, h_xmax);
  plot1D("h_cr2_met"+tag_selection+tag_njets+flav_tag_dl, min(sTree->t1met10_, x_ovflw), evtweight, h_1d, nbins, h_xmin, h_xmax);
  //pseudo-met
  plot1D("h_cr2_pseudomet"+tag_selection	  +flav_tag_dl, min(t1met10_lep, x_ovflw), evtweight, h_1d, nbins, h_xmin, h_xmax);
  plot1D("h_cr2_pseudomet"+tag_selection+tag_njets+flav_tag_dl, min(t1met10_lep, x_ovflw), evtweight, h_1d, nbins, h_xmin, h_xmax);
  //positive lepton pt - enters mT calculation
  float leppt = isfirstp ? sTree->lep1_.Pt() : sTree->lep2_.Pt();
  plot1D("h_cr2_leppt"+tag_selection          +flav_tag_dl, min(leppt, x_ovflw), evtweight, h_1d, nbins, h_xmin, h_xmax);
  plot1D("h_cr2_leppt"+tag_selection+tag_njets+flav_tag_dl, min(leppt, x_ovflw), evtweight, h_1d, nbins, h_xmin, h_xmax);
  //angle between pos-lep and pseudopseudomet
  float dphi_pseudometlep = isfirstp ?
    getdphi( sTree->lep2_.Phi() , t1met10phi_lep ) :
    getdphi( sTree->lep1_.Phi() , t1met10phi_lep );
  plot1D("h_cr2_dphi_pseudometlep"+tag_selection          +flav_tag_dl, dphi_pseudometlep, evtweight, h_1d, 15, 0., TMath::Pi());
  plot1D("h_cr2_dphi_pseudometlep"+tag_selection+tag_njets+flav_tag_dl, dphi_pseudometlep, evtweight, h_1d, 15, 0., TMath::Pi());
  //pseudo-mt
  plot1D("h_cr2_pseudomt"      +tag_selection          +flav_tag_dl, min(t1met10mt_lep, x_ovflw), evtweight, h_1d, nbins, h_xmin, h_xmax);
  plot1D("h_cr2_pseudomt"      +tag_selection+tag_njets+flav_tag_dl, min(t1met10mt_lep, x_ovflw), evtweight, h_1d, nbins, h_xmin, h_xmax);
  plot1D("h_cr2_pseudomt_count"+tag_selection          +flav_tag_dl, pseudomt_count, evtweight, h_1d, 2, 0, 2);
  plot1D("h_cr2_pseudomt_count"+tag_selection+tag_njets+flav_tag_dl, pseudomt_count, evtweight, h_1d, 2, 0, 2);


}

void StopTreeLooper::makeCR4Plots( const StopTree *sTree, float evtweight, std::map<std::string, TH1F*> &h_1d, 
				   string tag_selection, string tag_njets, string flav_tag_dl, float mtcut ) 
{

  //binning for mT plots
  int nbins = 50;
  float h_xmin = 0.;
  float h_xmax = 500.;
  float x_ovflw = h_xmax-0.001;
  
  float mt_count = -1.;
  if ( sTree->t1met10mt_ > 60. 
       && sTree->t1met10mt_ < 100. )    mt_count = 0.5;
  else if ( sTree->t1met10mt_ > mtcut ) mt_count = 1.5;
  
  //default met
  plot1D("h_cr4_met"+tag_selection          +flav_tag_dl, min(sTree->t1met10_, x_ovflw), evtweight, h_1d, nbins, h_xmin, h_xmax);
  plot1D("h_cr4_met"+tag_selection+tag_njets+flav_tag_dl, min(sTree->t1met10_, x_ovflw), evtweight, h_1d, nbins, h_xmin, h_xmax);
  //leading lepton pt - enters mT calculation
  plot1D("h_cr4_leppt"+tag_selection          +flav_tag_dl, min(sTree->lep1_.Pt(), x_ovflw), evtweight, h_1d, nbins, h_xmin, h_xmax);
  plot1D("h_cr4_leppt"+tag_selection+tag_njets+flav_tag_dl, min(sTree->lep1_.Pt(), x_ovflw), evtweight, h_1d, nbins, h_xmin, h_xmax);
  //subleading lepton pt
  plot1D("h_cr4_subleadleppt"+tag_selection          +flav_tag_dl, min(sTree->lep2_.Pt(), x_ovflw), evtweight, h_1d, nbins, h_xmin, h_xmax);
  plot1D("h_cr4_subleadleppt"+tag_selection+tag_njets+flav_tag_dl, min(sTree->lep2_.Pt(), x_ovflw), evtweight, h_1d, nbins, h_xmin, h_xmax);
  //angle between lead-lep and met
  float dphi_metlep = getdphi( sTree->lep1_.Phi() , sTree->t1met10phi_ );
  plot1D("h_cr4_dphi_metlep"+tag_selection          +flav_tag_dl, dphi_metlep, evtweight, h_1d, 15, 0., TMath::Pi());
  plot1D("h_cr4_dphi_metlep"+tag_selection+tag_njets+flav_tag_dl, dphi_metlep, evtweight, h_1d, 15, 0., TMath::Pi());
  //MT
  plot1D("h_cr4_mt"+tag_selection          +flav_tag_dl, min(sTree->t1met10mt_, x_ovflw), evtweight, h_1d, nbins, h_xmin, h_xmax);
  plot1D("h_cr4_mt"+tag_selection+tag_njets+flav_tag_dl, min(sTree->t1met10mt_, x_ovflw), evtweight, h_1d, nbins, h_xmin, h_xmax);
  plot1D("h_cr4_mt_count"+tag_selection          +flav_tag_dl, mt_count, evtweight, h_1d, 2, 0, 2);
  plot1D("h_cr4_mt_count"+tag_selection+tag_njets+flav_tag_dl, mt_count, evtweight, h_1d, 2, 0, 2);

  //Plot more angles between various objects
  //angle between 2 leptons
  float dphi_dilep = getdphi( sTree->lep1_.Phi() ,  sTree->lep2_.Phi() );
  plot1D("h_cr4_dphi_dilep"+tag_selection          +flav_tag_dl, dphi_dilep, evtweight, h_1d, 15, 0., TMath::Pi());
  plot1D("h_cr4_dphi_dilep"+tag_selection+tag_njets+flav_tag_dl, dphi_dilep, evtweight, h_1d, 15, 0., TMath::Pi());
  //dR between 2 leptons
  float dR_dilep = dRbetweenVectors( sTree->lep1_ ,  sTree->lep2_ );
  plot1D("h_cr4_dR_dilep"+tag_selection          +flav_tag_dl, dR_dilep, evtweight, h_1d, 15, 0., TMath::Pi());
  plot1D("h_cr4_dR_dilep"+tag_selection+tag_njets+flav_tag_dl, dR_dilep, evtweight, h_1d, 15, 0., TMath::Pi());

}

void StopTreeLooper::makeCR5Plots( const StopTree *sTree, float evtweight, std::map<std::string, TH1F*> &h_1d, 
				   string tag_selection, string tag_njets, string flav_tag, float mtcut ) 
{

  //binning for mT plots
  int nbins = 50;
  float h_xmin = 0.;
  float h_xmax = 500.;
  float x_ovflw = h_xmax-0.001;
  
  float mt_count = -1.;
  if ( sTree->t1met10mt_ > 60. 
       && sTree->t1met10mt_ < 100. )    mt_count = 0.5;
  else if ( sTree->t1met10mt_ > mtcut ) mt_count = 1.5;
  
  //default met
  plot1D("h_cr5_met"+tag_selection          +flav_tag, min(sTree->t1met10_, x_ovflw), evtweight, h_1d, nbins, h_xmin, h_xmax);
  plot1D("h_cr5_met"+tag_selection+tag_njets+flav_tag, min(sTree->t1met10_, x_ovflw), evtweight, h_1d, nbins, h_xmin, h_xmax);
  //leading lepton pt - enters mT calculation
  plot1D("h_cr5_leppt"+tag_selection          +flav_tag, min(sTree->lep1_.Pt(), x_ovflw), evtweight, h_1d, nbins, h_xmin, h_xmax);
  plot1D("h_cr5_leppt"+tag_selection+tag_njets+flav_tag, min(sTree->lep1_.Pt(), x_ovflw), evtweight, h_1d, nbins, h_xmin, h_xmax);
  //isolated track pt
  plot1D("h_cr5_isotrkpt"+tag_selection          +flav_tag, min(sTree->pfcand10_.Pt(), x_ovflw), evtweight, h_1d, nbins, h_xmin, h_xmax);
  plot1D("h_cr5_isotrkpt"+tag_selection+tag_njets+flav_tag, min(sTree->pfcand10_.Pt(), x_ovflw), evtweight, h_1d, nbins, h_xmin, h_xmax);
  //angle between lead-lep and met
  float dphi_metlep = getdphi( sTree->lep1_.Phi() , sTree->t1met10phi_ );
  plot1D("h_cr5_dphi_metlep"+tag_selection          +flav_tag, dphi_metlep, evtweight, h_1d, 15, 0., TMath::Pi());
  plot1D("h_cr5_dphi_metlep"+tag_selection+tag_njets+flav_tag, dphi_metlep, evtweight, h_1d, 15, 0., TMath::Pi());
  //MT
  plot1D("h_cr5_mt"+tag_selection          +flav_tag, min(sTree->t1met10mt_, x_ovflw), evtweight, h_1d, nbins, h_xmin, h_xmax);
  plot1D("h_cr5_mt"+tag_selection+tag_njets+flav_tag, min(sTree->t1met10mt_, x_ovflw), evtweight, h_1d, nbins, h_xmin, h_xmax);
  plot1D("h_cr5_mt_count"+tag_selection          +flav_tag, mt_count, evtweight, h_1d, 2, 0, 2);
  plot1D("h_cr5_mt_count"+tag_selection+tag_njets+flav_tag, mt_count, evtweight, h_1d, 2, 0, 2);

  //Plot more angles between various objects
  //angle between lepton and isolated track
  float dphi_leptrk = getdphi( sTree->lep1_.Phi() ,  sTree->pfcand10_.Phi() );
  plot1D("h_cr5_dphi_leptrk"+tag_selection          +flav_tag, dphi_leptrk, evtweight, h_1d, 15, 0., TMath::Pi());
  plot1D("h_cr5_dphi_leptrk"+tag_selection+tag_njets+flav_tag, dphi_leptrk, evtweight, h_1d, 15, 0., TMath::Pi());
  //dR between lepton and isolated track
  float dR_leptrk = dRbetweenVectors( sTree->lep1_ ,  sTree->pfcand10_ );
  plot1D("h_cr5_dR_leptrk"+tag_selection          +flav_tag, dR_leptrk, evtweight, h_1d, 15, 0., TMath::Pi());
  plot1D("h_cr5_dR_leptrk"+tag_selection+tag_njets+flav_tag, dR_leptrk, evtweight, h_1d, 15, 0., TMath::Pi());

}

void StopTreeLooper::makeCR1Plots( const StopTree *sTree, float evtweight, std::map<std::string, TH1F*> &h_1d, 
				   string tag_selection, string tag_njets, string flav_tag, float mtcut ) 
{

 //binning for mT plots
  int nbins = 150;
  int nbins2 = 200;
  float h_xmin = 0.;
  float h_xmax = 1500.;
  float x_ovflw = h_xmax-0.001;
  
  float mt_count = -1.;
  if ( sTree->t1met10mt_ > 60. 
       && sTree->t1met10mt_ < 100. )    mt_count = 0.5;
  else if ( sTree->t1met10mt_ > mtcut ) mt_count = 1.5;
  
  float dphi_metlep = getdphi( sTree->lep1_.Phi() , sTree->t1met10phi_ );
  
  //if(sTree->t1met10mt_>250 && sTree->t1met10mt_<450) {
  //if(sTree->t1met10mt_>450 ) {
  //if(sTree->t1met10mt_>250 ) {
  //if(sTree->t1met10mt_>250 && fabs(sTree->pflep1_.Pt() -sTree->lep1_.Pt() )<10) {
  //if(sTree->t1met10mt_>250 && dphi_metlep<3. && fabs(sTree->pflep1_.Pt() -sTree->lep1_.Pt() )<10 && fabs(sTree->lep1_.Eta())<1.4442 ) {
  if(sTree->t1met10mt_>250 && dphi_metlep<(TMath::Pi()*19./20.) && fabs(sTree->pflep1_.Pt() -sTree->lep1_.Pt() )<10 && fabs(sTree->lep1_.Eta())<1.4442 ) {
  //if(sTree->t1met10mt_>0) {
  
  plot1D("h_cr1_njets"    +tag_selection+flav_tag, min(sTree->npfjets30_,4),  evtweight, h_1d, 5,0,5);
  plot1D("h_cr1_njets_all"+tag_selection+flav_tag, min(sTree->npfjets30_,9),  evtweight, h_1d, 10, 0, 10);
  
  //default met
  plot1D("h_cr1_met"+tag_selection          +flav_tag, min(sTree->t1met10_, x_ovflw), evtweight, h_1d, nbins2, 0., 1000.);
  plot1D("h_cr1_met"+tag_selection+tag_njets+flav_tag, min(sTree->t1met10_, x_ovflw), evtweight, h_1d, nbins2, 0., 1000.);
  
  //lepton pt - enters mT calculation
  plot1D("h_cr1_leppt"+tag_selection          +flav_tag, min(sTree->lep1_.Pt(), x_ovflw), evtweight, h_1d, nbins, h_xmin, h_xmax);
  plot1D("h_cr1_leppt"+tag_selection+tag_njets+flav_tag, min(sTree->lep1_.Pt(), x_ovflw), evtweight, h_1d, nbins, h_xmin, h_xmax);
  
  //angle between lepton and met
  //float dphi_metlep = getdphi( sTree->lep1_.Phi() , sTree->t1met10phi_ );
  plot1D("h_cr1_dphi_metlep"+tag_selection          +flav_tag, dphi_metlep, evtweight, h_1d, 200, 0., TMath::Pi());
  plot1D("h_cr1_dphi_metlep"+tag_selection+tag_njets+flav_tag, dphi_metlep, evtweight, h_1d, 200, 0., TMath::Pi());
  
  //MT
  plot1D("h_cr1_mt"+tag_selection          +flav_tag, min(sTree->t1met10mt_, x_ovflw), evtweight, h_1d, nbins, h_xmin, h_xmax);
  plot1D("h_cr1_mt"+tag_selection+tag_njets+flav_tag, min(sTree->t1met10mt_, x_ovflw), evtweight, h_1d, nbins, h_xmin, h_xmax);
  
  //lepton eta
  float x_ovflw2 = 3-0.001;
  plot1D("h_cr1_lepeta"+tag_selection          +flav_tag, (sTree->lep1_.Eta()>0) ? min(sTree->lep1_.Eta(), x_ovflw2) : max(sTree->lep1_.Eta(), -x_ovflw2), evtweight, h_1d, nbins, -3, 3);
  plot1D("h_cr1_lepeta"+tag_selection+tag_njets+flav_tag, (sTree->lep1_.Eta()>0) ? min(sTree->lep1_.Eta(), x_ovflw2) : max(sTree->lep1_.Eta(), -x_ovflw2), evtweight, h_1d, nbins, -3, 3);
  
  //lepton phi
  plot1D("h_cr1_lepphi"+tag_selection          +flav_tag, sTree->lep1_.Phi(), evtweight, h_1d, nbins, -TMath::Pi(),TMath::Pi());
  plot1D("h_cr1_lepphi"+tag_selection+tag_njets+flav_tag, sTree->lep1_.Phi(), evtweight, h_1d, nbins, -TMath::Pi(),TMath::Pi()) ;
  
  //angle between met and high pT jet
  float dphi_metjet1 = getdphi( sTree->pfjet1_.Phi() , sTree->t1met10phi_ );
  plot1D("h_cr1_dphi_metjet1"+tag_selection          +flav_tag, dphi_metjet1, evtweight, h_1d, 200, 0., TMath::Pi());
  plot1D("h_cr1_dphi_metjet1"+tag_selection+tag_njets+flav_tag, dphi_metjet1, evtweight, h_1d, 200, 0., TMath::Pi());
  
  //angle between met and closest jet
  float dphi_metjet[] = {getdphi( sTree->pfjet1_.Phi() , sTree->t1met10phi_ ), getdphi( sTree->pfjet2_.Phi() , sTree->t1met10phi_ ), getdphi( sTree->pfjet3_.Phi() , sTree->t1met10phi_ ), getdphi( sTree->pfjet4_.Phi() , sTree->t1met10phi_ ) };
  //cout<<sizeof(dphi_metjet)/sizeof(dphi_metjet[0])<<endl;
  std::sort(dphi_metjet, dphi_metjet + sizeof(dphi_metjet)/sizeof(dphi_metjet[0]) );
  //cout<<dphi_metjet[0]<<" "<<dphi_metjet[1]<<" "<<dphi_metjet[2]<<" "<<dphi_metjet[3]<<endl;
  float dphi_metclosestjet = ( ( (TMath::Pi() - dphi_metjet[3]) < dphi_metjet[0] ) ? dphi_metjet[3] : dphi_metjet[0] );
  //cout<<dphi_metclosestjet <<endl;
  plot1D("h_cr1_dphi_metclosestjet"+tag_selection          +flav_tag, dphi_metclosestjet, evtweight, h_1d, 200, 0., TMath::Pi());
  plot1D("h_cr1_dphi_metclosestjet"+tag_selection+tag_njets+flav_tag, dphi_metclosestjet, evtweight, h_1d, 200, 0., TMath::Pi());
  
  //nvtx
  plot1D("h_cr1_nvtx"    +tag_selection          +flav_tag, sTree->nvtx_,  evtweight, h_1d, 40,0,40);
  plot1D("h_cr1_nvtx"    +tag_selection+tag_njets+flav_tag, sTree->nvtx_,  evtweight, h_1d, 40,0,40);
  
  //vtxw
  plot1D("h_cr1_nvtxw"    +tag_selection          +flav_tag, sTree-> nvtxweight_,  evtweight, h_1d, 41, -4., 4.);
  plot1D("h_cr1_nvtxw"    +tag_selection+tag_njets+flav_tag, sTree-> nvtxweight_,  evtweight, h_1d, 41, -4., 4.);

  //rho
  plot1D("h_cr1_rho"    +tag_selection          +flav_tag, sTree-> rhovor_,  evtweight, h_1d, 80,0,40);
  plot1D("h_cr1_rho"    +tag_selection+tag_njets+flav_tag, sTree-> rhovor_,  evtweight, h_1d, 80,0,40);

  
  //nbtagscsvl
  plot1D("h_cr1_nbtagscsvl"    +tag_selection          +flav_tag, sTree-> nbtagscsvl_,  evtweight, h_1d, 6,0,6);
  plot1D("h_cr1_nbtagscsvl"    +tag_selection+tag_njets+flav_tag, sTree-> nbtagscsvl_,  evtweight, h_1d, 6,0,6);
  
  //nbtagscsvm
  plot1D("h_cr1_nbtagscsvm"    +tag_selection          +flav_tag, sTree-> nbtagscsvm_,  evtweight, h_1d, 6,0,6);
  plot1D("h_cr1_nbtagscsvm"    +tag_selection+tag_njets+flav_tag, sTree-> nbtagscsvm_,  evtweight, h_1d, 6,0,6);
  
  float t1met10mtPFlep_     = getMT( sTree->pflep1_.Pt() , sTree->pflep1_.Phi() , sTree->t1met10_     , sTree->t1met10phi_ );
  //MT with pflep
  plot1D("h_cr1_mtpflep"+tag_selection          +flav_tag, min(t1met10mtPFlep_, x_ovflw), evtweight, h_1d, nbins, h_xmin, h_xmax);
  plot1D("h_cr1_mtpflep"+tag_selection+tag_njets+flav_tag, min(t1met10mtPFlep_, x_ovflw), evtweight, h_1d, nbins, h_xmin, h_xmax);
  
  //pf lepton pt - enters mT calculation
  plot1D("h_cr1_pfleppt"+tag_selection          +flav_tag, min(sTree->pflep1_.Pt(), x_ovflw), evtweight, h_1d, nbins, h_xmin, h_xmax);
  plot1D("h_cr1_pfleppt"+tag_selection+tag_njets+flav_tag, min(sTree->pflep1_.Pt(), x_ovflw), evtweight, h_1d, nbins, h_xmin, h_xmax);
  
  //angle between pf lepton and met
  float dphi_metpflep = getdphi( sTree->pflep1_.Phi() , sTree->t1met10phi_ );
  plot1D("h_cr1_dphi_metpflep"+tag_selection          +flav_tag, dphi_metpflep, evtweight, h_1d, 200, 0., TMath::Pi());
  plot1D("h_cr1_dphi_metpflep"+tag_selection+tag_njets+flav_tag, dphi_metpflep, evtweight, h_1d, 200, 0., TMath::Pi());
  
  //pf lepton eta
  plot1D("h_cr1_pflepeta"+tag_selection          +flav_tag, (sTree->pflep1_.Eta()>0) ? min(sTree->pflep1_.Eta(), x_ovflw2) : max(sTree->pflep1_.Eta(), -x_ovflw2), evtweight, h_1d, nbins, -3, 3);
  plot1D("h_cr1_pflepeta"+tag_selection+tag_njets+flav_tag, (sTree->pflep1_.Eta()>0) ? min(sTree->pflep1_.Eta(), x_ovflw2) : max(sTree->pflep1_.Eta(), -x_ovflw2), evtweight, h_1d, nbins, -3, 3);
  
  //pf lepton phi
  plot1D("h_cr1_pflepphi"+tag_selection          +flav_tag, sTree->pflep1_.Phi(), evtweight, h_1d, nbins, -TMath::Pi(),TMath::Pi());
  plot1D("h_cr1_pflepphi"+tag_selection+tag_njets+flav_tag, sTree->pflep1_.Phi(), evtweight, h_1d, nbins, -TMath::Pi(),TMath::Pi()) ;
  
   //phi between pf lepton and lepton
  float dphi_leppflep = getdphi( sTree->pflep1_.Phi() , sTree->lep1_.Phi());
  plot1D("h_cr1_dphi_leppflep"+tag_selection          +flav_tag, min(dphi_leppflep,float(TMath::Pi()/8 - 0.0001)), evtweight, h_1d, 200, 0., TMath::Pi()/8);
  plot1D("h_cr1_dphi_leppflep"+tag_selection+tag_njets+flav_tag, min(dphi_leppflep,float(TMath::Pi()/8 - 0.0001)), evtweight, h_1d, 200, 0., TMath::Pi()/8);
  
  //ratio of pT between pf lepton and lepton
  //plot1D("h_cr1_pTrat_leppflep"+tag_selection          +flav_tag, sTree->pflep1_.Pt() / sTree->lep1_.Pt(), evtweight, h_1d, 200, 0.5, 1.5);
  //plot1D("h_cr1_pTrat_leppflep"+tag_selection+tag_njets+flav_tag, sTree->pflep1_.Pt() / sTree->lep1_.Pt(), evtweight, h_1d, 200, 0.5, 1.5);
  x_ovflw2 = 1.5;
  plot1D("h_cr1_pTrat_leppflep"+tag_selection          +flav_tag, ((sTree->pflep1_.Pt() / sTree->lep1_.Pt())>1) ? min((sTree->pflep1_.Pt() / sTree->lep1_.Pt()), x_ovflw2) : max((sTree->pflep1_.Pt() / sTree->lep1_.Pt()), 2-x_ovflw2), evtweight, h_1d, 200, 2-x_ovflw2, x_ovflw2);
  plot1D("h_cr1_pTrat_leppflep"+tag_selection+tag_njets+flav_tag, ((sTree->pflep1_.Pt() / sTree->lep1_.Pt())>1) ? min((sTree->pflep1_.Pt() / sTree->lep1_.Pt()), x_ovflw2) : max((sTree->pflep1_.Pt() / sTree->lep1_.Pt()), 2-x_ovflw2), evtweight, h_1d, 200, 2-x_ovflw2, x_ovflw2);

  //difference of pT between pf lepton and lepton
  //plot1D("h_cr1_pTdiff_leppflep"+tag_selection          +flav_tag, sTree->pflep1_.Pt() - sTree->lep1_.Pt(), evtweight, h_1d, 200, -10, 10);
  //plot1D("h_cr1_pTdiff_leppflep"+tag_selection+tag_njets+flav_tag, sTree->pflep1_.Pt() - sTree->lep1_.Pt(), evtweight, h_1d, 200, -10, 10);
  x_ovflw2 = 40.;
  plot1D("h_cr1_pTdiff_leppflep"+tag_selection          +flav_tag, ((sTree->pflep1_.Pt() - sTree->lep1_.Pt())>0) ? min((sTree->pflep1_.Pt() - sTree->lep1_.Pt()), x_ovflw2) : max((sTree->pflep1_.Pt() - sTree->lep1_.Pt()), -x_ovflw2), evtweight, h_1d, 200, -x_ovflw2, x_ovflw2);
  plot1D("h_cr1_pTdiff_leppflep"+tag_selection+tag_njets+flav_tag, ((sTree->pflep1_.Pt() - sTree->lep1_.Pt())>0) ? min((sTree->pflep1_.Pt() - sTree->lep1_.Pt()), x_ovflw2) : max((sTree->pflep1_.Pt() - sTree->lep1_.Pt()), -x_ovflw2), evtweight, h_1d, 200, -x_ovflw2, x_ovflw2);

  
  h_xmin = 0;
  h_xmax = 10;
  plot1D("h_cr1_nwzpartons"+tag_selection          +flav_tag, (sTree->nwzpartons_>0) ? min(sTree->nwzpartons_, 10) : max(sTree->nwzpartons_, 0), evtweight, h_1d, nbins2, h_xmin, h_xmax);
  plot1D("h_cr1_nwzpartons"+tag_selection+tag_njets+flav_tag, (sTree->nwzpartons_>0) ? min(sTree->nwzpartons_, 10) : max(sTree->nwzpartons_, 0), evtweight, h_1d, nbins2, h_xmin, h_xmax);

  h_xmin = 0.;
  h_xmax = 1000.;
  plot1D("h_cr1_trkmet_nolepcorr"+tag_selection          +flav_tag, (sTree->trkmet_nolepcorr_>0) ? min(sTree->trkmet_nolepcorr_, float(h_xmax-0.0001)) : max(sTree->trkmet_nolepcorr_, float(h_xmin+0.0001)), evtweight, h_1d, nbins2, h_xmin, h_xmax);
  plot1D("h_cr1_trkmet_nolepcorr"+tag_selection+tag_njets+flav_tag, (sTree->trkmet_nolepcorr_>0) ? min(sTree->trkmet_nolepcorr_, float(h_xmax-0.0001)) : max(sTree->trkmet_nolepcorr_, float(h_xmin+0.0001)), evtweight, h_1d, nbins2, h_xmin, h_xmax);

  h_xmin = -TMath::Pi();
  h_xmax = TMath::Pi();
  plot1D("h_cr1_trkmetphi_nolepcorr"+tag_selection          +flav_tag, (sTree->trkmetphi_nolepcorr_>0) ? min(sTree->trkmetphi_nolepcorr_, float(h_xmax-0.0001)) : max(sTree->trkmetphi_nolepcorr_, float(h_xmin+0.0001)), evtweight, h_1d, nbins2, h_xmin, h_xmax);
  plot1D("h_cr1_trkmetphi_nolepcorr"+tag_selection+tag_njets+flav_tag, (sTree->trkmetphi_nolepcorr_>0) ? min(sTree->trkmetphi_nolepcorr_, float(h_xmax-0.0001)) : max(sTree->trkmetphi_nolepcorr_, float(h_xmin+0.0001)), evtweight, h_1d, nbins2, h_xmin, h_xmax);

  h_xmin = 0.;
  h_xmax = 1000.;
  plot1D("h_cr1_trkmet"+tag_selection          +flav_tag, (sTree->trkmet_>0) ? min(sTree->trkmet_, float(h_xmax-0.0001)) : max(sTree->trkmet_, float(h_xmin+0.0001)), evtweight, h_1d, nbins2, h_xmin, h_xmax);
  plot1D("h_cr1_trkmet"+tag_selection+tag_njets+flav_tag, (sTree->trkmet_>0) ? min(sTree->trkmet_, float(h_xmax-0.0001)) : max(sTree->trkmet_, float(h_xmin+0.0001)), evtweight, h_1d, nbins2, h_xmin, h_xmax);

  h_xmin = -TMath::Pi();
  h_xmax = TMath::Pi();
  plot1D("h_cr1_trkmetphi"+tag_selection          +flav_tag, (sTree->trkmetphi_>0) ? min(sTree->trkmetphi_, float(h_xmax-0.0001)) : max(sTree->trkmetphi_, float(h_xmin+0.0001)), evtweight, h_1d, nbins2, h_xmin, h_xmax);
  plot1D("h_cr1_trkmetphi"+tag_selection+tag_njets+flav_tag, (sTree->trkmetphi_>0) ? min(sTree->trkmetphi_, float(h_xmax-0.0001)) : max(sTree->trkmetphi_, float(h_xmin+0.0001)), evtweight, h_1d, nbins2, h_xmin, h_xmax);

  h_xmin = 0.;
  h_xmax = 0.5;
  plot1D("h_cr1_isopf1"+tag_selection          +flav_tag, (sTree->isopf1_>0) ? min(sTree->isopf1_, float(h_xmax-0.0001)) : max(sTree->isopf1_, float(h_xmin+0.0001)), evtweight, h_1d, nbins2, h_xmin, h_xmax);
  plot1D("h_cr1_isopf1"+tag_selection+tag_njets+flav_tag, (sTree->isopf1_>0) ? min(sTree->isopf1_, float(h_xmax-0.0001)) : max(sTree->isopf1_, float(h_xmin+0.0001)), evtweight, h_1d, nbins2, h_xmin, h_xmax);

  h_xmin = 0.;
  h_xmax = 0.5;
  plot1D("h_cr1_isopfold1"+tag_selection          +flav_tag, (sTree->isopfold1_>0) ? min(sTree->isopfold1_, float(h_xmax-0.0001)) : max(sTree->isopfold1_, float(h_xmin+0.0001)), evtweight, h_1d, nbins2, h_xmin, h_xmax);
  plot1D("h_cr1_isopfold1"+tag_selection+tag_njets+flav_tag, (sTree->isopfold1_>0) ? min(sTree->isopfold1_, float(h_xmax-0.0001)) : max(sTree->isopfold1_, float(h_xmin+0.0001)), evtweight, h_1d, nbins2, h_xmin, h_xmax);

  h_xmin = 0.;
  h_xmax = 80.;
  plot1D("h_cr1_eoverpin"+tag_selection          +flav_tag, (sTree->eoverpin_>0) ? min(sTree->eoverpin_, float(h_xmax-0.0001)) : max(sTree->eoverpin_, float(h_xmin+0.0001)), evtweight, h_1d, nbins2, h_xmin, h_xmax);
  plot1D("h_cr1_eoverpin"+tag_selection+tag_njets+flav_tag, (sTree->eoverpin_>0) ? min(sTree->eoverpin_, float(h_xmax-0.0001)) : max(sTree->eoverpin_, float(h_xmin+0.0001)), evtweight, h_1d, nbins2, h_xmin, h_xmax);

  h_xmin = 0.;
  h_xmax = 80.;
  plot1D("h_cr1_eoverpout"+tag_selection          +flav_tag, (sTree->eoverpout_>0) ? min(sTree->eoverpout_, float(h_xmax-0.0001)) : max(sTree->eoverpout_, float(h_xmin+0.0001)), evtweight, h_1d, nbins2, h_xmin, h_xmax);
  plot1D("h_cr1_eoverpout"+tag_selection+tag_njets+flav_tag, (sTree->eoverpout_>0) ? min(sTree->eoverpout_, float(h_xmax-0.0001)) : max(sTree->eoverpout_, float(h_xmin+0.0001)), evtweight, h_1d, nbins2, h_xmin, h_xmax);

  h_xmin = -0.02;
  h_xmax = 0.02;
  plot1D("h_cr1_dEtaIn"+tag_selection          +flav_tag, (sTree->dEtaIn_>0) ? min(sTree->dEtaIn_, float(h_xmax-0.0001)) : max(sTree->dEtaIn_, float(h_xmin+0.0001)), evtweight, h_1d, nbins2, h_xmin, h_xmax);
  plot1D("h_cr1_dEtaIn"+tag_selection+tag_njets+flav_tag, (sTree->dEtaIn_>0) ? min(sTree->dEtaIn_, float(h_xmax-0.0001)) : max(sTree->dEtaIn_, float(h_xmin+0.0001)), evtweight, h_1d, nbins2, h_xmin, h_xmax);

  h_xmin = -0.1;
  h_xmax = 0.1;
  plot1D("h_cr1_dPhiIn"+tag_selection          +flav_tag, (sTree->dPhiIn_>0) ? min(sTree->dPhiIn_, float(h_xmax-0.0001)) : max(sTree->dPhiIn_, float(h_xmin+0.0001)), evtweight, h_1d, nbins2, h_xmin, h_xmax);
  plot1D("h_cr1_dPhiIn"+tag_selection+tag_njets+flav_tag, (sTree->dPhiIn_>0) ? min(sTree->dPhiIn_, float(h_xmax-0.0001)) : max(sTree->dPhiIn_, float(h_xmin+0.0001)), evtweight, h_1d, nbins2, h_xmin, h_xmax);

  h_xmin = 0.;
  h_xmax = 0.1;
  plot1D("h_cr1_sigmaIEtaIEta"+tag_selection          +flav_tag, (sTree->sigmaIEtaIEta_>0) ? min(sTree->sigmaIEtaIEta_, float(h_xmax-0.0001)) : max(sTree->sigmaIEtaIEta_, float(h_xmin+0.0001)), evtweight, h_1d, nbins2, h_xmin, h_xmax);
  plot1D("h_cr1_sigmaIEtaIEta"+tag_selection+tag_njets+flav_tag, (sTree->sigmaIEtaIEta_>0) ? min(sTree->sigmaIEtaIEta_, float(h_xmax-0.0001)) : max(sTree->sigmaIEtaIEta_, float(h_xmin+0.0001)), evtweight, h_1d, nbins2, h_xmin, h_xmax);

  h_xmin = 0.;
  h_xmax = 0.2;
  plot1D("h_cr1_hOverE"+tag_selection          +flav_tag, (sTree->hOverE_>0) ? min(sTree->hOverE_, float(h_xmax-0.0001)) : max(sTree->hOverE_, float(h_xmin+0.0001)), evtweight, h_1d, nbins2, h_xmin, h_xmax);
  plot1D("h_cr1_hOverE"+tag_selection+tag_njets+flav_tag, (sTree->hOverE_>0) ? min(sTree->hOverE_, float(h_xmax-0.0001)) : max(sTree->hOverE_, float(h_xmin+0.0001)), evtweight, h_1d, nbins2, h_xmin, h_xmax);

  h_xmin = 0.;
  h_xmax = 0.1;
  plot1D("h_cr1_ooemoop"+tag_selection          +flav_tag, (sTree->ooemoop_>0) ? min(sTree->ooemoop_, float(h_xmax-0.0001)) : max(sTree->ooemoop_, float(h_xmin+0.0001)), evtweight, h_1d, nbins2, h_xmin, h_xmax);
  plot1D("h_cr1_ooemoop"+tag_selection+tag_njets+flav_tag, (sTree->ooemoop_>0) ? min(sTree->ooemoop_, float(h_xmax-0.0001)) : max(sTree->ooemoop_, float(h_xmin+0.0001)), evtweight, h_1d, nbins2, h_xmin, h_xmax);

  h_xmin = -0.1;
  h_xmax = 0.1;
  plot1D("h_cr1_d0vtx"+tag_selection          +flav_tag, (sTree->d0vtx_>0) ? min(sTree->d0vtx_, float(h_xmax-0.0001)) : max(sTree->d0vtx_, float(h_xmin+0.0001)), evtweight, h_1d, nbins2, h_xmin, h_xmax);
  plot1D("h_cr1_d0vtx"+tag_selection+tag_njets+flav_tag, (sTree->d0vtx_>0) ? min(sTree->d0vtx_, float(h_xmax-0.0001)) : max(sTree->d0vtx_, float(h_xmin+0.0001)), evtweight, h_1d, nbins2, h_xmin, h_xmax);

  h_xmin = -0.20;
  h_xmax = 0.20;
  plot1D("h_cr1_dzvtx"+tag_selection          +flav_tag, (sTree->dzvtx_>0) ? min(sTree->dzvtx_, float(h_xmax-0.0001)) : max(sTree->dzvtx_, float(h_xmin+0.0001)), evtweight, h_1d, nbins2, h_xmin, h_xmax);
  plot1D("h_cr1_dzvtx"+tag_selection+tag_njets+flav_tag, (sTree->dzvtx_>0) ? min(sTree->dzvtx_, float(h_xmax-0.0001)) : max(sTree->dzvtx_, float(h_xmin+0.0001)), evtweight, h_1d, nbins2, h_xmin, h_xmax);

  h_xmin = 0;
  h_xmax = 4;
  plot1D("h_cr1_expinnerlayers"+tag_selection          +flav_tag, (sTree->expinnerlayers_>0) ? min(sTree->expinnerlayers_, float(h_xmax-0.0001)) : max(sTree->expinnerlayers_, float(h_xmin+0.0001)), evtweight, h_1d, nbins2, h_xmin, h_xmax);
  plot1D("h_cr1_expinnerlayers"+tag_selection+tag_njets+flav_tag, (sTree->expinnerlayers_>0) ? min(sTree->expinnerlayers_, float(h_xmax-0.0001)) : max(sTree->expinnerlayers_, float(h_xmin+0.0001)), evtweight, h_1d, nbins2, h_xmin, h_xmax);

  h_xmin = -1.;
  h_xmax = 1.;
  plot1D("h_cr1_fbrem"+tag_selection          +flav_tag, (sTree->fbrem_>0) ? min(sTree->fbrem_, float(h_xmax-0.0001)) : max(sTree->fbrem_, float(h_xmin+0.0001)), evtweight, h_1d, nbins2, h_xmin, h_xmax);
  plot1D("h_cr1_fbrem"+tag_selection+tag_njets+flav_tag, (sTree->fbrem_>0) ? min(sTree->fbrem_, float(h_xmax-0.0001)) : max(sTree->fbrem_, float(h_xmin+0.0001)), evtweight, h_1d, nbins2, h_xmin, h_xmax);

  h_xmin = 0.;
  h_xmax = 50.;
  plot1D("h_cr1_pfisoch"+tag_selection          +flav_tag, (sTree->pfisoch_>0) ? min(sTree->pfisoch_, float(h_xmax-0.0001)) : max(sTree->pfisoch_, float(h_xmin+0.0001)), evtweight, h_1d, nbins2, h_xmin, h_xmax);
  plot1D("h_cr1_pfisoch"+tag_selection+tag_njets+flav_tag, (sTree->pfisoch_>0) ? min(sTree->pfisoch_, float(h_xmax-0.0001)) : max(sTree->pfisoch_, float(h_xmin+0.0001)), evtweight, h_1d, nbins2, h_xmin, h_xmax);

  h_xmin = 0.;
  h_xmax = 50.;
  plot1D("h_cr1_pfisoem"+tag_selection          +flav_tag, (sTree->pfisoem_>0) ? min(sTree->pfisoem_, float(h_xmax-0.0001)) : max(sTree->pfisoem_, float(h_xmin+0.0001)), evtweight, h_1d, nbins2, h_xmin, h_xmax);
  plot1D("h_cr1_pfisoem"+tag_selection+tag_njets+flav_tag, (sTree->pfisoem_>0) ? min(sTree->pfisoem_, float(h_xmax-0.0001)) : max(sTree->pfisoem_, float(h_xmin+0.0001)), evtweight, h_1d, nbins2, h_xmin, h_xmax);

  h_xmin = 0.;
  h_xmax = 50.;
  plot1D("h_cr1_pfisonh"+tag_selection          +flav_tag, (sTree->pfisonh_>0) ? min(sTree->pfisonh_, float(h_xmax-0.0001)) : max(sTree->pfisonh_, float(h_xmin+0.0001)), evtweight, h_1d, nbins2, h_xmin, h_xmax);
  plot1D("h_cr1_pfisonh"+tag_selection+tag_njets+flav_tag, (sTree->pfisonh_>0) ? min(sTree->pfisonh_, float(h_xmax-0.0001)) : max(sTree->pfisonh_, float(h_xmin+0.0001)), evtweight, h_1d, nbins2, h_xmin, h_xmax);

  h_xmin = 0.;
  h_xmax = 1500.;
  plot1D("h_cr1_eSC"+tag_selection          +flav_tag, (sTree->eSC_>0) ? min(sTree->eSC_, float(h_xmax-0.0001)) : max(sTree->eSC_, float(h_xmin+0.0001)), evtweight, h_1d, nbins2, h_xmin, h_xmax);
  plot1D("h_cr1_eSC"+tag_selection+tag_njets+flav_tag, (sTree->eSC_>0) ? min(sTree->eSC_, float(h_xmax-0.0001)) : max(sTree->eSC_, float(h_xmin+0.0001)), evtweight, h_1d, nbins2, h_xmin, h_xmax);

  h_xmin = -TMath::Pi();
  h_xmax = TMath::Pi();
  plot1D("h_cr1_phiSC"+tag_selection          +flav_tag, (sTree->phiSC_>0) ? min(sTree->phiSC_, float(h_xmax-0.0001)) : max(sTree->phiSC_, float(h_xmin+0.0001)), evtweight, h_1d, nbins2, h_xmin, h_xmax);
  plot1D("h_cr1_phiSC"+tag_selection+tag_njets+flav_tag, (sTree->phiSC_>0) ? min(sTree->phiSC_, float(h_xmax-0.0001)) : max(sTree->phiSC_, float(h_xmin+0.0001)), evtweight, h_1d, nbins2, h_xmin, h_xmax);

  h_xmin = 0.;
  h_xmax = 1500.;
  plot1D("h_cr1_eSCRaw"+tag_selection          +flav_tag, (sTree->eSCRaw_>0) ? min(sTree->eSCRaw_, float(h_xmax-0.0001)) : max(sTree->eSCRaw_, float(h_xmin+0.0001)), evtweight, h_1d, nbins2, h_xmin, h_xmax);
  plot1D("h_cr1_eSCRaw"+tag_selection+tag_njets+flav_tag, (sTree->eSCRaw_>0) ? min(sTree->eSCRaw_, float(h_xmax-0.0001)) : max(sTree->eSCRaw_, float(h_xmin+0.0001)), evtweight, h_1d, nbins2, h_xmin, h_xmax);

  h_xmin = 0.;
  h_xmax = 200.;
  plot1D("h_cr1_eSCPresh"+tag_selection          +flav_tag, (sTree->eSCPresh_>0) ? min(sTree->eSCPresh_, float(h_xmax-0.0001)) : max(sTree->eSCPresh_, float(h_xmin+0.0001)), evtweight, h_1d, nbins2, h_xmin, h_xmax);
  plot1D("h_cr1_eSCPresh"+tag_selection+tag_njets+flav_tag, (sTree->eSCPresh_>0) ? min(sTree->eSCPresh_, float(h_xmax-0.0001)) : max(sTree->eSCPresh_, float(h_xmin+0.0001)), evtweight, h_1d, nbins2, h_xmin, h_xmax);

  h_xmin = -3.;
  h_xmax = 3.;
  plot1D("h_cr1_etasc1"+tag_selection          +flav_tag, (sTree->etasc1_>0) ? min(sTree->etasc1_, float(h_xmax-0.0001)) : max(sTree->etasc1_, float(h_xmin+0.0001)), evtweight, h_1d, nbins2, h_xmin, h_xmax);
  plot1D("h_cr1_etasc1"+tag_selection+tag_njets+flav_tag, (sTree->etasc1_>0) ? min(sTree->etasc1_, float(h_xmax-0.0001)) : max(sTree->etasc1_, float(h_xmin+0.0001)), evtweight, h_1d, nbins2, h_xmin, h_xmax);

  h_xmin = 0.;
  h_xmax = 0.5;
  plot1D("h_cr1_iso1"+tag_selection          +flav_tag, (sTree->iso1_>0) ? min(sTree->iso1_, float(h_xmax-0.0001)) : max(sTree->iso1_, float(h_xmin+0.0001)), evtweight, h_1d, nbins2, h_xmin, h_xmax);
  plot1D("h_cr1_iso1"+tag_selection+tag_njets+flav_tag, (sTree->iso1_>0) ? min(sTree->iso1_, float(h_xmax-0.0001)) : max(sTree->iso1_, float(h_xmin+0.0001)), evtweight, h_1d, nbins2, h_xmin, h_xmax);

  h_xmin = 0.;
  h_xmax = 0.5;
  plot1D("h_cr1_isont1"+tag_selection          +flav_tag, (sTree->isont1_>0) ? min(sTree->isont1_, float(h_xmax-0.0001)) : max(sTree->isont1_, float(h_xmin+0.0001)), evtweight, h_1d, nbins2, h_xmin, h_xmax);
  plot1D("h_cr1_isont1"+tag_selection+tag_njets+flav_tag, (sTree->isont1_>0) ? min(sTree->isont1_, float(h_xmax-0.0001)) : max(sTree->isont1_, float(h_xmin+0.0001)), evtweight, h_1d, nbins2, h_xmin, h_xmax);
      
  }
 // else cout<<"*****PROBLEM WITH MT******** "<<sTree->t1met10mt_<<endl;
  
  plot1D("h_cr1_mt_count"+tag_selection          +flav_tag, mt_count, evtweight, h_1d, 2, 0, 2);
  plot1D("h_cr1_mt_count"+tag_selection+tag_njets+flav_tag, mt_count, evtweight, h_1d, 2, 0, 2);

}


void StopTreeLooper::makeNJPlots( const StopTree *sTree, float evtweight, std::map<std::string, TH1F*> &h_1d, 
				  string tag_selection, string flav_tag ) 
{

  plot1D("h_njets"    +tag_selection,          min(sTree->npfjets30_,4), evtweight, h_1d, 4,1,5);
  plot1D("h_njets"    +tag_selection+flav_tag, min(sTree->npfjets30_,4), evtweight, h_1d, 4,1,5);
  plot1D("h_njets_all"+tag_selection,          min(sTree->npfjets30_,8), evtweight, h_1d, 7,1,8);
  plot1D("h_njets_all"+tag_selection+flav_tag, min(sTree->npfjets30_,8), evtweight, h_1d, 7,1,8);

}

void StopTreeLooper::makeZPlots( const StopTree *sTree, float evtweight, std::map<std::string, TH1F*> &h_1d, 
				  string tag_selection, string tag_njets, string flav_tag ) 
{

  int nbins = 40;
  float h_xmin = 0.;
  float h_xmax = 400.;
  float x_ovflw = h_xmax-0.001;

  plot1D("h_z_met"+tag_selection+flav_tag, min(sTree->t1met10_,x_ovflw), evtweight, h_1d, nbins, h_xmin, h_xmax);
  plot1D("h_z_met"+tag_selection+tag_njets+flav_tag, min(sTree->t1met10_,x_ovflw), evtweight, h_1d, nbins, h_xmin, h_xmax);

  string lep1type =  abs(sTree->id1_)==13 ? "h_muo" : "h_ele";
  string lep2type =  abs(sTree->id2_)==13 ? "h_muo" : "h_ele";
  plot1D(lep1type+"pt"+tag_selection+flav_tag, min(sTree->lep1_.Pt(),(float)199.99), evtweight, h_1d, 40, 20, 200);
  plot1D(lep2type+"pt"+tag_selection+flav_tag, min(sTree->lep2_.Pt(),(float)199.99), evtweight, h_1d, 40, 20, 200);
  plot1D(lep1type+"pt"+tag_selection+tag_njets+flav_tag, min(sTree->lep1_.Pt(),(float)199.99), evtweight, h_1d, 40, 20, 200);
  plot1D(lep2type+"pt"+tag_selection+tag_njets+flav_tag, min(sTree->lep2_.Pt(),(float)199.99), evtweight, h_1d, 40, 20, 200);
  
  plot1D("h_z_leppt"  +tag_selection+flav_tag, min(sTree->lep1_.Pt(),(float)299.99), evtweight, h_1d, 50, 20., 300.);
  plot1D("h_z_lepeta" +tag_selection+flav_tag, sTree->lep1_.Eta(), evtweight, h_1d, 24, -2.4, 2.4);
  plot1D("h_z_lep2pt" +tag_selection+flav_tag, min(sTree->lep2_.Pt(),(float)199.99), evtweight, h_1d, 50, 20., 200.);
  plot1D("h_z_lep2eta"+tag_selection+flav_tag, sTree->lep2_.Eta(), evtweight, h_1d, 24, -2.4, 2.4);
  plot1D("h_z_leppt"  +tag_selection+tag_njets+flav_tag, min(sTree->lep1_.Pt(),(float)299.99), evtweight, h_1d, 50, 20., 300.);
  plot1D("h_z_lepeta" +tag_selection+tag_njets+flav_tag, sTree->lep1_.Eta(), evtweight, h_1d, 24, -2.4, 2.4);
  plot1D("h_z_lep2pt" +tag_selection+tag_njets+flav_tag, min(sTree->lep2_.Pt(),(float)199.99), evtweight, h_1d, 50, 20., 200.);
  plot1D("h_z_lep2eta"+tag_selection+tag_njets+flav_tag, sTree->lep2_.Eta(), evtweight, h_1d, 24, -2.4, 2.4);

  float dphi_metlep = getdphi(sTree->lep1_.Phi(), sTree->t1met10phi_);
  plot1D("h_z_dphi_metl"+tag_selection+flav_tag, dphi_metlep, evtweight, h_1d, 15, 0., 3.14159);
  plot1D("h_z_dphi_metl"+tag_selection+tag_njets+flav_tag, dphi_metlep, evtweight, h_1d, 15, 0., 3.14159);
  plot1D("h_z_mt"+tag_selection+flav_tag, min(sTree->t1met10mt_, x_ovflw), evtweight, h_1d, nbins, h_xmin, h_xmax);
  plot1D("h_z_mt"+tag_selection+tag_njets+flav_tag, min(sTree->t1met10mt_, x_ovflw), evtweight, h_1d, nbins, h_xmin, h_xmax);

  if ( sTree->npfjets30_<1 ) return;
  plot1D("h_z_j1pt" +tag_selection+flav_tag, min(sTree->pfjet1_.Pt(), (float)399.99), evtweight, h_1d, 20, 30., 400.);
  plot1D("h_z_j1eta"+tag_selection+flav_tag, sTree->pfjet1_.Eta(), evtweight, h_1d, 24, -2.4, 2.4);
  if ( sTree->npfjets30_<2 ) return;
  plot1D("h_z_j2pt" +tag_selection+flav_tag, min(sTree->pfjet2_.Pt(), (float)299.99), evtweight, h_1d, 20, 30., 300.);
  plot1D("h_z_j2eta"+tag_selection+flav_tag, sTree->pfjet2_.Eta(), evtweight, h_1d, 24, -2.4, 2.4);
  if ( sTree->npfjets30_<3 ) return;
  plot1D("h_z_j3pt" +tag_selection+flav_tag, min(sTree->pfjet3_.Pt(), (float)199.99), evtweight, h_1d, 20, 30., 200.);
  plot1D("h_z_j3eta"+tag_selection+flav_tag, sTree->pfjet3_.Eta(), evtweight, h_1d, 24, -2.4, 2.4);
  if ( sTree->npfjets30_<4 ) return;
  plot1D("h_z_j4pt" +tag_selection+flav_tag, min(sTree->pfjet4_.Pt(), (float)119.99), evtweight, h_1d, 20, 30., 120.);
  plot1D("h_z_j4eta"+tag_selection+flav_tag, sTree->pfjet4_.Eta(), evtweight, h_1d, 24, -2.4, 2.4);
  
}
