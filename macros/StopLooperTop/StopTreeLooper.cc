#include "../../CORE/jetSmearingTools.h"

#include "StopTreeLooper.h"
#include "../Plotting/PlotUtilities.h"
#include "../Core/MT2Utility.h"
#include "../Core/mt2bl_bisect.h"
#include "../Core/mt2w_bisect.h"

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
  t1metphicorr = -9999.;
  t1metphicorrphi = -9999.;
  t1metphicorrmt = -9999.;
  min_mtpeak = -9999.;
  max_mtpeak = -9999.; 
}

StopTreeLooper::~StopTreeLooper()
{
}

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

bool is_badLaserEvent (const DorkyEventIdentifier &id) {
  if (events_lasercalib.find(id) != events_lasercalib.end()) return true;
  return false;
}



bool compare_candidates( Candidate x, Candidate y ){
  return x.chi2 < y.chi2;
}

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


double fchi2 (double c1, double pt1, double sigma1, double pt2, double sigma2,
              double m12, double m22, double m02){
    double rat1 = pt1 * (1 - c1) / sigma1;
    double rat2 = pt2 * (1 - fc2(c1, m12, m22, m02)) / sigma2;

    return ( rat1 * rat1 + rat2 * rat2);
}

//void StopSelector::minuitFunction(int& npar, double *gout, double &result, double par[], int flg)
void minuitFunction(int&, double* , double &result, double par[], int){
  result=fchi2(par[0], par[1], par[2], par[3], par[4], par[5], par[6], par[7]);
}

float getDataMCRatio(float eta){
  if (eta >=0.0 && eta < 0.5) return 1.052;
  if (eta >=0.5 && eta < 1.1) return 1.057;
  if (eta >=1.1 && eta < 1.7) return 1.096;
  if (eta >=1.7 && eta < 2.3) return 1.134;
  if (eta >=2.3 && eta < 5.0) return 1.288;
  return 1.0;
}


/* Reconstruct the hadronic top candidates, select the best candidate and
 * store the chi2 =  (m_jj - m_W)^2/sigma_m_jj + (m_jjj - m_t)^2/sigma_m_jjj
 * return the number of candidates found.
 *
 * n_jets - number of jets.
 * jets - jets
 * btag - b-tagging information of the jets
 * mc - qgjet montecarlo match number for the jets
 * 
 * returns a list of candidates sorted by chi2 ( if __sort = true in .h );
*/ 

list<Candidate> StopTreeLooper::recoHadronicTop(StopTree* tree, bool isData){

  assert( tree->pfjets_->size() == tree->pfjets_csv_.size() );

  static JetSmearer* jetSmearer = 0;
  if (jetSmearer == 0 ){
    std::vector<std::string> list_of_file_names;
    list_of_file_names.push_back("../../CORE/jetsmear/data/Spring10_PtResolution_AK5PF.txt");
    list_of_file_names.push_back("../../CORE/jetsmear/data/Spring10_PhiResolution_AK5PF.txt");
    list_of_file_names.push_back("../../COREjetsmear/data/jet_resolutions.txt");
    jetSmearer = makeJetSmearer(list_of_file_names);
  }

  LorentzVector* lep = &tree->lep1_;
  double met = tree->t1metphicorr_;
  double metphi = tree->t1metphicorrphi_;

  vector<LorentzVector> jets;
  vector<float> btag;

  //cout << endl << "baby branches:" << endl;
  for( unsigned int i = 0 ; i < tree->pfjets_->size() ; ++i ){
    //cout << i << " pt csv " << tree->pfjets_->at(i).pt() << " " << tree->pfjets_csv_.at(i) << endl;
    jets.push_back( tree->pfjets_->at(i) );
    btag.push_back( tree->pfjets_csv_.at(i) );
  } 

  // cout << endl << "stored:" << endl;
  // for( unsigned int i = 0 ; i < jets.size() ; ++i ){
  //   cout << i << " pt csv " << jets.at(i).pt() << " " << btag.at(i) << endl;
  // } 

  assert( jets.size() == btag.size() );

  float metx = met * cos( metphi );
  float mety = met * sin( metphi );

  int n_jets = jets.size();
  double sigma_jets[n_jets];
  for (int i=0; i<n_jets; ++i)
    sigma_jets[i] = getJetResolution(jets[i], jetSmearer);

  if ( isData )
    for (int i=0; i<n_jets; ++i)
      sigma_jets[i] *= getDataMCRatio(jets[i].eta());

  /*
    // TO BE UPDATED
    // we don't currently store the qgjet variables, to be added

  vector<int> mc;

  if (!isData) {
    mc.push_back( tree->qgjet1_ );
    mc.push_back( tree->qgjet2_ );
    mc.push_back( tree->qgjet3_ );
    mc.push_back( tree->qgjet4_ );
    mc.push_back( tree->qgjet5_ );
    mc.push_back( tree->qgjet6_ );
  }

  int ibl[5];
  int iw1[5];
  int iw2[5];
  int ib[5];

  if ( !isData ){
     
     // Matching MC algoritm search over all conbinations  until the 
     // right combination is found. More than one candidate is suported 
     //  but later only the first is used.
     // 
    int match = 0;
    for (int jbl=0; jbl<n_jets; ++jbl )
      for (int jb=0; jb<n_jets; ++jb )
        for (int jw1=0; jw1<n_jets; ++jw1 )
          for (int jw2=jw1+1; jw2<n_jets; ++jw2 )
            if ( (mc.at(jw2)==2 && mc.at(jw1)==2 && mc.at(jb)==1 && mc.at(jbl)==-1) ||
                 (mc.at(jw2)==-2 && mc.at(jw1)==-2 && mc.at(jb)==-1 && mc.at(jbl)==1) ) {
                    ibl[match] = jbl;
                    iw1[match] = jw1;
                    iw2[match] = jw2;
                    ib[match] = jb;
                    match++;
            }
  }
  */
  
////////    * Combinatorics. j_1 Pt must be > PTMIN_W1 and so on.
  
  vector<int> v_i, v_j;
  vector<double> v_k1, v_k2;
  for ( int i=0; i<n_jets; ++i )
    for ( int j=i+1; j<n_jets; ++j ){
      double pt_w1 = jets[i].Pt();
      double pt_w2 = jets[j].Pt();
      if (pt_w1 < PTMIN_J1  || pt_w2 < PTMIN_J2)
        continue;

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
      double c2 = fc2(c1, jets[i].mass2(), jets[j].mass2(), hadW.mass2());
                
      delete minimizer;

     
  //     * W Mass check :)
  //     *  Never trust a computer you can't throw out a window. 
 //      *  - Steve Wozniak 

//      cout << "c1 = " <<  c1 << "  c1 = " << c2 << "   M_jj = " 
//           << ((jets[i] * c1) + (jets[j] * c2)).mass() << endl;
      
      v_i.push_back(i);
      v_j.push_back(j);
      v_k1.push_back(c1);
      v_k2.push_back(c2);
    }


  list<Candidate> chi2candidates;
        
  mt2_bisect::mt2 mt2_event;
  mt2bl_bisect::mt2bl mt2bl_event;
  mt2w_bisect::mt2w mt2w_event;
  
  for ( int b=0; b<n_jets; ++b )
    for (int o=0; o<n_jets; ++o){
      if ( b == o )
        continue;

      if ( btag[b] < BTAG_MIN && btag[o] < BTAG_MIN )
        continue;

      double pt_b = jets[b].Pt();
      if ( btag[b] >= BTAG_MIN && pt_b < PTMIN_BTAG )
        continue;

      if ( btag[b] < BTAG_MIN && pt_b < PTMIN_B )
        continue;

      double pt_o = jets[o].Pt();
      if ( btag[o] >= BTAG_MIN && pt_o < PTMIN_OTAG )
        continue;

      if ( btag[o] < BTAG_MIN && pt_o < PTMIN_O)
        continue;

      ///
      //  MT2 Variables
      ///
         
      double pl[4];     // Visible lepton
      double pb1[4];    // bottom on the same side as the visible lepton
      double pb2[4];    // other bottom, paired with the invisible W
      double pmiss[3];  // <unused>, pmx, pmy   missing pT
      pl[0]= lep->E(); pl[1]= lep->Px(); pl[2]= lep->Py(); pl[3]= lep->Pz();
      pb1[0] = jets[o].E();  pb1[1] = jets[o].Px(); 
      pb1[2] = jets[o].Py(); pb1[3] = jets[o].Pz();
      pb2[0] = jets[b].E();  pb2[1] = jets[b].Px(); 
      pb2[2] = jets[b].Py(); pb2[3] = jets[b].Pz();
      pmiss[0] = 0.; pmiss[1] = metx; pmiss[2] = mety;

      double pmiss_lep[3];
      pmiss_lep[0] = 0.;
      pmiss_lep[1] = pmiss[1]+pl[1]; pmiss_lep[2] = pmiss[2]+pl[2];

      mt2_event.set_momenta( pb1, pb2, pmiss_lep );
      mt2_event.set_mn( 0.0 );   // Invisible particle mass
      double c_mt2b = mt2_event.get_mt2();
         
      mt2bl_event.set_momenta(pl, pb1, pb2, pmiss); 
      double c_mt2bl = mt2bl_event.get_mt2bl();

      mt2w_event.set_momenta(pl, pb1, pb2, pmiss); 
      double c_mt2w = mt2w_event.get_mt2w();

//      cout << b << ":"<< btag[b] << " - " << o << ":" << btag[o] << " = " << c_mt2w << endl;

      for (unsigned int w = 0; w < v_i.size() ; ++w ){
        int i = v_i[w];
        int j = v_j[w];
        if ( i==o || i==b || j==o || j==b )
            continue;

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
        double sigma_w2 = pt_w1*sigma_jets[i] * pt_w1*sigma_jets[i]
                        + pt_w2*sigma_jets[j] * pt_w2*sigma_jets[j];
        double smw2 = (1.+2.*pt_w*pt_w/massW/massW)*sigma_w2;
        double pt_t = hadT.Pt();
        double sigma_t2 = c1*pt_w1*sigma_jets[i] * c1*pt_w1*sigma_jets[i]
                        + c2*pt_w2*sigma_jets[j] * c2*pt_w2*sigma_jets[j]
                        + pt_b*sigma_jets[b] * pt_b*sigma_jets[b];
        double smtop2 = (1.+2.*pt_t*pt_t/massT/massT)*sigma_t2;

        double c_chi2 = (massT-PDG_TOP_MASS)*(massT-PDG_TOP_MASS)/smtop2
                      + (massW-PDG_W_MASS)*(massW-PDG_W_MASS)/smw2;

        //bool c_match = ( !isData &&  iw1[0]==i && iw2[0]==j && ib[0]==b && ibl[0]==o );
        bool c_match = false;

        Candidate c;
        c.chi2  = c_chi2;
        c.mt2b  = c_mt2b;
        c.mt2w  = c_mt2w;
        c.mt2bl = c_mt2bl;
        c.j1 = i;
        c.j2 = j;
        c.bi = b;
        c.oi = o;
        c.k1 = c1;
        c.k2 = c2;
        c.match = c_match;

        chi2candidates.push_back(c);
      }
    }

   if (__SORT) 
     chi2candidates.sort(compare_candidates);

   return chi2candidates;
}


//--------------------------------------------------------------------

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
  std::map<std::string, TH1F*> h_1d;
  //also for control regions
  std::map<std::string, TH1F*> h_1d_cr1, h_1d_cr2, h_1d_cr4, h_1d_cr5;
  //for signal region 
  std::map<std::string, TH1F*> h_1d_sig;
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

      float evtweight = isData ? 1. : ( tree->weight_ * 9.708 * tree->nvtxweight_ * tree->mgcor_ );
      // to reweight from file - also need to comment stuff before
      //      float vtxweight = vtxweight_n( tree->nvtx_, h_vtx_wgt, isData );

      plot1D("h_vtx",       tree->nvtx_,       evtweight, h_1d_z, 40, 0, 40);
      plot1D("h_vtxweight", tree->nvtxweight_, evtweight, h_1d_z, 41, -4., 4.);

      // 
      // selection criteria
      // 

      //preselection
      if ( !passEvtSelection(tree, name) ) continue;

      list<Candidate> candidates = recoHadronicTop(tree, isData );

      if (candidates.size() < 1 ) continue;    

      //----------------------------------------------------------------------------------------
      // calculate the hadronic top and MT2 variables and compare to values stored in babies
      //----------------------------------------------------------------------------------------

      cout << endl << endl << endl;
      cout << "On-the-fly : " << candidates.size()         << " candidates" << endl;
      cout << "From baby  : " << tree->candidates_->size() << " candidates" << endl << endl;

      list<Candidate>::iterator candIter;

      int i = 0;

      for(candIter = candidates.begin() ; candIter != candidates.end() ; candIter++ ){
	cout << "chi2        " << (*candIter).chi2   << " " << (tree->candidates_->at(i)).chi2  << endl;
	cout << "mt2w        " << (*candIter).mt2w   << " " << (tree->candidates_->at(i)).mt2w  << endl;
	cout << "mt2bl       " << (*candIter).mt2bl  << " " << (tree->candidates_->at(i)).mt2bl << endl;
	cout << "mt2b        " << (*candIter).mt2b   << " " << (tree->candidates_->at(i)).mt2b  << endl;
	cout << endl;
	++i;
      }


      plot1D("h_chi2",  candidates.front().chi2,       evtweight, h_1d_z, 100, 0, 20);

    } // end event loop

    // delete tree;

  } // end file loop

    //
    // finish
    //

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
