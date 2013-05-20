#include <algorithm>
#include <iostream>
#include <vector>
#include <math.h>
#include <fstream>
#include "TCanvas.h"
#include "TLegend.h"
#include "TChain.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TROOT.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"
#include "TPad.h"
#include "TCut.h"
#include "TProfile.h"
#include "THStack.h"
#include "TLatex.h"
#include "TGraphErrors.h"
#include "TStyle.h"
#include "TLine.h"
#include "TMath.h"
#include <sstream>
#include <iomanip>

using namespace std;

double zbi(double n_on, double mu_b_hat, double sigma_b, bool verbose = false ){
 
  //double n_on     = 140.;                         // total events in signal region (S+B)
  //double mu_b_hat = 83.33;                        // mean num of BG events expected in sig. region
  //double sigma_b  = 8.333;                        // uncertainty of mu_b_hat
  
  double tau      = mu_b_hat / (sigma_b*sigma_b); // scale factor to corresp. Noff/Non              
  double n_off    = tau*mu_b_hat;
  double P_Bi     = TMath::BetaIncomplete(1./(1.+tau), n_on, n_off+1);
  double Z_Bi     = sqrt(2.)*TMath::ErfInverse(1 - 2.*P_Bi);           

  if( verbose ){
    cout  <<"  total events in signal region (S+B)               - n_on     " <<n_on      <<endl
	  <<"  mean num of BG events expected in sig. region     - mu_b_hat " <<mu_b_hat  <<endl
	  <<"  uncertainty of mu_b_hat                           - sigma_b  " <<sigma_b   <<endl
	  <<"  scale factor to corresp. Noff/Non                 - tau      " <<tau       <<endl
	  <<"  tau*mu_b_hat                                      - n_off    " <<n_off     <<endl
	  <<"  TMath::BetaIncomplete(1./(1.+tau), n_on, n_off+1) - P_Bi     " <<P_Bi      <<endl
	  <<"  sqrt(2.)*TMath::ErfInverse(1 - 2.*P_Bi)           - Z_Bi     " <<Z_Bi      <<endl;
  }

  return Z_Bi;
}

double simple_significance( double nsig , double nbkg , double sbkg ){
  float toterr       = sqrt( nbkg + sbkg * sbkg );
  float significance = nsig / toterr;
  return significance;
}

void SMS_projection(char* sample = "T2tt" , int x = 1, bool doBDT = true , char* pol = "" , bool print = false){

  gStyle->SetPaintTextFormat(".2f");

  //--------------------------------------------------
  // input parameters
  //--------------------------------------------------

  const float lumi      = 19500;
  char* suffix          = (char*) "";
  char* denomhistoname  = (char*) "masses";
  char* filename        = (char*) "";
  char* denomname       = (char*) "";
  char* BDTchar         = (char*) "";
  char* xchar           = (char*) "";
  char* label           = (char*)"";
  bool  doISRWeight     = true;

  if( doBDT ){
    BDTchar  = (char*) "_BDT";
    cout << "Doing BDT analysis" << endl;
  }else{
    cout << "Doing cut-based analysis" << endl;
  }

  //--------------------------------------------------
  // set up input files
  //--------------------------------------------------

  if( TString(sample).Contains("T2tt") ){    
    filename  = (char*) "/tas/cms2/stop/cms2V05-03-26_stoplooperV00-02-24/T2tt_mad/minibaby_V00-03-06/Skim_4jets_MET100_MT120/merged*root";
    denomname = (char*) "/tas/cms2/stop/cms2V05-03-26_stoplooperV00-02-24/T2tt_mad/minibaby_V00-03-06/Skim_4jets_MET100_MT120/myMassDB_T2tt_combined_25GeVbins.root";

    label     = (char*)"pp #rightarrow #tilde{t}#tilde{t}, #tilde{t} #rightarrow t #tilde{#chi}_{1}^{0}";
  }

  else if( TString(sample).Contains("T2bw_MG") ){

    if( x==25 ){
      xchar          = (char*) "_x25";
      denomhistoname = (char*) "masses25";
      filename       = (char*) "/tas/cms2/stop/cms2V05-03-26_stoplooperV00-02-25/T2bw_mad/minibaby_V00-03-10/Skim_4jets_MET100_MT120/merged*x025*root";
      denomname      = (char*) "/tas/cms2/stop/cms2V05-03-26_stoplooperV00-02-25/T2bw_mad/minibaby_V00-03-10/Skim_4jets_MET100_MT120/myMassDB_mStop_x25_25GeVbins.root";
      label          = (char*) "pp #rightarrow #tilde{t}#tilde{t}, #tilde{t} #rightarrow b+#tilde{#chi}_{1}^{#pm}, x=0.25";
    }

    else if( x==50 ){
      xchar          = (char*) "_x50";
      denomhistoname = (char*) "masses";
      filename       = (char*) "/tas/cms2/stop/cms2V05-03-26_stoplooperV00-02-25/T2bw_mad/minibaby_V00-03-10/Skim_4jets_MET100_MT120/merged*x050*root";
      denomname      = (char*) "/tas/cms2/stop/cms2V05-03-26_stoplooperV00-02-25/T2bw_mad/minibaby_V00-03-10/Skim_4jets_MET100_MT120/myMassDB_mStop_x50_25GeVbins.root";
      label          = (char*) "pp #rightarrow #tilde{t}#tilde{t}, #tilde{t} #rightarrow b+#tilde{#chi}_{1}^{#pm}, x=0.50";
    }

    else if( x==75 ){
      xchar          = (char*) "_x75";
      denomhistoname = (char*) "masses75";
      filename       = (char*) "/tas/cms2/stop/cms2V05-03-26_stoplooperV00-02-25/T2bw_mad/minibaby_V00-03-10/Skim_4jets_MET100_MT120/merged*x075*root";
      denomname      = (char*) "/tas/cms2/stop/cms2V05-03-26_stoplooperV00-02-25/T2bw_mad/minibaby_V00-03-10/Skim_4jets_MET100_MT120/myMassDB_mStop_x75_25GeVbins.root";
      label          = (char*) "pp #rightarrow #tilde{t}#tilde{t}, #tilde{t} #rightarrow b+#tilde{#chi}_{1}^{#pm}, x=0.75";
    }
    
    else{
      cout << "ERROR! unrecognized x value " << x << ", quitting!!!" << endl;
    }
  }

  else if( TString(sample).Contains("T2bw") ){

    denomname      = (char*) "/tas/cms2/stop/cms2V05-03-26_stoplooperV00-02-25/T2bw_pythiaCoarse/minibaby_V00-03-09/Skim_4jets_MET100_MT120/myMassDB_T2bw_coarse_25GeVbins.root";
    doISRWeight    = false;

    if( x==25 ){
      xchar          = (char*) "_x25";
      denomhistoname = (char*) "masses25";
      filename       = (char*) "/tas/cms2/stop/cms2V05-03-26_stoplooperV00-02-25/T2bw_pythiaCoarse/minibaby_V00-03-09/Skim_4jets_MET100_MT120/T2bw_x25.root";
      label          = (char*) "pp #rightarrow #tilde{t}#tilde{t}, #tilde{t} #rightarrow b+#tilde{#chi}_{1}^{#pm}, x=0.25";
    }

    else if( x==50 ){
      xchar          = (char*) "_x50";
      denomhistoname = (char*) "masses";
      filename       = (char*) "/tas/cms2/stop/cms2V05-03-26_stoplooperV00-02-25/T2bw_pythiaCoarse/minibaby_V00-03-09/Skim_4jets_MET100_MT120/T2bw_x50.root";
      label          = (char*) "pp #rightarrow #tilde{t}#tilde{t}, #tilde{t} #rightarrow b+#tilde{#chi}_{1}^{#pm}, x=0.50";
    }

    else if( x==75 ){
      xchar          = (char*) "_x75";
      denomhistoname = (char*) "masses75";
      filename       = (char*) "/tas/cms2/stop/cms2V05-03-26_stoplooperV00-02-25/T2bw_pythiaCoarse/minibaby_V00-03-09/Skim_4jets_MET100_MT120/T2bw_x75.root";
      label          = (char*) "pp #rightarrow #tilde{t}#tilde{t}, #tilde{t} #rightarrow b+#tilde{#chi}_{1}^{#pm}, x=0.75";
    }

    else{
      cout << "ERROR! unrecognized x value " << x << ", quitting!!!" << endl;
    }
  }

 
  else{
    cout << "ERROR! unrecognized sample " << sample << ", quitting!!!" << endl;
    exit(0);
  }

  //--------------------------------------------------
  // set up logfile
  //--------------------------------------------------

  char* logfilename = Form("logfiles/%s%s%s%s%s_projections.log",sample,xchar,suffix,pol,BDTchar);
  ofstream* logfile = new ofstream();
  logfile->open(logfilename,ios::trunc);

  cout << "----------------------------------------------------" << endl;
  cout << "Sample            " << sample          << endl;
  cout << "logfile           " << logfilename     << endl;
  cout << "x                 " << x               << endl;
  cout << "pol               " << pol             << endl;
  cout << "Using file        " << filename        << endl;
  cout << "Using denominator " << denomname       << endl;
  cout << "Denom histo       " << denomhistoname  << endl;
  cout << "Using lumi        " << lumi            << " pb-1" << endl;
  cout << "----------------------------------------------------" << endl;

  *logfile << "----------------------------------------------------" << endl;
  *logfile << "Sample            " << sample          << endl;
  *logfile << "logfile           " << logfilename     << endl;
  *logfile << "x                 " << x               << endl;
  *logfile << "pol               " << pol             << endl;
  *logfile << "Using file        " << filename        << endl;
  *logfile << "Using denominator " << denomname       << endl;
  *logfile << "Denom histo       " << denomhistoname  << endl;
  *logfile << "Using lumi        " << lumi            << " pb-1" << endl;
  *logfile << "----------------------------------------------------" << endl;

  TFile* fdenom = TFile::Open(denomname);
  TH2F*  hdenom = (TH2F*) fdenom->Get(denomhistoname);

  //--------------------------------------------------
  // read in TChain
  //--------------------------------------------------

  TChain *ch = new TChain("t");
  ch->Add(filename);

  //--------------------------------------------------
  // read in reference cross section
  //--------------------------------------------------

  TFile *xsecfile = TFile::Open("stop_xsec.root");
  TH1F* refxsec   = (TH1F*) xsecfile->Get("h_stop_xsec");

  //--------------------------------------------------
  // preselection
  //--------------------------------------------------

  TCut rho("rhovor>0 && rhovor<40");
  TCut filters("isdata==0 || (csc==0 && hbhe==1 && hcallaser==1 && ecaltp==1 && trkfail==1 && eebadsc==1 && hbhenew==1)");
  TCut goodlep("ngoodlep > 0 && abs( pflep1.Pt() - lep1.Pt() ) < 10.0 && abs(isopf1 * lep1.Pt() ) < 5.0");
  TCut el("leptype==0 && abs(lep1->Eta())<1.4442 && lep1->Pt()>30.0 && eoverpin < 4.0 && (isdata==0 || ele27wp80==1)");
  TCut mu("leptype==1 && abs(lep1->Eta())<2.1    && lep1->Pt()>25.0 && (isdata==0 || isomu24==1)");
  TCut njets4("mini_njets >= 4");
  TCut btag1("mini_nb >= 1");
  TCut passisotrk("mini_passisotrk==1");
  TCut tauveto("mini_passtauveto == 1");
  TCut mt120(" mini_mt > 120");
  TCut mt150(" mini_mt > 150");
  TCut dphi("mini_dphimjmin>0.8");
  TCut chi2("mini_chi2<5.0");
  TCut mt2w("mini_mt2w>200.0");
  TCut bpt100("mini_pt_b > 100.0");
  TCut met100("mini_met > 100");
  TCut met150("mini_met > 150.0");
  TCut met200("mini_met > 200.0");
  TCut met250("mini_met > 250.0");
  TCut met300("mini_met > 300.0");
  TCut testing("event%2==0");
  TCut isrweight("mini_isrweight");
  TCut x25("x < 0.3");
  TCut x50("x > 0.3 && x < 0.7");
  TCut x75("x > 0.7");

  TCut presel;

  //-------------------------------------------
  // THESE CUTS DEFINE PRESELECTION REGION
  //-------------------------------------------

  //presel += filters;
  presel += rho;
  presel += goodlep;
  presel += (el||mu);
  presel += njets4;
  presel += btag1;
  presel += passisotrk;
  presel += tauveto;
  presel += met100;
  presel += mt120;
  if( doBDT ){
    presel += testing;
  }
  if( !doISRWeight ){
    isrweight = TCut("1");
  }

  if( TString(sample).Contains("T2bw") && x==25 ) presel += x25;
  if( TString(sample).Contains("T2bw") && x==50 ) presel += x50;
  if( TString(sample).Contains("T2bw") && x==75 ) presel += x75;

  int BDTweight = 1;
  if( doBDT ) BDTweight = 2;

  TCut weight(Form("mini_sltrigeff * %i * mini_weight",BDTweight)); 
  if( TString(pol).Contains("left") )  weight = TCut(Form("mini_sltrigeff * weightleft  * %i",BDTweight));
  if( TString(pol).Contains("right") ) weight = TCut(Form("mini_sltrigeff * weightright * %i",BDTweight));

  if( TString(pol).Contains("T2BW_LR") ) weight = TCut(Form("mini_sltrigeff * mini_t2bwweight_lr * %i",BDTweight));
  if( TString(pol).Contains("T2BW_LS") ) weight = TCut(Form("mini_sltrigeff * mini_t2bwweight_ls * %i",BDTweight));
  if( TString(pol).Contains("T2BW_LL") ) weight = TCut(Form("mini_sltrigeff * mini_t2bwweight_ll * %i",BDTweight));

  if( TString(pol).Contains("T2BW_SR") ) weight = TCut(Form("mini_sltrigeff * mini_t2bwweight_sr * %i",BDTweight));
  if( TString(pol).Contains("T2BW_SS") ) weight = TCut(Form("mini_sltrigeff * mini_t2bwweight_ss * %i",BDTweight));
  if( TString(pol).Contains("T2BW_SL") ) weight = TCut(Form("mini_sltrigeff * mini_t2bwweight_sl * %i",BDTweight));

  if( TString(pol).Contains("T2BW_RR") ) weight = TCut(Form("mini_sltrigeff * mini_t2bwweight_rr * %i",BDTweight));
  if( TString(pol).Contains("T2BW_RS") ) weight = TCut(Form("mini_sltrigeff * mini_t2bwweight_rs * %i",BDTweight));
  if( TString(pol).Contains("T2BW_RL") ) weight = TCut(Form("mini_sltrigeff * mini_t2bwweight_rl * %i",BDTweight));

  cout << "Using pre-selection   " << presel.GetTitle()    << endl;
  cout << "Using weight          " << weight.GetTitle()    << endl;
  cout << "Using ISR weight      " << isrweight.GetTitle() << endl;

  *logfile << "Using pre-selection   " << presel.GetTitle()    << endl;
  *logfile << "Using weight          " << weight.GetTitle()    << endl;
  *logfile << "Using ISR weight      " << isrweight.GetTitle() << endl;

  //--------------------------------------------------
  // signal regions
  //--------------------------------------------------

  vector<TCut>    sigcuts;
  vector<string>  signames;
  vector<string>  labels;
  vector<float>   uls;

  vector<float>   nbkg;
  vector<float>   sbkg;

  if( doBDT ){

    //-----------------------------
    // T2tt BDT
    //-----------------------------

    if( TString(sample).Contains("T2tt") ){

      cout << "Doing T2tt BDT signal regions" << endl;  

      sigcuts.push_back(TCut(presel+"mini_bdt[1] > 0.30"));  signames.push_back("T2TT_BDT1L");  labels.push_back("T2TT_BDT1L");  uls.push_back(-1.0); nbkg.push_back(763.0); sbkg.push_back(102.0);
      sigcuts.push_back(TCut(presel+"mini_bdt[1] > 0.40"));  signames.push_back("T2TT_BDT1T");  labels.push_back("T2TT_BDT1T");  uls.push_back(-1.0); nbkg.push_back(124.0); sbkg.push_back( 21.0);
      sigcuts.push_back(TCut(presel+"mini_bdt[2] > 0.55"));  signames.push_back("T2TT_BDT2");   labels.push_back("T2TT_BDT2");   uls.push_back(-1.0); nbkg.push_back( 85.0); sbkg.push_back( 16.0);
      sigcuts.push_back(TCut(presel+"mini_bdt[3] > 0.65"));  signames.push_back("T2TT_BDT3");   labels.push_back("T2TT_BDT3");   uls.push_back(-1.0); nbkg.push_back( 13.0); sbkg.push_back(  4.0);
      sigcuts.push_back(TCut(presel+"mini_bdt[4] > 0.50"));  signames.push_back("T2TT_BDT4");   labels.push_back("T2TT_BDT4");   uls.push_back(-1.0); nbkg.push_back(  2.9); sbkg.push_back(  1.1);
      sigcuts.push_back(TCut(presel+"mini_bdt[5] > 0.30"));  signames.push_back("T2TT_BDT5");   labels.push_back("T2TT_BDT5");   uls.push_back(-1.0); nbkg.push_back( 87.0); sbkg.push_back( 18.0);
    }

    //-----------------------------
    // T2bw BDT
    //-----------------------------

    else if( TString(sample).Contains("T2bw") ){
      cout << "NOT SETUP FOR T2BW" << endl;
      exit(0);

      cout << "Doing T2bw BDT signal regions" << endl;  

      if( x==25 ){
      	sigcuts.push_back(TCut(presel+"mini_bdt[8]  > 0.30"));  signames.push_back("T2BW25_BDT2");   labels.push_back("T2BW25_BDT2");   uls.push_back(8.73);
      	sigcuts.push_back(TCut(presel+"mini_bdt[9]  > 0.45"));  signames.push_back("T2BW25_BDT3");   labels.push_back("T2BW25_BDT3");   uls.push_back(5.97);
      }

      else if( x==50 ){
      	sigcuts.push_back(TCut(presel+"mini_bdt[12] > 0.35"));  signames.push_back("T2BW50_BDT1");   labels.push_back("T2BW50_BDT1");   uls.push_back(28.3);
      	sigcuts.push_back(TCut(presel+"mini_bdt[13] > 0.45"));  signames.push_back("T2BW50_BDT2L");  labels.push_back("T2BW50_BDT2L");  uls.push_back(21.8);
      	sigcuts.push_back(TCut(presel+"mini_bdt[13] > 0.55"));  signames.push_back("T2BW50_BDT2T");  labels.push_back("T2BW50_BDT2T");  uls.push_back(10.4);
      	sigcuts.push_back(TCut(presel+"mini_bdt[14] > 0.35"));  signames.push_back("T2BW50_BDT3");   labels.push_back("T2BW50_BDT3");   uls.push_back(11.5);
      }

      else if( x==75 ){
      	sigcuts.push_back(TCut(presel+"mini_bdt[17] > 0.30"));  signames.push_back("T2BW75_BDT1");   labels.push_back("T2BW75_BDT1");   uls.push_back(24.2);
      	sigcuts.push_back(TCut(presel+"mini_bdt[18] > 0.55"));  signames.push_back("T2BW75_BDT2");   labels.push_back("T2BW75_BDT2");   uls.push_back(14.4);
      	sigcuts.push_back(TCut(presel+"mini_bdt[19] > 0.50"));  signames.push_back("T2BW75_BDT3");   labels.push_back("T2BW75_BDT3");   uls.push_back(7.96);
      	sigcuts.push_back(TCut(presel+"mini_bdt[20] > 0.25"));  signames.push_back("T2BW75_BDT4");   labels.push_back("T2BW75_BDT4");   uls.push_back(129.7);
      }

    }

  }

  else{

    cout << "NOT SETUP FOR CUT-BASED LIMITS" << endl;
    exit(0);

    //-----------------------------
    // T2tt cut-and-count
    //-----------------------------

    if( TString(sample).Contains("T2tt") ){

      cout << "Doing T2tt cut-based signal regions" << endl;  

      // // low-mass
      sigcuts.push_back(presel+met150+dphi+chi2);        signames.push_back("T2TT_LM150");   labels.push_back("T2TT_LM150");   uls.push_back(79.6);
      sigcuts.push_back(presel+met200+dphi+chi2);        signames.push_back("T2TT_LM200");   labels.push_back("T2TT_LM200");   uls.push_back(42.1);
      sigcuts.push_back(presel+met250+dphi+chi2);        signames.push_back("T2TT_LM250");   labels.push_back("T2TT_LM250");   uls.push_back(19.0);
      sigcuts.push_back(presel+met300+dphi+chi2);        signames.push_back("T2TT_LM300");   labels.push_back("T2TT_LM300");   uls.push_back( 9.9);

      // // high-mass
      sigcuts.push_back(presel+met150+dphi+chi2+mt2w);   signames.push_back("T2TT_HM150");   labels.push_back("T2TT_HM150");   uls.push_back(17.5);
      sigcuts.push_back(presel+met200+dphi+chi2+mt2w);   signames.push_back("T2TT_HM200");   labels.push_back("T2TT_HM200");   uls.push_back(12.3);
      sigcuts.push_back(presel+met250+dphi+chi2+mt2w);   signames.push_back("T2TT_HM250");   labels.push_back("T2TT_HM250");   uls.push_back( 9.0);
      sigcuts.push_back(presel+met300+dphi+chi2+mt2w);   signames.push_back("T2TT_HM300");   labels.push_back("T2TT_HM300");   uls.push_back( 6.4);

    }

    //-----------------------------
    // T2bw cut-and-count
    //-----------------------------

    else if( TString(sample).Contains("T2bw") ){

      cout << "Doing T2bw cut-based signal regions" << endl;  

      // low-mass
      sigcuts.push_back(presel+met100+dphi);               signames.push_back("T2BW_LM100");   labels.push_back("T2BW_LM100");   uls.push_back(79.6);
      sigcuts.push_back(presel+met150+dphi);               signames.push_back("T2BW_LM150");   labels.push_back("T2BW_LM150");   uls.push_back(42.1);
      sigcuts.push_back(presel+met200+dphi);               signames.push_back("T2BW_LM200");   labels.push_back("T2BW_LM200");   uls.push_back(19.0);
      sigcuts.push_back(presel+met250+dphi);               signames.push_back("T2BW_LM250");   labels.push_back("T2BW_LM250");   uls.push_back( 9.9);

      // high-mass
      sigcuts.push_back(presel+met100+dphi+bpt100+mt2w);   signames.push_back("T2BW_HM100");   labels.push_back("T2BW_HM100");   uls.push_back(17.5);
      sigcuts.push_back(presel+met150+dphi+bpt100+mt2w);   signames.push_back("T2BW_HM150");   labels.push_back("T2BW_HM150");   uls.push_back(12.3);
      sigcuts.push_back(presel+met200+dphi+bpt100+mt2w);   signames.push_back("T2BW_HM200");   labels.push_back("T2BW_HM200");   uls.push_back( 9.0);
      sigcuts.push_back(presel+met250+dphi+bpt100+mt2w);   signames.push_back("T2BW_HM250");   labels.push_back("T2BW_HM250");   uls.push_back( 6.4);

    }

  }

  const unsigned int nsig = sigcuts.size();

  //--------------------------------------------------
  // make efficiency and xsec TH2's
  //--------------------------------------------------
  
  TH2F* hyield[nsig];
  TH2F* hsignificance[nsig];
  TH2F* hdisc[nsig];
  
  TCanvas *ctemp = new TCanvas();
  ctemp->cd();

  for( unsigned int i = 0 ; i < nsig ; ++i ){

    cout << endl << endl << endl;
    cout << "Signal region       : " << labels.at(i)   << endl << endl;
    cout << "Selection           : " << sigcuts.at(i)  << endl << endl;

    *logfile << "Signal region       : " << labels.at(i)   << endl << endl;
    *logfile << "Selection           : " << sigcuts.at(i)  << endl << endl;

    int   nbinsx  =      41;
    float xmin    =   -12.5;
    float xmax    =  1012.5;
    int   nbinsy  =      41;
    float ymin    =   -12.5;
    float ymax    =  1012.5;

    hyield[i]      = new TH2F(Form("hyield_%i",i)        , Form("hyield_%i",i)       , nbinsx , xmin , xmax , nbinsy , ymin , ymax );
    hsignificance[i]        = new TH2F(Form("hsignificance_%i",i)          , Form("hsignificance_%i",i)         , nbinsx , xmin , xmax , nbinsy , ymin , ymax );
    hdisc[i]       = new TH2F(Form("hdisc_%i",i)         , Form("hdisc_%i",i)        , nbinsx , xmin , xmax , nbinsy , ymin , ymax );

    ch->Draw(Form("ml:mg>>hyield_%i",i),sigcuts.at(i)*weight*isrweight);

    for( int ibin = 1 ; ibin <= nbinsx ; ibin++ ){
      for( int jbin = 1 ; jbin <= nbinsy ; jbin++ ){

	float yield          = hyield[i]->GetBinContent(ibin,jbin);

	hsignificance[i]->SetBinContent(ibin,jbin,0);
	hdisc[i]->SetBinContent(ibin,jbin,0);

	if( yield   < 1e-20 ) continue;

	//float my_zbi = zbi( yield + nbkg.at(i) , nbkg.at(i) , sbkg.at(i) );
	float my_zbi = simple_significance( yield , nbkg.at(i) , sbkg.at(i) );

	hsignificance[i]->SetBinContent(ibin,jbin,my_zbi);

	if( my_zbi > 5.0 ) hdisc[i]->SetBinContent(ibin,jbin,1);	
      }
    }
  }

  delete ctemp;

  cout << endl << endl;

  //--------------------------------------------------
  // make pretty pictures
  //--------------------------------------------------
  
  TLatex *t = new TLatex();
  t->SetNDC();
  t->SetTextSize(0.04);

  TCanvas* can[nsig];
  
  for( unsigned int i = 0 ; i < nsig ; ++i ){

    can[i] = new TCanvas(Form("can_%i",i),Form("can_%i",i),1200,600);
    can[i]->Divide(2,1);

    //-------------------------------
    // Yield
    //-------------------------------
  
    can[i]->cd(1);
    gPad->SetTopMargin(0.1);
    gPad->SetRightMargin(0.2);
    hyield[i]->GetXaxis()->SetLabelSize(0.035);
    hyield[i]->GetYaxis()->SetLabelSize(0.035);
    hyield[i]->GetYaxis()->SetTitle("m_{#tilde{#chi}^{0}_{1}}  [GeV]");
    hyield[i]->GetXaxis()->SetTitle("m_{#tilde{t}}  [GeV]");
    hyield[i]->GetZaxis()->SetTitle("Signal Yield");
    hyield[i]->GetZaxis()->SetTitleOffset(1.2);
    hyield[i]->GetYaxis()->SetTitleOffset(1.1);
    hyield[i]->GetXaxis()->SetRangeUser(100,800);
    hyield[i]->GetYaxis()->SetRangeUser(0,600);
    hyield[i]->Draw("colz");
  
    t->DrawLatex(0.2,0.83,label);
    t->DrawLatex(0.2,0.78,signames.at(i).c_str());
    t->DrawLatex(0.2,0.68,Form("N_{bkg} = %.1f #pm %.1f",nbkg.at(i),sbkg.at(i)));
    t->DrawLatex(0.15,0.92,"CMS Preliminary  #sqrt{s} = 8 TeV, #scale[0.6]{#int}Ldt = 19.5 fb^{-1}");
  
    //-------------------------------
    // cross section
    //-------------------------------
    
    can[i]->cd(2);
    gPad->SetTopMargin(0.1);
    gPad->SetRightMargin(0.2);
  
    hsignificance[i]->GetXaxis()->SetLabelSize(0.035);
    hsignificance[i]->GetYaxis()->SetLabelSize(0.035);
    hsignificance[i]->GetYaxis()->SetTitle("m_{#tilde{#chi}^{0}_{1}}  [GeV]");
    hsignificance[i]->GetXaxis()->SetTitle("m_{#tilde{t}}  [GeV]");
    hsignificance[i]->GetZaxis()->SetTitle("Z_{bi} Significance");
    hsignificance[i]->GetZaxis()->SetTitleOffset(1.2);
    hsignificance[i]->GetYaxis()->SetTitleOffset(1.1);
    hsignificance[i]->Draw("colz");
    hsignificance[i]->SetMinimum(0.01);
    hsignificance[i]->SetMaximum(5);
    hsignificance[i]->GetXaxis()->SetRangeUser(100,800);
    hsignificance[i]->GetYaxis()->SetRangeUser(0,600);
    hdisc[i]->Draw("samebox");

    t->DrawLatex(0.2,0.83,label);
    t->DrawLatex(0.2,0.78,signames.at(i).c_str());
    t->DrawLatex(0.15,0.92,"CMS Preliminary  #sqrt{s} = 8 TeV, #scale[0.6]{#int}Ldt = 19.5 fb^{-1}");
  
    if( print ){
      can[i]->          Print(Form("plots/%s%s%s_%s%s%s_projection.pdf"       ,sample,xchar,suffix,labels.at(i).c_str(),pol,BDTchar));
    }
  
    int bin250 = hyield[i]->FindBin(250,50);
    int bin650 = hyield[i]->FindBin(650,50);

    cout << endl;

    cout << "Signal region           : " << labels.at(i)   << endl << endl;
    cout << "Background              : " << nbkg.at(i) << " +/- " << sbkg.at(i) << endl;

    cout << "Signal yield (250/50)   : " << hyield[i]->GetBinContent(bin250)  << endl;
    cout << "Significance (250/50)   : " << hsignificance[i]->GetBinContent(bin250)    << endl;

    cout << "Signal yield (650/50)   : " << hyield[i]->GetBinContent(bin650)  << endl;
    cout << "Significance (650/50)   : " << hsignificance[i]->GetBinContent(bin650)    << endl;

    cout << endl << endl;
  }

  
  TFile *outfile = TFile::Open(Form("rootfiles/%s%s%s%s%s_projection_histos.root",sample,xchar,suffix,pol,BDTchar),"RECREATE");

  outfile->cd();
  for( unsigned int i = 0 ; i < nsig ; ++i ){
    hyield[i]->Write();
    hsignificance[i]->Write();
  }
  outfile->Close();

}


void doAll(){

  SMS_projection("T2tt", 1,true,""     ,true);
  SMS_projection("T2tt", 1,false,""     ,true);

}
