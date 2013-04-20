#include "Utils/SMS_utils.C"
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

#include "limits_T2tt_LM150.C"
#include "limits_T2tt_LM200.C"
#include "limits_T2tt_LM250.C"
#include "limits_T2tt_LM300.C"
#include "limits_T2tt_HM150.C"
#include "limits_T2tt_HM200.C"
#include "limits_T2tt_HM250.C"
#include "limits_T2tt_HM300.C"

#include "limits_T2TT_BDT1L.C"
#include "limits_T2TT_BDT1T.C"
#include "limits_T2TT_BDT2.C"
#include "limits_T2TT_BDT3.C"
#include "limits_T2TT_BDT4.C"
#include "limits_T2TT_BDT5.C"

using namespace std;

void SMS(char* sample = "T2tt" , int x = 1, bool doBDT = false , char* pol = "" , bool print = false){

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

    label     = (char*)"pp #rightarrow #tilde{t}#tilde{t}, #tilde{t} #rightarrow t+#tilde{#chi}_{1}^{0}";
  }

  else if( TString(sample).Contains("T2bw_MG") ){

    if( x==25 ){
      xchar          = (char*) "_x25";
      denomhistoname = (char*) "masses25";
      filename       = (char*) "/tas/cms2/stop/cms2V05-03-26_stoplooperV00-02-25/T2bw_mad/minibaby_V00-03-08/Skim_4jets_MET100_MT120/merged*x025*root";
      denomname      = (char*) "/tas/cms2/stop/cms2V05-03-26_stoplooperV00-02-25/T2bw_mad/minibaby_V00-03-08/Skim_4jets_MET100_MT120/myMassDB_mStop_x25_25GeVbins.root";
      label          = (char*) "pp #rightarrow #tilde{t}#tilde{t}, #tilde{t} #rightarrow b+#tilde{#chi}_{1}^{#pm}, x=0.25";
    }

    else if( x==50 ){
      xchar          = (char*) "_x50";
      denomhistoname = (char*) "masses";
      filename       = (char*) "/tas/cms2/stop/cms2V05-03-26_stoplooperV00-02-25/T2bw_mad/minibaby_V00-03-08/Skim_4jets_MET100_MT120/merged*x050*root";
      denomname      = (char*) "/tas/cms2/stop/cms2V05-03-26_stoplooperV00-02-25/T2bw_mad/minibaby_V00-03-08/Skim_4jets_MET100_MT120/myMassDB_mStop_x50_25GeVbins.root";
      label          = (char*) "pp #rightarrow #tilde{t}#tilde{t}, #tilde{t} #rightarrow b+#tilde{#chi}_{1}^{#pm}, x=0.50";
    }

    else if( x==75 ){
      xchar          = (char*) "_x75";
      denomhistoname = (char*) "masses75";
      filename       = (char*) "/tas/cms2/stop/cms2V05-03-26_stoplooperV00-02-25/T2bw_mad/minibaby_V00-03-08/Skim_4jets_MET100_MT120/merged*x075*root";
      denomname      = (char*) "/tas/cms2/stop/cms2V05-03-26_stoplooperV00-02-25/T2bw_mad/minibaby_V00-03-08/Skim_4jets_MET100_MT120/myMassDB_mStop_x75_25GeVbins.root";
      label          = (char*) "pp #rightarrow #tilde{t}#tilde{t}, #tilde{t} #rightarrow b+#tilde{#chi}_{1}^{#pm}, x=0.75";
    }
    
    else{
      cout << "ERROR! unrecognized x value " << x << ", quitting!!!" << endl;
    }
  }

  else if( TString(sample).Contains("T2bw") ){

    if( x==25 ){
      xchar          = (char*) "_x25";
      denomhistoname = (char*) "masses25";
      filename       = (char*) "/tas/cms2/stop/cms2V05-03-25_stoplooperV00-02-18/T2bw/minibabyV00-03-03/Skim_4jets_MET100_MT120/T2bw_x25.root";
      denomname      = (char*) "/tas/cms2/stop/cms2V05-03-25_stoplooperV00-02-18/T2bw/minibabyV00-03-03/myMassDB_T2bw_25GeVbins.root";
      label          = (char*) "pp #rightarrow #tilde{t}#tilde{t}, #tilde{t} #rightarrow b+#tilde{#chi}_{1}^{#pm}, x=0.25";
    }

    else if( x==50 ){
      xchar          = (char*) "_x50";
      denomhistoname = (char*) "masses";
      filename       = (char*) "/tas/cms2/stop/cms2V05-03-25_stoplooperV00-02-18/T2bw/minibabyV00-03-03/Skim_4jets_MET100_MT120/T2bw_x50.root";
      denomname      = (char*) "/tas/cms2/stop/cms2V05-03-25_stoplooperV00-02-18/T2bw/minibabyV00-03-03/myMassDB_T2bw_25GeVbins.root";
      label          = (char*) "pp #rightarrow #tilde{t}#tilde{t}, #tilde{t} #rightarrow b+#tilde{#chi}_{1}^{#pm}, x=0.50";
    }

    else if( x==75 ){
      xchar          = (char*) "_x75";
      denomhistoname = (char*) "masses75";
      filename       = (char*) "/tas/cms2/stop/cms2V05-03-25_stoplooperV00-02-18/T2bw/minibabyV00-03-03/Skim_4jets_MET100_MT120/T2bw_x75.root";
      denomname      = (char*) "/tas/cms2/stop/cms2V05-03-25_stoplooperV00-02-18/T2bw/minibabyV00-03-03/myMassDB_T2bw_25GeVbins.root";
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

  char* logfilename = Form("%s%s%s%s%s.log",sample,xchar,suffix,pol,BDTchar);
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

  if( TString(sample).Contains("T2bw") && x==25 ) presel += x25;
  if( TString(sample).Contains("T2bw") && x==50 ) presel += x50;
  if( TString(sample).Contains("T2bw") && x==75 ) presel += x75;

  int BDTweight = 1;
  if( doBDT ) BDTweight = 2;

  TCut weight(Form("mini_sltrigeff * %i",BDTweight)); 
  if( TString(pol).Contains("left") )  weight = TCut(Form("mini_sltrigeff * weightleft  * %i",BDTweight));
  if( TString(pol).Contains("right") ) weight = TCut(Form("mini_sltrigeff * weightright * %i",BDTweight));

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

  if( doBDT ){

    //-----------------------------
    // T2tt BDT
    //-----------------------------

    if( TString(sample).Contains("T2tt") ){

      cout << "Doing T2tt BDT signal regions" << endl;  

      sigcuts.push_back(TCut(presel+"mini_bdt[1] > 0.30"));  signames.push_back("T2TT_BDT1L");  labels.push_back("T2TT_BDT1L");  uls.push_back(-1.0);
      sigcuts.push_back(TCut(presel+"mini_bdt[1] > 0.40"));  signames.push_back("T2TT_BDT1T");  labels.push_back("T2TT_BDT1T");  uls.push_back(-1.0);
      sigcuts.push_back(TCut(presel+"mini_bdt[2] > 0.55"));  signames.push_back("T2TT_BDT2");   labels.push_back("T2TT_BDT2");   uls.push_back(-1.0);
      sigcuts.push_back(TCut(presel+"mini_bdt[3] > 0.65"));  signames.push_back("T2TT_BDT3");   labels.push_back("T2TT_BDT3");   uls.push_back(-1.0);
      sigcuts.push_back(TCut(presel+"mini_bdt[4] > 0.50"));  signames.push_back("T2TT_BDT4");   labels.push_back("T2TT_BDT4");   uls.push_back(-1.0);
      sigcuts.push_back(TCut(presel+"mini_bdt[5] > 0.30"));  signames.push_back("T2TT_BDT5");   labels.push_back("T2TT_BDT5");   uls.push_back(-1.0);
    }

    //-----------------------------
    // T2bw BDT
    //-----------------------------

    else if( TString(sample).Contains("T2bw") ){

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
  
  TH2F* heff[nsig];
  TH2F* hnevents[nsig];
  TH2F* heff_noisr[nsig];
  TH2F* heffup[nsig];
  TH2F* heffdn[nsig];
  TH2F* heffbup[nsig];
  TH2F* heffbdn[nsig];
  TH2F* hxsec[nsig];
  TH2F* hxsec_exp[nsig];
  TH2F* hxsec_expp1[nsig];
  TH2F* hxsec_expm1[nsig];
  TH2F* hexcl[nsig];
  TH2F* hexcl_exp[nsig];
  TH2F* hexcl_expp1[nsig];
  TH2F* hexcl_expm1[nsig];
  TH2F* hexcl_obsp1[nsig];
  TH2F* hexcl_obsm1[nsig];
  TH2F* hjes[nsig];
  TH2F* hbtagerr[nsig];
  TH2F* hstaterr[nsig];
  TH2F* htoterr[nsig];
  TH2F* hisrerr[nsig];
  
  TCanvas *ctemp = new TCanvas();
  ctemp->cd();

  for( unsigned int i = 0 ; i < nsig ; ++i ){

    //------------------------------------------------
    // calculate point-by-point signal uncertainties
    //------------------------------------------------

    TString jesup(sigcuts.at(i));
    jesup.ReplaceAll("mini_njets"         , "mini_njetsup"    );
    jesup.ReplaceAll("mini_met"           , "mini_metup"      );
    jesup.ReplaceAll(" mini_mt"           , "mini_mtup"       );
    jesup.ReplaceAll("mini_chi2"          , "mini_chi2up"     );
    jesup.ReplaceAll("mini_mt2w"          , "mini_mt2wup"     );
    jesup.ReplaceAll("mini_pt_b"          , "mini_pt_b_up"     );
    jesup.ReplaceAll("mini_bdt"           , "mini_bdtup"      );

    TString jesdown(sigcuts.at(i));
    jesdown.ReplaceAll("mini_njets"       , "mini_njetsdown"  );
    jesdown.ReplaceAll("mini_met"         , "mini_metdown"    );
    jesdown.ReplaceAll(" mini_mt"         , "mini_mtdown"     );
    jesdown.ReplaceAll("mini_chi2"        , "mini_chi2down"   );
    jesdown.ReplaceAll("mini_mt2w"        , "mini_mt2wdown"   );
    jesdown.ReplaceAll("mini_pt_b"        , "mini_pt_b_down"   );
    jesdown.ReplaceAll("mini_bdt"         , "mini_bdtdown"    );

    TString btagup(sigcuts.at(i));
    btagup.ReplaceAll("mini_nb"           , "mini_nbupBC"     );
    btagup.ReplaceAll("mini_chi2"         , "mini_chi2bup"    );
    btagup.ReplaceAll("mini_mt2w"         , "mini_mt2wbup"    );
    btagup.ReplaceAll("mini_pt_b"         , "mini_pt_b_bup"   );
    btagup.ReplaceAll("mini_bdt"          , "mini_bdtbup"     );

    TString btagdn(sigcuts.at(i));
    btagdn.ReplaceAll("mini_nb"         , "mini_nbdownBC"   );
    btagdn.ReplaceAll("mini_chi2"       , "mini_chi2bdown"  );
    btagdn.ReplaceAll("mini_mt2w"       , "mini_mt2wbdown"  );
    btagdn.ReplaceAll("mini_pt_b"       , "mini_pt_b_bdown" );
    btagdn.ReplaceAll("mini_bdt"        , "mini_bdtbdown"   );

    TCut jesupcut(jesup);
    TCut jesdncut(jesdown);
    TCut btagupcut(btagup);
    TCut btagdncut(btagdn);

    cout << endl << endl << endl;
    cout << "Signal region       : " << labels.at(i)   << endl << endl;
    cout << "Selection           : " << sigcuts.at(i)  << endl << endl;
    cout << "Selection JES up    : " << jesupcut       << endl << endl;
    cout << "Selection JES down  : " << jesdncut       << endl << endl;
    cout << "Selection btag up   : " << btagupcut      << endl << endl;
    cout << "Selection btag down : " << btagdncut      << endl << endl;

    *logfile << "Signal region       : " << labels.at(i)   << endl << endl;
    *logfile << "Selection           : " << sigcuts.at(i)  << endl << endl;
    *logfile << "Selection JES up    : " << jesupcut       << endl << endl;
    *logfile << "Selection JES down  : " << jesdncut       << endl << endl;
    *logfile << "Selection btag up   : " << btagupcut      << endl << endl;
    *logfile << "Selection btag down : " << btagdncut      << endl << endl;

    int   nbinsx  =      41;
    float xmin    =   -12.5;
    float xmax    =  1012.5;
    int   nbinsy  =      41;
    float ymin    =   -12.5;
    float ymax    =  1012.5;

    heff[i]        = new TH2F(Form("heff_%i",i)          , Form("heff_%i",i)         , nbinsx , xmin , xmax , nbinsy , ymin , ymax );
    hnevents[i]    = new TH2F(Form("hnevents_%i",i)      , Form("hnevents_%i",i)     , nbinsx , xmin , xmax , nbinsy , ymin , ymax );
    heff_noisr[i]  = new TH2F(Form("heff_noisr_%i",i)    , Form("heff_noisr_%i",i)   , nbinsx , xmin , xmax , nbinsy , ymin , ymax );
    heffup[i]      = new TH2F(Form("heffup_%i",i)        , Form("heffup_%i",i)       , nbinsx , xmin , xmax , nbinsy , ymin , ymax );
    heffdn[i]      = new TH2F(Form("heffdn_%i",i)        , Form("heffdn_%i",i)       , nbinsx , xmin , xmax , nbinsy , ymin , ymax );
    heffbup[i]     = new TH2F(Form("heffbup_%i",i)       , Form("heffbup_%i",i)      , nbinsx , xmin , xmax , nbinsy , ymin , ymax );
    heffbdn[i]     = new TH2F(Form("heffbdn_%i",i)       , Form("heffbdn_%i",i)      , nbinsx , xmin , xmax , nbinsy , ymin , ymax );
    hxsec[i]       = new TH2F(Form("hxsec_%i",i)         , Form("hxsec_%i",i)        , nbinsx , xmin , xmax , nbinsy , ymin , ymax );
    hxsec_exp[i]   = new TH2F(Form("hxsec_exp_%i",i)     , Form("hxsec_exp_%i",i)    , nbinsx , xmin , xmax , nbinsy , ymin , ymax );
    hxsec_expp1[i] = new TH2F(Form("hxsec_expp1_%i",i)   , Form("hxsec_expp1_%i",i)  , nbinsx , xmin , xmax , nbinsy , ymin , ymax );
    hxsec_expm1[i] = new TH2F(Form("hxsec_expm1_%i",i)   , Form("hxsec_expm1_%i",i)  , nbinsx , xmin , xmax , nbinsy , ymin , ymax );
 
    hexcl[i]       = new TH2F(Form("hexcl_%i",i)         , Form("hexcl_%i",i)        , nbinsx , xmin , xmax , nbinsy , ymin , ymax );
    hexcl_obsp1[i] = new TH2F(Form("hexcl_obsp1_%i",i)   , Form("hexcl_obsp1_%i",i)  , nbinsx , xmin , xmax , nbinsy , ymin , ymax );
    hexcl_obsm1[i] = new TH2F(Form("hexcl_obsm1_%i",i)   , Form("hexcl_obsm1_%i",i)  , nbinsx , xmin , xmax , nbinsy , ymin , ymax );
    hexcl_exp[i]   = new TH2F(Form("hexcl_exp_%i",i)     , Form("hexcl_exp_%i",i)    , nbinsx , xmin , xmax , nbinsy , ymin , ymax );
    hexcl_expp1[i] = new TH2F(Form("hexcl_expp1_%i",i)   , Form("hexcl_expp1_%i",i)  , nbinsx , xmin , xmax , nbinsy , ymin , ymax );
    hexcl_expm1[i] = new TH2F(Form("hexcl_expm1_%i",i)   , Form("hexcl_expm1_%i",i)  , nbinsx , xmin , xmax , nbinsy , ymin , ymax );

    hjes[i]        = new TH2F(Form("hjes_%i",i)          , Form("hjes_%i",i)         , nbinsx , xmin , xmax , nbinsy , ymin , ymax );
    htoterr[i]     = new TH2F(Form("htoterr_%i",i)       , Form("htoterr_%i",i)      , nbinsx , xmin , xmax , nbinsy , ymin , ymax );
    hisrerr[i]     = new TH2F(Form("hisrerr_%i",i)       , Form("hisrerr_%i",i)      , nbinsx , xmin , xmax , nbinsy , ymin , ymax );
    hbtagerr[i]    = new TH2F(Form("hbtagerr_%i",i)      , Form("hbtagerr_%i",i)     , nbinsx , xmin , xmax , nbinsy , ymin , ymax );
    hstaterr[i]    = new TH2F(Form("hstaterr_%i",i)      , Form("hstaterr_%i",i)     , nbinsx , xmin , xmax , nbinsy , ymin , ymax );

    heff[i]->Sumw2();

    // nominal
    ch->Draw(Form("ml:mg>>heff_%i",i),sigcuts.at(i)*weight*isrweight);
    TH2F* hnevents = (TH2F*) heff[i]->Clone(Form("hnevents_%i",i));
    heff[i]->Divide(hdenom);

    // raw number of signal events
    ch->Draw(Form("ml:mg>>hnevents_%i",i),sigcuts.at(i));

    // JES up
    ch->Draw(Form("ml:mg>>heffup_%i",i),jesupcut*weight*isrweight);
    heffup[i]->Divide(hdenom);

    // JES down
    ch->Draw(Form("ml:mg>>heffdn_%i",i),jesdncut*weight*isrweight);
    heffdn[i]->Divide(hdenom);

    // btag up
    ch->Draw(Form("ml:mg>>heffbup_%i",i),btagupcut*weight*isrweight);
    heffbup[i]->Divide(hdenom);

    // btag down
    ch->Draw(Form("ml:mg>>heffbdn_%i",i),btagdncut*weight*isrweight);
    heffbdn[i]->Divide(hdenom);

    // remove ISR weight
    ch->Draw(Form("ml:mg>>heff_noisr_%i",i),sigcuts.at(i)*weight);
    heff_noisr[i]->Divide(hdenom);

    for( int ibin = 1 ; ibin <= nbinsx ; ibin++ ){
      for( int jbin = 1 ; jbin <= nbinsy ; jbin++ ){

	hjes[i]->SetBinContent(ibin,jbin,0.0);
	hbtagerr[i]->SetBinContent(ibin,jbin,0.0);
	htoterr[i]->SetBinContent(ibin,jbin,0.0);
	hisrerr[i]->SetBinContent(ibin,jbin,0.0);
	hstaterr[i]->SetBinContent(ibin,jbin,0.0);

	float mg = heff[i]->GetXaxis()->GetBinCenter(ibin);
	//float ml = heff[i]->GetYaxis()->GetBinCenter(jbin);

	// nominal
	float eff          = heff[i]->GetBinContent(ibin,jbin);
	float eff_staterr  = 1.0/sqrt(hnevents->GetBinContent(ibin,jbin));

	// JES up/down
	float effup      = heffup[i]->GetBinContent(ibin,jbin);
	float effdn      = heffdn[i]->GetBinContent(ibin,jbin);

	// btag up/down
	float effbup     = heffbup[i]->GetBinContent(ibin,jbin);
	float effbdn     = heffbdn[i]->GetBinContent(ibin,jbin);

	// no ISR weighting
	float eff_noisr  = heff_noisr[i]->GetBinContent(ibin,jbin);

	if( eff   < 1e-20 ) continue;

	// JES uncertainty
	float dup    = fabs(effup/eff-1);
	float ddn    = fabs(1-effdn/eff);
	float djes   = 0.5 * (dup+ddn);

	// btag uncertainty
	float dbup    = fabs(effbup/eff-1);
	float dbdn    = fabs(1-effbdn/eff);
	float dbtag   = 0.5 * (dbup+dbdn);

	// ISR uncertainty
	float disr   = fabs(1-eff/eff_noisr);

	// lumi (4.4%), trigger (3%), lepton selection (5%), JES, ISR, btagging
	//float toterr  = sqrt( 0.044*0.044 + 0.03*0.03 + 0.05*0.05 + djes*djes + disr*disr + dbtag * dbtag + eff_staterr*eff_staterr);
	float toterr  = sqrt( 0.044*0.044 + 0.03*0.03 + 0.05*0.05 + djes*djes + disr*disr + dbtag * dbtag );
	htoterr[i]->SetBinContent(ibin,jbin,toterr);
	hisrerr[i]->SetBinContent(ibin,jbin,disr);
	hbtagerr[i]->SetBinContent(ibin,jbin,dbtag);
	hjes[i]->SetBinContent(ibin,jbin,djes);
	hstaterr[i]->SetBinContent(ibin,jbin,eff_staterr);

	float this_ul;
	float this_ul_exp;
	float this_ul_expp1;
	float this_ul_expm1;

	// this_ul       = uls.at(i);
	// this_ul_exp   = uls.at(i);
	// this_ul_expp1 = uls.at(i);
	// this_ul_expm1 = uls.at(i);
	
	//------------------------------------------
	// T2tt cut-and-count
	//------------------------------------------

	if( TString(labels.at(i)).Contains("T2TT_LM150") ){
	  // this_ul         = 66.7991;
	  // this_ul_exp     = 76.7344;
	  // this_ul_expp1   = 104.545;
	  // this_ul_expm1   = 56.8783;

	  this_ul       = getUpperLimit_T2tt_LM150( toterr );
	  this_ul_exp   = getExpectedUpperLimit_T2tt_LM150( toterr );
	  this_ul_expp1 = getExpectedP1UpperLimit_T2tt_LM150( toterr );
	  this_ul_expm1 = getExpectedM1UpperLimit_T2tt_LM150( toterr );
	}

	else if( TString(labels.at(i)).Contains("T2TT_LM200") ){
	  // this_ul         = 25.714;
	  // this_ul_exp     = 32.2226;
	  // this_ul_expp1   = 43.9916;
	  // this_ul_expm1   = 23.6334;

	  this_ul       = getUpperLimit_T2tt_LM200( toterr );
	  this_ul_exp   = getExpectedUpperLimit_T2tt_LM200( toterr );
	  this_ul_expp1 = getExpectedP1UpperLimit_T2tt_LM200( toterr );
	  this_ul_expm1 = getExpectedM1UpperLimit_T2tt_LM200( toterr );
	}

	else if( TString(labels.at(i)).Contains("T2TT_LM250") ){
	  // this_ul         = 10.3561;
	  // this_ul_exp     = 15.1845;
	  // this_ul_expp1   = 23.3334;
	  // this_ul_expm1   = 10.9627;

	  this_ul       = getUpperLimit_T2tt_LM250( toterr );
	  this_ul_exp   = getExpectedUpperLimit_T2tt_LM250( toterr );
	  this_ul_expp1 = getExpectedP1UpperLimit_T2tt_LM250( toterr );
	  this_ul_expm1 = getExpectedM1UpperLimit_T2tt_LM250( toterr );
	}

	else if( TString(labels.at(i)).Contains("T2TT_LM300") ){
	  // this_ul         = 7.7604;
	  // this_ul_exp     = 9.05201;
	  // this_ul_expp1   = 12.7193;
	  // this_ul_expm1   = 6.55692;

	  this_ul       = getUpperLimit_T2tt_LM300( toterr );
	  this_ul_exp   = getExpectedUpperLimit_T2tt_LM300( toterr );
	  this_ul_expp1 = getExpectedP1UpperLimit_T2tt_LM300( toterr );
	  this_ul_expm1 = getExpectedM1UpperLimit_T2tt_LM300( toterr );
	}

	else if( TString(labels.at(i)).Contains("T2TT_HM150") ){
	  // this_ul         = 11.5352;
	  // this_ul_exp     = 14.5708;
	  // this_ul_expp1   = 20.3566;
	  // this_ul_expm1   = 10.4048;

	  this_ul       = getUpperLimit_T2tt_HM150( toterr );
	  this_ul_exp   = getExpectedUpperLimit_T2tt_HM150( toterr );
	  this_ul_expp1 = getExpectedP1UpperLimit_T2tt_HM150( toterr );
	  this_ul_expm1 = getExpectedM1UpperLimit_T2tt_HM150( toterr );
	}

	else if( TString(labels.at(i)).Contains("T2TT_HM200") ){
	  // this_ul         = 7.4327;
	  // this_ul_exp     = 10.2928;
	  // this_ul_expp1   = 14.6985;
	  // this_ul_expm1   = 7.41028;

	  this_ul       = getUpperLimit_T2tt_HM200( toterr );
	  this_ul_exp   = getExpectedUpperLimit_T2tt_HM200( toterr );
	  this_ul_expp1 = getExpectedP1UpperLimit_T2tt_HM200( toterr );
	  this_ul_expm1 = getExpectedM1UpperLimit_T2tt_HM200( toterr );
	}

	else if( TString(labels.at(i)).Contains("T2TT_HM250") ){
	  // this_ul         = 4.15662;
	  // this_ul_exp     = 7.37302;
	  // this_ul_expp1   = 10.4398;
	  // this_ul_expm1   = 5.09026;

	  this_ul       = getUpperLimit_T2tt_HM250( toterr );
	  this_ul_exp   = getExpectedUpperLimit_T2tt_HM250( toterr );
	  this_ul_expp1 = getExpectedP1UpperLimit_T2tt_HM250( toterr );
	  this_ul_expm1 = getExpectedM1UpperLimit_T2tt_HM250( toterr );
	}

	else if( TString(labels.at(i)).Contains("T2TT_HM300") ){
	  // this_ul         = 3.97749;
	  // this_ul_exp     = 5.52676;
	  // this_ul_expp1   = 7.95513;
	  // this_ul_expm1   = 4.05663;

	  this_ul       = getUpperLimit_T2tt_HM300( toterr );
	  this_ul_exp   = getExpectedUpperLimit_T2tt_HM300( toterr );
	  this_ul_expp1 = getExpectedP1UpperLimit_T2tt_HM300( toterr );
	  this_ul_expm1 = getExpectedM1UpperLimit_T2tt_HM300( toterr );
	}

	//------------------------------------------
	// T2bw cut-and-count
	//------------------------------------------

	else if( TString(labels.at(i)).Contains("T2BW_LM100") ){
	  this_ul         = 321.149;
	  this_ul_exp     = 340.097;
	  this_ul_expp1   = 476.991;
	  this_ul_expm1   = 248.095;
	}

	else if( TString(labels.at(i)).Contains("T2BW_LM150") ){
	  this_ul         = 97.8571;
	  this_ul_exp     = 121.299;
	  this_ul_expp1   = 166.086;
	  this_ul_expm1   = 89.0909;
	}

	else if( TString(labels.at(i)).Contains("T2BW_LM200") ){
	  this_ul         = 35.556;
	  this_ul_exp     = 48.3854;
	  this_ul_expp1   = 66.5798;
	  this_ul_expm1   = 35.6169;
	}

	else if( TString(labels.at(i)).Contains("T2BW_LM250") ){
	  this_ul         = 17.3482;
	  this_ul_exp     = 24.0754;
	  this_ul_expp1   = 33.372;
	  this_ul_expm1   = 17.1881;
	}

	else if( TString(labels.at(i)).Contains("T2BW_HM100") ){
	  this_ul         = 36.9713;
	  this_ul_exp     = 29.0217;
	  this_ul_expp1   = 40.1313;
	  this_ul_expm1   = 21.6054;
	}

	else if( TString(labels.at(i)).Contains("T2BW_HM150") ){
	  this_ul         = 19.3805;
	  this_ul_exp     = 18.3767;
	  this_ul_expp1   = 25.448;
	  this_ul_expm1   = 13.2759;
	}

	else if( TString(labels.at(i)).Contains("T2BW_HM200") ){
	  this_ul         = 11.4308;
	  this_ul_exp     = 11.9613;
	  this_ul_expp1   = 17.4688;
	  this_ul_expm1   = 8.99217;
	}

	else if( TString(labels.at(i)).Contains("T2BW_HM250") ){
	  this_ul         = 5.04346;
	  this_ul_exp     = 7.73816;
	  this_ul_expp1   = 11.1644;
	  this_ul_expm1   = 5.39447;
	}

	//------------------------------------------
	// T2tt BDT
	//------------------------------------------

	else if( TString(labels.at(i)).Contains("T2TT_BDT1L") ){
	  // this_ul         = 150.058;
	  // this_ul_exp     = 165.876;
	  // this_ul_expp1   = 228.763;
	  // this_ul_expm1   = 121.923;

	  this_ul       = getUpperLimit_T2TT_BDT1L( toterr );
	  this_ul_exp   = getExpectedUpperLimit_T2TT_BDT1L( toterr );
	  this_ul_expp1 = getExpectedP1UpperLimit_T2TT_BDT1L( toterr );
	  this_ul_expm1 = getExpectedM1UpperLimit_T2TT_BDT1L( toterr );
	}

	else if( TString(labels.at(i)).Contains("T2TT_BDT1T") ){
	  // this_ul         = 27.8416;
	  // this_ul_exp     = 36.9395;
	  // this_ul_expp1   = 50.9219;
	  // this_ul_expm1   = 27.1247;

	  this_ul       = getUpperLimit_T2TT_BDT1T( toterr );
	  this_ul_exp   = getExpectedUpperLimit_T2TT_BDT1T( toterr );
	  this_ul_expp1 = getExpectedP1UpperLimit_T2TT_BDT1T( toterr );
	  this_ul_expm1 = getExpectedM1UpperLimit_T2TT_BDT1T( toterr );
	}

	else if( TString(labels.at(i)).Contains("T2TT_BDT2") ){
	  // this_ul         = 15.5022;
	  // this_ul_exp     = 26.193;
	  // this_ul_expp1   = 35.9907;
	  // this_ul_expm1   = 18.5872;

	  this_ul       = getUpperLimit_T2TT_BDT2( toterr );
	  this_ul_exp   = getExpectedUpperLimit_T2TT_BDT2( toterr );
	  this_ul_expp1 = getExpectedP1UpperLimit_T2TT_BDT2( toterr );
	  this_ul_expm1 = getExpectedM1UpperLimit_T2TT_BDT2( toterr );
	}

	else if( TString(labels.at(i)).Contains("T2TT_BDT3") ){
	  // this_ul         = 6.42278;
	  // this_ul_exp     = 9.12343;
	  // this_ul_expp1   = 12.9316;
	  // this_ul_expm1   = 6.50077;

	  this_ul       = getUpperLimit_T2TT_BDT3( toterr );
	  this_ul_exp   = getExpectedUpperLimit_T2TT_BDT3( toterr );
	  this_ul_expp1 = getExpectedP1UpperLimit_T2TT_BDT3( toterr );
	  this_ul_expm1 = getExpectedM1UpperLimit_T2TT_BDT3( toterr );
	}

	else if( TString(labels.at(i)).Contains("T2TT_BDT4") ){
	  // this_ul         = 4.36022;
	  // this_ul_exp     = 4.87426;
	  // this_ul_expp1   = 7.392;
	  // this_ul_expm1   = 3.68327;

	  this_ul       = getUpperLimit_T2TT_BDT4( toterr );
	  this_ul_exp   = getExpectedUpperLimit_T2TT_BDT4( toterr );
	  this_ul_expp1 = getExpectedP1UpperLimit_T2TT_BDT4( toterr );
	  this_ul_expm1 = getExpectedM1UpperLimit_T2TT_BDT4( toterr );
	}

	else if( TString(labels.at(i)).Contains("T2TT_BDT5") ){
	  // this_ul         = 26.1597;
	  // this_ul_exp     = 31.2588;
	  // this_ul_expp1   = 42.6817;
	  // this_ul_expm1   = 22.7823;

	  this_ul       = getUpperLimit_T2TT_BDT5( toterr );
	  this_ul_exp   = getExpectedUpperLimit_T2TT_BDT5( toterr );
	  this_ul_expp1 = getExpectedP1UpperLimit_T2TT_BDT5( toterr );
	  this_ul_expm1 = getExpectedM1UpperLimit_T2TT_BDT5( toterr );
	}

	//------------------------------------------
	// T2bw BDT
	//------------------------------------------

	else if( TString(labels.at(i)).Contains("T2BW25_BDT2") ){
	  this_ul         = 6.36826;
	  this_ul_exp     = 8.05362;
	  this_ul_expp1   = 11.4484;
	  this_ul_expm1   = 5.6441;
	}

	else if( TString(labels.at(i)).Contains("T2BW25_BDT3") ){
	  this_ul         = 3.99136;
	  this_ul_exp     = 5.38962;
	  this_ul_expp1   = 8.48055;
	  this_ul_expm1   = 3.95123;
	}

	else if( TString(labels.at(i)).Contains("T2BW50_BDT1") ){
	  this_ul         = 20.6657;
	  this_ul_exp     = 26.0337;
	  this_ul_expp1   = 35.8818;
	  this_ul_expm1   = 18.7248;
	}

	else if( TString(labels.at(i)).Contains("T2BW50_BDT2L") ){
	  this_ul         = 11.6687;
	  this_ul_exp     = 18.6775;
	  this_ul_expp1   = 26.1434;
	  this_ul_expm1   = 13.0451;
	}

	else if( TString(labels.at(i)).Contains("T2BW50_BDT2T") ){
	  this_ul         = 9.68818;
	  this_ul_exp     = 9.84362;
	  this_ul_expp1   = 13.8548;
	  this_ul_expm1   = 7.31892;
	}

	else if( TString(labels.at(i)).Contains("T2BW50_BDT3") ){
	  this_ul         = 8.27178;
	  this_ul_exp     = 10.6125;
	  this_ul_expp1   = 16.1522;
	  this_ul_expm1   = 7.67213;
	}

	else if( TString(labels.at(i)).Contains("T2BW75_BDT1") ){
	  this_ul         = 16.1085;
	  this_ul_exp     = 22.0458;
	  this_ul_expp1   = 30.6627;
	  this_ul_expm1   = 15.5538;
	}

	else if( TString(labels.at(i)).Contains("T2BW75_BDT2") ){
	  this_ul         = 7.14169;
	  this_ul_exp     = 11.6002;
	  this_ul_expp1   = 16.4583;
	  this_ul_expm1   = 8.26508;
	}

	else if( TString(labels.at(i)).Contains("T2BW75_BDT3") ){
	  this_ul         = 5.45205;
	  this_ul_exp     = 7.23529;
	  this_ul_expp1   = 10.3813;
	  this_ul_expm1   = 5.19484;
	}

	else if( TString(labels.at(i)).Contains("T2BW75_BDT4") ){
	  this_ul         = 103.755;
	  this_ul_exp     = 121.37;
	  this_ul_expp1   = 164.736;
	  this_ul_expm1   = 88.7762;
	}
	else{
	  cout << "ERROR! UNRECOGNIZED SIGNAL REGION " << labels.at(i) << endl;
	  exit(0);
	}

	// cout << endl;
	// cout << "toterr " << toterr << endl;
	// cout << "this_ul        " << this_ul       << endl;
	// cout << "this_ul_exp    " << this_ul_exp   << endl;
	// cout << "this_ul_expp1  " << this_ul_expp1 << endl;
	// cout << "this_ul_expm1  " << this_ul_expm1 << endl;

	float xsecul        = this_ul       / ( lumi * eff );
	float xsecul_exp    = this_ul_exp   / ( lumi * eff );
	float xsecul_expp1  = this_ul_expp1 / ( lumi * eff );
	float xsecul_expm1  = this_ul_expm1 / ( lumi * eff );

	if( eff > 0 ){
	  hxsec[i]->SetBinContent(ibin,jbin, xsecul );
	  hxsec_exp[i]->SetBinContent(ibin,jbin, xsecul_exp );
	  hxsec_expp1[i]->SetBinContent(ibin,jbin, xsecul_expp1 );
	  hxsec_expm1[i]->SetBinContent(ibin,jbin, xsecul_expm1 );
	}

	int   bin     = refxsec->FindBin(mg);
	float xsec    = refxsec->GetBinContent(bin);
	float xsec_up = refxsec->GetBinContent(bin) + refxsec->GetBinError(bin);
	float xsec_dn = refxsec->GetBinContent(bin) - refxsec->GetBinError(bin);

	hexcl[i]->SetBinContent(ibin,jbin,0);
	if( xsec > xsecul )   hexcl[i]->SetBinContent(ibin,jbin,1);

	hexcl_exp[i]->SetBinContent(ibin,jbin,0);
	if( xsec > xsecul_exp )   hexcl_exp[i]->SetBinContent(ibin,jbin,1);

	hexcl_expp1[i]->SetBinContent(ibin,jbin,0);
	if( xsec > xsecul_expp1 )   hexcl_expp1[i]->SetBinContent(ibin,jbin,1);

	hexcl_expm1[i]->SetBinContent(ibin,jbin,0);
	if( xsec > xsecul_expm1 )   hexcl_expm1[i]->SetBinContent(ibin,jbin,1);

	hexcl_obsp1[i]->SetBinContent(ibin,jbin,0);
	if( xsec_up > xsecul )   hexcl_obsp1[i]->SetBinContent(ibin,jbin,1);

	hexcl_obsm1[i]->SetBinContent(ibin,jbin,0);
	if( xsec_dn > xsecul )   hexcl_obsm1[i]->SetBinContent(ibin,jbin,1);

	//cout << "ibin jbin mg xsec " << ibin << " " << jbin << " " << mg << " " << xsec << endl;
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
  TCanvas* can_exclusion[nsig];

  TGraph* gr[nsig];
  TGraph* gr_exp[nsig];
  TGraph* gr_expp1[nsig];
  TGraph* gr_expm1[nsig];


  for( unsigned int i = 0 ; i < nsig ; ++i ){

    //TGraph *gr     = getGraph( sample , "observed" , signames.at(i) );
    //TGraph *gr_exp = getGraph( sample , "expected" , signames.at(i) );

    gr[i]      = getRefXsecGraph(hxsec[i]     , (char*) "T2tt", 1.0);
    gr_exp[i]  = getRefXsecGraph(hxsec_exp[i] , (char*) "T2tt", 1.0);

    gr[i]->SetLineWidth(3);
    gr_exp[i]->SetLineWidth(3);
    gr_exp[i]->SetLineStyle(2);

    gr[i]->SetName(Form("gr_%i",i));
    gr_exp[i]->SetName(Form("gr_exp_%i",i));
    gr[i]->SetTitle(Form("gr_%i",i));
    gr_exp[i]->SetTitle(Form("gr_exp_%i",i));

    //can[i] = new TCanvas(Form("can_%i",i),Form("can_%i",i),1200,600);
    //can[i]->Divide(2,1);
    //can[i] = new TCanvas(Form("can_%i",i),Form("can_%i",i),1800,600);
    //can[i]->Divide(3,1);

    can[i] = new TCanvas(Form("can_%i",i),Form("can_%i",i),1200,600);
    can[i]->Divide(2,1);

    //-------------------------------
    // efficiency
    //-------------------------------

    can[i]->cd(1);
    gPad->SetTopMargin(0.1);
    gPad->SetRightMargin(0.2);
    heff[i]->Scale(100);
    heff[i]->GetXaxis()->SetLabelSize(0.035);
    heff[i]->GetYaxis()->SetLabelSize(0.035);
    heff[i]->GetYaxis()->SetTitle("#chi^{0}_{1} mass (GeV)");
    heff[i]->GetXaxis()->SetTitle("stop mass (GeV)");
    heff[i]->GetZaxis()->SetTitle("efficiency (%)");
    heff[i]->GetZaxis()->SetTitleOffset(1.2);
    heff[i]->GetXaxis()->SetRangeUser(100,800);
    heff[i]->GetYaxis()->SetRangeUser(0,600);
    heff[i]->Draw("colz");
    //heff[i]->Draw("sametext");

    t->DrawLatex(0.2,0.83,label);
    //t->DrawLatex(0.2,0.77,"m(#tilde{q}) >> m(#tilde{g})");
    t->DrawLatex(0.2,0.78,signames.at(i).c_str());
    t->DrawLatex(0.15,0.92,"CMS Preliminary  #sqrt{s} = 8 TeV, #scale[0.6]{#int}Ldt = 19.5 fb^{-1}");

    //-------------------------------
    // cross section
    //-------------------------------
  
    can[i]->cd(2);
    gPad->SetTopMargin(0.1);
    gPad->SetRightMargin(0.2);
    gPad->SetLogz();

    hxsec[i]->GetXaxis()->SetLabelSize(0.035);
    hxsec[i]->GetYaxis()->SetLabelSize(0.035);
    hxsec[i]->GetYaxis()->SetTitle("#chi^{0}_{1} mass (GeV)");
    hxsec[i]->GetXaxis()->SetTitle("stop mass (GeV)");
    hxsec[i]->GetZaxis()->SetTitle("#sigma upper limit [pb]");
    hxsec[i]->GetZaxis()->SetTitleOffset(1.2);
    hxsec[i]->Draw("colz");
    //hxsec[i]->Draw("sametext");
    hxsec[i]->SetMinimum(0.01);
    hxsec[i]->SetMaximum(100);
    hxsec[i]->GetXaxis()->SetRangeUser(100,800);
    hxsec[i]->GetYaxis()->SetRangeUser(0,600);
    hexcl[i]->Draw("samebox");

    gr[i]->Draw("same");
    gr_exp[i]->Draw("same");

    TLegend *leg = new TLegend(0.2,0.6,0.4,0.75);
    leg->AddEntry(gr[i],    "observed" ,"l");
    leg->AddEntry(gr_exp[i],"expected" ,"l");
    leg->SetFillColor(0);
    leg->SetBorderSize(0);
    leg->Draw();


    // TGraph* gr_excl      = getRefXsecGraph(hxsec[i], "T5zz", 1.0);
    // TGraph* gr_excl_down = getRefXsecGraph(hxsec[i], "T5zz", 1./3.);
    // TGraph* gr_excl_up   = getRefXsecGraph(hxsec[i], "T5zz", 3.);

    // gr_excl->SetLineWidth(2);
    // gr_excl_up->SetLineWidth(2);
    // gr_excl_down->SetLineWidth(2);
    // gr_excl_up->SetLineStyle(2);
    // gr_excl_down->SetLineStyle(3);
    // gr_excl->Draw("same");
    // gr_excl_up->Draw("same");
    // gr_excl_down->Draw("same");


    t->DrawLatex(0.2,0.83,label);
    //t->DrawLatex(0.2,0.77,"m(#tilde{q}) >> m(#tilde{g})");
    t->DrawLatex(0.2,0.78,signames.at(i).c_str());
    t->DrawLatex(0.15,0.92,"CMS Preliminary  #sqrt{s} = 8 TeV, #scale[0.6]{#int}Ldt = 19.5 fb^{-1}");

    //-------------------------------
    // excluded points
    //-------------------------------
    /*
    can_exclusion[i] = new TCanvas(Form("can_exclusion_%i",i),Form("can_exclusion_%i",i),1200,600);
    can_exclusion[i]->Divide(2,1);

    can_exclusion[i]->cd(1);    
    gPad->SetRightMargin(0.2);
    gPad->SetTopMargin(0.1);
    hexcl[i]->SetMinimum(0);
    hexcl[i]->SetMaximum(1);
    hexcl[i]->GetYaxis()->SetTitle("#chi^{0}_{1} mass (GeV)");
    hexcl[i]->GetXaxis()->SetTitle("#tilde{t} mass (GeV)");
    hexcl[i]->GetZaxis()->SetTitle("observed excluded points");
    hexcl[i]->Draw("colz");
    gr[i]->Draw("l");

    t->DrawLatex(0.2,0.83,label);
    //t->DrawLatex(0.2,0.77,"m(#tilde{q}) >> m(#tilde{g})");
    t->DrawLatex(0.2,0.71,signames.at(i).c_str());
    t->DrawLatex(0.15,0.92,"CMS Preliminary  #sqrt{s} = 8 TeV, #scale[0.6]{#int}Ldt = 19.5 fb^{-1}");

    //-------------------------------
    // JES uncertainty
    //-------------------------------

    can_exclusion[i]->cd(2);
    gPad->SetRightMargin(0.2);
    gPad->SetTopMargin(0.1);
    hexcl_exp[i]->SetMinimum(0);
    hexcl_exp[i]->SetMaximum(1);
    hexcl_exp[i]->GetYaxis()->SetTitle("#chi^{0}_{1} mass (GeV)");
    hexcl_exp[i]->GetXaxis()->SetTitle("#tilde{t} mass (GeV)");
    hexcl_exp[i]->GetZaxis()->SetTitle("expected excluded points");
    hexcl_exp[i]->Draw("colz");
    gr_exp[i]->Draw("l");

    t->DrawLatex(0.2,0.83,label);
    //t->DrawLatex(0.2,0.77,"m(#tilde{q}) >> m(#tilde{g})");
    t->DrawLatex(0.2,0.71,signames.at(i).c_str());
    t->DrawLatex(0.15,0.92,"CMS Preliminary   #sqrt{s} = 8 TeV, #scale[0.6]{#int}Ldt = 19.5 fb^{-1}");
    */
    
    /*
    hjes[i]->GetYaxis()->SetTitle("#chi^{0}_{1} mass (GeV)");
    hjes[i]->GetXaxis()->SetTitle("#tilde{t} mass (GeV)");
    hjes[i]->GetZaxis()->SetTitle("JES uncertainty");
    hjes[i]->Draw("colz");

    t->DrawLatex(0.2,0.83,label);
    //t->DrawLatex(0.2,0.77,"m(#tilde{q}) >> m(#tilde{g})");
    t->DrawLatex(0.2,0.71,signames.at(i).c_str());
    t->DrawLatex(0.18,0.92,"CMS Preliminary            #sqrt{s} = 7 TeV, #scale[0.6]{#int}Ldt = 4.98 fb^{-1}");
    */

    if( print ){
      can[i]->          Print(Form("plots/%s%s%s_%s%s%s.pdf"       ,sample,xchar,suffix,labels.at(i).c_str(),pol,BDTchar));
      //can_exclusion[i]->Print(Form("../../plots/%s%s_%s_points%s.pdf",sample,suffix,labels.at(i).c_str(),pol));
    }

    int bin = heff[i]->FindBin(300,50);

    cout << "efficiency (300,50)  " << heff[i]->GetBinContent(bin)      << endl;
    cout << "nevents              " << hnevents[i]->GetBinContent(bin)  << endl;
    cout << "xsec UL              " << hxsec[i]->GetBinContent(bin)     << endl;
    cout << "xsec UL exp          " << hxsec_exp[i]->GetBinContent(bin) << endl;
    cout << "JES                  " << hjes[i]->GetBinContent(bin)      << endl;
    cout << "tot err              " << htoterr[i]->GetBinContent(bin)   << endl;
    cout << "btag err             " << hbtagerr[i]->GetBinContent(bin)  << endl;
    cout << "ISR err              " << hisrerr[i]->GetBinContent(bin)   << endl;
    cout << "stat err             " << hstaterr[i]->GetBinContent(bin)  << endl;
    cout << endl << endl;

    *logfile << "efficiency (300,50)  " << heff[i]->GetBinContent(bin)      << endl;
    *logfile << "nevents              " << hnevents[i]->GetBinContent(bin)  << endl;
    *logfile << "xsec UL              " << hxsec[i]->GetBinContent(bin)     << endl;
    *logfile << "xsec UL exp          " << hxsec_exp[i]->GetBinContent(bin) << endl;
    *logfile << "JES                  " << hjes[i]->GetBinContent(bin)      << endl;
    *logfile << "tot err              " << htoterr[i]->GetBinContent(bin)   << endl;
    *logfile << "btag err             " << hbtagerr[i]->GetBinContent(bin)  << endl;
    *logfile << "ISR err              " << hisrerr[i]->GetBinContent(bin)   << endl;
    *logfile << "stat err             " << hstaterr[i]->GetBinContent(bin)  << endl;
    *logfile << endl << endl;

  }


  TFile *outfile = TFile::Open(Form("%s%s%s%s%s_histos.root",sample,xchar,suffix,pol,BDTchar),"RECREATE");

  outfile->cd();
  for( unsigned int i = 0 ; i < nsig ; ++i ){
    hxsec[i]->Write();
    hxsec_exp[i]->Write();
    hxsec_expp1[i]->Write();
    hxsec_expm1[i]->Write();
    hexcl[i]->Write();
    hexcl_obsp1[i]->Write();
    hexcl_obsm1[i]->Write();
    hexcl_exp[i]->Write();
    hexcl_expp1[i]->Write();
    hexcl_expm1[i]->Write();
    heff[i]->Write();
    hnevents[i]->Write();
    hjes[i]->Write();
    htoterr[i]->Write();
    hisrerr[i]->Write();
    hbtagerr[i]->Write();
    hstaterr[i]->Write();
    gr[i]->Write();
    gr_exp[i]->Write();
  }
  outfile->Close();

}


void doAll(){

  // C&C
  //SMS("T2tt", 1,false,""     ,true); //unpolarized
  SMS("T2tt", 1,false,"left" ,true); //topL
  SMS("T2tt", 1,false,"right",true); //topR
  // SMS("T2bw",25,false,"",true);
  // SMS("T2bw",50,false,"",true);
  // SMS("T2bw",75,false,"",true);

  // BDT
  //SMS("T2tt", 1,true,""     ,true);
  SMS("T2tt", 1,true,"left" ,true);
  SMS("T2tt", 1,true,"right",true);
  // SMS("T2bw",25,true,"",true);
  // SMS("T2bw",50,true,"",true);
  // SMS("T2bw",75,true,"",true);
}
