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
//#include "upperLimits.C"

using namespace std;

void SMS(char* sample = "T2tt" , int x = 1, bool print = false){

  gStyle->SetPaintTextFormat(".2f");

  //--------------------------------------------------
  // input parameters
  //--------------------------------------------------

  char* suffix         = "";
  char* denomhistoname = "masses";

  if( TString(sample).Contains("T2bw") && x==25 ){
    denomhistoname = "masses25";
    suffix         = "_x25";
  }
  if( TString(sample).Contains("T2bw") && x==50 ){
    denomhistoname = "masses";
    suffix         = "_x50";
  }
  if( TString(sample).Contains("T2bw") && x==75 ){
    denomhistoname = "masses75";
    suffix         = "_x75";
  }

  const float lumi      = 19500;
  const char* filename  = "/tas/cms2/stop/cms2V05-03-26_stoplooperV00-02-22/T2tt_mad/minibabyV00-03-03/Skim_4jets_MET100_MT120/T2tt*root";
  const char* denomname = "/tas/cms2/stop/cms2V05-03-26_stoplooperV00-02-22/T2tt_mad/minibabyV00-03-03/Skim_4jets_MET100_MT120/myMassDB_T2tt_MG_25GeVbins.root";
  const bool  doBDT     = true;

  //const char* filename  = Form("/tas/benhoob/testFiles/%s_8TeV/merged*njets4%s.root",sample,suffix);
  //const char* denomname = Form("/tas/benhoob/testFiles/%s_8TeV/myMassDB.root",sample);
  //const float btagerr   = 0.02;

  cout << "----------------------------------------------------" << endl;
  cout << "Sample            " << sample          << endl;
  cout << "x                 " << x               << endl;
  cout << "Using file        " << filename        << endl;
  cout << "Using denominator " << denomname       << endl;
  cout << "Denom histo       " << denomhistoname  << endl;
  cout << "Using lumi        " << lumi            << " pb-1" << endl;
  cout << "----------------------------------------------------" << endl;

  TFile* fdenom = TFile::Open(denomname);
  TH2F*  hdenom = (TH2F*) fdenom->Get(denomhistoname);

  char* label = (char*)"";

  // T2tt
  if     ( TString(sample).Contains("T2tt") ){
    label = (char*)"pp #rightarrow #tilde{t}#tilde{t}, #tilde{t} #rightarrow t+#chi_{1}^{0}";
  }

  // T2bw, various x slices
  else if( TString(sample).Contains("T2bw") ){
    if( x == 25 ){
      cout << "Doing x=0.25 slice" << endl;
      label = (char*)"pp #rightarrow #tilde{t}#tilde{t}, #tilde{t} #rightarrow b+#chi_{1}^{#pm}, x=0.25";
    }
    else if( x == 50 ){
      cout << "Doing x=0.50 slice" << endl;
      label = (char*)"pp #rightarrow #tilde{t}#tilde{t}, #tilde{t} #rightarrow b+#chi_{1}^{#pm}, x=0.50";
    }
    else if( x == 75 ){
      cout << "Doing x=0.75 slice" << endl;
      label = (char*)"pp #rightarrow #tilde{t}#tilde{t}, #tilde{t} #rightarrow b+#chi_{1}^{#pm}, x=0.75";
    }
    else{
      cout << "ERROR! unrecognized x value " << x << endl;
      exit(0);
    }
  }
  else{
    cout << "ERROR! unrecognized sample " << sample << endl;
    exit(0);
  }

  //if( TString(sample).Contains("T2tt") ) label = "";

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
  TCut met100("t1metphicorr > 100");
  TCut mt120("t1metphicorrmt > 120");
  TCut dphi("mini_dphimjmin>0.8");
  TCut chi2("mini_chi2<5.0");
  TCut mt2w("mini_mt2w>200.0");

  // TCut testing("event%2==0");
  // TCut HM("mini_t2ttHM==1");
  // TCut LM("mini_t2ttLM==1");
  // TCut met50("t1metphicorr > 50");
  // TCut met150("t1metphicorr > 150");
  // TCut met200("t1metphicorr > 200");
  // TCut met250("t1metphicorr > 250");
  // TCut met300("t1metphicorr > 300");
  // TCut SRA("t1metphicorr > 100 && t1metphicorrmt > 150");
  // TCut SRB("t1metphicorr > 150 && t1metphicorrmt > 120");

  TCut presel;

  //-------------------------------------------
  // THESE CUTS DEFINE PRESELECTION REGION
  //-------------------------------------------
  presel += rho;
  presel += filters;
  presel += goodlep;
  presel += (el||mu);
  presel += njets4;
  presel += btag1;
  presel += passisotrk;
  presel += tauveto;
  presel += met100;
  presel += mt120;
  //presel += dphi;
  //presel += chi2;
  //presel += mt2w;

  TCut weight("mini_sltrigeff"); // trigger efficiency X btagging SF

  cout << "Using pre-selection   " << presel.GetTitle() << endl;
  cout << "Using weight          " << weight.GetTitle() << endl;

  //------------------------------------------
  // cut-based signal region definitions
  //------------------------------------------

  TCut   SR[8];
  string SRname[8]={"LM150","LM200","LM250","LM300","HM150","HM200","HM250","HM300"};

  SR[0]=TCut("mini_met > 150.0")+dphi+chi2;
  SR[1]=TCut("mini_met > 200.0")+dphi+chi2;
  SR[2]=TCut("mini_met > 250.0")+dphi+chi2;
  SR[3]=TCut("mini_met > 300.0")+dphi+chi2;

  SR[4]=TCut("mini_met > 150.0")+dphi+chi2+mt2w;
  SR[5]=TCut("mini_met > 200.0")+dphi+chi2+mt2w;
  SR[6]=TCut("mini_met > 250.0")+dphi+chi2+mt2w;
  SR[7]=TCut("mini_met > 300.0")+dphi+chi2+mt2w;

  //------------------------------------------
  // BDT signal region definitions
  //------------------------------------------

  TCut   SR_BDT[5];
  string SRname_BDT[5];

  SR_BDT[0]=TCut("mini_bdt[1] > 0.30"); SRname_BDT[0] = "BDT1L";
  SR_BDT[1]=TCut("mini_bdt[1] > 0.40"); SRname_BDT[1] = "BDT1T";
  SR_BDT[2]=TCut("mini_bdt[2] > 0.55"); SRname_BDT[2] = "BDT2";
  SR_BDT[3]=TCut("mini_bdt[3] > 0.65"); SRname_BDT[3] = "BDT3";
  SR_BDT[4]=TCut("mini_bdt[4] > 0.50"); SRname_BDT[4] = "BDT4";


  //SR_BDT[3]=TCut("mini_bdt[2] > 0.65"); SRname_BDT[3] = "BDT2_65";
  //SR_BDT[4]=TCut("mini_bdt[3] > 0.55"); SRname_BDT[4] = "BDT3_55";
  //SR_BDT[5]=TCut("mini_bdt[3] > 0.60"); SRname_BDT[5] = "BDT3_60";
  //SR_BDT[7]=TCut("mini_bdt[3] > 0.70"); SRname_BDT[7] = "BDT3_70";
  //SR_BDT[8]=TCut("mini_bdt[3] > 0.75"); SRname_BDT[8] = "BDT3_75";

  //--------------------------------------------------
  // signal regions
  //--------------------------------------------------

  vector<TCut>    sigcuts;
  vector<string>  signames;
  vector<string>  labels;
  vector<int>     cuts;

  if( doBDT ){
    cout << "Doing BDT signal regions" << endl;
    suffix = "_BDT";

    sigcuts.push_back(TCut(presel+SR_BDT[0]));  signames.push_back(SRname_BDT[0]);  labels.push_back(SRname_BDT[0]);  cuts.push_back(1);
    sigcuts.push_back(TCut(presel+SR_BDT[1]));  signames.push_back(SRname_BDT[1]);  labels.push_back(SRname_BDT[1]);  cuts.push_back(1);
    sigcuts.push_back(TCut(presel+SR_BDT[2]));  signames.push_back(SRname_BDT[2]);  labels.push_back(SRname_BDT[2]);  cuts.push_back(1);
    sigcuts.push_back(TCut(presel+SR_BDT[3]));  signames.push_back(SRname_BDT[3]);  labels.push_back(SRname_BDT[3]);  cuts.push_back(1);
    sigcuts.push_back(TCut(presel+SR_BDT[4]));  signames.push_back(SRname_BDT[4]);  labels.push_back(SRname_BDT[4]);  cuts.push_back(1);

    // sigcuts.push_back(TCut(presel+SR_BDT[5]));  signames.push_back(SRname_BDT[5]);  labels.push_back(SRname_BDT[5]);  cuts.push_back(1);
    // sigcuts.push_back(TCut(presel+SR_BDT[6]));  signames.push_back(SRname_BDT[6]);  labels.push_back(SRname_BDT[6]);  cuts.push_back(1);
    // sigcuts.push_back(TCut(presel+SR_BDT[7]));  signames.push_back(SRname_BDT[7]);  labels.push_back(SRname_BDT[7]);  cuts.push_back(1);
    // sigcuts.push_back(TCut(presel+SR_BDT[8]));  signames.push_back(SRname_BDT[8]);  labels.push_back(SRname_BDT[8]);  cuts.push_back(1);
    // sigcuts.push_back(TCut(presel+SR_BDT[9]));  signames.push_back(SRname_BDT[9]);  labels.push_back(SRname_BDT[9]);  cuts.push_back(1);
  }

  else{
    cout << "Doing cut-based signal regions" << endl;  

    // low-mass
    sigcuts.push_back(TCut(presel+SR[0]));  signames.push_back(SRname[0]);  labels.push_back(SRname[0]);  cuts.push_back(1);
    sigcuts.push_back(TCut(presel+SR[1]));  signames.push_back(SRname[1]);  labels.push_back(SRname[1]);  cuts.push_back(1);
    sigcuts.push_back(TCut(presel+SR[2]));  signames.push_back(SRname[2]);  labels.push_back(SRname[2]);  cuts.push_back(1);
    sigcuts.push_back(TCut(presel+SR[3]));  signames.push_back(SRname[3]);  labels.push_back(SRname[3]);  cuts.push_back(1);

    // high-mass
    sigcuts.push_back(TCut(presel+SR[4]));  signames.push_back(SRname[4]);  labels.push_back(SRname[4]);  cuts.push_back(1);
    sigcuts.push_back(TCut(presel+SR[5]));  signames.push_back(SRname[5]);  labels.push_back(SRname[5]);  cuts.push_back(1);
    sigcuts.push_back(TCut(presel+SR[6]));  signames.push_back(SRname[6]);  labels.push_back(SRname[6]);  cuts.push_back(1);
    sigcuts.push_back(TCut(presel+SR[7]));  signames.push_back(SRname[7]);  labels.push_back(SRname[7]);  cuts.push_back(1);
  }


  const unsigned int nsig = sigcuts.size();

  //--------------------------------------------------
  // make efficiency and xsec TH2's
  //--------------------------------------------------
  
  TH2F* heff[nsig];
  TH2F* heffup[nsig];
  TH2F* heffdn[nsig];
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
  
  TCanvas *ctemp = new TCanvas();
  ctemp->cd();

  for( unsigned int i = 0 ; i < nsig ; ++i ){

    cout << endl ;
    cout << "Signal Region : " << i << " " << signames.at(i) << endl;
    cout << "Selection     : " << sigcuts.at(i) << endl;

    //------------------------------------------------
    // turn off point-by-point signal uncertainties
    //------------------------------------------------

    /*
    TString jesup(sigcuts.at(i));
    jesup.ReplaceAll("npfjets30"       , "njetsUp");
    jesup.ReplaceAll("t1metphicorr "   , "t1metphicorrup ");
    jesup.ReplaceAll("t1metphicorrmt " , "t1metphicorrmtup ");

    TString jesdn(sigcuts.at(i));
    jesdn.ReplaceAll("npfjets30"       , "njetsDown");
    jesdn.ReplaceAll("t1metphicorr "   , "t1metphicorrdn ");
    jesdn.ReplaceAll("t1metphicorrmt " , "t1metphicorrmtdn ");

    TCut jesupcut(jesup);
    TCut jesdncut(jesdn);

    cout << endl << endl;
    cout << "Signal region : " << labels.at(i)  << endl;
    cout << "Selection     : " << sigcuts.at(i) << endl;
    cout << "Selection up  : " << jesupcut      << endl;
    cout << "Selection dn  : " << jesdncut      << endl;
    */

    int   nbinsx  =      41;
    float xmin    =   -12.5;
    float xmax    =  1012.5;
    int   nbinsy  =      41;
    float ymin    =   -12.5;
    float ymax    =  1012.5;

    // if( TString(sample).Contains("T2bw") ){
    //   nbinsx =   9;
    //   xmin   = 175;
    //   xmax   = 625;
    //   nbinsy =   7;
    //   ymin   = -25;
    //   ymax   = 325;
    // }

    heff[i]        = new TH2F(Form("heff_%i",i)          , Form("heff_%i",i)         , nbinsx , xmin , xmax , nbinsy , ymin , ymax );
    heffup[i]      = new TH2F(Form("heffup_%i",i)        , Form("heffup_%i",i)       , nbinsx , xmin , xmax , nbinsy , ymin , ymax );
    heffdn[i]      = new TH2F(Form("heffdn_%i",i)        , Form("heffdn_%i",i)       , nbinsx , xmin , xmax , nbinsy , ymin , ymax );
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

    hjes[i]      = new TH2F(Form("hjes_%i",i)        , Form("hjes_%i",i)       , nbinsx , xmin , xmax , nbinsy , ymin , ymax );

    ch->Draw(Form("ml:mg>>heff_%i",i),sigcuts.at(i)*weight);
    //heff[i]->Scale(1.0/denom);
    heff[i]->Divide(hdenom);

    /*
    ch->Draw(Form("ml:mg>>heffup_%i",i),jesupcut*weight);
    //heffup[i]->Scale(1.0/denom);
    heffup[i]->Divide(hdenom);

    ch->Draw(Form("ml:mg>>heffdn_%i",i),jesdncut*weight);
    //heffdn[i]->Scale(1.0/denom);
    heffdn[i]->Divide(hdenom);
    */

    for( unsigned int ibin = 1 ; ibin <= nbinsx ; ibin++ ){
      for( unsigned int jbin = 1 ; jbin <= nbinsy ; jbin++ ){

	float mg = heff[i]->GetXaxis()->GetBinCenter(ibin);
	//float ml = heff[i]->GetYaxis()->GetBinCenter(jbin);

	float eff    = heff[i]->GetBinContent(ibin,jbin);
	//float effup  = heffup[i]->GetBinContent(ibin,jbin);
	//float effdn  = heffdn[i]->GetBinContent(ibin,jbin);

	if( eff   < 1e-20 ) continue;

	// float dup    = effup/eff-1;
	// float ddn    = 1-effdn/eff;
	// float djes   = 0.5 * (dup+ddn);

	float djes = 0.10;

	hjes[i]->SetBinContent(ibin,jbin,djes);

	// lumi, lepton selection, trigger, b-tagging, JES

	float toterr  = 0.1;
	//float toterr  = sqrt( 0.04*0.04 + 0.02*0.02 + 0.03*0.03 + btagerr*btagerr + djes*djes );

	// float this_ul     = getObservedLimit( toterr , labels.at(i) );
	// float this_ul_exp = getExpectedLimit( toterr , labels.at(i) );

	// float xsecul      = this_ul / ( lumi * eff );
	// float xsecul_exp  = this_ul_exp / ( lumi * eff );

	float this_ul;
	float this_ul_exp;
	float this_ul_expp1;
	float this_ul_expm1;

	if( TString(labels.at(i)).Contains("LM150") ){
	  this_ul       = 79.6;
	  this_ul_exp   = this_ul;
	  this_ul_expp1 = this_ul;
	  this_ul_expm1 = this_ul;
	  // this_ul       = getUpperLimit_SRA( toterr );
	  // this_ul_exp   = getExpectedUpperLimit_SRA( toterr );
	  // this_ul_expp1 = getExpectedP1UpperLimit_SRA( toterr );
	  // this_ul_expm1 = getExpectedM1UpperLimit_SRA( toterr );
	}

	else if( TString(labels.at(i)).Contains("LM200") ){
	  this_ul       = 42.1;
	  this_ul_exp   = this_ul;
	  this_ul_expp1 = this_ul;
	  this_ul_expm1 = this_ul;
	}

	else if( TString(labels.at(i)).Contains("LM250") ){
	  this_ul       = 19.0;
	  this_ul_exp   = this_ul;
	  this_ul_expp1 = this_ul;
	  this_ul_expm1 = this_ul;
	}

	else if( TString(labels.at(i)).Contains("LM300") ){
	  this_ul       = 9.9;
	  this_ul_exp   = this_ul;
	  this_ul_expp1 = this_ul;
	  this_ul_expm1 = this_ul;
	}

	else if( TString(labels.at(i)).Contains("HM150") ){
	  this_ul       = 17.5;
	  this_ul_exp   = this_ul;
	  this_ul_expp1 = this_ul;
	  this_ul_expm1 = this_ul;
	}

	else if( TString(labels.at(i)).Contains("HM200") ){
	  this_ul       = 12.2;
	  this_ul_exp   = this_ul;
	  this_ul_expp1 = this_ul;
	  this_ul_expm1 = this_ul;
	}

	else if( TString(labels.at(i)).Contains("HM250") ){
	  this_ul       = 9.9;
	  this_ul_exp   = this_ul;
	  this_ul_expp1 = this_ul;
	  this_ul_expm1 = this_ul;
	}

	else if( TString(labels.at(i)).Contains("HM300") ){
	  this_ul       = 6.4;
	  this_ul_exp   = this_ul;
	  this_ul_expp1 = this_ul;
	  this_ul_expm1 = this_ul;
	}

	else if( TString(labels.at(i)).Contains("BDT1L") ){
	  this_ul       = 162.0;
	  this_ul_exp   = 162.0;
	  this_ul_expp1 = 162.0;
	  this_ul_expm1 = 162.0;
	}

	else if( TString(labels.at(i)).Contains("BDT1_40") ){
	  this_ul       = 46.0;
	  this_ul_exp   = 46.0;
	  this_ul_expp1 = 46.0;
	  this_ul_expm1 = 46.0;
	}

	else if( TString(labels.at(i)).Contains("BDT2_55") ){
	  this_ul       = 38.0;
	  this_ul_exp   = 38.0;
	  this_ul_expp1 = 38.0;
	  this_ul_expm1 = 38.0;
	}

	else if( TString(labels.at(i)).Contains("BDT2_65") ){
	  this_ul       = 23.0;
	  this_ul_exp   = 23.0;
	  this_ul_expp1 = 23.0;
	  this_ul_expm1 = 23.0;
	}

	else if( TString(labels.at(i)).Contains("BDT3_55") ){
	  this_ul       = 17.0;
	  this_ul_exp   = 17.0;
	  this_ul_expp1 = 17.0;
	  this_ul_expm1 = 17.0;
	}

	else if( TString(labels.at(i)).Contains("BDT3_60") ){
	  this_ul       = 11.4;
	  this_ul_exp   = 11.4;
	  this_ul_expp1 = 11.4;
	  this_ul_expm1 = 11.4;
	}

	else if( TString(labels.at(i)).Contains("BDT3_65") ){
	  this_ul       = 8.0;
	  this_ul_exp   = 8.0;
	  this_ul_expp1 = 8.0;
	  this_ul_expm1 = 8.0;
	}

	else if( TString(labels.at(i)).Contains("BDT3_70") ){
	  this_ul       = 5.9;
	  this_ul_exp   = 5.9;
	  this_ul_expp1 = 5.9;
	  this_ul_expm1 = 5.9;
	}

	else if( TString(labels.at(i)).Contains("BDT3_75") ){
	  this_ul       = 5.0;
	  this_ul_exp   = 5.0;
	  this_ul_expp1 = 5.0;
	  this_ul_expm1 = 5.0;
	}
	
	else if( TString(labels.at(i)).Contains("BDT4_50") ){
	  this_ul       = 3.7;
	  this_ul_exp   = 3.7;
	  this_ul_expp1 = 3.7;
	  this_ul_expm1 = 3.7;
	}

	else{
	  cout << "WARNING UNRECOGNIZED SIGNAL REGION " << labels.at(i) << " QUITTING!!!" << endl;
	  exit(0);
	}


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

    gr[i]      = getRefXsecGraph(hxsec[i]     , "T2tt", 1.0);
    gr_exp[i]  = getRefXsecGraph(hxsec_exp[i] , "T2tt", 1.0);

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
    heff[i]->GetXaxis()->SetRangeUser(100,700);
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
    hxsec[i]->GetZaxis()->SetTitle("#sigma upper limit");
    hxsec[i]->GetZaxis()->SetTitleOffset(1.2);
    hxsec[i]->Draw("colz");
    //hxsec[i]->Draw("sametext");
    hxsec[i]->SetMinimum(0.01);
    hxsec[i]->SetMaximum(100);
    hxsec[i]->GetXaxis()->SetRangeUser(200,700);
    hxsec[i]->GetYaxis()->SetRangeUser(0,600);

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
      can[i]->          Print(Form("../../plots/%s%s_%s.pdf"       ,sample,suffix,labels.at(i).c_str()));
      can_exclusion[i]->Print(Form("../../plots/%s%s_%s_points.pdf",sample,suffix,labels.at(i).c_str()));
    }

    int bin = heff[i]->FindBin(300,50);

    float toterr = 0.1; //sqrt(pow(hjes[i]->GetBinContent(bin),2)+0.04*0.04 + 0.02*0.02 + 0.03*0.03 + btagerr*btagerr);
    //float toterr = sqrt(pow(hjes[i]->GetBinContent(bin),2)+0.04*0.04 + 0.02*0.02 + 0.03*0.03 + btagerr*btagerr);
    cout << "efficiency (300,50)  " << heff[i]->GetBinContent(bin) << endl;
    cout << "xsec UL              " << hxsec[i]->GetBinContent(bin) << endl;
    cout << "xsec UL exp          " << hxsec_exp[i]->GetBinContent(bin) << endl;
    cout << "JES                  " << hjes[i]->GetBinContent(bin) << endl;
    cout << "tot err              " << toterr << endl;
    //cout << "obs limit            " << getObservedLimit(toterr,labels.at(i)) << endl;
    //cout << "exp limit            " << getExpectedLimit(toterr,labels.at(i)) << endl;
    cout << endl << endl;
  }
  
  TFile *outfile = TFile::Open(Form("%s%s_histos.root",sample,suffix),"RECREATE");

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
    hjes[i]->Write();
    gr[i]->Write();
    gr_exp[i]->Write();
  }
  outfile->Close();

}
