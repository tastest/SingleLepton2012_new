#include <TH1F.h>
#include <TH2F.h>
#include <TLegend.h>
#include <TFile.h>
#include <TMath.h>
#include <TLine.h>
#include <TLatex.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <vector>
#include <iostream>
#include <algorithm>
#include <sstream>
#include <iomanip>
using namespace std;

float addSqr( float e1 , float e2 ){

  float err2 = pow(e1,2) + pow(e2,2);

  return sqrt(err2);
}

float subtractSqr( float e1 , float e2 ){

  float err2 = pow(e1,2) - pow(e2,2);

  return sqrt(err2);
}

float histError( TH1F* h , int minbin , int maxbin ){

  float err2 = 0;
  for( int i = minbin ; i <= maxbin ; ++i ){
    err2 += pow( h->GetBinError(i) , 2 );
  }

  return sqrt(err2);
}

void zeroHistError( TH1F* &h ) {

  for (int i=1; i<=h->GetNbinsX(); ++i) 
    h->SetBinError(i,0.); 
}

void zeroHist( TH1F* &h ) {

  for (int i=1; i<=h->GetNbinsX(); ++i) {
    h->SetBinContent(i,0.); 
    h->SetBinError(i,0.); 
  }
}

double error_poisson_up(double data) {
  double y1 = data + 1.0;
  double d = 1.0 - 1.0/(9.0*y1) + 1.0/(3*TMath::Sqrt(y1));
  return y1*d*d*d-data;
}

double error_poisson_down(double data) {
  double y = data;
  if (y == 0.0) return 0.0;
  double d = 1.0 - 1.0/(9.0*y) - 1.0/(3.0*TMath::Sqrt(y));
  return data-y*d*d*d;
}

//settings go here
const int NSAMPLE = 3;

char * selection[NSAMPLE] = {"_1l1or2b_mt150_g4j", "_1l1or2b_mt150_g5j", "_1l1or2b_mt150_g6j"};
const char* samplename[NSAMPLE] = {"1l + 2b + #geq 4j (M_{T}>150)", "1l + 2b + #geq 5j (M_{T}>150)", "1l + 2b + #geq 6j (M_{T}>150)"};
char * printtag = "_mt150";
int crbin = 1;
float errfactor = 0.1;

// char * selection[NSAMPLE] = {"_1l1or2b_mt120_g4j", "_1l1or2b_mt120_g5j", "_1l1or2b_mt120_g6j"};
// const char* samplename[NSAMPLE] = {"1l + 2b + #geq 4j (M_{T}>120)", "1l + 2b + #geq 5j (M_{T}>120)", "1l + 2b + #geq 6j (M_{T}>120)"};
// char * printtag = "_mt120";
// int crbin = 0;
// float errfactor = 0.05;

bool addhighmbb[NSAMPLE] = {false, false, false};
char * histotag = "";

bool printplot = false;

TFile *dt_dl;
TFile *dt_sl;

void plotComparison( TH1F* h_dt , TH1F* h_mc , TH1F *h_extra, char* label, bool dolog, bool drawbkg);
void runSimplePred(char* ttbar_tag = "", bool issyst = false);
//////////////////////////////////////////////////////////////////////////////

//ttbar prediction and uncertainty
float ttpred[NSAMPLE]; 
float ettpred[NSAMPLE];
float ttdlpred[NSAMPLE]; 
float ettdlpred[NSAMPLE];
//rare prediction and uncertainty
float rarepred[NSAMPLE]; 
float erarepred[NSAMPLE];
//total prediction and uncertainty
float pred[NSAMPLE];
float epred[NSAMPLE];
float sf_rmbb[NSAMPLE];
float esf_rmbb[NSAMPLE];

float xs_unc_rare = 0.;
float xs_unc_ttdl = 1.;

void doSimplePredNJCheckTTSplit(bool doaltttbar = false) {

  dt_dl = TFile::Open("SIGoutput/data_dl_histos.root");
  dt_sl = TFile::Open("SIGoutput/data_sl_histos.root");
  
  for (int isr=0; isr<NSAMPLE; ++isr) {
    ttpred[isr] = -999.;
    ettpred[isr] = -999.;
    ttdlpred[isr] = -999.;
    ettdlpred[isr] = -999.;
    rarepred[isr] = -999.;
    erarepred[isr] = -999.;
    pred[isr] = -999.;
    epred[isr] = -999.;
    sf_rmbb[isr] = -999.;
    esf_rmbb[isr] = -999.;
  }

  runSimplePred("_lmgtau", false);

  float sf_rmbb_default[NSAMPLE];
  float esf_rmbb_default[NSAMPLE]; 
  std::copy(&sf_rmbb[0], &sf_rmbb[NSAMPLE], sf_rmbb_default);
  std::copy(&esf_rmbb[0], &esf_rmbb[NSAMPLE], esf_rmbb_default);
  
  xs_unc_ttdl = 1.1;
  
  runSimplePred("_lmgtau", false);
  
  float sf_rmbb_ttdlup[NSAMPLE];
  float esf_rmbb_ttdlup[NSAMPLE]; 
  std::copy(&sf_rmbb[0], &sf_rmbb[NSAMPLE], sf_rmbb_ttdlup);
  std::copy(&esf_rmbb[0], &esf_rmbb[NSAMPLE], esf_rmbb_ttdlup);
  
  float aunc_ttdlup[NSAMPLE]; 
  float runc_ttdlup[NSAMPLE]; 
  for (int isr=0; isr<NSAMPLE; ++isr) {
    aunc_ttdlup[isr] = fabs(sf_rmbb_ttdlup[isr]-sf_rmbb_default[isr]);
    runc_ttdlup[isr] = aunc_ttdlup[isr]/sf_rmbb_default[isr]*100.;
  }
  
  xs_unc_ttdl = 0.9;
  
  runSimplePred("_lmgtau", false);
  
  float sf_rmbb_ttdldn[NSAMPLE];
  float esf_rmbb_ttdldn[NSAMPLE]; 
  std::copy(&sf_rmbb[0], &sf_rmbb[NSAMPLE], sf_rmbb_ttdldn);
  std::copy(&esf_rmbb[0], &esf_rmbb[NSAMPLE], esf_rmbb_ttdldn);
  
  float aunc_ttdldn[NSAMPLE]; 
  float runc_ttdldn[NSAMPLE]; 
  for (int isr=0; isr<NSAMPLE; ++isr) {
    aunc_ttdldn[isr] = fabs(sf_rmbb_ttdldn[isr]-sf_rmbb_default[isr]);
    runc_ttdldn[isr] = aunc_ttdldn[isr]/sf_rmbb_default[isr]*100.;
  }
  
  cout<<endl;
  cout<<endl;
  cout<<endl;
  cout<<"*************************************************"<<endl;
  cout<<"***** RESULTS OF BKG SUBTRACTION VARIATIONS *****"<<endl;
  cout<<"*************************************************"<<endl;

  printf("ttdl up \t ");
  for (int isr=0; isr<NSAMPLE; ++isr) 
    printf("\t & $%.1f$ %s", runc_ttdlup[isr], isr==NSAMPLE-1 ? " \\\\ \n":"");
  printf("ttdl down \t ");
  for (int isr=0; isr<NSAMPLE; ++isr) 
    printf("\t & $%.1f$ %s", runc_ttdldn[isr], isr==NSAMPLE-1 ? " \\\\ \n":"");
  
  cout<<endl;

  xs_unc_ttdl = 1.;
  printplot = true;

  runSimplePred("_lmgtau", false);

}

void runSimplePred(char* ttbar_tag, bool issyst){
  
  // this flag sets all uncertainties except for those
  // associated with the ttbar sample to zero
  // in this way it calculates the uncorrelated uncertainty
  //  bool issyst = true;
  //  bool issyst = (strlen(ttbar_tag) == 0) ? false : true;

  // to dump all the information about the yields for 
  // the various MC samples
  bool doverbose = false;

  //for single lepton, split single top
  const int MCID = 3; 
  const char* mcsample[MCID]={
    "ttsl",
    "ttdl",
    "rare"};

  //legends only needed for final predictions
  const char* legend[MCID] = {
    "t#bar{t}#rightarrow #font[12]{l} + jets",
    "t#bar{t}#rightarrow #font[12]{ll}",
    "Rare"};

  enum sample{TTSL=0, 
	      TTDL,
	      RARE};

  // bins for regions
  float i_ctr = 1;
  float i_sig = 2;

  //Store raw information
  //Open single lepton files
  TFile *mc_dl[MCID];
  for (int j=0;j<MCID;++j) {
    if (j<2)
      mc_dl[j] = TFile::Open(Form("SIGoutput/%s%s_histos.root",
				  mcsample[j], ttbar_tag));
    else 
      mc_dl[j] = TFile::Open(Form("SIGoutput/%s_histos.root",
				  mcsample[j]));
  }

  char * histoname = "h_sig_mt_count";

  //indices are for samples
  TH1F *h_dt[NSAMPLE];
  TH1F *h_mc[MCID][NSAMPLE];
  TH1F *h_mc_tot[NSAMPLE];
  TH1F *h_dt_cor[NSAMPLE];//data corrected for rare backgrounds 

  for (int isr=0; isr<NSAMPLE; ++isr) {

    //data
    h_dt[isr]= (TH1F*)dt_sl->Get(Form("%s%s%s",histoname,histotag,selection[isr]));
    if (doverbose) {
      cout<<"Data "<<Form("%s%s%s",histoname,histotag,selection[isr])
	  <<" "<<h_dt[isr]<<endl;
      cout<<"-------------------------------------------------"<<endl;
      printf("Data Yields for selection %s%s%s \n",histoname,histotag,selection[isr]);
      printf("Control: %.2f pm %.2f (%.1f %%) \n", 
	     h_dt[isr]->GetBinContent(i_ctr), 
	     h_dt[isr]->GetBinError(i_ctr), 
	     h_dt[isr]->GetBinError(i_ctr)/
	     h_dt[isr]->GetBinContent(i_ctr)*100.);
      printf("Signal: %.2f pm %.2f (%.1f %%) \n", 
      	     h_dt[isr]->GetBinContent(i_sig), 
      	     h_dt[isr]->GetBinError(i_sig), 
      	     h_dt[isr]->GetBinError(i_sig)/
      	     h_dt[isr]->GetBinContent(i_sig)*100.);
    }

    h_dt[isr]->SetName(Form("h_dt%s",selection[isr]));
    if (issyst)	zeroHistError(h_dt[isr]);
    h_dt_cor[isr] = (TH1F*)h_dt[isr]->Clone();
    h_dt_cor[isr]->SetName(Form("h_dt_corr%s",selection[isr]));

    if (addhighmbb[isr]) {
      float binval  =  h_dt[isr]->GetBinContent(1) + h_dt[isr]->GetBinContent(3);
      float ebinval =  addSqr(h_dt[isr]->GetBinError(1), h_dt[isr]->GetBinError(3));
      h_dt[isr]->SetBinContent(1, binval);
      h_dt[isr]->SetBinError(1, ebinval);
      h_dt[isr]->SetBinContent(3, 0.);
      h_dt[isr]->SetBinError(3, 0.);
    } 

    //mc
    for (int j=0;j<MCID;++j) {
      h_mc[j][isr] = (TH1F*)mc_dl[j]->Get(Form("%s%s%s",histoname,histotag,selection[isr]));
      if (doverbose) 
	cout<<"MC "<<mcsample[j]<<" "
	    <<Form("%s%s%s",histoname,histotag,selection[isr])
	    <<" "<<h_mc[j][isr]<<endl;
      if (h_mc[j][isr]==0) {
	h_mc[j][isr] = (TH1F*)mc_dl[0]->Get(Form("%s%s%s",histoname,histotag,selection[isr]));
	h_mc[j][isr]->SetName(Form("h_mc_%s%s", mcsample[j],selection[isr]));
	zeroHist(h_mc[j][isr]);
      }
	
      if (doverbose) 
	cout<<"MC "<<mcsample[j]<<" "
	    <<Form("%s%s%s",histoname,histotag,selection[isr])
	    <<" "<<h_mc[j][isr]<<endl;
      h_mc[j][isr]->SetName(Form("h_mc_%s%s", mcsample[j],selection[isr]));
	
      if (issyst && j>1)
	zeroHistError(h_mc[j][isr]);

      if (addhighmbb[isr]) {
	float bin  =  h_mc[j][isr]->GetBinContent(1) + h_mc[j][isr]->GetBinContent(3);
	float ebin =  addSqr(h_mc[j][isr]->GetBinError(1), h_mc[j][isr]->GetBinError(3));
	h_mc[j][isr]->SetBinContent(1, bin);
	h_mc[j][isr]->SetBinError(1, ebin);
	h_mc[j][isr]->SetBinContent(3, 0.);
	h_mc[j][isr]->SetBinError(3, 0.);
      } 
      
      if (j==TTDL||j==RARE) 
	h_mc[j][isr]->Scale(xs_unc_ttdl);

      if (j==0) {
	h_mc_tot[isr] = (TH1F*)h_mc[j][isr]->Clone();
	h_mc_tot[isr]->SetName(Form("h_mc_tot%s",selection[isr]));
      } else 
	h_mc_tot[isr]->Add(h_mc[j][isr]);
	
      if (j!=TTSL) 
	h_dt_cor[isr]->Add(h_mc[j][isr], -1.);

      if (doverbose) {
	cout<<"-------------------------------------------------"<<endl;
	printf("%s MC Yields for selection %s \n", mcsample[j], selection[isr]);
	printf("Control: %.2f pm %.2f (%.1f %%) \n", 
	       h_mc[j][isr]->GetBinContent(i_ctr), 
	       h_mc[j][isr]->GetBinError(i_ctr), 
	       h_mc[j][isr]->GetBinError(i_ctr)/
	       h_mc[j][isr]->GetBinContent(i_ctr)*100.);
	printf("Signal: %.2f pm %.2f (%.1f %%) \n", 
	       h_mc[j][isr]->GetBinContent(i_sig), 
	       h_mc[j][isr]->GetBinError(i_sig), 
	       h_mc[j][isr]->GetBinError(i_sig)/
	       h_mc[j][isr]->GetBinContent(i_sig)*100.);
      }
    }
    
  }


  //Calculate the SFs
  cout<<"-------------------------------------------------"<<endl;
  cout<<"\\hline"<<endl;  
  printf("\\multicolumn{%i}{c}{$M_{b\\bar{b}} \\le 100$~GeV} \\\\\n", NSAMPLE+1);
  cout<<"\\hline"<<endl;  
  cout<<"Sample  ";
  for (int isr=0; isr<NSAMPLE; ++isr) 
    printf("& %s ", samplename[isr]);
  printf(" \\\\\n");
  cout<<"\\hline"<<endl;
  cout<<"Data$\_{\\rm peak}$  ";
  for (int isr=0; isr<NSAMPLE; ++isr) {
    printf("& $%.2f \\pm %.2f$ ",
	   h_dt_cor[isr]->GetBinContent(i_ctr), 
	   h_dt_cor[isr]->GetBinError(i_ctr));
  }
  printf(" \\\\\n");
  cout<<"\\hline"<<endl;
  cout<<"MC$\_{\\rm peak}$  ";
  for (int isr=0; isr<NSAMPLE; ++isr) {
    printf("& $%.2f \\pm %.2f$ ",
	   h_mc[TTSL][isr]->GetBinContent(i_ctr), 
	   h_mc[TTSL][isr]->GetBinError(i_ctr));
  }
  printf(" \\\\\n");
  cout<<"\\hline"<<endl;
  cout<<"\\hline"<<endl;
  
  TH1F *h_sfout[NSAMPLE];
  cout<<"SF$\_{\\rm peak}$ ";
  for (int isr=0; isr<NSAMPLE; ++isr) {
    h_sfout[isr] = (TH1F*)h_dt_cor[isr]->Clone();
    h_sfout[isr]->SetName(Form("h_sfout%s",selection[isr]));
    h_sfout[isr]->Divide(h_mc[TTSL][isr]);
    printf("& $%.2f \\pm %.2f$ ",
	   h_sfout[isr]->GetBinContent(i_ctr), 
	   h_sfout[isr]->GetBinError(i_ctr));
  }
  printf(" \\\\\n");
  cout<<"\\hline"<<endl;
  cout<<"\\hline"<<endl;
  
  cout<<endl;
  cout<<endl;
  cout<<endl;
  
  cout<<"\\hline"<<endl;
  cout<<"MC$\_{\\rm tail}$  ";
  for (int isr=0; isr<NSAMPLE; ++isr) {
    printf("& $%.2f \\pm %.2f$ ",
	   h_mc[TTSL][isr]->GetBinContent(i_sig), 
	   h_mc[TTSL][isr]->GetBinError(i_sig));
  }
  printf(" \\\\\n");
  cout<<"MC$\_{\\rm peak}$  ";
  for (int isr=0; isr<NSAMPLE; ++isr) {
    printf("& $%.2f \\pm %.2f$ ",
	   h_mc[TTSL][isr]->GetBinContent(i_ctr), 
	   h_mc[TTSL][isr]->GetBinError(i_ctr));
  }
  printf(" \\\\\n");
  cout<<"\\hline"<<endl;
  
  float rmbb[NSAMPLE];
  float ermbb[NSAMPLE];
  cout<<"${\\rm R}_{{\\rm M}_{\\rm T}}^{\\rm MC}$  ";
  for (int isr=0; isr<NSAMPLE; ++isr) {
    rmbb[isr] = h_mc[TTSL][isr]->GetBinContent(i_sig)/h_mc[TTSL][isr]->GetBinContent(i_ctr);
    ermbb[isr] = addSqr( (h_mc[TTSL][isr]->GetBinError(i_sig)/h_mc[TTSL][isr]->GetBinContent(i_sig)), 
			 (h_mc[TTSL][isr]->GetBinError(i_ctr)/h_mc[TTSL][isr]->GetBinContent(i_ctr)) )*rmbb[isr];
    printf("& $%.3f \\pm %.3f$ ",rmbb[isr], ermbb[isr]);
  }
  printf(" \\\\\n");
  cout<<"\\hline"<<endl;
  cout<<"\\hline"<<endl;

  printf("Data$\_{\\rm tail}$ ");
  for (int isr=0; isr<NSAMPLE; ++isr) 
    printf("& $%.2f \\pm %.2f$ ",
	 h_dt_cor[isr]->GetBinContent(i_sig), 
	 h_dt_cor[isr]->GetBinError(i_sig));
  
  // printf("Data$\_{\\rm tail}$ & $%.2f \\pm %.2f$",
  // 	 h_dt_cor[0]->GetBinContent(i_sig), 
  // 	 h_dt_cor[0]->GetBinError(i_sig));
  // for (int isr=0; isr<NSAMPLE; ++isr) 
  //   printf("& - ");
  printf(" \\\\\n");
  cout<<"Data$\_{\\rm peak}$  ";
  for (int isr=0; isr<NSAMPLE; ++isr) {
    printf("& $%.2f \\pm %.2f$ ",
	   h_dt_cor[isr]->GetBinContent(i_ctr), 
	   h_dt_cor[isr]->GetBinError(i_ctr));
  }
  printf(" \\\\\n");
  cout<<"\\hline"<<endl;
  cout<<"${\\rm R}_{{\\rm M}_{\\rm T}}^{\\rm Data}$  ";
  float rmbb_dt[NSAMPLE];
  float ermbb_dt[NSAMPLE];
  for (int isr=0; isr<NSAMPLE; ++isr) {
    rmbb_dt[isr] = h_dt_cor[isr]->GetBinContent(i_sig)/h_dt_cor[isr]->GetBinContent(i_ctr);
    ermbb_dt[isr] = 
      addSqr( (h_dt_cor[isr]->GetBinError(i_sig)/h_dt_cor[isr]->GetBinContent(i_sig)),
	      (h_dt_cor[isr]->GetBinError(i_ctr)/h_dt_cor[isr]->GetBinContent(i_ctr)) )*rmbb_dt[isr];
    printf("& $%.3f \\pm %.3f$ ",rmbb_dt[isr], ermbb_dt[isr]);
  }
  // for (int isr=1; isr<NSAMPLE; ++isr) 
  //   printf("& - ");
  printf(" \\\\\n");
  cout<<"\\hline"<<endl;
  cout<<"\\hline"<<endl;
  
  cout<<"SF(${\\rm R}_{{\\rm M}_{\\rm T}}$)  ";
  
  for (int isr=0; isr<NSAMPLE; ++isr) {
    sf_rmbb[isr] = rmbb_dt[isr]/rmbb[isr];
    esf_rmbb[isr] = 
      addSqr( (ermbb_dt[isr]/rmbb_dt[isr]),(ermbb[isr]/rmbb[isr]) )*sf_rmbb[isr]; 
    printf("& $%.2f \\pm %.2f$ ",sf_rmbb[isr], esf_rmbb[isr]);
  }
  printf(" \\\\\n");
  cout<<"\\hline"<<endl;
  cout<<"\\hline"<<endl;
  
  cout<<"Top pred. = Data$\_{\\rm peak}$ * ${\\rm R}\_{{\\rm M}\_{\\rm T}}^{\\rm MC}$ ";
  
  for (int isr=0; isr<NSAMPLE; ++isr) {
    //    if (isr==0) {
      ttpred[isr]  = h_dt_cor[isr]->GetBinContent(i_ctr)*rmbb[isr];
      ettpred[isr] = sqrt( pow( (h_dt_cor[isr]->GetBinError(i_ctr)/h_dt_cor[isr]->GetBinContent(i_ctr)), 2) + 
			   pow( (ermbb[isr]/rmbb[isr]), 2) ) * ttpred[isr];
    // } else {
    //   ttpred[isr]  = h_dt_cor[isr]->GetBinContent(i_ctr)*rmbb[isr]*sf_rmbb;
    //   ettpred[isr] = sqrt( pow( (h_dt_cor[isr]->GetBinError(i_ctr)/h_dt_cor[isr]->GetBinContent(i_ctr)), 2) + 
    // 			   pow( (ermbb[isr]/rmbb[isr]), 2) + pow( (esf_rmbb/sf_rmbb),2) ) * ttpred[isr];
    // }
    printf("& $%.2f \\pm %.2f$ ",ttpred[isr], ettpred[isr]);    
  }
  printf(" \\\\\n");
  cout<<"Top Dilepton ";
  
  for (int isr=0; isr<NSAMPLE; ++isr) {
    ttdlpred[isr]  = h_mc[TTDL][isr]->GetBinContent(i_sig);
    ettdlpred[isr] = h_mc[TTDL][isr]->GetBinError(i_sig);
    //    ettdlpred[isr] = addSqr( h_mc[TTDL][isr]->GetBinError(i_sig), (0.1*ttdlpred[isr]) );
    printf("& $%.2f \\pm %.2f$ ",ttdlpred[isr], ettdlpred[isr]);    
  }
  printf(" \\\\\n");
  
  cout<<"\\hline"<<endl;
  cout<<"Rare ";
  
  for (int isr=0; isr<NSAMPLE; ++isr) {
    rarepred[isr]  = h_mc[RARE][isr]->GetBinContent(i_sig);
    erarepred[isr] = addSqr( h_mc[RARE][isr]->GetBinError(i_sig), (xs_unc_rare*rarepred[isr]) );
    printf("& $%.2f \\pm %.2f$ ",rarepred[isr], erarepred[isr]);    
  }
  printf(" \\\\\n");
  
  cout<<"\\hline"<<endl;
  cout<<"Total Pred. ";
  for (int isr=0; isr<NSAMPLE; ++isr) {
    pred[isr] = ttpred[isr] + ttdlpred[isr] + rarepred[isr];
    epred[isr] = addSqr( ettpred[isr], ettdlpred[isr]);
    epred[isr] = addSqr( epred[isr], erarepred[isr]);
    printf("& $%.2f \\pm %.2f$ ",pred[isr], epred[isr]);    
  }
  printf(" \\\\\n");
  cout<<"\\hline"<<endl;
  cout<<"Data ";
  for (int isr=0; isr<NSAMPLE; ++isr) {
    printf("& $%.0f$ ",h_dt[isr]->GetBinContent(i_sig));    
  }
  printf(" \\\\\n"); 
  cout<<"\\hline"<<endl;
  cout<<"Data/Pred ";
  float closure[NSAMPLE];
  float eclosure[NSAMPLE];  
  for (int isr=0; isr<NSAMPLE; ++isr) {
    closure[isr] = h_dt[isr]->GetBinContent(i_sig)/pred[isr];
    eclosure[isr] = h_dt[isr]->GetBinContent(i_sig)>0. ? 
      addSqr(h_dt[isr]->GetBinError(i_sig)/h_dt[isr]->GetBinContent(i_sig), 
	     epred[isr]/pred[isr]) * closure[isr] : 0.;
    printf("& $%.2f \\pm %.2f$ ",closure[isr], eclosure[isr]);    
  }
  printf(" \\\\\n"); 
  
  cout<<endl;
  cout<<endl;
  cout<<endl;

  //make plot of results
  float x[NSAMPLE];
  for (int i=0; i<NSAMPLE; ++i) 
    x[i] = i*1+0.5;
  float xerr[NSAMPLE];
  for (int i=0; i<NSAMPLE; ++i) 
    xerr[i] = 0.;
    
  float err[NSAMPLE];
  for (int i=0; i<NSAMPLE; ++i) 
    err[i]=sf_rmbb[crbin];

  TGraphErrors* gclose = new TGraphErrors(NSAMPLE,x,sf_rmbb,xerr,esf_rmbb);

  TCanvas *c1 = new TCanvas();
  c1->cd();
  gPad->SetGridx();
  gPad->SetGridy();

  gclose->SetLineColor(2);
  gclose->SetMarkerColor(2);

  gclose->SetMinimum(0);

  TH2F* hdummy = new TH2F("hdummy","",NSAMPLE,0,NSAMPLE,100,0.8,1.3);
  for (int i=0; i<NSAMPLE; ++i) 
    hdummy->GetXaxis()->SetBinLabel(i+1, samplename[i]);
  // hdummy->GetXaxis()->SetBinLabel(1,"BDT1 Loose");
  // hdummy->GetXaxis()->SetBinLabel(2,"BDT1 Tight");
  // hdummy->GetXaxis()->SetBinLabel(3,"BDT2");
  // hdummy->GetXaxis()->SetBinLabel(4,"BDT3");
  // hdummy->GetXaxis()->SetBinLabel(5,"BDT4");
  // hdummy->GetXaxis()->SetBinLabel(6,"BDT5");

  hdummy->GetXaxis()->SetTitle("control region");
  hdummy->GetYaxis()->SetTitle("SF(R_{M_{T}})");
  //  hdummy->GetYaxis()->SetTitle("data / MC (M_{T} tail)");

  hdummy->Draw();
  gclose->Draw("sameP");

  TH1F* hct = new TH1F("hct","",NSAMPLE,0,NSAMPLE);
  for (int i=0; i<NSAMPLE; ++i) 
    hct->SetBinContent(i+1, sf_rmbb[crbin]);

  TH1F* hup = new TH1F("hup","",NSAMPLE,0,NSAMPLE);
  for (int i=0; i<NSAMPLE; ++i) 
    hup->SetBinContent(i+1, sf_rmbb[crbin]+errfactor);
  //    hup->SetBinContent(i+1,1+err[i]);

  TH1F* hdn = new TH1F("hdn","",NSAMPLE,0,NSAMPLE);
  for (int i=0; i<NSAMPLE; ++i) 
    hdn->SetBinContent(i+1,sf_rmbb[crbin]-errfactor);

  hct->SetLineWidth(3);
  hct->SetLineColor(kGray);
  hup->SetLineWidth(3);
  hdn->SetLineWidth(3);


  hct->Draw("samehist");
  hup->Draw("samehist");
  hdn->Draw("samehist");

  gclose->Draw("sameP");

  if (printplot) {
    c1->Print(Form("plots/njetscheck%s.pdf",printtag));
    c1->Print(Form("plots/njetscheck%s.png",printtag));
  }

}

void plotComparison( TH1F* h_dt , TH1F* h_mc , TH1F *h_extra, char* label, bool dolog, bool drawbkg) {

  TPad* fullpad = new TPad();
  TPad* plotpad = new TPad();
  TPad* respad  = new TPad();
  fullpad = new TPad("fullpad","fullpad",0,0,1,1);
  fullpad->Draw();
  fullpad->cd();
  plotpad = new TPad("plotpad","plotpad",0,0,1,0.8);
  plotpad->Draw();
  plotpad->cd();
  if (dolog) plotpad->SetLogy();

  h_dt->GetYaxis()->SetTitle("Entries");
  h_dt->GetXaxis()->SetTitle(Form("%s", label));
  h_dt->GetYaxis()->SetTitleSize(0.05);
  h_dt->GetXaxis()->SetTitleSize(0.05);
  h_dt->GetYaxis()->SetTitleOffset(1.5);
  h_dt->GetXaxis()->SetTitleOffset(1.3);
  if (!dolog) h_dt->GetYaxis()->SetRangeUser(0., 1.4*h_dt->GetMaximum());
  h_dt->SetLineColor(kBlack);
  h_dt->SetMarkerColor(kBlack);
  h_mc->SetLineColor(kBlue);
  h_mc->SetMarkerColor(kBlue);
  h_extra->SetLineColor(kRed);
  h_extra->SetLineWidth(2);
  h_mc->SetLineWidth(2);
  h_dt->Draw();
  h_mc->Draw("HISTSAME");
  if (drawbkg) h_extra->Draw("HISTSAME");
  h_dt->Draw("ESAME");

  TLegend *legComp = new TLegend( 0.653, 0.663, 0.944, 0.870);
  legComp->AddEntry(h_dt, "Data", "lp");
  legComp->AddEntry(h_mc, "MC", "l");
  if (drawbkg) legComp->AddEntry(h_extra, "MC Bkg", "l");
  legComp->SetFillColor(0);
  legComp->SetBorderSize(0);
  legComp->Draw();
  
  TLatex *text = new TLatex();
  text->SetNDC();
  text->SetTextSize(0.04);
  //  float xtex = 0.65;
  //  text->DrawLatex(xtex,0.88,"1 lepton + jets Sample");
  
  fullpad->cd();
  
  respad = new TPad("respad","respad",0,0.8,1,1);
  respad->Draw();
  respad->cd();
  
  //gPad->SetGridy();
  
  TH1F* ratio = (TH1F*) h_dt->Clone("ratio");
  ratio->Divide(h_mc);

  ratio->GetYaxis()->SetTitleOffset(0.3);
  ratio->GetYaxis()->SetTitleSize(0.2);
  ratio->GetYaxis()->SetNdivisions(5);
  ratio->GetYaxis()->SetLabelSize(0.2);
  if (dolog) ratio->GetYaxis()->SetRangeUser(0.5,1.5);
  else ratio->GetYaxis()->SetRangeUser(0.7,1.3);
  ratio->GetYaxis()->SetTitle("Ratio    ");
  ratio->GetXaxis()->SetLabelSize(0);
  ratio->GetXaxis()->SetTitleSize(0);
  ratio->SetMarkerSize(1);
  ratio->SetLineWidth(2);
  ratio->SetLineColor(kBlue);
  ratio->SetMarkerColor(kBlue);
  ratio->SetFillColor(kBlue);
  ratio->SetFillStyle(3002);
  ratio->Draw("E2");
  
  TLine line;
  line.SetLineWidth(2);
  line.DrawLine(h_dt->GetXaxis()->GetXmin(),1,h_dt->GetXaxis()->GetXmax(),1);

}
