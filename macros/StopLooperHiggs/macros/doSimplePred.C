#include <TH1F.h>
#include <TH2F.h>
#include <TLegend.h>
#include <TFile.h>
#include <TMath.h>
#include <TLine.h>
#include <TLatex.h>
#include <TCanvas.h>
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
bool do1l = true;

const int NSAMPLE = 3;
string selection[NSAMPLE] = {Form("_%il2b", do1l?1:2), Form("_%il3b", do1l?1:2), Form("_%il4b", do1l?1:2)};    
string samplename[NSAMPLE] = {Form("%il + 2b", do1l?1:2), Form("%il + 3b", do1l?1:2), Form("%il + 4b", do1l?1:2)};
bool addhighmbb[NSAMPLE] = {false, false, true};
string histotag = "";
TFile *f_dt;

void plotComparison( TH1F* h_dt , TH1F* h_mc , TH1F *h_extra, char* label, bool dolog, bool drawbkg);
void runSimplePred(char* ttbar_tag = "", bool issyst = false);
//////////////////////////////////////////////////////////////////////////////

//ttbar prediction and uncertainty
float ttpred[NSAMPLE]; 
float ettpred[NSAMPLE];
//rare prediction and uncertainty
float rarepred[NSAMPLE]; 
float erarepred[NSAMPLE];
//total prediction and uncertainty
float pred[NSAMPLE];
float epred[NSAMPLE];

float xs_unc_rare = 1.;

bool doprintout = true;

void doSimplePred(bool doaltttbar = false) {

  cout<<"-------------------------------"<<endl;
  if (do1l) cout<<"1-LEPTON BACKGROUND PREDICTIONS"<<endl;
  else cout<<"2-LEPTON BACKGROUND PREDICTIONS"<<endl;
  cout<<"-------------------------------"<<endl;
  
  f_dt = TFile::Open(Form("SIGoutput/data_%s_histos.root", do1l ? "sl" : "dl"));

  if (do1l) {
    for (int isr=0; isr<NSAMPLE; ++isr) 
      addhighmbb[isr] = false;
  }

  for (int isr=0; isr<NSAMPLE; ++isr) {
    ttpred[isr] = -999.;
    ettpred[isr] = -999.;
    rarepred[isr] = -999.;
    erarepred[isr] = -999.;
    pred[isr] = -999.;
    epred[isr] = -999.;
  }

  runSimplePred("_lmgtau", false);

  //cout<<"----- SYSTEMATICS -----"<<endl;

    // printf("Prediction \t ");
    // for (int isr=1; isr<NSAMPLE; ++isr) 
    //   printf("\t & $%.1f \\pm %.1f$ %s",
    // 	     pred[isr],epred[isr], 
    // 	     isr==NSAMPLE-1 ? " \\\\ \n":"");
    // cout<<"\\hline"<<endl;
    // cout<<endl;

  float pred_default[NSAMPLE];
  std::copy(&pred[0], &pred[NSAMPLE], pred_default);
  float epred_default[NSAMPLE]; 
  std::copy(&epred[0], &epred[NSAMPLE], epred_default);

  //this is the statistical uncertainty
  float aunc_stat[NSAMPLE]; 
  float runc_stat[NSAMPLE]; 
  std::copy(&epred_default[0], &epred_default[NSAMPLE], aunc_stat);
  for (int isr=1; isr<NSAMPLE; ++isr) {
    runc_stat[isr] = aunc_stat[isr]/pred_default[isr]*100.;
  }

  doprintout = false;

  //cout<<"---------------------- Normalization of rare -----------------------"<<endl;
  xs_unc_rare = 1.5;

  runSimplePred("_lmgtau", false);

  float pred_xsrare[NSAMPLE];
  std::copy(&pred[0], &pred[NSAMPLE], pred_xsrare);
  float epred_xsrare[NSAMPLE]; 
  std::copy(&epred[0], &epred[NSAMPLE], epred_xsrare);

  float aunc_xsrare[NSAMPLE];
  float runc_xsrare[NSAMPLE];
  for (int isr=1; isr<NSAMPLE; ++isr) {
    aunc_xsrare[isr] = fabs(pred_xsrare[isr]-pred_default[isr]);
    //subtractSqr(epred_xsrare[isr], epred_default[isr]);
    runc_xsrare[isr] = aunc_xsrare[isr]/pred_default[isr]*100.;
    epred_default[isr] = addSqr(epred_default[isr], aunc_xsrare[isr]);
  }

  // printf("Rare Syst \t ");
  // for (int isr=1; isr<NSAMPLE; ++isr) 
  //   printf("\t & $%.1f \\pm %.1f$ %s",
  // 	   pred[isr],epred[isr], 
  // 	   isr==NSAMPLE-1 ? " \\\\ \n":"");
  // cout<<"\\hline"<<endl;
  // cout<<endl;

  xs_unc_rare = 1.;


  //cout<<"---------------------- Alternative tt MC -----------------------"<<endl;

  runSimplePred("_lpowheg", true);

  float pred_ttmc[NSAMPLE];
  std::copy(&pred[0], &pred[NSAMPLE], pred_ttmc);
  float epred_ttmc[NSAMPLE]; 
  std::copy(&epred[0], &epred[NSAMPLE], epred_ttmc);

  float aunc_ttmc[NSAMPLE];
  float runc_ttmc[NSAMPLE];
  for (int isr=1; isr<NSAMPLE; ++isr) {
    aunc_ttmc[isr] = fabs(pred_ttmc[isr]-pred_default[isr]);
    //    aunc_ttmc[isr] = addSqr(aunc_ttmc[isr], epred_ttmc[isr]);
    runc_ttmc[isr] = aunc_ttmc[isr]/pred_default[isr]*100.;
    //add this uncertainty
    epred_default[isr] = addSqr(epred_default[isr], aunc_ttmc[isr]);
  }

  cout<<endl;
  printf("Alternative MC \t ");
  for (int isr=1; isr<NSAMPLE; ++isr) 
    printf("\t & $%.2f \\pm %.2f$ %s",
  	   pred[isr],epred[isr], 
  	   isr==NSAMPLE-1 ? " \\\\ \n":"");
  cout<<"\\hline"<<endl;
  cout<<endl;
  
  //  cout<<"---------------------- Alternative control region range -----------------------"<<endl;
  histotag = "_alt";
  runSimplePred("_lmgtau");

  float pred_altrange[NSAMPLE];
  std::copy(&pred[0], &pred[NSAMPLE], pred_altrange);
  float epred_altrange[NSAMPLE]; 
  std::copy(&epred[0], &epred[NSAMPLE], epred_altrange);

  float aunc_altrange[NSAMPLE];
  float runc_altrange[NSAMPLE];
  for (int isr=1; isr<NSAMPLE; ++isr) {
    aunc_altrange[isr] = fabs(pred_altrange[isr]-pred_default[isr]);
    runc_altrange[isr] = aunc_altrange[isr]/pred_default[isr]*100.;
    //add this uncertainty
    epred_default[isr] = addSqr(epred_default[isr], aunc_altrange[isr]);
  }

  // printf("Alternative control region range \t ");
  // for (int isr=1; isr<NSAMPLE; ++isr) 
  //   printf("\t & $%.1f \\pm %.1f$ %s",
  //  	   pred[isr],epred[isr], 
  //  	   isr==NSAMPLE-1 ? " \\\\ \n":"");
  // cout<<"\\hline"<<endl;
  // cout<<endl;
  
  //  cout<<"---------------------- Top pt weight -----------------------"<<endl;
  histotag = "_topptwgt";
  runSimplePred("_lmgtau");

  float pred_toppt[NSAMPLE];
  std::copy(&pred[0], &pred[NSAMPLE], pred_toppt);
  float epred_toppt[NSAMPLE]; 
  std::copy(&epred[0], &epred[NSAMPLE], epred_toppt);

  float aunc_toppt[NSAMPLE];
  float runc_toppt[NSAMPLE];
  for (int isr=1; isr<NSAMPLE; ++isr) {
    aunc_toppt[isr] = fabs(pred_toppt[isr]-pred_default[isr]);
    runc_toppt[isr] = aunc_toppt[isr]/pred_default[isr]*100.;
    //add this uncertainty
    epred_default[isr] = addSqr(epred_default[isr], aunc_toppt[isr]);
  }

  // printf("Top pT weight \t ");
  // for (int isr=1; isr<NSAMPLE; ++isr) 
  //   printf("\t & $%.1f \\pm %.1f$ %s",
  //  	   pred[isr],epred[isr], 
  //  	   isr==NSAMPLE-1 ? " \\\\ \n":"");
  // cout<<"\\hline"<<endl;
  // cout<<endl;

  histotag = "";

  // cout<<endl;
  // cout<<endl;
  // cout<<endl;

  //Absolute systematic uncertainties

  printf("Stat. \t ");
  for (int isr=1; isr<NSAMPLE; ++isr) 
    printf("\t & $%.1f$ %s", aunc_stat[isr], isr==NSAMPLE-1 ? " \\\\ \n":"");
  printf("Alt. tt MC \t ");
  for (int isr=1; isr<NSAMPLE; ++isr) 
    printf("\t & $%.1f \\pm %.1f$ %s", aunc_ttmc[isr], epred_ttmc[isr], isr==NSAMPLE-1 ? " \\\\ \n":"");
  printf("Alt. CR Range \t ");
  for (int isr=1; isr<NSAMPLE; ++isr) 
    printf("\t & $%.1f$ %s", aunc_altrange[isr], isr==NSAMPLE-1 ? " \\\\ \n":"");
  printf("Top pT Modeling \t ");
  for (int isr=1; isr<NSAMPLE; ++isr) 
    printf("\t & $%.1f$ %s", aunc_toppt[isr], isr==NSAMPLE-1 ? " \\\\ \n":"");
  printf("rare bkg. Norm \t ");
  for (int isr=1; isr<NSAMPLE; ++isr) 
    printf("\t & $%.1f$ %s", aunc_xsrare[isr], isr==NSAMPLE-1 ? " \\\\ \n":"");
  cout<<"\\hline"<<endl;
  printf("Total \t ");
  for (int isr=1; isr<NSAMPLE; ++isr) 
    printf("\t & $%.1f$ %s", epred_default[isr], isr==NSAMPLE-1 ? " \\\\ \n":"");
  cout<<"\\hline"<<endl;
  cout<<endl;

  cout<<endl;
  cout<<endl;
  cout<<endl;

  //Relative systematic uncertainties
  printf("Stat. \t ");
  for (int isr=1; isr<NSAMPLE; ++isr) 
    printf("\t & $%.1f$ %s", runc_stat[isr], isr==NSAMPLE-1 ? " \\\\ \n":"");
  printf("Alt. tt MC \t ");
  for (int isr=1; isr<NSAMPLE; ++isr) 
    printf("\t & $%.1f \\pm %.1f$ %s", runc_ttmc[isr], epred_ttmc[isr]/pred_default[isr]*100., isr==NSAMPLE-1 ? " \\\\ \n":"");
  printf("Alt. CR Range \t ");
  for (int isr=1; isr<NSAMPLE; ++isr) 
    printf("\t & $%.1f$ %s", runc_altrange[isr], isr==NSAMPLE-1 ? " \\\\ \n":"");
  printf("Top pT Modeling \t ");
  for (int isr=1; isr<NSAMPLE; ++isr) 
    printf("\t & $%.1f$ %s", runc_toppt[isr], isr==NSAMPLE-1 ? " \\\\ \n":"");
  printf("rare bkg. Norm \t ");
  for (int isr=1; isr<NSAMPLE; ++isr) 
    printf("\t & $%.1f$ %s", runc_xsrare[isr], isr==NSAMPLE-1 ? " \\\\ \n":"");
  cout<<"\\hline"<<endl;
  printf("Total \t ");
  for (int isr=1; isr<NSAMPLE; ++isr) 
    printf("\t & $%.1f$ %s", 
	   (epred_default[isr]/pred_default[isr])*100., 
	   isr==NSAMPLE-1 ? " \\\\ \n":"");
  cout<<"\\hline"<<endl;
  cout<<endl;


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
  const int MCID = 2; 
  const char* mcsample[MCID]={
    "ttall",
    "rare"};

  //legends only needed for final predictions
  const char* legend[MCID] = {
    "$t\\bar{t}$",
    "Rare"};

  enum sample{TTBAR=0, 
	      RARE};

  // bins for regions
  float i_ctr = 1;
  float i_sig = 2;

  //Store raw information
  //Open single lepton files
  TFile *f_mc[MCID];
  for (int j=0;j<MCID;++j) {
    if (j<1)
      f_mc[j] = TFile::Open(Form("SIGoutput/%s%s_histos.root",
				  mcsample[j], ttbar_tag));
    else 
      f_mc[j] = TFile::Open(Form("SIGoutput/%s_histos.root",
				  mcsample[j]));
  }

  char * histoname = "h_sig_mbb_count";
  if (do1l) histoname = "h_sig_mt_count";

  //indices are for samples
  TH1F *h_dt[NSAMPLE];
  TH1F *h_mc[MCID][NSAMPLE];
  TH1F *h_mc_tot[NSAMPLE];
  TH1F *h_dt_cor[NSAMPLE];//data corrected for rare backgrounds 

  for (int isr=0; isr<NSAMPLE; ++isr) {

    //data
    h_dt[isr] = (TH1F*)f_dt->Get(Form("%s%s%s",histoname,histotag.c_str(),selection[isr].c_str()));
    if (doverbose) 
      cout<<"Data "<<Form("%s%s%s",histoname,histotag.c_str(),selection[isr].c_str())
	  <<" "<<h_dt[isr]<<endl;
    h_dt[isr]->SetName(Form("h_dt%s",selection[isr].c_str()));
    if (issyst)	zeroHistError(h_dt[isr]);
    h_dt_cor[isr] = (TH1F*)h_dt[isr]->Clone();
    h_dt_cor[isr]->SetName(Form("h_dt_corr%s",selection[isr].c_str()));

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
      h_mc[j][isr] = (TH1F*)f_mc[j]->Get(Form("%s%s%s",histoname,histotag.c_str(),selection[isr].c_str()));
      if (doverbose) 
	cout<<"MC "<<mcsample[j]<<" "
	    <<Form("%s%s%s",histoname,histotag.c_str(),selection[isr].c_str())
	    <<" "<<h_mc[j][isr]<<endl;
      if (h_mc[j][isr]==0) {
      h_mc[j][isr] = (TH1F*)f_mc[0]->Get(Form("%s%s%s",histoname,histotag.c_str(),selection[isr].c_str()));
	h_mc[j][isr]->SetName(Form("h_mc_%s%s", mcsample[j],selection[isr].c_str()));
	zeroHist(h_mc[j][isr]);
      }
	
      if (doverbose) 
	cout<<"MC "<<mcsample[j]<<" "
	    <<Form("%s%s%s",histoname,histotag.c_str(),selection[isr].c_str())
	    <<" "<<h_mc[j][isr]<<endl;
      h_mc[j][isr]->SetName(Form("h_mc_%s%s", mcsample[j],selection[isr].c_str()));
	
      if (issyst && j!=TTBAR)
	zeroHistError(h_mc[j][isr]);

      if (addhighmbb[isr]) {
	float bin  =  h_mc[j][isr]->GetBinContent(1) + h_mc[j][isr]->GetBinContent(3);
	float ebin =  addSqr(h_mc[j][isr]->GetBinError(1), h_mc[j][isr]->GetBinError(3));
	h_mc[j][isr]->SetBinContent(1, bin);
	h_mc[j][isr]->SetBinError(1, ebin);
	h_mc[j][isr]->SetBinContent(3, 0.);
	h_mc[j][isr]->SetBinError(3, 0.);
      } 
      
      if (j!=TTBAR) {
	h_mc[j][isr]->Scale(xs_unc_rare);
	h_dt_cor[isr]->Add(h_mc[j][isr], -1.);
      }

      if (j==0) {
	h_mc_tot[isr] = (TH1F*)h_mc[j][isr]->Clone();
	h_mc_tot[isr]->SetName(Form("h_mc_tot%s",selection[isr].c_str()));
      } else 
	h_mc_tot[isr]->Add(h_mc[j][isr]);
	
      if (doverbose) {
	cout<<"-------------------------------------------------"<<endl;
	printf("%s MC Yields for selection %s \n", mcsample[j], selection[isr].c_str());
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


  string printtag = do1l ? "{{\\rm M}_{\\rm T}}" : "{{\\rm M}_{{\\rm b}\\bar{{\\rm b}}}}";

  //Calculate the SFs
  if (doprintout) {
    cout<<"-------------------------------------------------"<<endl;
    cout<<"\\hline"<<endl;  
    cout<<"\\hline"<<endl;  
    cout<<"Sample  ";
    for (int isr=0; isr<NSAMPLE; ++isr) 
      printf("& %s ", samplename[isr].c_str());
    printf(" \\\\\n");
    cout<<"\\hline"<<endl;
    printf(Form("Data$\_{\\rm %s}$", do1l ? "peak": "out"));
    for (int isr=0; isr<NSAMPLE; ++isr) {
      printf("& $%.2f \\pm %.2f$ ",
	     h_dt_cor[isr]->GetBinContent(i_ctr), 
	     h_dt_cor[isr]->GetBinError(i_ctr));
    }
    printf(" \\\\\n");
    cout<<"\\hline"<<endl;
    printf(Form("MC$\_{\\rm %s}$", do1l ? "peak": "out"));
    for (int isr=0; isr<NSAMPLE; ++isr) {
      printf("& $%.2f \\pm %.2f$ ",
	     h_mc[TTBAR][isr]->GetBinContent(i_ctr), 
	     h_mc[TTBAR][isr]->GetBinError(i_ctr));
    }
    printf(" \\\\\n");
    cout<<"\\hline"<<endl;
    cout<<"\\hline"<<endl;
    
    TH1F *h_sfout[NSAMPLE];
    printf(Form("SF$\_{\\rm %s}$", do1l ? "peak": "out"));
    for (int isr=0; isr<NSAMPLE; ++isr) {
      h_sfout[isr] = (TH1F*)h_dt_cor[isr]->Clone();
      h_sfout[isr]->SetName(Form("h_sfout%s",selection[isr].c_str()));
      h_sfout[isr]->Divide(h_mc[TTBAR][isr]);
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
    printf(Form("MC$\_{\\rm %s}$", do1l ? "tail": "in"));
    for (int isr=0; isr<NSAMPLE; ++isr) {
      printf("& $%.2f \\pm %.2f$ ",
	     h_mc[TTBAR][isr]->GetBinContent(i_sig), 
	     h_mc[TTBAR][isr]->GetBinError(i_sig));
    }
    printf(" \\\\\n");
    printf(Form("MC$\_{\\rm %s}$", do1l ? "peak": "out"));
    for (int isr=0; isr<NSAMPLE; ++isr) {
      printf("& $%.2f \\pm %.2f$ ",
	     h_mc[TTBAR][isr]->GetBinContent(i_ctr), 
	     h_mc[TTBAR][isr]->GetBinError(i_ctr));
    }
    printf(" \\\\\n");
    cout<<"\\hline"<<endl;
  }

  float rmbb[NSAMPLE];
  float ermbb[NSAMPLE];
  if (doprintout) printf(Form("${\\rm R}_%s^{\\rm MC}$  ", printtag.c_str()));
  for (int isr=0; isr<NSAMPLE; ++isr) {
    rmbb[isr] = h_mc[TTBAR][isr]->GetBinContent(i_sig)/h_mc[TTBAR][isr]->GetBinContent(i_ctr);
    ermbb[isr] = addSqr( (h_mc[TTBAR][isr]->GetBinError(i_sig)/h_mc[TTBAR][isr]->GetBinContent(i_sig)), 
			 (h_mc[TTBAR][isr]->GetBinError(i_ctr)/h_mc[TTBAR][isr]->GetBinContent(i_ctr)) )*rmbb[isr];
    if (doprintout) printf("& $%.3f \\pm %.3f$ ",rmbb[isr], ermbb[isr]);
  }
  if (doprintout) {
    printf(" \\\\\n");
    cout<<"\\hline"<<endl;
    cout<<"\\hline"<<endl;

    printf(Form("Data$\_{\\rm %s}$ & $%.2f \\pm %.2f$ & - & - \\\\\n", 
		do1l ? "tail" : "in",
		h_dt_cor[0]->GetBinContent(i_sig), 
		h_dt_cor[0]->GetBinError(i_sig)));
    printf(Form("Data$\_{\\rm %s}$", do1l ? "peak": "out"));
    for (int isr=0; isr<NSAMPLE; ++isr) {
      printf("& $%.2f \\pm %.2f$ ",
	     h_dt_cor[isr]->GetBinContent(i_ctr), 
	     h_dt_cor[isr]->GetBinError(i_ctr));
    }
    printf(" \\\\\n");
    cout<<"\\hline"<<endl;
    printf(Form("${\\rm R}_%s^{\\rm Data}$  ", printtag.c_str()));
  }
  float rmbb_dt = h_dt_cor[0]->GetBinContent(i_sig)/h_dt_cor[0]->GetBinContent(i_ctr);
  float ermbb_dt = addSqr( (h_dt_cor[0]->GetBinError(i_sig)/h_dt_cor[0]->GetBinContent(i_sig)), 
			   (h_dt_cor[0]->GetBinError(i_ctr)/h_dt_cor[0]->GetBinContent(i_ctr)) )*rmbb_dt;
  if (doprintout) {
    printf("& $%.3f \\pm %.3f$ & - & - ",rmbb_dt, ermbb_dt);
    printf(" \\\\\n");
    cout<<"\\hline"<<endl;
    cout<<"\\hline"<<endl;
    
    printf(Form("SF(${\\rm R}_%s$)  ", printtag.c_str()));
  }
  float sf_rmbb = rmbb_dt/rmbb[0];
  float esf_rmbb = addSqr( (ermbb_dt/rmbb_dt), (ermbb[0]/rmbb[0]) )*sf_rmbb; 
  if (doprintout) {
    printf("& $%.2f \\pm %.2f$ & - & - ",sf_rmbb, esf_rmbb);
    printf(" \\\\\n");
    cout<<"\\hline"<<endl;
    cout<<"\\hline"<<endl;
    
    printf(Form("Top pred. = Data$\_{\\rm %s}$ * ${\\rm R}\_%s^{\\rm MC}$ * SF(${\\rm R}\_%s$, 2b) & -  ", 
		do1l ? "peak" : "out", printtag.c_str(), printtag.c_str()));
  }
  for (int isr=1; isr<NSAMPLE; ++isr) {
    ttpred[isr]  = h_dt_cor[isr]->GetBinContent(i_ctr)*rmbb[isr]*sf_rmbb;
    ettpred[isr] = sqrt( pow( (h_dt_cor[isr]->GetBinError(i_ctr)/h_dt_cor[isr]->GetBinContent(i_ctr)), 2) + 
			 pow( (ermbb[isr]/rmbb[isr]), 2) + pow( (esf_rmbb/sf_rmbb),2) ) * ttpred[isr];
    if (doprintout) printf("& $%.2f \\pm %.2f$ ",ttpred[isr], ettpred[isr]);    
  }
  if (doprintout) {
    printf(" \\\\\n");
    cout<<"Rare & - ";
  }
  for (int isr=1; isr<NSAMPLE; ++isr) {
    rarepred[isr]  = h_mc[RARE][isr]->GetBinContent(i_sig);
    erarepred[isr] = h_mc[RARE][isr]->GetBinError(i_sig);//addSqr( h_mc[RARE][isr]->GetBinError(i_sig), (xs_unc_rare*rarepred[isr]) );
    if (doprintout) printf("& $%.2f \\pm %.2f$ ",rarepred[isr], erarepred[isr]);    
  }
  if (doprintout) {
    printf(" \\\\\n");
    
    cout<<"\\hline"<<endl;
    cout<<"Total Pred. & -  ";
  }
  for (int isr=1; isr<NSAMPLE; ++isr) {
    pred[isr] = ttpred[isr] + rarepred[isr];
    epred[isr] = addSqr( ettpred[isr], erarepred[isr]);
    if (doprintout) printf("& $%.2f \\pm %.2f$ ",pred[isr], epred[isr]);    
  }
  if (doprintout) {
    printf(" \\\\\n"); 
    cout<<endl;
    cout<<endl;
    cout<<endl;
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
