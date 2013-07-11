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

//even entries are for SF derivation 
//odd entries are for bkg estimate
const int NSAMPLE = 8;
string selection[NSAMPLE] = 
  { "_1l1or2b_mt150_g5j", "_1l3b", "_1l1or2b_mt120_g4j", "_1l4b", "_2l2b", "_2l3b", "_2l2b", "_2l4b" }; 
string samplename[NSAMPLE] = 
  { "1l + 1 or 2 b (#geq 5j, M_{T}>150 GeV)", 
    "1l + 3 b (#geq 5j, M_{T}>150 GeV)", 
    "1l + 1 or 2 b (#geq 4j, M_{T}>120 GeV)", 
    "1l + 4 b (#geq 4j, M_{T}>120 GeV)", 
    "2l + 2b (#geq 4j)",
    "2l + 3b (#geq 5j)",
    "2l + 2b (#geq 4j)",
    "2l + 4b (#geq 4j)"};
bool addhighmbb[NSAMPLE] = {false, false, false, false, true, true, true, true};

char * histoname_1l = "h_sig_mt_count";
char * histoname_2l = "h_sig_mbb_count";

string histotag = "";
TFile *f_dt;

void plotComparison( TH1F* h_dt , TH1F* h_mc , TH1F *h_extra, char* label, bool dolog, bool drawbkg);
void runSimplePred(char* ttbar_tag = "", bool issyst = false);
//////////////////////////////////////////////////////////////////////////////

//ttbar prediction and uncertainty
float ttpred[NSAMPLE]; 
float ettpred[NSAMPLE];
float ttsmpred[NSAMPLE]; 
float ettsmpred[NSAMPLE];
//rare prediction and uncertainty
float rarepred[NSAMPLE]; 
float erarepred[NSAMPLE];
//total prediction and uncertainty
float pred[NSAMPLE];
float epred[NSAMPLE];

float xs_unc_rare = 1.;
float xs_unc_smtt = 1.;

bool doprintout = true;

void doPredAll(bool doaltttbar = false) {

  f_dt = TFile::Open("SIGoutput/data_histos.root");

  for (int isr=0; isr<NSAMPLE; ++isr) {
    ttpred[isr] = -999.;
    ettpred[isr] = -999.;
    ttsmpred[isr] = -999.;
    ettsmpred[isr] = -999.;
    rarepred[isr] = -999.;
    erarepred[isr] = -999.;
    pred[isr] = -999.;
    epred[isr] = -999.;
  }

  runSimplePred("_lmgtau", false);

  //cout<<"----- SYSTEMATICS -----"<<endl;

    // printf("Prediction \t ");
    // for (int isr=2; isr<NSAMPLE; ++isr) 
    //   printf("\t & $%.1f \\pm %.1f$ %s",
    // 	     pred[isr],epred[isr], 
    // 	     isr==NSAMPLE-1 ? " \\\\ \n":"");
    // cout<<"\\hline"<<endl;
    // cout<<endl;

  float ttpred_default[NSAMPLE];
  std::copy(&ttpred[0], &ttpred[NSAMPLE], ttpred_default);
  float ettpred_default[NSAMPLE]; 
  std::copy(&ettpred[0], &ettpred[NSAMPLE], ettpred_default);

  float pred_default[NSAMPLE];
  std::copy(&pred[0], &pred[NSAMPLE], pred_default);
  float epred_default[NSAMPLE]; 
  std::copy(&epred[0], &epred[NSAMPLE], epred_default);

  //this is the statistical uncertainty
  float aunc_stat[NSAMPLE]; 
  float runc_stat[NSAMPLE]; 
  std::copy(&epred_default[0], &epred_default[NSAMPLE], aunc_stat);
  for (int isr=0; isr<NSAMPLE; ++isr) {
    if (isr%2==0) continue;
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
  for (int isr=0; isr<NSAMPLE; ++isr) {
    if (isr%2==0) continue;
    aunc_xsrare[isr] = fabs(pred_xsrare[isr]-pred_default[isr]);
    //subtractSqr(epred_xsrare[isr], epred_default[isr]);
    runc_xsrare[isr] = aunc_xsrare[isr]/pred_default[isr]*100.;
    epred_default[isr] = addSqr(epred_default[isr], aunc_xsrare[isr]);
  }

  // printf("Rare Syst \t ");
  // for (int isr=2; isr<NSAMPLE; ++isr) 
  //   printf("\t & $%.1f \\pm %.1f$ %s",
  // 	   pred[isr],epred[isr], 
  // 	   isr==NSAMPLE-1 ? " \\\\ \n":"");
  // cout<<"\\hline"<<endl;
  // cout<<endl;

  xs_unc_rare = 1.;


  //cout<<"---------------------- Normalization of small ttbar background -----------------------"<<endl;
  xs_unc_smtt = 1.1;

  runSimplePred("_lmgtau", false);

  float pred_xssmtt[NSAMPLE];
  std::copy(&pred[0], &pred[NSAMPLE], pred_xssmtt);
  float epred_xssmtt[NSAMPLE]; 
  std::copy(&epred[0], &epred[NSAMPLE], epred_xssmtt);

  float aunc_xssmtt[NSAMPLE];
  float runc_xssmtt[NSAMPLE];
  for (int isr=0; isr<NSAMPLE; ++isr) {
    if (isr%2==0) continue;
    aunc_xssmtt[isr] = fabs(pred_xssmtt[isr]-pred_default[isr]);
    //subtractSqr(epred_xssmtt[isr], epred_default[isr]);
    runc_xssmtt[isr] = aunc_xssmtt[isr]/pred_default[isr]*100.;
    epred_default[isr] = addSqr(epred_default[isr], aunc_xssmtt[isr]);
  }

  printf("Subdominant ttbar Syst \t ");
  for (int isr=0; isr<NSAMPLE; ++isr) {
    if (isr%2==0) continue;
    printf("\t & $%.1f \\pm %.1f$ %s",
  	   pred[isr],epred[isr], 
  	   isr==NSAMPLE-1 ? " \\\\ \n":"");
  }
  cout<<"\\hline"<<endl;
  cout<<endl;

  xs_unc_smtt = 1.;

  //cout<<"---------------------- Alternative tt MC -----------------------"<<endl;

  runSimplePred("_lpowheg", true);

  float pred_ttmc[NSAMPLE];
  std::copy(&pred[0], &pred[NSAMPLE], pred_ttmc);
  float epred_ttmc[NSAMPLE]; 
  std::copy(&epred[0], &epred[NSAMPLE], epred_ttmc);

  float aunc_ttmc[NSAMPLE];
  float runc_ttmc[NSAMPLE];
  for (int isr=0; isr<NSAMPLE; ++isr) {
    if (isr%2==0) continue;
    aunc_ttmc[isr] = fabs(pred_ttmc[isr]-pred_default[isr]);
    //    aunc_ttmc[isr] = addSqr(aunc_ttmc[isr], epred_ttmc[isr]);
    runc_ttmc[isr] = aunc_ttmc[isr]/pred_default[isr]*100.;
    //add this uncertainty
    //    epred_default[isr] = addSqr(epred_default[isr], aunc_ttmc[isr]);
  }

  cout<<endl;
  printf("Alternative MC \t ");
  for (int isr=0; isr<NSAMPLE; ++isr) {
    if (isr%2==0) continue;
    printf("\t & $%.2f \\pm %.2f$ %s",
  	   pred[isr],epred[isr], 
  	   isr==NSAMPLE-1 ? " \\\\ \n":"");
  }
  cout<<"\\hline"<<endl;
  cout<<endl;

  //cout<<"---------------------- Alternative tt MC MCatNLO-----------------------"<<endl;

  runSimplePred("_mcatnlo", true);

  float pred_ttmc2[NSAMPLE];
  std::copy(&pred[0], &pred[NSAMPLE], pred_ttmc2);
  float epred_ttmc2[NSAMPLE]; 
  std::copy(&epred[0], &epred[NSAMPLE], epred_ttmc2);

  float aunc_ttmc2[NSAMPLE];
  float runc_ttmc2[NSAMPLE];
  for (int isr=0; isr<NSAMPLE; ++isr) {
    if (isr%2==0) continue;
    aunc_ttmc2[isr] = fabs(pred_ttmc2[isr]-pred_default[isr]);
    //    aunc_ttmc2[isr] = addSqr(aunc_ttmc2[isr], epred_ttmc2[isr]);
    runc_ttmc2[isr] = aunc_ttmc2[isr]/pred_default[isr]*100.;
    //add this uncertainty
    //    epred_default[isr] = addSqr(epred_default[isr], aunc_ttmc2[isr]);
  }

  cout<<endl;
  printf("Alternative MC \t ");
  for (int isr=0; isr<NSAMPLE; ++isr) {
    if (isr%2==0) continue;
    printf("\t & $%.2f \\pm %.2f$ %s",
  	   pred[isr],epred[isr], 
  	   isr==NSAMPLE-1 ? " \\\\ \n":"");
  }
  cout<<"\\hline"<<endl;
  cout<<endl;
  
  //  cout<<"---------------------- Alternative control region range -----------------------"<<endl;
  histotag = "_alt";
  //  addhighmbb[5] = true;
  runSimplePred("_lmgtau");

  float pred_altrange[NSAMPLE];
  std::copy(&pred[0], &pred[NSAMPLE], pred_altrange);
  float epred_altrange[NSAMPLE]; 
  std::copy(&epred[0], &epred[NSAMPLE], epred_altrange);

  float aunc_altrange[NSAMPLE];
  float runc_altrange[NSAMPLE];
  for (int isr=2; isr<NSAMPLE; ++isr) {
    aunc_altrange[isr] = fabs(pred_altrange[isr]-pred_default[isr]);
    runc_altrange[isr] = aunc_altrange[isr]/pred_default[isr]*100.;
    //add this uncertainty
    //epred_default[isr] = addSqr(epred_default[isr], aunc_altrange[isr]);
  }

  //  addhighmbb[5] = false;

  // printf("Alternative control region range \t ");
  // for (int isr=2; isr<NSAMPLE; ++isr) 
  //   printf("\t & $%.1f \\pm %.1f$ %s",
  //  	   pred[isr],epred[isr], 
  //  	   isr==NSAMPLE-1 ? " \\\\ \n":"");
  // cout<<"\\hline"<<endl;
  // cout<<endl;
  
  //  cout<<"---------------------- Top pt weight -----------------------"<<endl;
  histotag = "_topptwgt";
  //histotag = "";
  runSimplePred("_lmgtau");

  float pred_toppt[NSAMPLE];
  std::copy(&pred[0], &pred[NSAMPLE], pred_toppt);
  float epred_toppt[NSAMPLE]; 
  std::copy(&epred[0], &epred[NSAMPLE], epred_toppt);

  float aunc_toppt[NSAMPLE];
  float runc_toppt[NSAMPLE];
  for (int isr=0; isr<NSAMPLE; ++isr) {
    if (isr%2==0) continue;
    aunc_toppt[isr] = fabs(pred_toppt[isr]-pred_default[isr]);
    runc_toppt[isr] = aunc_toppt[isr]/pred_default[isr]*100.;
    //add this uncertainty
    epred_default[isr] = addSqr(epred_default[isr], aunc_toppt[isr]);
  }

  // printf("Top pT weight \t ");
  // for (int isr=2; isr<NSAMPLE; ++isr) 
  //   printf("\t & $%.1f \\pm %.1f$ %s",
  //  	   pred[isr],epred[isr], 
  //  	   isr==NSAMPLE-1 ? " \\\\ \n":"");
  // cout<<"\\hline"<<endl;
  // cout<<endl;

  histotag = "";

  // cout<<endl;
  // cout<<endl;
  // cout<<endl;

  //cout<<"---------------------- Uncertainty from the SFs -----------------------"<<endl;

  float aunc_sf[NSAMPLE];
  float runc_sf[NSAMPLE] = {0., 0.14, 0., 0.071, 0., 0.10, 0., 0.3};
  for (int isr=0; isr<NSAMPLE; ++isr) {
    if (isr%2==0) continue;
    if (isr==5 || isr==7) aunc_sf[isr] = runc_sf[isr]*pred_default[isr];
    else if (isr==1 || isr==3) aunc_sf[isr] = runc_sf[isr]*ttpred_default[isr];
    else aunc_sf[isr] = 0.;
    runc_sf[isr] = aunc_sf[isr]/pred_default[isr]*100.;
    //add this uncertainty
    epred_default[isr] = addSqr(epred_default[isr], aunc_sf[isr]);
  }

  //Absolute systematic uncertainties

  printf("Stat. \t ");
  for (int isr=0; isr<NSAMPLE; ++isr) 
    if (isr%2==1) printf("\t & $%.1f$ %s", aunc_stat[isr], isr==NSAMPLE-1 ? " \\\\ \n":"");
  printf("Alt. tt MC \t ");
  for (int isr=0; isr<NSAMPLE; ++isr) 
    if (isr%2==1) printf("\t & $%.1f \\pm %.1f$ %s", aunc_ttmc[isr], epred_ttmc[isr], isr==NSAMPLE-1 ? " \\\\ \n":"");
  printf("Alt. tt MC (MCatNLO) \t ");
  for (int isr=0; isr<NSAMPLE; ++isr) 
    if (isr%2==1) printf("\t & $%.1f \\pm %.1f$ %s", aunc_ttmc2[isr], epred_ttmc2[isr], isr==NSAMPLE-1 ? " \\\\ \n":"");
  printf("Alt. CR Range \t ");
  for (int isr=0; isr<NSAMPLE; ++isr) 
    if (isr%2==1) printf("\t & $%.1f$ %s", aunc_altrange[isr], isr==NSAMPLE-1 ? " \\\\ \n":"");
  printf("Top pT Modeling \t ");
  for (int isr=0; isr<NSAMPLE; ++isr) 
    if (isr%2==1) printf("\t & $%.1f$ %s", aunc_toppt[isr], isr==NSAMPLE-1 ? " \\\\ \n":"");
  printf("rare bkg. Norm \t ");
  for (int isr=0; isr<NSAMPLE; ++isr) 
    if (isr%2==1) printf("\t & $%.1f$ %s", aunc_xsrare[isr], isr==NSAMPLE-1 ? " \\\\ \n":"");
  printf("Subdominant ttbar Norm \t ");
  for (int isr=0; isr<NSAMPLE; ++isr) 
    if (isr%2==1) printf("\t & $%.1f$ %s", aunc_xssmtt[isr], isr==NSAMPLE-1 ? " \\\\ \n":"");
  printf("SF uncertainty \t ");
  for (int isr=0; isr<NSAMPLE; ++isr) 
    if (isr%2==1) printf("\t & $%.1f$ %s", aunc_sf[isr], isr==NSAMPLE-1 ? " \\\\ \n":"");
  cout<<"\\hline"<<endl;
  printf("Total \t ");
  for (int isr=0; isr<NSAMPLE; ++isr) 
    if (isr%2==1) printf("\t & $%.1f$ %s", epred_default[isr], isr==NSAMPLE-1 ? " \\\\ \n":"");
  cout<<"\\hline"<<endl;
  cout<<endl;

  cout<<endl;
  cout<<endl;
  cout<<endl;

  //Relative systematic uncertainties
  printf("Stat. \t ");
  for (int isr=0; isr<NSAMPLE; ++isr) 
    if (isr%2==1) printf("\t & $%.1f$ %s", runc_stat[isr], isr==NSAMPLE-1 ? " \\\\ \n":"");
  printf("Alt. tt MC \t ");
  for (int isr=0; isr<NSAMPLE; ++isr) 
    if (isr%2==1) printf("\t & $%.1f \\pm %.1f$ %s", 
	   runc_ttmc[isr], epred_ttmc[isr]/pred_default[isr]*100., isr==NSAMPLE-1 ? " \\\\ \n":"");
  printf("Alt. tt MC (MCatNLO) \t ");
  for (int isr=0; isr<NSAMPLE; ++isr) 
    if (isr%2==1) printf("\t & $%.1f \\pm %.1f$ %s", 
	   runc_ttmc2[isr], epred_ttmc2[isr]/pred_default[isr]*100., isr==NSAMPLE-1 ? " \\\\ \n":"");
  printf("Alt. CR Range \t ");
  for (int isr=0; isr<NSAMPLE; ++isr) 
    if (isr%2==1) printf("\t & $%.1f$ %s", runc_altrange[isr], isr==NSAMPLE-1 ? " \\\\ \n":"");
  printf("Top pT Modeling \t ");
  for (int isr=0; isr<NSAMPLE; ++isr) 
    if (isr%2==1) printf("\t & $%.1f$ %s", runc_toppt[isr], isr==NSAMPLE-1 ? " \\\\ \n":"");
  printf("rare bkg. Norm \t ");
  for (int isr=0; isr<NSAMPLE; ++isr) 
    if (isr%2==1) printf("\t & $%.1f$ %s", runc_xsrare[isr], isr==NSAMPLE-1 ? " \\\\ \n":"");
  printf("Subdominant ttbar Norm \t ");
  for (int isr=0; isr<NSAMPLE; ++isr) 
    if (isr%2==1) printf("\t & $%.1f$ %s", runc_xssmtt[isr], isr==NSAMPLE-1 ? " \\\\ \n":"");
  printf("SF uncertainty \t ");
  for (int isr=0; isr<NSAMPLE; ++isr) 
    if (isr%2==1) printf("\t & $%.1f$ %s", runc_sf[isr], isr==NSAMPLE-1 ? " \\\\ \n":"");
  cout<<"\\hline"<<endl;
  printf("Total \t ");
  for (int isr=0; isr<NSAMPLE; ++isr) 
    if (isr%2==1) printf("\t & $%.1f$ %s", 
			 (epred_default[isr]/pred_default[isr])*100., 
			 isr==NSAMPLE-1 ? " \\\\ \n":"");
  cout<<"\\hline"<<endl;
  cout<<endl;

  //sum all the signal regions
  float pred_total = 0.;
  float epred_total = 0.;
  for (int isr=0; isr<NSAMPLE; ++isr) {
    if (isr%2==0) continue;
    pred_total += pred[isr];
    epred_total = addSqr(epred_total, epred[isr]);
  }
  printf("TOTAL: \t & $%.2f \\pm %.2f (%.1f)$ \\\\ \n", pred_total, epred_total, epred_total/pred_total*100.);

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
    "ttdl",
    "ttsl",
    "rare"};

  //legends only needed for final predictions
  const char* legend[MCID] = {
    "t#bar{t}#rightarrow #font[12]{ll}",
    "t#bar{t}#rightarrow #font[12]{l} + jets",
    "Rare"};

  enum sample{TTDL=0, 
	      TTSL,
	      RARE};

  // bins for regions
  float i_ctr = 1;
  float i_sig = 2;

  //Store raw information
  //Open single lepton files
  TFile *f_mc[MCID];
  for (int j=0;j<MCID;++j) {
    if (j<2)
      f_mc[j] = TFile::Open(Form("SIGoutput/%s%s_histos.root",
				  mcsample[j], ttbar_tag));
    else 
      f_mc[j] = TFile::Open(Form("SIGoutput/%s_histos.root",
				  mcsample[j]));
  }

  //indices are for samples
  TH1F *h_dt[NSAMPLE];
  TH1F *h_mc[MCID][NSAMPLE];
  TH1F *h_mc_tot[NSAMPLE];
  TH1F *h_dt_cor[NSAMPLE];//data corrected for rare backgrounds 
  float sf_ttsm[NSAMPLE];
  float esf_ttsm[NSAMPLE];
  for (int isr=0; isr<NSAMPLE; ++isr) {

    //data
    h_dt[isr] = (TH1F*)f_dt->Get(Form("%s%s%s",isr<4 ? histoname_1l : histoname_2l,
				      histotag.c_str(),selection[isr].c_str()));
    if (doverbose) 
      cout<<"Data "<<Form("%s%s%s",isr<4 ? histoname_1l : histoname_2l,
			  histotag.c_str(),selection[isr].c_str())
	  <<" "<<h_dt[isr]<<endl;
    h_dt[isr]->SetName(Form("h_dt%s",selection[isr].c_str()));
    if (issyst)	zeroHistError(h_dt[isr]);

    if (addhighmbb[isr]) {
      float binval  =  h_dt[isr]->GetBinContent(1) + h_dt[isr]->GetBinContent(3);
      float ebinval =  addSqr(h_dt[isr]->GetBinError(1), h_dt[isr]->GetBinError(3));
      h_dt[isr]->SetBinContent(1, binval);
      h_dt[isr]->SetBinError(1, ebinval);
      h_dt[isr]->SetBinContent(3, 0.);
      h_dt[isr]->SetBinError(3, 0.);
    } 

    h_dt_cor[isr] = (TH1F*)h_dt[isr]->Clone();
    h_dt_cor[isr]->SetName(Form("h_dt_corr%s",selection[isr].c_str()));

    //mc
    for (int j=0;j<MCID;++j) {
      h_mc[j][isr] = (TH1F*)f_mc[j]->Get(Form("%s%s%s",isr<4 ? histoname_1l : histoname_2l,
					      histotag.c_str(),selection[isr].c_str()));
      if (doverbose) 
	cout<<"MC "<<mcsample[j]<<" "
	    <<Form("%s%s%s",isr<4 ? histoname_1l : histoname_2l,
		   histotag.c_str(),selection[isr].c_str())
	    <<" "<<h_mc[j][isr]<<endl;
      if (h_mc[j][isr]==0) {
      h_mc[j][isr] = (TH1F*)f_mc[0]->Get(Form("%s%s%s",isr<4 ? histoname_1l : histoname_2l,
					      histotag.c_str(),selection[isr].c_str()));
	h_mc[j][isr]->SetName(Form("h_mc_%s%s", mcsample[j],selection[isr].c_str()));
	zeroHist(h_mc[j][isr]);
      }
	
      if (doverbose) 
	cout<<"MC "<<mcsample[j]<<" "
	    <<Form("%s%s%s",isr<4 ? histoname_1l : histoname_2l,
		   histotag.c_str(),selection[isr].c_str())
	    <<" "<<h_mc[j][isr]<<endl;
      h_mc[j][isr]->SetName(Form("h_mc_%s%s", mcsample[j],selection[isr].c_str()));
	
      if (issyst && j!=TTSL && j!=TTDL)
	zeroHistError(h_mc[j][isr]);

      if (addhighmbb[isr]) {
	float bin  =  h_mc[j][isr]->GetBinContent(1) + h_mc[j][isr]->GetBinContent(3);
	float ebin =  addSqr(h_mc[j][isr]->GetBinError(1), h_mc[j][isr]->GetBinError(3));
	h_mc[j][isr]->SetBinContent(1, bin);
	h_mc[j][isr]->SetBinError(1, ebin);
	h_mc[j][isr]->SetBinContent(3, 0.);
	h_mc[j][isr]->SetBinError(3, 0.);
      } 
      
      if (j==RARE) h_mc[j][isr]->Scale(xs_unc_rare);
      //For 1l scale ttlj, for 2l scale ttll
      if (isr<4 && j!=TTSL) {
	//to calculate the uncertainty on the ttdl background
	if (j==TTDL) h_mc[j][isr]->Scale(xs_unc_smtt);
	h_dt_cor[isr]->Add(h_mc[j][isr], -1.);
      }
      if (isr>=4 && j!=TTDL) h_dt_cor[isr]->Add(h_mc[j][isr], -1.);

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

    //For the subdominant ttbar background, derive SFs from sideband regions
    int isample = isr<4 ? TTSL : TTDL;
    sf_ttsm[isr] = h_dt_cor[isr]->GetBinContent(i_ctr)/h_mc[isample][isr]->GetBinContent(i_ctr);
    esf_ttsm[isr] = addSqr((h_dt_cor[isr]->GetBinError(i_ctr)/h_dt_cor[isr]->GetBinContent(i_ctr)), 
			   (h_mc[isample][isr]->GetBinError(i_ctr)/h_mc[isample][isr]->GetBinContent(i_ctr)))
      *sf_ttsm[isr];
  }


  string printtag_1l = "\\mt";
  string printtag_2l = "\\mbb";

  //Calculate the SFs

  if (doprintout) {
    cout<<"-------------------------------------------------"<<endl;
    cout<<endl;
    cout<<endl;
    cout<<endl;
    
    cout<<"\\hline"<<endl;
    printf("${\\large \\textbf{N}\\Big|_{SR}^{\\textbf{MC}}}$");
    for (int isr=0; isr<NSAMPLE; ++isr) {
      printf("& $%.2f \\pm %.2f$ ",
	     isr<4 ? h_mc[TTSL][isr]->GetBinContent(i_sig) : h_mc[TTDL][isr]->GetBinContent(i_sig), 
	     isr<4 ? h_mc[TTSL][isr]->GetBinError(i_sig) : h_mc[TTDL][isr]->GetBinError(i_sig));
    }
    printf(" \\\\\n");
    printf("${\\large \\textbf{N}\\Big|_{CR}^{\\textbf{MC}}}$");
    for (int isr=0; isr<NSAMPLE; ++isr) {
      printf("& $%.2f \\pm %.2f$ ",
	     isr<4 ? h_mc[TTSL][isr]->GetBinContent(i_ctr) : h_mc[TTDL][isr]->GetBinContent(i_ctr), 
	     isr<4 ? h_mc[TTSL][isr]->GetBinError(i_ctr) : h_mc[TTDL][isr]->GetBinError(i_ctr));
    }
    printf(" \\\\\n");
    cout<<"\\hline"<<endl;
  }

  float rmbb[NSAMPLE];
  float ermbb[NSAMPLE];
 
  if (doprintout) printf("${\\large \\textbf{R}\\Big|^{\\textbf{MC}}}$");
  for (int isr=0; isr<NSAMPLE; ++isr) {
    int isample = isr<4 ? TTSL : TTDL;
    rmbb[isr] = h_mc[isample][isr]->GetBinContent(i_sig)/h_mc[isample][isr]->GetBinContent(i_ctr);
    ermbb[isr] = addSqr( (h_mc[isample][isr]->GetBinError(i_sig)/h_mc[isample][isr]->GetBinContent(i_sig)), 
			 (h_mc[isample][isr]->GetBinError(i_ctr)/h_mc[isample][isr]->GetBinContent(i_ctr)) )*rmbb[isr];
    if (doprintout) printf("& $%.3f \\pm %.3f$ ",rmbb[isr], ermbb[isr]);
  }
  if (doprintout) {
    printf(" \\\\\n");
    cout<<"\\hline"<<endl;
    cout<<"\\hline"<<endl;

    
    printf("${\\large \\textbf{N}\\Big|_{SR}^\\textbf{DATA}}$");
    for (int isr=0; isr<NSAMPLE; ++isr) {
      if (isr%2==0) 
	printf("& $%.2f \\pm %.2f$ ",
	       h_dt_cor[isr]->GetBinContent(i_sig), 
	       h_dt_cor[isr]->GetBinError(i_sig));
      else printf("& - ");
    }
    printf(" \\\\\n");
    //    cout<<"\\hline"<<endl;
    printf("${\\large \\textbf{N}\\Big|_{CR}^\\textbf{DATA}}$");
    for (int isr=0; isr<NSAMPLE; ++isr) {
      printf("& $%.2f \\pm %.2f$ ",
	     h_dt_cor[isr]->GetBinContent(i_ctr), 
	     h_dt_cor[isr]->GetBinError(i_ctr));
    }
    printf(" \\\\\n");
    cout<<"\\hline"<<endl;
    printf("${\\large \\textbf{R}\\Big|^{\\textbf{DATA}}}$  "); 
  }
  float rmbb_dt[NSAMPLE];
  float ermbb_dt[NSAMPLE];
  for (int isr=0; isr<NSAMPLE; ++isr) {
    if (isr%2==0) {
      rmbb_dt[isr] = h_dt_cor[isr]->GetBinContent(i_sig)/h_dt_cor[isr]->GetBinContent(i_ctr);
      ermbb_dt[isr] = addSqr( (h_dt_cor[isr]->GetBinError(i_sig)/h_dt_cor[isr]->GetBinContent(i_sig)), 
			      (h_dt_cor[isr]->GetBinError(i_ctr)/h_dt_cor[isr]->GetBinContent(i_ctr)) )*rmbb_dt[isr];
      if (doprintout) printf("& $%.3f \\pm %.3f$ ",rmbb_dt[isr],ermbb_dt[isr]);
    } else {
      rmbb_dt[isr] = -999.;
      ermbb_dt[isr] = -999.;
      if (doprintout) printf("& - ");
    }
  }
  if (doprintout) {
    printf(" \\\\\n");
    cout<<"\\hline"<<endl;
    cout<<"\\hline"<<endl;

    printf("${\\large \\textbf{SF}\\Big|_{R}}$  "); 
  }
  float sf_rmbb[NSAMPLE];
  float esf_rmbb[NSAMPLE];
  for (int isr=0; isr<NSAMPLE; ++isr) {
    if (isr%2==0) {
      sf_rmbb[isr] = rmbb_dt[isr]/rmbb[isr];
      esf_rmbb[isr] = addSqr( (ermbb_dt[isr]/rmbb_dt[isr]), (ermbb[isr]/rmbb[isr]) )*sf_rmbb[isr]; 
    } else {
      sf_rmbb[isr] = 0.;
      esf_rmbb[isr] = 0.; 
    }
  }
  if (doprintout) {
    for (int isr=0; isr<NSAMPLE; ++isr) 
      if (isr%2==0) printf("& $%.2f \\pm %.2f$ ",sf_rmbb[isr],esf_rmbb[isr]);
      else printf("& - ");
    printf(" \\\\\n");
    cout<<"\\hline"<<endl;
    printf("${\\large \\textbf{SF}\\Big|_{Norm. CR}}$  "); 
    for (int isr=0; isr<NSAMPLE; ++isr) 
      if (isr%2==1 && isr<4) printf("& $%.2f \\pm %.2f$ ",sf_ttsm[isr+4],esf_ttsm[isr+4]);
      else printf("& - ");
    printf(" \\\\\n");
    cout<<"\\hline"<<endl;
    cout<<"\\hline"<<endl;

    
    printf("Pred. ${\\large \\textbf{N}\\Big|_{SR}}$ ");
  }
  for (int isr=0; isr<NSAMPLE; ++isr) {
    if (isr%2==0) {
      ttpred[isr] = -999.;
      ettpred[isr] = -999.;
      if (doprintout) printf("& - ");
    } else {
      ttpred[isr]  = h_dt_cor[isr]->GetBinContent(i_ctr)*rmbb[isr]*sf_rmbb[isr-1];
      ettpred[isr] = sqrt( pow( (h_dt_cor[isr]->GetBinError(i_ctr)/h_dt_cor[isr]->GetBinContent(i_ctr)), 2) + 
			     pow( (ermbb[isr]/rmbb[isr]), 2) + pow( (esf_rmbb[isr-1]/sf_rmbb[isr-1]),2) ) * ttpred[isr];
      if (doprintout) printf("& $%.2f \\pm %.2f$ ",ttpred[isr], ettpred[isr]);    
    }
  }
  if (doprintout) 
    printf(" \\\\\n");

  cout<<endl;
  cout<<endl;
  cout<<endl;

  if (doprintout) {
    cout<<"\\hline"<<endl;
    cout<<"\\hline"<<endl;
    
    printf("Dominant \\ttbar\\ Prediction ");
    for (int isr=0; isr<NSAMPLE; ++isr) {
      if (isr%2==0) continue;
      printf("& $%.2f \\pm %.2f$ ",ttpred[isr], ettpred[isr]);    
    }
    printf(" \\\\\n");

    cout<<"Sub-dominant \\ttbar\\ ";
  }
  for (int isr=0; isr<NSAMPLE; ++isr) {
    if (isr%2==0) continue;
    int isample = isr<4 ? TTDL : TTSL;
    ttsmpred[isr]  = h_mc[isample][isr]->GetBinContent(i_sig);
    ettsmpred[isr] = h_mc[isample][isr]->GetBinError(i_sig);    
    //addSqr( h_mc[TTDL][isr]->GetBinError(i_sig), (xs_unc_ttdl*ttsmpred[isr]) );
    //for the single lepton SRs, scale the dilepton background
    //for the dilepton SRs, have no yields for single lepton ttbar, so this doesn't matter
    if (isr<4) {
      ttsmpred[isr] *= sf_ttsm[isr];
      ettsmpred[isr] = addSqr( (ettsmpred[isr]/ttsmpred[isr]), (esf_ttsm[isr]/sf_ttsm[isr]) );
    } 
    if (doprintout) printf("& $%.2f \\pm %.2f$ ",ttsmpred[isr], ettsmpred[isr]);    
  }
  if (doprintout) {
    printf(" \\\\\n");
    cout<<"Rare ";
  }

  for (int isr=0; isr<NSAMPLE; ++isr) {
    if (isr%2==0) continue;
    rarepred[isr]  = h_mc[RARE][isr]->GetBinContent(i_sig);
    erarepred[isr] = h_mc[RARE][isr]->GetBinError(i_sig);
    //addSqr( h_mc[RARE][isr]->GetBinError(i_sig), (xs_unc_rare*rarepred[isr]) );
    if (doprintout) printf("& $%.2f \\pm %.2f$ ",rarepred[isr], erarepred[isr]);    
  }
    
  if (doprintout) {
    printf(" \\\\\n");
    cout<<"\\hline"<<endl;
    cout<<"Total ";
  }
  for (int isr=0; isr<NSAMPLE; ++isr) {
    if (isr%2==0) continue;
    pred[isr] = ttpred[isr] + ttsmpred[isr] + rarepred[isr];
    epred[isr] = addSqr( ettpred[isr], ettsmpred[isr]);
    epred[isr] = addSqr( epred[isr], erarepred[isr]);
    if (doprintout) printf("& $%.2f \\pm %.2f$ ",pred[isr], epred[isr]);    
  }
  if (doprintout) {
    printf(" \\\\\n"); 
    cout<<endl;
    cout<<endl;
    cout<<endl;
  }


  // if (doprintout) {
  //   cout<<"------------------ SINGLE LEPTON TABLE -------------------------------"<<endl;
  //   cout<<endl;
  //   cout<<endl;
  //   cout<<endl;
    
  //   cout<<"\\hline"<<endl;
  //   printf("${\\large \\textbf{N}\\Big|_{SR}^{\\textbf{MC}}}$");
  //   for (int isr=0; isr<4; ++isr) 
  //     printf("& $%.2f \\pm %.2f$ ", 
  // 	     h_mc[TTSL][isr]->GetBinContent(i_sig), h_mc[TTSL][isr]->GetBinError(i_sig));
  //   printf(" \\\\\n");
  //   printf("${\\large \\textbf{N}\\Big|_{CR}^{\\textbf{MC}}}$");
  //   for (int isr=0; isr<4; ++isr) 
  //     printf("& $%.2f \\pm %.2f$ ", 
  // 	     h_mc[TTSL][isr]->GetBinContent(i_ctr), h_mc[TTSL][isr]->GetBinError(i_ctr));
  //   printf(" \\\\\n");
  //   cout<<"\\hline"<<endl;
  //   if (doprintout) printf("${\\large \\textbf{R}\\Big|^{\\textbf{MC}}}$");
  //   for (int isr=0; isr<4; ++isr) 
  //     printf("& $%.3f \\pm %.3f$ ",rmbb[isr], ermbb[isr]);
  //   printf(" \\\\\n");
  //   cout<<"\\hline"<<endl;
  //   cout<<"\\hline"<<endl;
  //   printf("${\\large \\textbf{N}\\Big|_{SR}^\\textbf{DATA}}$");
  //   for (int isr=0; isr<4; ++isr) {
  //     if (isr%2==0) 
  // 	printf("& $%.2f \\pm %.2f$ ",
  // 	       h_dt_cor[isr]->GetBinContent(i_sig), h_dt_cor[isr]->GetBinError(i_sig));
  //     else printf("& - ");
  //   }
  //   printf(" \\\\\n");
  //   printf("${\\large \\textbf{N}\\Big|_{CR}^\\textbf{DATA}}$");
  //   for (int isr=0; isr<4; ++isr) {
  //     printf("& $%.2f \\pm %.2f$ ",
  // 	     h_dt_cor[isr]->GetBinContent(i_ctr), h_dt_cor[isr]->GetBinError(i_ctr));
  //   }
  //   printf(" \\\\\n");
  //   cout<<"\\hline"<<endl;
  //   printf("${\\large \\textbf{R}\\Big|^{\\textbf{DATA}}}$  "); 
  //   for (int isr=0; isr<4; ++isr) {
  //     if (isr%2==0) printf("& $%.3f \\pm %.3f$ ",rmbb_dt[isr],ermbb_dt[isr]);
  //     else printf("& - ");
  //   }
  //   printf(" \\\\\n");
  //   cout<<"\\hline"<<endl;
  //   cout<<"\\hline"<<endl;
    
  //   printf("${\\large \\textbf{SF}\\Big|_{R}}$  "); 
  //   for (int isr=0; isr<4; ++isr) 
  //     if (isr%2==0) printf("& $%.2f \\pm %.2f$ ",sf_rmbb[isr],esf_rmbb[isr]);
  //     else printf("& - ");
  //   printf(" \\\\\n");
  //   cout<<"\\hline"<<endl;
  //   printf("${\\large \\textbf{SF}\\Big|_{Norm. CR}}$  "); 
  //   for (int isr=0; isr<4; ++isr) 
  //     if (isr%2==1 && isr<4) printf("& $%.2f \\pm %.2f$ ",sf_ttsm[isr+4],esf_ttsm[isr+4]);
  //     else printf("& - ");
  //   printf(" \\\\\n");
  //   cout<<"\\hline"<<endl;
  //   cout<<"\\hline"<<endl;
  //   printf("Pred. ${\\large \\textbf{N}\\Big|_{SR}}$ ");
  //   for (int isr=0; isr<4; ++isr) 
  //   if (isr%2==0) printf("& - ");
  //   else printf("& $%.2f \\pm %.2f$ ",ttpred[isr], ettpred[isr]);    
  //   printf(" \\\\\n");

  //   cout<<endl;
  //   cout<<endl;
  //   cout<<endl;
  // }


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
