//////////////////////////////////////////////////////////////////////////////
// Assumptions:
//
// - Data and MC files are in two directories output_sl and output_dl  
//   corresponding to the single lepton and dilepton selections. The
//   samples to be added are listed at the top of the file
// - Canvases are saved to a directory named with the variable outdirname
//
// Usage:
//
// - compile the macro in root:
//  [] .L doFullPred.C++
// - dump predictions and uncertainties for central values:
//  [] doFullPred()
// - for alternative samples, pass tag for sample, for example
//  [] doFullPred("_powheg")
//
//////////////////////////////////////////////////////////////////////////////
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

void plotComparison( TH1F* h_dt , TH1F* h_mc , TH1F *h_extra, char* label, bool dolog, bool drawbkg);

//////////////////////////////////////////////////////////////////////////////

void doPrintYields(char* ttbar_tag = "") {

  // this flag sets all uncertainties except for those
  // associated with the ttbar sample to zero
  // in this way it calculates the uncorrelated uncertainty
  //  bool issyst = true;
  bool issyst = (strlen(ttbar_tag) == 0) ? false : true;

  // to dump all the information about the yields for 
  // the various MC samples
  bool doverbose = true;

  //for single lepton, split single top
  const int MCID = 9; 
  const char* mcsample[MCID]={
    "ttdl_powheg",
    "ttsl_powheg",
    "w1to4jets",
    "tW_lepsl",
    "tW_lepdl",
    "ttV", 
    "diboson",
    "DY1to4Jtot", 
    "triboson"};

  //legends only needed for final predictions
  const char* legend[MCID] = {
    "$t\\bar{t} \\rightarrow l^{+}l^{-}$",
    "$t\\bar{t} \\rightarrow l^{\\pm}$ + jets",
    "W+jets",
    "single top (s/t-chan, 1-lep)",
    "single top (tW, 2-lep)",
    "ttW/Z/$\\gamma$", 
    "WW/WZ/ZZ",
    "DY+jets", 
    "triboson"};
  //other labels that might be useful in the future
    // "t#bar{t} #rightarrow #font[12]{l^{+}l^{-}} (e/#mu)",
    // "t#bar{t} #rightarrow #font[12]{l^{+}l^{-}} (lost)",
    // "t#bar{t} #rightarrow #font[12]{l^{+}l^{-}} (#tau_{had}#rightarrow1-prong)",
    // "t#bar{t} #rightarrow #font[12]{l^{+}l^{-}} (#tau_{had}#rightarrow3-prong)",
    // "t#bar{t} #rightarrow #font[12]{l^{+}l^{-}} (#tau_{lep})"
  //  "T2tt m(stop) = 250 m($\chi^0$) = 50",
  //  "T2tt m(stop) = 350 m($\chi^0$) = 50",
  //  "T2tt m(stop) = 400 m($\chi^0$) = 50",
  //  "T2tt m(stop) = 350 m($\chi^0$) = 100",
  //  "T2tt m(stop) = 400 m($\chi^0$) = 100",

  enum sample{TTDL=0, 
	      TTSL, 
	      WJETS, 
	      TWSL, 
	      TWDL, 
	      TTV, 
	      DIBO, 
	      DY, 
	      VVV};


  //-------------------------------
  // DILEPTON - N JET SCALING
  //-------------------------------

  //list of samples
  const int MCIDDL = 7; 
  const char* mcsampledl[MCIDDL]={
    "ttdl_powheg",
    "ttsl_powheg",
    //    "wstitchjets",
    "tWall_lep",
    "ttV", 
    "diboson",
    "DY1to4Jtot", 
    "triboson"};

  // the input n jets histogram contains
  // 4 bins for multiplicities 1, 2, 3 and >=4

  TFile *dt_dl[3];
  dt_dl[0] = TFile::Open("NJoutput/data_mueg_histos.root");
  dt_dl[1] = TFile::Open("NJoutput/data_dimu_histos.root");
  dt_dl[2] = TFile::Open("NJoutput/data_diel_histos.root");
  
  TH1F *h_dt_dl = (TH1F*)dt_dl[0]->Get("h_njets_met100_mueg");
  h_dt_dl->SetName("h_dt_dl");
  TH1F *h_dt_dl_tmp = (TH1F*)dt_dl[1]->Get("h_njets_met100_dimu");
  h_dt_dl->Add(h_dt_dl_tmp);
  h_dt_dl_tmp = (TH1F*)dt_dl[2]->Get("h_njets_met100_diel");
  h_dt_dl->Add(h_dt_dl_tmp);

  // Load histograms
  
  TFile *mc_dl[MCIDDL];
  TH1F *h_mc_dl[MCIDDL];
  TH1F *h_mc_tot_dl;
  TH1F *h_bg_dl;

  for (int j=0;j<MCIDDL;++j) {
    if (j<2)
      mc_dl[j] = TFile::Open(Form("NJoutput/%s%s_histos.root",
     				  mcsampledl[j], ttbar_tag));
    else 
      mc_dl[j] = TFile::Open(Form("NJoutput/%s_histos.root",
				  mcsampledl[j]));
    h_mc_dl[j] = (TH1F*)mc_dl[j]->Get("h_njets_met100");
    h_mc_dl[j]->SetName(Form("h_mc_dl_%s", mcsampledl[j]));
    //do not include dilepton in background
    if (j==TTDL) {
      h_mc_tot_dl = (TH1F*)h_mc_dl[j]->Clone("h_mc_tot_dl");
    } else {
      h_mc_tot_dl->Add(h_mc_dl[j]); 
      if (j==TTSL) {
	h_bg_dl = (TH1F*)h_mc_dl[j]->Clone("h_bg_dl");
      } else if (j>TTSL) 
	h_bg_dl->Add(h_mc_dl[j]);
    }
  }
  
  // normalize MC to data
  float norm = h_dt_dl->Integral()*1./h_mc_tot_dl->Integral();
  h_bg_dl->Scale(norm);
  h_mc_dl[TTDL]->Scale(norm);

  cout<<"-------------------------------------------------"<<endl;
  cout<<"*************************************************"<<endl;
  cout<<"-------------------------------------------------"<<endl;
  cout<<"2-LEPTON SAMPLE JETS SCALING"<<endl;
  cout<<"-------------------------------------------------"<<endl;

  int i_2j = 1;
  int i_max_2j = 2;
  float n_dt_2j = h_dt_dl->Integral(i_2j, i_max_2j);
  float e_n_dt_2j = histError(h_dt_dl, i_2j, i_max_2j);
  float n_bg_2j = h_bg_dl->Integral(i_2j, i_max_2j);
  float e_n_bg_2j = histError(h_bg_dl, i_2j, i_max_2j);
  float n_sg_2j = h_mc_dl[TTDL]->Integral(i_2j, i_max_2j);
  float e_n_sg_2j = histError(h_mc_dl[TTDL], i_2j, i_max_2j);
  cout<<"--- <=2 Jet Bin ---"<<endl;
  printf("Data: %.1f pm %.1f \n", n_dt_2j, e_n_dt_2j);
  printf("MC non-dil: %.1f pm %.1f \n", n_bg_2j, e_n_bg_2j);
  printf("MC dil: %.1f pm %.1f \n", n_sg_2j, e_n_sg_2j);

  int i_3j = 3;
  float n_dt_3j = h_dt_dl->GetBinContent(i_3j);
  float e_n_dt_3j = h_dt_dl->GetBinError(i_3j);
  float n_bg_3j = h_bg_dl->GetBinContent(i_3j);
  float e_n_bg_3j = h_bg_dl->GetBinError(i_3j);
  float n_sg_3j = h_mc_dl[TTDL]->GetBinContent(i_3j);
  float e_n_sg_3j = h_mc_dl[TTDL]->GetBinError(i_3j);
  cout<<"--- 3 Jet Bin ---"<<endl;
  printf("Data: %.1f pm %.1f \n", n_dt_3j, e_n_dt_3j);
  printf("MC non-dil: %.1f pm %.1f \n", n_bg_3j, e_n_bg_3j);
  printf("MC dil: %.1f pm %.1f \n", n_sg_3j, e_n_sg_3j);

  int i_4j = 4;
  float n_dt_4j = h_dt_dl->GetBinContent(i_4j);
  float e_n_dt_4j = h_dt_dl->GetBinError(i_4j);
  float n_bg_4j = h_bg_dl->GetBinContent(i_4j);
  float e_n_bg_4j = h_bg_dl->GetBinError(i_4j);
  float n_sg_4j = h_mc_dl[TTDL]->GetBinContent(i_4j);
  float e_n_sg_4j = h_mc_dl[TTDL]->GetBinError(i_4j);
  cout<<"--- >=4 Jet Bin ---"<<endl;
  printf("DATA: %.1f pm %.1f \n", n_dt_4j, e_n_dt_4j);
  printf("MC non-dil: %.1f pm %.1f \n", n_bg_4j, e_n_bg_4j);
  printf("MC dil: %.1f pm %.1f \n", n_sg_4j, e_n_sg_4j);


  cout<<"-------------------------------------------------"<<endl;
  cout<<"*************************************************"<<endl;
  cout<<"-------------------------------------------------"<<endl;

  // Calculate K Factors
  // these can be defined in terms of SFs in the <=2j, 3j and >=4j ranges
  // K3 = SF2*SF3 = M2/N2 * N3/M3 and K4 = SF2*SF4 = M2/N2 * N4/M4
  cout<<"MCSF"<<endl;
  //for N<=2, define the SF as MC/Data
  float mcsf_2j = n_sg_2j/(n_dt_2j-n_bg_2j);
  float den_2j = n_dt_2j-n_bg_2j;
  float e_den_2j = sqrt( pow(e_n_dt_2j,2) + pow(e_n_bg_2j, 2) );
  float e_mcsf_2j = sqrt( pow( (e_den_2j/den_2j), 2) 
			  + pow( (e_n_sg_2j/n_sg_2j), 2 )) * mcsf_2j;
  printf("1/SF2 = M2/D2: %.2f pm %.2f \n", mcsf_2j, e_mcsf_2j);

  //for rest, define as Data/MC
  float mcsf_3j = (n_dt_3j-n_bg_3j)/n_sg_3j;
  float num_3j = n_dt_3j-n_bg_3j;
  float e_num_3j = sqrt( pow(e_n_dt_3j,2) + pow(e_n_bg_3j, 2) );
  float e_mcsf_3j = sqrt( pow( (e_num_3j/num_3j), 2) 
			  + pow( (e_n_sg_3j/n_sg_3j), 2 )) * mcsf_3j;
  printf("SF3 = D3/M3: %.2f pm %.2f \n", mcsf_3j, e_mcsf_3j);

  float mcsf_4j = (n_dt_4j-n_bg_4j)/n_sg_4j;
  float num_4j = n_dt_4j-n_bg_4j;
  float e_num_4j = sqrt( pow(e_n_dt_4j,2) + pow(e_n_bg_4j, 2) );
  float e_mcsf_4j = sqrt( pow( (e_num_4j/num_4j), 2) 
			  + pow( (e_n_sg_4j/n_sg_4j), 2 )) * mcsf_4j;
  printf("SF4 = D4/M4: %.2f pm %.2f \n", mcsf_4j, e_mcsf_4j);

  //Calculate K factors
  cout<<"-------------------------------------------------"<<endl;
  cout<<"K FACTORS"<<endl;
  float k3 = mcsf_2j * mcsf_3j;
  float e_k3 = sqrt( pow( (e_mcsf_2j/mcsf_2j), 2) 
		     + pow( (e_mcsf_3j/mcsf_3j), 2 )) * k3;
  printf("K3 = SF3/SF2: %.2f pm %.2f \n", k3, e_k3);

  float k4 = mcsf_2j * mcsf_4j;
  float e_k4 = sqrt( pow( (e_mcsf_2j/mcsf_2j), 2) 
		     + pow( (e_mcsf_4j/mcsf_4j), 2 )) * k4;
  printf("K4 = SF4/SF2: %.2f pm %.2f \n", k4, e_k4);
  cout<<"-------------------------------------------------"<<endl;
  cout<<"*************************************************"<<endl;
  cout<<"-------------------------------------------------"<<endl;

  //-------------------------------
  // SINGLE LEPTON SAMPLE
  //-------------------------------

  //settings go here
  const int NSAMPLE = 4;
  char * metcut[NSAMPLE]  = {"_met150", "_met200", "_met250", "_met300"};

  // const int NSAMPLE = 7;
  // char * metcut[NSAMPLE]  = {"_met100", "_met150", "_met200", "_met250", "_met300", "_met350", "_met400"};

  // bins for regions
  float i_ctr = 1;
  float i_sig = 2;

  float ibin = i_sig;

  bool dokfscaling = false;
  bool dopreveto = false;

  //Store raw information for electron and muon
  enum leptontype{MUO = 0, ELE, COMB};
  char * leptype[2] = {"_muo", "_ele"};

  //Open single lepton files
  TFile *dt_sl[2];
  dt_sl[MUO] = TFile::Open("SIGoutput/data_muo_histos.root");
  dt_sl[ELE] = TFile::Open("SIGoutput/data_ele_histos.root");
  TFile *mc_sl[MCID];
  for (int j=0;j<MCID;++j) {
    if (j<2)
      mc_sl[j] = TFile::Open(Form("SIGoutput/%s%s_histos.root",
				  mcsample[j], ttbar_tag));
    else 
      mc_sl[j] = TFile::Open(Form("SIGoutput/%s_histos.root",
				  mcsample[j]));
  }

  char * histoname = "h_sig_mt_count";

  //-------------------------------
  // MEAT OF CALCULATION
  //-------------------------------
  cout<<"SINGLE LEPTON SAMPLE PREDICTION DERIVATION"<<endl;
  cout<<"-------------------------------------------------"<<endl;
  cout<<"MC CONTENT DUMP"<<endl;
  cout<<"Control = mT in range 60-100 GeV"<<endl;
  cout<<"Signal = mT>150 GeV"<<endl;

  //indices are for muon, electron and sum
  TH1F *h_dt[3][NSAMPLE];
  TH1F *h_mc[MCID][3][NSAMPLE];
  TH1F *h_mc_sl_tot[3][NSAMPLE];
  TH1F *h_mc_ot_tot[3][NSAMPLE];
  TH1F *h_mc_tot[3][NSAMPLE];

  for (int isr=0; isr<NSAMPLE; ++isr) {

    for (int itag=0; itag<2; ++itag) {
      h_dt[itag][isr]= (TH1F*)dt_sl[itag]->Get(Form("%s%s%s",histoname,metcut[isr],leptype[itag]));
      h_dt[itag][isr]->SetName(Form("h_dt%s%s",metcut[isr],leptype[itag]));
      if (dopreveto) {
      	TH1F *h_dt_tmp = (TH1F*)dt_sl[itag]->Get(Form("%s_wisotrk%s%s",histoname,metcut[isr],leptype[itag]));
      	h_dt_tmp->SetName("h_dt_tmp");
      	h_dt_tmp->Sumw2();
      	h_dt[itag][isr]->Add(h_dt_tmp);
      }
      if (issyst)	zeroHistError(h_dt[itag][isr]);

    for (int j=0;j<MCID;++j) {
      h_mc[j][itag][isr] = (TH1F*)mc_sl[j]->Get(Form("%s%s%s",
						     histoname,metcut[isr],leptype[itag]));
      if (h_mc[j][itag][isr]==0) {
	h_mc[j][itag][isr] = (TH1F*)mc_sl[0]->Get(Form("%s%s%s",
						       histoname,metcut[isr],leptype[itag]));
	h_mc[j][itag][isr]->SetName(Form("h_mc_%s%s%s", 
					 mcsample[j],metcut[isr],leptype[itag]));
	zeroHist(h_mc[j][itag][isr]);
      }
      h_mc[j][itag][isr]->Sumw2();
      if (dopreveto) {
      	TH1F *h_mc_tmp = (TH1F*)mc_sl[j]->Get(Form("%s_wisotrk%s%s",
      						     histoname,metcut[isr],leptype[itag]));      
      	if (h_mc_tmp!=0) h_mc[j][itag][isr]->Add(h_mc_tmp);
      }
      if (dokfscaling && j==TTDL) {
      	h_mc[j][itag][isr] = (TH1F*)mc_sl[j]->Get(Form("%s%s_K3%s",
      						       histoname,metcut[isr],leptype[itag]));
      	h_mc[j][itag][isr]->SetName(Form("h_mc_%s%s%s", 
      					 mcsample[j],metcut[isr],leptype[itag]));
      	h_mc[j][itag][isr]->Sumw2();
      	h_mc[j][itag][isr]->Scale(k3);
      	if (dopreveto) {
      	  TH1F *h_mc_tmp = (TH1F*)mc_sl[j]->Get(Form("%s_wisotrk%s_K3%s",
      						     histoname,metcut[isr],leptype[itag]));      
      	  if (h_mc_tmp!=0) {
      	    h_mc_tmp->Sumw2();
      	    h_mc_tmp->Scale(k3);
      	    h_mc[j][itag][isr]->Add(h_mc_tmp);
      	  }
      	}
      	TH1F *h_tmp = (TH1F*)mc_sl[j]->Get(Form("%s%s_K4%s",
      						histoname,metcut[isr],leptype[itag]));
      	h_tmp->SetName("h_tmp");
      	h_tmp->Sumw2();
      	h_tmp->Scale(k4);
      	if (dopreveto) {
      	  TH1F *h_mc_tmp = (TH1F*)mc_sl[j]->Get(Form("%s_wisotrk%s_K4%s",
      						     histoname,metcut[isr],leptype[itag]));      
      	  if (h_mc_tmp!=0) {
      	    h_mc_tmp->Scale(k4);
      	    h_tmp->Add(h_mc_tmp);
      	  }
      	}
      	h_mc[j][itag][isr]->Add(h_tmp);	
      }

      if (doverbose) 
	cout<<"MC "<<mcsample[j]<<" "
	    <<Form("%s%s%s",histoname,metcut[isr],leptype[itag])
	    <<" "<<h_mc[j][itag][isr]<<endl;
      h_mc[j][itag][isr]->SetName(Form("h_mc_%s%s%s", 
				       mcsample[j],metcut[isr],leptype[itag]));
      
      if (issyst && j>1)
	zeroHistError(h_mc[j][itag][isr]);

      if (j==TTDL) {
	h_mc_tot[itag][isr] = (TH1F*)h_mc[j][itag][isr]->Clone();
	h_mc_tot[itag][isr]->SetName(Form("h_mc_tot%s%s",metcut[isr],leptype[itag]));
      } else 
	h_mc_tot[itag][isr]->Add(h_mc[j][itag][isr]);

      //default combined single lepton top
      if (j==TTSL) {
	h_mc_sl_tot[itag][isr] = (TH1F*)h_mc[j][itag][isr]->Clone("h_mc_sl_tot");
	h_mc_sl_tot[itag][isr]->SetName(Form("h_mc_sl_tot%s%s",metcut[isr],leptype[itag]));
	cout<<"sample "<<j<<" should be TTSL content"<<h_mc[j][itag][isr]->GetBinContent(ibin)<<endl;
      } else if (j==TWSL) {
	h_mc_sl_tot[itag][isr]->Add(h_mc[j][itag][isr]);
	cout<<"sample "<<j<<" should be TWSL content"<<h_mc[j][itag][isr]->GetBinContent(ibin)<<endl;
	cout<<"gives total "<<h_mc_sl_tot[itag][isr]->GetBinContent(ibin)<<endl;
      }
      //make histogram for rare samples
      if (j==TWDL) {
	h_mc_ot_tot[itag][isr] = (TH1F*)h_mc[j][itag][isr]->Clone("h_mc_ot_tot");
	h_mc_ot_tot[itag][isr]->SetName(Form("h_mc_ot_tot%s%s",metcut[isr],leptype[itag]));
      } else if (j==TTV || j==DIBO || j==DY || j==VVV)
	h_mc_ot_tot[itag][isr]->Add(h_mc[j][itag][isr]);	
      
      if (doverbose) {
	cout<<"-------------------------------------------------"<<endl;
	printf("%s MC Yields for MET cut %s and %s \n", mcsample[j], metcut[isr],
	       itag==MUO ? "muon" : "electron");
	printf("Control: %.2f pm %.2f (%.1f %%) \n", 
	       h_mc[j][itag][isr]->GetBinContent(i_ctr), 
	       h_mc[j][itag][isr]->GetBinError(i_ctr), 
	       h_mc[j][itag][isr]->GetBinError(i_ctr)/
	       h_mc[j][itag][isr]->GetBinContent(i_ctr)*100.);
	printf("Signal: %.2f pm %.2f (%.1f %%) \n", 
	       h_mc[j][itag][isr]->GetBinContent(ibin), 
	       h_mc[j][itag][isr]->GetBinError(ibin), 
	       h_mc[j][itag][isr]->GetBinError(ibin)/
	       h_mc[j][itag][isr]->GetBinContent(ibin)*100.);
      }
    }
  }

  //Sum of two lepton flavors
    h_dt[COMB][isr]= (TH1F*)h_dt[MUO][isr]->Clone();
    h_dt[COMB][isr]->SetName(Form("h_dt_comb_%s", metcut[isr]));
    h_dt[COMB][isr]->Add(h_dt[ELE][isr]);

  for (int j=0;j<MCID;++j) {
    h_mc[j][COMB][isr]= (TH1F*)h_mc[j][MUO][isr]->Clone();
    h_mc[j][COMB][isr]->SetName(Form("h_mc_comb_%s%s", mcsample[j],metcut[isr]));
    h_mc[j][COMB][isr]->Add(h_mc[j][ELE][isr]);
  
    if (doverbose) {
      cout<<"-------------------------------------------------"<<endl;
      printf("%s COMBINED MC Yields MET cut %s \n", mcsample[j], metcut[isr] );
      printf("Control: %.2f pm %.2f (%.1f %%) \n", 
	     h_mc[j][COMB][isr]->GetBinContent(i_ctr), 
	     h_mc[j][COMB][isr]->GetBinError(i_ctr), 
	     h_mc[j][COMB][isr]->GetBinError(i_ctr)/
	     h_mc[j][COMB][isr]->GetBinContent(i_ctr)*100.);
      printf("Signal: %.2f pm %.2f (%.1f %%) \n", 
	     h_mc[j][COMB][isr]->GetBinContent(ibin), 
	     h_mc[j][COMB][isr]->GetBinError(ibin), 
	     h_mc[j][COMB][isr]->GetBinError(ibin)/
	     h_mc[j][COMB][isr]->GetBinContent(ibin)*100.);
    }

  }

  //Sum of two lepton flavors
  h_mc_sl_tot[COMB][isr]= (TH1F*)h_mc_sl_tot[MUO][isr]->Clone();
  h_mc_sl_tot[COMB][isr]->SetName("h_sl_tot_comb");
  h_mc_sl_tot[COMB][isr]->Add(h_mc_sl_tot[ELE][isr]);
  h_mc_ot_tot[COMB][isr]= (TH1F*)h_mc_ot_tot[MUO][isr]->Clone();
  h_mc_ot_tot[COMB][isr]->SetName("h_ot_tot_comb");
  h_mc_ot_tot[COMB][isr]->Add(h_mc_ot_tot[ELE][isr]);
  h_mc_tot[COMB][isr]= (TH1F*)h_mc_tot[MUO][isr]->Clone();
  h_mc_tot[COMB][isr]->SetName("h_mc_tot_comb");
  h_mc_tot[COMB][isr]->Add(h_mc_tot[ELE][isr]);

  cout<<"-------------------------------------------------"<<endl;
  cout<<"*************************************************"<<endl;
  cout<<"-------------------------------------------------"<<endl;
  cout<<"PRINT DETAILED YIELDS TABLE"<<endl;
  cout<<"\\hline"<<endl;
  cout<<"\\hline"<<endl;
  printf("\\ttdl\\ \t\t & $%.1f \\pm %.1f$ ($%.1f \\%%$) & $%.1f \\pm %.1f$ ($%.1f \\%%$) & $%.1f \\pm %.1f$ ($%.1f \\%%$) \\\\\n",
	 h_mc[TTDL][MUO][isr]->GetBinContent(ibin), 
	 h_mc[TTDL][MUO][isr]->GetBinError(ibin), 
	 h_mc[TTDL][MUO][isr]->GetBinError(ibin)/
	 h_mc[TTDL][MUO][isr]->GetBinContent(ibin)*100., 
	 h_mc[TTDL][ELE][isr]->GetBinContent(ibin), 
	 h_mc[TTDL][ELE][isr]->GetBinError(ibin), 
	 h_mc[TTDL][ELE][isr]->GetBinError(ibin)/
	 h_mc[TTDL][ELE][isr]->GetBinContent(ibin)*100., 
	 h_mc[TTDL][COMB][isr]->GetBinContent(ibin), 
	 h_mc[TTDL][COMB][isr]->GetBinError(ibin), 
	 h_mc[TTDL][COMB][isr]->GetBinError(ibin)/
	 h_mc[TTDL][COMB][isr]->GetBinContent(ibin)*100.);      
  printf("\\ttsl+ single top (1l) \t\t & $%.1f \\pm %.1f$ ($%.1f \\%%$) & $%.1f \\pm %.1f$ ($%.1f \\%%$) & $%.1f \\pm %.1f$ ($%.1f \\%%$) \\\\\n",
	 h_mc_sl_tot[MUO][isr]->GetBinContent(ibin), 
	 h_mc_sl_tot[MUO][isr]->GetBinError(ibin), 
	 h_mc_sl_tot[MUO][isr]->GetBinError(ibin)/
	 h_mc_sl_tot[MUO][isr]->GetBinContent(ibin)*100., 
	 h_mc_sl_tot[ELE][isr]->GetBinContent(ibin), 
	 h_mc_sl_tot[ELE][isr]->GetBinError(ibin), 
	 h_mc_sl_tot[ELE][isr]->GetBinError(ibin)/
	 h_mc_sl_tot[ELE][isr]->GetBinContent(ibin)*100., 
	 h_mc_sl_tot[COMB][isr]->GetBinContent(ibin), 
	 h_mc_sl_tot[COMB][isr]->GetBinError(ibin), 
	 h_mc_sl_tot[COMB][isr]->GetBinError(ibin)/
	 h_mc_sl_tot[COMB][isr]->GetBinContent(ibin)*100.);      
  printf("\\wjets\\ \t\t & $%.1f \\pm %.1f$ ($%.1f \\%%$) & $%.1f \\pm %.1f$ ($%.1f \\%%$) & $%.1f \\pm %.1f$ ($%.1f \\%%$) \\\\\n",
  	 h_mc[WJETS][MUO][isr]->GetBinContent(ibin), 
  	 h_mc[WJETS][MUO][isr]->GetBinError(ibin), 
  	 h_mc[WJETS][MUO][isr]->GetBinError(ibin)/
  	 h_mc[WJETS][MUO][isr]->GetBinContent(ibin)*100., 
  	 h_mc[WJETS][ELE][isr]->GetBinContent(ibin), 
  	 h_mc[WJETS][ELE][isr]->GetBinError(ibin), 
  	 h_mc[WJETS][ELE][isr]->GetBinError(ibin)/
  	 h_mc[WJETS][ELE][isr]->GetBinContent(ibin)*100., 
  	 h_mc[WJETS][COMB][isr]->GetBinContent(ibin), 
  	 h_mc[WJETS][COMB][isr]->GetBinError(ibin), 
  	 h_mc[WJETS][COMB][isr]->GetBinError(ibin)/
  	 h_mc[WJETS][COMB][isr]->GetBinContent(ibin)*100.);
  printf("rare \t\t & $%.1f \\pm %.1f$ ($%.1f \\%%$) & $%.1f \\pm %.1f$ ($%.1f \\%%$) & $%.1f \\pm %.1f$ ($%.1f \\%%$) \\\\\n",
  	 h_mc_ot_tot[MUO][isr]->GetBinContent(ibin), 
  	 h_mc_ot_tot[MUO][isr]->GetBinError(ibin), 
  	 h_mc_ot_tot[MUO][isr]->GetBinError(ibin)/
  	 h_mc_ot_tot[MUO][isr]->GetBinContent(ibin)*100., 
  	 h_mc_ot_tot[ELE][isr]->GetBinContent(ibin), 
  	 h_mc_ot_tot[ELE][isr]->GetBinError(ibin), 
  	 h_mc_ot_tot[ELE][isr]->GetBinError(ibin)/
  	 h_mc_ot_tot[ELE][isr]->GetBinContent(ibin)*100., 
  	 h_mc_ot_tot[COMB][isr]->GetBinContent(ibin), 
  	 h_mc_ot_tot[COMB][isr]->GetBinError(ibin), 
  	 h_mc_ot_tot[COMB][isr]->GetBinError(ibin)/
	 h_mc_ot_tot[COMB][isr]->GetBinContent(ibin)*100.);     
  
  cout<<"\\hline"<<endl;
  cout<<"\\hline"<<endl;
  printf("Total \t\t & $%.1f \\pm %.1f$ ($%.1f \\%%$) & $%.1f \\pm %.1f$ ($%.1f \\%%$) & $%.1f \\pm %.1f$ ($%.1f \\%%$) \\\\\n",
  	 h_mc_tot[MUO][isr]->GetBinContent(ibin), 
  	 h_mc_tot[MUO][isr]->GetBinError(ibin), 
  	 h_mc_tot[MUO][isr]->GetBinError(ibin)/
  	 h_mc_tot[MUO][isr]->GetBinContent(ibin)*100., 
  	 h_mc_tot[ELE][isr]->GetBinContent(ibin), 
  	 h_mc_tot[ELE][isr]->GetBinError(ibin), 
  	 h_mc_tot[ELE][isr]->GetBinError(ibin)/
  	 h_mc_tot[ELE][isr]->GetBinContent(ibin)*100., 
  	 h_mc_tot[COMB][isr]->GetBinContent(ibin), 
  	 h_mc_tot[COMB][isr]->GetBinError(ibin), 
  	 h_mc_tot[COMB][isr]->GetBinError(ibin)/
  	 h_mc_tot[COMB][isr]->GetBinContent(ibin)*100.);      
  cout<<"\\hline"<<endl;
  cout<<"\\hline"<<endl;

  cout<<"-------------------------------------------------"<<endl;
  cout<<"*************************************************"<<endl;
  cout<<"-------------------------------------------------"<<endl;

  }


  cout<<"-------------------------------------------------"<<endl;
  cout<<"*************************************************"<<endl;
  cout<<"-------------------------------------------------"<<endl;
  cout<<"PRINT SUMMARY YIELDS TABLE"<<endl;
  if (dopreveto) cout<<"BEFORE ISO. TRK VETO!!!!"<<endl;
  if (dokfscaling) cout<<"SCALING TTDL BY K-FACTORS!!!!"<<endl;
  for (int isample=0; isample<3; ++isample) {
    cout<<"\\hline"<<endl;
    cout<<"\\hline"<<endl;  
    printf("\\multicolumn{%i}{c}{%s} \\\\\n", NSAMPLE+1, isample==MUO ? "Muon":
	   (isample==ELE ? "Electron" : "Muon+Electron Combined"));
    cout<<"\\hline"<<endl;  
    printf("\\ttdl\\ \t\t ");
    for (int isr=0; isr<NSAMPLE; ++isr) {
      printf("& $%.1f \\pm %.1f$",// ($%.1f \\%%$)",
	     h_mc[TTDL][isample][isr]->GetBinContent(ibin), 
	     h_mc[TTDL][isample][isr]->GetBinError(ibin));
      // h_mc[TTDL][isample][isr]->GetBinError(ibin)/
      // h_mc[TTDL][isample][isr]->GetBinContent(ibin)*100.);     
    }
    printf(" \\\\\n");
    
    printf("\\ttsl\\ \\& single top (1\\Lep) \t\t ");
    for (int isr=0; isr<NSAMPLE; ++isr) {
      printf("& $%.1f \\pm %.1f$",// ($%.1f \\%%$)",
	     h_mc_sl_tot[isample][isr]->GetBinContent(ibin), 
	     h_mc_sl_tot[isample][isr]->GetBinError(ibin));
      // h_mc_sl_tot[isample][isr]->GetBinError(ibin)/
      // h_mc_sl_tot[isample][isr]->GetBinContent(ibin)*100.);     
    }
    printf(" \\\\\n");
    
    printf("\\wjets\\ \t\t ");
    for (int isr=0; isr<NSAMPLE; ++isr) {
      printf("& $%.1f \\pm %.1f$",// ($%.1f \\%%$)",
	     h_mc[WJETS][isample][isr]->GetBinContent(ibin), 
	     h_mc[WJETS][isample][isr]->GetBinError(ibin));
      // h_mc[WJETS][isample][isr]->GetBinError(ibin)/
      // h_mc[WJETS][isample][isr]->GetBinContent(ibin)*100.);     
    }
    printf(" \\\\\n");
    
    printf("Rare \t\t ");
    for (int isr=0; isr<NSAMPLE; ++isr) {
      printf("& $%.1f \\pm %.1f$",// ($%.1f \\%%$)",
	     h_mc_ot_tot[isample][isr]->GetBinContent(ibin), 
	     h_mc_ot_tot[isample][isr]->GetBinError(ibin));
      // h_mc_ot_tot[isample][isr]->GetBinError(ibin)/
      // h_mc_ot_tot[isample][isr]->GetBinContent(ibin)*100.);     
    }
    printf(" \\\\\n");
    
    cout<<"\\hline"<<endl;
    
    printf("Total \t\t ");
    for (int isr=0; isr<NSAMPLE; ++isr) {
      printf("& $%.1f \\pm %.1f$",// ($%.1f \\%%$)",
	     h_mc_tot[isample][isr]->GetBinContent(ibin), 
	     h_mc_tot[isample][isr]->GetBinError(ibin));
      // h_mc_tot[isample][isr]->GetBinError(ibin)/
      // h_mc_tot[isample][isr]->GetBinContent(ibin)*100.);     
    }
    printf(" \\\\\n");
    cout<<"\\hline"<<endl;
    if (ibin!=i_sig) {
      cout<<"\\hline"<<endl;
      printf("Data \t\t ");
      for (int isr=0; isr<NSAMPLE; ++isr) {
	printf("& $%.0f$",// ($%.1f \\%%$)",
	       h_dt[isample][isr]->GetBinContent(ibin));
      }
      printf(" \\\\\n");
    }
    cout<<"\\hline"<<endl;
  }
  cout<<"-------------------------------------------------"<<endl;
  cout<<"*************************************************"<<endl;
  cout<<"-------------------------------------------------"<<endl;
  


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
