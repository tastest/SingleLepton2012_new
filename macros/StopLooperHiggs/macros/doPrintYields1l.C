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

void doPrintYields1l(char* ttbar_tag = "") {

  // this flag sets all uncertainties except for those
  // associated with the ttbar sample to zero
  // in this way it calculates the uncorrelated uncertainty
  //  bool issyst = true;
  bool issyst = (strlen(ttbar_tag) == 0) ? false : true;

  // to dump all the information about the yields for 
  // the various MC samples
  bool doverbose = false;

  //for single lepton, split single top
  const int MCID = 8; 
  const char* mcsample[MCID]={
    "ttsl_lmgtau",
    "ttdl_lmgtau",
    "extratop_sl",
    "extratop_dl",
    "ttH",
    "bosons", 
    "T6tthh_350_smallTree",
    "T6tthh_450_smallTree"};

  //legends only needed for final predictions
  const char* legend[MCID] = {
    "$t\\bar{t} \\rightarrow l^{\\pm}$ + jets",
    "$t\\bar{t} \\rightarrow l^{+}l^{-}$",
    "ttV + single top (s/t-chan) 1-lep",
    "ttV + single top (tW) 2-lep",
    "ttH",
    "V+jets + diboson + triboson", 
    "T6tthh 350-175-0", 
    "T6tthh 450-275-100"};
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

  enum sample{TTSL=0, 
	      TTDL, 
	      EXTSL, 
	      EXTDL, 
	      RARE};

  //-------------------------------
  // SINGLE LEPTON SAMPLE
  //-------------------------------

  //settings go here
  const int NSAMPLE = 3;
  char * selection[NSAMPLE] = {"_1l2b", "_1l3b", "_1l4b"};    
  const char* srname[NSAMPLE] = {"1l + 2b", "1l + 3b", "1l + 4b"};

  // bins for regions
  float i_ctr = 1;
  float i_sig = 2;

  float ibin = i_sig;

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

  //indices are for muon, electron and sum
  TH1F *h_dt[3][NSAMPLE];
  TH1F *h_mc[MCID][3][NSAMPLE];
  TH1F *h_mc_sl_tot[3][NSAMPLE];
  TH1F *h_mc_ot_tot[3][NSAMPLE];
  TH1F *h_mc_tot[3][NSAMPLE];

  for (int isr=0; isr<NSAMPLE; ++isr) {

    for (int itag=0; itag<2; ++itag) {
      h_dt[itag][isr]= (TH1F*)dt_sl[itag]->Get(Form("%s%s%s",histoname,selection[isr],leptype[itag]));
      h_dt[itag][isr]->SetName(Form("h_dt%s%s",selection[isr],leptype[itag]));
      if (issyst)	zeroHistError(h_dt[itag][isr]);
      
      for (int j=0;j<MCID;++j) {
	h_mc[j][itag][isr] = (TH1F*)mc_sl[j]->Get(Form("%s%s%s",
						       histoname,selection[isr],leptype[itag]));
	if (h_mc[j][itag][isr]==0) {
	  h_mc[j][itag][isr] = (TH1F*)mc_sl[0]->Get(Form("%s%s%s",
							 histoname,selection[isr],leptype[itag]));
	  h_mc[j][itag][isr]->SetName(Form("h_mc_%s%s%s", 
					   mcsample[j],selection[isr],leptype[itag]));
	  zeroHist(h_mc[j][itag][isr]);
	}
	//	h_mc[j][itag][isr]->Sumw2();
	
	if (doverbose) 
	  cout<<"MC "<<mcsample[j]<<" "
	      <<Form("%s%s%s",histoname,selection[isr],leptype[itag])
	      <<" "<<h_mc[j][itag][isr]<<endl;
	h_mc[j][itag][isr]->SetName(Form("h_mc_%s%s%s", 
					 mcsample[j],selection[isr],leptype[itag]));
	
	if (issyst && j>1)
	  zeroHistError(h_mc[j][itag][isr]);
	
	if (j==0) {
	  h_mc_tot[itag][isr] = (TH1F*)h_mc[j][itag][isr]->Clone();
	  h_mc_tot[itag][isr]->SetName(Form("h_mc_tot%s%s",selection[isr],leptype[itag]));
	} else if (j<6)
	  h_mc_tot[itag][isr]->Add(h_mc[j][itag][isr]);
	
	if (doverbose) {
	  cout<<"-------------------------------------------------"<<endl;
	  printf("%s MC Yields for MET cut %s and %s \n", mcsample[j], selection[isr],
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
    h_dt[COMB][isr]->SetName(Form("h_dt_comb_%s", selection[isr]));
    h_dt[COMB][isr]->Add(h_dt[ELE][isr]);
    
    for (int j=0;j<MCID;++j) {
      h_mc[j][COMB][isr]= (TH1F*)h_mc[j][MUO][isr]->Clone();
      h_mc[j][COMB][isr]->SetName(Form("h_mc_comb_%s%s", mcsample[j],selection[isr]));
      h_mc[j][COMB][isr]->Add(h_mc[j][ELE][isr]);
  
      if (doverbose) {
	cout<<"-------------------------------------------------"<<endl;
	printf("%s COMBINED MC Yields MET cut %s \n", mcsample[j], selection[isr] );
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

    h_mc_tot[COMB][isr]= (TH1F*)h_mc_tot[MUO][isr]->Clone();
    h_mc_tot[COMB][isr]->SetName(Form("h_mc_tot_comb_%s", selection[isr]));
    h_mc_tot[COMB][isr]->Add(h_mc_tot[ELE][isr]);
    
  }

  
  // cout<<"-------------------------------------------------"<<endl;
  // cout<<"*************************************************"<<endl;
  cout<<"-------------------------------------------------"<<endl;
  //  cout<<"PRINT SUMMARY YIELDS TABLE"<<endl;
  for (int isample=0; isample<3; ++isample) {
    cout<<"\\hline"<<endl;
    cout<<"\\hline"<<endl;  
    printf("\\multicolumn{%i}{c}{%s} \\\\\n", NSAMPLE+1, isample==MUO ? "Muon":
   	   (isample==ELE ? "Electron" : "Muon+Electron Combined"));
    cout<<"\\hline"<<endl;  
    for (int j=0;j<MCID;++j) {
      printf("%s \t\t ", legend[j]);
      for (int isr=0; isr<NSAMPLE; ++isr) {
	printf("& $%.2f \\pm %.2f$",// ($%.2f \\%%$)",
	       h_mc[j][isample][isr]->GetBinContent(ibin), 
	       h_mc[j][isample][isr]->GetBinError(ibin));
      }
      printf(" \\\\\n");
      if (j==5) {
	cout<<"\\hline"<<endl;  
	printf("Total \t\t ");
	for (int isr=0; isr<NSAMPLE; ++isr) {
	  printf("& $%.2f \\pm %.2f$",// ($%.2f \\%%$)",
		 h_mc_tot[isample][isr]->GetBinContent(ibin), 
		 h_mc_tot[isample][isr]->GetBinError(ibin));
	}
	printf(" \\\\\n");
	cout<<"\\hline"<<endl;  
      }
    }
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
