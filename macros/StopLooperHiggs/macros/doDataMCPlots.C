//////////////////////////////////////////////////////////////////////////////
// Usage:
//
// - compile the macro in root:
//  [] .L doDataMCPlots.C++
// - launch functions for 1D and 2D plotting:
//  [] doDataMCPlots()
//
//////////////////////////////////////////////////////////////////////////////
#include <TH1F.h>
#include <TH2F.h>
#include <THStack.h>
#include <TLegend.h>
#include <TFile.h>
#include <TMath.h>
#include <TF1.h>
#include <TLine.h>
#include <TList.h>
#include <TLatex.h>
#include <TGraphAsymmErrors.h>
#include <TCanvas.h>
#include <vector>
#include <iostream>
#include <algorithm>
#include <sstream>
#include <iomanip>

using namespace std;

// This is meant to be passed as the third argument, the predicate, of the standard library sort algorithm
inline bool sortByIntegral(TH1F *h1, TH1F *h2 ) {
  return h1->Integral() < h2->Integral();
}

double error_poisson_up(double data);
double error_poisson_down(double data);

double compatibilityTest(TH1F* hist, THStack* hstack);

double straightline ( double *x, double *parm) {
  return parm[0] + parm[1]*x[0];
}

float addSqr( float e1 , float e2 ){
  float err2 = pow(e1,2) + pow(e2,2);
  return sqrt(err2);
}

void zeroHist( TH1F* &h ) {

  for (int i=1; i<=h->GetNbinsX(); ++i) {
    h->SetBinContent(i,0.); 
    h->SetBinError(i,0.); 
  }
}

TGraphAsymmErrors* GetPoissonizedGraph(TH1F* histo);

void myText(Double_t x,Double_t y,Color_t color,char *text)
{

  TLatex l;
  l.SetNDC();
  l.SetTextColor(color);
  l.SetTextSize(0.05);
  l.DrawLatex(x,y,text);
}

void myOtherText(Double_t x,Double_t y,double tsize,char *text)
{
  TLatex l;
  l.SetTextSize(tsize);
  l.SetNDC();
  l.DrawLatex(x,y,text);
}

void myCMSLabel(Double_t x,Double_t y, int type) {

  // CMS Label types:
  // 0: no label
  // 1: CMS

  if (type==0)
    return;
  
  TLatex l;
  l.SetNDC();
  l.SetTextFont(72);
  l.SetTextSize(0.04); 
  l.DrawLatex(x,y,"CMS");
  l.SetTextFont(42);
  l.SetTextSize(0.04);
  l.DrawLatex(x+0.13,y,"Preliminary");

}

void printYields( vector<TH1F*> mc , const char* labels[] , TH1F* chdata , bool latex = false );
void initSymbols(bool);
void printLine(bool);
void printHeader();
void print( TH1F* h , string label , bool correlatedError = false );
char* pm;         
char* percent;         
char* delim;      
char* delimstart; 
char* delimend;   
char* e;         
char* m;         
int   width1;
int   width2;
int   linelength;

void scaleHistError( TH1F* &h, float sf, float esf ) {

  float ierr = 0.;
  for (int i=1; i<=h->GetNbinsX(); ++i) {
    ierr = h->GetBinError(i)/h->GetBinContent(i); 
    ierr = addSqr( ierr, esf/sf) * h->GetBinContent(i);
    h->SetBinError(i,ierr); 
  }

}

void zeroHistError( TH1F* &h ) {

  for (int i=1; i<=h->GetNbinsX(); ++i) 
    h->SetBinError(i,0.); 
}

float histError( TH1F* h , int minbin , int maxbin ){

  float err2 = 0;
  for( int i = minbin ; i <= maxbin ; ++i ){
    err2 += pow( h->GetBinError(i) , 2 );
  }

  return sqrt(err2);
}

pair<float, float> datamcsf( TH1F* h_dt, vector<TH1F*> h_mc, int minbin, int maxbin ) {
  
  float n_dt = h_dt->Integral(minbin, maxbin+1);
  float e_dt = histError( h_dt, minbin, maxbin+1);
  //Add up MC samples
  if (h_mc.size()<1) return make_pair(0.,0.);

  int nbins = h_mc.at(0)->GetNbinsX();
  float xmin = h_mc.at(0)->GetBinLowEdge(1);
  float xmax = h_mc.at(0)->GetBinLowEdge(nbins);
  cout<<"Going from "<<minbin<<" to "<<maxbin
      <<" nbin range "<<h_mc.at(0)->GetXaxis()->GetBinLowEdge(minbin)
      <<" - "<<h_mc.at(0)->GetXaxis()->GetBinUpEdge(maxbin)<<endl;
  TH1F* hmctot = new TH1F("hmctot","hmctot", nbins, xmin, xmax);
  hmctot->Sumw2();

  for(unsigned int imc = 0 ; imc < h_mc.size() ; imc++){
    h_mc.at(imc)->Sumw2();
    if( imc == 0 ) hmctot = (TH1F*) h_mc.at(imc)->Clone();
    else           hmctot->Add(h_mc.at(imc));
  }

  float n_mc = hmctot->Integral(minbin, maxbin);
  float e_mc = histError( hmctot, minbin, maxbin );

  float sf_dtmc = n_dt/n_mc;
  float e_sf_dtmc = sqrt(pow( (e_dt/n_dt), 2) + pow( (e_mc/n_mc), 2 )) * sf_dtmc;

  return make_pair(sf_dtmc, e_sf_dtmc);
}


//////////////////////////////////////////////////////////////////////////////
void doDataMCPlots(const char* ttbar_tag = "_lmgtau") {

  //for single lepton, split single top
  const int MCID = 3; 
  const char* mcsample[MCID]={
    "ttdl",
    "ttsl",
    "rare"};
    // "T6tthh_350_smallTree",
    // "T6tthh_450_smallTree"};

  //legends only needed for final predictions
  const char* legend[MCID] = {
    "t#bar{t}#rightarrow #font[12]{ll}",
    "t#bar{t}#rightarrow #font[12]{l}+jets",
    "rare"};
    // "T6tthh 350-175-0", 
    // "T6tthh 450-275-100"};

  enum sample{TTDL=0, 
	      TTSL, 
	      RARE};

  const int mccolor[]={kAzure+1,2,kGreen+1,6,4,kOrange,8,kAzure-9, kGreen+1, kViolet,kGreen+1, 15,12,13,27};

  const int NSAMPLE = 1;
  char * selection[NSAMPLE] = {"_2l2b"};    
  const char* samplename[NSAMPLE] = {"2l + 2b"};

  TFile *dt_sl = TFile::Open("SIGoutput/data_dl_histos.root");
  
  TFile *mc_sl[MCID];   
  for (int j=0;j<MCID;++j) {
    if (j<2)
      mc_sl[j] = TFile::Open(Form("SIGoutput/%s%s_histos.root",
				  mcsample[j], ttbar_tag));
    else 
      mc_sl[j] = TFile::Open(Form("SIGoutput/%s_histos.root",
				  mcsample[j]));
  }

  double lumi = 19.5;

  for (int isr=0; isr<NSAMPLE; ++isr) {
    
    const int N1DHISTS = 7;
    
    TH1F *h_dt1d[N1DHISTS];
    TH1F *h_mc1d[N1DHISTS][MCID];
    vector<TH1F*> sorted_mc1d[N1DHISTS];
    TH1F *h_mc1d_tot[N1DHISTS];
    
    // Output file names
    const char* file1dname[N1DHISTS]= {
      "met",
      "mbb", 
      "leppt", 
      "mbb_count",
      "dphi_metlep", 
      "bpt1",
      "bpt2"};
    
    // Histogram names
    const char* hist1dname[N1DHISTS] ={
      "h_sig_met",
      "h_sig_mbb", 
      "h_sig_leppt", 
      "h_sig_mbb_count",
      "h_sig_dphi_metlep", 
      "h_sig_bpt1", 
      "h_sig_bpt2"};
    
    // List of Log scale plots:
    vector<int> logScale;
    logScale.push_back(0);
    logScale.push_back(1);
    logScale.push_back(2);
    logScale.push_back(3);
    logScale.push_back(4);
    logScale.push_back(5);
    logScale.push_back(6);
    logScale.push_back(29);
    
    // End MC Samples definition and open files ///////////////////////////////
    // End getting all info needed: K3, K4 and mtsf ///////////////////////////
    // Open and book histograms ///////////////////////////////////////////////
    
    // List of rebin factors:
    int rebinFactor[N1DHISTS] = {3,2,2,1,1,5,5};
    if (isr>4) 
      rebinFactor[1] = 6;
    
    const char* xtitle1d[N1DHISTS]= {
      //    "N Jets",
      "ME_{T} [GeV]", 
      "M_{bb} [GeV]", 
      "Lepton p_{T} [GeV]", 
      "M_{bb} Range", 
      "#Delta #phi(ME_{T}, lep) [rad]", 
      "Lead b p_{T} [GeV]", 
      "Sub-Lead b p_{T} [GeV]"}; 
    
      const char* ytitle1d[N1DHISTS]= {
	//    "",
	"GeV",
	"GeV",
	"GeV",
	"",
	// "GeV",
	// "GeV",
	// "GeV", 
	// "GeV", 
	"rad", 
	"GeV", 
	"GeV"};

      // Book histograms
      for (int i=0;i<N1DHISTS;++i) {
	h_dt1d[i] = (TH1F*)dt_sl->Get(Form("%s%s", hist1dname[i], selection[isr]));
	cout << "DT: " << h_dt1d[i] << " " << Form("%s%s", hist1dname[i], selection[isr]) << endl;
	h_dt1d[i]->SetName(Form("%s",hist1dname[i]));
	h_dt1d[i]->Rebin(rebinFactor[i]);

	bool doinit = false;
	for (int j=0;j<MCID;++j) {

	  h_mc1d[i][j] =
	    (TH1F*)mc_sl[j]->Get(Form("%s%s",hist1dname[i], selection[isr]));
	  if (h_mc1d[i][j]==0) {
	    h_mc1d[i][j] =
	      (TH1F*)mc_sl[0]->Get(Form("%s%s",hist1dname[i], selection[isr]));
	    h_mc1d[i][j]->SetName(Form("%s_%s",mcsample[j],hist1dname[i]));
	    zeroHist(h_mc1d[i][j]);
	  }
	  h_mc1d[i][j]->SetName(Form("%s_%s",mcsample[j],hist1dname[i]));
      
	  cout << "MC" << j << " "<<mcsample[j]<<": " <<h_mc1d[i][j] << " " 
	       << Form("%s",hist1dname[i]) << endl;

	  h_mc1d[i][j]->Sumw2();

	  // Rebin, set limits
	  h_mc1d[i][j]->Rebin(rebinFactor[i]);
      
	  //get total
	  if (!doinit) {
	    h_mc1d_tot[i] = (TH1F*)h_mc1d[i][j]->Clone(Form("mctot_%s",hist1dname[i]));
	    doinit = true;
	  } else h_mc1d_tot[i]->Add(h_mc1d[i][j]);

	}
      }

      //after applying all the scalings, fill vector for sorting
      for (int i=0;i<N1DHISTS;++i) {
	cout<<"-------------------------------------------------"<<endl;
	float mcall = h_mc1d_tot[i]->Integral();
	float dtall = h_dt1d[i]->Integral();
	cout << file1dname[i] << " MC ALL " << mcall <<endl;
	cout << file1dname[i] << "Data ALL " << dtall <<endl;
	float mcsf_all = dtall/mcall;
	cout << "RATIO ALL " << mcsf_all <<endl;
	cout<<"-------------------------------------------------"<<endl;
	
	for (int j=0;j<MCID;++j) 
	  sorted_mc1d[i].push_back( h_mc1d[i][j] );
	//sort mc histograms
	sort( sorted_mc1d[i].begin(), sorted_mc1d[i].end(), sortByIntegral );
      }

      // End of open and book histograms ////////////////////////////////////////

      // Style up histograms ////////////////////////////////////////////////////

      cout << "Styling 1D" << endl;
      for (int i=0;i<N1DHISTS;++i) {
	h_dt1d[i]->SetLineColor(kBlack);
	h_dt1d[i]->SetMarkerColor(kBlack);
	h_dt1d[i]->SetFillColor(0);
	h_dt1d[i]->SetFillStyle(0);
	for (int j=0;j<MCID;++j) {
	  h_mc1d[i][j]->SetLineColor(kBlack);
	  h_mc1d[i][j]->SetMarkerColor(mccolor[j]);
	  h_mc1d[i][j]->SetFillColor(mccolor[j]);
	  //      h_mc1d[i][j]->SetFillColor(j>0?mccolor[j]:10);
	  h_mc1d[i][j]->SetFillStyle(1001);
	}
    
	h_dt1d[i]->SetLineColor(kBlack);
	h_dt1d[i]->SetMarkerColor(kBlack);
	h_dt1d[i]->SetFillColor(0);
	h_dt1d[i]->SetFillStyle(0);
	for (int j=0;j<MCID;++j) {
	  h_mc1d[i][j]->SetLineColor(kBlack);
	  h_mc1d[i][j]->SetMarkerColor(mccolor[j]);
	  h_mc1d[i][j]->SetFillColor(mccolor[j]);
	  //      h_mc1d[i][j]->SetFillColor(j>0?mccolor[j]:10);
	  h_mc1d[i][j]->SetFillStyle(1001);
	}
      }
      // End of style up histograms /////////////////////////////////////////////

      // Add titles to axes /////////////////////////////////////////////////////

      cout << "Setting Titles" << endl;
      for (int i=0;i<N1DHISTS;++i) {
	if (i==3) {
	  h_dt1d[i]->GetXaxis()->SetBinLabel(1,"Out Low");	
	  h_dt1d[i]->GetXaxis()->SetBinLabel(2,"In");	
	  h_dt1d[i]->GetXaxis()->SetBinLabel(3,"Out High");	
	  h_dt1d[i]->GetXaxis()->SetTitle("M_{bb} Region");
	  h_dt1d[i]->GetYaxis()->SetTitle("Entries");
	} else {
	  h_dt1d[i]->GetXaxis()->SetTitle(Form("%s",
					       xtitle1d[i]));
	  h_dt1d[i]->GetXaxis()->SetTitleOffset(1.);
	  h_dt1d[i]->GetYaxis()->SetTitle(Form("Entries / %3.0f %s",
					       h_dt1d[i]->GetBinWidth(1),
					       ytitle1d[i]));
	  h_dt1d[i]->GetXaxis()->SetTitle(Form("%s",
						    xtitle1d[i]));
	  h_dt1d[i]->GetYaxis()->SetTitle(Form("Entries / %3.0f %s",
					       h_dt1d[i]->GetBinWidth(1),
						    ytitle1d[i]));
	}

	if (find(logScale.begin(),logScale.end(),i)!=logScale.end()) {
	  h_dt1d[i]->GetYaxis()->SetTitleOffset(1.1);
	}
	else {
	  h_dt1d[i]->GetYaxis()->SetTitleOffset(1.1);
	}

	for (int j=0;j<MCID;++j) {
      
	  h_mc1d[i][j]->GetXaxis()->SetTitle(Form("%s",
						  xtitle1d[i]));
	  h_mc1d[i][j]->GetXaxis()->SetTitleOffset(1.);
	  h_mc1d[i][j]->GetYaxis()->SetTitle(Form("Entries / %3.0f %s",
						  h_dt1d[i]->GetBinWidth(1),
						  ytitle1d[i]));
	  if (find(logScale.begin(),logScale.end(),i)!=logScale.end()) {
	    h_mc1d[i][j]->GetYaxis()->SetTitleOffset(1.1);
	  }
	  else {
	    h_mc1d[i][j]->GetYaxis()->SetTitleOffset(1.1);
	  }
      
	}
      }
      // End of add titles to axes /////////////////////////////////////////////

      // Print expected number of MC events, after all normalizations
      for (int j=0;j<MCID;++j) {
	cout << endl << "Expected " << legend[j] 
	     << ": " << h_mc1d[1][j]->Integral(1,h_mc1d[1][j]->GetNbinsX())
	     << endl;
      }
      cout << endl << "Expected Data " << h_dt1d[1]->Integral(1,h_dt1d[1]->GetNbinsX())
	   << endl;

      // Stack histograms ////////////////////////////////////////////////////
      THStack *s_mc1d[N1DHISTS];

      for (int i=0;i<N1DHISTS;++i) {
	s_mc1d[i] = new THStack(hist1dname[i],hist1dname[i]);
	s_mc1d[i]->SetName(Form("stack_%s",hist1dname[i]));
	bool doinit = false;
	for (int j=0;j<MCID;++j) {
	  s_mc1d[i]->Add(sorted_mc1d[i].at(j));
	}
	cout << endl << "MC " << h_mc1d_tot[i]->Integral(1,h_mc1d_tot[i]->GetNbinsX())
	     << endl;
	cout << endl << "Data " << h_dt1d[1]->Integral(1,h_dt1d[1]->GetNbinsX())
	     << endl;
	cout<<endl;
	// if (i==0) {
	//   TFile *outfile_sl=new TFile("mt_leadmuo_njge4.root","RECREATE"); 
	//   s_mc1d[i]->Write();
	//   h_dt1d[i]->Write();
	//   outfile_sl->Write();
	//   outfile_sl->Close();
	// }


      }
      // End of stack histograms /////////////////////////////////////////////

      // Legends!
      cout << "Doing Legends" << endl;
      TLegend *leg1d[N1DHISTS];

      // Will make them in same position, put "if(i==..)" to change them
      for (int i=0;i<N1DHISTS;++i) {
	leg1d[i] = new TLegend(0.71, 0.60, 0.93, 0.931);//0.56,0.64,0.92,0.915);
	leg1d[i]->SetName(Form("leg_%s",h_dt1d[i]->GetName()));
	leg1d[i]->SetFillColor(0);
	leg1d[i]->
	  AddEntry(h_dt1d[i],
		   "data ","lp");

	for (int j=0;j<MCID;++j) {
	  leg1d[i]->AddEntry(h_mc1d[i][j],legend[j],"f");
	}
      }

      // Canvases, at last!
      TCanvas *canv1d[N1DHISTS];

      for (int i=0;i<N1DHISTS;++i) {

	canv1d[i] = new TCanvas(Form("c_%s",file1dname[i]),
				Form("c_%s",file1dname[i]),
				600,700);
	canv1d[i]->SetLeftMargin(0.18);
	canv1d[i]->SetRightMargin(0.05);

	TPad* fullpad = new TPad();
	TPad* plotpad = new TPad();
	TPad* respad  = new TPad();
    
	fullpad = new TPad("fullpad","fullpad",0,0,1,1);
	fullpad->Draw();
	fullpad->cd();
    
	plotpad = new TPad("plotpad","plotpad",0,0,1,0.8);
	plotpad->Draw();
	plotpad->cd();

	if (find(logScale.begin(),logScale.end(),i)!=logScale.end()) {
	  plotpad->SetLogy();
	  float maxval = h_dt1d[i]->GetMaximum()*100.;
	  if (i==4 || i==5 || i==6) maxval = h_dt1d[i]->GetMaximum()*5000.;
	  h_dt1d[i]->SetMaximum(maxval);
	  float minval = 5e-2;
	  if (i==0 || i==2 || isr>2) minval = 5e-2;
	  h_dt1d[i]->SetMinimum(minval);//5e-1);
	  s_mc1d[i]->SetMinimum(minval);//5e-1);
	} else {
	  h_dt1d[i]->SetMaximum(h_dt1d[i]->GetMaximum()*2.);      
	  s_mc1d[i]->SetMaximum(h_dt1d[i]->GetMaximum()*2.);      
	  s_mc1d[i]->SetMinimum(0.);
	  h_dt1d[i]->SetMinimum(0.);
	}

	TH1F *h_basic = (TH1F*)h_dt1d[i]->Clone();
	h_basic->SetName("h_basic");
	h_basic->SetMarkerColor(kWhite);
	h_basic->SetMarkerStyle(1);
	h_basic->SetMarkerSize(0.0000000001);
	h_basic->Draw("e3");
	s_mc1d[i]->Draw("hist,same");

	TGraphAsymmErrors* g_data = GetPoissonizedGraph(h_dt1d[i]);
	// g_data->SetMarkerStyle(20);
	// g_data->SetMarkerSize(0.9);
	// g_data->SetLineWidth(2);
	g_data->Draw("e1,p,z,same");
	// h_dt1d[i]->SetMarkerStyle(20);
	// h_dt1d[i]->SetMarkerSize(0.9);
	// h_dt1d[i]->Draw("p,same");
	h_dt1d[i]->Draw("axis,same");

	leg1d[i]->Draw("same");

	// Labels! 
	//    myCMSLabel(0.2,0.88,CMSLabel);
	TLatex *text = new TLatex();
	text->SetNDC();
	text->SetTextSize(0.04);
	float xtex = 0.2;//0.16;//used to be 0.2
	text->DrawLatex(xtex,0.88,"CMS Preliminary");
	text->DrawLatex(xtex,0.83,Form("#sqrt{s} = 8 TeV, #scale[0.6]{#int}Ldt = %.1f fb^{-1}", lumi));
	text->DrawLatex(xtex,0.78,Form("%s",samplename[isr]));

	fullpad->cd();

	respad = new TPad("respad","respad",0,0.8,1,1);
	respad->Draw();
	respad->cd();
    
	//gPad->SetGridy();
	cout<<"----------------------------------------------"<<endl;
	TH1F* ratio = (TH1F*) h_dt1d[i]->Clone(Form("%s_ratio",h_dt1d[i]->GetName()));
	ratio->Divide(h_mc1d_tot[i]);
    
	ratio->GetYaxis()->SetTitleOffset(0.3);
	ratio->GetYaxis()->SetTitleSize(0.2);
	ratio->GetYaxis()->SetNdivisions(5);
	ratio->GetYaxis()->SetLabelSize(0.2);
	//    if (dozoom) ratio->GetYaxis()->SetRangeUser(0.8,1.2);
	//    if (dozoom) ratio->GetYaxis()->SetRangeUser(0.7,1.3);
	//    else 
	//    ratio->GetYaxis()->SetRangeUser(0.7,1.3);
	ratio->GetYaxis()->SetRangeUser(0.,2.);
	if (i==3) ratio->GetYaxis()->SetRangeUser(0.5,1.5);
	//ratio->GetYaxis()->SetRangeUser(0.5,1.5);
	ratio->GetYaxis()->SetTitle("data/SM  ");
	ratio->GetXaxis()->SetLabelSize(0);
	ratio->GetXaxis()->SetTitleSize(0);
	//	ratio->SetMarkerSize(1);
	ratio->Draw();

	// float xmin = ratio->GetXaxis()->GetBinLowEdge(1);
	// float xmax = ratio->GetXaxis()->GetBinUpEdge(ratio->GetNbinsX());

	// TF1 *f_line = new TF1("f_line", straightline, xmin, xmax, 2);
	// ratio->Fit("f_line");
	// f_line->SetLineColor(kRed);
	// f_line->SetLineWidth(2);
	// f_line->Draw("SAME");

	TLine line;
	line.SetLineWidth(1);
	line.DrawLine(h_dt1d[i]->GetXaxis()->GetXmin(),1,h_dt1d[i]->GetXaxis()->GetXmax(),1);

	canv1d[i]->Print(Form("SIGplots/%s%s%s.C",
			      file1dname[i], selection[isr],
			      ttbar_tag));
	canv1d[i]->Print(Form("SIGplots/%s%s%s.png",
			      file1dname[i], selection[isr],
			      ttbar_tag));
	canv1d[i]->Print(Form("SIGplots/%s%s%s.pdf",
			      file1dname[i], selection[isr],
			      ttbar_tag));
	canv1d[i]->Print(Form("SIGplots/%s%s%s.pdf",
			      file1dname[i], selection[isr],
			      ttbar_tag));
	
	cout << " Probability " << file1dname[i] << " : "
	     << compatibilityTest(h_dt1d[i],s_mc1d[i]) << endl;
	
	//	if (i==6) break;
	//delete canv1d[i];
      }
  }
}

//////////////////////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////////////////////
  TGraphAsymmErrors* GetPoissonizedGraph(TH1F* histo) {

    TGraphAsymmErrors* graph = new TGraphAsymmErrors();
    graph->SetName(Form("g_%s",histo->GetName()));

    int j=1;
    for (int i=1;i<=histo->GetNbinsX();++i) {

      if (histo->GetBinContent(i)!=0) {
  
	graph->SetPoint(j,histo->GetBinCenter(i),histo->GetBinContent(i));

	graph->SetPointError(j,
			     histo->GetBinWidth(i)/2.,
			     histo->GetBinWidth(i)/2.,
			     error_poisson_down(histo->GetBinContent(i)),
			     error_poisson_up(histo->GetBinContent(i)));

	++j;
      }
    }
    return graph;

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
  //////////////////////////////////////////////////////////////////////////////


  // New function for comparing errors

  double compatibilityTest(TH1F* hist, THStack* hstack) {

    // A bit of a disaster: need to loop on list of histograms in stack

    TList* stack_hists = hstack->GetHists();

    int ndf = 0;
    double chisq = 0.;
    double h_entries = 0.; // entries in histo
    double s_entries = 0.; // entries in stack
    double s_sqerror = 0.;
    for (int i = 1; i<=hist->GetNbinsX(); ++i) {
    
      h_entries += hist->GetBinContent(i);
    
      TIter iter(stack_hists->MakeIterator());
      while (TH1F* stack_hist = dynamic_cast<TH1F*>(iter())) {
	s_entries += stack_hist->GetBinContent(i);
	s_sqerror += (stack_hist->GetBinError(i)*
		      stack_hist->GetBinError(i));
      }
      //minimum requirement used to be 20, now lowered to 10
      if (h_entries<10 && i<=hist->GetNbinsX()) 
	continue;
      else {
	ndf++;
      
	chisq += (h_entries-s_entries)*(h_entries-s_entries)/
	  (s_sqerror + h_entries);
      

	h_entries = s_entries = s_sqerror = 0.;
      }
    
    } 
    cout << "Chisq, ndf: " << chisq << ", " << ndf << endl; 
    return TMath::Prob(chisq,ndf);
  }


  //------------------------------------
  // print yield table
  //------------------------------------
  void printYields( vector<TH1F*> h_mc , const char* labels[] , TH1F* h_data , bool latex ){

    initSymbols( latex );//for latex plots

    TCanvas *ctemp = new TCanvas();

    printLine(latex);
    printHeader();
    printLine(latex);

    TH1F* hmctot = new TH1F("hmctot","hmctot",2,0,2);
    hmctot->Sumw2();

    //----------------------
    // print SM MC samples
    //----------------------

    for(unsigned int imc = 0 ; imc < h_mc.size() ; imc++){

      print( h_mc.at(imc) , labels[imc] , false );

      if( imc == 0 ) hmctot = (TH1F*) h_mc.at(imc)->Clone();
      else           hmctot->Add(h_mc.at(imc));
    }

    printLine(latex);

    //-------------------------------
    // print sum of SM MC samples
    //-------------------------------

    print( hmctot , "total SM MC" );

    printLine(latex);
 
    print( h_data , "data" );
    
    printLine(latex);

    delete ctemp;
  }

  void initSymbols( bool latex ){

    //-------------------------------------------------------
    // table format
    //-------------------------------------------------------

    width1      = 20;
    width2      = 4;
    linelength  = (width1+width2)*4+1;

    //-------------------------------------------------------
    // symbols
    //-------------------------------------------------------
  
    if( latex ){
      percent    = " \\% ";
      pm         = " $\\pm$ ";
      delim      = "&";
      delimstart = "";
      delimend   = "\\\\";
      e          = "$e$";
      m          = "$\\mu$";
    }else{
      percent    = " % ";
      pm         = " +/- ";
      delim      = "|";
      delimstart = "|";
      delimend   = "|";
      e          = "e";
      m          = "m";
    }

  }

  void printLine( bool latex ){

    if( latex ){
      cout << "\\hline" << endl;
    }
    else{
      for( int i = 0 ; i < linelength ; ++i ) cout << "-";
      cout << endl;
    }
  }

  void printHeader(){

    cout << delimstart << setw(width1) << "Sample"    << setw(width2)
	 << delim      << setw(width1) << e           << setw(width2)
	 << delim      << setw(width1) << m           << setw(width2)
	 << delim      << setw(width1) << "total"     << setw(width2) 
	 << delimend   << endl;

  }

  void print( TH1F* h , string label , bool correlatedError ){

    stringstream se;
    stringstream sm;
    stringstream stot;

    if( label == "data" ){
      se   << Form( "%.0f" , h->GetBinContent(1) );
      sm   << Form( "%.0f" , h->GetBinContent(2) );
      stot << Form( "%.0f" , h->GetBinContent(3) );
    }else{
      //see  << Form( "%.1f" , h->GetBinContent(1) );
      //smm  << Form( "%.1f" , h->GetBinContent(2) );
      //stot << Form( "%.1f" , h->Integral()       );

      se   << Form( "%.1f" , h->GetBinContent(1) ) << pm << Form( "%.1f" , h->GetBinError(1) );
      sm   << Form( "%.1f" , h->GetBinContent(2) ) << pm << Form( "%.1f" , h->GetBinError(2) );
      stot << Form( "%.1f" , h->GetBinContent(3) ) << pm << Form( "%.1f" , h->GetBinError(3) );

      float error = 0;
      if( correlatedError ) error = h->GetBinError(1) + h->GetBinError(2) + h->GetBinError(3);
      else                  error = histError(h,1,4);
    
      //    stot << Form( "%.1f" , h->Integral()       ) << pm << Form( "%.1f" , error  );
    }

    cout << delimstart << setw(width1) << label      << setw(width2)
	 << delim      << setw(width1) << se.str()   << setw(width2)
	 << delim      << setw(width1) << sm.str()   << setw(width2)
	 << delim      << setw(width1) << stot.str() << setw(width2)
	 << delimend   << endl;
  
  
  }
