#include <algorithm>
#include <iostream>
#include <iomanip>
#include <vector>
#include <math.h>
#include <fstream>
#include <sstream>

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

void compareHisto( TCut sel, TChain* chdl, TChain* chsl, TChain* chsig , char* var, char* name , int nbins , float xmin, float xmax ){

  TH1F* hdl  = new TH1F(Form("%s_dl",var) ,Form("%s_dl",var) ,nbins,xmin,xmax);
  TH1F* hsl  = new TH1F(Form("%s_sl",var) ,Form("%s_sl",var) ,nbins,xmin,xmax);
  TH1F* hsig = new TH1F(Form("%s_sig",var),Form("%s_sig",var),nbins,xmin,xmax);

  hdl->Sumw2();
  hsl->Sumw2();
  hsig->Sumw2();

  TCanvas *ctemp = new TCanvas();
  ctemp->cd();
  chdl ->Draw(Form("min(%s,%f)>>%s_dl" ,var,xmax-0.0001,var),sel);
  chsl ->Draw(Form("min(%s,%f)>>%s_sl" ,var,xmax-0.0001,var),sel);
  chsig->Draw(Form("min(%s,%f)>>%s_sig",var,xmax-0.0001,var),sel);
  delete ctemp;

  TCanvas *c1 = new TCanvas();
  c1->cd();

  hdl->Scale(1.0/hdl->Integral());
  hsl->Scale(1.0/hsl->Integral());
  hsig->Scale(1.0/hsig->Integral());

  hsl->SetLineColor(2);
  hsig->SetLineColor(4);
  
  float max = hdl->GetMaximum();
  if( hsl->GetMaximum() > max ) max = hsl->GetMaximum();
  if( hsig->GetMaximum() > max ) max = hsig->GetMaximum();
  hdl->SetMaximum(1.1*max);

  hdl->GetXaxis()->SetTitle(name);
  hdl->Draw("hist");
  hsl->Draw("samehist");
  hsig->Draw("samehist");

  cout << endl;
  cout << "Variable       : " << var << endl;
  cout << "tt->ll         : " << hdl->GetBinContent(2) << endl;
  cout << "tt->l+jets     : " << hsl->GetBinContent(2) << endl;
  cout << "T2tt (450,0)   : " << hsig->GetBinContent(2) << endl;

}


void compareVariables(){

  //-----------------------------------------------------
  // define event selection (store in TCut sel)
  //-----------------------------------------------------

  TCut njets4("njets>=4");
  TCut met100("met>=100");
  TCut mt120("mt>=120");
  TCut nb1("nb>=1");
  TCut isotrk("passisotrk==1");
  TCut lep_pt30("nlep>=1 && lep1pt>30.0");
   
  //TCut  sel      = lep_pt30 + isotrk + njets4 + met100 + mt120 + nb1;
  TCut  sel      = lep_pt30 + njets4 + met100 + mt120 + nb1;

  TChain* chdl  = new TChain("t");
  TChain* chsl  = new TChain("t");
  TChain* chsig = new TChain("t");

  chdl->Add("output/ttdl_powheg_mini.root");
  chsl->Add("output/ttsl_powheg_mini.root");
  chsig->Add("output/T2tt_450_0_mini.root");

  vector<char*> variables;
  vector<char*> names;
  vector<int>   nbins;
  vector<float> xmin;
  vector<float> xmax;

  //variables.push_back("met");                   names.push_back("E_{T}^{miss} [GeV]"); nbins.push_back(50); xmin.push_back(0.0); xmax.push_back(500.0);
  variables.push_back("passisotrk");            names.push_back("iso-track (old)");    nbins.push_back(2);  xmin.push_back(0);   xmax.push_back(2);
  variables.push_back("passisotrkv2");          names.push_back("iso-track (new)");    nbins.push_back(2);  xmin.push_back(0);   xmax.push_back(2);
  // variables.push_back("passisotrk");                    
  // variables.push_back("passisotrkv2");                    
  // variables.push_back("mt");                    
  // variables.push_back("mt2w");                    
  // variables.push_back("mt2bl");                    
  // variables.push_back("mt2b");                    
  // variables.push_back("min(chi2,100)");                    
  // variables.push_back("lep1pt");                                     
  // variables.push_back("lep1eta");                                    
  // variables.push_back("thrjetlm");                                   
  // variables.push_back("apljetlm");                                   
  // variables.push_back("sphjetlm");                                   
  // variables.push_back("cirjetlm");                                   
  // variables.push_back("min(chi2min,100)");                           
  // variables.push_back("chi2min_mt2b");                               
  // variables.push_back("chi2min_mt2bl");                              
  // variables.push_back("chi2min_mt2w");                               
  // variables.push_back("mt2bmin");                                    
  // variables.push_back("mt2blmin");                                   
  // variables.push_back("mt2wmin");                                    
  // variables.push_back("mt2bmin_chi2");                               
  // variables.push_back("mt2blmin_chi2");                              
  // variables.push_back("mt2wmin_chi2");                               
  // variables.push_back("htssl/(htssl+htosl)");                       
  // variables.push_back("htssm/(htssm+htosm)");                       
  // variables.push_back("dphimj1");                                    
  // variables.push_back("dphimj2");                                    
  // variables.push_back("rand");                                       

  for( int i = 0 ; i < variables.size() ; ++i){
    compareHisto( sel , chdl , chsl , chsig , variables.at(i) , names.at(i) , nbins.at(i) , xmin.at(i) , xmax.at(i) );
  }


}
