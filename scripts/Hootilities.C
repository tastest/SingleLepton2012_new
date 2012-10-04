#include "Hootilities.h"
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


float histError( TH1F* h , int minbin , int maxbin ){

  float err2 = 0;
  for( int i = minbin ; i <= maxbin ; ++i ){
    err2 += pow( h->GetBinError(i) , 2 );
  }

  return sqrt(err2);
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


void print( TH1F* h , string label , bool correlatedError = false ){

  stringstream se;
  stringstream sm;
  stringstream stot;

  if( label == "data" ){
    se   << Form( "%.0f" , h->GetBinContent(1) );
    sm   << Form( "%.0f" , h->GetBinContent(2) );
    stot << Form( "%.0f" , h->Integral()       );
  }else{
    //see  << Form( "%.1f" , h->GetBinContent(1) );
    //smm  << Form( "%.1f" , h->GetBinContent(2) );
    //stot << Form( "%.1f" , h->Integral()       );

    se   << Form( "%.1f" , h->GetBinContent(1) ) << pm << Form( "%.1f" , h->GetBinError(1) );
    sm   << Form( "%.1f" , h->GetBinContent(2) ) << pm << Form( "%.1f" , h->GetBinError(2) );

    float error = 0;
    if( correlatedError ) error = h->GetBinError(1) + h->GetBinError(2);
    else                  error = histError(h,1,4);
    
    stot << Form( "%.1f" , h->Integral()       ) << pm << Form( "%.1f" , error  );
  }

  cout << delimstart << setw(width1) << label      << setw(width2)
       << delim      << setw(width1) << se.str()   << setw(width2)
       << delim      << setw(width1) << sm.str()   << setw(width2)
       << delim      << setw(width1) << stot.str() << setw(width2)
       << delimend   << endl;
  
  
}


#include <TList.h>
#include <TIterator.h>

void deleteHistos() {
   // Delete all existing histograms in memory
   TObject* obj;
   TList* list = gDirectory->GetList() ;
   TIterator* iter = list->MakeIterator();
   while ((obj=iter->Next())) {
     if (obj->IsA()->InheritsFrom(TH1::Class()) ||
         obj->IsA()->InheritsFrom(TH2::Class()) ) {delete obj;}
   }
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
    pm         = " $\\pm$ ";
    delim      = "&";
    delimstart = "";
    delimend   = "\\\\";
    e          = "$e$";
    m          = "$\\mu$";
  }else{
    pm         = " +/- ";
    delim      = "|";
    delimstart = "|";
    delimend   = "|";
    e          = "e";
    m          = "m";
  }

}

//TH1F* doDYestimate( TChain *ch , TCut sel ){
void doDYestimate( TChain *ch , TCut sel , TH1F* hyield ){

  //-------------------
  // set parameters
  //-------------------

  bool  verbose = true;
  float R       = 0.13;
  float R_err   = 0.07;
  float k       = 1.11;

  if( verbose ){
    cout << "R : " << R << " +/- " << R_err << endl;
    cout << "k : " << k << endl;
  }

  //-------------------
  // invert Z-veto
  //-------------------

  TString selstring(sel.GetTitle());
  selstring.ReplaceAll("passz == 0","dilmass>76&&dilmass<106");
  TCut newsel = TCut(selstring);

  if( verbose ){
    cout << "Pre  : " << sel.GetTitle() << endl;
    cout << "Post : " << newsel.GetTitle() << endl;
  }

  //-------------------
  // get data yields
  //-------------------

  float ndata_ee_in = ch->GetEntries(newsel+"leptype==0");
  float ndata_mm_in = ch->GetEntries(newsel+"leptype==1");
  float ndata_em_in = ch->GetEntries(newsel+"leptype==2");

  if( verbose ){
    cout << "nee     " << ndata_ee_in << endl;
    cout << "nmm     " << ndata_mm_in << endl;
    cout << "nem     " << ndata_em_in << endl;
  }

  //-------------------
  // do DY estimate
  //-------------------

  float neepred     = R     * ( ndata_ee_in - (0.5/k)  * ndata_em_in );
  float nmmpred     = R     * ( ndata_mm_in - k/2.     * ndata_em_in );

  float neeprederr  = 0.;
  float nmmprederr  = 0.;

  neeprederr  += pow( R_err * ( ndata_ee_in - (0.5/k) * ndata_em_in ) , 2 );
  nmmprederr  += pow( R_err * ( ndata_mm_in - k/2.    * ndata_em_in ) , 2 );

  neeprederr  += pow( R * sqrt( ndata_ee_in + pow(0.5/k,2)  * ndata_em_in ) , 2 );
  nmmprederr  += pow( R * sqrt( ndata_mm_in + pow(k/2.,2)   * ndata_em_in ) , 2 );

  neeprederr = sqrt(neeprederr);
  nmmprederr = sqrt(nmmprederr);

  float ntotpred    = neepred + nmmpred;
  float ntotprederr = sqrt( pow(neeprederr,2) + pow(nmmprederr,2) );

  if( verbose ){
    cout << "nee  pred " << neepred  << " +/- " << neeprederr  << endl;
    cout << "nmm  pred " << nmmpred  << " +/- " << nmmprederr  << endl;
    cout << "ntot pred " << ntotpred << " +/- " << ntotprederr << endl;
  }

  //-------------------
  // set hist contents
  //-------------------

  hyield->SetBinContent( 1 , neepred );
  hyield->SetBinContent( 2 , nmmpred );
  hyield->SetBinContent( 3 , 0 );
  hyield->SetBinContent( 4 , 0 );

  hyield->SetBinError( 1 , neeprederr );
  hyield->SetBinError( 2 , nmmprederr );
  hyield->SetBinError( 3 , 0 );
  hyield->SetBinError( 4 , 0 );

}

void printYields( vector<TChain*> chmc , vector<char*> labels , TChain* chdata , TCut sel , TCut weight , bool latex ){

  initSymbols( latex );

  TCanvas *ctemp = new TCanvas();

  printLine(latex);
  printHeader();
  printLine(latex);

  TH1F* hyield = new TH1F("hyield","yield",4,0,4);
  TH1F* hmctot = new TH1F("hmctot","hmctot",4,0,4);
  hyield->Sumw2();
  hmctot->Sumw2();

  //----------------------
  // print SM MC samples
  //----------------------

  for(unsigned int imc = 0 ; imc < chmc.size() ; imc++){

    bool correlatedError = false;

    if( TString(labels[imc]).Contains("T2") ) continue;

    // data-driven DY estimate
    if( strcmp(labels[imc],"DYdata")   == 0 ){
      //hyield = doDYestimate( chmc[imc] , sel );
      doDYestimate( chmc[imc] , sel , hyield );
    }

    // fake estimate
    else if( strcmp(labels[imc],"single fakes")   == 0 || strcmp(labels[imc],"double fakes") == 0){

      //correlatedError = true;

      TString weightstring(weight.GetTitle());
      weightstring.ReplaceAll("ndavtxweight","1");
      TCut newweight = TCut(weightstring);

      chmc[imc]->Draw("leptype>>hyield",sel*newweight);

      // SF --> SF - 2 X DF
      if( strcmp(labels[imc],"single fakes")   == 0 ){
	
	TH1F *hyielddf = new TH1F("hyielddf","hyielddf",4,0,4);

	chmc[imc+1]->Draw("leptype>>hyielddf",sel*newweight);
	hyield->Add(hyielddf,-2);
      }      

      hyield->SetBinError(1,0.5*hyield->GetBinContent(1));
      hyield->SetBinError(2,0.5*hyield->GetBinContent(2));
      
    }

    //vanilla MC
    else{
      chmc[imc]->Draw("leptype>>hyield",sel*weight);

      //do efficiency correction
      //hyield->SetBinContent  ( 2 , hyield->GetBinContent(2) * 0.90);
      //hyield->SetBinContent  ( 3 , hyield->GetBinContent(3) * 0.95);
      //hyield->SetBinError    ( 2 , hyield->GetBinError(2)   * 0.90);
      //hyield->SetBinError    ( 3 , hyield->GetBinError(3)   * 0.95);

    }

    if( imc == 0 ) hmctot = (TH1F*) hyield->Clone();
    else           hmctot->Add(hyield);
    
    print( hyield , labels[imc] , correlatedError );

    //hyield->Reset();
  }

  printLine(latex);

  //-------------------------------
  // print sum of SM MC samples
  //-------------------------------

  print( hmctot , "total SM MC" );

  printLine(latex);
 
  chdata->Draw("leptype>>hyield",sel);

  print( hyield , "data" );
    
  printLine(latex);

  //----------------------
  // print T2tt MC samples
  //----------------------

  for(unsigned int imc = 0 ; imc < chmc.size() ; imc++){

    if( !TString(labels[imc]).Contains("T2") ) continue;

    chmc[imc]->Draw("leptype>>hyield",sel*weight);

    if( TString(labels[imc]).Contains("X5")  ) hyield->Scale(5);
    if( TString(labels[imc]).Contains("X6")  ) hyield->Scale(6);
    if( TString(labels[imc]).Contains("X10") ) hyield->Scale(10);

    //do efficiency correction
    //hyield->SetBinContent  ( 2 , hyield->GetBinContent(2) * 0.90);
    //hyield->SetBinContent  ( 3 , hyield->GetBinContent(3) * 0.95);
    //hyield->SetBinError    ( 2 , hyield->GetBinError(2)   * 0.90);
    //hyield->SetBinError    ( 3 , hyield->GetBinError(3)   * 0.95);
    
    print( hyield , labels[imc] );

  }

  printLine(latex);


  
  delete ctemp;
}


TLegend *getLegend( vector<TChain*> chmc , vector<char*> labels , bool overlayData, float x1, float y1, float x2, float y2){

  int colors[]={6,2,7,4,5,8,9,15,12};
  int sigcolors[]={1,1,7,4,5,8,9,15,12};
  int isigmc = 0;

  TLegend *leg = new TLegend(x1,y1,x2,y2);

  TH1F*    datahist = new TH1F("datahist","datahist",1,0,1);
  datahist->Sumw2();

  if( overlayData ) leg->AddEntry(datahist,"data");

  const int nmc = chmc.size();
  TH1F*    mchist[nmc];

  //-----------------
  // SM samples
  //-----------------

  //for( unsigned int imc = 0 ; imc < nmc ; imc++ ){
  for( int imc = nmc - 1 ; imc >= 0 ; imc-- ){

    char* t = labels.at(imc);

    if( TString(t).Contains("T2") ) continue;

    mchist[imc] = new TH1F(Form("mc_%i",imc),Form("mc_%i",imc),1,0,1);

    if( TString( labels.at(imc) ).Contains("T2") ){
      mchist[imc]->SetFillColor( 0 );
      mchist[imc]->SetLineStyle(2);
      mchist[imc]->SetLineColor( sigcolors[imc] );
    }else{
      mchist[imc]->SetFillColor( colors[imc] );
    }


    if( strcmp("ttall",t)    == 0 ) t = "t#bar{t}";
    if( strcmp("ttl",t)      == 0 ) t = "t#bar{t} #rightarrow #font[12]{l}+jets";
    if( strcmp("tt1l",t)     == 0 ) t = "t#bar{t} #rightarrow #font[12]{l^{#pm}} + jets";
    if( strcmp("ttsl",t)     == 0 ) t = "t#bar{t} #rightarrow #font[12]{l^{#pm}} + jets";
    if( strcmp("tt2l",t)     == 0 ) t = "t#bar{t} #rightarrow #font[12]{l^{+}l^{-}}";
    if( strcmp("ttdl",t)     == 0 ) t = "t#bar{t} #rightarrow #font[12]{l^{+}l^{-}}";
    if( strcmp("ttdl_lost",t)  == 0 ) t = "t#bar{t} #rightarrow #font[12]{l^{+}l^{-}} (lost)";
    if( strcmp("ttdl_lep",t)   == 0 ) t = "t#bar{t} #rightarrow #font[12]{l^{+}l^{-}} (e/#mu)";
    if( strcmp("ttdl_tauh",t)  == 0 ) t = "t#bar{t} #rightarrow #font[12]{l^{+}l^{-}} (#tau_{had})";
    if( strcmp("ttdl_tauh1",t) == 0 ) t = "t#bar{t} #rightarrow #font[12]{l^{+}l^{-}} (#tau_{had}#rightarrow1-prong)";
    if( strcmp("ttdl_tauhm",t) == 0 ) t = "t#bar{t} #rightarrow #font[12]{l^{+}l^{-}} (#tau_{had}#rightarrow3-prong)";
    if( strcmp("ttdl_taul",t)  == 0 ) t = "t#bar{t} #rightarrow #font[12]{l^{+}l^{-}} (#tau_{lep})";
    if( strcmp("ttfake",t)   == 0 ) t = "t#bar{t} #rightarrow hadrons";
    if( strcmp("ttll",t)     == 0 ) t = "t#bar{t} #rightarrow #font[12]{l}#font[12]{l}";
    if( strcmp("ttltau",t)   == 0 ) t = "t#bar{t} #rightarrow #font[1]{l}#font[12]{#tau}";
    if( strcmp("tttau",t)    == 0 ) t = "t#bar{t} #rightarrow #font[12]{#tau}+jets";
    if( strcmp("tttautau",t) == 0 ) t = "t#bar{t} #rightarrow #font[12]{#tau}#font[12]{#tau}";
    if( strcmp("ttotr",t)    == 0 ) t = "t#bar{t} #rightarrow other";
    if( strcmp("t",t)        == 0 ) t = "single top";
    if( strcmp("qcd",t)      == 0 ) t = "QCD";
    if( strcmp("wjets",t)    == 0 ) t = "W+jets";
    if( strcmp("WW",t)       == 0 ) t = "W^{+}W^{-}";
    if( strcmp("WZ",t)       == 0 ) t = "W^{#pm}Z^{0}";
    if( strcmp("ZZ",t)       == 0 ) t = "Z^{0}Z^{0}";

    //leg->AddEntry(mchist[imc],labels.at(imc),"f");
    leg->AddEntry(mchist[imc],t,"f");
    
  }

  //-----------------
  // T2tt samples
  //-----------------

  for( unsigned int imc = 0 ; imc < nmc ; imc++ ){
  //for( int imc = nmc - 1 ; imc >= 0 ; imc-- ){

    char* t = labels.at(imc);

    if( !TString(t).Contains("T2") ) continue;

    mchist[imc] = new TH1F(Form("mc_%i",imc),Form("mc_%i",imc),1,0,1);

    if( TString( labels.at(imc) ).Contains("T2") ){
      mchist[imc]->SetFillColor( 0 );
      //mchist[imc]->SetLineStyle(2);
      mchist[imc]->SetLineWidth(2);
      mchist[imc]->SetLineColor( sigcolors[isigmc++] );
    }else{
      mchist[imc]->SetFillColor( colors[imc] );
    }

    if( strcmp("ttall",t) == 0 ) t = "t#bar{t}";
    if( strcmp("t",t)     == 0 ) t = "single top";
    if( strcmp("wjets",t) == 0 ) t = "W+jets";
    if( strcmp("WW",t)    == 0 ) t = "W^{+}W^{-}";
    if( strcmp("WZ",t)    == 0 ) t = "W^{#pm}Z^{0}";
    if( strcmp("ZZ",t)    == 0 ) t = "Z^{0}Z^{0}";

    //leg->AddEntry(mchist[imc],labels.at(imc),"f");
    leg->AddEntry(mchist[imc],t,"f");
    
  }

  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  
  return leg;

}

void compareDataMC( vector<TChain*> chmc , vector<char*> labels , TChain* chdata , char* var , 
		    TCut sel , TCut weight , int nbins ,  float xmin , float xmax ,  
		    char* xtitle , bool overlayData , bool residual , bool drawLegend , bool log , char* flavor ){

  TPad* fullpad = new TPad();
  TPad* plotpad = new TPad();
  TPad* respad  = new TPad();

  if( residual ){
    fullpad = new TPad("fullpad","fullpad",0,0,1,1);
    fullpad->Draw();
    fullpad->cd();

    plotpad = new TPad("plotpad","plotpad",0,0,1,0.8);
    plotpad->Draw();
    plotpad->cd();
    if( log ) plotpad->SetLogy();
  }
  else{
    if( log ) gPad->SetLogy();
  }

  TString tvar(var);
  tvar.ReplaceAll("()","");
  tvar.ReplaceAll(".","");
  tvar.ReplaceAll("(","");
  tvar.ReplaceAll(")","");
  const char* myvar = tvar;

  cout << "Plotting var " << myvar << " flavor " << flavor << endl;

  int colors[]={6,2,7,4,5,8,9,15,12};
  int sigcolors[]={1,1,7,4,5,8,9,15,12};
  int isigmc = 0;

  assert( chmc.size() == labels.size() );
  const unsigned int nmc = chmc.size();

  THStack* mcstack = new THStack("mcstack","mcstack");
  TH1F*    mctothist = new TH1F();
  TH1F*    smtothist = new TH1F();
  TH1F*    mchist[nmc];
  TH1F*    datahist = new TH1F(Form("%s_datahist_%s",myvar,flavor),Form("%s_datahist_%s",myvar,flavor),nbins,xmin,xmax);

  vector<TH1F*> sighist;

  float trigeff = 1.0;
  //if     ( TString(flavor).Contains("ee")  ) trigeff = 1.00;
  //else if( TString(flavor).Contains("mm")  ) trigeff = 0.90;
  //else if( TString(flavor).Contains("em")  ) trigeff = 0.95;
  //else if( TString(flavor).Contains("all") ) trigeff = 0.95;

  TCut trigweight(Form("%.2f",trigeff));

  for( unsigned int imc = 0 ; imc < nmc ; imc++ ){
  //for( int imc = nmc-1 ; imc > -1 ; imc-- ){

    bool isSignal = TString( labels.at(imc) ).Contains("T2");

    mchist[imc] = new TH1F(Form("%s_mc_%i_%s",myvar,imc,flavor),Form("%s_mc_%i_%s",myvar,imc,flavor),nbins,xmin,xmax);
    mchist[imc]->Sumw2();

    chmc.at(imc)->Draw(Form("TMath::Min(%s,%f)>>%s_mc_%i_%s",var,xmax-0.01,myvar,imc,flavor),sel*weight*trigweight);

    if( isSignal ){
      mchist[imc]->SetFillColor( 0 );
      mchist[imc]->SetLineStyle(2);
      mchist[imc]->SetLineWidth(2);
      mchist[imc]->SetLineColor( sigcolors[isigmc++] );
      if( TString( labels.at(imc) ).Contains("X5") ){
	mchist[imc]->Scale(5);
	cout << "Scaling signal MC by 5" << endl;
      }
      if( TString( labels.at(imc) ).Contains("X6") ){
	mchist[imc]->Scale(6);
	cout << "Scaling signal MC by 6" << endl;
      }
      if( TString( labels.at(imc) ).Contains("X10") ){
	mchist[imc]->Scale(10);
	cout << "Scaling signal MC by 10" << endl;
      }
    }else{
      mchist[imc]->SetLineWidth(1);
      mchist[imc]->SetFillColor( colors[imc] );
      mchist[imc]->SetLineColor( 1 );
    }

    //mcstack->Add( mchist[imc] );

    if( !isSignal ){
      mcstack->Add( mchist[imc] );

      if( imc == 0 ){
	mctothist = (TH1F*) mchist[imc]->Clone();
	smtothist = (TH1F*) mchist[imc]->Clone();
      }
      else{
	mctothist->Add(mchist[imc]);
	smtothist->Add(mchist[imc]);
      }
    }
    else{
      mctothist->Add(mchist[imc]);
      sighist.push_back( mchist[imc] );
    }

    cout << "MC yield " << labels[imc] << " " << Form("%.2f",mchist[imc]->Integral()) << endl;
  }

  chdata->Draw(Form("TMath::Min(%s,%f)>>%s_datahist_%s",var,xmax-0.01,myvar,flavor),sel);

  if( overlayData ){

    float max = datahist->GetMaximum() + datahist->GetBinError(datahist->GetMaximumBin());
    if( mctothist->GetMaximum() > max ) max = mctothist->GetMaximum();
    if( log ) datahist->SetMaximum( 15 * max );
    else      datahist->SetMaximum( 1.4 * max );

    datahist->GetXaxis()->SetTitle(xtitle);
    datahist->Draw("E1");
    mcstack->Draw("samehist");
    
    for( unsigned int isig = 0 ; isig < sighist.size() ; isig++ ){
      sighist.at(isig)->Add(smtothist);
      sighist.at(isig)->Draw("samehist");
    } 

    datahist->Draw("sameE1");
    datahist->Draw("sameaxis");

    if(!log) datahist->GetYaxis()->SetRangeUser(0.,1.4*max);
    
  }
  else{
    float max = mctothist->GetMaximum();
    if( log ) mctothist->SetMaximum( 15 * max );
    else      mctothist->SetMaximum( 1.4 * max );

    mctothist->SetLineColor(0);
    mctothist->SetFillColor(0);

    mctothist->GetXaxis()->SetTitle(xtitle);
    mctothist->Draw("hist");
    mcstack->Draw("samehist");

    for( unsigned int isig = 0 ; isig < sighist.size() ; isig++ ){
      //sighist.at(isig)->Add(smtothist);
      //sighist.at(isig)->Scale( smtothist->Integral() / sighist.at(isig)->Integral() );
      sighist.at(isig)->Draw("samehist");
    } 

    mctothist->Draw("sameaxis");
  }

  if( drawLegend ){
    TLegend* myleg = getLegend( chmc , labels , overlayData );
    myleg->Draw();
  }

  TLatex *text = new TLatex();
  text->SetNDC();
  text->SetTextSize(0.04);
  //text->DrawLatex(0.2,0.88,"CMS Preliminary");
  //text->DrawLatex(0.2,0.83,"0.98 fb^{-1} at #sqrt{s} = 7 TeV");
  //text->DrawLatex(0.2,0.83,"#sqrt{s} = 7 TeV, #scale[0.6]{#int}Ldt = 5.0 fb^{-1}");

  // if     ( TString(flavor).Contains("e")  )  text->DrawLatex(0.2,0.78,"e-channel");
  // else if( TString(flavor).Contains("m")  )  text->DrawLatex(0.2,0.78,"#mu-channel");
  // else if( TString(flavor).Contains("all") ) text->DrawLatex(0.2,0.78,"e/#mu-channel");

  if( residual ){
    fullpad->cd();

    respad = new TPad("respad","respad",0,0.8,1,1);
    respad->Draw();
    respad->cd();

    //gPad->SetGridy();

    TH1F* ratio = (TH1F*) datahist->Clone(Form("%s_ratio",datahist->GetName()));
    ratio->Divide(smtothist);

    ratio->GetYaxis()->SetTitleOffset(0.3);
    ratio->GetYaxis()->SetTitleSize(0.2);
    ratio->GetYaxis()->SetNdivisions(5);
    ratio->GetYaxis()->SetLabelSize(0.2);
    ratio->GetYaxis()->SetRangeUser(0.5,1.5);
    ratio->GetYaxis()->SetTitle("data/SM  ");
    ratio->GetXaxis()->SetLabelSize(0);
    ratio->GetXaxis()->SetTitleSize(0);
    ratio->SetMarkerSize(1);
    ratio->Draw();

    TLine line;
    line.SetLineWidth(1);
    line.DrawLine(datahist->GetXaxis()->GetXmin(),1,datahist->GetXaxis()->GetXmax(),1);

  }






}
