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

//Main part of the macro
void doSimplePredUnc() {

  bool doverbose = false;

  //apply k-factor corrections - to be added for systematic samples
  //  bool dokfactscale = false;

  const int NLEP = 2;
  enum leptontype{MUO = 0, ELE};
  string leptag[NLEP] = {"_muo", "_ele"};

  enum hbin{PEAK = 1, TAIL};

  //list of sample names
  const int MCID = 5; 
  const char* mcsample[MCID]={
    "ttdl_lpowheg",
    "ttsl_lpowheg",
    "w1to4jets",
    "tW_lepsl",
    "rare"};
  
  enum sample{TTDL=0, 
	      TTSL, 
	      WJETS, 
	      TWSL, 
	      RARE};
  
  //input histogram name
  string histoname = "h_sig_mt_count";
  const int NSR = 4;
  string selection[NSR] = {"_met100", "_met150", "_met200", "_met250"};    
  string trkiso_tag = "_wisotrk";
  string ttbar_tag = "";
  //low mass
  //string samplename[NSR] = {"LM SRB", "LM SRC", "LM SRD", "LM SRF"};
  //high mass
  string samplename[NSR] = {"HM SRA", "HM SRB", "HM SRC", "HM SRD"};

  //Inputs
  //uncertainty on the w+jets and rare background normalization - the same for all signal regions
  float xs_unc_wjets = 0.5;
  float xs_unc_rare = 0.5;

  //Dilepton uncertainties
  //rel unc. on K3, K4 	
  float e_kfact = 0.03;
  //unc on 2nd lep veto for dil -from the T&P studies (applies only to those where it matter)
  float e_2ndlepveto = 0.06;
  float frac_2ndlepveto = 0.38;
  //unc on tau veto - guess for now
  float e_2ndtauveto = 0.07;
  float frac_2ndtauveto = 0.20;

  //dil uncertainty from MC variations - to cover the range of variations from scale up/down
  //low mass values
  //float e_dl_vars[NSR] = {0.10, 0.15, 0.25, 0.40};
  //high mass values
  float e_dl_vars[NSR] = {0.15, 0.20, 0.30, 0.45};

  //from closure tests in CR1 - W+Jets control sample
  //THESE NUMBERS NEED TO BE UPDATED
  //constant numbers
  float in_ttp_wj_sf[NSR] = { 1.2, 1.2, 1.2, 1.2 };
  float in_ttp_wjets_aunc[NSR] = { 0.3, 0.3, 0.3, 0.3 };
  //low mass numbers
  // float in_ttp_wj_sf[NSR] = { 1.102,  1.442,  1.439,  0.630};
  // float in_ttp_wjets_aunc[NSR] = { 0.280,  0.508,  0.673,  0.583};
  //high mass numbers
  // float in_ttp_wj_sf[NSR] = { 1.277,  1.046,  1.166,  1.134};
  // float in_ttp_wjets_aunc[NSR] = { 0.264,  0.275,  0.402,  0.461}; 

  //ftop/fW numbers used to obtain ttp_ttsl_sf from in_ttp_wj_sf
  //note last number is cooked up -- no entries in MC
  //low mass values
  // float in_ftopW[NSR] = { 3.894,  5.501,  9.873,  20.};
  // float in_ftopW_aunc[NSR] = { 0.591,  1.537,  4.942,  20.};
  //High mass values
  float in_ftopW[NSR] = { 1.857,  2.460,  5.234,  9.003};
  float in_ftopW_aunc[NSR] = { 0.140,  0.392,  1.696,  4.631}; 

  //calculate SFRtop and its uncertainties (correlated part is 100% correlated with SFRwjets uncertainty)
  float ttp_ttsl_sf[NSR];
  float ttp_ttsl_sf_aunc[NSR];
  float ttp_ttsl_sf_runc[NSR];
  for (int isr = 0; isr < NSR; ++isr)
  {
  	ttp_ttsl_sf[isr] = 1.0 + (in_ttp_wj_sf[isr]-1.0)*in_ftopW[isr];
  	float ttp_ttsl_sf_aunc_correlated = in_ftopW[isr]*in_ttp_wjets_aunc[isr];
  	float ttp_ttsl_sf_aunc_uncorrelated = (in_ttp_wj_sf[isr]-1.0)*in_ftopW_aunc[isr];
  	ttp_ttsl_sf_aunc[isr] = sqrt( pow(ttp_ttsl_sf_aunc_correlated,2) + pow(ttp_ttsl_sf_aunc_uncorrelated,2) );
  	if(true) cout<< "SFRtop " <<samplename[isr]<<": "
		     <<ttp_ttsl_sf[isr] << " +/- " <<ttp_ttsl_sf_aunc_correlated 
		     << " +/- " <<ttp_ttsl_sf_aunc_uncorrelated << " ( +/- " <<ttp_ttsl_sf_aunc[isr] 
		     << " total)" <<endl;
	ttp_ttsl_sf_runc[isr] = ttp_ttsl_sf_aunc[isr]/ttp_ttsl_sf[isr];
  }
  cout<<endl;

  //Fill SFRwjets to be used in later calculations. For SRs with SFRwjets relative uncertainty > 50%, we take SFRwjets from the tightest SR with relative uncertainty <= 50%.
  float ttp_wj_sf_constantabove50pcunc[NSR];
  float ttp_wjets_aunc_constantabove50pcunc[NSR]; 
  int SR_l50pcunc = 0;
  //flag to set which ttbar prediction is used, based on the same criterion as above.
  //doav true --> average between optimistic and pessimistic
  //doav false --> decomposition method and CR1 SF
  //we are now using the average method everywhere, with the decomposition method used as a cross-check
  bool doav[NSR];

  for (int isr = 0; isr < NSR; ++isr)
  {
  	if(in_ttp_wjets_aunc[isr]/in_ttp_wj_sf[isr] <= 0.4) {
  		SR_l50pcunc = isr;
  		doav[isr] = true;
  		cout<<"INFO::using Av. Estimate for 1-lepton ttbar in "<<samplename[isr]<<endl;
  		//cout<<"INFO::using Decomposition method for 1-lepton ttbar in "<<samplename[isr]<<endl;
  		ttp_wj_sf_constantabove50pcunc[isr] = in_ttp_wj_sf[isr];
  		ttp_wjets_aunc_constantabove50pcunc[isr] = in_ttp_wjets_aunc[isr];
  	}
  	else {
  		doav[isr] = true;
  		cout<<"INFO::using Av. Estimate for 1-lepton ttbar in "<<samplename[isr]<<endl;
  		ttp_wj_sf_constantabove50pcunc[isr] = in_ttp_wj_sf[SR_l50pcunc];
  		ttp_wjets_aunc_constantabove50pcunc[isr] = in_ttp_wjets_aunc[SR_l50pcunc];		
  	}
  }
  cout<<endl;


  float ttp_tt_av[NLEP][NSR];
  float ttp_tt_av_runc[NLEP][NSR];
  float ttp_tt_corr[NLEP][NSR];
  float ttp_tt_sf_runc[NLEP][NSR];
  float ttp_wjets_comb[NSR];
  float e_ttp_wjets_comb[NSR];
  float ttp_ttsl_comb[NSR];
  float e_ttp_ttsl_comb[NSR];

  //Open files
  TFile *dt_sl[2];
  TFile *mc_sl[MCID];
  dt_sl[MUO] = TFile::Open("SIGoutput/data_muo_histos.root");
  dt_sl[ELE] = TFile::Open("SIGoutput/data_ele_histos.root");
  
  for (int j=0;j<MCID;++j) {
    if (j<2)
      mc_sl[j] = TFile::Open(Form("SIGoutput/%s%s_histos.root",
  				  mcsample[j], ttbar_tag.c_str()));
    else 
      mc_sl[j] = TFile::Open(Form("SIGoutput/%s_histos.root",
  				  mcsample[j]));
  }
  TFile *mc_sl_alt[2];
  mc_sl_alt[TTDL] = TFile::Open("SIGoutput/ttdl_lmg_histos.root");
  mc_sl_alt[TTSL] = TFile::Open("SIGoutput/ttsl_lmg_histos.root");

  float pred_ttdl[NLEP][NSR];
  float e_pred_ttdl[NLEP][NSR];
  float pred_top1l[NLEP][NSR];
  float e_pred_top1l[NLEP][NSR];
  float pred_wjets[NLEP][NSR];
  float e_pred_wjets[NLEP][NSR];
  float pred_rare[NLEP][NSR];
  float e_pred_rare[NLEP][NSR];
  float pred_total_final[NLEP][NSR];
  float e_pred_total_final[NLEP][NSR];
  float pred_ttdl_comb[NSR];
  float e_pred_ttdl_comb[NSR];
  float e_mc_pred_ttdl_comb[NSR];
  float pred_top1l_comb[NSR];
  float e_pred_top1l_comb[NSR];
  float pred_wjets_comb[NSR];
  float e_pred_wjets_comb[NSR];
  float pred_rare_comb[NSR];
  float e_pred_rare_comb[NSR];
  float pred_total_final_comb[NSR];
  float e_pred_total_final_comb[NSR];

  //scale factors
  float pv_sf[NLEP][NSR];
  float e_pv_sf[NLEP][NSR];
  float sf[NLEP][NSR];
  float e_sf[NLEP][NSR];

  //truth fractions
  float frac_dl_lepveto[NSR];
  float e_frac_dl_lepveto[NSR];
  float frac_dl_tauveto[NSR];
  float e_frac_dl_tauveto[NSR];
  float frac_dl_outacc[NSR];
  float e_frac_dl_outacc[NSR];
  //alternative samples
  float pred_ttdl_alt[NLEP][NSR];
  float e_pred_ttdl_alt[NLEP][NSR];
  float frac_dl_lepveto_alt[NSR];
  float e_frac_dl_lepveto_alt[NSR];
  float frac_dl_tauveto_alt[NSR];
  float e_frac_dl_tauveto_alt[NSR];
  float frac_dl_outacc_alt[NSR];
  float e_frac_dl_outacc_alt[NSR];

  for (int isr=0; isr<NSR; ++isr) {
    //postveto
    TH1F *h_dt[NLEP];
    TH1F *h_mc[MCID][NLEP];
    TH1F *h_mc_prebtag[MCID][NLEP];
    //preveto
    TH1F *h_pv_dt[NLEP];
    TH1F *h_pv_mc[MCID][NLEP];
    //alternative sample
    TH1F *h_pv_mc_alt[2][NLEP];
    TH1F *h_mc_alt[2][NLEP];

    for (int il=0; il<NLEP; ++il) {
      
      //pass veto
      h_dt[il]= (TH1F*)dt_sl[il]->Get(Form("%s%s%s",histoname.c_str(),
					   selection[isr].c_str(),
					   leptag[il].c_str()));
      h_dt[il]->SetName(Form("h_dt_%s",leptag[il].c_str()));
      
      //preveto is sum of pass and fail veto
      h_pv_dt[il] = (TH1F*)h_dt[il]->Clone(Form("h_pv_dt_%s",leptag[il].c_str()));
      TH1F *h_tmp = (TH1F*)dt_sl[il]->Get(Form("%s%s%s%s",histoname.c_str(),
					       trkiso_tag.c_str(),selection[isr].c_str(),
					       leptag[il].c_str()));
      h_tmp->SetName("h_tmp");
      h_pv_dt[il]->Add(h_tmp);
      
      for (int imc=0;imc<MCID;++imc) {
	
	//pass veto
	h_mc[imc][il] = (TH1F*)mc_sl[imc]->Get(Form("%s%s%s",histoname.c_str(),
						    selection[isr].c_str(),
						    leptag[il].c_str()));

	if (h_mc[imc][il]==0) {
	  h_mc[imc][il] = (TH1F*)mc_sl[0]->Get(Form("%s%s%s",histoname.c_str(),
						    selection[isr].c_str(),
						    leptag[il].c_str()));
	  h_mc[imc][il]->SetName(Form("h_mc_%s%s", mcsample[imc],leptag[il].c_str()));
	  zeroHist(h_mc[imc][il]);
	}

	//alternative MC
	if (imc<2) {
	  h_mc_alt[imc][il] = (TH1F*)mc_sl_alt[imc]->Get(Form("%s%s%s",histoname.c_str(),
							  selection[isr].c_str(),
							  leptag[il].c_str()));
	  if (h_mc_alt[imc][il]==0) {
	    h_mc_alt[imc][il] = 
	      (TH1F*)h_mc[imc][il]->Clone(Form("h_mc_alt_%s%s", mcsample[imc],leptag[il].c_str()));
	    h_mc_alt[imc][il]->SetName(Form("h_mc_alt_%s%s", mcsample[imc],leptag[il].c_str()));
	    zeroHist(h_mc_alt[imc][il]);
	  }
	}

	//pass veto, prebtag
	h_mc_prebtag[imc][il] = (TH1F*)mc_sl[imc]->Get(Form("%s_prebtag%s%s",histoname.c_str(),
						    selection[isr].c_str(),
						    leptag[il].c_str()));
	if (h_mc_prebtag[imc][il]==0) {
	  h_mc_prebtag[imc][il] = (TH1F*)mc_sl[0]->Get(Form("%s_prebtag%s%s",histoname.c_str(),
						    selection[isr].c_str(),
						    leptag[il].c_str()));
	  h_mc_prebtag[imc][il]->SetName(Form("h_mc_prebtag_%s%s", mcsample[imc],leptag[il].c_str()));
	  zeroHist(h_mc[imc][il]);
	}

	//preveto is sum of pass and fail veto
	h_pv_mc[imc][il] = 
	  (TH1F*)h_mc[imc][il]->Clone(Form("h_pv_mc_%s%s", mcsample[imc],leptag[il].c_str()));
	TH1F *h_tmp_mc = (TH1F*)mc_sl[imc]->Get(Form("%s%s%s%s",histoname.c_str(),
						     trkiso_tag.c_str(),selection[isr].c_str(),
						     leptag[il].c_str()));
	if (h_tmp_mc!=0) {
	  h_tmp_mc->SetName("h_tmp_mc");
	  h_pv_mc[imc][il]->Add(h_tmp_mc);
	}

	//alternative MC
	if (imc<2) {
	  h_pv_mc_alt[imc][il] = 
	    (TH1F*)h_mc_alt[imc][il]->Clone(Form("h_pv_mc_alt_%s%s", 
						 mcsample[imc],leptag[il].c_str()));
	  TH1F *h_tmp_mc_alt = 
	    (TH1F*)mc_sl_alt[imc]->Get(Form("%s%s%s%s",histoname.c_str(),
					    trkiso_tag.c_str(),selection[isr].c_str(),
					    leptag[il].c_str()));
	  if (h_tmp_mc_alt!=0) {
	    h_tmp_mc_alt->SetName("h_tmp_mc_alt");
	    h_pv_mc_alt[imc][il]->Add(h_tmp_mc_alt);
	  }
	}
	
      }//end loop over samples
    }//end loop over leptons	
    
    //raw ttp from MC
    ttp_wjets_comb[isr] = 
      (h_mc_prebtag[WJETS][MUO]->GetBinContent(TAIL)+h_mc_prebtag[WJETS][ELE]->GetBinContent(TAIL))/
      (h_mc_prebtag[WJETS][MUO]->GetBinContent(PEAK)+h_mc_prebtag[WJETS][ELE]->GetBinContent(PEAK));
    TH1F *h_wjets_tmp = (TH1F*)h_mc_prebtag[WJETS][MUO]->Clone("h_wjets_tmp");
    h_wjets_tmp->Add(h_mc_prebtag[WJETS][ELE]);
    e_ttp_wjets_comb[isr] = sqrt( pow( (h_wjets_tmp->GetBinError(TAIL)/h_wjets_tmp->GetBinContent(TAIL)), 2) +
				  pow( (h_wjets_tmp->GetBinError(PEAK)/h_wjets_tmp->GetBinContent(PEAK)), 2) ) * 
      ttp_wjets_comb[isr];
    //to check what we get before the b-tagging requirement
    // TH1F *h_mc_wj_prebtag[NLEP];

    // for (int il=0; il<NLEP; ++il) {
    //   h_mc_wj_prebtag[il] = (TH1F*)mc_sl[WJETS]->Get(Form("%s_prebtag%s%s",histoname.c_str(),
    // 							selection[isr].c_str(),
    // 							leptag[il].c_str()));
    //   if (h_mc_wj_prebtag[il]==0) {
    // 	h_mc_wj_prebtag[il] = (TH1F*)mc_sl[0]->Get(Form("%s_prebtag%s%s",histoname.c_str(),
    // 							selection[isr].c_str(),
    // 							leptag[il].c_str()));
    // 	h_mc_wj_prebtag[il]->SetName(Form("h_mc_wj_prebtag%s", leptag[il].c_str()));
    // 	zeroHist(h_mc_wj_prebtag[il]);
    //   }    
    // }

    // float n_wj = h_mc_wj_prebtag[MUO]->GetBinContent(TAIL) + h_mc_wj_prebtag[ELE]->GetBinContent(TAIL) 
    //   + h_mc_wj_prebtag[MUO]->GetBinContent(TAIL) + h_mc_wj_prebtag[ELE]->GetBinContent(TAIL);
    // float d_wj = h_mc_wj_prebtag[MUO]->GetBinContent(PEAK) + h_mc_wj_prebtag[ELE]->GetBinContent(PEAK) 
    //   + h_mc_wj_prebtag[MUO]->GetBinContent(PEAK) + h_mc_wj_prebtag[ELE]->GetBinContent(PEAK);
    // float ttp_prebtag_wj = n_wj/d_wj;
    // cout<<ttp_prebtag_wj<<" "<<ttp_wjets_comb[isr]<<endl;
    float ttp_wjets[NLEP] = {ttp_wjets_comb[isr],ttp_wjets_comb[isr]};
    if (doverbose) cout<<"WJETS TTP "<<ttp_wjets_comb[isr]<<endl;
    //from b-veto control sample - fully correlated with electron 
    float ttp_wjets_sf[NLEP] = {ttp_wj_sf_constantabove50pcunc[isr], ttp_wj_sf_constantabove50pcunc[isr]};
    float ttp_wjets_aunc[NLEP] = {ttp_wjets_aunc_constantabove50pcunc[isr], ttp_wjets_aunc_constantabove50pcunc[isr]};
    float ttp_wjets_runc[NLEP];
    float ttp_wjets_corr[NLEP];
    for (int il = 0; il<NLEP; il++) {
      ttp_wjets_runc[il] = ttp_wjets_aunc[il]/ttp_wjets_sf[il];
      ttp_wjets_corr[il] = ttp_wjets[il] * ttp_wjets_sf[il];
    }
    
    //raw ttp from MC
    ttp_ttsl_comb[isr] = 
      (h_mc_prebtag[TTSL][MUO]->GetBinContent(TAIL)+h_mc_prebtag[TTSL][ELE]->GetBinContent(TAIL))/
      (h_mc_prebtag[TTSL][MUO]->GetBinContent(PEAK)+h_mc_prebtag[TTSL][ELE]->GetBinContent(PEAK));
    TH1F *h_ttsl_tmp = (TH1F*)h_mc_prebtag[TTSL][MUO]->Clone("h_ttsl_tmp");
    h_ttsl_tmp->Add(h_mc_prebtag[TTSL][ELE]);
    e_ttp_ttsl_comb[isr] = sqrt( pow( (h_ttsl_tmp->GetBinError(TAIL)/h_ttsl_tmp->GetBinContent(TAIL)), 2) +
				  pow( (h_ttsl_tmp->GetBinError(PEAK)/h_ttsl_tmp->GetBinContent(PEAK)), 2) ) * 
      ttp_ttsl_comb[isr];
    float ttp_tt_mc[NLEP] = {ttp_ttsl_comb[isr],ttp_ttsl_comb[isr]};
    if (doverbose) cout<<"TTSL TTP "<<ttp_ttsl_comb[isr]<<endl;
    //in between optimistic and pessimistic
    //half difference between optimistic and pessimistic + SF uncertainty from b-veto
    for (int il = 0; il<NLEP; il++) {
      ttp_tt_av[il][isr] = (ttp_tt_mc[il]+ttp_wjets[il])/2.*ttp_wjets_sf[il];
      float e_ttpsum = addSqr(e_ttp_wjets_comb[isr],e_ttp_ttsl_comb[isr])/(ttp_wjets_comb[isr]+ttp_ttsl_comb[isr]);
      ttp_tt_av_runc[il][isr] = sqrt( pow(ttp_wjets[il]-ttp_tt_mc[il],2)/
				      pow(ttp_wjets[il]+ttp_tt_mc[il],2) +
				      pow(ttp_wjets_runc[il],2) +
				      pow(e_ttpsum, 2) );
    }
    //alternative calculation based on decomposition of component with true MT passing requirement
    float ttp_tt_sf[NLEP]   = {ttp_ttsl_sf[isr], ttp_ttsl_sf[isr]};
    float ttp_tt_sf_aunc[NLEP] = {ttp_ttsl_sf_aunc[isr], ttp_ttsl_sf_aunc[isr]};
    for (int il = 0; il<NLEP; il++) {
      ttp_tt_sf_runc[il][isr] = ttp_tt_sf_aunc[il]/ttp_tt_sf[il];
      ttp_tt_corr[il][isr] = ttp_tt_mc[il] * ttp_tt_sf[il];
      //cout<< (1. +ttp_wjets[il]/ttp_tt_mc[il])/2.<<endl;
      if(il==MUO && true) cout<<samplename[isr]
			      <<":  SFRtop*ttp_tt: "<<ttp_tt_corr[il][isr]
			      <<" +/- "<<100.*ttp_tt_sf_runc[il][isr]
			      <<"%  optimistic and pessimistic average: "<<ttp_tt_av[il][isr]
			      <<" +/- "<<100.*ttp_tt_av_runc[il][isr]<<"%"<<endl;
    }
    float ttp_tt[NLEP];
    float ttp_tt_runc[NLEP];
    for (int il = 0; il<NLEP; il++) {
      ttp_tt[il] = doav[isr] ? ttp_tt_av[il][isr] : ttp_tt_corr[il][isr];
      ttp_tt_runc[il] = doav[isr] ? ttp_tt_av_runc[il][isr] : ttp_tt_sf_runc[il][isr];
    }

    //dead reckoning dil BG - must be already already corrected for K3 etc
    float mc_dl[NLEP] = {h_mc[TTDL][MUO]->GetBinContent(TAIL), 
			 h_mc[TTDL][ELE]->GetBinContent(TAIL)};
    //stat. error on dil BG in number of events
    float e_mc_dl[NLEP] = {h_mc[TTDL][MUO]->GetBinError(TAIL), 
			   h_mc[TTDL][ELE]->GetBinError(TAIL)};
    
    //dead reckoning rare MC - rare MC = tWdl + ttV + diboson + triboson + DY		
    float mc_rare[NLEP] = {h_mc[RARE][MUO]->GetBinContent(TAIL), 
			   h_mc[RARE][ELE]->GetBinContent(TAIL)};
    
    //frac of dil where lep veto matters - 1/3 of sample
    //    float frac_dl_lepveto = 0.33;
    //TO UPDATE check if has a truth e/mu down to 5 GeV
    //Add tau uncertainty
    string wtrk_tag = "_wtruetrk";
    string wtau_tag = "_notruetrk_wtruetau";
    TH1F *h_ttdl_wtrk[NLEP];
    TH1F *h_ttdl_wtau[NLEP];
    //alternative sample - note also need deak reckoning
    TH1F *h_ttdl_wtrk_alt[NLEP];
    TH1F *h_ttdl_wtau_alt[NLEP];
    TH1F *h_ttdl_alt[NLEP];
    for (int il = 0; il<NLEP; il++) {
      h_ttdl_wtrk[il] = (TH1F*)mc_sl[TTDL]->Get(Form("%s%s%s%s",histoname.c_str(),
						     wtrk_tag.c_str(),
						     selection[isr].c_str(),
						     leptag[il].c_str()));
      h_ttdl_wtrk[il]->SetName(Form("h_mc_ttdl%s%s", 
				    wtrk_tag.c_str(),leptag[il].c_str()));
      h_ttdl_wtau[il] = (TH1F*)mc_sl[TTDL]->Get(Form("%s%s%s%s",histoname.c_str(),
						     wtau_tag.c_str(),
						     selection[isr].c_str(),
						     leptag[il].c_str()));
      h_ttdl_wtau[il]->SetName(Form("h_mc_ttdl%s%s", 
				    wtau_tag.c_str(),leptag[il].c_str()));
      //alternative MC sample
      h_ttdl_alt[il] = (TH1F*)mc_sl_alt[TTDL]->Get(Form("%s%s%s",histoname.c_str(),
							selection[isr].c_str(),
							leptag[il].c_str()));
      h_ttdl_alt[il]->SetName(Form("h_mc_ttdl_alt%s",leptag[il].c_str()));
      h_ttdl_wtrk_alt[il] = (TH1F*)mc_sl_alt[TTDL]->Get(Form("%s%s%s%s",histoname.c_str(),
							     wtrk_tag.c_str(),
							     selection[isr].c_str(),
							     leptag[il].c_str()));
      h_ttdl_wtrk_alt[il]->SetName(Form("h_mc_ttdl_alt%s%s", 
					wtrk_tag.c_str(),leptag[il].c_str()));
      h_ttdl_wtau_alt[il] = (TH1F*)mc_sl_alt[TTDL]->Get(Form("%s%s%s%s",
							     histoname.c_str(),
							     wtau_tag.c_str(),
							     selection[isr].c_str(),
							     leptag[il].c_str()));
      h_ttdl_wtau_alt[il]->SetName(Form("h_mc_ttdl_alt%s%s", 
					wtau_tag.c_str(),leptag[il].c_str()));
    }
    float mc_dl_tot = mc_dl[MUO]+mc_dl[ELE];
    float e_mc_dl_tot = addSqr(e_mc_dl[MUO], e_mc_dl[ELE]);
    float invsqrtn = e_mc_dl_tot/mc_dl_tot;
    frac_dl_lepveto[isr] = ( h_ttdl_wtrk[MUO]->GetBinContent(TAIL)
			     +h_ttdl_wtrk[ELE]->GetBinContent(TAIL) )/mc_dl_tot;
    e_frac_dl_lepveto[isr] = sqrt(frac_dl_lepveto[isr]*(1.-frac_dl_lepveto[isr]))*invsqrtn;
    frac_dl_tauveto[isr] = ( h_ttdl_wtau[MUO]->GetBinContent(TAIL)
			     +h_ttdl_wtau[ELE]->GetBinContent(TAIL) )/mc_dl_tot;
    e_frac_dl_tauveto[isr] = sqrt(frac_dl_tauveto[isr]*(1.-frac_dl_tauveto[isr]))*invsqrtn;
    
    frac_dl_outacc[isr] = 1.-frac_dl_lepveto[isr]-frac_dl_tauveto[isr];
    e_frac_dl_outacc[isr] =  sqrt(frac_dl_outacc[isr]*(1.-frac_dl_outacc[isr]))*invsqrtn;

    //alternative sample
    float mc_dl_alt[NLEP] = { h_ttdl_alt[MUO]->GetBinContent(TAIL),
			      h_ttdl_alt[ELE]->GetBinContent(TAIL) }; 
    float e_mc_dl_alt[NLEP] = { h_ttdl_alt[MUO]->GetBinError(TAIL),
				h_ttdl_alt[ELE]->GetBinError(TAIL) }; 
    float mc_dl_comb_alt = mc_dl_alt[MUO] + mc_dl_alt[ELE];
    float e_mc_dl_comb_alt = addSqr(e_mc_dl_alt[MUO],e_mc_dl_alt[ELE]);
    float invsqrtn_alt = e_mc_dl_comb_alt/mc_dl_comb_alt;
    frac_dl_lepveto_alt[isr] = ( h_ttdl_wtrk_alt[MUO]->GetBinContent(TAIL)
			     +h_ttdl_wtrk_alt[ELE]->GetBinContent(TAIL) )/mc_dl_comb_alt;
    e_frac_dl_lepveto_alt[isr] = sqrt(frac_dl_lepveto_alt[isr]*
				      (1.-frac_dl_lepveto_alt[isr]))*invsqrtn_alt;
    frac_dl_tauveto_alt[isr] = ( h_ttdl_wtau_alt[MUO]->GetBinContent(TAIL)
				 +h_ttdl_wtau_alt[ELE]->GetBinContent(TAIL) )/mc_dl_comb_alt;
    e_frac_dl_tauveto_alt[isr] = sqrt(frac_dl_tauveto_alt[isr]*
				      (1.-frac_dl_tauveto_alt[isr]))*invsqrtn_alt;
    
    frac_dl_outacc_alt[isr] = 1.-frac_dl_lepveto_alt[isr]-frac_dl_tauveto_alt[isr];
    e_frac_dl_outacc_alt[isr] =  sqrt(frac_dl_outacc_alt[isr]*
				      (1.-frac_dl_outacc_alt[isr]))*invsqrtn_alt;

    if (doverbose) {
      cout<<"-------------------------------------------------"<<endl;
      cout<<"Fraction of events with a true e/mu "<<endl;
      cout<<"or charged track from single prong tau decay "<<endl;
      cout<<"with pT>(5)10 for (leptons)hadrons in |eta|<2.5 is "
      	  <<frac_dl_lepveto[isr]<<" pm "<<e_frac_dl_lepveto[isr]
      	  <<endl;
      cout<<"Fraction of events with a true tau "<<endl;
      cout<<"with pT>20 in|eta|<2.5 (and not with single particle) is "
      	  <<frac_dl_tauveto[isr]<<" pm "<<e_frac_dl_tauveto[isr]
      	  <<endl;
      cout<<"-------------------------------------------------"<<endl;
    }

    //Here we enter the number of events in the control region in data and various MC.  
    //Also the stat uncertainties	
    //Yields pre-veto in control region
    //sample	N	unc
    float pv_data[NLEP] = {h_pv_dt[MUO]->GetBinContent(PEAK), 
			   h_pv_dt[ELE]->GetBinContent(PEAK)};
    float pv_ttsl[NLEP] = {h_pv_mc[TTSL][MUO]->GetBinContent(PEAK), 
			   h_pv_mc[TTSL][ELE]->GetBinContent(PEAK)};
    float e_pv_ttsl[NLEP] = {h_pv_mc[TTSL][MUO]->GetBinError(PEAK), 
			     h_pv_mc[TTSL][ELE]->GetBinError(PEAK)};
    float pv_wjets[NLEP] = {h_pv_mc[WJETS][MUO]->GetBinContent(PEAK), 
			    h_pv_mc[WJETS][ELE]->GetBinContent(PEAK)};
    float e_pv_wjets[NLEP] = {h_pv_mc[WJETS][MUO]->GetBinError(PEAK), 
			      h_pv_mc[WJETS][ELE]->GetBinError(PEAK)};
    float pv_st[NLEP] = {h_pv_mc[TWSL][MUO]->GetBinContent(PEAK), 
			 h_pv_mc[TWSL][ELE]->GetBinContent(PEAK)};
    float e_pv_st[NLEP] = {h_pv_mc[TWSL][MUO]->GetBinError(PEAK), 
			   h_pv_mc[TWSL][ELE]->GetBinError(PEAK)};
    float pv_rare[NLEP] = {h_pv_mc[RARE][MUO]->GetBinContent(PEAK), 
			   h_pv_mc[RARE][ELE]->GetBinContent(PEAK)};
    float e_pv_rare[NLEP] = {h_pv_mc[RARE][MUO]->GetBinError(PEAK), 
			     h_pv_mc[RARE][ELE]->GetBinError(PEAK)};
    float pv_ttdl[NLEP] = {h_pv_mc[TTDL][MUO]->GetBinContent(PEAK), 
			   h_pv_mc[TTDL][ELE]->GetBinContent(PEAK)};
    float e_pv_ttdl[NLEP] = {h_pv_mc[TTDL][MUO]->GetBinError(PEAK), 
			     h_pv_mc[TTDL][ELE]->GetBinError(PEAK)};
    float pv_total[NLEP];
    float e_pv_total[NLEP];
    for (int il = 0; il<NLEP; il++) {
      pv_total[il] = pv_ttsl[il] + pv_wjets[il] + pv_st[il] + pv_rare[il] + pv_ttdl[il];
      e_pv_total[il] = sqrt( pow(e_pv_ttsl[il],2) + pow(e_pv_wjets[il],2) + 
			     pow(e_pv_st[il],2) + pow(e_pv_rare[il],2) + 
			     pow(e_pv_ttdl[il],2) );
    }

    //Yields post veto in control region
    //sample	N	unc
    float data[NLEP] = {h_dt[MUO]->GetBinContent(PEAK), 
			h_dt[ELE]->GetBinContent(PEAK)};
    float ttsl[NLEP] = {h_mc[TTSL][MUO]->GetBinContent(PEAK), 
			h_mc[TTSL][ELE]->GetBinContent(PEAK)};
    float e_ttsl[NLEP] = {h_mc[TTSL][MUO]->GetBinError(PEAK), 
			  h_mc[TTSL][ELE]->GetBinError(PEAK)};
    float wjets[NLEP] = {h_mc[WJETS][MUO]->GetBinContent(PEAK), 
			 h_mc[WJETS][ELE]->GetBinContent(PEAK)};
    float e_wjets[NLEP] = {h_mc[WJETS][MUO]->GetBinError(PEAK), 
			   h_mc[WJETS][ELE]->GetBinError(PEAK)};
    float st[NLEP] = {h_mc[TWSL][MUO]->GetBinContent(PEAK), 
		      h_mc[TWSL][ELE]->GetBinContent(PEAK)};
    float e_st[NLEP] = {h_mc[TWSL][MUO]->GetBinError(PEAK), 
			h_mc[TWSL][ELE]->GetBinError(PEAK)};
    float rare[NLEP] = {h_mc[RARE][MUO]->GetBinContent(PEAK), 
			h_mc[RARE][ELE]->GetBinContent(PEAK)};
    float e_rare[NLEP] = {h_mc[RARE][MUO]->GetBinError(PEAK), 
			  h_mc[RARE][ELE]->GetBinError(PEAK)};
    float ttdl[NLEP] = {h_mc[TTDL][MUO]->GetBinContent(PEAK), 
			h_mc[TTDL][ELE]->GetBinContent(PEAK)};
    float e_ttdl[NLEP] = {h_mc[TTDL][MUO]->GetBinError(PEAK), 
			  h_mc[TTDL][ELE]->GetBinError(PEAK)};
    float total[NLEP];
    float e_total[NLEP];
    for (int il = 0; il<NLEP; il++) {
      total[il] = ttsl[il] + wjets[il] + st[il] + rare[il] + ttdl[il];
      e_total[il] = sqrt( pow(e_ttsl[il],2) + pow(e_wjets[il],2) + 
			  pow(e_st[il],2) + pow(e_rare[il],2) + 
			  pow(e_ttdl[il],2) );
    }

    //The scale factors below do not include the systematics on wjets and rare.
    //This is because there are correlations everywhere.
    //What you want to do for these guys is to scale them up and down using the SF at the top of the macro and then see how the bottom line changes

    //PreVeto ScaleFactor (nominal)
    float e_pv_sf_rel[NLEP];
    for (int il = 0; il<NLEP; il++) {
      pv_sf[il][isr] = (pv_data[il]-pv_rare[il])/(pv_ttsl[il] + pv_wjets[il] + pv_st[il] + pv_ttdl[il]);
      e_pv_sf[il][isr] = 1./(pv_ttsl[il] + pv_wjets[il] + pv_st[il] + pv_ttdl[il]) *
	sqrt( pv_data[il] + pow(e_pv_rare[il],2) +
	      pow(pv_sf[il][isr],2)*( pow(e_pv_ttsl[il],2) + pow(e_pv_wjets[il],2) + 
				 pow(e_pv_st[il],2) + pow(e_pv_ttdl[il],2) ) );
      e_pv_sf_rel[il] = e_pv_sf[il][isr]/pv_sf[il][isr];//fractional error
      if (doverbose) 
	cout<<"DEFAULT PREVETO "<<pv_sf[il][isr]<<" pm "<< e_pv_sf[il][isr]<<endl;
    }

    //PreVeto ScaleFactor for alternative sample
    float e_pv_sf_alt_rel[NLEP];
    float pv_sf_alt[NLEP];
    float e_pv_sf_alt[NLEP];
    float pv_ttdl_alt[NLEP];
    float pv_ttsl_alt[NLEP];    
    float e_pv_ttdl_alt[NLEP];
    float e_pv_ttsl_alt[NLEP];    
    for (int il = 0; il<NLEP; il++) {
      pv_ttdl_alt[il] = h_pv_mc_alt[TTDL][il]->GetBinContent(PEAK);
      e_pv_ttdl_alt[il] = h_pv_mc_alt[TTDL][il]->GetBinError(PEAK);
      pv_ttsl_alt[il] = h_pv_mc_alt[TTSL][il]->GetBinContent(PEAK);
      e_pv_ttsl_alt[il] = h_pv_mc_alt[TTSL][il]->GetBinError(PEAK);
      pv_sf_alt[il] = (pv_data[il]-pv_rare[il])/
	(pv_ttsl_alt[il] + pv_wjets[il] + pv_st[il] + pv_ttdl_alt[il]);
      e_pv_sf_alt[il] = 1./(pv_ttsl_alt[il] + pv_wjets[il] + pv_st[il] + pv_ttdl_alt[il]) *
    	sqrt( pv_data[il] + pow(e_pv_rare[il],2) +
    	      pow(pv_sf_alt[il],2)*( pow(e_pv_ttsl_alt[il],2) + pow(e_pv_wjets[il],2) + 
    					  pow(e_pv_st[il],2) + pow(e_pv_ttdl_alt[il],2) ) );
      e_pv_sf_alt_rel[il] = e_pv_sf_alt[il]/pv_sf_alt[il];//fractional error
      if (doverbose) 
	cout<<"ALTERNATIVE PREVETO "<<pv_sf_alt[il]<<" pm "<< e_pv_sf_alt[il]<<endl;
    }

    //PostVeto ScaleFactor
    //in the uncertainty we have neglected the error on the PreVeto SF - think it is OK
    float e_sf_rel[NLEP];
    for (int il = 0; il<NLEP; il++) {
      sf[il][isr] = (data[il]-rare[il]-pv_sf[il][isr]*ttdl[il])/(ttsl[il]+wjets[il]+st[il] );
      e_sf[il][isr] = 1./(ttsl[il]+wjets[il]+st[il]) *
	sqrt( data[il]+pow(e_rare[il],2)+pow(pv_sf[il][isr]*e_ttdl[il],2) +
	      pow(sf[il][isr],2)*( pow(e_ttsl[il],2)+pow(e_wjets[il],2)+pow(e_st[il],2) ) );
      e_sf_rel[il] = e_sf[il][isr]/sf[il][isr];//fractional error
      if (doverbose) cout<<"POSTVETO "<<sf[il][isr]<<" pm "<<e_sf[il][isr]<<endl;
    }

    //BG prediction
    //tt and single top in l+jets	
    float pred_total[NLEP];
    float e_pred_total[NLEP];
    for (int il = 0; il<NLEP; il++) {
      pred_top1l[il][isr] = sf[il][isr]*ttp_tt[il]*(ttsl[il]+st[il]);
      //sig/ctrl dominates but include scale factor unc because it matters in low stat regions
      e_pred_top1l[il][isr] = pred_top1l[il][isr]*sqrt( pow(ttp_tt_runc[il],2) + pow(e_sf_rel[il],2) );
      pred_wjets[il][isr] = sf[il][isr]*ttp_wjets_corr[il]*wjets[il];
      e_pred_wjets[il][isr] = pred_wjets[il][isr]>0. ? 
	pred_wjets[il][isr]*sqrt( pow(ttp_wjets_runc[il],2) + 
				  pow(e_sf_rel[il],2) +
				  pow( (e_ttp_wjets_comb[isr]/ttp_wjets_comb[isr]), 2) )
	: 0.4;
      pred_rare[il][isr] = mc_rare[il];
      e_pred_rare[il][isr] = pred_rare[il][isr]*xs_unc_rare;
      pred_ttdl[il][isr] = pv_sf[il][isr]*mc_dl[il];
      e_pred_ttdl[il][isr] = pred_ttdl[il][isr] * 
	sqrt( pow(e_kfact,2) +
	      pow((e_2ndlepveto*frac_2ndlepveto),2) +
	      pow((e_2ndtauveto*frac_2ndtauveto),2) +
	      // pow((e_2ndlepveto*frac_dl_lepveto[isr]),2) +
	      // pow((e_2ndtauveto*frac_dl_tauveto[isr]),2) +
	      pow(e_dl_vars[isr],2) +
	      pow((e_pv_sf[il][isr]/pv_sf[il][isr]),2) +
	      pow((e_mc_dl[il]/mc_dl[il]),2) );
      //alternative ttdl 
      pred_ttdl_alt[il][isr] = pv_sf_alt[il]*mc_dl_alt[il];
      e_pred_ttdl_alt[il][isr] = pv_sf_alt[il]*e_mc_dl_alt[il];
      //Missing the Wjets systematic and a piece of the rare MC
      pred_total[il] = pred_top1l[il][isr]+pred_wjets[il][isr]+pred_rare[il][isr]+pred_ttdl[il][isr];
      e_pred_total[il] = sqrt( pow(e_pred_top1l[il][isr],2)+
			       pow(e_pred_wjets[il][isr],2)+
			       pow(e_pred_rare[il][isr],2)+
			       pow(e_pred_ttdl[il][isr],2) );
    }

    //cross section variations in the ScaleFactors
    float pv_sf_raredw[NLEP];
    float pv_sf_rareup[NLEP];
    float pv_sf_wjetsdw[NLEP];
    float pv_sf_wjetsup[NLEP];
    float sf_raredw[NLEP];
    float sf_rareup[NLEP];
    float sf_wjetsdw[NLEP];
    float sf_wjetsup[NLEP];
    for (int il = 0; il<NLEP; il++) {
      pv_sf_raredw[il] = 
	(pv_data[il]-(1.-xs_unc_rare)*pv_rare[il])/(pv_ttsl[il]+pv_wjets[il]+pv_st[il]+pv_ttdl[il]);
      pv_sf_rareup[il] = 
	(pv_data[il]-(1.+xs_unc_rare)*pv_rare[il])/(pv_ttsl[il]+pv_wjets[il]+pv_st[il]+pv_ttdl[il]);
      pv_sf_wjetsdw[il] = 
	(pv_data[il]-pv_rare[il])/(pv_ttsl[il]+(1.-xs_unc_wjets)*pv_wjets[il]+pv_st[il]+pv_ttdl[il]);
      pv_sf_wjetsup[il] = 
	(pv_data[il]-pv_rare[il])/(pv_ttsl[il]+(1.+xs_unc_wjets)*pv_wjets[il]+pv_st[il]+pv_ttdl[il]);
      sf_raredw[il] = 
	(data[il]-(1.-xs_unc_rare)*rare[il]-pv_sf[il][isr]*ttdl[il])/(ttsl[il]+wjets[il]+st[il] );
      sf_rareup[il] = 
	(data[il]-(1.+xs_unc_rare)*rare[il]-pv_sf[il][isr]*ttdl[il])/(ttsl[il]+wjets[il]+st[il] );
      sf_wjetsdw[il] = 
	(data[il]-rare[il]-pv_sf[il][isr]*ttdl[il])/(ttsl[il]+(1.-xs_unc_wjets)*wjets[il]+st[il] );
      sf_wjetsup[il] = 
	(data[il]-rare[il]-pv_sf[il][isr]*ttdl[il])/(ttsl[il]+(1.+xs_unc_wjets)*wjets[il]+st[il] );
    }

    //rare down
    float pred_top1l_raredw[NLEP];
    float pred_wjets_raredw[NLEP];
    float pred_ttdl_raredw[NLEP];
    float pred_total_raredw_norare[NLEP];
    float pred_rare_raredw[NLEP];
    float pred_total_raredw[NLEP];
    for (int il = 0; il<NLEP; il++) {
      pred_top1l_raredw[il] = sf_raredw[il]*ttp_tt[il]*(ttsl[il]+st[il]);
      pred_wjets_raredw[il] = sf_raredw[il]*ttp_wjets_corr[il]*wjets[il];
      pred_ttdl_raredw[il] = pv_sf_raredw[il]*mc_dl[il];
      //total without the rare component
      pred_total_raredw_norare[il] = pred_top1l_raredw[il]+pred_wjets_raredw[il]+pred_ttdl_raredw[il];
      pred_rare_raredw[il] = (1.-xs_unc_rare)*mc_rare[il];
      pred_total_raredw[il] = pred_total_raredw_norare[il]+pred_rare_raredw[il];
    }

    //rare up
    float pred_top1l_rareup[NLEP];
    float pred_wjets_rareup[NLEP];
    float pred_ttdl_rareup[NLEP];
    float pred_total_rareup_norare[NLEP];
    float pred_rare_rareup[NLEP];
    float pred_total_rareup[NLEP];
    for (int il = 0; il<NLEP; il++) {
      pred_top1l_rareup[il] = sf_rareup[il]*ttp_tt[il]*(ttsl[il]+st[il]);
      pred_wjets_rareup[il] = sf_rareup[il]*ttp_wjets_corr[il]*wjets[il];
      pred_ttdl_rareup[il] = pv_sf_rareup[il]*mc_dl[il];
      //total without the rare component
      pred_total_rareup_norare[il] = 
	pred_top1l_rareup[il]+pred_wjets_rareup[il]+pred_ttdl_rareup[il];
      pred_rare_rareup[il] = (1.+xs_unc_rare)*mc_rare[il];
      pred_total_rareup[il] =pred_total_rareup_norare[il]+pred_rare_rareup[il];
    }
    //wjets down
    float pred_top1l_wjetsdw[NLEP];
    float pred_wjets_wjetsdw[NLEP];
    float pred_ttdl_wjetsdw[NLEP];
    float pred_total_wjetsdw_norare[NLEP];
    float pred_rare_wjetsdw[NLEP];
    float pred_total_wjetsdw[NLEP];
    for (int il = 0; il<NLEP; il++) {
      pred_top1l_wjetsdw[il] = sf_wjetsdw[il]*ttp_tt[il]*(ttsl[il]+st[il]);
      pred_wjets_wjetsdw[il] = sf_wjetsdw[il]*ttp_wjets_corr[il]*(1.-xs_unc_wjets)*wjets[il];
      pred_ttdl_wjetsdw[il] = pv_sf_wjetsdw[il]*mc_dl[il];
      //total without the wjets component
      pred_total_wjetsdw_norare[il] = 
	pred_top1l_wjetsdw[il]+pred_wjets_wjetsdw[il]+pred_ttdl_wjetsdw[il];
      pred_rare_wjetsdw[il] = mc_rare[il];
      pred_total_wjetsdw[il] = pred_total_wjetsdw_norare[il]+pred_rare_wjetsdw[il];
    }
    //wjets up
    float pred_top1l_wjetsup[NLEP];
    float pred_wjets_wjetsup[NLEP];
    float pred_ttdl_wjetsup[NLEP];
    float pred_total_wjetsup_norare[NLEP];
    float pred_rare_wjetsup[NLEP];
    float pred_total_wjetsup[NLEP];
    for (int il = 0; il<NLEP; il++) {
      pred_top1l_wjetsup[il] = sf_wjetsup[il]*ttp_tt[il]*(ttsl[il]+st[il]);
      pred_wjets_wjetsup[il] = sf_wjetsup[il]*ttp_wjets_corr[il]*(1.+xs_unc_wjets)*wjets[il];
      pred_ttdl_wjetsup[il] = pv_sf_wjetsup[il]*mc_dl[il];
      //total without the wjets component
      pred_total_wjetsup_norare[il] = 
	pred_top1l_wjetsup[il]+pred_wjets_wjetsup[il]+pred_ttdl_wjetsup[il];
      pred_rare_wjetsup[il] = mc_rare[il];
      pred_total_wjetsup[il] = pred_total_wjetsup_norare[il]+pred_rare_wjetsup[il];
    }

    //systematic uncertainties on prediction
    float pred_wjets_syst[NLEP];
    float pred_rare_syst[NLEP];
    for (int il = 0; il<NLEP; il++) {
      pred_wjets_syst[il] = fabs(pred_total_wjetsup_norare[il]-pred_total_wjetsdw_norare[il])/2.;
      pred_rare_syst[il] = (pred_total_rareup_norare[il]-pred_total_raredw_norare[il])/2.;
      if (doverbose) {
	cout<<pred_wjets_syst[il]<<endl;
	cout<<pred_rare_syst[il]<<endl;
      }
    }
    //FINAL ANSWER BG prediction
    for (int il = 0; il<NLEP; il++) {
      pred_total_final[il][isr] = pred_total[il];
      e_pred_total_final[il][isr] = sqrt( pow(e_pred_top1l[il][isr],2)+
					  pow(e_pred_wjets[il][isr],2)+
					  pow((e_pred_rare[il][isr]+pred_rare_syst[il]),2) +
					  pow(e_pred_ttdl[il][isr],2) +
					  pow(pred_wjets_syst[il],2) );
      if (doverbose) cout<<pred_total_final[il][isr]<<" "<<e_pred_total_final[il][isr]<<endl;
    }


    //tt and single top in l+jets	
    pred_top1l_comb[isr] = pred_top1l[MUO][isr] + pred_top1l[ELE][isr];
    e_pred_top1l_comb[isr] = sqrt( pow( ((pred_top1l[ELE][isr]*ttp_tt_runc[ELE])+
    				  (pred_top1l[MUO][isr]*ttp_tt_runc[MUO])), 2) +
    			    pow( (pred_top1l[ELE][isr]*e_sf_rel[ELE]), 2) + 
    			    pow( (pred_top1l[MUO][isr]*e_sf_rel[MUO]), 2) ); 
    pred_wjets_comb[isr] = pred_wjets[MUO][isr] + pred_wjets[ELE][isr];
    e_pred_wjets_comb[isr] = pred_wjets_comb[isr]>0. ? 
      sqrt( pow( ((pred_wjets[ELE][isr]*ttp_wjets_runc[ELE])+
		  (pred_wjets[MUO][isr]*ttp_wjets_runc[MUO])), 2) +
	    pow( (pred_wjets[ELE][isr]*e_sf_rel[ELE]), 2) + 
	    pow( (pred_wjets[MUO][isr]*e_sf_rel[MUO]), 2) +
	    pow( (pred_wjets_comb[isr]*e_ttp_wjets_comb[isr]/ttp_wjets_comb[isr]), 2) ) 
      : 0.6; 
    pred_rare_comb[isr] = pred_rare[MUO][isr] + pred_rare[ELE][isr];
    e_pred_rare_comb[isr] = e_pred_rare[MUO][isr] + e_pred_rare[ELE][isr];
    pred_ttdl_comb[isr] = pred_ttdl[MUO][isr] + pred_ttdl[ELE][isr];
    float frac_ttdl[NLEP];
    frac_ttdl[MUO] = pred_ttdl[MUO][isr]/(pred_ttdl_comb[isr]);
    frac_ttdl[ELE] = 1.-frac_ttdl[MUO];
    e_pred_ttdl_comb[isr] = pred_ttdl_comb[isr] * 
      sqrt( pow(e_kfact,2) +
	    pow((e_2ndlepveto*frac_2ndlepveto),2) +
	    pow((e_2ndtauveto*frac_2ndtauveto),2) +
	    // pow((e_2ndlepveto*frac_dl_lepveto[isr]),2) +
	    // pow((e_2ndtauveto*frac_dl_tauveto[isr]),2) +
	    pow(e_dl_vars[isr],2) +
	    pow((frac_ttdl[MUO]*e_pv_sf[MUO][isr]/pv_sf[MUO][isr]),2) +
	    pow((frac_ttdl[ELE]*e_pv_sf[ELE][isr]/pv_sf[ELE][isr]),2) +
	    pow( (addSqr(e_mc_dl[MUO],e_mc_dl[ELE])/(mc_dl[MUO]+mc_dl[ELE])),2) );
    e_mc_pred_ttdl_comb[isr] = pred_ttdl_comb[isr] * 
      addSqr(e_mc_dl[MUO],e_mc_dl[ELE])/(mc_dl[MUO]+mc_dl[ELE]);

    //Missing the Wjets systematic and a piece of the rare MC
    float pred_total_comb = pred_top1l_comb[isr]+
      pred_wjets_comb[isr]+pred_rare_comb[isr]+pred_ttdl_comb[isr];
    float e_pred_total_comb = sqrt( pow(e_pred_top1l_comb[isr],2)+
				    pow(e_pred_wjets_comb[isr],2)+
				    pow(e_pred_rare_comb[isr],2)+
				    pow(e_pred_ttdl_comb[isr],2) );
    if (doverbose) cout<<"total pred "<<pred_total_comb<<" pm "<<e_pred_total_comb<<endl;

    //rare down
    //total without the rare component
    float pred_total_raredw_norare_comb = 
      pred_total_raredw_norare[MUO] + pred_total_raredw_norare[ELE];

    //rare up
    //total without the rare component
    float pred_total_rareup_norare_comb = 
      pred_total_rareup_norare[MUO] + pred_total_rareup_norare[ELE];

    //wjets down
    //total without the rare component
    float pred_total_wjetsdw_norare_comb = 
      pred_total_wjetsdw_norare[MUO] + pred_total_wjetsdw_norare[ELE];

    //wjets up
    //total without the rare component
    float pred_total_wjetsup_norare_comb = 
      pred_total_wjetsup_norare[MUO] + pred_total_wjetsup_norare[ELE];

    //systematic uncertainties on prediction
    float pred_wjets_syst_comb = fabs(pred_total_wjetsup_norare_comb-pred_total_wjetsdw_norare_comb)/2.;
    float pred_rare_syst_comb = (pred_total_rareup_norare_comb-pred_total_raredw_norare_comb)/2.;
    if (doverbose) {
      cout<<pred_wjets_syst_comb<<endl;
      cout<<pred_rare_syst_comb<<endl;
    }
    //FINAL ANSWER BG prediction (COMBINATION)
    pred_total_final_comb[isr] = pred_total_comb;
    e_pred_total_final_comb[isr] = sqrt( pow(e_pred_top1l_comb[isr],2)+
					 pow(e_pred_wjets_comb[isr],2)+
					 pow((e_pred_rare_comb[isr]+pred_rare_syst_comb),2) +
					 pow(e_pred_ttdl_comb[isr],2) +
					 pow(pred_wjets_syst_comb,2) );
    if (doverbose) cout<<pred_total_final_comb[isr]<<" "<<e_pred_total_final_comb[isr]<<endl;

  }


  //dump tables with information
  cout<<endl;
  cout<<endl;
  cout<<endl;

  cout<<"-------------------------------------------------"<<endl;
  cout<<"*************************************************"<<endl;
  cout<<"-------------------------------------------------"<<endl;
  // cout<<"PRINT SUMMARY YIELDS TABLE"<<endl;
  // cout<<"-------------------------------------------------"<<endl;

  cout<<endl;
  cout<<endl;

 cout<<"\\hline"<<endl;
  printf("Sample \t ");
  for (int isr=0; isr<NSR; ++isr) {
    if (isr==NSR-1) printf("\t & %s \\\\ \n",samplename[isr].c_str());
    else printf("\t & %s",samplename[isr].c_str());
  }
  cout<<"\\hline"<<endl;
  cout<<"\\hline"<<endl;
  printf("Frac. \\ttdl\\ with true iso. trk. \t ");
  for (int isr=0; isr<NSR; ++isr) {
    printf(" & $%.2f \\pm %.2f$",frac_dl_lepveto[isr], e_frac_dl_lepveto[isr]);
  }
  printf(" \\\\\n");
  printf("Frac. \\ttdl\\ with true tau.\t ");
  for (int isr=0; isr<NSR; ++isr) {
    printf(" & $%.2f \\pm %.2f$",frac_dl_tauveto[isr], e_frac_dl_tauveto[isr]);
  }
  printf(" \\\\\n");
 cout<<"\\hline"<<endl;
  printf("Frac. \\ttdl\\ out of acc. \t ");
  for (int isr=0; isr<NSR; ++isr) {
    printf(" & $%.2f \\pm %.2f$",frac_dl_outacc[isr], e_frac_dl_outacc[isr]);
  }
  printf(" \\\\\n");
 cout<<"\\hline"<<endl;
  printf("Pred. \\ttdl\\ out of acc. \t ");
  for (int isr=0; isr<NSR; ++isr) {
    printf(" & $%.2f$",frac_dl_outacc[isr]*pred_ttdl_comb[isr]);
  }
  printf(" \\\\\n");
 cout<<"\\hline"<<endl;

  cout<<endl;
  cout<<endl;
  cout<<endl;

 cout<<"\\hline"<<endl;
  printf("Sample \t ");
  for (int isr=0; isr<NSR; ++isr) {
    if (isr==NSR-1) printf("\t & %s \\\\ \n",samplename[isr].c_str());
    else printf("\t & %s",samplename[isr].c_str());
  }
  cout<<"\\hline"<<endl;
  cout<<"\\hline"<<endl;
  printf("Frac. \\ttdl\\ with true iso. trk. \t ");
  for (int isr=0; isr<NSR; ++isr) {
    printf(" & $%.2f \\pm %.2f$",frac_dl_lepveto_alt[isr], e_frac_dl_lepveto_alt[isr]);
  }
  printf(" \\\\\n");
  printf("Frac. \\ttdl\\ with true tau.\t ");
  for (int isr=0; isr<NSR; ++isr) {
    printf(" & $%.2f \\pm %.2f$",frac_dl_tauveto_alt[isr], e_frac_dl_tauveto_alt[isr]);
  }
  printf(" \\\\\n");
 cout<<"\\hline"<<endl;
  printf("Frac. \\ttdl\\ out of acc. \t ");
  for (int isr=0; isr<NSR; ++isr) {
    printf(" & $%.2f \\pm %.2f$",frac_dl_outacc_alt[isr], e_frac_dl_outacc_alt[isr]);
  }
  printf(" \\\\\n");
 cout<<"\\hline"<<endl;
  printf("Pred. \\ttdl\\ out of acc. \t ");
  for (int isr=0; isr<NSR; ++isr) {
    printf(" & $%.2f$",frac_dl_outacc_alt[isr]*(pred_ttdl_alt[MUO][isr]+pred_ttdl_alt[ELE][isr]));
  }
  printf(" \\\\\n");
  cout<<"\\hline"<<endl;
  cout<<"\\hline"<<endl;
  
  cout<<endl;
  cout<<endl;
  cout<<endl;
  cout<<"\\hline"<<endl;
  printf("Powheg \t ");
  for (int isr=0; isr<NSR; ++isr) {
    printf(" & $%.1f \\pm %.1f$",pred_ttdl_comb[isr], e_mc_pred_ttdl_comb[isr]);
  }
  printf(" \\\\\n");
  printf("Madgraph \t ");
  for (int isr=0; isr<NSR; ++isr) {
    float pred_ttdl_alt_comb = pred_ttdl_alt[MUO][isr]+pred_ttdl_alt[ELE][isr];
    float e_pred_ttdl_alt_comb = addSqr(e_pred_ttdl_alt[MUO][isr], e_pred_ttdl_alt[ELE][isr]);
    printf(" & $%.1f \\pm %.1f$",pred_ttdl_alt_comb, e_pred_ttdl_alt_comb);
  }
  printf(" \\\\\n");
  cout<<"\\hline"<<endl;
  printf("Difference (rel.) \t ");
  for (int isr=0; isr<NSR; ++isr) {
    float absdiff = fabs(pred_ttdl_comb[isr]-pred_ttdl_alt[MUO][isr]-pred_ttdl_alt[ELE][isr]);
    float fracdiff = absdiff/pred_ttdl_comb[isr];
    printf(" & $%.1f$ ($%.0f$\\%%)", absdiff, fracdiff*100.);
  }
  printf(" \\\\\n");
  cout<<"\\hline"<<endl;
  
  cout<<endl;
  cout<<endl;
  cout<<endl;


  cout<<"\\hline"<<endl;
  printf("Sample \t ");
  for (int isr=0; isr<NSR; ++isr) {
    if (isr==NSR-1) printf("\t & %s \\\\ \n",samplename[isr].c_str());
    else printf("\t & %s",samplename[isr].c_str());
  }
  for (int il=0; il<2; ++il) {
    //    printf("\\multicolumn{%i}{c}{%s} \\\\\n", NSR+1, il==MUO ? "Muon": "Electron");
    cout<<"\\hline"<<endl;
    cout<<"\\hline"<<endl;
    printf("%s pre-veto \\mt-SF \t ", il==MUO ? "$\\mu$": "e");
    for (int isr=0; isr<NSR; ++isr) {
      printf(" & $%.2f \\pm %.2f$",pv_sf[il][isr], e_pv_sf[il][isr]);
    }
    printf(" \\\\\n");
    printf("%s post-veto \\mt-SF \t ", il==MUO ? "$\\mu$": "e");
    for (int isr=0; isr<NSR; ++isr) {
      printf(" & $%.2f \\pm %.2f$",sf[il][isr], e_sf[il][isr]);
    }
    printf(" \\\\\n");
    cout<<"\\hline"<<endl;
    printf("%s veto \\mt-SF \t ", il==MUO ? "$\\mu$": "e");
    for (int isr=0; isr<NSR; ++isr) {
      printf(" & $%.2f$ ($%.2f$)",sf[il][isr]/pv_sf[il][isr], fabs(1.-(sf[il][isr]/pv_sf[il][isr])));
    }
    printf(" \\\\\n");

  }

  cout<<endl;
  cout<<endl;
  cout<<"\\hline"<<endl;
  cout<<"\\hline"<<endl;
  printf("%s \t ","$R_{wjet}$");
  for (int isr=0; isr<NSR; ++isr) 
    printf(" & $%.3f \\pm %.3f$ ", ttp_wjets_comb[isr], e_ttp_wjets_comb[isr]);
  printf(" \\\\\n");
  cout<<"\\hline"<<endl;
  printf("%s \t ","$R_{top}$");
  for (int isr=0; isr<NSR; ++isr) 
    printf(" & $%.3f \\pm %.3f$ ", ttp_ttsl_comb[isr], e_ttp_ttsl_comb[isr]);
  printf(" \\\\\n");

  cout<<"\\hline"<<endl;
  cout<<"\\hline"<<endl;
  printf("%s \t ","$SF_{wjet}$");
  for (int isr=0; isr<NSR; ++isr) 
    printf(" & $%.3f \\pm %.3f$ ", ttp_wj_sf_constantabove50pcunc[isr], ttp_wjets_aunc_constantabove50pcunc[isr]);
  printf(" \\\\\n");
  cout<<"\\hline"<<endl;
  cout<<"\\hline"<<endl;
  
  printf("%s \t ","$SF_{top}$ (from decomp. method)");
  for (int isr=0; isr<NSR; ++isr) 
    printf(" & $%.3f \\pm %.3f$ ", ttp_ttsl_sf[isr], ttp_ttsl_sf_aunc[isr]);
  printf(" \\\\\n");
  cout<<"\\hline"<<endl;
  printf("%s \t ","$SF_{top}*R_{top}$ (from decomp. method)");
  for (int isr=0; isr<NSR; ++isr) 
    printf(" & $%.3f \\pm %.3f$ ", ttp_tt_corr[MUO][isr], ttp_tt_sf_runc[MUO][isr]*ttp_tt_corr[MUO][isr]);
  printf(" \\\\\n");
  cout<<"\\hline"<<endl;
  printf("%s \t ","$SF_{top}*R_{top}$ (from pess-opt av.)");
  for (int isr=0; isr<NSR; ++isr) 
    printf(" & $%.3f \\pm %.3f$ ", ttp_tt_av[MUO][isr], ttp_tt_av_runc[MUO][isr]*ttp_tt_av[MUO][isr]);
  printf(" \\\\\n");
  cout<<"\\hline"<<endl;

  // for (int il=0; il<2; ++il) {
  //   //    printf("\\multicolumn{%i}{c}{%s} \\\\\n", NSR+1, il==MUO ? "Muon": "Electron");
  //   cout<<"\\hline"<<endl;
  //   cout<<"\\hline"<<endl;
  //   // printf("%s TTP SL Top %s \t ", il==MUO ? "$\\mu$": "e", 
  //   // 	   il==MUO ? "($R^{\\mu}\_{top}$)": "($R^e\_{top}$)");
  //   printf("%s \t ", il==MUO ? "$R^{\\mu}\_{top}$": "$R^e\_{top}$");
  //   for (int isr=0; isr<NSR; ++isr) {
  //     printf(" & $%.3f \\pm %.3f$ ",
  // 	     ttp_tt[isr][il], e_ttsl_ttp[isr][il]);
  //     // printf(" & %.3f * %.3f (%.2f %%) ",
  //     // 	     ttp_tt[isr][il],ettp_tt[isr][il],ettp_tt[isr][il]/ttp_tt[isr][il]*100.);
  //   }
  //   printf(" \\\\\n");
  //   // printf("%s TTP W+Jets %s \t ", il==MUO ? "$\\mu$": "e", 
  //   // 	   il==MUO ? "($R^{\\mu}\_{wjet}$)": "($R^e\_{wjet}$)");
  //   printf("%s \t ", il==MUO ? "$R^{\\mu}\_{wjet}$": "$R^e\_{wjet}$");
  //   for (int isr=0; isr<NSR; ++isr) {
  //     printf(" & $%.3f \\pm %.3f$ ",
  // 	     ttp_wj[isr][il], e_wjets_ttp[isr][il]);
  //     // printf(" & %.3f * %.3f (%.2f %%) ",
  //     // 	     ttp_wj[isr][il],ettp_wj[isr][il],ettp_wj[isr][il]/ttp_wj[isr][il]*100.);
  //   }
  //   printf(" \\\\\n");
  // }

  // for (int il=0; il<2; ++il) {
  //   //    printf("\\multicolumn{%i}{c}{%s} \\\\\n", NSR+1, il==MUO ? "Muon": "Electron");
  //   cout<<"\\hline"<<endl;
  //   cout<<"\\hline"<<endl;
  //   printf("%s Frac. \\ttdl\\ with true iso. trk. \t ", il==MUO ? "$\\mu$": "e");
  //   for (int isr=0; isr<NSR; ++isr) {
  //     printf("& $%.2f \\pm %.2f$ ",
  // 	     mplls_wtrk_frac[isr][il],emplls_wtrk_frac[isr][il]);
  //   }
  //   printf(" \\\\\n");
  // }


  //dump table with results
  cout<<endl;
  cout<<endl;
  cout<<endl;
  cout<<"\\hline"<<endl;
  printf("Sample \t ");
  for (int isr=0; isr<NSR; ++isr) {
    if (isr==NSR-1) printf("\t & %s \\\\ \n",samplename[isr].c_str());
    else printf("\t & %s",samplename[isr].c_str());
  }
  cout<<"\\hline"<<endl;
  cout<<"\\hline"<<endl;
  printf("\\multicolumn{%i}{c}{Muon} \\\\ \n", NSR+1);
  cout<<"\\hline"<<endl;
  printf("\\ttdl\\ \t ");
  for (int isr=0; isr<NSR; ++isr) {
    if (isr==NSR-1) printf("\t & $%.1f \\pm %.1f$ ($%.0f$\\%%) \\\\ \n",
		pred_ttdl[MUO][isr],e_pred_ttdl[MUO][isr],
		e_pred_ttdl[MUO][isr]/pred_ttdl[MUO][isr]*100.);
    else printf("\t & $%.1f \\pm %.1f$ ($%.0f$\\%%) ",
		pred_ttdl[MUO][isr],e_pred_ttdl[MUO][isr],
		e_pred_ttdl[MUO][isr]/pred_ttdl[MUO][isr]*100.);
  }
  printf("\\ttsl\\ \\& single top ($1\\ell$) \t ");
  for (int isr=0; isr<NSR; ++isr) {
    if (isr==NSR-1) printf("\t & $%.1f \\pm %.1f$ ($%.0f$\\%%) \\\\ \n",
			   pred_top1l[MUO][isr],e_pred_top1l[MUO][isr],
			   e_pred_top1l[MUO][isr]/pred_top1l[MUO][isr]*100.);
    else printf("\t & $%.1f \\pm %.1f$ ($%.0f$\\%%) ",
		pred_top1l[MUO][isr],e_pred_top1l[MUO][isr],
		e_pred_top1l[MUO][isr]/pred_top1l[MUO][isr]*100.);
  }
  printf("\\wjets\\  ");
  for (int isr=0; isr<NSR; ++isr) {
    if (isr==NSR-1) printf("\t & $%.1f \\pm %.1f$ ($%.0f$\\%%) \\\\ \n",
			   pred_wjets[MUO][isr],e_pred_wjets[MUO][isr],
			   pred_wjets[MUO][isr]>0. ? e_pred_wjets[MUO][isr]/pred_wjets[MUO][isr]*100. : 100.);
    else printf("\t & $%.1f \\pm %.1f$ ($%.0f$\\%%) ",
		pred_wjets[MUO][isr],e_pred_wjets[MUO][isr],
		pred_wjets[MUO][isr]>0. ? e_pred_wjets[MUO][isr]/pred_wjets[MUO][isr]*100. : 100.);
  }
  printf("Rare \t ");
  for (int isr=0; isr<NSR; ++isr) {
    if (isr==NSR-1) printf("\t & $%.1f \\pm %.1f$ ($%.0f$\\%%) \\\\ \n",
			   pred_rare[MUO][isr],e_pred_rare[MUO][isr],
			   e_pred_rare[MUO][isr]/pred_rare[MUO][isr]*100.);
    else printf("\t & $%.1f \\pm %.1f$ ($%.0f$\\%%) ",
		pred_rare[MUO][isr],e_pred_rare[MUO][isr],
		e_pred_rare[MUO][isr]/pred_rare[MUO][isr]*100.);
  }
  cout<<"\\hline"<<endl;
  printf("Total \t ");
  for (int isr=0; isr<NSR; ++isr) {
    if (isr==NSR-1) printf("\t & $%.1f \\pm %.1f$ ($%.0f$\\%%) \\\\ \n",
			   pred_total_final[MUO][isr],e_pred_total_final[MUO][isr],
			   e_pred_total_final[MUO][isr]/pred_total_final[MUO][isr]*100.);
    else printf("\t & $%.1f \\pm %.1f$ ($%.0f$\\%%) ",
		pred_total_final[MUO][isr],e_pred_total_final[MUO][isr],
		e_pred_total_final[MUO][isr]/pred_total_final[MUO][isr]*100.);
  }

  cout<<"\\hline"<<endl;
  cout<<"\\hline"<<endl;
  // printf("Data \t --- \n");
  // cout<<"\\hline"<<endl;
  // cout<<"\\hline"<<endl;
  cout<<"\\hline"<<endl;
  printf("\\multicolumn{%i}{c}{Electron} \\\\ \n", NSR+1);
  cout<<"\\hline"<<endl;
  printf("\\ttdl\\ \t ");
  for (int isr=0; isr<NSR; ++isr) {
    if (isr==NSR-1) printf("\t & $%.1f \\pm %.1f$ ($%.0f$\\%%) \\\\ \n",
		pred_ttdl[ELE][isr],e_pred_ttdl[ELE][isr],
		e_pred_ttdl[ELE][isr]/pred_ttdl[ELE][isr]*100.);
    else printf("\t & $%.1f \\pm %.1f$ ($%.0f$\\%%) ",
		pred_ttdl[ELE][isr],e_pred_ttdl[ELE][isr],
		e_pred_ttdl[ELE][isr]/pred_ttdl[ELE][isr]*100.);
  }
  printf("\\ttsl\\ \\& single top ($1\\ell$) \t ");
  for (int isr=0; isr<NSR; ++isr) {
    if (isr==NSR-1) printf("\t & $%.1f \\pm %.1f$ ($%.0f$\\%%) \\\\ \n",
			   pred_top1l[ELE][isr],e_pred_top1l[ELE][isr],
			   e_pred_top1l[ELE][isr]/pred_top1l[ELE][isr]*100.);
    else printf("\t & $%.1f \\pm %.1f$ ($%.0f$\\%%) ",
		pred_top1l[ELE][isr],e_pred_top1l[ELE][isr],
		e_pred_top1l[ELE][isr]/pred_top1l[ELE][isr]*100.);
  }
  printf("\\wjets\\  ");
  for (int isr=0; isr<NSR; ++isr) {
    if (isr==NSR-1) printf("\t & $%.1f \\pm %.1f$ ($%.0f$\\%%) \\\\ \n",
			   pred_wjets[ELE][isr],e_pred_wjets[ELE][isr],
			   pred_wjets[ELE][isr]>0. ? e_pred_wjets[ELE][isr]/pred_wjets[ELE][isr]*100. : 100.);
    else printf("\t & $%.1f \\pm %.1f$ ($%.0f$\\%%) ",
		pred_wjets[ELE][isr],e_pred_wjets[ELE][isr],
		pred_wjets[ELE][isr]>0. ? e_pred_wjets[ELE][isr]/pred_wjets[ELE][isr]*100. : 100.);
  }
  printf("Rare \t ");
  for (int isr=0; isr<NSR; ++isr) {
    if (isr==NSR-1) printf("\t & $%.1f \\pm %.1f$ ($%.0f$\\%%) \\\\ \n",
			   pred_rare[ELE][isr],e_pred_rare[ELE][isr],
			   e_pred_rare[ELE][isr]/pred_rare[ELE][isr]*100.);
    else printf("\t & $%.1f \\pm %.1f$ ($%.0f$\\%%) ",
		pred_rare[ELE][isr],e_pred_rare[ELE][isr],
		e_pred_rare[ELE][isr]/pred_rare[ELE][isr]*100.);
  }
  cout<<"\\hline"<<endl;
  printf("Total \t ");
  for (int isr=0; isr<NSR; ++isr) {
    if (isr==NSR-1) printf("\t & $%.1f \\pm %.1f$ ($%.0f$\\%%) \\\\ \n",
			   pred_total_final[ELE][isr],e_pred_total_final[ELE][isr],
			   e_pred_total_final[ELE][isr]/pred_total_final[ELE][isr]*100.);
    else printf("\t & $%.1f \\pm %.1f$ ($%.0f$\\%%) ",
		pred_total_final[ELE][isr],e_pred_total_final[ELE][isr],
		e_pred_total_final[ELE][isr]/pred_total_final[ELE][isr]*100.);
  }
  cout<<"\\hline"<<endl;
  cout<<"\\hline"<<endl;
  // printf("Data \t --- \n");
  // cout<<"\\hline"<<endl;
  // cout<<"\\hline"<<endl;
  cout<<"\\hline"<<endl;
  printf("\\multicolumn{%i}{c}{Muon+Electron Combined} \\\\ \n", NSR+1);
  cout<<"\\hline"<<endl;
  printf("\\ttdl\\ \t ");
  for (int isr=0; isr<NSR; ++isr) {
    if (isr==NSR-1) printf("\t & $%.1f \\pm %.1f$ ($%.0f$\\%%) \\\\ \n",
		pred_ttdl_comb[isr],e_pred_ttdl_comb[isr],
		e_pred_ttdl_comb[isr]/pred_ttdl_comb[isr]*100.);
    else printf("\t & $%.1f \\pm %.1f$ ($%.0f$\\%%) ",
		pred_ttdl_comb[isr],e_pred_ttdl_comb[isr],
		e_pred_ttdl_comb[isr]/pred_ttdl_comb[isr]*100.);
  }
  printf("\\ttsl\\ \\& single top ($1\\ell$) \t ");
  for (int isr=0; isr<NSR; ++isr) {
    if (isr==NSR-1) printf("\t & $%.1f \\pm %.1f$ ($%.0f$\\%%) \\\\ \n",
			   pred_top1l_comb[isr],e_pred_top1l_comb[isr],
			   e_pred_top1l_comb[isr]/pred_top1l_comb[isr]*100.);
    else printf("\t & $%.1f \\pm %.1f$ ($%.0f$\\%%) ",
		pred_top1l_comb[isr],e_pred_top1l_comb[isr],
		e_pred_top1l_comb[isr]/pred_top1l_comb[isr]*100.);
  }
  printf("\\wjets\\  ");
  for (int isr=0; isr<NSR; ++isr) {
    if (isr==NSR-1) printf("\t & $%.1f \\pm %.1f$ ($%.0f$\\%%) \\\\ \n",
			   pred_wjets_comb[isr],e_pred_wjets_comb[isr],
			   pred_wjets_comb[isr]>0. ? e_pred_wjets_comb[isr]/pred_wjets_comb[isr]*100. : 100.);
    else printf("\t & $%.1f \\pm %.1f$ ($%.0f$\\%%) ",
		pred_wjets_comb[isr],e_pred_wjets_comb[isr],
		pred_wjets_comb[isr]>0. ? e_pred_wjets_comb[isr]/pred_wjets_comb[isr]*100. : 100.);
  }
  printf("Rare \t ");
  for (int isr=0; isr<NSR; ++isr) {
    if (isr==NSR-1) printf("\t & $%.1f \\pm %.1f$ ($%.0f$\\%%) \\\\ \n",
			   pred_rare_comb[isr],e_pred_rare_comb[isr],
			   e_pred_rare_comb[isr]/pred_rare_comb[isr]*100.);
    else printf("\t & $%.1f \\pm %.1f$ ($%.0f$\\%%) ",
		pred_rare_comb[isr],e_pred_rare_comb[isr],
		e_pred_rare_comb[isr]/pred_rare_comb[isr]*100.);
  }
  cout<<"\\hline"<<endl;
  printf("Total \t ");
  for (int isr=0; isr<NSR; ++isr) {
    if (isr==NSR-1) printf("\t & $%.1f \\pm %.1f$ ($%.0f$\\%%) \\\\ \n",
			   pred_total_final_comb[isr],e_pred_total_final_comb[isr],
			   e_pred_total_final_comb[isr]/pred_total_final_comb[isr]*100.);
    else printf("\t & $%.1f \\pm %.1f$ ($%.0f$\\%%) ",
		pred_total_final_comb[isr],e_pred_total_final_comb[isr],
		e_pred_total_final_comb[isr]/pred_total_final_comb[isr]*100.);
  }
  cout<<"\\hline"<<endl;
  cout<<"\\hline"<<endl;
  // printf("Data \t --- \n");
  // cout<<"\\hline"<<endl;
  // cout<<"\\hline"<<endl;
  // cout<<"\\hline"<<endl;
  cout<<endl;
  cout<<endl;
  




  cout<<"\\hline"<<endl;
  printf("\\multicolumn{%i}{c}{Muon+Electron Combined} \\\\ \n", NSR+1);
  cout<<"\\hline"<<endl;
  printf("\\ttdl\\ \t ");
  for (int isr=0; isr<NSR; ++isr) {
    if (isr==NSR-1) printf("\t & $%.1f \\pm %.1f$ \\\\ \n",
			   pred_ttdl_comb[isr],e_pred_ttdl_comb[isr]);
    else printf("\t & $%.1f \\pm %.1f$ ",
		pred_ttdl_comb[isr],e_pred_ttdl_comb[isr]);
  }
  printf("\\ttsl\\ \\& single top ($1\\ell$) \t ");
  for (int isr=0; isr<NSR; ++isr) {
    if (isr==NSR-1) printf("\t & $%.1f \\pm %.1f$ \\\\ \n",
			   pred_top1l_comb[isr],e_pred_top1l_comb[isr]);
    else printf("\t & $%.1f \\pm %.1f$ ",
		pred_top1l_comb[isr],e_pred_top1l_comb[isr]);
  }
  printf("\\wjets\\  ");
  for (int isr=0; isr<NSR; ++isr) {
    if (isr==NSR-1) printf("\t & $%.1f \\pm %.1f$ \\\\ \n",
			   pred_wjets_comb[isr],e_pred_wjets_comb[isr]);
    else printf("\t & $%.1f \\pm %.1f$ ",
		pred_wjets_comb[isr],e_pred_wjets_comb[isr]);
  }
  printf("Rare \t ");
  for (int isr=0; isr<NSR; ++isr) {
    if (isr==NSR-1) printf("\t & $%.1f \\pm %.1f$ \\\\ \n",
			   pred_rare_comb[isr],e_pred_rare_comb[isr]);
    else printf("\t & $%.1f \\pm %.1f$ ",
		pred_rare_comb[isr],e_pred_rare_comb[isr]);
  }
  cout<<"\\hline"<<endl;
  printf("Total \t ");
  for (int isr=0; isr<NSR; ++isr) {
    if (isr==NSR-1) printf("\t & $%.1f \\pm %.1f$ ($%.0f$\\%%) \\\\ \n",
			   pred_total_final_comb[isr],e_pred_total_final_comb[isr],
			   e_pred_total_final_comb[isr]/pred_total_final_comb[isr]*100.);
    else printf("\t & $%.1f \\pm %.1f$ ($%.0f$\\%%) ",
		pred_total_final_comb[isr],e_pred_total_final_comb[isr],
		e_pred_total_final_comb[isr]/pred_total_final_comb[isr]*100.);
  }
  cout<<"\\hline"<<endl;
  cout<<"\\hline"<<endl;
  // printf("Data \t --- \n");
  // cout<<"\\hline"<<endl;
  // cout<<"\\hline"<<endl;
  // cout<<"\\hline"<<endl;
  cout<<endl;
  cout<<endl;

}
