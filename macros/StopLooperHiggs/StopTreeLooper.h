#ifndef STOPTREELOOPER_H
#define STOPTREELOOPER_H

#include "TChain.h"
#include "TFile.h"
#include "TString.h"
#include "TH1F.h"
#include "TH2F.h"

#include <iostream>
#include "Math/LorentzVector.h"
 
#include <cmath>
#include <map>

using namespace std;

class StopTreeLooper {

    public:
  typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > LorentzVector;

        StopTreeLooper();
        ~StopTreeLooper();
 
        void setOutFileName(string filename); 
        void loop(TChain *chain, TString name);

	//plotting
	void makeSIGPlots(float evtweight, std::map<std::string, TH1F*> &h_1d, 
			  string tag_selection, string flav_tag, float mtcut ); 
	void makeMinSIGPlots( float evtweight, std::map<std::string, TH1F*> &h_1d, 
			      string tag_selection, string flav_tag, 
			      float mtcut, float mbbval );
	//calculate mbb
	float getMbb(vector<LorentzVector> &bjets);

    private:

	string m_outfilename_;
	// njets requirement
	int min_njets;
	int min_nbjets;
	//for phi corrected met
	float t1metphicorr;
	float t1metphicorrphi;
	float t1metphicorrmt;
	//subtract lepton
	float t1metphicorr_lep;
	float t1metphicorrphi_lep;
	float t1metphicorrmt_lep;
	float dphi_pseudometlep;
	float leppt;
	//min dphi {j1,j2} 
	float dphimj12min;
	float dphimj1;
	float dphimj2;
	float dphimjmin;
	float dphimjmax;
	//b pt
	float pt_b;
	float dRleptB1;
	//mbb
	float mbb40;
	float mbb;
	//maria variables
	float htssl;
	float htosl;
	float htratiol;
	float htssm;
	float htosm;
	float htratiom;
	//for mt peak definition
	float min_mtpeak;
	float max_mtpeak; 
	//for mbb peak definition
	float min_mbb;
	float max_mbb; 
	//jets information
	int n_jets;
	int n_bjets40;
	int n_bjets30;
	int n_ljets;
	vector<LorentzVector> bjets40;
	vector<LorentzVector> bjets30;
	vector<float> btag;
	vector<float> sigma_jets;
	vector<int> mc;
	//systematics
	int n_bjets40_upBCShape;
	int n_bjets30_upBCShape;
	int n_ljets_upBCShape;
	vector<LorentzVector> bjets40_upBCShape;
	vector<LorentzVector> bjets30_upBCShape;
	int n_bjets40_downBCShape;
	int n_bjets30_downBCShape;
	int n_ljets_downBCShape;
	vector<LorentzVector> bjets40_downBCShape;
	vector<LorentzVector> bjets30_downBCShape;
	int n_bjets40_upLShape;
	int n_bjets30_upLShape;
	int n_ljets_upLShape;
	vector<LorentzVector> bjets40_upLShape;
	vector<LorentzVector> bjets30_upLShape;
	int n_bjets40_downLShape;
	int n_bjets30_downLShape;
	int n_ljets_downLShape;
	vector<LorentzVector> bjets40_downLShape;
	vector<LorentzVector> bjets30_downLShape;
	float mbb40_upBCShape;
	float mbb_upBCShape;
	float mbb40_downBCShape;
	float mbb_downBCShape;
	float mbb40_upLShape;
	float mbb_upLShape;
	float mbb40_downLShape;
	float mbb_downLShape;

	float chi2min;
	float chi2minprob;
	float mt2bmin;
	float mt2blmin;
	float mt2wmin;

	float pfcalo_metratio;
	float pfcalo_metdphi;
	
	//flags to blind signal region
	bool issigmbb;
	bool issigmbb40;
	bool issigmt120;
	bool issigmt150;


};

#endif
