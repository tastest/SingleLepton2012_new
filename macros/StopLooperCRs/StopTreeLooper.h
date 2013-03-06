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
			   string tag_selection, string tag_kbin, string flav_tag, float mtcut ); 
	void makeCR1Plots(float evtweight, std::map<std::string, TH1F*> &h_1d, 
			   string tag_selection, string tag_njets, string tag_kbin, string flav_tag, float mtcut ); 
	void makeCR2Plots(float evtweight, std::map<std::string, TH1F*> &h_1d, 
			   string tag_selection, string tag_njets, string tag_kbin, string flav_tag_dl, float mtcut );
	void makeCR4Plots(float evtweight, std::map<std::string, TH1F*> &h_1d, 
			   string tag_selection, string tag_njets, string tag_kbin, string flav_tag_dl, float mtcut );
	void makeCR5Plots(float evtweight, std::map<std::string, TH1F*> &h_1d, 
			   string tag_selection, string tag_njets, string tag_kbin, string flav_tag_dl, float mtcut );
	void makeNJPlots( float evtweight, std::map<std::string, TH1F*> &h_1d, 
			   string tag_selection, string flav_tag ); 
	void makeZPlots(  float evtweight, std::map<std::string, TH1F*> &h_1d, 
			   string tag_selection, string tag_njets, string flav_tag );

    private:

	string m_outfilename_;
	// njets requirement
	int min_njets;
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
	//jets information
	int n_jets;
	int n_bjets;
	int n_ljets;
	vector<LorentzVector> jets;
	vector<float> btag;
	vector<float> sigma_jets;
	vector<int> mc;

	float chi2min_;
	float chi2minprob_;
	float mt2bmin_;
	float mt2blmin_;
	float mt2wmin_;

};

#endif
