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

	void getYield(float evtweight, string tag_selection, std::map<std::string, TH1F*> &h_1d, std::map<std::string, TH2F*> &h_2d, TString name, bool isData);
	void studyIsoTrack(float evtweight, string tag_selection, std::map<std::string, TH1F*> &h_1d, std::map<std::string, TH2F*> &h_2d, bool isData);
	void studyIsoTaus(float evtweight, string tag_selection, std::map<std::string, TH1F*> &h_1d, std::map<std::string, TH2F*> &h_2d, bool isData);
	void studyGenJet(float evtweight, string tag_selection, std::map<std::string, TH1F*> &h_1d,std::map<std::string, TH2F*> &h_2d, bool isData);
	void classify3B(float evtweight, string tag_selection, std::map<std::string, TH1F*> &h_1d, std::map<std::string, TH2F*> &h_2d, bool isData, bool doLRM);
	void studyAcceptance(float evtweight, string tag_selection, std::map<std::string, TH1F*> &h_1d);

	void fillSignalYieldTable(float evtweight, string tag_selection , std::map<std::string, TH1F*> &h_1d, std::map<std::string, TH1F*> &h_1d_qg, std::map<std::string, TH2F*> &h_2d_qg, bool isData);

	void studyMassbb(float evtweight, string tag_selection , std::map<std::string, TH1F*> &h_1d, std::map<std::string, TH2F*> &h_2d, bool isData) ;

	void plotCR4(float evtweight, string tag_selection , std::map<std::string, TH1F*> &h_1d, bool isData, bool doLRM, bool doTAU);
        void plotCR1(float evtweight, string tag_selection , std::map<std::string, TH1F*> &h_1d, bool isData);

        void setOutFileName(string filename); 
        void loop(TChain *chain, TString name);

    private:

	string m_outfilename_;
	//for phi corrected met
	float t1metphicorr;
	float t1metphicorrphi;
	float t1metphicorrmt;
	//for mt peak definition
	float min_mtpeak;
	float max_mtpeak; 

	vector<LorentzVector> jetsPUP4;
	vector<LorentzVector> jetsP4;
	vector<LorentzVector> jetsP4csv;
	vector<LorentzVector> genJetsP4;
        vector<float> csv;
        vector<float> mc3;
        vector<float> mc1;
        vector<float> sigma;
        vector<float> lrm;
        vector<float> chm;
        vector<float> neu;

	bool dataset_1l;
	bool dataset_CR4;
	float detaFB;
	int n_taujets;
	int n_ljets;

	float MHT;
	float MHTphi;

	float tightDiscr;
	float looseDiscr;

	float massHbb;
	float massHbbTau;

	float massHbb30;
	float massHbbTau30;

	float massHbbTight;
	float massHbbTauTight;

	vector<int> indexJetsP4; // this is for the 2l+3b  

	vector<int> indexLooseB40; // this is for the 2l+3b  
	vector<int> indexLooseB30; // this is for the 1l+4b and also the 2l +4b
	vector<int> indexTightB40; // this is for the 1l + 3b  

	vector<int> indexLooseB40Tau; //+3b  
	vector<int> indexLooseB30Tau; // +4b 
	vector<int> indexTightB40Tau; // this is for the 1l + 3b  (used for the 2l prediction)

	float minMassbbCut;
	float maxMassbbCut;

	float csvNonbMax;
	float csvNonbTauMax;

};

#endif
