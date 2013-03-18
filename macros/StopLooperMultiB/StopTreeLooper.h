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

	void studyIsoTrack(float evtweight, string tag_selection, std::map<std::string, TH1F*> &h_1d, std::map<std::string, TH2F*> &h_2d, bool isData);
	void studyIsoTaus(float evtweight, string tag_selection, std::map<std::string, TH1F*> &h_1d, std::map<std::string, TH2F*> &h_2d, bool isData);
	void classify3B(float evtweight, string tag_selection, std::map<std::string, TH1F*> &h_1d, std::map<std::string, TH2F*> &h_2d, bool isData, bool doLRM);

	void plotCR2(float evtweight, string tag_selection , std::map<std::string, TH1F*> &h_1d, std::map<std::string, TH2F*> &h_2d, bool isData);
        void plotCR4(float evtweight, string tag_selection , std::map<std::string, TH1F*> &h_1d, bool isData);
        void plotCR5(float evtweight, string tag_selection , std::map<std::string, TH1F*> &h_1d, bool isData);
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

	vector<LorentzVector> jetsP4;
        vector<float> csv;
        vector<float> mc3;
        vector<float> mc1;
        vector<float> sigma;
        vector<float> lrm;

};

#endif
