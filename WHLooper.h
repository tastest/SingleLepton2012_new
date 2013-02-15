#ifndef WHLOOPER_H
#define WHLOOPER_H

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

class WHLooper {

    public:
  typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > LorentzVector;

        WHLooper();
        ~WHLooper();

        void setOutFileName(string filename); 
        void loop(TChain *chain, TString name);

    private:

        void fillHists1D(std::map<std::string, TH1F*>& h_1d, const float evtweight = 1., const std::string& suffix = "");

	string m_outfilename_;
	//for phi corrected met
	float t1metphicorr;
	float t1metphicorrphi;
	float t1metphicorrmt;
	//for mt peak definition
	float min_mtpeak;
	float max_mtpeak; 

	std::vector<LorentzVector> * myPfJets;

};

#endif
