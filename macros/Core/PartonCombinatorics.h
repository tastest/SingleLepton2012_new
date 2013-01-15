#ifndef PARTONCOMBINATORICS_H
#define PARTONCOMBINATORICS_H

#include "../Core/stopUtils.h"



#include <list>
#include <vector>



        typedef struct {
                float one_chi2;
                float two_mt2b, two_mt2bl, two_mt2w;
                float three_mt2b, three_mt2bl, three_mt2w;
                float four_chi2b, four_chi2bl, four_chi2w;
        } MT2CHI2;

using namespace std;

class PartonCombinatorics {
	typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > LorentzVector;

public:
	PartonCombinatorics(vector<LorentzVector> jets, vector<float> btag, vector<float> sigma_jets,
			vector<float> mc, LorentzVector lep, float met, float metphi, bool isData);
	PartonCombinatorics(vector<LorentzVector> jets, vector<float> btag, vector<float> sigma_jets,
			vector<int> mc, LorentzVector lep, float met, float metphi, bool isData);

	MT2CHI2 getMt2Chi2();

	bool isData_;
	int n_jets;
	vector<LorentzVector> jets_;
	vector<float> btag_;
	vector<float> sigma_jets_;
	vector<int> mc_;
	LorentzVector lep_;
	double met_;
	double metphi_;
        int n_jets_;

	list<Candidate> candidates_;
	list<Candidate> b_candidates_;
	MT2CHI2 mt2chi2_;

	static bool compare_in_chi2( Candidate &x, Candidate &y );
	static bool compare_in_mt2b( Candidate &x, Candidate &y );
	static bool compare_in_mt2bl( Candidate &x, Candidate &y );
	static bool compare_in_mt2w( Candidate &x, Candidate &y );

	static double fc2 (double c1, double m12, double m22, double m02, bool verbose);
	static double fchi2 (double c1, double pt1, double sigma1, double pt2, double sigma2, double m12, double m22, double m02);
	static void minuitFunction(int&, double* , double &result, double par[], int);


	static float min_with_value(list<Candidate> &candidates, float value, const char* fix, const char* var);

	void recoHadronicTop();
	void applyBConsistency(float btagcut);
	void MT2CHI2Calculator();

	static const float PTMIN_J1 = 30;
	static const float PTMIN_J2 = 30;
	static const float PTMIN_B = 30;
	static const float PTMIN_O = 30;
	static const float JET_PT = 30.;
	static const float JET_ETA = 3.0;
	static const float PDG_TOP_MASS = 173.5;
	static const float PDG_W_MASS = 80.385;
	static const float BTAG_MED = 0.679;
	static const float BTAG_LOW = 0.244;

	static const int N_JETS_TO_CONSIDER = 6;
	static const int NUM_LEAD_JETS_0B = 3;
	static const int NUM_LEAD_JETS_1B = 3;

	static const bool __SORT = true;
	static const bool __debug = false;

};

#endif
