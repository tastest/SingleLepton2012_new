#ifndef PARTONCOMBINATORICS_H
#define PARTONCOMBINATORICS_H

#include <list>
#include <vector>

class PartonCombinatorics {
	typedef struct {
		float one_chi2;
		float two_mt2b, two_mt2bl, two_mt2w;
		float three_mt2b, three_mt2bl, three_mt2w;
		float four_chi2b, four_chi2bl, four_chi2w;
	} MT2CHI2;

public:
	PartonCombinatorics(vector<LorentzVector> jets, vector<float> btag, vector<float> sigma_jets,
			vector<int> mc, LorentzVector lep, double met, double metphi, bool isData);

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

	list<Candidate> candidates_;
	list<Candidate> b_candidates_;
	MT2CHI2 mt2chi2_;

	static bool compare_in_chi2( Candidate &x, Candidate &y );
	static bool compare_in_mt2b( Candidate &x, Candidate &y );
	static bool compare_in_mt2bl( Candidate &x, Candidate &y );
	static bool compare_in_mt2w( Candidate &x, Candidate &y );
	static float min_with_value(list<Candidate> &candidates, float value, const char* fix, const char* var);

	void recoHadronicTop(StopTree*, bool);
	void applyBConsistency(list<Candidate> &candidates, StopTree* tree, float btagcut);
	void MT2CHI2Calculator(list<Candidate> &candidates, StopTree* tree);

	const float PTMIN_J1 = 30;
	const float PTMIN_J2 = 30;
	const float PTMIN_B = 30;
	const float PTMIN_O = 30;
	const float JET_PT = 30.;
	const float JET_ETA = 3.0;
	const float PDG_TOP_MASS = 173.5;
	const float PDG_W_MASS = 80.385;
	const float BTAG_MED = 0.679;
	const float BTAG_LOW = 0.244;

	const int N_JETS_TO_CONSIDER = 6;
	const int NUM_LEAD_JETS_0B = 3;
	const int NUM_LEAD_JETS_1B = 3;

	const bool __SORT = true;

};

#endif
