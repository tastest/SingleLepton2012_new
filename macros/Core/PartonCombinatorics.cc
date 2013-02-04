#include "PartonCombinatorics.h"

#include "TFitter.h"
#include "../Core/MT2Utility.h"
#include "../Core/mt2bl_bisect.h"
#include "../Core/mt2w_bisect.h"

PartonCombinatorics::PartonCombinatorics(vector<LorentzVector> jets, vector<float> btag, vector<float> sigma_jets,
		vector<float> mc, LorentzVector lep, float met, float metphi, bool isData){
        vector<int> new_mc;
        for (unsigned int i =0; i < mc.size(); ++i)
           new_mc.push_back((int) mc.at(i));

        PartonCombinatorics(jets, btag, sigma_jets, new_mc, lep, met, metphi, isData);
}

PartonCombinatorics::PartonCombinatorics(vector<LorentzVector> jets, vector<float> btag, vector<float> sigma_jets,
		vector<int> mc, LorentzVector lep, float met, float metphi, bool isData){
        if ( __debug ) cout << "PartonCombinatorics::Constructor " << endl;
	isData_ = isData;
	jets_ = jets;
	btag_ = btag;
	sigma_jets_ = sigma_jets;
	mc_ = mc;
	lep_ = lep;
	met_ = met;
	metphi_ = metphi;
        n_jets_ = jets.size();

        assert( jets_.size() == btag_.size() );
        assert( jets_.size() == mc_.size() );
        assert( jets_.size() == sigma_jets_.size() );

        if ( __debug ) cout << "PartonCombinatorics::Constructor isData = " << isData << endl;

	recoHadronicTop();
	applyBConsistency(BTAG_MED);
	MT2CHI2Calculator();
}

MT2CHI2 PartonCombinatorics::getMt2Chi2(){
	return mt2chi2_;
}

//--------------------------------------------------------------------
bool PartonCombinatorics::compare_in_chi2( Candidate &x, Candidate &y ){
  return x.chi2 < y.chi2;
}

//--------------------------------------------------------------------
bool PartonCombinatorics::compare_in_mt2b( Candidate &x, Candidate &y ){
  return x.mt2b < y.mt2b;
}

//--------------------------------------------------------------------
bool PartonCombinatorics::compare_in_mt2bl( Candidate &x, Candidate &y ){
  return x.mt2bl < y.mt2bl;
}

//--------------------------------------------------------------------
bool PartonCombinatorics::compare_in_mt2w( Candidate &x, Candidate &y ){
  return x.mt2w < y.mt2w;
}

float PartonCombinatorics::min_with_value(list<Candidate> &candidates, float value, const char* fix, const char* var){
  float min_value = 9999.0;

  for(list<Candidate>::iterator it = candidates.begin(); it != candidates.end(); it++){
	  if (strcmp(fix, "chi2") == 0 && it->chi2 == value ){
		     if ( strcmp(var, "mt2b")  == 0  && it->mt2b  < min_value ) min_value = it->mt2b;
		     if ( strcmp(var, "mt2bl") == 0  && it->mt2bl < min_value ) min_value = it->mt2bl;
		     if ( strcmp(var, "mt2w")  == 0  && it->mt2w  < min_value ) min_value = it->mt2w;
	  }
	  if (strcmp(fix, "mt2b") == 0 && it->mt2b == value ){
		     if ( strcmp(var, "chi2")  == 0  && it->chi2  < min_value ) min_value = it->chi2;
		     if ( strcmp(var, "mt2bl") == 0  && it->mt2bl < min_value ) min_value = it->mt2bl;
		     if ( strcmp(var, "mt2w")  == 0  && it->mt2w  < min_value ) min_value = it->mt2w;
	  }
	  if (strcmp(fix, "mt2bl") == 0 && it->mt2bl == value ){
		     if ( strcmp(var, "chi2")  == 0  && it->chi2  < min_value ) min_value = it->chi2;
		     if ( strcmp(var, "mt2b")  == 0  && it->mt2b  < min_value ) min_value = it->mt2b;
		     if ( strcmp(var, "mt2w")  == 0  && it->mt2w  < min_value ) min_value = it->mt2w;
	  }
	  if (strcmp(fix, "mt2w") == 0 && it->mt2w == value ){
		     if ( strcmp(var, "chi2")  == 0  && it->chi2  < min_value ) min_value = it->chi2;
		     if ( strcmp(var, "mt2b")  == 0  && it->mt2b  < min_value ) min_value = it->mt2b;
		     if ( strcmp(var, "mt2bl") == 0  && it->mt2bl < min_value ) min_value = it->mt2bl;
	  }
  }

  if ( min_value > 9998.0 ) return -0.999;
  return min_value;
}
//--------------------------------------------------------------------
double PartonCombinatorics::fc2 (double c1, double m12, double m22, double m02, bool verbose = false)
{
  if (verbose) {
    printf("c1: %4.2f\n", c1);
    printf("m12: %4.2f\n", m12);
    printf("m22: %4.2f\n", m22);
    printf("m02: %4.2f\n", m02);
  }

  double a = m22;
  double b = (m02 - m12 - m22) * c1;
  double c = m12 * c1 * c1 - PDG_W_MASS * PDG_W_MASS;

  if (verbose) {
    printf("a: %4.2f\n", a);
    printf("b: %4.2f\n", b);
    printf("c: %4.2f\n", c);
  }

  double num = -1. * b + sqrt(b * b - 4 * a * c);
  double den = 2 * a;

  if (verbose) {
    printf("num: %4.2f\n", num);
    printf("den: %4.2f\n", den);
    printf("num/den: %4.2f\n", num/den);
  }

  return (num/den);
}

//--------------------------------------------------------------------
double PartonCombinatorics::fchi2 (double c1, double pt1, double sigma1, double pt2, double sigma2,
              double m12, double m22, double m02){
  double rat1 = pt1 * (1 - c1) / sigma1;
  double rat2 = pt2 * (1 - fc2(c1, m12, m22, m02)) / sigma2;

  return ( rat1 * rat1 + rat2 * rat2);
}

//--------------------------------------------------------------------
void PartonCombinatorics::minuitFunction(int&, double* , double &result, double par[], int){
  result=fchi2(par[0], par[1], par[2], par[3], par[4], par[5], par[6], par[7]);
}

/* Reconstruct the hadronic top candidates, select the best candidate and
 * store the chi2 =  (m_jj - m_W)^2/sigma_m_jj + (m_jjj - m_t)^2/sigma_m_jjj
 * return the number of candidates found.
 *
 * n_jets - number of jets.
 * jets - jets
 * btag - b-tagging information of the jets
 * jets_ - qgjet montecarlo match number for the jets
 *
 * returns a list of candidates sorted by chi2 ( if __sort = true in .h );
 */
void PartonCombinatorics::recoHadronicTop(){
  if ( __debug ) cout << "PartonCombinatorics::recoHadronicTop " << endl;

  float metx = met_ * cos( metphi_ );
  float mety = met_ * sin( metphi_ );

  int ibl[5];
  int iw1[5];
  int iw2[5];
  int ib[5];

  if ( !isData_ ){
    // Matching MC algoritm search over all conbinations  until the
    // right combination is found. More than one candidate is suported
    //  but later only the first is used.

    int match = 0;
    for (int jbl=0; jbl<n_jets_; ++jbl )
      for (int jb=0; jb<n_jets_; ++jb )
        for (int jw1=0; jw1<n_jets_; ++jw1 )
          for (int jw2=jw1+1; jw2<n_jets_; ++jw2 )
            if ( (mc_.at(jw2)==2 && mc_.at(jw1)==2 && mc_.at(jb)==1 && mc_.at(jbl)==-1) ||
                 (mc_.at(jw2)==-2 && mc_.at(jw1)==-2 && mc_.at(jb)==-1 && mc_.at(jbl)==1) ||
		 (mc_.at(jw2)==5 && mc_.at(jw1)==5 && mc_.at(jb)==1 && mc_.at(jbl)==-1) ||
                 (mc_.at(jw2)==-5 && mc_.at(jw1)==-5 && mc_.at(jb)==-1 && mc_.at(jbl)==1) ) {
	      if ( match == 5 ) break;
              if ( __debug ) cout << "PartonCombinatorics::recoHadronicTop MC found:" << match << endl;
	      ibl[match] = jbl;
	      iw1[match] = jw1;
	      iw2[match] = jw2;
	      ib[match] = jb;
	      match++;
            }
  }

  if ( __debug ) cout << "PartonCombinatorics::recoHadronicTop MC done" << endl;
  ////////    * Combinatorics. j_1 Pt must be > PTMIN_W1 and so on.

  vector<int> v_i, v_j;
  vector<double> v_k1, v_k2;
  for ( int i=0; i<n_jets_; ++i )
    for ( int j=i+1; j<n_jets_; ++j ){
      double pt_w1 = jets_[i].Pt();
      double pt_w2 = jets_[j].Pt();
      if ( pt_w1 < PTMIN_J1 || fabs(jets_[i].Eta()) > JET_ETA ) continue;
      if ( pt_w2 < PTMIN_J2 || fabs(jets_[j].Eta()) > JET_ETA ) continue;

      //
      //  W
      //
      LorentzVector hadW = jets_[i] + jets_[j];

      //
      //  W Mass Constraint.
      //
      TFitter *minimizer = new TFitter();
      double p1 = -1;

      minimizer->ExecuteCommand("SET PRINTOUT", &p1, 1);
      minimizer->SetFCN(minuitFunction);
      minimizer->SetParameter(0 , "c1"     , 1.1             , 1 , 0 , 0);
      minimizer->SetParameter(1 , "pt1"    , 1.0             , 1 , 0 , 0);
      minimizer->SetParameter(2 , "sigma1" , sigma_jets_[i]   , 1 , 0 , 0);
      minimizer->SetParameter(3 , "pt2"    , 1.0             , 1 , 0 , 0);
      minimizer->SetParameter(4 , "sigma2" , sigma_jets_[j]   , 1 , 0 , 0);
      minimizer->SetParameter(5 , "m12"    , jets_[i].mass2() , 1 , 0 , 0);
      minimizer->SetParameter(6 , "m22"    , jets_[j].mass2() , 1 , 0 , 0);
      minimizer->SetParameter(7 , "m02"    , hadW.mass2()    , 1 , 0 , 0);

      for (unsigned int k = 1; k < 8; k++)
        minimizer->FixParameter(k);

      minimizer->ExecuteCommand("SIMPLEX", 0, 0);
      minimizer->ExecuteCommand("MIGRAD", 0, 0);

      double c1 = minimizer->GetParameter(0);
      if (c1!=c1) {
        cout<<"[PartonCombinatorics::recoHadronicTop] ERROR: c1 parameter is NAN! Skipping this parton combination"
	    <<endl;
        continue;
      }
      double c2 = fc2(c1, jets_[i].mass2(), jets_[j].mass2(), hadW.mass2());

      delete minimizer;


      //     * W Mass check :)
      //     *  Never trust a computer you can't throw out a window.
      //      *  - Steve Wozniak

      //      cout << "c1 = " <<  c1 << "  c1 = " << c2 << "   M_jj = "
      //           << ((jets_[i] * c1) + (jets_[j] * c2)).mass() << endl;

      v_i.push_back(i);
      v_j.push_back(j);
      v_k1.push_back(c1);
      v_k2.push_back(c2);
    }

  if ( __debug ) cout << "PartonCombinatorics::recoHadronicTop W done" << endl;

  list<Candidate> chi2candidates;

  mt2_bisect::mt2 mt2_event;
  mt2bl_bisect::mt2bl mt2bl_event;
  mt2w_bisect::mt2w mt2w_event;

  for ( int b=0; b<n_jets_; ++b )
    for (int o=0; o<n_jets_; ++o){

      if ( b == o ) continue;
      //apply pt and eta requirements
      double pt_b = jets_[b].Pt();
      double pt_o = jets_[o].Pt();
      if ( pt_b < PTMIN_B || fabs(jets_[b].Eta()) > JET_ETA ) continue;
      if ( pt_o < PTMIN_O || fabs(jets_[o].Eta()) > JET_ETA ) continue;

      ///
      //  MT2 Variables
      ///

      double pl[4];     // Visible lepton
      double pb1[4];    // bottom on the same side as the visible lepton
      double pb2[4];    // other bottom, paired with the invisible W
      double pmiss[3];  // <unused>, pmx, pmy   missing pT
      pl[0]= lep_.E(); pl[1]= lep_.Px(); pl[2]= lep_.Py(); pl[3]= lep_.Pz();
      pb1[1] = jets_[o].Px();  pb1[2] = jets_[o].Py();   pb1[3] = jets_[o].Pz();
      pb2[1] = jets_[b].Px();  pb2[2] = jets_[b].Py();   pb2[3] = jets_[b].Pz();
      pmiss[0] = 0.; pmiss[1] = metx; pmiss[2] = mety;

      double pmiss_lep[3];
      pmiss_lep[0] = 0.;
      pmiss_lep[1] = pmiss[1]+pl[1]; pmiss_lep[2] = pmiss[2]+pl[2];

      pb1[0] = jets_[o].mass();
      pb2[0] = jets_[b].mass();
      mt2_event.set_momenta( pb1, pb2, pmiss_lep );
      mt2_event.set_mn( 80.385 );   // Invisible particle mass
      double c_mt2b = mt2_event.get_mt2();

      pb1[0] = jets_[o].E();
      pb2[0] = jets_[b].E();
      mt2bl_event.set_momenta(pl, pb1, pb2, pmiss);
      double c_mt2bl = mt2bl_event.get_mt2bl();

      mt2w_event.set_momenta(pl, pb1, pb2, pmiss);
      double c_mt2w = mt2w_event.get_mt2w();

      //      cout << b << ":"<< btag[b] << " - " << o << ":" << btag[o] << " = " << c_mt2w << endl;

      for (unsigned int w = 0; w < v_i.size() ; ++w ){
        int i = v_i[w];
        int j = v_j[w];
        if ( i==o || i==b || j==o || j==b ) continue;

        double pt_w1 = jets_[i].Pt();
        double pt_w2 = jets_[j].Pt();

	///
	//  W Mass.
	///
	LorentzVector hadW = jets_[i] + jets_[j];
	double massW = hadW.mass();

	double c1 = v_k1[w];
	double c2 = v_k2[w];

	///
	// Top Mass.
	///
        LorentzVector hadT = (jets_[i] * c1) + (jets_[j] * c2) + jets_[b];
        double massT = hadT.mass();

        double pt_w = hadW.Pt();
        double sigma_w2 = pt_w1*sigma_jets_[i] * pt_w1*sigma_jets_[i]
	  + pt_w2*sigma_jets_[j] * pt_w2*sigma_jets_[j];
        double smw2 = (1.+2.*pt_w*pt_w/massW/massW)*sigma_w2;
        double pt_t = hadT.Pt();
        double sigma_t2 = c1*pt_w1*sigma_jets_[i] * c1*pt_w1*sigma_jets_[i]
	  + c2*pt_w2*sigma_jets_[j] * c2*pt_w2*sigma_jets_[j]
	  + pt_b*sigma_jets_[b] * pt_b*sigma_jets_[b];
        double smtop2 = (1.+2.*pt_t*pt_t/massT/massT)*sigma_t2;

        double c_chi2 = (massT-PDG_TOP_MASS)*(massT-PDG_TOP_MASS)/smtop2
	  + (massW-PDG_W_MASS)*(massW-PDG_W_MASS)/smw2;

        bool c_match = ( !isData_ &&  iw1[0]==i && iw2[0]==j && ib[0]==b && ibl[0]==o );

        Candidate c;
        c.chi2  = c_chi2;
        c.mt2b  = c_mt2b;
        c.mt2w  = c_mt2w;
        c.mt2bl = c_mt2bl;
        c.j1 = i;
        c.j2 = j;
        c.bi = b;
        c.oi = o;
        c.k1 = c1;
        c.k2 = c2;
        c.match = c_match;

        candidates_.push_back(c);
      }
    }

  if (__SORT)
    candidates_.sort(compare_in_chi2);

}

//--------------------------------------------------------------------
void PartonCombinatorics::applyBConsistency(float btagcut){
	int n_btag = 0;
	vector<int> non_bjets;
	for( int i = 0 ; i < n_jets_ ; i++ ){
		if( btag_.at(i) > btagcut )
			n_btag++;
		else
			non_bjets.push_back(i);
	}

	for(list<Candidate>::iterator c_it = candidates_.begin() ; c_it != candidates_.end() ; c_it++ ){
		int bi = c_it->bi;
		int oi = c_it->oi;

		bool b_btag  = (btag_.at(bi)  > btagcut);
		bool o_btag  = (btag_.at(oi)  > btagcut);

		if (n_btag == 0){
			if ( bi > (NUM_LEAD_JETS_0B - 1) ||
			     oi > (NUM_LEAD_JETS_0B - 1) )
				continue;
		} else if (n_btag == 1){
			if ( b_btag ){
				bool is_lead = false;
				for (int i=0; i < NUM_LEAD_JETS_1B ; i++)
					if (oi == non_bjets.at(i))
						is_lead = true;

				if ( not is_lead )
					continue;

			} else if ( o_btag ){
				bool is_lead = false;
				for (int i=0; i < NUM_LEAD_JETS_1B ; i++)
					if (bi == non_bjets.at(i))
						is_lead = true;

				if ( not is_lead )
					continue;

			} else {
				continue;
			}
		} else if (n_btag >= 2){
			if ( not (b_btag and o_btag) )
				continue;
		}

		b_candidates_.push_back(*c_it);
	}
}

//--------------------------------------------------------------------
void PartonCombinatorics::MT2CHI2Calculator(){
	mt2chi2_.one_chi2    = -0.999;
	mt2chi2_.two_mt2b    = -0.999;
	mt2chi2_.two_mt2bl   = -0.999;
	mt2chi2_.two_mt2w    = -0.999;
	mt2chi2_.three_mt2b  = -0.999;
	mt2chi2_.three_mt2bl = -0.999;
	mt2chi2_.three_mt2w  = -0.999;
	mt2chi2_.four_chi2b  = -0.999;
	mt2chi2_.four_chi2bl = -0.999;
	mt2chi2_.four_chi2w  = -0.999;

  if (b_candidates_.size() == 0)
    return;

  // Calculate Variable 1
  b_candidates_.sort(compare_in_chi2);
  mt2chi2_.one_chi2 = b_candidates_.front().chi2;

  //Calculate Variable 2b, 2bl, 2bw
  mt2chi2_.two_mt2b = min_with_value(b_candidates_, mt2chi2_.one_chi2, "chi2", "mt2b");
  mt2chi2_.two_mt2bl = min_with_value(b_candidates_, mt2chi2_.one_chi2, "chi2", "mt2bl");
  mt2chi2_.two_mt2w = min_with_value(b_candidates_, mt2chi2_.one_chi2, "chi2", "mt2w");

  //Calculate Variable 3b, 3bl, 3bw
  b_candidates_.sort(compare_in_mt2b);
  mt2chi2_.three_mt2b = b_candidates_.front().mt2b;

  b_candidates_.sort(compare_in_mt2bl);
  mt2chi2_.three_mt2bl = b_candidates_.front().mt2bl;

  b_candidates_.sort(compare_in_mt2w);
  mt2chi2_.three_mt2w = b_candidates_.front().mt2w;

  //Calculate Variable 4b, 4bl, 4w
  mt2chi2_.four_chi2b  = min_with_value(b_candidates_, mt2chi2_.three_mt2b, "mt2b", "chi2");
  mt2chi2_.four_chi2bl = min_with_value(b_candidates_, mt2chi2_.three_mt2bl, "mt2bl", "chi2");
  mt2chi2_.four_chi2w  = min_with_value(b_candidates_, mt2chi2_.three_mt2w,  "mt2w", "chi2");
}
