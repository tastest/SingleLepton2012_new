#include "HadronicTop.h"

#include "TROOT.h"
#include "TFile.h"
#include "TList.h"
#include "TObject.h"
#include "TObjArray.h"
#include "TH2.h"
#include "TKey.h"
#include "TRegexp.h"
#include <TIterator.h>
#include "TClass.h"
#include "TLegend.h"

#include <iostream>
#include <cmath>

list<Candidate> recoHadronicTop(JetSmearer* jetSmearer, bool isData,
                                 LorentzVector* lep, double met, double metphi,
                                 VofP4 jets, std::vector<float> btag){

  assert( jets.size() == btag.size() );

  float metx = met * cos( metphi );
  float mety = met * sin( metphi );

  int n_jets = jets.size();
  double sigma_jets[n_jets];
  for (int i=0; i<n_jets; ++i)
    sigma_jets[i] = getJetResolution(jets[i], jetSmearer);

  if ( isData )
    for (int i=0; i<n_jets; ++i)
      sigma_jets[i] *= getDataMCRatio(jets[i].eta());

  
  vector<int> mc;
  for (unsigned int i=0; i<jets.size(); i++)
    mc.push_back( isGenQGMatched( jets.at(i), 0.4 ) );

  int ibl[5];
  int iw1[5];
  int iw2[5];
  int ib[5];

  if ( !isData ){
     
     // Matching MC algoritm search over all conbinations  until the 
     // right combination is found. More than one candidate is suported 
     //  but later only the first is used.
     // 
    int match = 0;
    for (int jbl=0; jbl<n_jets; ++jbl )
      for (int jb=0; jb<n_jets; ++jb )
        for (int jw1=0; jw1<n_jets; ++jw1 )
          for (int jw2=jw1+1; jw2<n_jets; ++jw2 )
            if ( (mc.at(jw2)==2 && mc.at(jw1)==2 && mc.at(jb)==1 && mc.at(jbl)==-1) ||
                 (mc.at(jw2)==-2 && mc.at(jw1)==-2 && mc.at(jb)==-1 && mc.at(jbl)==1) ) {
                    ibl[match] = jbl;
                    iw1[match] = jw1;
                    iw2[match] = jw2;
                    ib[match] = jb;
                    match++;
            }
  }


  
////////    * Combinatorics. j_1 Pt must be > PTMIN_W1 and so on.
  
  vector<int> v_i, v_j;
  vector<double> v_k1, v_k2;
  for ( int i=0; i<n_jets; ++i )
    for ( int j=i+1; j<n_jets; ++j ){
      double pt_w1 = jets[i].Pt();
      double pt_w2 = jets[j].Pt();
      if (pt_w1 < PTMIN_J1  || pt_w2 < PTMIN_J2)
        continue;

      //
      //  W
      //
      LorentzVector hadW = jets[i] + jets[j];

      //
      //  W Mass Constraint.
      //
      TFitter *minimizer = new TFitter();
      double p1 = -1;

      minimizer->ExecuteCommand("SET PRINTOUT", &p1, 1);
      minimizer->SetFCN(minuitFunction);
      minimizer->SetParameter(0 , "c1"     , 1.1             , 1 , 0 , 0);
      minimizer->SetParameter(1 , "pt1"    , 1.0             , 1 , 0 , 0);
      minimizer->SetParameter(2 , "sigma1" , sigma_jets[i]   , 1 , 0 , 0);
      minimizer->SetParameter(3 , "pt2"    , 1.0             , 1 , 0 , 0);
      minimizer->SetParameter(4 , "sigma2" , sigma_jets[j]   , 1 , 0 , 0);
      minimizer->SetParameter(5 , "m12"    , jets[i].mass2() , 1 , 0 , 0);
      minimizer->SetParameter(6 , "m22"    , jets[j].mass2() , 1 , 0 , 0);
      minimizer->SetParameter(7 , "m02"    , hadW.mass2()    , 1 , 0 , 0);

      for (unsigned int k = 1; k < 8; k++)
        minimizer->FixParameter(k);

      minimizer->ExecuteCommand("SIMPLEX", 0, 0);
      minimizer->ExecuteCommand("MIGRAD", 0, 0);

      double c1 = minimizer->GetParameter(0);
      double c2 = fc2(c1, jets[i].mass2(), jets[j].mass2(), hadW.mass2());
                
      delete minimizer;

     
  //     * W Mass check :)
  //     *  Never trust a computer you can't throw out a window. 
 //      *  - Steve Wozniak 

//      cout << "c1 = " <<  c1 << "  c1 = " << c2 << "   M_jj = " 
//           << ((jets[i] * c1) + (jets[j] * c2)).mass() << endl;
      
      v_i.push_back(i);
      v_j.push_back(j);
      v_k1.push_back(c1);
      v_k2.push_back(c2);
    }


  list<Candidate> chi2candidates;
        
  mt2_bisect::mt2 mt2_event;
  mt2bl_bisect::mt2bl mt2bl_event;
  mt2w_bisect::mt2w mt2w_event;
  
  for ( int b=0; b<n_jets; ++b )
    for (int o=0; o<n_jets; ++o){
      if ( b == o )
        continue;

      if ( btag[b] < BTAG_MIN && btag[o] < BTAG_MIN )
        continue;

      double pt_b = jets[b].Pt();
      if ( btag[b] >= BTAG_MIN && pt_b < PTMIN_BTAG )
        continue;

      if ( btag[b] < BTAG_MIN && pt_b < PTMIN_B )
        continue;

      double pt_o = jets[o].Pt();
      if ( btag[o] >= BTAG_MIN && pt_o < PTMIN_OTAG )
        continue;

      if ( btag[o] < BTAG_MIN && pt_o < PTMIN_O)
        continue;

      ///
      //  MT2 Variables
      ///
         
      double pl[4];     // Visible lepton
      double pb1[4];    // bottom on the same side as the visible lepton
      double pb2[4];    // other bottom, paired with the invisible W
      double pmiss[3];  // <unused>, pmx, pmy   missing pT
      pl[0]= lep->E(); pl[1]= lep->Px(); pl[2]= lep->Py(); pl[3]= lep->Pz();
      pb1[0] = jets[o].E();  pb1[1] = jets[o].Px(); 
      pb1[2] = jets[o].Py(); pb1[3] = jets[o].Pz();
      pb2[0] = jets[b].E();  pb2[1] = jets[b].Px(); 
      pb2[2] = jets[b].Py(); pb2[3] = jets[b].Pz();
      pmiss[0] = 0.; pmiss[1] = metx; pmiss[2] = mety;

      double pmiss_lep[3];
      pmiss_lep[0] = 0.;
      pmiss_lep[1] = pmiss[1]+pl[1]; pmiss_lep[2] = pmiss[2]+pl[2];

      mt2_event.set_momenta( pb1, pb2, pmiss_lep );
      mt2_event.set_mn( 0.0 );   // Invisible particle mass
      double c_mt2b = mt2_event.get_mt2();
         
      mt2bl_event.set_momenta(pl, pb1, pb2, pmiss); 
      double c_mt2bl = mt2bl_event.get_mt2bl();

      mt2w_event.set_momenta(pl, pb1, pb2, pmiss); 
      double c_mt2w = mt2w_event.get_mt2w();

//      cout << b << ":"<< btag[b] << " - " << o << ":" << btag[o] << " = " << c_mt2w << endl;

      for (unsigned int w = 0; w < v_i.size() ; ++w ){
        int i = v_i[w];
        int j = v_j[w];
        if ( i==o || i==b || j==o || j==b )
            continue;

        double pt_w1 = jets[i].Pt();
        double pt_w2 = jets[j].Pt();

	///
	//  W Mass.
	///
	LorentzVector hadW = jets[i] + jets[j];
	double massW = hadW.mass();

	double c1 = v_k1[w];
	double c2 = v_k2[w];

	///
	// Top Mass.
	///
        LorentzVector hadT = (jets[i] * c1) + (jets[j] * c2) + jets[b];
        double massT = hadT.mass();

        double pt_w = hadW.Pt();
        double sigma_w2 = pt_w1*sigma_jets[i] * pt_w1*sigma_jets[i]
                        + pt_w2*sigma_jets[j] * pt_w2*sigma_jets[j];
        double smw2 = (1.+2.*pt_w*pt_w/massW/massW)*sigma_w2;
        double pt_t = hadT.Pt();
        double sigma_t2 = c1*pt_w1*sigma_jets[i] * c1*pt_w1*sigma_jets[i]
                        + c2*pt_w2*sigma_jets[j] * c2*pt_w2*sigma_jets[j]
                        + pt_b*sigma_jets[b] * pt_b*sigma_jets[b];
        double smtop2 = (1.+2.*pt_t*pt_t/massT/massT)*sigma_t2;

        double c_chi2 = (massT-PDG_TOP_MASS)*(massT-PDG_TOP_MASS)/smtop2
                      + (massW-PDG_W_MASS)*(massW-PDG_W_MASS)/smw2;

        bool c_match = ( !isData &&  iw1[0]==i && iw2[0]==j && ib[0]==b && ibl[0]==o );

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

        chi2candidates.push_back(c);
      }
    }

   if (__SORT) 
     chi2candidates.sort(compare_candidates);

   return chi2candidates;
}
