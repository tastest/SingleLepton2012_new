#include "stopUtils.h"

#include "TLorentzVector.h"
#include "Math/LorentzVector.h"
#include "Math/VectorUtil.h"
#include "Math/Point3D.h"
#include "TMath.h"
#include "TBranch.h"
#include "TTree.h"
#include "TH1F.h"
#include "TFile.h"
#include "TBits.h"

#include <vector> 
#include <list> 
#include <iostream>
#include <fstream>

#include "TFitter.h"
#include "MT2Utility.h"
#include "mt2bl_bisect.h"
#include "mt2w_bisect.h"

using namespace Stop;
using namespace std;

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > LorentzVector;

float getThetaStar(LorentzVector& vec1, LorentzVector& vec2 )
{

  float costheta=9999;

  TLorentzVector pair;
  pair.SetPtEtaPhiE((vec1+vec2).pt(), (vec1+vec2).eta(), (vec1+vec2).phi(), (vec1+vec2).energy());

  TLorentzVector child;
  child.SetPtEtaPhiE(vec1.pt(), vec1.eta(), vec1.phi(), vec1.energy());
  child.Boost(-pair.BoostVector());

  TLorentzVector child2;
  child2.SetPtEtaPhiE(vec2.pt(), vec2.eta(), vec2.phi(), vec2.energy());
  child2.Boost(-pair.BoostVector());

  // Do not reweight if by any reason top/fermion directions are undefined 
  // This should be pathological if things are fine                                                                                                                                                   
  if (pair.P()<=0 || child.P()<=0) {
    printf("Warning: particles at rest, no weight applied: pPAIR: %.3e, pf: %.3e\n", pair.P(), child.P());

    return costheta;

  }

  costheta =cos(deltaPhi(child.Phi(),pair.Phi()));

  return costheta;

}




float passLRM(double pt) {

  float b=0.1;
  float a=(-0.1)/(150.);

  float lrmCut = a * pt +b;

  return lrmCut;

}


unsigned int getNJets(const float etacut){

  unsigned int njets=0;

  for ( unsigned int i=0; i<stopt.pfjets().size() ; i++) {
    
    if( stopt.pfjets().at(i).pt()<30 )  continue;
    if( fabs(stopt.pfjets().at(i).eta())>etacut )  continue;
    if ( (fabs(stopt.pfjets().at(i).eta()) < 2.5) 
	&& (stopt.pfjets_beta2_0p5().at(i)<0.2) ) continue;

    njets++;

  }

  return njets;

}

Float_t getCSVNonb(vector<LorentzVector> jets, vector<float> csv, float discr, double ptTH, double etaTH, bool doTAU) {

  Float_t max=-1.;

  for( unsigned int i = 0 ; i < jets.size() ; ++i ){

    if( jets.at(i).pt() < ptTH)  continue;
    if( fabs(jets.at(i).eta()) > etaTH )  continue;

    if(csv.at(i)>=discr) continue;

    if( doTAU && stopt.pfTau_leadPtcandID()!=(-1) ) {

      if(deltaR(jets.at(i).eta() , jets.at(i).phi() , stopt.pfTau().eta(), stopt.pfTau().phi())<0.4) continue;

    }

    if(csv.at(i)>max) max=csv.at(i);

  }

  return max;

}



vector<int> getBJetIndex(double discr, int iskip1, int iskip2, vector<LorentzVector>jets, vector<float>csv, vector<float> lrm, double ptTH, double etaTH, bool doLRM, bool doTAU)
{

  vector<int> btagJets;
  
  for ( int i=0; i< (int) jets.size() ; i++) {
    
    if( jets.at(i).pt() < ptTH)  continue;
    if( fabs(jets.at(i).eta()) > etaTH )  continue;
    if( csv.at(i)      < discr   ) continue;
    if( doLRM && lrm.at(i) < passLRM(jets.at(i).pt())) continue;
    
    if( doTAU && stopt.pfTau_leadPtcandID()!=(-1) ) {
      
      if(deltaR(jets.at(i).eta() , jets.at(i).phi() , stopt.pfTau().eta(), stopt.pfTau().phi())<0.4) continue;
      
    }
    
    // skip these indices			\
      
    if( int(i) == iskip1 ) continue;
    if( int(i) == iskip2 ) continue;
      
    //    nbtags++;                                                                                                                                                                 
    btagJets.push_back(i);

  }

  return btagJets;
  
}

int leadingJetIndex(vector<LorentzVector> jets, int iskip1 = -1, int iskip2 = -1){

  int imaxpt  = -1;
  float maxpt = -1.0;

  // loop over jets
  for ( int i = 0 ; i < (int) jets.size() ; ++i ){      
   
    // skip these indices
    if( i == iskip1 ) continue;
    if( i == iskip2 ) continue;

    // store index of max pt jet
    if( jets.at(i).pt() > maxpt ){
      maxpt  = jets.at(i).pt();
      imaxpt = i;
    }
  }

  return imaxpt;
}

//------------------------------------------------------------------------------------------------
//this is for the Mbb
//------------------------------------------------------------------------------------------------

pair<int, int> getIndexPair(vector<int> listBJetIndex,vector<LorentzVector> jets, bool docleaning,int skip1, int skip2)
{

  bool doMin=false;
  bool doDR=false;
  bool doMass=true;

  Float_t maxptbb=-1;
  Float_t minDRbb=10.;
  Float_t minptbb=999999.;
  Float_t minMassbb=999999.;

  //  vector<pair<int, int>> index_pair_Vec;                                                                                                                                                   

  pair<int, int> index_pair=make_pair(-1.,-1.);

  for ( int k=0; k < (int) listBJetIndex.size() ; k++) {

    for ( int l=0; l < (int) listBJetIndex.size() ; l++) {

      if(l>=k) continue;

      //      if(i==skip1 || j==skip2 || i==skip2 || j==skip1) continue;

      int i=listBJetIndex.at(k);
      int j=listBJetIndex.at(l);

      double deltaRHbb = deltaR(jets.at(i).eta(),jets.at(i).phi(), jets.at(j).eta(),jets.at(j).phi());
      double massHbb = (jets.at(i)+jets.at(j)).M();
      double ptHbb = (jets.at(i)+jets.at(j)).pt();
      double deltaPhibb=deltaPhi(jets.at(i).phi(),jets.at(j).phi());

      //      float condition = (massHbb/ptHbb)*deltaRHbb;                                                                                                                                     
      //      if(condition>1.5 && docleaning) continue;                                                                                                                                        
      float conditionAngle = (massHbb/ptHbb)/deltaRHbb;
      //      float conditionAngle = (massHbb/ptHbb)/deltaPhibb;                                                                                                                               
      float conditionPt = min(jets.at(i).pt(),jets.at(j).pt())/((jets.at(i)+jets.at(j)).pt());
      float conditionRapidity = fabs(jets.at(i).Rapidity()-jets.at(j).Rapidity());

      float dPhiHL=deltaPhi((jets.at(i)+jets.at(j)).phi(),stopt.lep1().phi());

      //      if((conditionAngle>0.65 || conditionPt>0.6 ) && docleaning) continue;                                                                                                            
      //      if((conditionAngle>0.65 || conditionPt>0.6 || conditionRapidity>1.5) && docleaning) continue;                                                                                    
      //      if((conditionAngle>0.65) && docleaning) continue;                                                                                                                                
      //      if((conditionAngle>0.65 || conditionPt>0.6 || deltaPhibb>=2) && docleaning) continue;                                                                                            
      //      if((deltaPhibb>(TMath::Pi()/2)|| conditionRapidity>1.5) && docleaning) continue;                                                                                                 
      //      if(dPhiHL>(2*TMath::Pi()/3) && (listBJetIndex.size()==3) && docleaning) continue;                                                                                                
      ///      if((deltaRHbb>(2.*TMath::Pi()/3.) || conditionAngle>0.65) && docleaning) continue;
      if((deltaRHbb>(2*TMath::Pi()/3) || conditionAngle>0.65 || conditionRapidity>1.2) && docleaning) continue;
      //      if((deltaRHbb>(2*TMath::Pi()/3) || conditionAngle>0.65 || conditionPt>0.6) && docleaning) continue;

      if(!doMin && !doMass && !doDR && (jets.at(i) + jets.at(j)).pt()>maxptbb) {

        maxptbb=(jets.at(i) + jets.at(j)).pt();
        index_pair = make_pair(i,j);

      }

      if(doMin && !doMass && !doDR && (jets.at(i) + jets.at(j)).pt()<minptbb) {

        minptbb=(jets.at(i) + jets.at(j)).pt();
        index_pair = make_pair(i,j);

      }

      if(doMass && fabs((jets.at(i) + jets.at(j)).M()-125.)<minMassbb) {

        minMassbb=fabs((jets.at(i) + jets.at(j)).M()-125.);
        index_pair = make_pair(i,j);

      }

      /*                                                                                                                                                                                       
      if(doDR && deltaRHbb < minDRbb) {                                                                                                                                                        
        minDRbb= deltaRHbb;                                                                                                                                                                    
        //        minMassbb=fabs((jets.at(i) + jets.at(j)).M()-125);                                                                                                                           
        index_pair = make_pair(i,j);                                                                                                                                                           
      }                                                                                                                                                                                        
      */

    }

  }

  return index_pair;

}




//------------------------------------------------------------------------------------------------
//this is for the jetResolutions
//------------------------------------------------------------------------------------------------

float getDataMCRatio(float eta){
  if (fabs(eta) >=0.0 && fabs(eta) < 0.5) return 1.052;
  if (fabs(eta) >=0.5 && fabs(eta) < 1.1) return 1.057;
  if (fabs(eta) >=1.1 && fabs(eta) < 1.7) return 1.096;
  if (fabs(eta) >=1.7 && fabs(eta) < 2.3) return 1.134;
  if (fabs(eta) >=2.3 && fabs(eta) < 5.0) return 1.288;
  return 1.0;
}

//------------------------------------------------------------------------------------------------
// Fix function to partial correction used in looper before 1.16
//------------------------------------------------------------------------------------------------

float getDataMCRatioFix(float eta){
  if ( eta > 0 ) 
	return 1.0;
  else
     	return getDataMCRatio(eta);
}


//-------------------------------------------
// efficiency for dilepton trigger
//-------------------------------------------

float getdltrigweight(int id1, int id2){ 
  if (abs(id1)==11 && abs(id2)==11) return 0.95;
  if (abs(id1)==13 && abs(id2)==13) return 0.88;
  if (abs(id1)!=abs(id2)) return 0.92;
  return -999.;

}

float getsltrigweight(int id1, float pt, float eta) 
{
  //--------------------------------------------------------------------------------------
  // This function returns the trigger efficiencies for HLT_IsoMu24 and HLT_Ele27_WP80
  // Efficiencies are calculated with the full 2012 data sample using Z tag-and-probe
  // 3 arguments:
  // id1 : lepton PDG ID, abs(id1)=11 (electron),  abs(id1)=13 (muon)  
  // pt  : lepton pt
  // eta : lepton eta
  //--------------------------------------------------------------------------------------
  //electron efficiencies
  if ( abs(id1)==11 ) {
    
    if ( fabs(eta)<1.5) {
      if ( pt>=20  && pt<22  ) return 0.00;
      if ( pt>=22  && pt<24  ) return 0.00;
      if ( pt>=24  && pt<26  ) return 0.00;
      if ( pt>=26  && pt<28  ) return 0.07;
      if ( pt>=28  && pt<30  ) return 0.57;
      if ( pt>=30  && pt<32  ) return 0.85;
      if ( pt>=32  && pt<34  ) return 0.88;
      if ( pt>=34  && pt<36  ) return 0.89;
      if ( pt>=36  && pt<38  ) return 0.91;
      if ( pt>=38  && pt<40  ) return 0.92;
      if ( pt>=40  && pt<50  ) return 0.94;
      if ( pt>=50  && pt<60  ) return 0.95;
      if ( pt>=60  && pt<80  ) return 0.96;
      if ( pt>=80  && pt<100 ) return 0.96;
      if ( pt>=100 && pt<150 ) return 0.97;
      if ( pt>=150 && pt<200 ) return 0.97;
      if ( pt>=200           ) return 0.97;
    } 
    else if ( fabs(eta)>=1.5 && fabs(eta)<2.1) {
      if ( pt>=20  && pt<22  ) return 0.00; 
      if ( pt>=22  && pt<24  ) return 0.00; 
      if ( pt>=24  && pt<26  ) return 0.03; 
      if ( pt>=26  && pt<28  ) return 0.22; 
      if ( pt>=28  && pt<30  ) return 0.52; 
      if ( pt>=30  && pt<32  ) return 0.65; 
      if ( pt>=32  && pt<34  ) return 0.70; 
      if ( pt>=34  && pt<36  ) return 0.72; 
      if ( pt>=36  && pt<38  ) return 0.74; 
      if ( pt>=38  && pt<40  ) return 0.75; 
      if ( pt>=40  && pt<50  ) return 0.77; 
      if ( pt>=50  && pt<60  ) return 0.79; 
      if ( pt>=60  && pt<80  ) return 0.79; 
      if ( pt>=80  && pt<100 ) return 0.80; 
      if ( pt>=100 && pt<150 ) return 0.82; 
      if ( pt>=150 && pt<200 ) return 0.83; 
      if ( pt>=200           ) return 0.85; 
    }
//     else{
//       std::cout << "WARNING: electron eta " << eta << " is out of range. Return trigger efficiency = 1" << std::endl;
//     }
  } 
  //muon efficiencies
  else if ( abs(id1)==13 ) {
    if ( fabs(eta)<0.8 ) {
      if (pt>=20  && pt<22  ) return 0.00;
      if (pt>=22  && pt<24  ) return 0.02;
      if (pt>=24  && pt<26  ) return 0.87;
      if (pt>=26  && pt<28  ) return 0.90;
      if (pt>=28  && pt<30  ) return 0.91;
      if (pt>=30  && pt<32  ) return 0.91;
      if (pt>=32  && pt<34  ) return 0.92;
      if (pt>=34  && pt<36  ) return 0.93;
      if (pt>=36  && pt<38  ) return 0.93;
      if (pt>=38  && pt<40  ) return 0.93;
      if (pt>=40  && pt<50  ) return 0.94;
      if (pt>=50  && pt<60  ) return 0.95;
      if (pt>=60  && pt<80  ) return 0.95;
      if (pt>=80  && pt<100 ) return 0.94;
      if (pt>=100 && pt<150 ) return 0.94;
      if (pt>=150 && pt<200 ) return 0.93;
      if (pt>=200           ) return 0.92;
    } 
    else if ( fabs(eta)>=0.8 && fabs(eta)<1.5 ) {
      if (pt>=20  && pt<22  ) return 0.00;
      if (pt>=22  && pt<24  ) return 0.05;
      if (pt>=24  && pt<26  ) return 0.78;
      if (pt>=26  && pt<28  ) return 0.80;
      if (pt>=28  && pt<30  ) return 0.81;
      if (pt>=30  && pt<32  ) return 0.82;
      if (pt>=32  && pt<34  ) return 0.82;
      if (pt>=34  && pt<36  ) return 0.82;
      if (pt>=36  && pt<38  ) return 0.83;
      if (pt>=38  && pt<40  ) return 0.83;
      if (pt>=40  && pt<50  ) return 0.84;
      if (pt>=50  && pt<60  ) return 0.84;
      if (pt>=60  && pt<80  ) return 0.84;
      if (pt>=80  && pt<100 ) return 0.84;
      if (pt>=100 && pt<150 ) return 0.84;
      if (pt>=150 && pt<200 ) return 0.83;
      if (pt>=200           ) return 0.83;
    } 
    else if ( fabs(eta)>=1.5 && fabs(eta)<2.1 ) {
      if (pt>=20 && pt<22   ) return 0.00;
      if (pt>=22 && pt<24   ) return 0.10;
      if (pt>=24 && pt<26   ) return 0.76;
      if (pt>=26 && pt<28   ) return 0.78;
      if (pt>=28 && pt<30   ) return 0.79;
      if (pt>=30 && pt<32   ) return 0.80;
      if (pt>=32 && pt<34   ) return 0.81;
      if (pt>=34 && pt<36   ) return 0.81;
      if (pt>=36 && pt<38   ) return 0.82;
      if (pt>=38 && pt<40   ) return 0.82;
      if (pt>=40 && pt<50   ) return 0.83;
      if (pt>=50 && pt<60   ) return 0.83;
      if (pt>=60 && pt<80   ) return 0.84;
      if (pt>=80 && pt<100  ) return 0.84;
      if (pt>=100 && pt<150 ) return 0.84;
      if (pt>=150 && pt<200 ) return 0.82;
      if (pt>=200           ) return 0.83;
    }
//     else{
//       std::cout << "WARNING: muon eta " << eta << " is out of range. Return trigger efficiency = 1" << std::endl;
//     }
  }
  else{
    std::cout << "WARNING: unrecognized lepton id " << id1 << ". Return trigger efficiency = 1" << std::endl;
  }
  return 1.;
}

float getideffweight(int id1, float pt, float eta) 
{
  //--------------------------------------------------------------------------------------
  // This function returns the ID efficiencies
  // Efficiencies are calculated with the full 2012 data sample using Z tag-and-probe
  // 3 arguments:
  // id1 : lepton PDG ID, abs(id1)=11 (electron),  abs(id1)=13 (muon)  
  // pt  : lepton pt
  // eta : lepton eta
  //--------------------------------------------------------------------------------------
  //muon efficiencies
  if ( abs(id1)==13 ) {

    if ( fabs(eta)<0.8 ) {
      if ( pt<20. ) 			return 0.;
      else if (pt>=20. && pt<30.) 	return 0.9839;
      else if (pt>=30. && pt<40.) 	return 0.9850;
      else if (pt>=40. && pt<50.) 	return 0.9865;
      else if (pt>=50. && pt<60.) 	return 0.9829;
      else if (pt>=60. && pt<80.) 	return 0.9835;
      else if (pt>=80. && pt<100.) 	return 0.9785;
      else if (pt>=100. && pt<150.) 	return 0.9847;
      else if (pt>=150. && pt<200.) 	return 0.9958;
      else if (pt>=200. && pt<300.) 	return 0.9937;
      else if (pt>=300. ) 		return 0.9754;
    } else if ( fabs(eta)<1.5) {
      if ( pt<20. ) 			return 0.;
      else if (pt>=20. && pt<30.) 	return 0.9850;
      else if (pt>=30. && pt<40.) 	return 0.9846;
      else if (pt>=40. && pt<50.) 	return 0.9866;
      else if (pt>=50. && pt<60.) 	return 0.9834;
      else if (pt>=60. && pt<80.) 	return 0.9818;
      else if (pt>=80. && pt<100.) 	return 0.9803;
      else if (pt>=100. && pt<150.) 	return 0.9765;
      else if (pt>=150. && pt<200.) 	return 1.0064;
      else if (pt>=200. && pt<300.) 	return 0.9867;
      else if (pt>=300. ) 		return 1.0348; 
    } else if ( fabs(eta)<2.1) {
      if ( pt<20. ) 			return 0.;
      else if (pt>=20. && pt<30.) 	return 0.9876;
      else if (pt>=30. && pt<40.) 	return 0.9890;
      else if (pt>=40. && pt<50.) 	return 0.9902;
      else if (pt>=50. && pt<60.) 	return 0.9864;
      else if (pt>=60. && pt<80.) 	return 0.9909;
      else if (pt>=80. && pt<100.) 	return 0.9995;
      else if (pt>=100. && pt<150.) 	return 0.9884;
      else if (pt>=150. && pt<200.) 	return 0.9613;
      else if (pt>=200. && pt<300.) 	return 0.9652;
      else if (pt>=300. ) 		return 0.4286; 
    } else return 1.;
  } 
  //electron efficiencies
  else if ( abs(id1)==11 ) {
    if ( fabs(eta)<0.8 ) {
      if ( pt<20. ) 			return 0.;
      else if (pt>=20. && pt<30.) 	return 0.9923;
      else if (pt>=30. && pt<40.) 	return 0.9883;
      else if (pt>=40. && pt<50.) 	return 0.9900;
      else if (pt>=50. && pt<60.) 	return 0.9880;
      else if (pt>=60. && pt<80.) 	return 0.9847;
      else if (pt>=80. && pt<100.) 	return 0.9924;
      else if (pt>=100. && pt<150.) 	return 0.9892;
      else if (pt>=150. && pt<200.) 	return 1.0216;
      else if (pt>=200. && pt<300.) 	return 0.9869;
      else if (pt>=300. ) 		return 1.0789; 
    } else if ( fabs(eta)<1.4442 ) {
      if ( pt<20. ) 			return 0.;
      else if (pt>=20. && pt<30.) 	return 0.9632;
      else if (pt>=30. && pt<40.) 	return 0.9707;
      else if (pt>=40. && pt<50.) 	return 0.9755;
      else if (pt>=50. && pt<60.) 	return 0.9777;
      else if (pt>=60. && pt<80.) 	return 0.9797;
      else if (pt>=80. && pt<100.) 	return 0.9687;
      else if (pt>=100. && pt<150.) 	return 0.9813;
      else if (pt>=150. && pt<200.) 	return 0.9940;
      else if (pt>=200. && pt<300.) 	return 0.8853;
      else if (pt>=300. ) 		return 1.0286; 
    } else return 1.;
  }

  return 1.;

}


float getisoeffweight(int id1, float pt, float eta) 
{
  //--------------------------------------------------------------------------------------
  // This function returns the ISO efficiencies
  // Efficiencies are calculated with the full 2012 data sample using Z tag-and-probe
  // 3 arguments:
  // id1 : lepton PDG ID, abs(id1)=11 (electron),  abs(id1)=13 (muon)  
  // pt  : lepton pt
  // eta : lepton eta
  //--------------------------------------------------------------------------------------
  //muon efficiencies
  if ( abs(id1)==13 ) {

    if ( fabs(eta)<0.8 ) {
      if ( pt<20. ) 			return 0.;
      else if (pt>=20. && pt<30.) 	return 0.9934;
      else if (pt>=30. && pt<40.) 	return 0.9969;
      else if (pt>=40. && pt<50.) 	return 0.9979;
      else if (pt>=50. && pt<60.) 	return 0.9985;
      else if (pt>=60. && pt<80.) 	return 0.9989;
      else if (pt>=80. && pt<100.) 	return 0.9999;
      else if (pt>=100. && pt<150.) 	return 1.0014;
      else if (pt>=150. && pt<200.) 	return 0.9802;
      else if (pt>=200. && pt<300.) 	return 1.0016;
      else if (pt>=300. ) 		return 0.9923;
    } else if ( fabs(eta)<1.5) {
      if ( pt<20. ) 			return 0.;
      else if (pt>=20. && pt<30.) 	return 0.9974;
      else if (pt>=30. && pt<40.) 	return 1.0004;
      else if (pt>=40. && pt<50.) 	return 1.0001;
      else if (pt>=50. && pt<60.) 	return 1.0007;
      else if (pt>=60. && pt<80.) 	return 0.9997;
      else if (pt>=80. && pt<100.) 	return 1.0075;
      else if (pt>=100. && pt<150.) 	return 1.0056;
      else if (pt>=150. && pt<200.) 	return 1.0203;
      else if (pt>=200. && pt<300.) 	return 1.0059;
      else if (pt>=300. ) 		return 0.9822; 
    } else if ( fabs(eta)<2.1) {
      if ( pt<20. ) 			return 0.;
      else if (pt>=20. && pt<30.) 	return 1.0068;
      else if (pt>=30. && pt<40.) 	return 1.0039;
      else if (pt>=40. && pt<50.) 	return 1.0023;
      else if (pt>=50. && pt<60.) 	return 1.0042;
      else if (pt>=60. && pt<80.) 	return 1.0046;
      else if (pt>=80. && pt<100.) 	return 1.0086;
      else if (pt>=100. && pt<150.) 	return 1.0071;
      else if (pt>=150. && pt<200.) 	return 0.9582;
      else if (pt>=200. && pt<300.) 	return 1.0261;
      else if (pt>=300. ) 		return 1.0000; 
    } else return 1.;

  } 
  //electron efficiencies
  else if ( abs(id1)==11 ) {

    if ( fabs(eta)<0.8 ) {
      if ( pt<20. ) 			return 0.;
      else if (pt>=20. && pt<30.) 	return 0.9938;
      else if (pt>=30. && pt<40.) 	return 0.9968;
      else if (pt>=40. && pt<50.) 	return 0.9973;
      else if (pt>=50. && pt<60.) 	return 0.9957;
      else if (pt>=60. && pt<80.) 	return 0.9962;
      else if (pt>=80. && pt<100.) 	return 0.9992;
      else if (pt>=100. && pt<150.) 	return 0.9964;
      else if (pt>=150. && pt<200.) 	return 0.9861;
      else if (pt>=200. && pt<300.) 	return 1.0025;
      else if (pt>=300. ) 		return 1.1525; 
    } else if ( fabs(eta)<1.4442 ) {
      if ( pt<20. ) 			return 0.;
      else if (pt>=20. && pt<30.) 	return 0.9939;
      else if (pt>=30. && pt<40.) 	return 0.9963;
      else if (pt>=40. && pt<50.) 	return 0.9965;
      else if (pt>=50. && pt<60.) 	return 0.9963;
      else if (pt>=60. && pt<80.) 	return 0.9952;
      else if (pt>=80. && pt<100.) 	return 1.0013;
      else if (pt>=100. && pt<150.) 	return 0.9882;
      else if (pt>=150. && pt<200.) 	return 1.0068;
      else if (pt>=200. && pt<300.) 	return 1.0076;
      else if (pt>=300. ) 		return 1.0084; 
    } else return 1.;

  }

  return 1.;

 }

//-----------------------------------------------------------------------
// basic event selection:
// >=1 good lepton, rho cut, MET filters, remove 2 nearby lepton events
//-----------------------------------------------------------------------

bool passEvtSelection(TString name, bool dometdphi) 
{

  if (!name.Contains("T2")  && !name.Contains("TChiwh") && !name.Contains("T6tt")) {

    //rho requirement
    if ( stopt.rhovor()<0. || stopt.rhovor()>=40. ) return false;

    //met filters
    if ( stopt.csc()      != 0 ) return false;
    if ( stopt.hbhe()     != 1 ) return false;
    if ( stopt.hcallaser()!= 1 ) return false;
    if ( stopt.ecaltp()   != 1 ) return false;
    if ( stopt.trkfail()  != 1 ) return false;
    if ( stopt.eebadsc()  != 1 ) return false;
    if ( stopt.hbhenew()  != 1 ) return false;
    /// to comment for now for the met>50 search
    if (dometdphi && (getdphi(stopt.t1metphicorrphi(), stopt.calometphi())>1.5)) return false;

  }

  //at least 1 lepton
  if ( stopt.ngoodlep() < 1 ) return false;

  /*
  //if have more than 1 lepton, remove cases where have 2 close together
  // moved to the passDilepton functions
  if ( stopt.ngoodlep() > 1 && 
       dRbetweenVectors( stopt.lep1() ,  stopt.lep2() )<0.1 ) return false;
  */
  return true;

}


// tightness : 2=loose 1=medium 0=tight
bool passMVAJetId(double corjetpt, double jeteta, double mvavalue, unsigned int tightness)         
{
  if(tightness<0 || tightness>2)
    {
      cout << "ERROR : tightness should be 0, 1, or 2. " << endl;
      return false;
    }


  double fMVACut[3][4][4];

  /*
  // working points from full_53x_wp defined in 
  // http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/UserCode/CMG/CMGTools/External/python/JetIdParams_cfi.py?revision=1.12&view=markup

  //Tight Id
  fMVACut[0][0][0] = -0.83; fMVACut[0][0][1] = -0.81; fMVACut[0][0][2] = -0.74; fMVACut[0][0][3] = -0.81;
  fMVACut[0][1][0] = -0.83; fMVACut[0][1][1] = -0.81; fMVACut[0][1][2] = -0.74; fMVACut[0][1][3] = -0.81;
  fMVACut[0][2][0] = -0.38; fMVACut[0][2][1] = -0.32; fMVACut[0][2][2] = -0.14; fMVACut[0][2][3] = -0.48;
  fMVACut[0][3][0] = -0.38; fMVACut[0][3][1] = -0.32; fMVACut[0][3][2] = -0.14; fMVACut[0][3][3] = -0.48;
  //Medium id
  fMVACut[1][0][0] = -0.83; fMVACut[1][0][1] = -0.92; fMVACut[1][0][2] = -0.90; fMVACut[1][0][3] = -0.92;
  fMVACut[1][1][0] = -0.83; fMVACut[1][1][1] = -0.92; fMVACut[1][1][2] = -0.90; fMVACut[1][1][3] = -0.92;
  fMVACut[1][2][0] = -0.40; fMVACut[1][2][1] = -0.49; fMVACut[1][2][2] = -0.50; fMVACut[1][2][3] = -0.65;
  fMVACut[1][3][0] = -0.40; fMVACut[1][3][1] = -0.49; fMVACut[1][3][2] = -0.50; fMVACut[1][3][3] = -0.65;
  //Loose Id 
  fMVACut[2][0][0] = -0.95; fMVACut[2][0][1] = -0.96; fMVACut[2][0][2] = -0.94; fMVACut[2][0][3] = -0.95;
  fMVACut[2][1][0] = -0.95; fMVACut[2][1][1] = -0.96; fMVACut[2][1][2] = -0.94; fMVACut[2][1][3] = -0.95;
  fMVACut[2][2][0] = -0.80; fMVACut[2][2][1] = -0.74; fMVACut[2][2][2] = -0.68; fMVACut[2][2][3] = -0.77;
  fMVACut[2][3][0] = -0.80; fMVACut[2][3][1] = -0.74; fMVACut[2][3][2] = -0.68; fMVACut[2][3][3] = -0.77;
  */

  // working points from full_5x_wp defined in 
  // http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/UserCode/CMG/CMGTools/External/python/JetIdParams_cfi.py?revision=1.12&view=markup
  //Tight Id                                                                                                                                                                       
  fMVACut[0][0][0] = -0.47; fMVACut[0][0][1] = -0.92; fMVACut[0][0][2] = -0.92; fMVACut[0][0][3] = -0.94;
  fMVACut[0][1][0] = -0.47; fMVACut[0][1][1] = -0.92; fMVACut[0][1][2] = -0.92; fMVACut[0][1][3] = -0.94;
  fMVACut[0][2][0] = +0.32; fMVACut[0][2][1] = -0.49; fMVACut[0][2][2] = -0.61; fMVACut[0][2][3] = -0.74;
  fMVACut[0][3][0] = +0.32; fMVACut[0][3][1] = -0.49; fMVACut[0][3][2] = -0.61; fMVACut[0][3][3] = -0.74;
  //Medium id
  fMVACut[1][0][0] = -0.83; fMVACut[1][0][1] = -0.96; fMVACut[1][0][2] = -0.95; fMVACut[1][0][3] = -0.96;
  fMVACut[1][1][0] = -0.83; fMVACut[1][1][1] = -0.96; fMVACut[1][1][2] = -0.95; fMVACut[1][1][3] = -0.96;
  fMVACut[1][2][0] = -0.40; fMVACut[1][2][1] = -0.74; fMVACut[1][2][2] = -0.76; fMVACut[1][2][3] = -0.81;
  fMVACut[1][3][0] = -0.40; fMVACut[1][3][1] = -0.74; fMVACut[1][3][2] = -0.76; fMVACut[1][3][3] = -0.81;
  //Loose Id 
  fMVACut[2][0][0] = -0.95; fMVACut[2][0][1] = -0.97; fMVACut[2][0][2] = -0.97; fMVACut[2][0][3] = -0.97;
  fMVACut[2][1][0] = -0.95; fMVACut[2][1][1] = -0.97; fMVACut[2][1][2] = -0.97; fMVACut[2][1][3] = -0.97;
  fMVACut[2][2][0] = -0.80; fMVACut[2][2][1] = -0.85; fMVACut[2][2][2] = -0.84; fMVACut[2][2][3] = -0.85;
  fMVACut[2][3][0] = -0.80; fMVACut[2][3][1] = -0.85; fMVACut[2][3][2] = -0.84; fMVACut[2][3][3] = -0.85;


  // pT categorization
  int ptId = 0;
  if( corjetpt > 10 && corjetpt < 20 ) ptId = 1;
  if( corjetpt > 20 && corjetpt < 30 ) ptId = 2;
  if( corjetpt > 30                  ) ptId = 3;

  // eta categorization
  int etaId = 0;
  if( fabs(jeteta) > 2.5  && fabs(jeteta) < 2.75 ) etaId = 1;
  if( fabs(jeteta) > 2.75 && fabs(jeteta) < 3.0  ) etaId = 2;
  if( fabs(jeteta) > 3.0  && fabs(jeteta) < 5.0  ) etaId = 3;

  // return  
  if( mvavalue > fMVACut[tightness][ptId][etaId] ) return true;
  return false;
}

//-------------------------------------------
// the isolated track veto used at HCP
//-------------------------------------------

bool passTauVeto() {

  if(stopt.pfTau_leadPtcandID()!=(-1)) return false;

  return true;

}

//-------------------------------------------
// the isolated track veto used at HCP
//-------------------------------------------

bool passIsoTrkVeto() 
{

  //pass isolated track veto
  //unfortunately changed default value to 9999.
  if ( stopt.pfcandpt10() <9998. && stopt.pfcandiso10() < 0.1 ) return false;

  return true;

}

//-------------------------------------------
// the isolated track veto looser E+mu
//-------------------------------------------

bool passIsoTrkVeto_v2() 
{

  //pass isolated track veto
  //unfortunately changed default value to 9999.
  if ( stopt.pfcandpt10() <9998. && stopt.pfcandiso10() < 0.1 ) return false;
  if ( stopt.pfcandpt5()  <9998. && abs(stopt.pfcandid5())==13 && stopt.pfcandiso5() < 0.2) return false;
  if ( stopt.pfcandpt5()  <9998. && abs(stopt.pfcandid5())==11 && stopt.pfcandiso5() < 0.2) return false;

  return true;

}


//-------------------------------------------
// the isolated track veto looser E+mu and OSTrack
//-------------------------------------------

bool passIsoTrkVeto_v3() 
{

  //pass isolated track veto
  //unfortunately changed default value to 9999.
  if ( stopt.pfcandptOS10() <9998. && abs(stopt.pfcandidOS10())!=13 && abs(stopt.pfcandidOS10())!=11 && stopt.pfcandisoOS10() < 0.1 ) return false;
  if ( stopt.pfcandpt5()  <9998. && abs(stopt.pfcandid5())==13 && stopt.pfcandiso5() < 0.2) return false;
  if ( stopt.pfcandpt5()  <9998. && abs(stopt.pfcandid5())==11 && stopt.pfcandiso5() < 0.2) return false;

  return true;

}


//-------------------------------------------
// the isolated track veto with looser dz
//-------------------------------------------

bool passIsoTrkVeto_v4()
{

  //pass isolated track veto
  //unfortunately changed default value to 9999.
  // We want to check for the generic track only there is now good e/mu candidate
  if ( stopt.pfcandptOS10looseZ() <9998. && abs(stopt.pfcandid5looseZ())!=13 && abs(stopt.pfcandid5looseZ())!=11 && stopt.pfcandisoOS10looseZ() < 0.1 ) return false;

  if ( stopt.pfcandpt5looseZ()  <9998. && abs(stopt.pfcandid5looseZ())==13 && stopt.pfcandiso5looseZ() < 0.2) return false;
  if ( stopt.pfcandpt5looseZ()  <9998. && abs(stopt.pfcandid5looseZ())==11 && stopt.pfcandiso5looseZ() < 0.2) return false;

  return true;

}


//-------------------------------------------
// >=1 selected lepton and trigger
//-------------------------------------------

bool passSingleLeptonSelection(bool isData) 
{
  //single lepton selection for 8 TeV 53 analysis

  //at least one lepton
  if ( stopt.ngoodlep() < 1 ) return false;

  //lepton flavor - trigger, pt and eta requirements
  if ( stopt.leptype() == 0 && stopt.lep1().Pt() < 30 )          return false;
  if ( stopt.leptype() == 1 && stopt.lep1().Pt() < 25 )          return false;

  if ( fabs( stopt.pflep1().Pt() - stopt.lep1().Pt() ) > 10. )  return false;
  if ( ( stopt.isopf1() * stopt.lep1().Pt() ) > 5. )  return false; 

  if ( stopt.leptype() == 0 ) {

    //pass trigger if data - single electron
    if ( isData && stopt.ele27wp80() != 1 ) return false;
    //    if ( isData && stopt.trgel1() != 1 )  return false;
 
    //barrel only electrons
    if ( fabs(stopt.lep1().Eta() ) > 1.4442) return false;
    if ( stopt.eoverpin() > 4. ) return false;
 

  } else if ( stopt.leptype() == 1 ) {

    //pass trigger if data - single muon
    if ( isData && stopt.isomu24() != 1 ) return false;
    //    if ( isData && stopt.trgmu1() != 1 )  return false;
 
    if ( fabs(stopt.lep1().Eta() ) > 2.1)  return false;

  }

  return true;

}

//-------------------------------------------
// the dilepton OS selection
//-------------------------------------------

bool passDileptonSelection(bool isData) 
{
  //two lepton selection for 8 TeV 53 analysis

  //exactly 2 leptons
  if ( stopt.ngoodlep() != 2 ) return false;

  //opposite sign
  if ( stopt.id1()*stopt.id2()>0 ) return false;

  //pass trigger if data - dilepton
  if ( isData && stopt.mm() != 1 && stopt.me() != 1 
       && stopt.em() != 1 && stopt.ee() != 1 ) return false;

  //passes pt and eta requirements
  if ( stopt.lep1().Pt() < 20 )          return false;
  if ( stopt.lep2().Pt() < 20 )          return false;
  if ( fabs(stopt.lep1().Eta() ) > 2.4)  return false;
  if ( fabs(stopt.lep2().Eta() ) > 2.4)  return false;

  //consistency with pf leptons
  if ( fabs( stopt.pflep1().Pt() - stopt.lep1().Pt() ) > 10. )  return false;
  if ( fabs( stopt.pflep2().Pt() - stopt.lep2().Pt() ) > 10. )  return false;

  //requirement for leading lepton
  if ( ( stopt.isopf1() * stopt.lep1().Pt() ) > 5. )  return false; 
  if ( fabs(stopt.id1())==11 && stopt.eoverpin() > 4. ) return false;

  //requirement for trailing lepton
  if ( ( stopt.isopf2() * stopt.lep2().Pt() ) > 5. )  return false; 
  if ( fabs(stopt.id2())==11 && stopt.eoverpin2() > 4. ) return false;

  //barrel only electrons
  if (fabs(stopt.id1())==11 && fabs(stopt.lep1().Eta() ) > 1.4442) return false;
  if (fabs(stopt.id2())==11 && fabs(stopt.lep2().Eta() ) > 1.4442) return false;

  //if have more than 1 lepton, remove cases where have 2 close together
  if ( stopt.ngoodlep() > 1 && 
       dRbetweenVectors( stopt.lep1() ,  stopt.lep2() )<0.1 ) return false;

  return true;

}

//-------------------------------------------
// the dilepton SS selection
//-------------------------------------------

bool passDileptonSSSelection(bool isData) 
{
  //two lepton selection for 8 TeV 53 analysis

  //exactly 2 leptons
  if ( stopt.ngoodlep() != 2 ) return false;

  //same sign
  if ( stopt.id1()*stopt.id2()<0 ) return false;

  //pass trigger if data - dilepton
  if ( isData && stopt.mm() != 1 && stopt.me() != 1 
       && stopt.em() != 1 && stopt.ee() != 1 ) return false;

  //passes pt and eta requirements
  if ( stopt.lep1().Pt() < 20 )          return false;
  if ( stopt.lep2().Pt() < 20 )          return false;
  if ( fabs(stopt.lep1().Eta() ) > 2.4)  return false;
  if ( fabs(stopt.lep2().Eta() ) > 2.4)  return false;

  //consistency with pf leptons
  if ( fabs( stopt.pflep1().Pt() - stopt.lep1().Pt() ) > 10. )  return false;
  if ( fabs( stopt.pflep2().Pt() - stopt.lep2().Pt() ) > 10. )  return false;

  //requirement for leading lepton
  if ( ( stopt.isopf1() * stopt.lep1().Pt() ) > 5. )  return false; 
  if ( fabs(stopt.id1())==11 && stopt.eoverpin() > 4. ) return false;

  //requirement for trailing lepton
  if ( ( stopt.isopf2() * stopt.lep2().Pt() ) > 5. )  return false; 
  if ( fabs(stopt.id2())==11 && stopt.eoverpin2() > 4. ) return false;

  //barrel only electrons
  if (fabs(stopt.id1())==11 && fabs(stopt.lep1().Eta() ) > 1.4442) return false;
  if (fabs(stopt.id2())==11 && fabs(stopt.lep2().Eta() ) > 1.4442) return false;

  //if have more than 1 lepton, remove cases where have 2 close together
  if ( stopt.ngoodlep() > 1 && 
       dRbetweenVectors( stopt.lep1() ,  stopt.lep2() )<0.1 ) return false;

  return true;

}

//-------------------------------------------
// good lepton + veto isolated track
//-------------------------------------------

bool passOneLeptonSelection( bool isData) 
{
  //single lepton selection for 8 TeV 53 analysis
  if ( !passSingleLeptonSelection(isData) ) return false;

  //pass isolated track veto
  //unfortunately changed default value to 9999.
  //  if ( stopt.pfcandpt10() <9998. && stopt.pfcandiso10() < 0.1 ) return false;
  if ( !passIsoTrkVeto_v4() ) return false;

  return true;

}

//-------------------------------------------
// dilepton selection
//-------------------------------------------

bool passTwoLeptonSelection(bool isData) 
{
  //single lepton selection for 8 TeV 53 analysis
  if ( !passDileptonSelection(isData) ) return false;

  //apply isolated track veto in addition to 2 leptons
  //default value for this one is -999
  if ( stopt.trkpt10loose() >0. && stopt.trkreliso10loose() < 0.1 ) return false;

  return true;

}


bool passLepPlusSSPionTrkSelection(bool isData)
{
  //single lepton plus iso trk selection for 8 TeV 53 analysis                                                                                                                                                         

  //at least one lepton                                                                                                                                                                                                
  if ( !passSingleLeptonSelection(isData) ) return false;

  //pass isolated track requirement                                                                                                                                                                                    
  //unfortunately changed default value to 9999.                                                                                                                                                                       
  //  if ( pfcandpt10() > 9990. || pfcandiso10() > 0.1 ) return false;                                                                                                                                                 
  if ( stopt.pfcandpt5looseZ()  > 9990. || stopt.pfcandpt5looseZ() < 10.) return false;
  if ( abs(stopt.pfcandid5looseZ())==13 || abs(stopt.pfcandid5looseZ())==11) return false;
  if ( stopt.pfcandiso5looseZ() > 0.1) return false;

  float charge=(stopt.pfcandid5looseZ()*stopt.id1());
  // charge < 0 is a SS , charge > 0 is a OS for e/mu; need to flip for pions                                                                                                                                           if((abs(stopt.pfcandidOS10looseZ())!=11) && (abs(stopt.pfcandidOS10looseZ())!=13)) charge*=(-1);                                                                                                                   

  if(charge > 0) return false;

  return true;

}



//-------------------------------------------
// lepton + isolated track selection CR5
//-------------------------------------------

bool passLepPlusIsoTrkSelection(bool isData) 
{
  //single lepton plus iso trk selection for 8 TeV 53 analysis

  //at least one lepton
  if ( !passSingleLeptonSelection(isData) ) return false;

  //pass isolated track requirement
  //unfortunately changed default value to 9999.
  //  if ( pfcandpt10() > 9990. || pfcandiso10() > 0.1 ) return false;
  if ( passIsoTrkVeto_v4() ) return false;
  if ( stopt.pfcandpt5looseZ()  > 9990.) return false;
  if ( stopt.pfcandptOS10looseZ()  > 9990.) return false;

  return true;

}

//-------------------------------------------
// lepton + isolated track selection CR5
//-------------------------------------------

bool passLepPlusIsoTrkSelectionWHMet(bool isData) 
{
  //single lepton plus iso trk selection for 8 TeV WH+MET analysis

  //at least one lepton
  if ( !passSingleLeptonSelection(isData) ) return false;

  // exactly one good lepton
  if ( stopt.ngoodlep() != 1 ) return false;

  //pass isolated track requirement
  if ( !passIsoTrkVeto_v4() ) return true;

  //unfortunately changed default value to 9999.
  //  if ( pfcandpt10() > 9990. || pfcandiso10() > 0.1 ) return false;

  // // if we have an isolated track cand, check eta for electron and muon cands
  // if ( !passIsoTrkVeto_v4() ) {
  //   if ( stopt.pfcandpt5looseZ()  <9998. && abs(stopt.pfcandid5looseZ())==13 && stopt.pfcandiso5looseZ() < 0.2) {
  //     if (fabs(stopt.pfcand5looseZ().eta()) < 2.1) return true;
  //     else return false;
  //   }
  //   else if ( stopt.pfcandpt5looseZ()  <9998. && abs(stopt.pfcandid5looseZ())==11 && stopt.pfcandiso5looseZ() < 0.2) {
  //     if (fabs(stopt.pfcand5looseZ().eta()) < 1.4442) return true;
  //     else return false;
  //   }

  //   // isolated track with pt > 10: don't restrict eta?
  //   //    if ( stopt.pfcandptOS10looseZ() <9998. && abs(stopt.pfcandid5looseZ())!=13 && abs(stopt.pfcandid5looseZ())!=11 && stopt.pfcandisoOS10looseZ() < 0.1 ) return false;
  //   else return true;

  // }

  // allow isolated tau cands, if they weren't already covered by previous cases
  if ( !passTauVeto() ) return true;

  return false;
}

//-------------------------------------------
// lepton + taus selection CR5
//-------------------------------------------

bool passLepPlusTauSelection(bool isData)
{

  //at least one lepton                                                                                                                                      
  if ( !passSingleLeptonSelection(isData) ) return false;

  //pass tauVeto                                                                                                                                             
  if ( passTauVeto() ) return false;

  return true;

}


bool passLepPlusTauSelection_v2(bool isData) 
{
  //single lepton plus iso trk selection for 8 TeV 53 analysis
  
  //at least one lepton
  if ( !passSingleLeptonSelection(isData) ) return false;
  
  // veto the second lepton 
  if(stopt.ngoodlep()>1) return false;
  
  // veto looser e/mu veto
  if ( stopt.pfcandpt5looseZ()  <9998. && abs(stopt.pfcandid5looseZ())==13 && stopt.pfcandiso5looseZ() < 0.2) return false;
  if ( stopt.pfcandpt5looseZ()  <9998. && abs(stopt.pfcandid5looseZ())==11 && stopt.pfcandiso5looseZ() < 0.2) return false;

  /*
  // this should be redundant at this stage
  if( stopt.id2()!=(-1) ) {
    if(stopt.lep2().pt() > 20) return false;
    if(deltaR(stopt.lep2().eta() , stopt.lep2().phi() , stopt.pfTau().eta(), stopt.pfTau().phi())<0.4) return false;
  }
  */
  
  //pass tauVeto
  if ( passTauVeto() ) return false;

  if(fabs(stopt.pfTau().eta())> 2.3) return false;

  return true;

}



bool pass_T2tt_LM(bool isData, TString name){

  if ( !passSingleLeptonSelection(isData) ) return false;

  if ( !passIsoTrkVeto_v4() ) return false;
  if ( !passTauVeto() ) return false;

  // this is just a preselection
  if(stopt.t1metphicorr()<100) return false;

  if(stopt.t1metphicorrmt()<120) return false;

  vector<LorentzVector> myJets;
  vector<float> myJetsTag;
  vector<int> myJetsMC;
  vector<float> myJetsSigma;
  int btag=0;

  for (int ijet =0; ijet<(int)stopt.pfjets().size(); ijet++){

    if ( stopt.pfjets().at(ijet).Pt() < 30 ) continue;
    if ( fabs(stopt.pfjets().at(ijet).eta()) > 2.4 ) continue;
    // if( stopt.pfjets_beta2_0p5().at(ijet)<0.2 )  continue;
    // bool passMediumPUid = passMVAJetId(stopt.pfjets().at(ijet).pt(), stopt.pfjets().at(ijet).eta(),stopt.pfjets_mvaPUid().at(ijet),1);
    bool passTightPUid = passMVAJetId(stopt.pfjets().at(ijet).pt(), stopt.pfjets().at(ijet).eta(),stopt.pfjets_mva5xPUid().at(ijet),0);
    if(!passTightPUid) continue;

    myJets.push_back(stopt.pfjets().at(ijet));
    myJetsTag.push_back(stopt.pfjets_csv().at(ijet));

    if(stopt.pfjets_csv().at(ijet) > 0.679) btag++;

    myJetsSigma.push_back(stopt.pfjets_sigma().at(ijet));

  }

  if(myJets.size()<4) return false;
  if(btag==0) return false;

  // mindPhi
  float dphimjmin=getMinDphi(stopt.t1metphicorrphi(), myJets.at(0),myJets.at(1));
  if(dphimjmin<0.8) return false;

  // chi2
  double chi2 = calculateChi2SNT(myJets, myJetsSigma, myJetsTag);
  if(chi2>5) return false;

  // mt2w
  //  double x_mt2w = calculateMT2w(myJets, myJetsTag, stopt.lep1(), stopt.t1metphicorr(), stopt.t1metphicorrphi());
  return true;

}


bool pass_T2tt_HM(bool isData,TString name){

  if ( !passSingleLeptonSelection(isData) ) return false;

  if ( !passIsoTrkVeto_v4() ) return false;
  if ( !passTauVeto() ) return false;
  
  // this is just a preselection                                                                                                                                                     
  if(stopt.t1metphicorr()<100) return false;

  if(stopt.t1metphicorrmt()<120) return false;

  vector<LorentzVector> myJets;
  vector<float> myJetsTag;
  vector<int> myJetsMC;
  vector<float> myJetsSigma;
  int btag;

  for (int ijet =0; ijet<(int)stopt.pfjets().size(); ijet++){

      if ( stopt.pfjets().at(ijet).Pt() < 30 ) continue;
      if ( fabs(stopt.pfjets().at(ijet).eta()) > 2.4 ) continue;
      // if( stopt.pfjets_beta2_0p5().at(ijet)<0.2 )  continue;
      // bool passMediumPUid = passMVAJetId(stopt.pfjets().at(ijet).pt(), stopt.pfjets().at(ijet).eta(),stopt.pfjets_mvaPUid().at(ijet),1);
      bool passTightPUid = passMVAJetId(stopt.pfjets().at(ijet).pt(), stopt.pfjets().at(ijet).eta(),stopt.pfjets_mva5xPUid().at(ijet),0);
      if(!passTightPUid) continue;

      myJets.push_back(stopt.pfjets().at(ijet));
      myJetsTag.push_back(stopt.pfjets_csv().at(ijet));

      if(stopt.pfjets_csv().at(ijet) > 0.679) btag++;
      myJetsSigma.push_back(stopt.pfjets_sigma().at(ijet));
  }

  if(myJets.size()<4) return false;
  if(btag==0) return false;

  // mindPhi                                                                                                                                                                         
  float dphimjmin=getMinDphi(stopt.t1metphicorrphi(), myJets.at(0),myJets.at(1));
  if(dphimjmin<0.8) return false;

  // chi2                                                                                                                                                                            
  double chi2 = calculateChi2SNT(myJets, myJetsSigma, myJetsTag);

  if(chi2>5) return false;

  // mt2w                                                                                                                                                                            
  double x_mt2w = calculateMT2w(myJets, myJetsTag, stopt.lep1(), stopt.t1metphicorr(), stopt.t1metphicorrphi());
  // cut a the top mass + some resoultion effect
  if(x_mt2w<200) return false;

  return true;


}



bool pass_T2bw_HM(bool isData,TString name){

  if ( !passSingleLeptonSelection(isData) ) return false;

  if ( !passIsoTrkVeto_v4() ) return false;
  if ( !passTauVeto() ) return false;
  
  // this is just a preselection                                                                                                                                                     
  if(stopt.t1metphicorr()<100) return false;

  if(stopt.t1metphicorrmt()<120) return false;

  vector<LorentzVector> myJets;
  vector<float> myJetsTag;
  vector<int> myJetsMC;
  vector<float> myJetsSigma;
  int btag;

  float bpt=-1;

  for (int ijet =0; ijet<(int)stopt.pfjets().size(); ijet++){

      if ( stopt.pfjets().at(ijet).Pt() < 30 ) continue;
      if ( fabs(stopt.pfjets().at(ijet).eta()) > 2.4 ) continue;
      // if( stopt.pfjets_beta2_0p5().at(ijet)<0.2 )  continue;
      // bool passMediumPUid = passMVAJetId(stopt.pfjets().at(ijet).pt(), stopt.pfjets().at(ijet).eta(),stopt.pfjets_mvaPUid().at(ijet),1);
      bool passTightPUid = passMVAJetId(stopt.pfjets().at(ijet).pt(), stopt.pfjets().at(ijet).eta(),stopt.pfjets_mva5xPUid().at(ijet),0);
      if(!passTightPUid) continue;

      myJets.push_back(stopt.pfjets().at(ijet));
      myJetsTag.push_back(stopt.pfjets_csv().at(ijet));

      if(stopt.pfjets_csv().at(ijet) > 0.679) btag++;
      if(btag==1) bpt=stopt.pfjets().at(ijet).Pt();

      myJetsSigma.push_back(stopt.pfjets_sigma().at(ijet));
  }

  if(myJets.size()<4) return false;
  if(btag==0) return false;

  // mindPhi                                                                                                                                                                       
  float dphimjmin=getMinDphi(stopt.t1metphicorrphi(), myJets.at(0),myJets.at(1));
  if(dphimjmin<0.8) return false;

  // b jets pt
  if(bpt<100) return false;

  // mt2w                                                                                                                                                                          
  double x_mt2w = calculateMT2w(myJets, myJetsTag, stopt.lep1(), stopt.t1metphicorr(), stopt.t1metphicorrphi());
  // cut a the top mass + some resoultion effect
  if(x_mt2w<200) return false;

  return true;


}


//-------------------------------------------
// on-the-fly MET phi corrections
//-------------------------------------------

bool passLepPlusIsoTrkSelection_noEMu(bool isData, bool pickMuE) 
{
  //single lepton plus iso trk selection for 8 TeV 53 analysis

  //at least one lepton
  if ( !passSingleLeptonSelection(isData) ) return false;

  //pass isolated track requirement
  //unfortunately changed default value to 9999.
  //////  if ( pfcandpt10() > 9990. || pfcandiso10() > 0.1 ) return false;
  bool foundIsoTrack_noEMU=(stopt.pfcandptOS10looseZ() <9998. && (abs(stopt.pfcandidOS10looseZ())!=13 && abs(stopt.pfcandidOS10looseZ())!=11 ) && stopt.pfcandisoOS10looseZ() < 0.1); 
  bool foundIsoTrack_EMU= (stopt.pfcandpt5looseZ() <9998. && (abs(stopt.pfcandid5looseZ())==13 || abs(stopt.pfcandid5looseZ())==11 ) && stopt.pfcandiso5looseZ() < 0.2);

  ///  if(foundIsoTrack_EMU) cout << "EMU " << foundIsoTrack_EMU << endl;
  ///  if(foundIsoTrack_noEMU) cout << "noEMU " << foundIsoTrack_noEMU << endl;

  bool muE=(pickMuE && foundIsoTrack_EMU );
  bool noMuE=(!pickMuE && foundIsoTrack_noEMU);

  if(muE) return true;
  if(noMuE) return true;

  return false;

}


//-------------------------------------------
// on-the-fly MET phi corrections
//-------------------------------------------

pair<float,float> getPhiCorrMET( float met, float metphi, int nvtx, bool ismc){

  //using met phi corrections from C. Veelken (revision 1.6)
  //functions are available here:                                                                                                            
  //http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/CMSSW/JetMETCorrections/Type1MET/python/pfMETsysShiftCorrections_cfi.py                               

  float metx = met * cos( metphi );
  float mety = met * sin( metphi );

  float shiftx = 0.;
  float shifty = 0.;

  //use correction for data vs. mc                                                                                                                  
  shiftx = ismc ? (0.1166 + 0.0200*nvtx)
    : (0.2661 + 0.3217*nvtx);
  shifty = ismc ? (0.2764 - 0.1280*nvtx)
    : (-0.2251 - 0.1747*nvtx);

  metx -= shiftx;
  mety -= shifty;

  pair<float, float> phicorrmet = make_pair( sqrt( metx*metx + mety*mety ), atan2( mety , metx ) );
  return phicorrmet;

}

//---------------------------------------
// utility functions
//---------------------------------------

float vtxweight_n( const int nvertices, TH1F *hist, bool isData ) {

  if( isData ) return 1;

  int nvtx = nvertices;
  if( nvtx > hist->GetNbinsX() )
    nvtx = hist->GetNbinsX();

  float weight = 0;
  weight = hist->GetBinContent( hist->FindBin(nvtx) );
  if( weight <= 0 ) //we don't want to kill events bc they have no weight
    weight = 1.;
  //  cout << "nvtx " << nvtx << " weight " << weight << endl;
  return weight;

}

float getdphi( float phi1 , float phi2 ){
  float dphi = fabs( phi1 - phi2 );
  if( dphi > TMath::Pi() ) dphi = TMath::TwoPi() - dphi;
  return dphi;
}

float getMT( float pt1 , float phi1 , float pt2 , float phi2 ){

  float dphi = getdphi(phi1, phi2);
  return sqrt( 2 * ( pt1 * pt2 * (1 - cos( dphi ) ) ) );

}

float dRbetweenVectors(LorentzVector& vec1,LorentzVector& vec2 ){ 

  float dphi = std::min(::fabs(vec1.Phi() - vec2.Phi()), 2 * M_PI - fabs(vec1.Phi() - vec2.Phi()));
  float deta = vec1.Eta() - vec2.Eta();

  return sqrt(dphi*dphi + deta*deta);
}

float getMinDphi(float metPhi, LorentzVector& vec1, LorentzVector& vec2 ) {

  float dphimj1_    = getdphi(metPhi, vec1.phi() );
  float dphimj2_    = getdphi(metPhi, vec2.phi() );
  float dphimjmin_  = TMath::Min( dphimj1_ , dphimj2_ );

  return dphimjmin_;

}


//---------------------------------------
// duplicates and filters
//---------------------------------------

bool DorkyEventIdentifier::operator < (const DorkyEventIdentifier &other) const
{
  if (run != other.run)
    return run < other.run;
  if (event != other.event)
    return event < other.event;
  if(lumi != other.lumi)
    return lumi < other.lumi;
  return false;
}

//--------------------------------------------------------------------

bool DorkyEventIdentifier::operator == (const DorkyEventIdentifier &other) const
{
  if (run != other.run)
    return false;
  if (event != other.event)
    return false;
  return true;
}

//--------------------------------------------------------------------

bool is_duplicate (const DorkyEventIdentifier &id, std::set<DorkyEventIdentifier> &already_seen) {
  std::pair<std::set<DorkyEventIdentifier>::const_iterator, bool> ret =
    already_seen.insert(id);
  return !ret.second;
}

//--------------------------------------------------------------------

int load_badlaserevents  (char* filename, std::set<DorkyEventIdentifier> &events_lasercalib) {

  ifstream in;
  in.open(filename);

   unsigned long int run, event, lumi;
   int nlines = 0;

   while (1) {
      in >> run >> event >> lumi;
      if (!in.good()) break;
      nlines++;
      DorkyEventIdentifier id = {run, event, lumi };
      events_lasercalib.insert(id);
   }
   printf(" found %d bad events \n",nlines);

   in.close();

   return 0;

}

//--------------------------------------------------------------------

bool is_badLaserEvent (const DorkyEventIdentifier &id, std::set<DorkyEventIdentifier> &events_lasercalib) {
  if (events_lasercalib.find(id) != events_lasercalib.end()) return true;
  return false;
}


//--------------------------------------------------------------------
double calculateMT2w(vector<LorentzVector>& jets, vector<float>& btag, LorentzVector& lep, float met, float metphi, MT2Type mt2type){

	// I am asumming that jets is sorted by Pt
	assert ( jets.size() == btag.size() );
	// require at least 2 jets
	if ( jets.size()<2 ) return 99999.; 

	// First we count the number of b-tagged jets, and separate those non b-tagged
	std::vector<int> bjets;
	std::vector<int> non_bjets;
	for( unsigned int i = 0 ; i < jets.size() ; i++ ){
	  if( btag.at(i) > BTAG_MED ) {
	    bjets.push_back(i);
	  } else {
	    non_bjets.push_back(i);
	  }
	}	

	int n_btag = (int) bjets.size();
	//	cout << "n_btag = " << n_btag << endl;

	// We do different things depending on the number of b-tagged jets
	// arXiv:1203.4813 recipe

	int nMax=-1;
	if(jets.size()<=3) nMax=non_bjets.size();
	else nMax=3;

	if (n_btag == 0){
	  // If no b-jets select the minimum of the mt2w from all combinations with 
	  // the three leading jets
	  float min_mt2w = 9999;

	  for (int i=0; i<nMax; i++)
	    for (int j=0; j<nMax; j++){
	      if (i == j) continue;
	      float c_mt2w = mt2wWrapper(lep, 
					 jets[non_bjets[i]],
					 jets[non_bjets[j]], met, metphi, mt2type);
	      if (c_mt2w < min_mt2w)
		min_mt2w = c_mt2w;
	    }
	  return min_mt2w;
	} else if (n_btag == 1 ){
	  // if only one b-jet choose the three non-b leading jets and choose the smaller
	  float min_mt2w = 9999;

	  for (int i=0; i<nMax; i++){
	    float c_mt2w = mt2wWrapper(lep, jets[bjets[0]], jets[non_bjets[i]], met, metphi, mt2type);
	    if (c_mt2w < min_mt2w)
	      min_mt2w = c_mt2w;
	  }
	  for (int i=0; i<nMax; i++){
	    float c_mt2w = mt2wWrapper(lep, jets[non_bjets[i]], jets[bjets[0]], met, metphi, mt2type);
	    if (c_mt2w < min_mt2w)
	      min_mt2w = c_mt2w;
	  }
	  return min_mt2w;
	} else if (n_btag >= 2) {
	  // if 3 or more b-jets the paper says ignore b-tag and do like 0-bjets 
	  // but we are going to make the combinations with the b-jets
	  float min_mt2w = 9999;
	  for (int i=0; i<n_btag; i++)
	    for (int j=0; j<n_btag; j++){
	      if (i == j) continue;
	      float c_mt2w = mt2wWrapper(lep, 
					 jets[bjets[i]],
					 jets[bjets[j]], met, metphi, mt2type);
	      if (c_mt2w < min_mt2w)
		min_mt2w = c_mt2w;
	    }
	  return min_mt2w;
	}

	return -1.;
}


//---------------------------------------------------------------------


// This funcion is a wrapper for mt2w_bisect etc that takes LorentzVectors instead of doubles
double mt2wWrapper(LorentzVector& lep, LorentzVector& jet_o, LorentzVector& jet_b, float met, float metphi, MT2Type mt2type){

	// same for all MT2x variables
	float metx = met * cos( metphi );
	float mety = met * sin( metphi );

	double pl[4];     // Visible lepton
	double pb1[4];    // bottom on the same side as the visible lepton
	double pb2[4];    // other bottom, paired with the invisible W
	double pmiss[3];  // <unused>, pmx, pmy   missing pT
	pl[0]= lep.E(); pl[1]= lep.Px(); pl[2]= lep.Py(); pl[3]= lep.Pz();
	pb1[1] = jet_o.Px();  pb1[2] = jet_o.Py();   pb1[3] = jet_o.Pz();
	pb2[1] = jet_b.Px();  pb2[2] = jet_b.Py();   pb2[3] = jet_b.Pz();
	pmiss[0] = 0.; pmiss[1] = metx; pmiss[2] = mety;

	// specifics for each variable
	if (mt2type == MT2b) {
	  double pmiss_lep[3];
	  pmiss_lep[0] = 0.;
	  pmiss_lep[1] = pmiss[1]+pl[1]; pmiss_lep[2] = pmiss[2]+pl[2];

	  pb1[0] = jet_o.mass();
	  pb2[0] = jet_b.mass();

	  mt2_bisect::mt2 mt2_event;
	  mt2_event.set_momenta( pb1, pb2, pmiss_lep );
	  mt2_event.set_mn( 80.385 );   // Invisible particle mass == W mass
	  return mt2_event.get_mt2();
	}

	else {
	  pb1[0] = jet_o.E();
	  pb2[0] = jet_b.E();

	  if (mt2type == MT2bl) {
	    mt2bl_bisect::mt2bl mt2bl_event;
	    mt2bl_event.set_momenta(pl, pb1, pb2, pmiss);
	    return mt2bl_event.get_mt2bl();
	  }
	  else if (mt2type == MT2w) {
	    mt2w_bisect::mt2w mt2w_event;
	    mt2w_event.set_momenta(pl, pb1, pb2, pmiss);
	    return mt2w_event.get_mt2w();
	  }
	}

	// SHOULDN'T GET HERE
	std::cout << "ERROR: " << __FILE__ << " " << __LINE__ << ": mt2wWrapper: shouldn't get here! mt2type == " 
		  << mt2type << std::endl;
	return -1.;
}


//--------------------------------------------------------------------
double fc2 (double c1, double m12, double m22, double m02, bool verbose = false)
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
double fchi2 (double c1, double pt1, double sigma1, double pt2, double sigma2,
              double m12, double m22, double m02){
  double rat1 = pt1 * (1 - c1) / sigma1;
  double rat2 = pt2 * (1 - fc2(c1, m12, m22, m02)) / sigma2;

  return ( rat1 * rat1 + rat2 * rat2);
}

//--------------------------------------------------------------------
void minuitFunction(int&, double* , double &result, double par[], int){
  result=fchi2(par[0], par[1], par[2], par[3], par[4], par[5], par[6], par[7]);
}

// This function calculates the hadronic chi2 - atlas version
double calculateChi2(vector<LorentzVector>& jets, vector<float>& sigma_jets){

	assert(jets.size() == sigma_jets.size());

	int n_jets = jets.size();

	int j1=-1;
	int j2=-1;
	int bi=-1;
	double min_dR = 9999;
	for (int i=0; i < n_jets; i++)
		for (int j=0; j < n_jets; j++){
			double dR =  ROOT::Math::VectorUtil::DeltaR(jets.at(i),jets.at(j));
			LorentzVector hadW = jets.at(i) + jets.at(j);
			if ( dR < min_dR && hadW.mass() > 60){
				min_dR = dR;
				j1 = i;
				j2 = j;
			}
		}
	if ( j1 < 0 || j2 < 0 ) return -1;
	
	LorentzVector hadW = jets.at(j1) + jets.at(j2);

	min_dR = 9999;	
	for (int j=0; j < n_jets; j++){
		double dR =  ROOT::Math::VectorUtil::DeltaR(hadW,jets.at(j));
		LorentzVector hadT = hadW + jets.at(j);
		if ( dR < min_dR && hadT.mass() > 130 && hadT.mass()< 205 ){
			min_dR = dR;
			bi = j;
		}
	}

	if ( bi <  0 ) return -1;

/*
	LorentzVector c1_direction = jets.at(0)/jets.at(0).P();
	LorentzVector c2_direction = jets.at(1)/jets.at(1).P();

	int cluster1[15];
	int cluster2[15];
	int n_cluster1, n_cluster2;

	float diff2;
	do {
		n_cluster1 = 0;
		n_cluster2 = 0;
		for(int i=0; i < n_jets; i++){
			double dRc1 = (c1_direction - jets.at(i)).P();
			double dRc2 = (c2_direction - jets.at(i)).P();
//			cout << "RR: " << dRc1 << "  " << dRc2 <<  "   " << jets.at(i) << endl;
			if ( dRc1 < dRc2 )
				cluster1[n_cluster1++] = i;
			else
				cluster2[n_cluster2++] = i;
		}

		LorentzVector c1_new;
		for ( int i=0; i<n_cluster1; i++)
			c1_new += jets.at(cluster1[i]);
		c1_new = c1_new/c1_new.P();

		LorentzVector c2_new;
		for ( int i=0; i<n_cluster2; i++)
			c2_new += jets.at(cluster2[i]);
		c2_new = c2_new/c2_new.P();

//		cout << c1_new << "    " << c2_new << endl;

		diff2 = ( c1_new - c1_direction ).P2() + ( c2_new - c2_direction ).P2();
//		cout << n_cluster1 << "   " << n_cluster2 << "   " << diff2 << endl;

		c1_direction = c1_new;
		c2_direction = c2_new;

		assert ( n_cluster1 > 0 ); 
		assert ( n_cluster2 > 0 ); 
	} while ( diff2 > 0.001 );

	cout << n_cluster1 << "   " << n_cluster2 << "   " << diff2 << endl;

	

	double min_chi2 = 9999;
	for ( int b=0; b<n_jets; ++b )
		for (int j1=0; j1<n_jets; ++j1)
			for (int j2=0; j2<n_jets; ++j2){
				if (b == j1 || b == j2 || j1 == j2) continue;
				
				double c_chi2 = getChi2(jets.at(b), jets.at(j1), jets.at(j2),
						sigma_jets.at(b), sigma_jets.at(j1), sigma_jets.at(j2));
				
				if (c_chi2 < min_chi2) min_chi2 = c_chi2;

			}
*/

	double c_chi2 =  getChi2(jets.at(bi), jets.at(j1), jets.at(j2), 
		sigma_jets.at(bi), sigma_jets.at(j1), sigma_jets.at(j2) ); 
	return c_chi2;
}

// This function calculates the hadronic chi2 - SNT version
double calculateChi2SNT(vector<LorentzVector>& jets, vector<float>& sigma_jets, vector<float>& btag){

  assert(jets.size() == sigma_jets.size());
  assert(jets.size() == btag.size());

  //check at most first 6 jets
  int n_jets = jets.size();
  if (n_jets>6) n_jets = 6;
  //consider at least 3 jets
  if (n_jets<3) return 999999.;
  
  vector<int> v_i, v_j;
  vector<double> v_k1, v_k2;
  for ( int i=0; i<n_jets; ++i )
    for ( int j=i+1; j<n_jets; ++j ){

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
      if (c1!=c1) {
        cout<<"[PartonCombinatorics::recoHadronicTop] ERROR: c1 parameter is NAN! Skipping this parton combination"
	    <<endl;
        continue;
      }
      double c2 = fc2(c1, jets[i].mass2(), jets[j].mass2(), hadW.mass2());

      delete minimizer;


      //     * W Mass check :)
      //     *  Never trust a computer you can't throw out a window.
      //      *  - Steve Wozniak

      // cout << "c1 = " <<  c1 << "  c1 = " << c2 << "   M_jj = "
      // 	   << ((jets[i] * c1) + (jets[j] * c2)).mass() << endl;

      v_i.push_back(i);
      v_j.push_back(j);
      v_k1.push_back(c1);
      v_k2.push_back(c2);
    }
  
  //Apply b-consistency requirement
  int n_btag = 0;
  for( int i = 0 ; i < n_jets ; i++ )
    if( btag.at(i) > BTAG_MED ) n_btag++;

  double chi2min = 99999.;

  //consider b-jet in leading 3 jets
  for ( int b=0; b<n_jets; ++b ) {    

    //if not tagged, consider only 3 leading jets
    if( btag.at(b)<BTAG_MED && b>2 ) continue;

    //require b-tagging if have more than 1 b-tag
    if( n_btag>1 && btag.at(b) < BTAG_MED ) continue;
    double pt_b = jets[b].Pt();
      
    for (unsigned int w = 0; w < v_i.size() ; ++w ) {
      int i = v_i[w];
      int j = v_j[w];
      if ( i==b || j==b ) continue;
      //count number of b-tagged Ws
      int nwb = 0;
      if (btag.at(i) > BTAG_MED) nwb++;
      if (btag.at(j) > BTAG_MED) nwb++;
      //no btagged jets in W if have few btags
      if ( n_btag<3  && nwb>0 ) continue;
      //In 3 b-tag case, allow for 1 W jet to be tagged
      // If have more b-tags then btagging information not useful
      if ( n_btag==3 && nwb>1 ) continue;
 
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
      double sigma_w2 = pow(pt_w1*sigma_jets[i], 2)
	+ pow(pt_w2*sigma_jets[j], 2);
      double smw2 = (1. + 2.*pow(pt_w,2)/pow(massW,2))*sigma_w2;
      double pt_t = hadT.Pt();
      double sigma_t2 = pow(c1*pt_w1*sigma_jets[i],2)
	+ pow(c2*pt_w2*sigma_jets[j],2)
	+ pow(pt_b*sigma_jets[b],2);
      double smtop2 = (1. + 2.*pow(pt_t,2)/pow(massT,2))*sigma_t2;
      
      double c_chi2 = pow(massT-PDG_TOP_MASS, 2)/smtop2
	+ pow(massW-PDG_W_MASS, 2)/smw2;
      if (c_chi2<chi2min) chi2min = c_chi2;

    }
  }

  return chi2min;
}


double getChi2(LorentzVector& jets_b, LorentzVector& jets_j1, LorentzVector& jets_j2,
		float sigma_b, float sigma_j1, float sigma_j2){

	LorentzVector hadW = jets_j1 + jets_j2;
	LorentzVector hadT = jets_b + hadW;

	double massT = hadT.mass();
	double massW = hadW.mass();

	double pt_w1 = jets_j1.Pt();
	double pt_w2 = jets_j2.Pt();
	double pt_b = jets_b.Pt();

	double pt_w = hadW.Pt();

	double sigma_w2 = pt_w1*sigma_j1 * pt_w1*sigma_j1
		+ pt_w2*sigma_j2 * pt_w2*sigma_j2;

	double smw2 = (1.+2.*pt_w*pt_w/massW/massW)*sigma_w2;

	double pt_t = hadT.Pt();

	double sigma_t2 = pt_w1*sigma_j1 * pt_w1*sigma_j1
		+ pt_w2*sigma_j2 * pt_w2*sigma_j2
		+ pt_b*sigma_b * pt_b*sigma_b;

	double smtop2 = (1.+2.*pt_t*pt_t/massT/massT)*sigma_t2;

	double c_chi2 = (massT-PDG_TOP_MASS)*(massT-PDG_TOP_MASS)/smtop2
		+ (massW-PDG_W_MASS)*(massW-PDG_W_MASS)/smw2;

	return c_chi2;
}


int getRegionNumber(float mstop, float mlsp) {
    float delta = mstop - mlsp;

    if ( delta < 173.5 ) return 5;
    if ( mstop < 350 ) return 1;
    if ( delta < 300 ) return 2;
    if ( delta < 500 ) return 3;

    return 4;
}

//--------------------------------------------------------------------

double TopPtWeight(double topPt){
  if( topPt<0 ) return 1;

  double p0 = 1.18246e+00;
  double p1 = 4.63312e+02;
  double p2 = 2.10061e-06;

  if( topPt>p1 ) topPt = p1;

  double result = p0 + p2 * topPt * ( topPt - 2 * p1 );
  return result;
}

//--------------------------------------------------------------------
//updated version of top pt reweighting function from June 2013
//see slide 12 of presentation: https://indico.cern.ch/getFile.py/access?contribId=2&resId=0&materialId=slides&confId=252018

double TopPtWeight_v2(double topPt){
  if( topPt<0 ) return 1.;
  double result = exp(0.156 - 0.00137 * topPt);
  return result;
}

//--------------------------------------------------------------------

// btag POG recommended SFs from Moriond 2013 recommendations, located here:
//  https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagPOG#Recommendation_for_b_c_tagging_a

// in particular, using (flat) ttbar results for b and c
// using functional form for light from code linked to twiki above
//  https://twiki.cern.ch/twiki/pub/CMS/BtagPOG/SFlightFuncs_Moriond2013.C

// !!! for CVSM only!
float getBtagSF(float pt, float eta, int id){

  // same SF for b and c
  if ((abs(id) == 4) || (abs(id) == 5)) return 0.963;

  // otherwise use functional forms for light
  if( fabs(eta) < 0.8 )
    return ((1.06238+(0.00198635*pt))+(-4.89082e-06*(pt*pt)))+(3.29312e-09*(pt*(pt*pt)));

  if( fabs(eta) >= 0.8 && fabs(eta) < 1.6 )
    return ((1.08048+(0.00110831*pt))+(-2.96189e-06*(pt*pt)))+(2.16266e-09*(pt*(pt*pt)));

  if( fabs(eta) >= 1.6 && fabs(eta) < 2.4 )
    return ((1.09145+(0.000687171*pt))+(-2.45054e-06*(pt*pt)))+(1.7844e-09*(pt*(pt*pt)));

  std::cout << __FILE__ << " " << __LINE__ << ":getBtagSF: invalid eta: " << eta << ", or invalid id: " << id << std::endl;
  return 1.;

}

