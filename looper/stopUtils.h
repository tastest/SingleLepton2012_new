#ifndef stopUtils_h
#define stopUtils_h

#include <algorithm>
#include <iostream>
#include <map>
#include <vector>
#include <set>
#include <sstream>
#include "TChain.h"
#include "TString.h"
#include "TChainElement.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TProfile.h"
#include "TTree.h"
#include "TVector2.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"
#include "TRandom3.h"
#include "Math/LorentzVector.h"



#include "../CORE/MT2/MT2Utility.h"                                                                                                                                                                                  
#include "mt2w_bisect.h"                                                                                                                                                                                             
#include "mt2bl_bisect.h"  

#include "../CORE/CMS2.h"

#include "../CORE/utilities.h"                                                                                                                                                                                       
#include "../CORE/ssSelections.h"                                                                                                                                                                                    
#include "../CORE/electronSelections.h"                                                                                                                                                                              
#include "../CORE/electronSelectionsParameters.h"                                                                                                                                                                    
#include "../CORE/MITConversionUtilities.h"                                                                                                                                                                          
#include "../CORE/muonSelections.h"                                                                                                                                                                                  
#include "../CORE/eventSelections.h"                                                                                                                                                                                 
#include "../CORE/trackSelections.h"                                                                                                                                                                                 
#include "../CORE/metSelections.h"                                                                                                                                                                                   
#include "../CORE/jetcorr/FactorizedJetCorrector.h"                                                                                                                                                                  
#include "../CORE/jetcorr/JetCorrectionUncertainty.h"                                                                                                                                                                
#include "../CORE/jetSelections.h"                                                                                                                                                                                   
#include "../CORE/jetSmearingTools.h"                                                                                                                                                                                
#include "../CORE/photonSelections.h"                                                                                                                                                                                
#include "../CORE/triggerUtils.h"                                                                                                                                                                                    
#include "../CORE/triggerSuperModel.h"                                                                                                                                                                               
#include "../CORE/mcSelections.h"                                                                                                                                                                                    
#include "../CORE/susySelections.h"                                                                                                                                                                                  
#include "../CORE/mcSUSYkfactor.h"                                                                                                                                                                                   
#include "../CORE/SimpleFakeRate.h"              

using namespace tas;
using namespace std;


typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > LorentzVector;
typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > P4;
typedef vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > VofP4;
typedef map<unsigned int, unsigned int> m_uiui;

struct indP4{
  LorentzVector p4obj;
  int p4ind;
};

typedef vector< indP4 > VofiP4;

inline bool sortIP4ByPt(indP4 iP41, indP4 iP42) {
  return iP41.p4obj.pt() > iP42.p4obj.pt();
}

float dz_trk_vtx( const unsigned int trkidx, const unsigned int vtxidx );

//--------------------------------------------------------------------                                                                                                                                                   
pair<float, float> ScaleMET( pair<float, float> p_met, LorentzVector p4_dilep, double rescale = 1.0);

pair<float,float> getPhiCorrMET( float met, float metphi, int nvtx, bool ismc );

pair<float,float> getTrackerMET( P4 *lep, double deltaZCut, bool dolepcorr );

pair<float,float> getTrackerMET( P4 *lep, double deltaZCut = 0.1, bool dolepcorr = true );

pair<float,float> Type1PFMET( VofP4 jets_p4 , vector<float> cors , vector<float> l1cors , float minpt );

float getDataMCRatio(float eta);

pair<float,float> Type1PFMETSmear(JetSmearer* jetSmearer, bool isData,
				  VofP4 jets_p4 , float met, float metphi);

float getMT( float leppt , float lepphi , float met , float metphi );

float trackIso( int thisPf , float coneR , float dz_thresh , bool dovtxcut , float pt_thresh );
float trackIso( int thisPf , float coneR = 0.3 , float dz_thresh = 0.05 , bool dovtxcut = false , float pt_thresh = 0.0); 

//--------------------------------------------------------------------                                                                                                                                                   
double dRbetweenVectors(const LorentzVector &vec1,
                        const LorentzVector &vec2 );

bool isGenBMatched ( LorentzVector p4, float dR );

int isGenQGMatched ( LorentzVector p4, float dR );

float dRGenJet ( LorentzVector p4, float genminpt=20. );

float getminjdr( VofiP4 jets, LorentzVector *particle );

int isSSVMTagged ( int ijet );

int isCSVTagged ( int ijet );

bool isBTagged ( LorentzVector p4, VofP4 bJets );

int getLeptonMatchIndex ( LorentzVector *jet, LorentzVector *lep1, LorentzVector *lep2, float dR );

int findTriggerIndex(TString trigName);

bool objectPassTrigger(const LorentzVector &obj, const std::vector<LorentzVector> &trigObjs, float pt);

TString triggerName(TString triggerPattern);

bool objectPassTrigger(const LorentzVector &obj, char* trigname, float drmax ); 

#endif
