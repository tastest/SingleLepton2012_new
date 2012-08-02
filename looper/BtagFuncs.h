#ifndef btagFuncs_h
#define btagFuncs_h

#include"TF1.h"

float getBtagSF(float jetpt, float jeteta, string tagger) {

  //Functions for b-tagging SFs obtained from 
  //https://twiki.cern.ch/twiki/pub/CMS/BtagPOG/MistagFuncs.C
  //SFb is provided for each tagger/op in the pt range 30-670 GeV and abs(eta) range 0-2.4 
  //for pt > 670 GeV: use the SFb value at 670 GeV with twice the quoted uncertainty 
  //for pt < 30 GeV: use the SFb value at 30 GeV with a +_0.12 absolute uncertainty 
  //i.e SFb(pt<30) = SFb(pt=30) +_ 0.12, so the relative uncertainty is 0.12/SFb(pt=30)

  //no correction outside acceptance
  if ( abs(jeteta)>2.4 ) 
    return 1.;
  
  float pt = jetpt; 
  if ( jetpt>670 ) pt = 670.;
  if ( jetpt<30  ) pt = 30.;

  //Tagger: SSVHEM within 30 < pt < 670 GeV, abs(eta) < 2.4, x = pt
  if (tagger=="SSVHEM") {
    return 0.896462*((1.+(0.00957275*pt))/(1.+(0.00837582*pt)));
  }
  // Tagger: TCHEL within 30 < pt < 670 GeV, abs(eta) < 2.4, x = pt
  else if (tagger=="TCHEL") 
    return 0.603913*((1.+(0.286361*pt))/(1.+(0.170474*pt)));
  
  // Tagger: TCHEM within 30 < pt < 670 GeV, abs(eta) < 2.4, x = pt
  else if (tagger=="TCHEM") 
    return 0.932251*((1.+(0.00335634*pt))/(1.+(0.00305994*pt)));
  
  //Tagger: CSVL within 30 < pt < 670 GeV, abs(eta) < 2.4, x = pt
  else if (tagger=="CSVL") 
    return 1.02658*((1.+(0.0195388*pt))/(1.+(0.0209145*pt)));
  
  //Tagger: CSVM within 30 < pt < 670 GeV, abs(eta) < 2.4, x = pt
  else if (tagger=="CSVM") 
    return 0.6981*((1.+(0.414063*pt))/(1.+(0.300155*pt)));
  
  // Tagger: CSVT within 30 < pt < 670 GeV, abs(eta) < 2.4, x = pt
  else if (tagger=="CSVT") 
    return 0.901615*((1.+(0.552628*pt))/(1.+(0.547195*pt)));
  
  //no correction otherwise
  return 1.;
  
}

//--------------------------------------------------------------------

float getBtagEff(float x, string tagger) {
  
  //Functions for b-tagging efficiencies obtained from 
  //https://twiki.cern.ch/twiki/pub/CMS/BtagPOG/eff_b_c-ttbar_payload.txt
  //NOTE x = discriminator value

  // Tagger: SSVHE
  if (tagger=="SSVHE") 
    return -0.00592795459517*x*x*x*x +  0.100276502585*x*x*x +  -0.562858452018*x*x +  1.02676379935*x +  0.0810613917651;

  // Tagger: TCHE
  else if (tagger=="TCHE") 
    return -3.67153247396e-07*x*x*x*x +  -2.81599797034e-05*x*x*x +  0.00293190163243*x*x +  -0.0849600849778*x +  0.928524440715;

  //Tagger: CSV
  else if (tagger=="CSV") 
    return -6.41591823466*x*x*x*x +  11.5812173893*x*x*x +  -6.94406444589*x*x +  1.13278339944*x +  0.889359753365;

  return -999.;
}

//--------------------------------------------------------------------

float getMistagSF(float jetpt, float jeteta, string tagger)
{

  //Functions for mistag efficiency obtained from 
  //https://twiki.cern.ch/twiki/pub/CMS/BtagPOG/SFlightFuncs.C
  //SFlight is provided for each tagger/op in the pt range 20-670 GeV in several abs(eta) bins: 
  //0-0.5, 0.5-1, 1-1.5, 1.5-2.4, 0-2.4 for Loose
  //0-0.8, 0.8-1.6, 1.6-2.4, 0-2.4 for Medium
  //0-2.4 for Tight operating points. 
  //for pt > 670 GeV and 0<abs(eta)<2.4: use the SFlight value at 670 GeV (integrated on all eta) 
  //with twice the quoted uncertainty

  //no correction outside acceptance
  if ( abs(jeteta)>2.4 ) 
    return 1.;

  //no correction for jets below 20 GeV
  if ( jetpt < 20 )
    return 1.;

  // Tagger: SSVHEM
  if (tagger=="SSVHEM") {
    
    //evaluate function averaged over eta at 670
    //((0.890254+(0.000553319*jetpt))+(-1.29993e-06*(jetpt*jetpt)))+(4.19294e-10*(jetpt*(jetpt*jetpt)))
    if ( jetpt>670. )
      return 0.803547274322000016;
      
    //otherwise do the eta dependent corrections
    if ( abs(jeteta) < 0.8 )
      return ((0.86318+(0.000801639*jetpt))+(-1.64119e-06*(jetpt*jetpt)))+(2.59121e-10*(jetpt*(jetpt*jetpt)));
    else if ( abs(jeteta) < 1.6 )
      return ((0.958973+(-0.000269555*jetpt))+(1.381e-06*(jetpt*jetpt)))+(-1.87744e-09*(jetpt*(jetpt*jetpt)));
    else //upto 2.4 in eta
      return ((0.923033+(-0.000898227*jetpt))+(4.74565e-06*(jetpt*jetpt)))+(-6.11053e-09*(jetpt*(jetpt*jetpt)));

  }

  // Tagger: TCHEL
  if (tagger=="TCHEL") {
    
    //evaluate function averaged over eta at 670
    //(1.10649*((1+(-9.00297e-05*jetpt))+(2.32185e-07*(jetpt*jetpt))))+(-4.04925e-10*(jetpt*(jetpt*(jetpt/(1+(-0.00051036*jetpt))))))
    if ( jetpt>670. )
      return 0.970004440844050575;
      
    //otherwise do the eta dependent corrections
    if ( abs(jeteta) < 0.5 )
      return (1.13615*((1+(-0.00119852*jetpt))+(1.17888e-05*(jetpt*jetpt))))+(-9.8581e-08*(jetpt*(jetpt*(jetpt/(1+(0.00689317*jetpt))))));
    else if ( abs(jeteta) < 1.0 )
      return (1.13277*((1+(-0.00084146*jetpt))+(3.80313e-06*(jetpt*jetpt))))+(-8.75061e-09*(jetpt*(jetpt*(jetpt/(1+(0.00118695*jetpt))))));
    else if ( abs(jeteta) < 1.5 )
      return (1.17163*((1+(-0.000828475*jetpt))+(3.0769e-06*(jetpt*jetpt))))+(-4.692e-09*(jetpt*(jetpt*(jetpt/(1+(0.000337759*jetpt))))));
    else //upto 2.4 in eta
      return (1.14554*((1+(-0.000128043*jetpt))+(4.10899e-07*(jetpt*jetpt))))+(-2.07565e-10*(jetpt*(jetpt*(jetpt/(1+(-0.00118618*jetpt))))));

  }

  // Tagger: TCHEM
  if (tagger=="TCHEM") {
    
    //evaluate function averaged over eta at 670
    //(1.06268*((1+(0.00390509*jetpt))+(-5.85405e-05*(jetpt*jetpt))))+(7.87135e-07*(jetpt*(jetpt*(jetpt/(1+(0.01259*jetpt))))))
    if ( jetpt>670. )
      return 1.00809635615323856;
      
    //otherwise do the eta dependent corrections
    if ( abs(jeteta) < 0.8 )
      return (1.2875*((1+(-0.000356371*jetpt))+(1.08081e-07*(jetpt*jetpt))))+(-6.89998e-11*(jetpt*(jetpt*(jetpt/(1+(-0.0012139*jetpt))))));
    else if ( abs(jeteta) < 1.6 )
      return (1.24986*((1+(-0.00039734*jetpt))+(5.37486e-07*(jetpt*jetpt))))+(-1.74023e-10*(jetpt*(jetpt*(jetpt/(1+(-0.00112954*jetpt))))));
    else //upto 2.4 in eta
      return (1.10763*((1+(-0.000105805*jetpt))+(7.11718e-07*(jetpt*jetpt))))+(-5.3001e-10*(jetpt*(jetpt*(jetpt/(1+(-0.000821215*jetpt))))));

  } 

  // Tagger: CSVL
  if (tagger=="CSVL") {

    //evaluate function averaged over eta at 670
    //((1.0344+(0.000962994*jetpt))+(-3.65392e-06*(jetpt*jetpt)))+(3.23525e-09*(jetpt*(jetpt*jetpt))
    if ( jetpt>670. )
      return 1.01240478774999998;
      
    //otherwise do the eta dependent corrections
    if ( abs(jeteta) < 0.5 )
      return ((1.07536+(0.000175506*jetpt))+(-8.63317e-07*(jetpt*jetpt)))+(3.27516e-10*(jetpt*(jetpt*jetpt)));
    else if ( abs(jeteta) < 1.0 )
      return ((1.07846+(0.00032458*jetpt))+(-1.30258e-06*(jetpt*jetpt)))+(8.50608e-10*(jetpt*(jetpt*jetpt)));
    else if ( abs(jeteta) < 1.5 )
      return ((1.08294+(0.000474818*jetpt))+(-1.43857e-06*(jetpt*jetpt)))+(1.13308e-09*(jetpt*(jetpt*jetpt)));
    else //upto 2.4 in eta
      return ((1.0617+(0.000173654*jetpt))+(-5.29009e-07*(jetpt*jetpt)))+(5.55931e-10*(jetpt*(jetpt*jetpt)));

  }

  // Tagger: CSVM
  if (tagger=="CSVM") {

    //evaluate function averaged over eta at 670
    //((1.04318+(0.000848162*jetpt))+(-2.5795e-06*(jetpt*jetpt)))+(1.64156e-09*(jetpt*(jetpt*jetpt)))
    if ( jetpt>670. )
      return 0.947231500280000027;
      
    //otherwise do the eta dependent corrections
    if ( abs(jeteta) < 0.8 )
      return ((1.06182+(0.000617034*jetpt))+(-1.5732e-06*(jetpt*jetpt)))+(3.02909e-10*(jetpt*(jetpt*jetpt)));
    else if ( abs(jeteta) < 1.6 )
      return ((1.111+(-9.64191e-06*jetpt))+(1.80811e-07*(jetpt*jetpt)))+(-5.44868e-10*(jetpt*(jetpt*jetpt)));
    else //upto 2.4 in eta
      return ((1.08498+(-0.000701422*jetpt))+(3.43612e-06*(jetpt*jetpt)))+(-4.11794e-09*(jetpt*(jetpt*jetpt)));

  }

  // Tagger: CSVT
  if (tagger=="CSVT") {

    //no eta dependent corrections available
    //evaluate function averaged over eta at 670
    //((0.948463+(0.00288102*jetpt))+(-7.98091e-06*(jetpt*jetpt)))+(5.50157e-09*(jetpt*(jetpt*jetpt)))
    if ( jetpt>670. )
      return 0.950784598910000112;
    else 
      return ((0.948463+(0.00288102*jetpt))+(-7.98091e-06*(jetpt*jetpt)))+(5.50157e-09*(jetpt*(jetpt*jetpt)));

  }

  //no correction otherwise
  return 1.;  

}

//--------------------------------------------------------------------

float getMistags(float jetpt, float jeteta, string tagger)
{

  //Functions for mistag efficiency obtained from 
  //https://twiki.cern.ch/twiki/pub/CMS/BtagPOG/MistagFuncs.C
  //eff_l is provided for each tagger/op in the pt range 20-670 GeV in several abs(eta) bins: 
  //0-0.5, 0.5-1, 1-1.5, 1.5-2.4, 0-2.4 for Loose 
  //0-0.8, 0.8-1.6, 1.6-2.4, 0-2.4 for Medium 
  //0-2.4 for Tight operating points. 
  //no uncertainty provided yet

  //no efficiency available outside acceptance
  if ( abs(jeteta)>2.4 ) 
    return 0.;
  
  //no efficiency available for jets with too low pT
  if ( jetpt<20  ) 
    return 0.;

  // Tagger: SSVHEM
  if (tagger=="SSVHEM") {

    //evaluate function averaged over eta at 670
    //(((-0.000420178+(0.00029105*jetpt))+(-8.9398e-07*(jetpt*jetpt)))+(1.35401e-09*(jetpt*(jetpt*jetpt))))+(-7.93826e-13*(jetpt*(jetpt*(jetpt*jetpt))))
    if ( jetpt>670. )
      return 0.0405469718405399859;
      
    //otherwise do the eta dependent corrections
    if ( abs(jeteta) < 0.8 )
      return (((0.000547883+(0.00023023*jetpt))+(-7.31792e-07*(jetpt*jetpt)))+(1.15659e-09*(jetpt*(jetpt*jetpt))))+(-7.00641e-13*(jetpt*(jetpt*(jetpt*jetpt))));
    else if ( abs(jeteta) < 1.6 )
      return (((0.000615562+(0.000240254*jetpt))+(-7.00237e-07*(jetpt*jetpt)))+(1.2566e-09*(jetpt*(jetpt*jetpt))))+(-8.59011e-13*(jetpt*(jetpt*(jetpt*jetpt))));
    else //upto 2.4 in eta
      return (((0.000372388+(0.000309735*jetpt))+(-4.35952e-07*(jetpt*jetpt)))+(3.63763e-10*(jetpt*(jetpt*jetpt))))+(-2.11993e-13*(jetpt*(jetpt*(jetpt*jetpt))));

  }

  // Tagger: TCHEL
  if (tagger=="TCHEL") {
    
    //evaluate function averaged over eta at 670
    //(((-0.0276197+(0.00291907*jetpt))+(-7.51594e-06*(jetpt*jetpt)))+(9.82128e-09*(jetpt*(jetpt*jetpt))))+(-5.33759e-12*(jetpt*(jetpt*(jetpt*jetpt))))
    if ( jetpt>670. )
      return 0.432545151256100491;
      
    //otherwise do the eta dependent corrections
    if ( abs(jeteta) < 0.5 )
      return (((-0.0235318+(0.00268868*jetpt))+(-6.47688e-06*(jetpt*jetpt)))+(7.92087e-09*(jetpt*(jetpt*jetpt))))+(-4.06519e-12*(jetpt*(jetpt*(jetpt*jetpt))));
    else if ( abs(jeteta) < 1.0 )
      return (((-0.0257274+(0.00289337*jetpt))+(-7.48879e-06*(jetpt*jetpt)))+(9.84928e-09*(jetpt*(jetpt*jetpt))))+(-5.40844e-12*(jetpt*(jetpt*(jetpt*jetpt))));
    else if ( abs(jeteta) < 1.5 )
      return (((-0.0310046+(0.00307803*jetpt))+(-7.94145e-06*(jetpt*jetpt)))+(1.06889e-08*(jetpt*(jetpt*jetpt))))+(-6.08971e-12*(jetpt*(jetpt*(jetpt*jetpt))));
    else //upto 2.4 in eta
      return (((-0.0274561+(0.00301096*jetpt))+(-8.89588e-06*(jetpt*jetpt)))+(1.40142e-08*(jetpt*(jetpt*jetpt))))+(-8.95723e-12*(jetpt*(jetpt*(jetpt*jetpt))));
    
  }

  // Tagger: TCHEM
  if (tagger=="TCHEM") {

    //evaluate function averaged over eta at 670
    //(-0.00256163+(0.000332759*jetpt))+(-2.39887e-07*(jetpt*jetpt))
    if ( jetpt>670. )
      return 0.112701625700000016;
      
    //otherwise do the eta dependent corrections
    if ( abs(jeteta) < 0.8 )
      return (0.000919586+(0.00026266*jetpt))+(-1.75723e-07*(jetpt*jetpt));
    else if ( abs(jeteta) < 1.6 )
      return (-0.00364137+(0.000350371*jetpt))+(-1.89967e-07*(jetpt*jetpt));
    else //upto 2.4 in eta
      return (-0.00483904+(0.000367751*jetpt))+(-1.36152e-07*(jetpt*jetpt));

  }

  // Tagger: CSVL
  if (tagger=="CSVL") {

    //evaluate function averaged over eta at 670
    //18168.8*(((1+(0.020356*jetpt))+(2.73475e-05*(jetpt*jetpt)))/(1+(5239.42*jetpt)))
    if ( jetpt>670. )
      return 0.139302678480796138;
      
    //otherwise do the eta dependent corrections
    if ( abs(jeteta) < 0.5 )
      return 242534*(((1+(0.0182863*jetpt))+(4.50105e-05*(jetpt*jetpt)))/(1+(108569*jetpt)));
    else if ( abs(jeteta) < 1.0 )
      return 129.938*(((1+(0.0197657*jetpt))+(4.73472e-05*(jetpt*jetpt)))/(1+(55.2415*jetpt)));
    else if ( abs(jeteta) < 1.5 )
      return 592.214*(((1+(0.00671207*jetpt))+(6.46109e-05*(jetpt*jetpt)))/(1+(134.318*jetpt)));
    else //upto 2.4 in eta
      return 93329*(((1+(0.0219705*jetpt))+(3.76566e-05*(jetpt*jetpt)))/(1+(18245.1*jetpt)));

  }

  // Tagger: CSVM
  if (tagger=="CSVM") {
    
    //evaluate function averaged over eta at 670
    //(0.0113428+(5.18983e-05*jetpt))+(-2.59881e-08*(jetpt*jetpt))
    if ( jetpt>670. )
      return 0.0344486029100000007;
      
    //otherwise do the eta dependent corrections
    if ( abs(jeteta) < 0.8 )
      return (0.00967751+(2.54564e-05*jetpt))+(-6.92256e-10*(jetpt*jetpt));
    else if ( abs(jeteta) < 1.6 )
      return (0.00974141+(5.09503e-05*jetpt))+(2.0641e-08*(jetpt*jetpt));
    else //upto 2.4 in eta
      return (0.013595+(0.000104538*jetpt))+(-1.36087e-08*(jetpt*jetpt));

  }

  // Tagger: CSVT
  if (tagger=="CSVT") {
    
    //no eta dependent corrections available
    //evaluate function averaged over eta at 670
    //0.00315116*(((1+(-0.00769281*jetpt))+(2.58066e-05*(jetpt*jetpt)))+(-2.02149e-08*(jetpt*(jetpt*jetpt))))
    if ( jetpt>670. )
      return 0.00425566071163770761;
    else 
      return 0.00315116*(((1+(-0.00769281*jetpt))+(2.58066e-05*(jetpt*jetpt)))+(-2.02149e-08*(jetpt*(jetpt*jetpt))));
  }

  //otherwise no efficiency available here
  return 0.;

}

//--------------------------------------------------------------------

bool getCorrBtag(bool isbtag, int pdgid, float btagSF, float mistagSF, float mistagEff, float rand) {
  //b-tagging efficiency
  //if have a true b and it is tagged, downgrade to untagged based on SF
  if ( pdgid==5 && isbtag  && rand > btagSF ) return false;
  //Correct mistags
  //if have a mistag, downgrade to untagged based on SF
  if (mistagSF<1 && pdgid!=5 && isbtag && rand > mistagSF) return false;
  //if not true b and not tagged, upgrade to tagged based on mistag difference
  if ( mistagSF>1 && pdgid!=5 && !isbtag ) {
    //fraction of jets that need to be upgraded
    float mistagFrac = (1.0 - mistagSF) / (1.0 - (mistagSF/mistagEff) );
    //upgrade to tagged
    if( rand < mistagFrac ) return true;
  } 
  return isbtag;

}

//--------------------------------------------------------------------

#endif
