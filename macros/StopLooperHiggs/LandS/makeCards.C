#include <algorithm>
#include <iostream>
#include <vector>
#include <math.h>
#include <fstream>
#include <sstream>
#include <iomanip>

#include "TCanvas.h"
#include "TLegend.h"
#include "TChain.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TROOT.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"
#include "TPad.h"
#include "TCut.h"
#include "TProfile.h"
#include "THStack.h"
#include "TLatex.h"
#include "TGraphErrors.h"
#include "TStyle.h"
#include "TLine.h"
#include "TMath.h"

using namespace std;

void makeThisCard(ofstream* dofile, char* cardname, int nobs, float nbkg, float sbkg, float seff = 1.00){

  ofstream* ofile = new ofstream();

  ofile->open(Form("%s.txt",cardname));

  *ofile <<  "imax 1  number of channels"                                                    << endl;
  *ofile <<  "jmax 1  number of backgrounds"                                                 << endl;
  *ofile <<  "kmax 2  number of nuisance parameters (sources of systematical uncertainties)" << endl;
  *ofile <<  "------------"                                                                  << endl;
  *ofile <<  "bin         1"                                                                 << endl;
  *ofile <<  Form("observation %i",nobs)                                                     << endl;
  *ofile <<  "------------"                                                                  << endl;
  *ofile <<  "bin             1      1"                                                      << endl;
  *ofile <<  "process       signal  background"                                              << endl; 
  *ofile <<  "process         0      1"                                                      << endl;
  *ofile <<  Form("rate            1    %.1f",nbkg)                                          << endl;
  *ofile <<  "------------"                                                                  << endl;
  *ofile <<  Form("deltaS  lnN   %.2f    -     uncertainty on signal",seff)                  << endl;
  *ofile <<  Form("deltaB  lnN     -   %.2f    uncertainty on background",sbkg)              << endl; 

  ofile->close();

  *dofile << Form("../test/lands.exe -d %s.txt -M Hybrid --freq --ExpectationHints Asymptotic --scanRs 1 --freq --nToysForCLsb 10000 --nToysForCLb 5000 --seed 1234 -n %s",cardname,cardname) << endl;

}

void makeCards(){

  ofstream* dofile = new ofstream();
  
  dofile->open("doLimits.sh");

  makeThisCard(dofile,"all" , 53 , 52.85 , 1.07 , 1.10);
  makeThisCard(dofile,"1l3b" , 12 , 11.79 , 1.13 , 1.10);
  makeThisCard(dofile,"1l4b" , 24 , 24.42 , 1.097 , 1.10);
  makeThisCard(dofile,"2l3b" , 9 , 8.55 , 1.232 , 1.10);
  makeThisCard(dofile,"2l4b" , 8 , 8.08 , 1.412 , 1.10);

}
