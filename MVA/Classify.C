/**********************************************************************************
 * Project   : TMVA - a Root-integrated toolkit for multivariate data analysis    *
 * Package   : TMVA                                                               *
 * Exectuable: TMVAClassificationApplication                                      *
 *                                                                                *
 * This macro provides a simple example on how to use the trained classifiers     *
 * within an analysis module                                                      *
 **********************************************************************************/

#include <cstdlib>
#include <vector>
#include <iostream>
#include <map>
#include <string>

#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TStopwatch.h"
#include "TChain.h"

#include "TMVAGui.C"

#if not defined(__CINT__) || defined(__MAKECINT__)
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"
#endif

using namespace TMVA;

//--------------------------------------------------------------------

void fillUnderOverFlow(TH1F *h1, float value, float weight = 1)
{
  float min = h1->GetXaxis()->GetXmin();
  float max = h1->GetXaxis()->GetXmax();

  if (value > max) value = h1->GetBinCenter(h1->GetNbinsX());
  if (value < min) value = h1->GetBinCenter(1);

  h1->Fill(value, weight);
}

//--------------------------------------------------------------------

void Classify_HWW( TString myMethodList = "" ) 
{   
#ifdef __CINT__
  gROOT->ProcessLine( ".O0" ); // turn off optimization in CINT
#endif

  //--------------------------------------------------------------------
  // path to weights dir (this is where MVA training info is stored)
  // output root file will be stored at [path]/output
  //--------------------------------------------------------------------

  TString path   = "Trainings/v5/H160_WW_10vars_dphi10/";
  //TString path   = "./";

  //-----------------------------------
  // select samples to run over
  //-----------------------------------

  char* babyPath = "/tas/cerati/HtoWWmvaBabies/latest";
  int mH         = 160;  // choose Higgs mass

  vector<char*> samples;
  samples.push_back("WWTo2L2Nu");
  samples.push_back("GluGluToWWTo4L");
  samples.push_back("WZ");
  samples.push_back("ZZ");
  samples.push_back("TTJets");
  samples.push_back("tW");
  samples.push_back("WJetsToLNu");
  samples.push_back("DY");
  //samples.push_back("WJetsFO3");

  if     ( mH == 130 ) samples.push_back("Higgs130");
  else if( mH == 160 ) samples.push_back("Higgs160");
  else if( mH == 200 ) samples.push_back("Higgs200");
  else{
    cout << "Error, unrecognized Higgs mass " << mH << " GeV, quitting" << endl;
    exit(0);
  }

  //--------------------------------------------------------------------------------
  // IMPORTANT: set the following variables to the same set used for MVA training!!!
  //--------------------------------------------------------------------------------
  
  std::map<std::string,int> mvaVar;
  mvaVar[ "lephard_pt" ]        = 1;
  mvaVar[ "lepsoft_pt" ]        = 1;
  mvaVar[ "dil_dphi" ]          = 1;
  mvaVar[ "dil_mass" ]          = 1;
  mvaVar[ "event_type" ]        = 0;
  mvaVar[ "met_projpt" ]        = 1;
  mvaVar[ "met_pt" ]            = 0;
  mvaVar[ "mt_lephardmet" ]     = 1;
  mvaVar[ "mt_lepsoftmet" ]     = 1;
  mvaVar[ "mthiggs" ]           = 1;
  mvaVar[ "dphi_lephardmet" ]   = 1;
  mvaVar[ "dphi_lepsoftmet" ]   = 1;
  mvaVar[ "lepsoft_fbrem" ]     = 0;
  mvaVar[ "lepsoft_eOverPIn" ]  = 0;
  mvaVar[ "lepsoft_qdphi" ]     = 0;

  //---------------------------------------------------------------

  // This loads the library
  TMVA::Tools::Instance();

  // Default MVA methods to be trained + tested
  std::map<std::string,int> Use;

  // --- Cut optimisation
  Use["Cuts"]            = 1;
  Use["CutsD"]           = 1;
  Use["CutsPCA"]         = 0;
  Use["CutsGA"]          = 0;
  Use["CutsSA"]          = 0;
  // 
  // --- 1-dimensional likelihood ("naive Bayes estimator")
  Use["Likelihood"]      = 1;
  Use["LikelihoodD"]     = 0; // the "D" extension indicates decorrelated input variables (see option strings)
  Use["LikelihoodPCA"]   = 1; // the "PCA" extension indicates PCA-transformed input variables (see option strings)
  Use["LikelihoodKDE"]   = 0;
  Use["LikelihoodMIX"]   = 0;
  //
  // --- Mutidimensional likelihood and Nearest-Neighbour methods
  Use["PDERS"]           = 1;
  Use["PDERSD"]          = 0;
  Use["PDERSPCA"]        = 0;
  Use["PDEFoam"]         = 1;
  Use["PDEFoamBoost"]    = 0; // uses generalised MVA method boosting
  Use["KNN"]             = 1; // k-nearest neighbour method
  //
  // --- Linear Discriminant Analysis
  Use["LD"]              = 1; // Linear Discriminant identical to Fisher
  Use["Fisher"]          = 0;
  Use["FisherG"]         = 0;
  Use["BoostedFisher"]   = 0; // uses generalised MVA method boosting
  Use["HMatrix"]         = 0;
  //
  // --- Function Discriminant analysis
  Use["FDA_GA"]          = 1; // minimisation of user-defined function using Genetics Algorithm
  Use["FDA_SA"]          = 0;
  Use["FDA_MC"]          = 0;
  Use["FDA_MT"]          = 0;
  Use["FDA_GAMT"]        = 0;
  Use["FDA_MCMT"]        = 0;
  //
  // --- Neural Networks (all are feed-forward Multilayer Perceptrons)
  Use["MLP"]             = 0; // Recommended ANN
  Use["MLPBFGS"]         = 0; // Recommended ANN with optional training method
  Use["MLPBNN"]          = 1; // Recommended ANN with BFGS training method and bayesian regulator
  Use["CFMlpANN"]        = 0; // Depreciated ANN from ALEPH
  Use["TMlpANN"]         = 0; // ROOT's own ANN
  //
  // --- Support Vector Machine 
  Use["SVM"]             = 1;
  // 
  // --- Boosted Decision Trees
  Use["BDT"]             = 1; // uses Adaptive Boost
  Use["BDTG"]            = 0; // uses Gradient Boost
  Use["BDTB"]            = 0; // uses Bagging
  Use["BDTD"]            = 0; // decorrelation + Adaptive Boost
  // 
  // --- Friedman's RuleFit method, ie, an optimised series of cuts ("rules")
  Use["RuleFit"]         = 1;
  // ---------------------------------------------------------------
  Use["Plugin"]          = 0;
  Use["Category"]        = 0;
  Use["SVM_Gauss"]       = 0;
  Use["SVM_Poly"]        = 0;
  Use["SVM_Lin"]         = 0;

  std::cout << std::endl;
  std::cout << "==> Start TMVAClassificationApplication" << std::endl;

  // Select methods (don't look at this code - not of interest)
  if (myMethodList != "") {
    for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) it->second = 0;

    std::vector<TString> mlist = gTools().SplitString( myMethodList, ',' );
    for (UInt_t i=0; i<mlist.size(); i++) {
      std::string regMethod(mlist[i]);

      if (Use.find(regMethod) == Use.end()) {
        std::cout << "Method \"" << regMethod 
                  << "\" not known in TMVA under this name. Choose among the following:" << std::endl;
        for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) {
          std::cout << it->first << " ";
        }
        std::cout << std::endl;
        return;
      }
      Use[regMethod] = 1;
    }
  }

  // --------------------------------------------------------------------------------------------------

  const unsigned int nsamples = samples.size();
  
  for( unsigned int i = 0 ; i < nsamples ; ++i ){

    // --- Create the Reader object

    TMVA::Reader *reader = new TMVA::Reader( "!Color:!Silent" );    

    // Create a set of variables and declare them to the reader
    // - the variable names MUST corresponds in name and type to those given in the weight file(s) used
    //    Float_t var1, var2;
    //    Float_t var3, var4;
    //    reader->AddVariable( "myvar1 := var1+var2", &var1 );
    //    reader->AddVariable( "myvar2 := var1-var2", &var2 );
    //    reader->AddVariable( "var3",                &var3 );
    //    reader->AddVariable( "var4",                &var4 );

    Float_t lephard_pt;
    Float_t lepsoft_pt;
    Float_t dil_dphi;
    Float_t dil_mass;
    Float_t event_type;
    Float_t met_projpt;
    Float_t met_pt;
    Float_t mt_lephardmet;
    Float_t mt_lepsoftmet;
    Float_t mthiggs;
    Float_t dphi_lephardmet;
    Float_t dphi_lepsoftmet;
    Float_t lepsoft_fbrem;
    Float_t lepsoft_eOverPIn;
    Float_t lepsoft_qdphi;

    if( mvaVar["lephard_pt"])       reader->AddVariable( "lephard_pt"                  ,   &lephard_pt        ); 
    if( mvaVar["lepsoft_pt"])       reader->AddVariable( "lepsoft_pt"                  ,   &lepsoft_pt        ); 
    if( mvaVar["dil_dphi"])         reader->AddVariable( "dil_dphi"                    ,   &dil_dphi          ); 
    if( mvaVar["dil_mass"])         reader->AddVariable( "dil_mass"                    ,   &dil_mass          ); 
    if( mvaVar["event_type"])       reader->AddVariable( "event_type"                  ,   &event_type        );
    if( mvaVar["met_projpt"])       reader->AddVariable( "met_projpt"                  ,   &met_pt            );
    if( mvaVar["met_pt"])           reader->AddVariable( "met_pt"                      ,   &met_pt            );
    if( mvaVar["mt_lephardmet"])    reader->AddVariable( "mt_lephardmet"               ,   &mt_lephardmet     );
    if( mvaVar["mt_lepsoftmet"])    reader->AddVariable( "mt_lepsoftmet"               ,   &mt_lepsoftmet     );
    if( mvaVar["mthiggs"])          reader->AddVariable( "mthiggs"                     ,   &mthiggs           );  
    if( mvaVar["dphi_lephardmet"])  reader->AddVariable( "dphi_lephardmet"             ,   &dphi_lephardmet   );
    if( mvaVar["dphi_lepsoftmet"])  reader->AddVariable( "dphi_lepsoftmet"             ,   &dphi_lepsoftmet   );
    if( mvaVar["lepsoft_fbrem"])    reader->AddVariable( "lepsoft_fbrem"               ,   &lepsoft_fbrem     );
    if( mvaVar["lepsoft_eOverPIn"]) reader->AddVariable( "lepsoft_eOverPIn"            ,   &lepsoft_eOverPIn  );
    if( mvaVar["lepsoft_qdphi"])    reader->AddVariable( "lepsoft_q * lepsoft_dPhiIn"  ,   &lepsoft_qdphi     );
 

    // Spectator variables declared in the training have to be added to the reader, too
    //    Float_t spec1,spec2;
    //    reader->AddSpectator( "spec1 := var1*2",   &spec1 );
    //    reader->AddSpectator( "spec2 := var1*3",   &spec2 );

    Float_t Category_cat1, Category_cat2, Category_cat3;
    if (Use["Category"]){
      // Add artificial spectators for distinguishing categories
      //       reader->AddSpectator( "Category_cat1 := var3<=0",             &Category_cat1 );
      //       reader->AddSpectator( "Category_cat2 := (var3>0)&&(var4<0)",  &Category_cat2 );
      //       reader->AddSpectator( "Category_cat3 := (var3>0)&&(var4>=0)", &Category_cat3 );
    }

    // --- Book the MVA methods

    //--------------------------------------------------------------------------------------
    // tell Classify_HWW where to find the weights dir, which contains the trained MVA's. 
    // In this example, the weights dir is located at [path]/[dir]
    // and the output root file is written to [path]/[output]
    //--------------------------------------------------------------------------------------

    TString dir    = path + "weights/";
    TString outdir = path + "output/";
    TString prefix = "TMVAClassification";

    // Book method(s)
    for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) {
      if (it->second) {
        TString methodName = TString(it->first) + TString(" method");
        TString weightfile = dir + prefix + TString("_") + TString(it->first) + TString(".weights.xml");
        reader->BookMVA( methodName, weightfile ); 
      }
    }
   
    // Book output histograms
    UInt_t nbin = 1000;
    TH1F   *histLk(0), *histLkD(0), *histLkPCA(0), *histLkKDE(0), *histLkMIX(0), *histPD(0), *histPDD(0);
    TH1F   *histPDPCA(0), *histPDEFoam(0), *histPDEFoamErr(0), *histPDEFoamSig(0), *histKNN(0), *histHm(0);
    TH1F   *histFi(0), *histFiG(0), *histFiB(0), *histLD(0), *histNn(0),*histNnbfgs(0),*histNnbnn(0);
    TH1F   *histNnC(0), *histNnT(0), *histBdt(0), *histBdtG(0), *histBdtD(0), *histRf(0), *histSVMG(0);
    TH1F   *histSVMP(0), *histSVML(0), *histFDAMT(0), *histFDAGA(0), *histCat(0), *histPBdt(0);

    if (Use["Likelihood"])    histLk      = new TH1F( "MVA_Likelihood",    "MVA_Likelihood",    nbin, -1, 1 );               
    if (Use["LikelihoodD"])   histLkD     = new TH1F( "MVA_LikelihoodD",   "MVA_LikelihoodD",   nbin, -1, 0.9999 );
    if (Use["LikelihoodPCA"]) histLkPCA   = new TH1F( "MVA_LikelihoodPCA", "MVA_LikelihoodPCA", nbin, -1, 1 );
    if (Use["LikelihoodKDE"]) histLkKDE   = new TH1F( "MVA_LikelihoodKDE", "MVA_LikelihoodKDE", nbin,  -0.00001, 0.99999 );
    if (Use["LikelihoodMIX"]) histLkMIX   = new TH1F( "MVA_LikelihoodMIX", "MVA_LikelihoodMIX", nbin,  0, 1 );
    if (Use["PDERS"])         histPD      = new TH1F( "MVA_PDERS",         "MVA_PDERS",         nbin,  0, 1 );
    if (Use["PDERSD"])        histPDD     = new TH1F( "MVA_PDERSD",        "MVA_PDERSD",        nbin,  0, 1 );
    if (Use["PDERSPCA"])      histPDPCA   = new TH1F( "MVA_PDERSPCA",      "MVA_PDERSPCA",      nbin,  0, 1 );
    if (Use["KNN"])           histKNN     = new TH1F( "MVA_KNN",           "MVA_KNN",           nbin,  0, 1 );
    if (Use["HMatrix"])       histHm      = new TH1F( "MVA_HMatrix",       "MVA_HMatrix",       nbin, -0.95, 1.55 );
    if (Use["Fisher"])        histFi      = new TH1F( "MVA_Fisher",        "MVA_Fisher",        nbin, -4, 4 );
    if (Use["FisherG"])       histFiG     = new TH1F( "MVA_FisherG",       "MVA_FisherG",       nbin, -1, 1 );
    if (Use["BoostedFisher"]) histFiB     = new TH1F( "MVA_BoostedFisher", "MVA_BoostedFisher", nbin, -2, 2 );
    if (Use["LD"])            histLD      = new TH1F( "MVA_LD",            "MVA_LD",            nbin, -2, 2 );
    if (Use["MLP"])           histNn      = new TH1F( "MVA_MLP",           "MVA_MLP",           nbin, -1.25, 1.5 );
    if (Use["MLPBFGS"])       histNnbfgs  = new TH1F( "MVA_MLPBFGS",       "MVA_MLPBFGS",       nbin, -1.25, 1.5 );
    if (Use["MLPBNN"])        histNnbnn   = new TH1F( "MVA_MLPBNN",        "MVA_MLPBNN",        nbin, -1.25, 1.5 );
    if (Use["CFMlpANN"])      histNnC     = new TH1F( "MVA_CFMlpANN",      "MVA_CFMlpANN",      nbin,  0, 1 );
    if (Use["TMlpANN"])       histNnT     = new TH1F( "MVA_TMlpANN",       "MVA_TMlpANN",       nbin, -1.3, 1.3 );
    if (Use["BDT"])           histBdt     = new TH1F( "MVA_BDT",           "MVA_BDT",           nbin, -1. , 1. );
    if (Use["BDTD"])          histBdtD    = new TH1F( "MVA_BDTD",          "MVA_BDTD",          nbin, -0.8, 0.8 );
    if (Use["BDTG"])          histBdtG    = new TH1F( "MVA_BDTG",          "MVA_BDTG",          nbin, -1.0, 1.0 );
    if (Use["RuleFit"])       histRf      = new TH1F( "MVA_RuleFit",       "MVA_RuleFit",       nbin, -2.0, 2.0 );
    if (Use["SVM_Gauss"])     histSVMG    = new TH1F( "MVA_SVM_Gauss",     "MVA_SVM_Gauss",     nbin,  0.0, 1.0 );
    if (Use["SVM_Poly"])      histSVMP    = new TH1F( "MVA_SVM_Poly",      "MVA_SVM_Poly",      nbin,  0.0, 1.0 );
    if (Use["SVM_Lin"])       histSVML    = new TH1F( "MVA_SVM_Lin",       "MVA_SVM_Lin",       nbin,  0.0, 1.0 );
    if (Use["FDA_MT"])        histFDAMT   = new TH1F( "MVA_FDA_MT",        "MVA_FDA_MT",        nbin, -2.0, 3.0 );
    if (Use["FDA_GA"])        histFDAGA   = new TH1F( "MVA_FDA_GA",        "MVA_FDA_GA",        nbin, -2.0, 3.0 );
    if (Use["Category"])      histCat     = new TH1F( "MVA_Category",      "MVA_Category",      nbin, -2., 2. );
    if (Use["Plugin"])        histPBdt    = new TH1F( "MVA_PBDT",          "MVA_BDT",           nbin, -0.8, 0.8 );

    if (Use["Likelihood"])    histLk      ->Sumw2();
    if (Use["LikelihoodD"])   histLkD     ->Sumw2();
    if (Use["LikelihoodPCA"]) histLkPCA   ->Sumw2();
    if (Use["LikelihoodKDE"]) histLkKDE   ->Sumw2();
    if (Use["LikelihoodMIX"]) histLkMIX   ->Sumw2();
    if (Use["PDERS"])         histPD      ->Sumw2();
    if (Use["PDERSD"])        histPDD     ->Sumw2();
    if (Use["PDERSPCA"])      histPDPCA   ->Sumw2();
    if (Use["KNN"])           histKNN     ->Sumw2();
    if (Use["HMatrix"])       histHm      ->Sumw2();
    if (Use["Fisher"])        histFi      ->Sumw2();
    if (Use["FisherG"])       histFiG     ->Sumw2();
    if (Use["BoostedFisher"]) histFiB     ->Sumw2();
    if (Use["LD"])            histLD      ->Sumw2();
    if (Use["MLP"])           histNn      ->Sumw2();
    if (Use["MLPBFGS"])       histNnbfgs  ->Sumw2();
    if (Use["MLPBNN"])        histNnbnn   ->Sumw2();
    if (Use["CFMlpANN"])      histNnC     ->Sumw2();
    if (Use["TMlpANN"])       histNnT     ->Sumw2();
    if (Use["BDT"])           histBdt     ->Sumw2();
    if (Use["BDTD"])          histBdtD    ->Sumw2();
    if (Use["BDTG"])          histBdtG    ->Sumw2();
    if (Use["RuleFit"])       histRf      ->Sumw2();
    if (Use["SVM_Gauss"])     histSVMG    ->Sumw2();
    if (Use["SVM_Poly"])      histSVMP    ->Sumw2();
    if (Use["SVM_Lin"])       histSVML    ->Sumw2();
    if (Use["FDA_MT"])        histFDAMT   ->Sumw2();
    if (Use["FDA_GA"])        histFDAGA   ->Sumw2();
    if (Use["Category"])      histCat     ->Sumw2();
    if (Use["Plugin"])        histPBdt    ->Sumw2();

    // PDEFoam also returns per-event error, fill in histogram, and also fill significance
    if (Use["PDEFoam"]) {
      histPDEFoam    = new TH1F( "MVA_PDEFoam",       "MVA_PDEFoam",              nbin,  0, 1 );
      histPDEFoamErr = new TH1F( "MVA_PDEFoamErr",    "MVA_PDEFoam error",        nbin,  0, 1 );
      histPDEFoamSig = new TH1F( "MVA_PDEFoamSig",    "MVA_PDEFoam significance", nbin,  0, 10 );
    }

    // Book example histogram for probability (the other methods are done similarly)
    TH1F *probHistFi(0), *rarityHistFi(0);
    if (Use["Fisher"]) {
      probHistFi   = new TH1F( "MVA_Fisher_Proba",  "MVA_Fisher_Proba",  nbin, 0, 1 );
      rarityHistFi = new TH1F( "MVA_Fisher_Rarity", "MVA_Fisher_Rarity", nbin, 0, 1 );
    }

    // Prepare input tree (this must be replaced by your data source)
    // in this example, there is a toy tree with signal and one with background events
    // we'll later on use only the "signal" events for the test in this example.
    //   

 
    TChain *ch = new TChain("Events");

    if( strcmp( samples.at(i) , "DY" ) == 0 ){
      ch -> Add( Form("%s/DYToMuMuM20_PU_testFinal_baby.root",babyPath) );
      ch -> Add( Form("%s/DYToMuMuM10To20_PU_testFinal_baby.root",babyPath) );
      ch -> Add( Form("%s/DYToEEM20_PU_testFinal_baby.root",babyPath) );
      ch -> Add( Form("%s/DYToEEM10To20_PU_testFinal_baby.root",babyPath) );
      ch -> Add( Form("%s/DYToTauTauM20_PU_testFinal_baby.root",babyPath) );
      ch -> Add( Form("%s/DYToTauTauM10To20_PU_testFinal_baby.root",babyPath) );
    }
    if( strcmp( samples.at(i) , "WJetsFO3" ) == 0 ){
      ch -> Add( Form("%s/WJetsToLNu_FOv3_PU_testFinal_baby.root",babyPath) );
      ch -> Add( Form("%s/WToLNu_FOv3_testFinal_baby.root",babyPath) );
    }
    else if( strcmp( samples.at(i) , "Higgs130" ) == 0 ){
      ch -> Add( Form("%s/HToWWTo2L2NuM130_PU_testFinal_baby.root",babyPath) );
      ch -> Add( Form("%s/HToWWToLNuTauNuM130_PU_testFinal_baby.root",babyPath) );
      ch -> Add( Form("%s/HToWWTo2Tau2NuM130_PU_testFinal_baby.root",babyPath) );
    }
    else if( strcmp( samples.at(i) , "Higgs160" ) == 0 ){
      ch -> Add( Form("%s/HToWWTo2L2NuM160_PU_testFinal_baby.root",babyPath) );
      ch -> Add( Form("%s/HToWWToLNuTauNuM160_PU_testFinal_baby.root",babyPath) );
      ch -> Add( Form("%s/HToWWTo2Tau2NuM160_PU_testFinal_baby.root",babyPath) );
    }
    else if( strcmp( samples.at(i) , "Higgs200" ) == 0 ){
      ch -> Add( Form("%s/HToWWTo2L2NuM200_PU_testFinal_baby.root",babyPath) );
      ch -> Add( Form("%s/HToWWToLNuTauNuM200_PU_testFinal_baby.root",babyPath) );
      ch -> Add( Form("%s/HToWWTo2Tau2NuM200_PU_testFinal_baby.root",babyPath) );
    }
    else{
      ch -> Add( Form("%s/%s_PU_testFinal_baby.root",babyPath,samples.at(i)) );
    }

    // --- Event loop

    // Prepare the event tree
    // - here the variable names have to corresponds to your tree
    // - you can use the same variables as above which is slightly faster,
    //   but of course you can use different ones and copy the values inside the event loop
    //
  
    TTree *theTree     = (TTree*) ch;

    std::cout << "--- Using input files: -------------------" <<  std::endl;

    TObjArray *listOfFiles = ch->GetListOfFiles();
    TIter fileIter(listOfFiles);
    TChainElement* currentFile = 0;
    
    while((currentFile = (TChainElement*)fileIter.Next())) {
      std::cout << currentFile->GetTitle() << std::endl;
    }

    Float_t lephard_pt_;
    Float_t lepsoft_pt_;
    Float_t lepsoft_fr_;
    Float_t dil_dphi_;
    Float_t dil_mass_;
    Float_t event_type_;
    Float_t met_projpt_;
    Int_t   jets_num_;
    Int_t   extralep_num_;
    Int_t   lowptbtags_num_;
    Int_t   softmu_num_;
    Float_t event_scale1fb_;
    Float_t met_pt_;
    Int_t   lepsoft_passTighterId_;
    Float_t mt_lephardmet_;
    Float_t mt_lepsoftmet_;
    Float_t mthiggs_;
    Float_t dphi_lephardmet_;
    Float_t dphi_lepsoftmet_;
    Float_t lepsoft_fbrem_;
    Float_t lepsoft_eOverPIn_;
    Float_t lepsoft_q_;
    Float_t lepsoft_dPhiIn_;

    theTree->SetBranchAddress( "lephard_pt_"             ,   &lephard_pt_              ); 
    theTree->SetBranchAddress( "lepsoft_pt_"             ,   &lepsoft_pt_              ); 
    theTree->SetBranchAddress( "lepsoft_fr_"             ,   &lepsoft_fr_              ); 
    theTree->SetBranchAddress( "dil_dphi_"               ,   &dil_dphi_                ); 
    theTree->SetBranchAddress( "dil_mass_"               ,   &dil_mass_                ); 
    theTree->SetBranchAddress( "event_type_"             ,   &event_type_              ); 
    theTree->SetBranchAddress( "met_projpt_"             ,   &met_projpt_              ); 
    theTree->SetBranchAddress( "jets_num_"               ,   &jets_num_                ); 
    theTree->SetBranchAddress( "extralep_num_"           ,   &extralep_num_            ); 
    theTree->SetBranchAddress( "lowptbtags_num_"         ,   &lowptbtags_num_          ); 
    theTree->SetBranchAddress( "softmu_num_"             ,   &softmu_num_              ); 
    theTree->SetBranchAddress( "event_scale1fb_"         ,   &event_scale1fb_          ); 
    theTree->SetBranchAddress( "lepsoft_passTighterId_"  ,   &lepsoft_passTighterId_   );
    theTree->SetBranchAddress( "met_pt_"                 ,   &met_pt_                  );
    theTree->SetBranchAddress( "mt_lephardmet_"          ,   &mt_lephardmet_           );
    theTree->SetBranchAddress( "mt_lepsoftmet_"          ,   &mt_lepsoftmet_           );
    theTree->SetBranchAddress( "mthiggs_"                ,   &mthiggs_                 );
    theTree->SetBranchAddress( "dphi_lephardmet_"        ,   &dphi_lephardmet_         );
    theTree->SetBranchAddress( "dphi_lepsoftmet_"        ,   &dphi_lepsoftmet_         );
    theTree->SetBranchAddress( "lepsoft_fbrem_"          ,   &lepsoft_fbrem_           );
    theTree->SetBranchAddress( "lepsoft_eOverPIn_"       ,   &lepsoft_eOverPIn_        );
    theTree->SetBranchAddress( "lepsoft_q_"              ,   &lepsoft_q_               );
    theTree->SetBranchAddress( "lepsoft_dPhiIn_"         ,   &lepsoft_dPhiIn_          );

    // Efficiency calculator for cut method
    Int_t    nSelCutsGA = 0;
    Double_t effS       = 0.7;

    std::vector<Float_t> vecVar(4); // vector for EvaluateMVA tests

    std::cout << "--- Processing: " << theTree->GetEntries() << " events" << std::endl;
    TStopwatch sw;
    sw.Start();

    int npass   = 0;
    float yield = 0.;
    
    for (Long64_t ievt=0; ievt<theTree->GetEntries();ievt++) {

      if (ievt%1000 == 0) std::cout << "--- ... Processing event: " << ievt << std::endl;

      theTree->GetEntry(ievt);

      //-------------------------------------------------------
      // event selection
      //-------------------------------------------------------

      if( dil_dphi_ > 1. ) continue;

      //em
      if( event_type_ > 0.5 && event_type_ < 2.5 ){
        if( met_projpt_ < 20. )   continue;
      }
      //ee/mm
      if( event_type_ < 0.5 || event_type_ > 2.5 ){
        if( met_projpt_ < 35. )   continue;
      }
      if( lephard_pt_ < 20.           )             continue;
      if( jets_num_ > 0               )             continue;
      if( extralep_num_ > 0           )             continue;
      if( lowptbtags_num_ > 0         )             continue;
      if( softmu_num_ > 0             )             continue;
      if( dil_mass_ < 12.             )             continue;
      if( lepsoft_passTighterId_ == 0 )             continue;
      //if( event_type_ < 1.5    )                    continue;
      //if( event_type > 1.5 && lepsoft_pt_ < 15. )   continue;

      //mH-dependent selection
      if( mH == 130 ){
        if( lepsoft_pt_ < 10.    )                  continue;      
        if( dil_mass_   > 90.    )                  continue;     
      }
      else if( mH == 160 ){
        if( lepsoft_pt_ < 20.    )                  continue;      
        if( dil_mass_   > 100.   )                  continue;     
      }
      else if( mH == 200 ){
        if( lepsoft_pt_ < 20.    )                  continue;      
        if( dil_mass_   > 130.   )                  continue;     
      }

      float weight = event_scale1fb_ * lepsoft_fr_ * 0.5;

      //--------------------------------------------------------
      // important: here we associate branches to MVA variables
      //--------------------------------------------------------

      lephard_pt        = lephard_pt_;
      lepsoft_pt        = lepsoft_pt_;
      dil_mass          = dil_mass_;
      dil_dphi          = dil_dphi_;
      event_type        = event_type_;
      met_pt            = met_pt_;
      met_projpt        = met_projpt_;
      mt_lephardmet     = mt_lephardmet_;
      mt_lepsoftmet     = mt_lepsoftmet_;
      mthiggs           = mthiggs_;
      dphi_lephardmet   = dphi_lephardmet_;
      dphi_lepsoftmet   = dphi_lepsoftmet_;
      lepsoft_fbrem     = lepsoft_fbrem_;
      lepsoft_eOverPIn  = lepsoft_eOverPIn_;
      lepsoft_qdphi     = lepsoft_q_ * lepsoft_dPhiIn_;

      npass++;
      yield+=weight;

      //       var1 = userVar1 + userVar2;
      //       var2 = userVar1 - userVar2;

      // --- Return the MVA outputs and fill into histograms

      if (Use["CutsGA"]) {
        // Cuts is a special case: give the desired signal efficienciy
        Bool_t passed = reader->EvaluateMVA( "CutsGA method", effS );
        if (passed) nSelCutsGA++;
      }

      if (Use["Likelihood"   ])   histLk     ->Fill( reader->EvaluateMVA( "Likelihood method"    ) , weight);
      if (Use["LikelihoodD"  ])   histLkD    ->Fill( reader->EvaluateMVA( "LikelihoodD method"   ) , weight);
      if (Use["LikelihoodPCA"])   histLkPCA  ->Fill( reader->EvaluateMVA( "LikelihoodPCA method" ) , weight);
      if (Use["LikelihoodKDE"])   histLkKDE  ->Fill( reader->EvaluateMVA( "LikelihoodKDE method" ) , weight);
      if (Use["LikelihoodMIX"])   histLkMIX  ->Fill( reader->EvaluateMVA( "LikelihoodMIX method" ) , weight);
      if (Use["PDERS"        ])   histPD     ->Fill( reader->EvaluateMVA( "PDERS method"         ) , weight);
      if (Use["PDERSD"       ])   histPDD    ->Fill( reader->EvaluateMVA( "PDERSD method"        ) , weight);
      if (Use["PDERSPCA"     ])   histPDPCA  ->Fill( reader->EvaluateMVA( "PDERSPCA method"      ) , weight);
      if (Use["KNN"          ])   histKNN    ->Fill( reader->EvaluateMVA( "KNN method"           ) , weight);
      if (Use["HMatrix"      ])   histHm     ->Fill( reader->EvaluateMVA( "HMatrix method"       ) , weight);
      if (Use["Fisher"       ])   histFi     ->Fill( reader->EvaluateMVA( "Fisher method"        ) , weight);
      if (Use["FisherG"      ])   histFiG    ->Fill( reader->EvaluateMVA( "FisherG method"       ) , weight);
      if (Use["BoostedFisher"])   histFiB    ->Fill( reader->EvaluateMVA( "BoostedFisher method" ) , weight);
      if (Use["LD"           ])   histLD     ->Fill( reader->EvaluateMVA( "LD method"            ) , weight);
      if (Use["MLP"          ])   histNn     ->Fill( reader->EvaluateMVA( "MLP method"           ) , weight);
      if (Use["MLPBFGS"      ])   histNnbfgs ->Fill( reader->EvaluateMVA( "MLPBFGS method"       ) , weight);
      if (Use["MLPBNN"       ])   histNnbnn  ->Fill( reader->EvaluateMVA( "MLPBNN method"        ) , weight);
      if (Use["CFMlpANN"     ])   histNnC    ->Fill( reader->EvaluateMVA( "CFMlpANN method"      ) , weight);
      if (Use["TMlpANN"      ])   histNnT    ->Fill( reader->EvaluateMVA( "TMlpANN method"       ) , weight);
      if (Use["BDT"          ])   histBdt    ->Fill( reader->EvaluateMVA( "BDT method"           ) , weight);
      if (Use["BDTD"         ])   histBdtD   ->Fill( reader->EvaluateMVA( "BDTD method"          ) , weight);
      if (Use["BDTG"         ])   histBdtG   ->Fill( reader->EvaluateMVA( "BDTG method"          ) , weight);
      if (Use["RuleFit"      ])   histRf     ->Fill( reader->EvaluateMVA( "RuleFit method"       ) , weight);
      if (Use["SVM_Gauss"    ])   histSVMG   ->Fill( reader->EvaluateMVA( "SVM_Gauss method"     ) , weight);
      if (Use["SVM_Poly"     ])   histSVMP   ->Fill( reader->EvaluateMVA( "SVM_Poly method"      ) , weight);
      if (Use["SVM_Lin"      ])   histSVML   ->Fill( reader->EvaluateMVA( "SVM_Lin method"       ) , weight);
      if (Use["FDA_MT"       ])   histFDAMT  ->Fill( reader->EvaluateMVA( "FDA_MT method"        ) , weight);
      if (Use["FDA_GA"       ])   histFDAGA  ->Fill( reader->EvaluateMVA( "FDA_GA method"        ) , weight);
      if (Use["Category"     ])   histCat    ->Fill( reader->EvaluateMVA( "Category method"      ) , weight);
      if (Use["Plugin"       ])   histPBdt   ->Fill( reader->EvaluateMVA( "P_BDT method"         ) , weight);

      // Retrieve also per-event error
      if (Use["PDEFoam"]) {
        Double_t val = reader->EvaluateMVA( "PDEFoam method" );
        Double_t err = reader->GetMVAError();
        histPDEFoam   ->Fill( val );
        histPDEFoamErr->Fill( err );         
        if (err>1.e-50) histPDEFoamSig->Fill( val/err , weight);
      }         

      // Retrieve probability instead of MVA output
      if (Use["Fisher"])   {
        probHistFi  ->Fill( reader->GetProba ( "Fisher method" ) , weight);
        rarityHistFi->Fill( reader->GetRarity( "Fisher method" ) , weight);
      }
    }

    std::cout << npass << " events passing selection, yield " << yield << std::endl;
 
    // Get elapsed time
    sw.Stop();
    std::cout << "--- End of event loop: "; sw.Print();

    // Get efficiency for cuts classifier
    if (Use["CutsGA"]) std::cout << "--- Efficiency for CutsGA method: " << double(nSelCutsGA)/theTree->GetEntries()
                                 << " (for a required signal efficiency of " << effS << ")" << std::endl;

    if (Use["CutsGA"]) {

      // test: retrieve cuts for particular signal efficiency
      // CINT ignores dynamic_casts so we have to use a cuts-secific Reader function to acces the pointer  
      TMVA::MethodCuts* mcuts = reader->FindCutsMVA( "CutsGA method" ) ;

      if (mcuts) {      
        std::vector<Double_t> cutsMin;
        std::vector<Double_t> cutsMax;
        mcuts->GetCuts( 0.7, cutsMin, cutsMax );
        std::cout << "--- -------------------------------------------------------------" << std::endl;
        std::cout << "--- Retrieve cut values for signal efficiency of 0.7 from Reader" << std::endl;
        for (UInt_t ivar=0; ivar<cutsMin.size(); ivar++) {
          std::cout << "... Cut: " 
                    << cutsMin[ivar] 
                    << " < \"" 
                    << mcuts->GetInputVar(ivar)
                    << "\" <= " 
                    << cutsMax[ivar] << std::endl;
        }
        std::cout << "--- -------------------------------------------------------------" << std::endl;
      }
    }

    // --- Write histograms
    cout << "dir " << dir << endl;
    char* mydir = outdir;
    TFile *target  = new TFile( Form("%s/%s.root",mydir,samples.at(i) ) ,"RECREATE" );
    cout << "Writing to file " << Form("%s/%s.root",mydir,samples.at(i) ) << endl;

    if (Use["Likelihood"   ])   histLk     ->Write();
    if (Use["LikelihoodD"  ])   histLkD    ->Write();
    if (Use["LikelihoodPCA"])   histLkPCA  ->Write();
    if (Use["LikelihoodKDE"])   histLkKDE  ->Write();
    if (Use["LikelihoodMIX"])   histLkMIX  ->Write();
    if (Use["PDERS"        ])   histPD     ->Write();
    if (Use["PDERSD"       ])   histPDD    ->Write();
    if (Use["PDERSPCA"     ])   histPDPCA  ->Write();
    if (Use["KNN"          ])   histKNN    ->Write();
    if (Use["HMatrix"      ])   histHm     ->Write();
    if (Use["Fisher"       ])   histFi     ->Write();
    if (Use["FisherG"      ])   histFiG    ->Write();
    if (Use["BoostedFisher"])   histFiB    ->Write();
    if (Use["LD"           ])   histLD     ->Write();
    if (Use["MLP"          ])   histNn     ->Write();
    if (Use["MLPBFGS"      ])   histNnbfgs ->Write();
    if (Use["MLPBNN"       ])   histNnbnn  ->Write();
    if (Use["CFMlpANN"     ])   histNnC    ->Write();
    if (Use["TMlpANN"      ])   histNnT    ->Write();
    if (Use["BDT"          ])   histBdt    ->Write();
    if (Use["BDTD"         ])   histBdtD   ->Write();
    if (Use["BDTG"         ])   histBdtG   ->Write(); 
    if (Use["RuleFit"      ])   histRf     ->Write();
    if (Use["SVM_Gauss"    ])   histSVMG   ->Write();
    if (Use["SVM_Poly"     ])   histSVMP   ->Write();
    if (Use["SVM_Lin"      ])   histSVML   ->Write();
    if (Use["FDA_MT"       ])   histFDAMT  ->Write();
    if (Use["FDA_GA"       ])   histFDAGA  ->Write();
    if (Use["Category"     ])   histCat    ->Write();
    if (Use["Plugin"       ])   histPBdt   ->Write();

    // Write also error and significance histos
    if (Use["PDEFoam"]) { histPDEFoam->Write(); histPDEFoamErr->Write(); histPDEFoamSig->Write(); }

    // Write also probability hists
    if (Use["Fisher"]) { if (probHistFi != 0) probHistFi->Write(); if (rarityHistFi != 0) rarityHistFi->Write(); }
    target->Close();

    delete reader;
    
    std::cout << "==> TMVAClassificationApplication is done with sample " << samples.at(i) << endl << std::endl;
  } 
}
