// @(#)root/tmva $Id: TMVA.C,v 1.1 2012/12/11 11:48:17 benhoob Exp $
/**********************************************************************************
 * Project   : TMVA - a ROOT-integrated toolkit for multivariate data analysis    *
 * Package   : TMVA                                                               *
 * Root Macro: TMVAClassification                                                 *
 *                                                                                *
 * This macro provides examples for the training and testing of the               *
 * TMVA classifiers.                                                              *
 *                                                                                *
 * As input data is used a toy-MC sample consisting of four Gaussian-distributed  *
 * and linearly correlated input variables.                                       *
 *                                                                                *
 * The methods to be used can be switched on and off by means of booleans, or     *
 * via the prompt command, for example:                                           *
 *                                                                                *
 *    root -l ./TMVAClassification.C\(\"Fisher,Likelihood\"\)                     *
 *                                                                                *
 * (note that the backslashes are mandatory)                                      *
 * If no method given, a default set of classifiers is used.                      *
 *                                                                                *
 * The output file "TMVA.root" can be analysed with the use of dedicated          *
 * macros (simply say: root -l <macro.C>), which can be conveniently              *
 * invoked through a GUI that will appear at the end of the run of this macro.    *
 * Launch the GUI via the command:                                                *
 *                                                                                *
 *    root -l ./TMVAGui.C                                                         *
 *                                                                                *
 **********************************************************************************/

#include <cstdlib>
#include <iostream>
#include <map>
#include <string>

#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"

#include "TMVAGui.C"

#if not defined(__CINT__) || defined(__MAKECINT__)
// needs to be included when makecint runs (ACLIC)
#include "TMVA/Factory.h"
#include "TMVA/Tools.h"
#endif


void TMVA_HWW( TString myMethodList = "" )
{
   // The explicit loading of the shared libTMVA is done in TMVAlogon.C, defined in .rootrc
   // if you use your private .rootrc, or run from a different directory, please copy the
   // corresponding lines from .rootrc

   // methods to be processed can be given as an argument; use format:
   //
   // mylinux~> root -l TMVAClassification.C\(\"myMethod1,myMethod2,myMethod3\"\)
   //
   // if you like to use a method via the plugin mechanism, we recommend using
   //
   // mylinux~> root -l TMVAClassification.C\(\"P_myMethod\"\)
   // (an example is given for using the BDT as plugin (see below),
   // but of course the real application is when you write your own
   // method based)


  //-----------------------------------------------------
  // define event selection (store in TCut sel)
  //-----------------------------------------------------
   
  TCut met_projpt = "((event_type<0.5||event_type>2.5)&met_projpt>35)|((event_type>0.5&&event_type<2.5)&met_projpt>20)";
  TCut pt2020      = "lephard_pt > 20 && lepsoft_pt > 20";
  TCut pt2010      = "lephard_pt > 20 && lepsoft_pt > 10";
  TCut jetveto     = "jets_num==0 && extralep_num==0 && lowptbtags_num==0 && softmu_num==0";
  TCut mll12       = "dil_mass > 12.";
  TCut mll90       = "dil_mass < 90.";
  TCut mll100      = "dil_mass < 100.";
  TCut mll130      = "dil_mass < 130.";
  TCut tight_el_ID = "lepsoft_passTighterId==1";
  TCut mutype      = "event_type < 1.5";
  TCut eltype      = "event_type > 1.5";
  TCut elsoft15    = "event_type < 1.5 || lepsoft_pt > 15.";
 
  //------------
  //Higgs130
  //------------
  
  //TCut sel     = pt2010 + met_projpt + jetveto + mll12 + mll90 + tight_el_ID;
  //int  mH      = 130;
 
  //------------
  //Higgs160
  //------------
  
  TCut sel     = pt2020 + met_projpt + jetveto + mll12 + mll100 + "dil_dphi<2.5";
  int  mH      = 160;

  //------------
  //Higgs200
  //------------

  //TCut sel     = pt2020 + met_projpt + jetveto + mll12 + mll130;
  //int  mH      = 200;

  //-----------------------------------------------------
  // choose which variables to include in MVA training
  //-----------------------------------------------------
  
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

  //---------------------------------
  //choose bkg samples to include
  //---------------------------------
  
  const char* babyPath = "/tas/cerati/HtoWWmvaBabies/latest";
  
  TChain *chbackground = new TChain("Events");
  chbackground->Add(Form("%s/WWTo2L2Nu_PU_testFinal_baby.root",babyPath));
  chbackground->Add(Form("%s/GluGluToWWTo4L_PU_testFinal_baby.root",babyPath));
//   chbackground->Add(Form("%s/WZ_PU_testFinal_baby.root",babyPath));
//   chbackground->Add(Form("%s/ZZ_PU_testFinal_baby.root",babyPath));
//   chbackground->Add(Form("%s/TTJets_PU_testFinal_baby.root",babyPath));
//   chbackground->Add(Form("%s/tW_PU_testFinal_baby.root",babyPath));
//   chbackground->Add(Form("%s/WJetsToLNu_PU_testFinal_baby.root",babyPath));
//   chbackground->Add(Form("%s/DYToMuMuM20_PU_testFinal_baby.root",babyPath) );
//   chbackground->Add(Form("%s/DYToMuMuM10To20_PU_testFinal_baby.root",babyPath) );
//   chbackground->Add(Form("%s/DYToEEM20_PU_testFinal_baby.root",babyPath) );
//   chbackground->Add(Form("%s/DYToEEM10To20_PU_testFinal_baby.root",babyPath) );
//   chbackground->Add(Form("%s/DYToTauTauM20_PU_testFinal_baby.root",babyPath) );
//   chbackground->Add(Form("%s/DYToTauTauM10To20_PU_testFinal_baby.root",babyPath) );
  
  //   chbackground->Add(Form("%s/WJetsToLNu_PU_testFinal_baby_FO1.root",babyPath));
  //   chbackground->Add(Form("%s/WJetsToLNu_PU_testFinal_baby_FO2.root",babyPath));
  //   chbackground->Add(Form("%s/WJetsToLNu_PU_testFinal_baby_FO3.root",babyPath));
  //   chbackground->Add(Form("%s/WToLNu_FOv3_testFinal_baby.root",babyPath));
  //   chbackground->Add(Form("%s/WJetsToLNu_PU_testFinal_baby_FO4.root",babyPath));

  //---------------------------------
  //choose signal sample to include
  //---------------------------------

  TChain *chsignal = new TChain("Events");

  if( mH == 130 ){
    chsignal->Add(Form("%s/HToWWTo2L2NuM130_PU_testFinal_baby.root",babyPath));
    chsignal->Add(Form("%s/HToWWToLNuTauNuM130_PU_testFinal_baby.root",babyPath));
    chsignal->Add(Form("%s/HToWWTo2Tau2NuM130_PU_testFinal_baby.root",babyPath));
  }
  else if( mH == 160 ){
    chsignal->Add(Form("%s/HToWWTo2L2NuM160_PU_testFinal_baby.root",babyPath));
    chsignal->Add(Form("%s/HToWWToLNuTauNuM160_PU_testFinal_baby.root",babyPath));
    chsignal->Add(Form("%s/HToWWTo2Tau2NuM160_PU_testFinal_baby.root",babyPath));
  }
  else if( mH == 200 ){
    chsignal->Add(Form("%s/HToWWTo2L2NuM200_PU_testFinal_baby.root",babyPath));
    chsignal->Add(Form("%s/HToWWToLNuTauNuM200_PU_testFinal_baby.root",babyPath));
    chsignal->Add(Form("%s/HToWWTo2Tau2NuM200_PU_testFinal_baby.root",babyPath));
  }
  else{
    std::cout << "Error, unrecognized higgs mass " << mH << " GeV, quitting" << std::endl;
    exit(0);
  }

  //-----------------------------------------------------
  // choose backgrounds to include for multiple outputs
  //-----------------------------------------------------
  
  bool doMultipleOutputs = false;

   TChain *chww = new TChain("Events");
   chww->Add(Form("%s/WWTo2L2Nu_PU_testFinal_baby.root",babyPath));
   chww->Add(Form("%s/GluGluToWWTo4L_PU_testFinal_baby.root",babyPath));

   TChain *chwjets = new TChain("Events");
   chwjets->Add(Form("%s/WJetsToLNu_PU_testFinal_baby.root",babyPath));

   TChain *chtt = new TChain("Events");
   chtt->Add(Form("%s/TTJets_PU_testFinal_baby.root",babyPath));

   std::map<std::string,int> includeBkg;
   includeBkg["ww"]      = 1;
   includeBkg["wjets"]   = 0;
   includeBkg["tt"]      = 0;

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
   //
   // --- multi-output MVA's
   Use["multi_BDTG"]      = 1;
   Use["multi_MLP"]       = 1;
   Use["multi_FDA_GA"]    = 0;
   //
   // ---------------------------------------------------------------



   std::cout << std::endl;
   std::cout << "==> Start TMVAClassification" << std::endl;

   // Select methods (don't look at this code - not of interest)
   if (myMethodList != "") {
      for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) it->second = 0;

      std::vector<TString> mlist = TMVA::gTools().SplitString( myMethodList, ',' );
      for (UInt_t i=0; i<mlist.size(); i++) {
         std::string regMethod(mlist[i]);

         if (Use.find(regMethod) == Use.end()) {
            std::cout << "Method \"" << regMethod << "\" not known in TMVA under this name. Choose among the following:" << std::endl;
            for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) std::cout << it->first << " ";
            std::cout << std::endl;
            return;
         }
         Use[regMethod] = 1;
      }
   }

   // --------------------------------------------------------------------------------------------------

   // --- Here the preparation phase begins

   // Create a ROOT output file where TMVA will store ntuples, histograms, etc.
   TString outfileName( "TMVA_HWW.root" );
   TFile* outputFile = TFile::Open( outfileName, "RECREATE" );

   TString multioutfileName( "TMVA_HWW_multi.root" );
   TFile* multioutputFile;
   if( doMultipleOutputs )
     multioutputFile = TFile::Open( multioutfileName, "RECREATE" );

   // Create the factory object. Later you can choose the methods
   // whose performance you'd like to investigate. The factory is 
   // the only TMVA object you have to interact with
   //
   // The first argument is the base of the name of all the
   // weightfiles in the directory weight/
   //
   // The second argument is the output file for the training results
   // All TMVA output can be suppressed by removing the "!" (not) in
   // front of the "Silent" argument in the option string
   TMVA::Factory *factory = new TMVA::Factory( "TMVAClassification", outputFile,
                                               "!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification" );
   
   TMVA::Factory *multifactory;
   if( doMultipleOutputs )
     multifactory= new TMVA::Factory( "TMVAMulticlass", multioutputFile,
                                      "!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=multiclass" );

   
   // If you wish to modify default settings
   // (please check "src/Config.h" to see all available global options)
   //    (TMVA::gConfig().GetVariablePlotting()).fTimesRMS = 8.0;
   //    (TMVA::gConfig().GetIONames()).fWeightFileDir = "myWeightDirectory";

   // Define the input variables that shall be used for the MVA training
   // note that you may also use variable expressions, such as: "3*var1/var2*abs(var3)"
   // [all types of expressions that can also be parsed by TTree::Draw( "expression" )]
   //factory->AddVariable( "myvar1 := var1+var2", 'F' );
   //factory->AddVariable( "myvar2 := var1-var2", "Expression 2", "", 'F' );
   //factory->AddVariable( "var3",                "Variable 3", "units", 'F' );
   //factory->AddVariable( "var4",                "Variable 4", "units", 'F' );

   //--------------------------------------------------------
   // choose which variables to include in training
   //--------------------------------------------------------
   
   if (mvaVar["lephard_pt"])       factory->AddVariable( "lephard_pt",                 "1st lepton pt",                "GeV", 'F' );
   if (mvaVar["lepsoft_pt"])       factory->AddVariable( "lepsoft_pt",                 "2nd lepton pt",                "GeV", 'F' );
   if (mvaVar["dil_dphi"])         factory->AddVariable( "dil_dphi",                   "dphi(ll)",                     "",    'F' );
   if (mvaVar["dil_mass"])         factory->AddVariable( "dil_mass",                   "M(ll)",                        "GeV", 'F' );
   if (mvaVar["event_type"])       factory->AddVariable( "event_type",                 "Dil Flavor Type",              "",    'F' );
   if (mvaVar["met_projpt"])       factory->AddVariable( "met_projpt",                 "Proj. MET",                    "GeV", 'F' );
   if (mvaVar["met_pt"])           factory->AddVariable( "met_pt",                     "MET",                          "GeV", 'F' );
   if (mvaVar["mt_lephardmet"])    factory->AddVariable( "mt_lephardmet",              "MT(lep1,MET)",                 "GeV", 'F' );
   if (mvaVar["mt_lepsoftmet"])    factory->AddVariable( "mt_lepsoftmet",              "MT(lep2,MET)",                 "GeV", 'F' );
   if (mvaVar["mthiggs"])          factory->AddVariable( "mthiggs",                    "MT(Higgs)",                    "GeV", 'F' );
   if (mvaVar["dphi_lephardmet"])  factory->AddVariable( "dphi_lephardmet",            "dphi(lep1,MET)",               "GeV", 'F' );
   if (mvaVar["dphi_lepsoftmet"])  factory->AddVariable( "dphi_lepsoftmet",            "dphi(lep2,MET)",               "GeV", 'F' );
   if (mvaVar["lepsoft_fbrem"])    factory->AddVariable( "lepsoft_fbrem",              "2nd lepton f_{brem}",          "",    'F' );
   if (mvaVar["lepsoft_eOverPIn"]) factory->AddVariable( "lepsoft_eOverPIn",           "2nd lepton E/p",               "",    'F' );
   if (mvaVar["lepsoft_qdphi"])    factory->AddVariable( "lepsoft_q * lepsoft_dPhiIn", "2nd lepton q#times#Delta#phi", "",    'F' );

   if( doMultipleOutputs ){
     if (mvaVar["lephard_pt"])       multifactory->AddVariable( "lephard_pt",                 "1st lepton pt",                "GeV", 'F' );
     if (mvaVar["lepsoft_pt"])       multifactory->AddVariable( "lepsoft_pt",                 "2nd lepton pt",                "GeV", 'F' );
     if (mvaVar["dil_dphi"])         multifactory->AddVariable( "dil_dphi",                   "dphi(ll)",                     "",    'F' );
     if (mvaVar["dil_mass"])         multifactory->AddVariable( "dil_mass",                   "M(ll)",                        "GeV", 'F' );
     if (mvaVar["event_type"])       multifactory->AddVariable( "event_type",                 "Dil Flavor Type",              "",    'F' );
     if (mvaVar["met_projpt"])       multifactory->AddVariable( "met_projpt",                 "Proj. MET",                    "GeV", 'F' );
     if (mvaVar["met_pt"])           multifactory->AddVariable( "met_pt",                     "MET",                          "GeV", 'F' );
     if (mvaVar["mt_lephardmet"])    multifactory->AddVariable( "mt_lephardmet",              "MT(lep1,MET)",                 "GeV", 'F' );
     if (mvaVar["mt_lepsoftmet"])    multifactory->AddVariable( "mt_lepsoftmet",              "MT(lep2,MET)",                 "GeV", 'F' );
     if (mvaVar["mthiggs"])          multifactory->AddVariable( "mthiggs",                    "MT(Higgs)",                    "GeV", 'F' );
     if (mvaVar["dphi_lephardmet"])  multifactory->AddVariable( "dphi_lephardmet",            "dphi(lep1,MET)",               "GeV", 'F' );
     if (mvaVar["dphi_lepsoftmet"])  multifactory->AddVariable( "dphi_lepsoftmet",            "dphi(lep2,MET)",               "GeV", 'F' );
     if (mvaVar["lepsoft_fbrem"])    multifactory->AddVariable( "lepsoft_fbrem",              "2nd lepton f_{brem}",          "",    'F' );
     if (mvaVar["lepsoft_eOverPIn"]) multifactory->AddVariable( "lepsoft_eOverPIn",           "2nd lepton E/p",               "",    'F' );
     if (mvaVar["lepsoft_qdphi"])    multifactory->AddVariable( "lepsoft_q * lepsoft_dPhiIn", "2nd lepton q#times#Delta#phi", "",    'F' );
   }

   if (mvaVar["lephard_pt"])       cout << "Adding variable to MVA training: lephard_pt"      << endl;
   if (mvaVar["lepsoft_pt"])       cout << "Adding variable to MVA training: lepsoft_pt"      << endl;
   if (mvaVar["dil_dphi"])         cout << "Adding variable to MVA training: dil_dphi"        << endl;
   if (mvaVar["dil_mass"])         cout << "Adding variable to MVA training: dil_mass"        << endl;
   if (mvaVar["event_type"])       cout << "Adding variable to MVA training: event_type"      << endl;
   if (mvaVar["met_projpt"])       cout << "Adding variable to MVA training: met_projpt"      << endl;
   if (mvaVar["met_pt"])           cout << "Adding variable to MVA training: met_pt"          << endl;
   if (mvaVar["mt_lephardmet"])    cout << "Adding variable to MVA training: mt_lephardmet"   << endl;
   if (mvaVar["mt_lepsoftmet"])    cout << "Adding variable to MVA training: mt_lepsoftmet"   << endl;
   if (mvaVar["mthiggs"])          cout << "Adding variable to MVA training: mthiggs"       << endl;
   if (mvaVar["dphi_lephardmet"])  cout << "Adding variable to MVA training: dphi_lephardmet" << endl;
   if (mvaVar["dphi_lepsoftmet"])  cout << "Adding variable to MVA training: dphi_lepsoftmet" << endl;
   if (mvaVar["lepsoft_fbrem"])    cout << "Adding variable to MVA training: lepsoft_fbrem"       << endl;
   if (mvaVar["lepsoft_eOverPIn"]) cout << "Adding variable to MVA training: lepsoft_eOverPIn" << endl;
   if (mvaVar["lepsoft_qdphi"])    cout << "Adding variable to MVA training: lepsoft_qdphi" << endl;

   // You can add so-called "Spectator variables", which are not used in the MVA training,
   // but will appear in the final "TestTree" produced by TMVA. This TestTree will contain the
   // input variables, the response values of all trained MVAs, and the spectator variables
   //factory->AddSpectator( "spec1 := var1*2",  "Spectator 1", "units", 'F' );
   //factory->AddSpectator( "spec2 := var1*3",  "Spectator 2", "units", 'F' );

   TTree *signal     = (TTree*) chsignal;
   TTree *background = (TTree*) chbackground;
   
   std::cout << "--- TMVAClassification       : Using bkg input files: -------------------" <<  std::endl;

   TObjArray *listOfBkgFiles = chbackground->GetListOfFiles();
   TIter bkgFileIter(listOfBkgFiles);
   TChainElement* currentBkgFile = 0;

   while((currentBkgFile = (TChainElement*)bkgFileIter.Next())) {
     std::cout << currentBkgFile->GetTitle() << std::endl;
   }

   std::cout << "--- TMVAClassification       : Using sig input files: -------------------" <<  std::endl;
   
   TObjArray *listOfSigFiles = chsignal->GetListOfFiles();
   TIter sigFileIter(listOfSigFiles);
   TChainElement* currentSigFile = 0;

   while((currentSigFile = (TChainElement*)sigFileIter.Next())) {
     std::cout << currentSigFile->GetTitle() << std::endl;
   }

   // global event weights per tree (see below for setting event-wise weights)
   Double_t signalWeight     = 1.0;
   Double_t backgroundWeight = 1.0;
   
   // You can add an arbitrary number of signal or background trees
   factory->AddSignalTree    ( signal,     signalWeight     );
   factory->AddBackgroundTree( background, backgroundWeight );
      
   // To give different trees for training and testing, do as follows:
   //    factory->AddSignalTree( signalTrainingTree, signalTrainWeight, "Training" );
   //    factory->AddSignalTree( signalTestTree,     signalTestWeight,  "Test" );
   
   // Use the following code instead of the above two or four lines to add signal and background
   // training and test events "by hand"
   // NOTE that in this case one should not give expressions (such as "var1+var2") in the input
   //      variable definition, but simply compute the expression before adding the event
   //
   //     // --- begin ----------------------------------------------------------
   //     std::vector<Double_t> vars( 4 ); // vector has size of number of input variables
   //     Float_t  treevars[4], weight;
   //     
   //     // Signal
   //     for (UInt_t ivar=0; ivar<4; ivar++) signal->SetBranchAddress( Form( "var%i", ivar+1 ), &(treevars[ivar]) );
   //     for (UInt_t i=0; i<signal->GetEntries(); i++) {
   //        signal->GetEntry(i);
   //        for (UInt_t ivar=0; ivar<4; ivar++) vars[ivar] = treevars[ivar];
   //        // add training and test events; here: first half is training, second is testing
   //        // note that the weight can also be event-wise
   //        if (i < signal->GetEntries()/2.0) factory->AddSignalTrainingEvent( vars, signalWeight );
   //        else                              factory->AddSignalTestEvent    ( vars, signalWeight );
   //     }
   //   
   //     // Background (has event weights)
   //     background->SetBranchAddress( "weight", &weight );
   //     for (UInt_t ivar=0; ivar<4; ivar++) background->SetBranchAddress( Form( "var%i", ivar+1 ), &(treevars[ivar]) );
   //     for (UInt_t i=0; i<background->GetEntries(); i++) {
   //        background->GetEntry(i);
   //        for (UInt_t ivar=0; ivar<4; ivar++) vars[ivar] = treevars[ivar];
   //        // add training and test events; here: first half is training, second is testing
   //        // note that the weight can also be event-wise
   //        if (i < background->GetEntries()/2) factory->AddBackgroundTrainingEvent( vars, backgroundWeight*weight );
   //        else                                factory->AddBackgroundTestEvent    ( vars, backgroundWeight*weight );
   //     }
         // --- end ------------------------------------------------------------
   //
   // --- end of tree registration 
   
   // Set individual event weights (the variables must exist in the original TTree)
   factory->SetSignalWeightExpression    ("event_scale1fb * lepsoft_fr");
   factory->SetBackgroundWeightExpression("event_scale1fb * lepsoft_fr");
   
   if( doMultipleOutputs ){
     multifactory->AddTree(signal,"Signal");
     multifactory->SetSignalWeightExpression    ("event_scale1fb");
     multifactory->SetBackgroundWeightExpression("event_scale1fb");
     multifactory->SetWeightExpression("event_scale1fb");
     
     if( includeBkg["ww"] ){
       TTree* ww = (TTree*) chww;
       multifactory->AddTree(ww,"WW");
       cout << "Added WW to multi-MVA" << endl;
     }
     if( includeBkg["wjets"] ){
       TTree* wjets = (TTree*) chwjets;
       multifactory->AddTree(wjets,"WJets");
       cout << "Added W+jets to multi-MVA" << endl;
     }
     if( includeBkg["tt"] ){
       TTree* tt = (TTree*) chtt;
       multifactory->AddTree(tt,"tt");
       cout << "Added ttbar multi-MVA" << endl;
     }
   }

   // Apply additional cuts on the signal and background samples (can be different)
   TCut mycuts = sel; // for example: TCut mycuts = "abs(var1)<0.5 && abs(var2-0.5)<1";
   TCut mycutb = sel; // for example: TCut mycutb = "abs(var1)<0.5";

   // Tell the factory how to use the training and testing events
   //
   // If no numbers of events are given, half of the events in the tree are used 
   // for training, and the other half for testing:
   //    factory->PrepareTrainingAndTestTree( mycut, "SplitMode=random:!V" );
   // To also specify the number of testing events, use:
   //    factory->PrepareTrainingAndTestTree( mycut,
   //                                         "NSigTrain=3000:NBkgTrain=3000:NSigTest=3000:NBkgTest=3000:SplitMode=Random:!V" );
   
   //Use random splitting
   factory->PrepareTrainingAndTestTree( mycuts, mycutb,
                                        "nTrain_Signal=0:nTrain_Background=0:SplitMode=Random:NormMode=NumEvents:!V" );

   if( doMultipleOutputs ){
     multifactory->PrepareTrainingAndTestTree( mycuts, mycutb,
                                               "nTrain_Signal=0:nTrain_Background=0:SplitMode=Random:NormMode=NumEvents:!V" );
   }

   //Use alternate splitting 
   //(this is preferable since its easier to track which events were used for training, but the job crashes! need to fix this...)
   //factory->PrepareTrainingAndTestTree( mycuts, mycutb,
   //                                     "nTrain_Signal=0:nTrain_Background=0:SplitMode=Alternate:NormMode=NumEvents:!V" );

   // ---- Book MVA methods
   //
   // Please lookup the various method configuration options in the corresponding cxx files, eg:
   // src/MethoCuts.cxx, etc, or here: http://tmva.sourceforge.net/optionRef.html
   // it is possible to preset ranges in the option string in which the cut optimisation should be done:
   // "...:CutRangeMin[2]=-1:CutRangeMax[2]=1"...", where [2] is the third input variable

   // Cut optimisation
   if (Use["Cuts"])
      factory->BookMethod( TMVA::Types::kCuts, "Cuts",
                           "!H:!V:FitMethod=MC:EffSel:SampleSize=200000:VarProp=FSmart" );

   if (Use["CutsD"])
      factory->BookMethod( TMVA::Types::kCuts, "CutsD",
                           "!H:!V:FitMethod=MC:EffSel:SampleSize=200000:VarProp=FSmart:VarTransform=Decorrelate" );

   if (Use["CutsPCA"])
      factory->BookMethod( TMVA::Types::kCuts, "CutsPCA",
                           "!H:!V:FitMethod=MC:EffSel:SampleSize=200000:VarProp=FSmart:VarTransform=PCA" );

   if (Use["CutsGA"])
      factory->BookMethod( TMVA::Types::kCuts, "CutsGA",
                           "H:!V:FitMethod=GA:CutRangeMin[0]=-10:CutRangeMax[0]=10:VarProp[1]=FMax:EffSel:Steps=30:Cycles=3:PopSize=400:SC_steps=10:SC_rate=5:SC_factor=0.95" );

   if (Use["CutsSA"])
      factory->BookMethod( TMVA::Types::kCuts, "CutsSA",
                           "!H:!V:FitMethod=SA:EffSel:MaxCalls=150000:KernelTemp=IncAdaptive:InitialTemp=1e+6:MinTemp=1e-6:Eps=1e-10:UseDefaultScale" );

   // Likelihood ("naive Bayes estimator")
   if (Use["Likelihood"])
      factory->BookMethod( TMVA::Types::kLikelihood, "Likelihood",
                           "H:!V:!TransformOutput:PDFInterpol=Spline2:NSmoothSig[0]=20:NSmoothBkg[0]=20:NSmoothBkg[1]=10:NSmooth=1:NAvEvtPerBin=50" );

   // Decorrelated likelihood
   if (Use["LikelihoodD"])
      factory->BookMethod( TMVA::Types::kLikelihood, "LikelihoodD",
                           "!H:!V:TransformOutput:PDFInterpol=Spline2:NSmoothSig[0]=20:NSmoothBkg[0]=20:NSmooth=5:NAvEvtPerBin=50:VarTransform=Decorrelate" );

   // PCA-transformed likelihood
   if (Use["LikelihoodPCA"])
      factory->BookMethod( TMVA::Types::kLikelihood, "LikelihoodPCA",
                           "!H:!V:!TransformOutput:PDFInterpol=Spline2:NSmoothSig[0]=20:NSmoothBkg[0]=20:NSmooth=5:NAvEvtPerBin=50:VarTransform=PCA" ); 

   // Use a kernel density estimator to approximate the PDFs
   if (Use["LikelihoodKDE"])
      factory->BookMethod( TMVA::Types::kLikelihood, "LikelihoodKDE",
                           "!H:!V:!TransformOutput:PDFInterpol=KDE:KDEtype=Gauss:KDEiter=Adaptive:KDEFineFactor=0.3:KDEborder=None:NAvEvtPerBin=50" ); 

   // Use a variable-dependent mix of splines and kernel density estimator
   if (Use["LikelihoodMIX"])
      factory->BookMethod( TMVA::Types::kLikelihood, "LikelihoodMIX",
                           "!H:!V:!TransformOutput:PDFInterpolSig[0]=KDE:PDFInterpolBkg[0]=KDE:PDFInterpolSig[1]=KDE:PDFInterpolBkg[1]=KDE:PDFInterpolSig[2]=Spline2:PDFInterpolBkg[2]=Spline2:PDFInterpolSig[3]=Spline2:PDFInterpolBkg[3]=Spline2:KDEtype=Gauss:KDEiter=Nonadaptive:KDEborder=None:NAvEvtPerBin=50" ); 

   // Test the multi-dimensional probability density estimator
   // here are the options strings for the MinMax and RMS methods, respectively:
   //      "!H:!V:VolumeRangeMode=MinMax:DeltaFrac=0.2:KernelEstimator=Gauss:GaussSigma=0.3" );
   //      "!H:!V:VolumeRangeMode=RMS:DeltaFrac=3:KernelEstimator=Gauss:GaussSigma=0.3" );
   if (Use["PDERS"])
      factory->BookMethod( TMVA::Types::kPDERS, "PDERS",
                           "!H:!V:NormTree=T:VolumeRangeMode=Adaptive:KernelEstimator=Gauss:GaussSigma=0.3:NEventsMin=400:NEventsMax=600" );

   if (Use["PDERSD"])
      factory->BookMethod( TMVA::Types::kPDERS, "PDERSD",
                           "!H:!V:VolumeRangeMode=Adaptive:KernelEstimator=Gauss:GaussSigma=0.3:NEventsMin=400:NEventsMax=600:VarTransform=Decorrelate" );

   if (Use["PDERSPCA"])
      factory->BookMethod( TMVA::Types::kPDERS, "PDERSPCA",
                           "!H:!V:VolumeRangeMode=Adaptive:KernelEstimator=Gauss:GaussSigma=0.3:NEventsMin=400:NEventsMax=600:VarTransform=PCA" );

   // Multi-dimensional likelihood estimator using self-adapting phase-space binning
   if (Use["PDEFoam"])
      factory->BookMethod( TMVA::Types::kPDEFoam, "PDEFoam",
                           "H:!V:SigBgSeparate=F:TailCut=0.001:VolFrac=0.0333:nActiveCells=500:nSampl=2000:nBin=5:Nmin=100:Kernel=None:Compress=T" );

   if (Use["PDEFoamBoost"])
      factory->BookMethod( TMVA::Types::kPDEFoam, "PDEFoamBoost",
                           "!H:!V:Boost_Num=30:Boost_Transform=linear:SigBgSeparate=F:MaxDepth=4:UseYesNoCell=T:DTLogic=MisClassificationError:FillFoamWithOrigWeights=F:TailCut=0:nActiveCells=500:nBin=20:Nmin=400:Kernel=None:Compress=T" );

   // K-Nearest Neighbour classifier (KNN)
   if (Use["KNN"])
      factory->BookMethod( TMVA::Types::kKNN, "KNN",
                           "H:nkNN=20:ScaleFrac=0.8:SigmaFact=1.0:Kernel=Gaus:UseKernel=F:UseWeight=T:!Trim" );

   // H-Matrix (chi2-squared) method
   if (Use["HMatrix"])
      factory->BookMethod( TMVA::Types::kHMatrix, "HMatrix", "!H:!V" );

   // Linear discriminant (same as Fisher discriminant)
   if (Use["LD"])
      factory->BookMethod( TMVA::Types::kLD, "LD", "H:!V:VarTransform=None:CreateMVAPdfs:PDFInterpolMVAPdf=Spline2:NbinsMVAPdf=50:NsmoothMVAPdf=10" );

   // Fisher discriminant (same as LD)
   if (Use["Fisher"])
      factory->BookMethod( TMVA::Types::kFisher, "Fisher", "H:!V:Fisher:CreateMVAPdfs:PDFInterpolMVAPdf=Spline2:NbinsMVAPdf=50:NsmoothMVAPdf=10" );

   // Fisher with Gauss-transformed input variables
   if (Use["FisherG"])
      factory->BookMethod( TMVA::Types::kFisher, "FisherG", "H:!V:VarTransform=Gauss" );

   // Composite classifier: ensemble (tree) of boosted Fisher classifiers
   if (Use["BoostedFisher"])
      factory->BookMethod( TMVA::Types::kFisher, "BoostedFisher", 
                           "H:!V:Boost_Num=20:Boost_Transform=log:Boost_Type=AdaBoost:Boost_AdaBoostBeta=0.2" );

   // Function discrimination analysis (FDA) -- test of various fitters - the recommended one is Minuit (or GA or SA)
   if (Use["FDA_MC"])
      factory->BookMethod( TMVA::Types::kFDA, "FDA_MC",
                           "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=MC:SampleSize=100000:Sigma=0.1" );

   if (Use["FDA_GA"]) // can also use Simulated Annealing (SA) algorithm (see Cuts_SA options])
      factory->BookMethod( TMVA::Types::kFDA, "FDA_GA",
                           "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=GA:PopSize=300:Cycles=3:Steps=20:Trim=True:SaveBestGen=1" );

   if (Use["FDA_SA"]) // can also use Simulated Annealing (SA) algorithm (see Cuts_SA options])
      factory->BookMethod( TMVA::Types::kFDA, "FDA_SA",
                           "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=SA:MaxCalls=15000:KernelTemp=IncAdaptive:InitialTemp=1e+6:MinTemp=1e-6:Eps=1e-10:UseDefaultScale" );

   if (Use["FDA_MT"])
      factory->BookMethod( TMVA::Types::kFDA, "FDA_MT",
                           "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=MINUIT:ErrorLevel=1:PrintLevel=-1:FitStrategy=2:UseImprove:UseMinos:SetBatch" );

   if (Use["FDA_GAMT"])
      factory->BookMethod( TMVA::Types::kFDA, "FDA_GAMT",
                           "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=GA:Converger=MINUIT:ErrorLevel=1:PrintLevel=-1:FitStrategy=0:!UseImprove:!UseMinos:SetBatch:Cycles=1:PopSize=5:Steps=5:Trim" );

   if (Use["FDA_MCMT"])
      factory->BookMethod( TMVA::Types::kFDA, "FDA_MCMT",
                           "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=MC:Converger=MINUIT:ErrorLevel=1:PrintLevel=-1:FitStrategy=0:!UseImprove:!UseMinos:SetBatch:SampleSize=20" );

   // TMVA ANN: MLP (recommended ANN) -- all ANNs in TMVA are Multilayer Perceptrons
   if (Use["MLP"])
      factory->BookMethod( TMVA::Types::kMLP, "MLP", "H:!V:NeuronType=tanh:VarTransform=N:NCycles=600:HiddenLayers=N+5:TestRate=5:!UseRegulator" );

   if (Use["MLPBFGS"])
      factory->BookMethod( TMVA::Types::kMLP, "MLPBFGS", "H:!V:NeuronType=tanh:VarTransform=N:NCycles=600:HiddenLayers=N+5:TestRate=5:TrainingMethod=BFGS:!UseRegulator" );

   if (Use["MLPBNN"])
      factory->BookMethod( TMVA::Types::kMLP, "MLPBNN", "H:!V:NeuronType=tanh:VarTransform=N:NCycles=600:HiddenLayers=N+5:TestRate=5:TrainingMethod=BFGS:UseRegulator" ); // BFGS training with bayesian regulators

   // CF(Clermont-Ferrand)ANN
   if (Use["CFMlpANN"])
      factory->BookMethod( TMVA::Types::kCFMlpANN, "CFMlpANN", "!H:!V:NCycles=2000:HiddenLayers=N+1,N"  ); // n_cycles:#nodes:#nodes:...  

   // Tmlp(Root)ANN
   if (Use["TMlpANN"])
      factory->BookMethod( TMVA::Types::kTMlpANN, "TMlpANN", "!H:!V:NCycles=200:HiddenLayers=N+1,N:LearningMethod=BFGS:ValidationFraction=0.3"  ); // n_cycles:#nodes:#nodes:...

   // Support Vector Machine
   if (Use["SVM"])
      factory->BookMethod( TMVA::Types::kSVM, "SVM", "Gamma=0.25:Tol=0.001:VarTransform=Norm" );

   // Boosted Decision Trees
   if (Use["BDTG"]) // Gradient Boost
      factory->BookMethod( TMVA::Types::kBDT, "BDTG",
                           "!H:!V:NTrees=1000:BoostType=Grad:Shrinkage=0.10:UseBaggedGrad:GradBaggingFraction=0.5:nCuts=20:NNodesMax=5" );

   if (Use["BDT"])  // Adaptive Boost
      factory->BookMethod( TMVA::Types::kBDT, "BDT",
                           "!H:!V:NTrees=850:nEventsMin=150:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:SeparationType=GiniIndex:nCuts=20:PruneMethod=NoPruning" );

   if (Use["BDTB"]) // Bagging
      factory->BookMethod( TMVA::Types::kBDT, "BDTB",
                           "!H:!V:NTrees=400:BoostType=Bagging:SeparationType=GiniIndex:nCuts=20:PruneMethod=NoPruning" );

   if (Use["BDTD"]) // Decorrelation + Adaptive Boost
      factory->BookMethod( TMVA::Types::kBDT, "BDTD",
                           "!H:!V:NTrees=400:nEventsMin=400:MaxDepth=3:BoostType=AdaBoost:SeparationType=GiniIndex:nCuts=20:PruneMethod=NoPruning:VarTransform=Decorrelate" );

   // RuleFit -- TMVA implementation of Friedman's method
   if (Use["RuleFit"])
      factory->BookMethod( TMVA::Types::kRuleFit, "RuleFit",
                           "H:!V:RuleFitModule=RFTMVA:Model=ModRuleLinear:MinImp=0.001:RuleMinDist=0.001:NTrees=20:fEventsMin=0.01:fEventsMax=0.5:GDTau=-1.0:GDTauPrec=0.01:GDStep=0.01:GDNSteps=10000:GDErrScale=1.02" );

   if( doMultipleOutputs ){
     if (Use["multi_BDTG"]) // gradient boosted decision trees
       multifactory->BookMethod( TMVA::Types::kBDT, "BDTG", "!H:!V:NTrees=1000:BoostType=Grad:Shrinkage=0.10:UseBaggedGrad:GradBaggingFraction=0.50:nCuts=20:NNodesMax=8");
     if (Use["multi_MLP"]) // neural network
       multifactory->BookMethod( TMVA::Types::kMLP, "MLP", "!H:!V:NeuronType=tanh:NCycles=1000:HiddenLayers=N+5,5:TestRate=5:EstimatorType=MSE");
     if (Use["multi_FDA_GA"]) // functional discriminant with GA minimizer
       multifactory->BookMethod( TMVA::Types::kFDA, "FDA_GA", "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=GA:PopSize=300:Cycles=3:Steps=20:Trim=True:SaveBestGen=1" );
   }
   
   // For an example of the category classifier usage, see: TMVAClassificationCategory

   // --------------------------------------------------------------------------------------------------

   // ---- Now you can optimize the setting (configuration) of the MVAs using the set of training events

   // factory->OptimizeAllMethods("SigEffAt001","Scan");
   // factory->OptimizeAllMethods("ROCIntegral","GA");

   // --------------------------------------------------------------------------------------------------

   // ---- Now you can tell the factory to train, test, and evaluate the MVAs
  
   // Train MVAs using the set of training events
   factory->TrainAllMethods();
  
   // ---- Evaluate all MVAs using the set of test events
   factory->TestAllMethods();
  
   // ----- Evaluate and compare performance of all configured MVAs
   factory->EvaluateAllMethods();
  
   if( doMultipleOutputs ){
     // Train nulti-MVAs using the set of training events
     multifactory->TrainAllMethods();
     
     // ---- Evaluate all multi-MVAs using the set of test events
     multifactory->TestAllMethods();
     
     // ----- Evaluate and compare performance of all configured multi-MVAs
     multifactory->EvaluateAllMethods();
   }
   
   // --------------------------------------------------------------

   // Save the output
   outputFile->Close();
   if( doMultipleOutputs )  multioutputFile->Close();

   std::cout << "==> Wrote root file: " << outputFile->GetName() << std::endl;
   std::cout << "==> TMVAClassification is done!" << std::endl;
  
   delete factory;

   // Launch the GUI for the root macros
   if (!gROOT->IsBatch()) TMVAGui( outfileName );
}
