{
   // --------- S t y l e ---------------------------
   const Bool_t UsePaperStyle = 0;
   // -----------------------------------------------
   
   TString curDynamicPath( gSystem->GetDynamicPath() );
   gSystem->SetDynamicPath( "/tas/home/benhoob/tmva-V04-01-00/lib:" + curDynamicPath );


   TString curIncludePath(gSystem->GetIncludePath());
   gSystem->SetIncludePath( " -I/tas/home/benhoob/tmva-V04-01-00/include " + curIncludePath );

   // load TMVA shared library created in local release 
   // (not required anymore with the use of rootmaps, but problems with MAC OSX)
   if (TString(gSystem->GetBuildArch()).Contains("macosx") ) gSystem->Load( "libTMVA.1" );

   TMVA::Tools::Instance();

   // welcome the user
   TMVA::gTools().TMVAWelcomeMessage();
   
#include "tmvaglob.C"
   
   TMVAGlob::SetTMVAStyle();
   cout << "TMVAlogon: use \"" << gStyle->GetName() << "\" style [" << gStyle->GetTitle() << "]" << endl;
   cout << endl;

   gROOT->ProcessLine(".L /tas03/home/benhoob/.root/tdrstyle_SUSY.C");
   setTDRStyle();

}
