// old, not sure where this came from...
//const float sf = 8.6e-4;
// from inc_2j region in wjets_nobb in V24_crs
const float sf_2med = 9.2e-4;
const float sf_1med = 3.0e-2;
const float sf_1med1loose = 4.9e-3;
const float sf_2loose = 2.6e-2;

void makeBgDir(TString input_region, TString output_region, TFile* fin, TFile* fout, float sf = sf_2med) {

  std::cout << "scaling hists from region " << input_region << " to " << output_region << std::endl;

  TDirectoryFile* dirf_out = fout->mkdir(output_region);

  TDirectoryFile* dirf_in = (TDirectoryFile*) fin->Get(input_region);
  TIter nextkey(dirf_in->GetListOfKeys());
  TKey *key;
  while (key = (TKey*)nextkey()) {
    TObject* obj = key->ReadObj();
    if (TString(obj->ClassName()) != "TH1F") continue;
    TH1F* h_orig = (TH1F*)obj;
    dirf_out->cd();
    TH1F* h_scaled = h_orig->Clone(h_orig->GetName());
    h_scaled->Scale(sf);
    //    h_scaled->Write();
  }

  dirf_out->Write();

}

// // run after making all other hists
// void makeMbbHists(TFile* fout) {

//   TIter nextkey(fout->GetListOfKeys());
//   TKey *key;
//   while (key = (TKey*)nextkey()) {
//     TObject* obj = key->ReadObj();
//     std::cout << "------ read object" << std::endl;
//     if (TString(obj->ClassName()) != "TDirectoryFile") continue;
//     TDirectoryFile* dirf = (TDirectoryFile*) obj;
//     TString name(dirf->GetName());
//     if (!name.Contains("cr14")) continue;
//     TString compname(name);
//     compname.ReplaceAll("cr14","sig_metlast");
//     TDirectoryFile* dirf_mbb = (TDirectoryFile*) fout->Get(compname);
//     std::cout << "------ got 2nd dir" << std::endl;
//     if (!dirf_mbb) continue;
//     std::cout << "making full bbmass hist for region: " << name << std::endl;
//     TString newname(name);
//     newname.ReplaceAll("cr14","sig_bbmasslast");
//     TDirectoryFile* dirf_out = fout->mkdir(newname);
//     TH1F* h_inv = (TH1F*) dirf->Get("h_bbmass");
//     if (!h_inv) {
//       std::cout << "problem with getting hist in dir: " << name << std::endl;
//       continue;
//     }
//     TH1F* h_mbb = (TH1F*) dirf_mbb->Get("h_bbmass");
//     if (!h_mbb) {
//       std::cout << "problem with getting hist in dir: " << compname << std::endl;
//       continue;
//     }
//     dirf_out->cd();
//     TH1F* h_full = h_inv->Clone(h_inv->GetName());
//     h_full->Add(h_mbb);
//     //    h_full->Write();
//     dirf_out->Write();
//   }

// }

void makeWLightBg(TString dir = "") {

  TFile* fin = new TFile(dir+"/data_histos.root");
  TFile* fout = new TFile(dir+"/wlight_dd_histos.root","RECREATE");

  TIter nextkey(fin->GetListOfKeys());
  TKey *key;
  while (key = (TKey*)nextkey()) {
    TObject* obj = key->ReadObj();
    if (TString(obj->ClassName()) != "TDirectoryFile") continue;
    TDirectoryFile* dirf = (TDirectoryFile*) obj;
    TString name(dirf->GetName());
    TString newname(name);
    TString newname2(name);
    TString newname3(name);
    TString newname4(name);
    if (name.Contains("cr5_metlast")) {
      makeBgDir(name,newname.ReplaceAll("cr5","sig"),fin,fout);
      // makeBgDir(name,newname2.ReplaceAll("cr5_metlast","sig_1med"),fin,fout,sf_1med);
      // makeBgDir(name,newname3.ReplaceAll("cr5_metlast","sig_2loose"),fin,fout,sf_2loose);
      // makeBgDir(name,newname4.ReplaceAll("cr5_metlast","sig_1med1loose"),fin,fout,sf_1med1loose);
    } else if (name.Contains("cr5_invmass")) {
      makeBgDir(name,newname.ReplaceAll("cr5_invmass","cr14"),fin,fout);
    } else if (name.Contains("cr5_bbmasslast")) {
      makeBgDir(name,newname.ReplaceAll("cr5_bbmasslast","sig_bbmasslast"),fin,fout);
    } else if (name.Contains("cr5_highmass")) {
      makeBgDir(name,newname.ReplaceAll("cr5_highmass","cr1_metlast"),fin,fout);
    } else if (name.Contains("cr5_lowmass")) {
      makeBgDir(name,newname.ReplaceAll("cr5_lowmass","cr8_metlast"),fin,fout);
    } 
  }

  //  makeMbbHists(fout);

  fin->Close();
  //  fout->Write();
  fout->Close();

}
