{

  float xsec250 = 5.57596;
  float xsec450 = 0.169668;
  float xsec600 = 0.0248009;
  float xsec650 = 0.0139566;
  float lumi    = 19500.0;

  // T2tt cut-based
  // TFile* f = TFile::Open("rootfiles/T2tt_histos.root");
  // int n = 8;

  // T2tt cut-based
  // TFile* f = TFile::Open("rootfiles/T2tt_BDT_histos.root");
  // int n = 6;

  // T2bw cut-based
  //TFile* f = TFile::Open("rootfiles/T2bw_MG_x50_histos.root");
  // TFile* f = TFile::Open("rootfiles/T2bw_MG_x50T2BW_SS_histos.root");
  // int n = 8;

  // T2bw x=0.5 BDT
  TFile* f = TFile::Open("rootfiles/T2bw_MG_x50T2BW_SS_BDT_histos.root");
  int n = 4;

  // // T2bw x=0.5 BDT
  // TFile* f = TFile::Open("rootfiles/T2bw_MG_x75T2BW_SS_BDT_histos.root");
  // int n = 4;

  // T2bw x=0.5 BDT
  // TFile* f = TFile::Open("rootfiles/T2bw_MG_x25T2BW_SS_BDT_histos.root");
  // int n = 2;



  TH2F* h[n];


  for( int i = 0 ; i < n ; i++ ){
    h[i] = (TH2F*) f->Get(Form("heff_%i",i));
    int bin250 = h[i]->FindBin(250,50);

    float eff250 = h[i]->GetBinContent(bin250)/100.0;
    float err250 = h[i]->GetBinError(bin250)/100.0;

    cout << Form("$%.1f",eff250 * lumi * xsec250) << " \\pm " << Form("%.1f$",err250 * lumi * xsec250) << "   &   ";
  }

  cout << endl;


  for( int i = 0 ; i < n ; i++ ){
    int bin650 = h[i]->FindBin(650,50);

    float eff650 = h[i]->GetBinContent(bin650)/100.0;
    float err650 = h[i]->GetBinError(bin650)/100.0;

    cout << Form("$%.1f",eff650 * lumi * xsec650) << " \\pm " << Form("%.1f$",err650 * lumi * xsec650) << "   &   ";
  }

  cout << endl;




  /*
  for( int i = 0 ; i < n ; i++ ){
    h[i] = (TH2F*) f->Get(Form("heff_%i",i));
    int bin450 = h[i]->FindBin(450,50);

    float eff450 = h[i]->GetBinContent(bin450)/100.0;
    float err450 = h[i]->GetBinError(bin450)/100.0;

    cout << Form("$%.1f",eff450 * lumi * xsec450) << " \\pm " << Form("%.1f$",err450 * lumi * xsec450) << "   &   ";
  }

  cout << endl;


  for( int i = 0 ; i < n ; i++ ){
    int bin600 = h[i]->FindBin(600,100);

    float eff600 = h[i]->GetBinContent(bin600)/100.0;
    float err600 = h[i]->GetBinError(bin600)/100.0;

    cout << Form("$%.1f",eff600 * lumi * xsec600) << " \\pm " << Form("%.1f$",err600 * lumi * xsec600) << "   &   ";
  }

  cout << endl;
  */


















}
