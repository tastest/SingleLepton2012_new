#include <iomanip>

Double_t count_chain(TChain* ch, TCut cut, bool signal = false){
    TH1D *h_weights = new TH1D("h_weights", "h_weights", 1, 0., 2.);

    if ( signal )
        ch->Draw("1>>h_weights", "mini_weight"*cut );    
    else
        ch->Draw("1>>h_weights", "mini_weight"*cut );

    double x =  h_weights->GetSumOfWeights();
    //    double x =  h_weights->Integral();
    delete h_weights;
    return x ;
}

printCount(const char* name, TChain& ch, TCut extra = "", bool isSignal = false){
    TCut rho("rhovor>0 && rhovor<40");
    TCut filters("isdata==0 || (csc==0 && hbhe==1 && hcallaser==1 && ecaltp==1 && trkfail==1 && eebadsc==1 && hbhenew==1)");
    TCut goodlep("ngoodlep > 0 && abs( pflep1.Pt() - lep1.Pt() ) < 10.0 && abs(isopf1 * lep1.Pt() ) < 5.0");
    TCut el("leptype==0 && abs(lep1->Eta())<1.4442 && lep1->Pt()>30.0 && eoverpin < 4.0 && (isdata==0 || ele27wp80==1)");
    TCut mu("leptype==1 && abs(lep1->Eta())<2.1 && lep1->Pt()>25.0 && (isdata==0 || isomu24==1)");
    TCut njets4("mini_njets >= 4");
    TCut btag1("mini_nb >= 1");
    TCut passisotrk("mini_passisotrk==1");
    TCut passtauveto("mini_passtauveto==1");
    TCut met100("t1metphicorr > 100");
    TCut mt120("t1metphicorrmt > 120");

    TCut sel;
    sel += rho;
    sel += filters;
    sel += goodlep;
    sel += (el||mu);
    sel += njets4;
    sel += btag1;
    sel += passisotrk;
    sel += met100;
    sel += mt120;

    if (!isSignal)
    sel += passtauveto;

    //-------------------------------------------
    // THESE CUTS DEFINE PRESELECTION REGION
    //-------------------------------------------

    TCut sel0 = sel+extra;
    TCut sel1 = sel0+btag1;
    TCut sel2 = sel1+passisotrk;
    TCut sel3 = sel2+met100;
    TCut sel4 = sel3+mt120;
    TCut sel5 = sel4+passtauveto;

    //double scale = isSignal? 19.5*0.73:1.;
    double scale = 1 ;

    cout << "|  " << name << "  | " << fixed  << setprecision(3)
         <<  count_chain(&ch, sel0,isSignal)*scale << " | " 
//         <<  count_chain(&ch, sel1,isSignal)*scale << " | " 
//         <<  count_chain(&ch, sel2,isSignal)*scale << " | " 
//         <<  count_chain(&ch, sel3,isSignal)*scale << " | " 
//         <<  count_chain(&ch, sel4,isSignal)*scale << " | " 
//         <<  count_chain(&ch, sel5,isSignal)*scale << " |" 
         <<  endl;
}


tableBen(){
//    TString base_path = "/nfs-3/userdata/stop/MiniBabies/V00-02-18__V00-03-00_3jetsMET50_bkg/";
    TString base_path = "/nfs-3/userdata/stop/MiniBabies/V00-02-20__V00-03-01__BDT004__4jetsMET100MT120fix_bkg/";

    TChain ch_ttsl("t");
    ch_ttsl.Add(base_path+"ttsl_powheg.root");
    ch_ttsl.Add(base_path+"tW_lepsl.root");

    TChain ch_ttdl("t");
    ch_ttdl.Add(base_path+"ttdl_powheg.root");

    TChain ch_rare("t");
    ch_rare.Add(base_path+"DY1to4Jtot*root");
    ch_rare.Add(base_path+"ttV.root");
    ch_rare.Add(base_path+"tW_lepdl.root");
    ch_rare.Add(base_path+"diboson*root");
    ch_rare.Add(base_path+"triboson.root");

    TChain ch_wjets("t");
    ch_wjets.Add(base_path+"w1to4jets*root");

    TChain ch_signal("t");
//    ch_signal.Add("/nfs-3/userdata/stop/MiniBabies/V00-02-18__V00-03-00_3jetsMET50_T2tt/merged*.root");
    ch_signal.Add("/nfs-3/userdata/stop/MiniBabies/V00-02-18__V00-03-01__BDT004__4jetsMET100MT120fix_T2tt/T2tt*.root");

    printCount("wjets", ch_wjets);
    printCount("rare", ch_rare);
    printCount("ttsl", ch_ttsl);
    printCount("ttdl", ch_ttdl);
    printCount("T2tt 650,50", ch_signal,  "mg == 650 && ml == 50", true);
    printCount("T2tt 400,200", ch_signal, "mg == 400 && ml == 200", true);
    printCount("T2tt 250,25", ch_signal,  "mg == 250 && ml == 25", true);
//    printCount("Others", ch_kids);
}
