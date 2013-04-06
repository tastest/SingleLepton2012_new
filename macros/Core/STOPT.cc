#include "STOPT.h"
STOPT stopt;
namespace Stop {
	int &acc_2010() { return stopt.acc_2010(); }
	int &acc_highmet() { return stopt.acc_highmet(); }
	int &acc_highht() { return stopt.acc_highht(); }
	int &eldup() { return stopt.eldup(); }
	int &csc() { return stopt.csc(); }
	int &hbhe() { return stopt.hbhe(); }
	int &hbhenew() { return stopt.hbhenew(); }
	int &hcallaser() { return stopt.hcallaser(); }
	int &ecaltp() { return stopt.ecaltp(); }
	int &trkfail() { return stopt.trkfail(); }
	int &eebadsc() { return stopt.eebadsc(); }
	int &lep1_badecallaser() { return stopt.lep1_badecallaser(); }
	int &lep2_badecallaser() { return stopt.lep2_badecallaser(); }
	int &isdata() { return stopt.isdata(); }
	int &jetid() { return stopt.jetid(); }
	int &jetid30() { return stopt.jetid30(); }
	int &json() { return stopt.json(); }
	float &htoffset() { return stopt.htoffset(); }
	float &htuncor() { return stopt.htuncor(); }
	float &ptt() { return stopt.ptt(); }
	float &pttbar() { return stopt.pttbar(); }
	float &ptttbar() { return stopt.ptttbar(); }
	float &mttbar() { return stopt.mttbar(); }
	int &npartons() { return stopt.npartons(); }
	int &nwzpartons() { return stopt.nwzpartons(); }
	int &hyptype() { return stopt.hyptype(); }
	float &maxpartonpt() { return stopt.maxpartonpt(); }
	float &etattbar() { return stopt.etattbar(); }
	int &njetsoffset() { return stopt.njetsoffset(); }
	int &njetsuncor() { return stopt.njetsuncor(); }
	float &costhetaweight() { return stopt.costhetaweight(); }
	float &weight() { return stopt.weight(); }
	float &weightleft() { return stopt.weightleft(); }
	float &weightright() { return stopt.weightright(); }
	float &mutrigweight() { return stopt.mutrigweight(); }
	float &mutrigweight2() { return stopt.mutrigweight2(); }
	float &sltrigweight() { return stopt.sltrigweight(); }
	float &dltrigweight() { return stopt.dltrigweight(); }
	float &trgeff() { return stopt.trgeff(); }
	float &pthat() { return stopt.pthat(); }
	float &qscale() { return stopt.qscale(); }
	float &mgcor() { return stopt.mgcor(); }
	int &wflav() { return stopt.wflav(); }
	float &ksusy() { return stopt.ksusy(); }
	float &ksusyup() { return stopt.ksusyup(); }
	float &ksusydn() { return stopt.ksusydn(); }
	float &xsecsusy() { return stopt.xsecsusy(); }
	float &xsecsusy2() { return stopt.xsecsusy2(); }
	float &smeff() { return stopt.smeff(); }
	float &k() { return stopt.k(); }
	float &mllgen() { return stopt.mllgen(); }
	float &ptwgen() { return stopt.ptwgen(); }
	float &ptzgen() { return stopt.ptzgen(); }
	int &nlep() { return stopt.nlep(); }
	int &nosel() { return stopt.nosel(); }
	int &ngoodlep() { return stopt.ngoodlep(); }
	int &ngoodel() { return stopt.ngoodel(); }
	int &ngoodmu() { return stopt.ngoodmu(); }
	int &mull() { return stopt.mull(); }
	int &mult() { return stopt.mult(); }
	int &mullgen() { return stopt.mullgen(); }
	int &multgen() { return stopt.multgen(); }
	int &proc() { return stopt.proc(); }
	int &leptype() { return stopt.leptype(); }
	float &topmass() { return stopt.topmass(); }
	float &dilmass() { return stopt.dilmass(); }
	float &dilrecoil() { return stopt.dilrecoil(); }
	float &dilrecoilparl() { return stopt.dilrecoilparl(); }
	float &dilrecoilperp() { return stopt.dilrecoilperp(); }
	float &tcmet() { return stopt.tcmet(); }
	float &genmet() { return stopt.genmet(); }
	float &gensumet() { return stopt.gensumet(); }
	float &genmetphi() { return stopt.genmetphi(); }
	float &calomet() { return stopt.calomet(); }
	float &calometphi() { return stopt.calometphi(); }
	float &trkmet() { return stopt.trkmet(); }
	float &trkmetphi() { return stopt.trkmetphi(); }
	float &trkmet_nolepcorr() { return stopt.trkmet_nolepcorr(); }
	float &trkmetphi_nolepcorr() { return stopt.trkmetphi_nolepcorr(); }
	float &pfmet() { return stopt.pfmet(); }
	float &pfmetveto() { return stopt.pfmetveto(); }
	float &pfmetsig() { return stopt.pfmetsig(); }
	float &pfmetphi() { return stopt.pfmetphi(); }
	float &pfsumet() { return stopt.pfsumet(); }
	float &mucormet() { return stopt.mucormet(); }
	float &mucorjesmet() { return stopt.mucorjesmet(); }
	float &tcmet35X() { return stopt.tcmet35X(); }
	float &tcmetevent() { return stopt.tcmetevent(); }
	float &tcmetlooper() { return stopt.tcmetlooper(); }
	float &tcmetphi() { return stopt.tcmetphi(); }
	float &tcsumet() { return stopt.tcsumet(); }
	float &tcmetUp() { return stopt.tcmetUp(); }
	float &tcmetDown() { return stopt.tcmetDown(); }
	float &tcmetTest() { return stopt.tcmetTest(); }
	float &pfmetUp() { return stopt.pfmetUp(); }
	float &pfmetDown() { return stopt.pfmetDown(); }
	float &pfmetTest() { return stopt.pfmetTest(); }
	float &sumjetpt() { return stopt.sumjetpt(); }
	float &dileta() { return stopt.dileta(); }
	float &dilpt() { return stopt.dilpt(); }
	float &dildphi() { return stopt.dildphi(); }
	int &ngenjets() { return stopt.ngenjets(); }
	int &njpt() { return stopt.njpt(); }
	int &trgmu1() { return stopt.trgmu1(); }
	int &trgmu2() { return stopt.trgmu2(); }
	int &trgel1() { return stopt.trgel1(); }
	int &trgel2() { return stopt.trgel2(); }
	int &isomu24() { return stopt.isomu24(); }
	int &ele27wp80() { return stopt.ele27wp80(); }
	int &mm() { return stopt.mm(); }
	int &mmtk() { return stopt.mmtk(); }
	int &me() { return stopt.me(); }
	int &em() { return stopt.em(); }
	int &mu() { return stopt.mu(); }
	int &ee() { return stopt.ee(); }
	int &npfjets30() { return stopt.npfjets30(); }
	int &npfjets30lepcorr() { return stopt.npfjets30lepcorr(); }
	float &knjets() { return stopt.knjets(); }
	float &rhovor() { return stopt.rhovor(); }
	float &htpf30() { return stopt.htpf30(); }
	float &t1met10() { return stopt.t1met10(); }
	float &t1met20() { return stopt.t1met20(); }
	float &t1met30() { return stopt.t1met30(); }
	float &t1met10phi() { return stopt.t1met10phi(); }
	float &t1met20phi() { return stopt.t1met20phi(); }
	float &t1met30phi() { return stopt.t1met30phi(); }
	float &t1met10mt() { return stopt.t1met10mt(); }
	float &t1met20mt() { return stopt.t1met20mt(); }
	float &t1met30mt() { return stopt.t1met30mt(); }
	float &lepmetpt() { return stopt.lepmetpt(); }
	float &lept1met10pt() { return stopt.lept1met10pt(); }
	float &t1met10s() { return stopt.t1met10s(); }
	float &t1met10sphi() { return stopt.t1met10sphi(); }
	float &t1met10smt() { return stopt.t1met10smt(); }
	float &t1metphicorr() { return stopt.t1metphicorr(); }
	float &t1metphicorrup() { return stopt.t1metphicorrup(); }
	float &t1metphicorrdn() { return stopt.t1metphicorrdn(); }
	float &t1metphicorrphi() { return stopt.t1metphicorrphi(); }
	float &t1metphicorrphiup() { return stopt.t1metphicorrphiup(); }
	float &t1metphicorrphidn() { return stopt.t1metphicorrphidn(); }
	float &t1metphicorrlep() { return stopt.t1metphicorrlep(); }
	float &t1metphicorrlepphi() { return stopt.t1metphicorrlepphi(); }
	float &t1metphicorrmt() { return stopt.t1metphicorrmt(); }
	float &t1metphicorrmtup() { return stopt.t1metphicorrmtup(); }
	float &t1metphicorrmtdn() { return stopt.t1metphicorrmtdn(); }
	float &t1metphicorrlepmt() { return stopt.t1metphicorrlepmt(); }
	float &t1met_off() { return stopt.t1met_off(); }
	float &t1metphi_off() { return stopt.t1metphi_off(); }
	float &t1metmt_off() { return stopt.t1metmt_off(); }
	float &t1metphicorr_off() { return stopt.t1metphicorr_off(); }
	float &t1metphicorrphi_off() { return stopt.t1metphicorrphi_off(); }
	float &t1metphicorrmt_off() { return stopt.t1metphicorrmt_off(); }
	float &mt2bmin() { return stopt.mt2bmin(); }
	float &mt2blmin() { return stopt.mt2blmin(); }
	float &mt2wmin() { return stopt.mt2wmin(); }
	float &chi2min() { return stopt.chi2min(); }
	float &chi2minprob() { return stopt.chi2minprob(); }
	int &nbtagsssv() { return stopt.nbtagsssv(); }
	int &nbtagstcl() { return stopt.nbtagstcl(); }
	int &nbtagstcm() { return stopt.nbtagstcm(); }
	int &nbtagscsvl() { return stopt.nbtagscsvl(); }
	int &nbtagscsvm() { return stopt.nbtagscsvm(); }
	int &nbtagscsvt() { return stopt.nbtagscsvt(); }
	int &nbtagsssvcorr() { return stopt.nbtagsssvcorr(); }
	int &nbtagstclcorr() { return stopt.nbtagstclcorr(); }
	int &nbtagstcmcorr() { return stopt.nbtagstcmcorr(); }
	int &nbtagscsvlcorr() { return stopt.nbtagscsvlcorr(); }
	int &nbtagscsvmcorr() { return stopt.nbtagscsvmcorr(); }
	int &nbtagscsvtcott() { return stopt.nbtagscsvtcott(); }
	int &njetsUp() { return stopt.njetsUp(); }
	int &njetsDown() { return stopt.njetsDown(); }
	float &htUp() { return stopt.htUp(); }
	float &htDown() { return stopt.htDown(); }
	int &ntruepu() { return stopt.ntruepu(); }
	int &npu() { return stopt.npu(); }
	int &npuMinusOne() { return stopt.npuMinusOne(); }
	int &npuPlusOne() { return stopt.npuPlusOne(); }
	int &nvtx() { return stopt.nvtx(); }
	int &indexfirstGoodVertex_() { return stopt.indexfirstGoodVertex_(); }
	float &nvtxweight() { return stopt.nvtxweight(); }
	float &n3dvtxweight() { return stopt.n3dvtxweight(); }
	int &pdfid1() { return stopt.pdfid1(); }
	int &pdfid2() { return stopt.pdfid2(); }
	float &pdfx1() { return stopt.pdfx1(); }
	float &pdfx2() { return stopt.pdfx2(); }
	float &pdfQ() { return stopt.pdfQ(); }
	float &vecjetpt() { return stopt.vecjetpt(); }
	int &pass() { return stopt.pass(); }
	int &passz() { return stopt.passz(); }
	float &m0() { return stopt.m0(); }
	float &mg() { return stopt.mg(); }
	float &ml() { return stopt.ml(); }
	float &x() { return stopt.x(); }
	float &m12() { return stopt.m12(); }
	float &lep1chi2ndf() { return stopt.lep1chi2ndf(); }
	float &lep2chi2ndf() { return stopt.lep2chi2ndf(); }
	float &lep1dpt() { return stopt.lep1dpt(); }
	float &lep2dpt() { return stopt.lep2dpt(); }
	int &id1() { return stopt.id1(); }
	int &id2() { return stopt.id2(); }
	int &leptype1() { return stopt.leptype1(); }
	int &leptype2() { return stopt.leptype2(); }
	int &w1() { return stopt.w1(); }
	int &w2() { return stopt.w2(); }
	float &iso1() { return stopt.iso1(); }
	float &isont1() { return stopt.isont1(); }
	float &isopfold1() { return stopt.isopfold1(); }
	float &isopf1() { return stopt.isopf1(); }
	float &etasc1() { return stopt.etasc1(); }
	float &etasc2() { return stopt.etasc2(); }
	float &eoverpin() { return stopt.eoverpin(); }
	float &eoverpout() { return stopt.eoverpout(); }
	float &dEtaIn() { return stopt.dEtaIn(); }
	float &dPhiIn() { return stopt.dPhiIn(); }
	float &sigmaIEtaIEta() { return stopt.sigmaIEtaIEta(); }
	float &hOverE() { return stopt.hOverE(); }
	float &ooemoop() { return stopt.ooemoop(); }
	float &d0vtx() { return stopt.d0vtx(); }
	float &dzvtx() { return stopt.dzvtx(); }
	float &expinnerlayers() { return stopt.expinnerlayers(); }
	float &fbrem() { return stopt.fbrem(); }
	float &pfisoch() { return stopt.pfisoch(); }
	float &pfisoem() { return stopt.pfisoem(); }
	float &pfisonh() { return stopt.pfisonh(); }
	float &eSC() { return stopt.eSC(); }
	float &phiSC() { return stopt.phiSC(); }
	float &eSCRaw() { return stopt.eSCRaw(); }
	float &eSCPresh() { return stopt.eSCPresh(); }
	float &lep1_scslasercormean() { return stopt.lep1_scslasercormean(); }
	float &lep1_scslasercormax() { return stopt.lep1_scslasercormax(); }
	float &eoverpin2() { return stopt.eoverpin2(); }
	float &eoverpout2() { return stopt.eoverpout2(); }
	float &dEtaIn2() { return stopt.dEtaIn2(); }
	float &dPhiIn2() { return stopt.dPhiIn2(); }
	float &sigmaIEtaIEta2() { return stopt.sigmaIEtaIEta2(); }
	float &hOverE2() { return stopt.hOverE2(); }
	float &ooemoop2() { return stopt.ooemoop2(); }
	float &d0vtx2() { return stopt.d0vtx2(); }
	float &dzvtx2() { return stopt.dzvtx2(); }
	float &expinnerlayers2() { return stopt.expinnerlayers2(); }
	float &fbrem2() { return stopt.fbrem2(); }
	float &pfisoch2() { return stopt.pfisoch2(); }
	float &pfisoem2() { return stopt.pfisoem2(); }
	float &pfisonh2() { return stopt.pfisonh2(); }
	float &eSC2() { return stopt.eSC2(); }
	float &phiSC2() { return stopt.phiSC2(); }
	float &eSCRaw2() { return stopt.eSCRaw2(); }
	float &eSCPresh2() { return stopt.eSCPresh2(); }
	float &lep2_scslasercormean() { return stopt.lep2_scslasercormean(); }
	float &lep2_scslasercormax() { return stopt.lep2_scslasercormax(); }
	float &scslasercormax() { return stopt.scslasercormax(); }
	float &scslasercormax_pt() { return stopt.scslasercormax_pt(); }
	float &scslasercormax_eta() { return stopt.scslasercormax_eta(); }
	float &iso2() { return stopt.iso2(); }
	float &ecalveto1() { return stopt.ecalveto1(); }
	float &ecalveto2() { return stopt.ecalveto2(); }
	float &hcalveto1() { return stopt.hcalveto1(); }
	float &hcalveto2() { return stopt.hcalveto2(); }
	float &isont2() { return stopt.isont2(); }
	float &isopf2() { return stopt.isopf2(); }
	float &ptl1() { return stopt.ptl1(); }
	float &ptl2() { return stopt.ptl2(); }
	float &etal1() { return stopt.etal1(); }
	float &etal2() { return stopt.etal2(); }
	float &phil1() { return stopt.phil1(); }
	float &phil2() { return stopt.phil2(); }
	float &meff() { return stopt.meff(); }
	float &mt() { return stopt.mt(); }
	int &run() { return stopt.run(); }
	int &lumi() { return stopt.lumi(); }
	int &event() { return stopt.event(); }
	float &y() { return stopt.y(); }
	float &ht() { return stopt.ht(); }
	float &htgen() { return stopt.htgen(); }
	float &htjpt() { return stopt.htjpt(); }
	int &nels() { return stopt.nels(); }
	int &nmus() { return stopt.nmus(); }
	int &ntaus() { return stopt.ntaus(); }
	int &nleps() { return stopt.nleps(); }
	int &nbs() { return stopt.nbs(); }
	float &dphijm() { return stopt.dphijm(); }
	float &ptjetraw() { return stopt.ptjetraw(); }
	float &ptjet23() { return stopt.ptjet23(); }
	float &ptjetF23() { return stopt.ptjetF23(); }
	float &ptjetO23() { return stopt.ptjetO23(); }
	int &mcid1() { return stopt.mcid1(); }
	float &mcdr1() { return stopt.mcdr1(); }
	int &mcdecay1() { return stopt.mcdecay1(); }
	int &mcndec1() { return stopt.mcndec1(); }
	int &mcndec2() { return stopt.mcndec2(); }
	int &mcndeckls1() { return stopt.mcndeckls1(); }
	int &mcndeckls2() { return stopt.mcndeckls2(); }
	int &mcndecem1() { return stopt.mcndecem1(); }
	int &mcndecem2() { return stopt.mcndecem2(); }
	int &mcid2() { return stopt.mcid2(); }
	float &mcdr2() { return stopt.mcdr2(); }
	int &mcdecay2() { return stopt.mcdecay2(); }
	float &mctaudpt1() { return stopt.mctaudpt1(); }
	float &mctaudpt2() { return stopt.mctaudpt2(); }
	int &mctaudid1() { return stopt.mctaudid1(); }
	int &mctaudid2() { return stopt.mctaudid2(); }
	int &mlepid() { return stopt.mlepid(); }
	int &mleppassid() { return stopt.mleppassid(); }
	int &mleppassiso() { return stopt.mleppassiso(); }
	float &mlepiso() { return stopt.mlepiso(); }
	float &mlepdr() { return stopt.mlepdr(); }
	float &pflepiso() { return stopt.pflepiso(); }
	float &pflepdr() { return stopt.pflepdr(); }
	float &pfleppt() { return stopt.pfleppt(); }
	float &pflepmindrj() { return stopt.pflepmindrj(); }
	float &pftaudiso() { return stopt.pftaudiso(); }
	float &pftauddr() { return stopt.pftauddr(); }
	float &pftaudpt() { return stopt.pftaudpt(); }
	float &pftaudmindrj() { return stopt.pftaudmindrj(); }
	int &pfcandid5() { return stopt.pfcandid5(); }
	float &pfcandiso5() { return stopt.pfcandiso5(); }
	float &pfcandpt5() { return stopt.pfcandpt5(); }
	float &pfcanddz5() { return stopt.pfcanddz5(); }
	float &pfcandmindrj5() { return stopt.pfcandmindrj5(); }
	int &pfcandid10() { return stopt.pfcandid10(); }
	float &pfcandiso10() { return stopt.pfcandiso10(); }
	float &pfcandpt10() { return stopt.pfcandpt10(); }
	float &pfcanddz10() { return stopt.pfcanddz10(); }
	float &pfcandmindrj10() { return stopt.pfcandmindrj10(); }
	int &pfcandidOS10() { return stopt.pfcandidOS10(); }
	float &pfcandisoOS10() { return stopt.pfcandisoOS10(); }
	float &pfcandptOS10() { return stopt.pfcandptOS10(); }
	float &pfcanddzOS10() { return stopt.pfcanddzOS10(); }
	int &pfcandid5looseZ() { return stopt.pfcandid5looseZ(); }
	float &pfcandiso5looseZ() { return stopt.pfcandiso5looseZ(); }
	float &pfcandpt5looseZ() { return stopt.pfcandpt5looseZ(); }
	float &pfcanddz5looseZ() { return stopt.pfcanddz5looseZ(); }
	int &pfcandidOS10looseZ() { return stopt.pfcandidOS10looseZ(); }
	float &pfcandisoOS10looseZ() { return stopt.pfcandisoOS10looseZ(); }
	float &pfcandptOS10looseZ() { return stopt.pfcandptOS10looseZ(); }
	float &pfcanddzOS10looseZ() { return stopt.pfcanddzOS10looseZ(); }
	int &pfcanddirid10() { return stopt.pfcanddirid10(); }
	float &pfcanddiriso10() { return stopt.pfcanddiriso10(); }
	float &pfcanddirpt10() { return stopt.pfcanddirpt10(); }
	float &pfcanddirmindrj10() { return stopt.pfcanddirmindrj10(); }
	int &pfcandvetoid10() { return stopt.pfcandvetoid10(); }
	float &pfcandvetoiso10() { return stopt.pfcandvetoiso10(); }
	float &pfcandvetopt10() { return stopt.pfcandvetopt10(); }
	float &pfcandvetomindrj10() { return stopt.pfcandvetomindrj10(); }
	int &pfcandvetoLid10() { return stopt.pfcandvetoLid10(); }
	float &pfcandvetoLiso10() { return stopt.pfcandvetoLiso10(); }
	float &pfcandvetoLpt10() { return stopt.pfcandvetoLpt10(); }
	float &pfcandvetoLmindrj10() { return stopt.pfcandvetoLmindrj10(); }
	float &emjet10() { return stopt.emjet10(); }
	float &mjj() { return stopt.mjj(); }
	float &emjet20() { return stopt.emjet20(); }
	float &trkpt5() { return stopt.trkpt5(); }
	float &trkpt10() { return stopt.trkpt10(); }
	float &mleptrk5() { return stopt.mleptrk5(); }
	float &mleptrk10() { return stopt.mleptrk10(); }
	float &trkreliso5() { return stopt.trkreliso5(); }
	float &trkreliso10() { return stopt.trkreliso10(); }
	float &trkpt5loose() { return stopt.trkpt5loose(); }
	float &trkpt10loose() { return stopt.trkpt10loose(); }
	float &trkreliso5loose() { return stopt.trkreliso5loose(); }
	float &trkreliso10loose() { return stopt.trkreliso10loose(); }
	float &trkpt10pt0p1() { return stopt.trkpt10pt0p1(); }
	float &trkpt10pt0p2() { return stopt.trkpt10pt0p2(); }
	float &trkpt10pt0p3() { return stopt.trkpt10pt0p3(); }
	float &trkpt10pt0p4() { return stopt.trkpt10pt0p4(); }
	float &trkpt10pt0p5() { return stopt.trkpt10pt0p5(); }
	float &trkpt10pt0p6() { return stopt.trkpt10pt0p6(); }
	float &trkpt10pt0p7() { return stopt.trkpt10pt0p7(); }
	float &trkpt10pt0p8() { return stopt.trkpt10pt0p8(); }
	float &trkpt10pt0p9() { return stopt.trkpt10pt0p9(); }
	float &trkpt10pt1p0() { return stopt.trkpt10pt1p0(); }
	float &trkreliso10pt0p1() { return stopt.trkreliso10pt0p1(); }
	float &trkreliso10pt0p2() { return stopt.trkreliso10pt0p2(); }
	float &trkreliso10pt0p3() { return stopt.trkreliso10pt0p3(); }
	float &trkreliso10pt0p4() { return stopt.trkreliso10pt0p4(); }
	float &trkreliso10pt0p5() { return stopt.trkreliso10pt0p5(); }
	float &trkreliso10pt0p6() { return stopt.trkreliso10pt0p6(); }
	float &trkreliso10pt0p7() { return stopt.trkreliso10pt0p7(); }
	float &trkreliso10pt0p8() { return stopt.trkreliso10pt0p8(); }
	float &trkreliso10pt0p9() { return stopt.trkreliso10pt0p9(); }
	float &trkreliso10pt1p0() { return stopt.trkreliso10pt1p0(); }
	float &pfcandpt10pt0p1() { return stopt.pfcandpt10pt0p1(); }
	float &pfcandpt10pt0p2() { return stopt.pfcandpt10pt0p2(); }
	float &pfcandpt10pt0p3() { return stopt.pfcandpt10pt0p3(); }
	float &pfcandpt10pt0p4() { return stopt.pfcandpt10pt0p4(); }
	float &pfcandpt10pt0p5() { return stopt.pfcandpt10pt0p5(); }
	float &pfcandpt10pt0p6() { return stopt.pfcandpt10pt0p6(); }
	float &pfcandpt10pt0p7() { return stopt.pfcandpt10pt0p7(); }
	float &pfcandpt10pt0p8() { return stopt.pfcandpt10pt0p8(); }
	float &pfcandpt10pt0p9() { return stopt.pfcandpt10pt0p9(); }
	float &pfcandpt10pt1p0() { return stopt.pfcandpt10pt1p0(); }
	float &pfcandiso10pt0p1() { return stopt.pfcandiso10pt0p1(); }
	float &pfcandiso10pt0p2() { return stopt.pfcandiso10pt0p2(); }
	float &pfcandiso10pt0p3() { return stopt.pfcandiso10pt0p3(); }
	float &pfcandiso10pt0p4() { return stopt.pfcandiso10pt0p4(); }
	float &pfcandiso10pt0p5() { return stopt.pfcandiso10pt0p5(); }
	float &pfcandiso10pt0p6() { return stopt.pfcandiso10pt0p6(); }
	float &pfcandiso10pt0p7() { return stopt.pfcandiso10pt0p7(); }
	float &pfcandiso10pt0p8() { return stopt.pfcandiso10pt0p8(); }
	float &pfcandiso10pt0p9() { return stopt.pfcandiso10pt0p9(); }
	float &pfcandiso10pt1p0() { return stopt.pfcandiso10pt1p0(); }
	float &mbb() { return stopt.mbb(); }
	float &lep1pfjetdr() { return stopt.lep1pfjetdr(); }
	float &lep2pfjetdr() { return stopt.lep2pfjetdr(); }
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &mclep() { return stopt.mclep(); }
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &mcnu() { return stopt.mcnu(); }
	float &mcmln() { return stopt.mcmln(); }
	float &mcmtln() { return stopt.mcmtln(); }
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &mlep() { return stopt.mlep(); }
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &lep1() { return stopt.lep1(); }
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &lep2() { return stopt.lep2(); }
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &trklep1() { return stopt.trklep1(); }
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &trklep2() { return stopt.trklep2(); }
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &gfitlep1() { return stopt.gfitlep1(); }
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &gfitlep2() { return stopt.gfitlep2(); }
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &lepp() { return stopt.lepp(); }
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &lepm() { return stopt.lepm(); }
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &pflep1() { return stopt.pflep1(); }
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &pflep2() { return stopt.pflep2(); }
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &leppfjet1() { return stopt.leppfjet1(); }
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &leppfjet2() { return stopt.leppfjet2(); }
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &mclep1() { return stopt.mclep1(); }
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &mclep2() { return stopt.mclep2(); }
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &mctaud1() { return stopt.mctaud1(); }
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &mctaud2() { return stopt.mctaud2(); }
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &mctaudvis1() { return stopt.mctaudvis1(); }
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &mctaudvis2() { return stopt.mctaudvis2(); }
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &pflep() { return stopt.pflep(); }
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &pftaud() { return stopt.pftaud(); }
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &pfcand5() { return stopt.pfcand5(); }
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &pfcand10() { return stopt.pfcand10(); }
	int &pfTau15_leadPtcandID() { return stopt.pfTau15_leadPtcandID(); }
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &pfTau15() { return stopt.pfTau15(); }
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &pfTau15_leadPtcand() { return stopt.pfTau15_leadPtcand(); }
	int &pfTau_leadPtcandID() { return stopt.pfTau_leadPtcandID(); }
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &pfTau() { return stopt.pfTau(); }
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &pfTau_leadPtcand() { return stopt.pfTau_leadPtcand(); }
	int &pfTauLoose_leadPtcandID() { return stopt.pfTauLoose_leadPtcandID(); }
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &pfTauLoose() { return stopt.pfTauLoose(); }
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &pfTauLoose_leadPtcand() { return stopt.pfTauLoose_leadPtcand(); }
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &pfcandOS10() { return stopt.pfcandOS10(); }
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &pfcandOS10looseZ() { return stopt.pfcandOS10looseZ(); }
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &pfcand5looseZ() { return stopt.pfcand5looseZ(); }
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &pfcanddir10() { return stopt.pfcanddir10(); }
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &pfcandveto10() { return stopt.pfcandveto10(); }
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &pfcandvetoL10() { return stopt.pfcandvetoL10(); }
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &jet() { return stopt.jet(); }
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &nonisoel() { return stopt.nonisoel(); }
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &nonisomu() { return stopt.nonisomu(); }
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &t() { return stopt.t(); }
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &tbar() { return stopt.tbar(); }
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &ttbar() { return stopt.ttbar(); }
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &lep_t() { return stopt.lep_t(); }
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &lep_tbar() { return stopt.lep_tbar(); }
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &stop_t() { return stopt.stop_t(); }
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &stop_tbar() { return stopt.stop_tbar(); }
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &neutralino_t() { return stopt.neutralino_t(); }
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &neutralino_tbar() { return stopt.neutralino_tbar(); }
	int &lep_t_id() { return stopt.lep_t_id(); }
	int &lep_tbar_id() { return stopt.lep_tbar_id(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &pfjets() { return stopt.pfjets(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &pfjets_genJet_() { return stopt.pfjets_genJet_(); }
	vector<float> &pfjets_csv() { return stopt.pfjets_csv(); }
	vector<float> &pfjets_chm() { return stopt.pfjets_chm(); }
	vector<float> &pfjets_neu() { return stopt.pfjets_neu(); }
	vector<float> &pfjets_l1corr() { return stopt.pfjets_l1corr(); }
	vector<float> &pfjets_corr() { return stopt.pfjets_corr(); }
	vector<int> &pfjets_mc3() { return stopt.pfjets_mc3(); }
	vector<int> &pfjets_mcflavorAlgo() { return stopt.pfjets_mcflavorAlgo(); }
	vector<int> &pfjets_mcflavorPhys() { return stopt.pfjets_mcflavorPhys(); }
	vector<float> &pfjets_uncertainty() { return stopt.pfjets_uncertainty(); }
	vector<int> &pfjets_flav() { return stopt.pfjets_flav(); }
	vector<float> &pfjets_lrm() { return stopt.pfjets_lrm(); }
	vector<float> &pfjets_lrm2() { return stopt.pfjets_lrm2(); }
	vector<float> &pfjets_qgtag() { return stopt.pfjets_qgtag(); }
	vector<float> &pfjets_genJetDr() { return stopt.pfjets_genJetDr(); }
	vector<float> &pfjets_sigma() { return stopt.pfjets_sigma(); }
	vector<int> &pfjets_lepjet() { return stopt.pfjets_lepjet(); }
	vector<float> &pfjets_beta() { return stopt.pfjets_beta(); }
	vector<float> &pfjets_beta2() { return stopt.pfjets_beta2(); }
	vector<float> &pfjets_beta_0p1() { return stopt.pfjets_beta_0p1(); }
	vector<float> &pfjets_beta_0p2() { return stopt.pfjets_beta_0p2(); }
	vector<float> &pfjets_beta2_0p1() { return stopt.pfjets_beta2_0p1(); }
	vector<float> &pfjets_beta2_0p5() { return stopt.pfjets_beta2_0p5(); }
	vector<float> &pfjets_mvaPUid() { return stopt.pfjets_mvaPUid(); }
	vector<float> &pfjets_mva5xPUid() { return stopt.pfjets_mva5xPUid(); }
	vector<float> &pfjets_mvaBeta() { return stopt.pfjets_mvaBeta(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &genbs() { return stopt.genbs(); }
	vector<int> &genps_pdgId() { return stopt.genps_pdgId(); }
	vector<int> &genps_firstMother() { return stopt.genps_firstMother(); }
	vector<float> &genps_energy() { return stopt.genps_energy(); }
	vector<float> &genps_pt() { return stopt.genps_pt(); }
	vector<float> &genps_eta() { return stopt.genps_eta(); }
	vector<float> &genps_phi() { return stopt.genps_phi(); }
}
