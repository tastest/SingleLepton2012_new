// -*- C++ -*-
#ifndef STOPT_H
#define STOPT_H
#include "Math/LorentzVector.h"
#include "Math/Point3D.h"
#include "TMath.h"
#include "TBranch.h"
#include "TTree.h"
#include "TH1F.h"
#include "TFile.h"
#include "TBits.h"
#include <vector> 
typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > LorentzVector;

using namespace std; 
class STOPT {
private: 
protected: 
	unsigned int index;
	int	acc_2010_;
	TBranch *acc_2010_branch;
	bool acc_2010_isLoaded;
	int	acc_highmet_;
	TBranch *acc_highmet_branch;
	bool acc_highmet_isLoaded;
	int	acc_highht_;
	TBranch *acc_highht_branch;
	bool acc_highht_isLoaded;
	int	eldup_;
	TBranch *eldup_branch;
	bool eldup_isLoaded;
	int	csc_;
	TBranch *csc_branch;
	bool csc_isLoaded;
	int	hbhe_;
	TBranch *hbhe_branch;
	bool hbhe_isLoaded;
	int	hbhenew_;
	TBranch *hbhenew_branch;
	bool hbhenew_isLoaded;
	int	hcallaser_;
	TBranch *hcallaser_branch;
	bool hcallaser_isLoaded;
	int	ecaltp_;
	TBranch *ecaltp_branch;
	bool ecaltp_isLoaded;
	int	trkfail_;
	TBranch *trkfail_branch;
	bool trkfail_isLoaded;
	int	eebadsc_;
	TBranch *eebadsc_branch;
	bool eebadsc_isLoaded;
	int	lep1_badecallaser_;
	TBranch *lep1_badecallaser_branch;
	bool lep1_badecallaser_isLoaded;
	int	lep2_badecallaser_;
	TBranch *lep2_badecallaser_branch;
	bool lep2_badecallaser_isLoaded;
	int	isdata_;
	TBranch *isdata_branch;
	bool isdata_isLoaded;
	int	jetid_;
	TBranch *jetid_branch;
	bool jetid_isLoaded;
	int	jetid30_;
	TBranch *jetid30_branch;
	bool jetid30_isLoaded;
	int	json_;
	TBranch *json_branch;
	bool json_isLoaded;
	float	htoffset_;
	TBranch *htoffset_branch;
	bool htoffset_isLoaded;
	float	htuncor_;
	TBranch *htuncor_branch;
	bool htuncor_isLoaded;
	float	ptt_;
	TBranch *ptt_branch;
	bool ptt_isLoaded;
	float	pttbar_;
	TBranch *pttbar_branch;
	bool pttbar_isLoaded;
	float	ptttbar_;
	TBranch *ptttbar_branch;
	bool ptttbar_isLoaded;
	float	mttbar_;
	TBranch *mttbar_branch;
	bool mttbar_isLoaded;
	int	npartons_;
	TBranch *npartons_branch;
	bool npartons_isLoaded;
	int	nwzpartons_;
	TBranch *nwzpartons_branch;
	bool nwzpartons_isLoaded;
	int	hyptype_;
	TBranch *hyptype_branch;
	bool hyptype_isLoaded;
	float	maxpartonpt_;
	TBranch *maxpartonpt_branch;
	bool maxpartonpt_isLoaded;
	float	etattbar_;
	TBranch *etattbar_branch;
	bool etattbar_isLoaded;
	int	njetsoffset_;
	TBranch *njetsoffset_branch;
	bool njetsoffset_isLoaded;
	int	njetsuncor_;
	TBranch *njetsuncor_branch;
	bool njetsuncor_isLoaded;
	float	costhetaweight_;
	TBranch *costhetaweight_branch;
	bool costhetaweight_isLoaded;
	float	weight_;
	TBranch *weight_branch;
	bool weight_isLoaded;
	float	weightleft_;
	TBranch *weightleft_branch;
	bool weightleft_isLoaded;
	float	weightright_;
	TBranch *weightright_branch;
	bool weightright_isLoaded;
	float	mutrigweight_;
	TBranch *mutrigweight_branch;
	bool mutrigweight_isLoaded;
	float	mutrigweight2_;
	TBranch *mutrigweight2_branch;
	bool mutrigweight2_isLoaded;
	float	sltrigweight_;
	TBranch *sltrigweight_branch;
	bool sltrigweight_isLoaded;
	float	dltrigweight_;
	TBranch *dltrigweight_branch;
	bool dltrigweight_isLoaded;
	float	trgeff_;
	TBranch *trgeff_branch;
	bool trgeff_isLoaded;
	float	pthat_;
	TBranch *pthat_branch;
	bool pthat_isLoaded;
	float	qscale_;
	TBranch *qscale_branch;
	bool qscale_isLoaded;
	float	mgcor_;
	TBranch *mgcor_branch;
	bool mgcor_isLoaded;
	int	wflav_;
	TBranch *wflav_branch;
	bool wflav_isLoaded;
	float	ksusy_;
	TBranch *ksusy_branch;
	bool ksusy_isLoaded;
	float	ksusyup_;
	TBranch *ksusyup_branch;
	bool ksusyup_isLoaded;
	float	ksusydn_;
	TBranch *ksusydn_branch;
	bool ksusydn_isLoaded;
	float	xsecsusy_;
	TBranch *xsecsusy_branch;
	bool xsecsusy_isLoaded;
	float	xsecsusy2_;
	TBranch *xsecsusy2_branch;
	bool xsecsusy2_isLoaded;
	float	smeff_;
	TBranch *smeff_branch;
	bool smeff_isLoaded;
	float	k_;
	TBranch *k_branch;
	bool k_isLoaded;
	float	mllgen_;
	TBranch *mllgen_branch;
	bool mllgen_isLoaded;
	float	ptwgen_;
	TBranch *ptwgen_branch;
	bool ptwgen_isLoaded;
	float	ptzgen_;
	TBranch *ptzgen_branch;
	bool ptzgen_isLoaded;
	int	nlep_;
	TBranch *nlep_branch;
	bool nlep_isLoaded;
	int	nosel_;
	TBranch *nosel_branch;
	bool nosel_isLoaded;
	int	ngoodlep_;
	TBranch *ngoodlep_branch;
	bool ngoodlep_isLoaded;
	int	ngoodel_;
	TBranch *ngoodel_branch;
	bool ngoodel_isLoaded;
	int	ngoodmu_;
	TBranch *ngoodmu_branch;
	bool ngoodmu_isLoaded;
	int	mull_;
	TBranch *mull_branch;
	bool mull_isLoaded;
	int	mult_;
	TBranch *mult_branch;
	bool mult_isLoaded;
	int	mullgen_;
	TBranch *mullgen_branch;
	bool mullgen_isLoaded;
	int	multgen_;
	TBranch *multgen_branch;
	bool multgen_isLoaded;
	int	proc_;
	TBranch *proc_branch;
	bool proc_isLoaded;
	int	leptype_;
	TBranch *leptype_branch;
	bool leptype_isLoaded;
	float	topmass_;
	TBranch *topmass_branch;
	bool topmass_isLoaded;
	float	dilmass_;
	TBranch *dilmass_branch;
	bool dilmass_isLoaded;
	float	dilrecoil_;
	TBranch *dilrecoil_branch;
	bool dilrecoil_isLoaded;
	float	dilrecoilparl_;
	TBranch *dilrecoilparl_branch;
	bool dilrecoilparl_isLoaded;
	float	dilrecoilperp_;
	TBranch *dilrecoilperp_branch;
	bool dilrecoilperp_isLoaded;
	float	tcmet_;
	TBranch *tcmet_branch;
	bool tcmet_isLoaded;
	float	genmet_;
	TBranch *genmet_branch;
	bool genmet_isLoaded;
	float	gensumet_;
	TBranch *gensumet_branch;
	bool gensumet_isLoaded;
	float	genmetphi_;
	TBranch *genmetphi_branch;
	bool genmetphi_isLoaded;
	float	calomet_;
	TBranch *calomet_branch;
	bool calomet_isLoaded;
	float	calometphi_;
	TBranch *calometphi_branch;
	bool calometphi_isLoaded;
	float	trkmet_;
	TBranch *trkmet_branch;
	bool trkmet_isLoaded;
	float	trkmetphi_;
	TBranch *trkmetphi_branch;
	bool trkmetphi_isLoaded;
	float	pfmet_;
	TBranch *pfmet_branch;
	bool pfmet_isLoaded;
	float	pfmetveto_;
	TBranch *pfmetveto_branch;
	bool pfmetveto_isLoaded;
	float	pfmetsig_;
	TBranch *pfmetsig_branch;
	bool pfmetsig_isLoaded;
	float	pfmetphi_;
	TBranch *pfmetphi_branch;
	bool pfmetphi_isLoaded;
	float	pfsumet_;
	TBranch *pfsumet_branch;
	bool pfsumet_isLoaded;
	float	mucormet_;
	TBranch *mucormet_branch;
	bool mucormet_isLoaded;
	float	mucorjesmet_;
	TBranch *mucorjesmet_branch;
	bool mucorjesmet_isLoaded;
	float	tcmet35X_;
	TBranch *tcmet35X_branch;
	bool tcmet35X_isLoaded;
	float	tcmetevent_;
	TBranch *tcmetevent_branch;
	bool tcmetevent_isLoaded;
	float	tcmetlooper_;
	TBranch *tcmetlooper_branch;
	bool tcmetlooper_isLoaded;
	float	tcmetphi_;
	TBranch *tcmetphi_branch;
	bool tcmetphi_isLoaded;
	float	tcsumet_;
	TBranch *tcsumet_branch;
	bool tcsumet_isLoaded;
	float	tcmetUp_;
	TBranch *tcmetUp_branch;
	bool tcmetUp_isLoaded;
	float	tcmetDown_;
	TBranch *tcmetDown_branch;
	bool tcmetDown_isLoaded;
	float	tcmetTest_;
	TBranch *tcmetTest_branch;
	bool tcmetTest_isLoaded;
	float	pfmetUp_;
	TBranch *pfmetUp_branch;
	bool pfmetUp_isLoaded;
	float	pfmetDown_;
	TBranch *pfmetDown_branch;
	bool pfmetDown_isLoaded;
	float	pfmetTest_;
	TBranch *pfmetTest_branch;
	bool pfmetTest_isLoaded;
	float	sumjetpt_;
	TBranch *sumjetpt_branch;
	bool sumjetpt_isLoaded;
	float	dileta_;
	TBranch *dileta_branch;
	bool dileta_isLoaded;
	float	dilpt_;
	TBranch *dilpt_branch;
	bool dilpt_isLoaded;
	float	dildphi_;
	TBranch *dildphi_branch;
	bool dildphi_isLoaded;
	int	ngenjets_;
	TBranch *ngenjets_branch;
	bool ngenjets_isLoaded;
	int	njpt_;
	TBranch *njpt_branch;
	bool njpt_isLoaded;
	int	trgmu1_;
	TBranch *trgmu1_branch;
	bool trgmu1_isLoaded;
	int	trgmu2_;
	TBranch *trgmu2_branch;
	bool trgmu2_isLoaded;
	int	trgel1_;
	TBranch *trgel1_branch;
	bool trgel1_isLoaded;
	int	trgel2_;
	TBranch *trgel2_branch;
	bool trgel2_isLoaded;
	int	isomu24_;
	TBranch *isomu24_branch;
	bool isomu24_isLoaded;
	int	ele27wp80_;
	TBranch *ele27wp80_branch;
	bool ele27wp80_isLoaded;
	int	mm_;
	TBranch *mm_branch;
	bool mm_isLoaded;
	int	mmtk_;
	TBranch *mmtk_branch;
	bool mmtk_isLoaded;
	int	me_;
	TBranch *me_branch;
	bool me_isLoaded;
	int	em_;
	TBranch *em_branch;
	bool em_isLoaded;
	int	mu_;
	TBranch *mu_branch;
	bool mu_isLoaded;
	int	ee_;
	TBranch *ee_branch;
	bool ee_isLoaded;
	int	npfjets30_;
	TBranch *npfjets30_branch;
	bool npfjets30_isLoaded;
	int	npfjets30lepcorr_;
	TBranch *npfjets30lepcorr_branch;
	bool npfjets30lepcorr_isLoaded;
	float	knjets_;
	TBranch *knjets_branch;
	bool knjets_isLoaded;
	float	rhovor_;
	TBranch *rhovor_branch;
	bool rhovor_isLoaded;
	float	htpf30_;
	TBranch *htpf30_branch;
	bool htpf30_isLoaded;
	float	t1met10_;
	TBranch *t1met10_branch;
	bool t1met10_isLoaded;
	float	t1met20_;
	TBranch *t1met20_branch;
	bool t1met20_isLoaded;
	float	t1met30_;
	TBranch *t1met30_branch;
	bool t1met30_isLoaded;
	float	t1met10phi_;
	TBranch *t1met10phi_branch;
	bool t1met10phi_isLoaded;
	float	t1met20phi_;
	TBranch *t1met20phi_branch;
	bool t1met20phi_isLoaded;
	float	t1met30phi_;
	TBranch *t1met30phi_branch;
	bool t1met30phi_isLoaded;
	float	t1met10mt_;
	TBranch *t1met10mt_branch;
	bool t1met10mt_isLoaded;
	float	t1met20mt_;
	TBranch *t1met20mt_branch;
	bool t1met20mt_isLoaded;
	float	t1met30mt_;
	TBranch *t1met30mt_branch;
	bool t1met30mt_isLoaded;
	float	lepmetpt_;
	TBranch *lepmetpt_branch;
	bool lepmetpt_isLoaded;
	float	lept1met10pt_;
	TBranch *lept1met10pt_branch;
	bool lept1met10pt_isLoaded;
	float	t1met10s_;
	TBranch *t1met10s_branch;
	bool t1met10s_isLoaded;
	float	t1met10sphi_;
	TBranch *t1met10sphi_branch;
	bool t1met10sphi_isLoaded;
	float	t1met10smt_;
	TBranch *t1met10smt_branch;
	bool t1met10smt_isLoaded;
	float	t1metphicorr_;
	TBranch *t1metphicorr_branch;
	bool t1metphicorr_isLoaded;
	float	t1metphicorrup_;
	TBranch *t1metphicorrup_branch;
	bool t1metphicorrup_isLoaded;
	float	t1metphicorrdn_;
	TBranch *t1metphicorrdn_branch;
	bool t1metphicorrdn_isLoaded;
	float	t1metphicorrphi_;
	TBranch *t1metphicorrphi_branch;
	bool t1metphicorrphi_isLoaded;
	float	t1metphicorrphiup_;
	TBranch *t1metphicorrphiup_branch;
	bool t1metphicorrphiup_isLoaded;
	float	t1metphicorrphidn_;
	TBranch *t1metphicorrphidn_branch;
	bool t1metphicorrphidn_isLoaded;
	float	t1metphicorrlep_;
	TBranch *t1metphicorrlep_branch;
	bool t1metphicorrlep_isLoaded;
	float	t1metphicorrlepphi_;
	TBranch *t1metphicorrlepphi_branch;
	bool t1metphicorrlepphi_isLoaded;
	float	t1metphicorrmt_;
	TBranch *t1metphicorrmt_branch;
	bool t1metphicorrmt_isLoaded;
	float	t1metphicorrmtup_;
	TBranch *t1metphicorrmtup_branch;
	bool t1metphicorrmtup_isLoaded;
	float	t1metphicorrmtdn_;
	TBranch *t1metphicorrmtdn_branch;
	bool t1metphicorrmtdn_isLoaded;
	float	t1metphicorrlepmt_;
	TBranch *t1metphicorrlepmt_branch;
	bool t1metphicorrlepmt_isLoaded;
	float	t1met_off_;
	TBranch *t1met_off_branch;
	bool t1met_off_isLoaded;
	float	t1metphi_off_;
	TBranch *t1metphi_off_branch;
	bool t1metphi_off_isLoaded;
	float	t1metmt_off_;
	TBranch *t1metmt_off_branch;
	bool t1metmt_off_isLoaded;
	float	t1metphicorr_off_;
	TBranch *t1metphicorr_off_branch;
	bool t1metphicorr_off_isLoaded;
	float	t1metphicorrphi_off_;
	TBranch *t1metphicorrphi_off_branch;
	bool t1metphicorrphi_off_isLoaded;
	float	t1metphicorrmt_off_;
	TBranch *t1metphicorrmt_off_branch;
	bool t1metphicorrmt_off_isLoaded;
	float	mht15_;
	TBranch *mht15_branch;
	bool mht15_isLoaded;
	float	mht15phi_;
	TBranch *mht15phi_branch;
	bool mht15phi_isLoaded;
	float	trkmet_mht15_;
	TBranch *trkmet_mht15_branch;
	bool trkmet_mht15_isLoaded;
	float	trkmetphi_mht15_;
	TBranch *trkmetphi_mht15_branch;
	bool trkmetphi_mht15_isLoaded;
	float	mettlj15_;
	TBranch *mettlj15_branch;
	bool mettlj15_isLoaded;
	float	mettlj15phi_;
	TBranch *mettlj15phi_branch;
	bool mettlj15phi_isLoaded;
	float	mt2bmin_;
	TBranch *mt2bmin_branch;
	bool mt2bmin_isLoaded;
	float	mt2blmin_;
	TBranch *mt2blmin_branch;
	bool mt2blmin_isLoaded;
	float	mt2wmin_;
	TBranch *mt2wmin_branch;
	bool mt2wmin_isLoaded;
	float	chi2min_;
	TBranch *chi2min_branch;
	bool chi2min_isLoaded;
	float	chi2minprob_;
	TBranch *chi2minprob_branch;
	bool chi2minprob_isLoaded;
	int	nbtagsssv_;
	TBranch *nbtagsssv_branch;
	bool nbtagsssv_isLoaded;
	int	nbtagstcl_;
	TBranch *nbtagstcl_branch;
	bool nbtagstcl_isLoaded;
	int	nbtagstcm_;
	TBranch *nbtagstcm_branch;
	bool nbtagstcm_isLoaded;
	int	nbtagscsvl_;
	TBranch *nbtagscsvl_branch;
	bool nbtagscsvl_isLoaded;
	int	nbtagscsvm_;
	TBranch *nbtagscsvm_branch;
	bool nbtagscsvm_isLoaded;
	int	nbtagscsvt_;
	TBranch *nbtagscsvt_branch;
	bool nbtagscsvt_isLoaded;
	int	nbtagsssvcorr_;
	TBranch *nbtagsssvcorr_branch;
	bool nbtagsssvcorr_isLoaded;
	int	nbtagstclcorr_;
	TBranch *nbtagstclcorr_branch;
	bool nbtagstclcorr_isLoaded;
	int	nbtagstcmcorr_;
	TBranch *nbtagstcmcorr_branch;
	bool nbtagstcmcorr_isLoaded;
	int	nbtagscsvlcorr_;
	TBranch *nbtagscsvlcorr_branch;
	bool nbtagscsvlcorr_isLoaded;
	int	nbtagscsvmcorr_;
	TBranch *nbtagscsvmcorr_branch;
	bool nbtagscsvmcorr_isLoaded;
	int	nbtagscsvtcott_;
	TBranch *nbtagscsvtcott_branch;
	bool nbtagscsvtcott_isLoaded;
	int	njetsUp_;
	TBranch *njetsUp_branch;
	bool njetsUp_isLoaded;
	int	njetsDown_;
	TBranch *njetsDown_branch;
	bool njetsDown_isLoaded;
	float	htUp_;
	TBranch *htUp_branch;
	bool htUp_isLoaded;
	float	htDown_;
	TBranch *htDown_branch;
	bool htDown_isLoaded;
	int	ntruepu_;
	TBranch *ntruepu_branch;
	bool ntruepu_isLoaded;
	int	npu_;
	TBranch *npu_branch;
	bool npu_isLoaded;
	int	npuMinusOne_;
	TBranch *npuMinusOne_branch;
	bool npuMinusOne_isLoaded;
	int	npuPlusOne_;
	TBranch *npuPlusOne_branch;
	bool npuPlusOne_isLoaded;
	int	nvtx_;
	TBranch *nvtx_branch;
	bool nvtx_isLoaded;
	int	indexfirstGoodVertex__;
	TBranch *indexfirstGoodVertex__branch;
	bool indexfirstGoodVertex__isLoaded;
	float	nvtxweight_;
	TBranch *nvtxweight_branch;
	bool nvtxweight_isLoaded;
	float	n3dvtxweight_;
	TBranch *n3dvtxweight_branch;
	bool n3dvtxweight_isLoaded;
	int	pdfid1_;
	TBranch *pdfid1_branch;
	bool pdfid1_isLoaded;
	int	pdfid2_;
	TBranch *pdfid2_branch;
	bool pdfid2_isLoaded;
	float	pdfx1_;
	TBranch *pdfx1_branch;
	bool pdfx1_isLoaded;
	float	pdfx2_;
	TBranch *pdfx2_branch;
	bool pdfx2_isLoaded;
	float	pdfQ_;
	TBranch *pdfQ_branch;
	bool pdfQ_isLoaded;
	float	vecjetpt_;
	TBranch *vecjetpt_branch;
	bool vecjetpt_isLoaded;
	int	pass_;
	TBranch *pass_branch;
	bool pass_isLoaded;
	int	passz_;
	TBranch *passz_branch;
	bool passz_isLoaded;
	float	m0_;
	TBranch *m0_branch;
	bool m0_isLoaded;
	float	mg_;
	TBranch *mg_branch;
	bool mg_isLoaded;
	float	ml_;
	TBranch *ml_branch;
	bool ml_isLoaded;
	float	x_;
	TBranch *x_branch;
	bool x_isLoaded;
	float	m12_;
	TBranch *m12_branch;
	bool m12_isLoaded;
	float	lep1chi2ndf_;
	TBranch *lep1chi2ndf_branch;
	bool lep1chi2ndf_isLoaded;
	float	lep2chi2ndf_;
	TBranch *lep2chi2ndf_branch;
	bool lep2chi2ndf_isLoaded;
	float	lep1dpt_;
	TBranch *lep1dpt_branch;
	bool lep1dpt_isLoaded;
	float	lep2dpt_;
	TBranch *lep2dpt_branch;
	bool lep2dpt_isLoaded;
	int	id1_;
	TBranch *id1_branch;
	bool id1_isLoaded;
	int	id2_;
	TBranch *id2_branch;
	bool id2_isLoaded;
	int	leptype1_;
	TBranch *leptype1_branch;
	bool leptype1_isLoaded;
	int	leptype2_;
	TBranch *leptype2_branch;
	bool leptype2_isLoaded;
	int	w1_;
	TBranch *w1_branch;
	bool w1_isLoaded;
	int	w2_;
	TBranch *w2_branch;
	bool w2_isLoaded;
	float	iso1_;
	TBranch *iso1_branch;
	bool iso1_isLoaded;
	float	isont1_;
	TBranch *isont1_branch;
	bool isont1_isLoaded;
	float	isopfold1_;
	TBranch *isopfold1_branch;
	bool isopfold1_isLoaded;
	float	isopf1_;
	TBranch *isopf1_branch;
	bool isopf1_isLoaded;
	float	etasc1_;
	TBranch *etasc1_branch;
	bool etasc1_isLoaded;
	float	etasc2_;
	TBranch *etasc2_branch;
	bool etasc2_isLoaded;
	float	eoverpin_;
	TBranch *eoverpin_branch;
	bool eoverpin_isLoaded;
	float	eoverpout_;
	TBranch *eoverpout_branch;
	bool eoverpout_isLoaded;
	float	dEtaIn_;
	TBranch *dEtaIn_branch;
	bool dEtaIn_isLoaded;
	float	dPhiIn_;
	TBranch *dPhiIn_branch;
	bool dPhiIn_isLoaded;
	float	sigmaIEtaIEta_;
	TBranch *sigmaIEtaIEta_branch;
	bool sigmaIEtaIEta_isLoaded;
	float	hOverE_;
	TBranch *hOverE_branch;
	bool hOverE_isLoaded;
	float	ooemoop_;
	TBranch *ooemoop_branch;
	bool ooemoop_isLoaded;
	float	d0vtx_;
	TBranch *d0vtx_branch;
	bool d0vtx_isLoaded;
	float	dzvtx_;
	TBranch *dzvtx_branch;
	bool dzvtx_isLoaded;
	float	expinnerlayers_;
	TBranch *expinnerlayers_branch;
	bool expinnerlayers_isLoaded;
	float	fbrem_;
	TBranch *fbrem_branch;
	bool fbrem_isLoaded;
	float	pfisoch_;
	TBranch *pfisoch_branch;
	bool pfisoch_isLoaded;
	float	pfisoem_;
	TBranch *pfisoem_branch;
	bool pfisoem_isLoaded;
	float	pfisonh_;
	TBranch *pfisonh_branch;
	bool pfisonh_isLoaded;
	float	eSC_;
	TBranch *eSC_branch;
	bool eSC_isLoaded;
	float	phiSC_;
	TBranch *phiSC_branch;
	bool phiSC_isLoaded;
	float	eSCRaw_;
	TBranch *eSCRaw_branch;
	bool eSCRaw_isLoaded;
	float	eSCPresh_;
	TBranch *eSCPresh_branch;
	bool eSCPresh_isLoaded;
	float	lep1_scslasercormean_;
	TBranch *lep1_scslasercormean_branch;
	bool lep1_scslasercormean_isLoaded;
	float	lep1_scslasercormax_;
	TBranch *lep1_scslasercormax_branch;
	bool lep1_scslasercormax_isLoaded;
	float	eoverpin2_;
	TBranch *eoverpin2_branch;
	bool eoverpin2_isLoaded;
	float	eoverpout2_;
	TBranch *eoverpout2_branch;
	bool eoverpout2_isLoaded;
	float	dEtaIn2_;
	TBranch *dEtaIn2_branch;
	bool dEtaIn2_isLoaded;
	float	dPhiIn2_;
	TBranch *dPhiIn2_branch;
	bool dPhiIn2_isLoaded;
	float	sigmaIEtaIEta2_;
	TBranch *sigmaIEtaIEta2_branch;
	bool sigmaIEtaIEta2_isLoaded;
	float	hOverE2_;
	TBranch *hOverE2_branch;
	bool hOverE2_isLoaded;
	float	ooemoop2_;
	TBranch *ooemoop2_branch;
	bool ooemoop2_isLoaded;
	float	d0vtx2_;
	TBranch *d0vtx2_branch;
	bool d0vtx2_isLoaded;
	float	dzvtx2_;
	TBranch *dzvtx2_branch;
	bool dzvtx2_isLoaded;
	float	expinnerlayers2_;
	TBranch *expinnerlayers2_branch;
	bool expinnerlayers2_isLoaded;
	float	fbrem2_;
	TBranch *fbrem2_branch;
	bool fbrem2_isLoaded;
	float	pfisoch2_;
	TBranch *pfisoch2_branch;
	bool pfisoch2_isLoaded;
	float	pfisoem2_;
	TBranch *pfisoem2_branch;
	bool pfisoem2_isLoaded;
	float	pfisonh2_;
	TBranch *pfisonh2_branch;
	bool pfisonh2_isLoaded;
	float	eSC2_;
	TBranch *eSC2_branch;
	bool eSC2_isLoaded;
	float	phiSC2_;
	TBranch *phiSC2_branch;
	bool phiSC2_isLoaded;
	float	eSCRaw2_;
	TBranch *eSCRaw2_branch;
	bool eSCRaw2_isLoaded;
	float	eSCPresh2_;
	TBranch *eSCPresh2_branch;
	bool eSCPresh2_isLoaded;
	float	lep2_scslasercormean_;
	TBranch *lep2_scslasercormean_branch;
	bool lep2_scslasercormean_isLoaded;
	float	lep2_scslasercormax_;
	TBranch *lep2_scslasercormax_branch;
	bool lep2_scslasercormax_isLoaded;
	float	scslasercormax_;
	TBranch *scslasercormax_branch;
	bool scslasercormax_isLoaded;
	float	scslasercormax_pt_;
	TBranch *scslasercormax_pt_branch;
	bool scslasercormax_pt_isLoaded;
	float	scslasercormax_eta_;
	TBranch *scslasercormax_eta_branch;
	bool scslasercormax_eta_isLoaded;
	float	iso2_;
	TBranch *iso2_branch;
	bool iso2_isLoaded;
	float	ecalveto1_;
	TBranch *ecalveto1_branch;
	bool ecalveto1_isLoaded;
	float	ecalveto2_;
	TBranch *ecalveto2_branch;
	bool ecalveto2_isLoaded;
	float	hcalveto1_;
	TBranch *hcalveto1_branch;
	bool hcalveto1_isLoaded;
	float	hcalveto2_;
	TBranch *hcalveto2_branch;
	bool hcalveto2_isLoaded;
	float	isont2_;
	TBranch *isont2_branch;
	bool isont2_isLoaded;
	float	isopf2_;
	TBranch *isopf2_branch;
	bool isopf2_isLoaded;
	float	ptl1_;
	TBranch *ptl1_branch;
	bool ptl1_isLoaded;
	float	ptl2_;
	TBranch *ptl2_branch;
	bool ptl2_isLoaded;
	float	etal1_;
	TBranch *etal1_branch;
	bool etal1_isLoaded;
	float	etal2_;
	TBranch *etal2_branch;
	bool etal2_isLoaded;
	float	phil1_;
	TBranch *phil1_branch;
	bool phil1_isLoaded;
	float	phil2_;
	TBranch *phil2_branch;
	bool phil2_isLoaded;
	float	meff_;
	TBranch *meff_branch;
	bool meff_isLoaded;
	float	mt_;
	TBranch *mt_branch;
	bool mt_isLoaded;
	unsigned int	run_;
	TBranch *run_branch;
	bool run_isLoaded;
	unsigned int	lumi_;
	TBranch *lumi_branch;
	bool lumi_isLoaded;
	unsigned int	event_;
	TBranch *event_branch;
	bool event_isLoaded;
	float	y_;
	TBranch *y_branch;
	bool y_isLoaded;
	float	ht_;
	TBranch *ht_branch;
	bool ht_isLoaded;
	float	htgen_;
	TBranch *htgen_branch;
	bool htgen_isLoaded;
	float	htjpt_;
	TBranch *htjpt_branch;
	bool htjpt_isLoaded;
	int	nels_;
	TBranch *nels_branch;
	bool nels_isLoaded;
	int	nmus_;
	TBranch *nmus_branch;
	bool nmus_isLoaded;
	int	ntaus_;
	TBranch *ntaus_branch;
	bool ntaus_isLoaded;
	int	nleps_;
	TBranch *nleps_branch;
	bool nleps_isLoaded;
	int	nbs_;
	TBranch *nbs_branch;
	bool nbs_isLoaded;
	float	dphijm_;
	TBranch *dphijm_branch;
	bool dphijm_isLoaded;
	float	ptjetraw_;
	TBranch *ptjetraw_branch;
	bool ptjetraw_isLoaded;
	float	ptjet23_;
	TBranch *ptjet23_branch;
	bool ptjet23_isLoaded;
	float	ptjetF23_;
	TBranch *ptjetF23_branch;
	bool ptjetF23_isLoaded;
	float	ptjetO23_;
	TBranch *ptjetO23_branch;
	bool ptjetO23_isLoaded;
	int	mcid1_;
	TBranch *mcid1_branch;
	bool mcid1_isLoaded;
	float	mcdr1_;
	TBranch *mcdr1_branch;
	bool mcdr1_isLoaded;
	int	mcdecay1_;
	TBranch *mcdecay1_branch;
	bool mcdecay1_isLoaded;
	int	mcndec1_;
	TBranch *mcndec1_branch;
	bool mcndec1_isLoaded;
	int	mcndec2_;
	TBranch *mcndec2_branch;
	bool mcndec2_isLoaded;
	int	mcndeckls1_;
	TBranch *mcndeckls1_branch;
	bool mcndeckls1_isLoaded;
	int	mcndeckls2_;
	TBranch *mcndeckls2_branch;
	bool mcndeckls2_isLoaded;
	int	mcndecem1_;
	TBranch *mcndecem1_branch;
	bool mcndecem1_isLoaded;
	int	mcndecem2_;
	TBranch *mcndecem2_branch;
	bool mcndecem2_isLoaded;
	int	mcid2_;
	TBranch *mcid2_branch;
	bool mcid2_isLoaded;
	float	mcdr2_;
	TBranch *mcdr2_branch;
	bool mcdr2_isLoaded;
	int	mcdecay2_;
	TBranch *mcdecay2_branch;
	bool mcdecay2_isLoaded;
	float	mctaudpt1_;
	TBranch *mctaudpt1_branch;
	bool mctaudpt1_isLoaded;
	float	mctaudpt2_;
	TBranch *mctaudpt2_branch;
	bool mctaudpt2_isLoaded;
	int	mctaudid1_;
	TBranch *mctaudid1_branch;
	bool mctaudid1_isLoaded;
	int	mctaudid2_;
	TBranch *mctaudid2_branch;
	bool mctaudid2_isLoaded;
	int	mlepid_;
	TBranch *mlepid_branch;
	bool mlepid_isLoaded;
	int	mleppassid_;
	TBranch *mleppassid_branch;
	bool mleppassid_isLoaded;
	int	mleppassiso_;
	TBranch *mleppassiso_branch;
	bool mleppassiso_isLoaded;
	float	mlepiso_;
	TBranch *mlepiso_branch;
	bool mlepiso_isLoaded;
	float	mlepdr_;
	TBranch *mlepdr_branch;
	bool mlepdr_isLoaded;
	float	pflepiso_;
	TBranch *pflepiso_branch;
	bool pflepiso_isLoaded;
	float	pflepdr_;
	TBranch *pflepdr_branch;
	bool pflepdr_isLoaded;
	float	pfleppt_;
	TBranch *pfleppt_branch;
	bool pfleppt_isLoaded;
	float	pflepmindrj_;
	TBranch *pflepmindrj_branch;
	bool pflepmindrj_isLoaded;
	float	pftaudiso_;
	TBranch *pftaudiso_branch;
	bool pftaudiso_isLoaded;
	float	pftauddr_;
	TBranch *pftauddr_branch;
	bool pftauddr_isLoaded;
	float	pftaudpt_;
	TBranch *pftaudpt_branch;
	bool pftaudpt_isLoaded;
	float	pftaudmindrj_;
	TBranch *pftaudmindrj_branch;
	bool pftaudmindrj_isLoaded;
	int	pfcandid5_;
	TBranch *pfcandid5_branch;
	bool pfcandid5_isLoaded;
	float	pfcandiso5_;
	TBranch *pfcandiso5_branch;
	bool pfcandiso5_isLoaded;
	float	pfcandpt5_;
	TBranch *pfcandpt5_branch;
	bool pfcandpt5_isLoaded;
	float	pfcanddz5_;
	TBranch *pfcanddz5_branch;
	bool pfcanddz5_isLoaded;
	float	pfcandmindrj5_;
	TBranch *pfcandmindrj5_branch;
	bool pfcandmindrj5_isLoaded;
	int	pfcandid10_;
	TBranch *pfcandid10_branch;
	bool pfcandid10_isLoaded;
	float	pfcandiso10_;
	TBranch *pfcandiso10_branch;
	bool pfcandiso10_isLoaded;
	float	pfcandpt10_;
	TBranch *pfcandpt10_branch;
	bool pfcandpt10_isLoaded;
	float	pfcanddz10_;
	TBranch *pfcanddz10_branch;
	bool pfcanddz10_isLoaded;
	float	pfcandmindrj10_;
	TBranch *pfcandmindrj10_branch;
	bool pfcandmindrj10_isLoaded;
	int	pfcandidOS10_;
	TBranch *pfcandidOS10_branch;
	bool pfcandidOS10_isLoaded;
	float	pfcandisoOS10_;
	TBranch *pfcandisoOS10_branch;
	bool pfcandisoOS10_isLoaded;
	float	pfcandptOS10_;
	TBranch *pfcandptOS10_branch;
	bool pfcandptOS10_isLoaded;
	float	pfcanddzOS10_;
	TBranch *pfcanddzOS10_branch;
	bool pfcanddzOS10_isLoaded;
	int	pfcandid5looseZ_;
	TBranch *pfcandid5looseZ_branch;
	bool pfcandid5looseZ_isLoaded;
	float	pfcandiso5looseZ_;
	TBranch *pfcandiso5looseZ_branch;
	bool pfcandiso5looseZ_isLoaded;
	float	pfcandpt5looseZ_;
	TBranch *pfcandpt5looseZ_branch;
	bool pfcandpt5looseZ_isLoaded;
	float	pfcanddz5looseZ_;
	TBranch *pfcanddz5looseZ_branch;
	bool pfcanddz5looseZ_isLoaded;
	int	pfcandidOS10looseZ_;
	TBranch *pfcandidOS10looseZ_branch;
	bool pfcandidOS10looseZ_isLoaded;
	float	pfcandisoOS10looseZ_;
	TBranch *pfcandisoOS10looseZ_branch;
	bool pfcandisoOS10looseZ_isLoaded;
	float	pfcandptOS10looseZ_;
	TBranch *pfcandptOS10looseZ_branch;
	bool pfcandptOS10looseZ_isLoaded;
	float	pfcanddzOS10looseZ_;
	TBranch *pfcanddzOS10looseZ_branch;
	bool pfcanddzOS10looseZ_isLoaded;
	int	pfcanddirid10_;
	TBranch *pfcanddirid10_branch;
	bool pfcanddirid10_isLoaded;
	float	pfcanddiriso10_;
	TBranch *pfcanddiriso10_branch;
	bool pfcanddiriso10_isLoaded;
	float	pfcanddirpt10_;
	TBranch *pfcanddirpt10_branch;
	bool pfcanddirpt10_isLoaded;
	float	pfcanddirmindrj10_;
	TBranch *pfcanddirmindrj10_branch;
	bool pfcanddirmindrj10_isLoaded;
	int	pfcandvetoid10_;
	TBranch *pfcandvetoid10_branch;
	bool pfcandvetoid10_isLoaded;
	float	pfcandvetoiso10_;
	TBranch *pfcandvetoiso10_branch;
	bool pfcandvetoiso10_isLoaded;
	float	pfcandvetopt10_;
	TBranch *pfcandvetopt10_branch;
	bool pfcandvetopt10_isLoaded;
	float	pfcandvetomindrj10_;
	TBranch *pfcandvetomindrj10_branch;
	bool pfcandvetomindrj10_isLoaded;
	int	pfcandvetoLid10_;
	TBranch *pfcandvetoLid10_branch;
	bool pfcandvetoLid10_isLoaded;
	float	pfcandvetoLiso10_;
	TBranch *pfcandvetoLiso10_branch;
	bool pfcandvetoLiso10_isLoaded;
	float	pfcandvetoLpt10_;
	TBranch *pfcandvetoLpt10_branch;
	bool pfcandvetoLpt10_isLoaded;
	float	pfcandvetoLmindrj10_;
	TBranch *pfcandvetoLmindrj10_branch;
	bool pfcandvetoLmindrj10_isLoaded;
	float	emjet10_;
	TBranch *emjet10_branch;
	bool emjet10_isLoaded;
	float	mjj_;
	TBranch *mjj_branch;
	bool mjj_isLoaded;
	float	emjet20_;
	TBranch *emjet20_branch;
	bool emjet20_isLoaded;
	float	trkpt5_;
	TBranch *trkpt5_branch;
	bool trkpt5_isLoaded;
	float	trkpt10_;
	TBranch *trkpt10_branch;
	bool trkpt10_isLoaded;
	float	mleptrk5_;
	TBranch *mleptrk5_branch;
	bool mleptrk5_isLoaded;
	float	mleptrk10_;
	TBranch *mleptrk10_branch;
	bool mleptrk10_isLoaded;
	float	trkreliso5_;
	TBranch *trkreliso5_branch;
	bool trkreliso5_isLoaded;
	float	trkreliso10_;
	TBranch *trkreliso10_branch;
	bool trkreliso10_isLoaded;
	float	trkpt5loose_;
	TBranch *trkpt5loose_branch;
	bool trkpt5loose_isLoaded;
	float	trkpt10loose_;
	TBranch *trkpt10loose_branch;
	bool trkpt10loose_isLoaded;
	float	trkreliso5loose_;
	TBranch *trkreliso5loose_branch;
	bool trkreliso5loose_isLoaded;
	float	trkreliso10loose_;
	TBranch *trkreliso10loose_branch;
	bool trkreliso10loose_isLoaded;
	float	trkpt10pt0p1_;
	TBranch *trkpt10pt0p1_branch;
	bool trkpt10pt0p1_isLoaded;
	float	trkpt10pt0p2_;
	TBranch *trkpt10pt0p2_branch;
	bool trkpt10pt0p2_isLoaded;
	float	trkpt10pt0p3_;
	TBranch *trkpt10pt0p3_branch;
	bool trkpt10pt0p3_isLoaded;
	float	trkpt10pt0p4_;
	TBranch *trkpt10pt0p4_branch;
	bool trkpt10pt0p4_isLoaded;
	float	trkpt10pt0p5_;
	TBranch *trkpt10pt0p5_branch;
	bool trkpt10pt0p5_isLoaded;
	float	trkpt10pt0p6_;
	TBranch *trkpt10pt0p6_branch;
	bool trkpt10pt0p6_isLoaded;
	float	trkpt10pt0p7_;
	TBranch *trkpt10pt0p7_branch;
	bool trkpt10pt0p7_isLoaded;
	float	trkpt10pt0p8_;
	TBranch *trkpt10pt0p8_branch;
	bool trkpt10pt0p8_isLoaded;
	float	trkpt10pt0p9_;
	TBranch *trkpt10pt0p9_branch;
	bool trkpt10pt0p9_isLoaded;
	float	trkpt10pt1p0_;
	TBranch *trkpt10pt1p0_branch;
	bool trkpt10pt1p0_isLoaded;
	float	trkreliso10pt0p1_;
	TBranch *trkreliso10pt0p1_branch;
	bool trkreliso10pt0p1_isLoaded;
	float	trkreliso10pt0p2_;
	TBranch *trkreliso10pt0p2_branch;
	bool trkreliso10pt0p2_isLoaded;
	float	trkreliso10pt0p3_;
	TBranch *trkreliso10pt0p3_branch;
	bool trkreliso10pt0p3_isLoaded;
	float	trkreliso10pt0p4_;
	TBranch *trkreliso10pt0p4_branch;
	bool trkreliso10pt0p4_isLoaded;
	float	trkreliso10pt0p5_;
	TBranch *trkreliso10pt0p5_branch;
	bool trkreliso10pt0p5_isLoaded;
	float	trkreliso10pt0p6_;
	TBranch *trkreliso10pt0p6_branch;
	bool trkreliso10pt0p6_isLoaded;
	float	trkreliso10pt0p7_;
	TBranch *trkreliso10pt0p7_branch;
	bool trkreliso10pt0p7_isLoaded;
	float	trkreliso10pt0p8_;
	TBranch *trkreliso10pt0p8_branch;
	bool trkreliso10pt0p8_isLoaded;
	float	trkreliso10pt0p9_;
	TBranch *trkreliso10pt0p9_branch;
	bool trkreliso10pt0p9_isLoaded;
	float	trkreliso10pt1p0_;
	TBranch *trkreliso10pt1p0_branch;
	bool trkreliso10pt1p0_isLoaded;
	float	pfcandpt10pt0p1_;
	TBranch *pfcandpt10pt0p1_branch;
	bool pfcandpt10pt0p1_isLoaded;
	float	pfcandpt10pt0p2_;
	TBranch *pfcandpt10pt0p2_branch;
	bool pfcandpt10pt0p2_isLoaded;
	float	pfcandpt10pt0p3_;
	TBranch *pfcandpt10pt0p3_branch;
	bool pfcandpt10pt0p3_isLoaded;
	float	pfcandpt10pt0p4_;
	TBranch *pfcandpt10pt0p4_branch;
	bool pfcandpt10pt0p4_isLoaded;
	float	pfcandpt10pt0p5_;
	TBranch *pfcandpt10pt0p5_branch;
	bool pfcandpt10pt0p5_isLoaded;
	float	pfcandpt10pt0p6_;
	TBranch *pfcandpt10pt0p6_branch;
	bool pfcandpt10pt0p6_isLoaded;
	float	pfcandpt10pt0p7_;
	TBranch *pfcandpt10pt0p7_branch;
	bool pfcandpt10pt0p7_isLoaded;
	float	pfcandpt10pt0p8_;
	TBranch *pfcandpt10pt0p8_branch;
	bool pfcandpt10pt0p8_isLoaded;
	float	pfcandpt10pt0p9_;
	TBranch *pfcandpt10pt0p9_branch;
	bool pfcandpt10pt0p9_isLoaded;
	float	pfcandpt10pt1p0_;
	TBranch *pfcandpt10pt1p0_branch;
	bool pfcandpt10pt1p0_isLoaded;
	float	pfcandiso10pt0p1_;
	TBranch *pfcandiso10pt0p1_branch;
	bool pfcandiso10pt0p1_isLoaded;
	float	pfcandiso10pt0p2_;
	TBranch *pfcandiso10pt0p2_branch;
	bool pfcandiso10pt0p2_isLoaded;
	float	pfcandiso10pt0p3_;
	TBranch *pfcandiso10pt0p3_branch;
	bool pfcandiso10pt0p3_isLoaded;
	float	pfcandiso10pt0p4_;
	TBranch *pfcandiso10pt0p4_branch;
	bool pfcandiso10pt0p4_isLoaded;
	float	pfcandiso10pt0p5_;
	TBranch *pfcandiso10pt0p5_branch;
	bool pfcandiso10pt0p5_isLoaded;
	float	pfcandiso10pt0p6_;
	TBranch *pfcandiso10pt0p6_branch;
	bool pfcandiso10pt0p6_isLoaded;
	float	pfcandiso10pt0p7_;
	TBranch *pfcandiso10pt0p7_branch;
	bool pfcandiso10pt0p7_isLoaded;
	float	pfcandiso10pt0p8_;
	TBranch *pfcandiso10pt0p8_branch;
	bool pfcandiso10pt0p8_isLoaded;
	float	pfcandiso10pt0p9_;
	TBranch *pfcandiso10pt0p9_branch;
	bool pfcandiso10pt0p9_isLoaded;
	float	pfcandiso10pt1p0_;
	TBranch *pfcandiso10pt1p0_branch;
	bool pfcandiso10pt1p0_isLoaded;
	float	mbb_;
	TBranch *mbb_branch;
	bool mbb_isLoaded;
	float	lep1pfjetdr_;
	TBranch *lep1pfjetdr_branch;
	bool lep1pfjetdr_isLoaded;
	float	lep2pfjetdr_;
	TBranch *lep2pfjetdr_branch;
	bool lep2pfjetdr_isLoaded;
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > *mclep_;
	TBranch *mclep_branch;
	bool mclep_isLoaded;
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > *mcnu_;
	TBranch *mcnu_branch;
	bool mcnu_isLoaded;
	float	mcmln_;
	TBranch *mcmln_branch;
	bool mcmln_isLoaded;
	float	mcmtln_;
	TBranch *mcmtln_branch;
	bool mcmtln_isLoaded;
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > *mlep_;
	TBranch *mlep_branch;
	bool mlep_isLoaded;
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > *lep1_;
	TBranch *lep1_branch;
	bool lep1_isLoaded;
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > *lep2_;
	TBranch *lep2_branch;
	bool lep2_isLoaded;
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > *trklep1_;
	TBranch *trklep1_branch;
	bool trklep1_isLoaded;
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > *trklep2_;
	TBranch *trklep2_branch;
	bool trklep2_isLoaded;
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > *gfitlep1_;
	TBranch *gfitlep1_branch;
	bool gfitlep1_isLoaded;
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > *gfitlep2_;
	TBranch *gfitlep2_branch;
	bool gfitlep2_isLoaded;
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > *lepp_;
	TBranch *lepp_branch;
	bool lepp_isLoaded;
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > *lepm_;
	TBranch *lepm_branch;
	bool lepm_isLoaded;
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > *pflep1_;
	TBranch *pflep1_branch;
	bool pflep1_isLoaded;
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > *pflep2_;
	TBranch *pflep2_branch;
	bool pflep2_isLoaded;
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > *leppfjet1_;
	TBranch *leppfjet1_branch;
	bool leppfjet1_isLoaded;
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > *leppfjet2_;
	TBranch *leppfjet2_branch;
	bool leppfjet2_isLoaded;
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > *mclep1_;
	TBranch *mclep1_branch;
	bool mclep1_isLoaded;
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > *mclep2_;
	TBranch *mclep2_branch;
	bool mclep2_isLoaded;
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > *mctaud1_;
	TBranch *mctaud1_branch;
	bool mctaud1_isLoaded;
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > *mctaud2_;
	TBranch *mctaud2_branch;
	bool mctaud2_isLoaded;
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > *mctaudvis1_;
	TBranch *mctaudvis1_branch;
	bool mctaudvis1_isLoaded;
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > *mctaudvis2_;
	TBranch *mctaudvis2_branch;
	bool mctaudvis2_isLoaded;
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > *pflep_;
	TBranch *pflep_branch;
	bool pflep_isLoaded;
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > *pftaud_;
	TBranch *pftaud_branch;
	bool pftaud_isLoaded;
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > *pfcand5_;
	TBranch *pfcand5_branch;
	bool pfcand5_isLoaded;
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > *pfcand10_;
	TBranch *pfcand10_branch;
	bool pfcand10_isLoaded;
	int	pfTau15_leadPtcandID_;
	TBranch *pfTau15_leadPtcandID_branch;
	bool pfTau15_leadPtcandID_isLoaded;
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > *pfTau15_;
	TBranch *pfTau15_branch;
	bool pfTau15_isLoaded;
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > *pfTau15_leadPtcand_;
	TBranch *pfTau15_leadPtcand_branch;
	bool pfTau15_leadPtcand_isLoaded;
	int	pfTau_leadPtcandID_;
	TBranch *pfTau_leadPtcandID_branch;
	bool pfTau_leadPtcandID_isLoaded;
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > *pfTau_;
	TBranch *pfTau_branch;
	bool pfTau_isLoaded;
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > *pfTau_leadPtcand_;
	TBranch *pfTau_leadPtcand_branch;
	bool pfTau_leadPtcand_isLoaded;
	int	pfTauLoose_leadPtcandID_;
	TBranch *pfTauLoose_leadPtcandID_branch;
	bool pfTauLoose_leadPtcandID_isLoaded;
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > *pfTauLoose_;
	TBranch *pfTauLoose_branch;
	bool pfTauLoose_isLoaded;
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > *pfTauLoose_leadPtcand_;
	TBranch *pfTauLoose_leadPtcand_branch;
	bool pfTauLoose_leadPtcand_isLoaded;
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > *pfcandOS10_;
	TBranch *pfcandOS10_branch;
	bool pfcandOS10_isLoaded;
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > *pfcandOS10looseZ_;
	TBranch *pfcandOS10looseZ_branch;
	bool pfcandOS10looseZ_isLoaded;
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > *pfcand5looseZ_;
	TBranch *pfcand5looseZ_branch;
	bool pfcand5looseZ_isLoaded;
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > *pfcanddir10_;
	TBranch *pfcanddir10_branch;
	bool pfcanddir10_isLoaded;
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > *pfcandveto10_;
	TBranch *pfcandveto10_branch;
	bool pfcandveto10_isLoaded;
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > *pfcandvetoL10_;
	TBranch *pfcandvetoL10_branch;
	bool pfcandvetoL10_isLoaded;
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > *jet_;
	TBranch *jet_branch;
	bool jet_isLoaded;
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > *nonisoel_;
	TBranch *nonisoel_branch;
	bool nonisoel_isLoaded;
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > *nonisomu_;
	TBranch *nonisomu_branch;
	bool nonisomu_isLoaded;
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > *t_;
	TBranch *t_branch;
	bool t_isLoaded;
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > *tbar_;
	TBranch *tbar_branch;
	bool tbar_isLoaded;
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > *ttbar_;
	TBranch *ttbar_branch;
	bool ttbar_isLoaded;
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > *lep_t_;
	TBranch *lep_t_branch;
	bool lep_t_isLoaded;
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > *lep_tbar_;
	TBranch *lep_tbar_branch;
	bool lep_tbar_isLoaded;
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > *stop_t_;
	TBranch *stop_t_branch;
	bool stop_t_isLoaded;
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > *stop_tbar_;
	TBranch *stop_tbar_branch;
	bool stop_tbar_isLoaded;
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > *neutralino_t_;
	TBranch *neutralino_t_branch;
	bool neutralino_t_isLoaded;
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > *neutralino_tbar_;
	TBranch *neutralino_tbar_branch;
	bool neutralino_tbar_isLoaded;
	int	lep_t_id_;
	TBranch *lep_t_id_branch;
	bool lep_t_id_isLoaded;
	int	lep_tbar_id_;
	TBranch *lep_tbar_id_branch;
	bool lep_tbar_id_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > *pfjets_;
	TBranch *pfjets_branch;
	bool pfjets_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > *pfjets_genJet__;
	TBranch *pfjets_genJet__branch;
	bool pfjets_genJet__isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > *pfjets_failjetid_;
	TBranch *pfjets_failjetid_branch;
	bool pfjets_failjetid_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > *pfjets_faillepolap_;
	TBranch *pfjets_faillepolap_branch;
	bool pfjets_faillepolap_isLoaded;
	vector<float> *pfjets_csv_;
	TBranch *pfjets_csv_branch;
	bool pfjets_csv_isLoaded;
	vector<float> *pfjets_chEfrac_;
	TBranch *pfjets_chEfrac_branch;
	bool pfjets_chEfrac_isLoaded;
	vector<float> *pfjets_chm_;
	TBranch *pfjets_chm_branch;
	bool pfjets_chm_isLoaded;
	vector<float> *pfjets_neu_;
	TBranch *pfjets_neu_branch;
	bool pfjets_neu_isLoaded;
	vector<float> *pfjets_l1corr_;
	TBranch *pfjets_l1corr_branch;
	bool pfjets_l1corr_isLoaded;
	vector<float> *pfjets_corr_;
	TBranch *pfjets_corr_branch;
	bool pfjets_corr_isLoaded;
	vector<int> *pfjets_mc3_;
	TBranch *pfjets_mc3_branch;
	bool pfjets_mc3_isLoaded;
	vector<int> *pfjets_mcflavorAlgo_;
	TBranch *pfjets_mcflavorAlgo_branch;
	bool pfjets_mcflavorAlgo_isLoaded;
	vector<int> *pfjets_mcflavorPhys_;
	TBranch *pfjets_mcflavorPhys_branch;
	bool pfjets_mcflavorPhys_isLoaded;
	vector<float> *pfjets_uncertainty_;
	TBranch *pfjets_uncertainty_branch;
	bool pfjets_uncertainty_isLoaded;
	vector<int> *pfjets_flav_;
	TBranch *pfjets_flav_branch;
	bool pfjets_flav_isLoaded;
	vector<float> *pfjets_lrm_;
	TBranch *pfjets_lrm_branch;
	bool pfjets_lrm_isLoaded;
	vector<float> *pfjets_lrm2_;
	TBranch *pfjets_lrm2_branch;
	bool pfjets_lrm2_isLoaded;
	vector<float> *pfjets_qgtag_;
	TBranch *pfjets_qgtag_branch;
	bool pfjets_qgtag_isLoaded;
	vector<float> *pfjets_genJetDr_;
	TBranch *pfjets_genJetDr_branch;
	bool pfjets_genJetDr_isLoaded;
	vector<float> *pfjets_sigma_;
	TBranch *pfjets_sigma_branch;
	bool pfjets_sigma_isLoaded;
	vector<int> *pfjets_lepjet_;
	TBranch *pfjets_lepjet_branch;
	bool pfjets_lepjet_isLoaded;
	vector<float> *pfjets_tobtecmult_;
	TBranch *pfjets_tobtecmult_branch;
	bool pfjets_tobtecmult_isLoaded;
	vector<float> *pfjets_tobtecfrac_;
	TBranch *pfjets_tobtecfrac_branch;
	bool pfjets_tobtecfrac_isLoaded;
	vector<float> *pfjets_beta_;
	TBranch *pfjets_beta_branch;
	bool pfjets_beta_isLoaded;
	vector<float> *pfjets_beta2_;
	TBranch *pfjets_beta2_branch;
	bool pfjets_beta2_isLoaded;
	vector<float> *pfjets_beta_0p1_;
	TBranch *pfjets_beta_0p1_branch;
	bool pfjets_beta_0p1_isLoaded;
	vector<float> *pfjets_beta_0p2_;
	TBranch *pfjets_beta_0p2_branch;
	bool pfjets_beta_0p2_isLoaded;
	vector<float> *pfjets_beta2_0p1_;
	TBranch *pfjets_beta2_0p1_branch;
	bool pfjets_beta2_0p1_isLoaded;
	vector<float> *pfjets_beta2_0p5_;
	TBranch *pfjets_beta2_0p5_branch;
	bool pfjets_beta2_0p5_isLoaded;
	vector<float> *pfjets_mvaPUid_;
	TBranch *pfjets_mvaPUid_branch;
	bool pfjets_mvaPUid_isLoaded;
	vector<float> *pfjets_mva5xPUid_;
	TBranch *pfjets_mva5xPUid_branch;
	bool pfjets_mva5xPUid_isLoaded;
	vector<float> *pfjets_mvaBeta_;
	TBranch *pfjets_mvaBeta_branch;
	bool pfjets_mvaBeta_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > *genjets_;
	TBranch *genjets_branch;
	bool genjets_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > *genqgs_;
	TBranch *genqgs_branch;
	bool genqgs_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > *genbs_;
	TBranch *genbs_branch;
	bool genbs_isLoaded;
	vector<int> *genps_pdgId_;
	TBranch *genps_pdgId_branch;
	bool genps_pdgId_isLoaded;
	vector<int> *genps_firstMother_;
	TBranch *genps_firstMother_branch;
	bool genps_firstMother_isLoaded;
	vector<float> *genps_energy_;
	TBranch *genps_energy_branch;
	bool genps_energy_isLoaded;
	vector<float> *genps_pt_;
	TBranch *genps_pt_branch;
	bool genps_pt_isLoaded;
	vector<float> *genps_eta_;
	TBranch *genps_eta_branch;
	bool genps_eta_isLoaded;
	vector<float> *genps_phi_;
	TBranch *genps_phi_branch;
	bool genps_phi_isLoaded;
	vector<float> *genps_mass_;
	TBranch *genps_mass_branch;
	bool genps_mass_isLoaded;
public: 
void Init(TTree *tree) {
	mclep_branch = 0;
	if (tree->GetBranch("mclep") != 0) {
		mclep_branch = tree->GetBranch("mclep");
		if (mclep_branch) {mclep_branch->SetAddress(&mclep_);}
	}
	mcnu_branch = 0;
	if (tree->GetBranch("mcnu") != 0) {
		mcnu_branch = tree->GetBranch("mcnu");
		if (mcnu_branch) {mcnu_branch->SetAddress(&mcnu_);}
	}
	mlep_branch = 0;
	if (tree->GetBranch("mlep") != 0) {
		mlep_branch = tree->GetBranch("mlep");
		if (mlep_branch) {mlep_branch->SetAddress(&mlep_);}
	}
	lep1_branch = 0;
	if (tree->GetBranch("lep1") != 0) {
		lep1_branch = tree->GetBranch("lep1");
		if (lep1_branch) {lep1_branch->SetAddress(&lep1_);}
	}
	lep2_branch = 0;
	if (tree->GetBranch("lep2") != 0) {
		lep2_branch = tree->GetBranch("lep2");
		if (lep2_branch) {lep2_branch->SetAddress(&lep2_);}
	}
	trklep1_branch = 0;
	if (tree->GetBranch("trklep1") != 0) {
		trklep1_branch = tree->GetBranch("trklep1");
		if (trklep1_branch) {trklep1_branch->SetAddress(&trklep1_);}
	}
	trklep2_branch = 0;
	if (tree->GetBranch("trklep2") != 0) {
		trklep2_branch = tree->GetBranch("trklep2");
		if (trklep2_branch) {trklep2_branch->SetAddress(&trklep2_);}
	}
	gfitlep1_branch = 0;
	if (tree->GetBranch("gfitlep1") != 0) {
		gfitlep1_branch = tree->GetBranch("gfitlep1");
		if (gfitlep1_branch) {gfitlep1_branch->SetAddress(&gfitlep1_);}
	}
	gfitlep2_branch = 0;
	if (tree->GetBranch("gfitlep2") != 0) {
		gfitlep2_branch = tree->GetBranch("gfitlep2");
		if (gfitlep2_branch) {gfitlep2_branch->SetAddress(&gfitlep2_);}
	}
	lepp_branch = 0;
	if (tree->GetBranch("lepp") != 0) {
		lepp_branch = tree->GetBranch("lepp");
		if (lepp_branch) {lepp_branch->SetAddress(&lepp_);}
	}
	lepm_branch = 0;
	if (tree->GetBranch("lepm") != 0) {
		lepm_branch = tree->GetBranch("lepm");
		if (lepm_branch) {lepm_branch->SetAddress(&lepm_);}
	}
	pflep1_branch = 0;
	if (tree->GetBranch("pflep1") != 0) {
		pflep1_branch = tree->GetBranch("pflep1");
		if (pflep1_branch) {pflep1_branch->SetAddress(&pflep1_);}
	}
	pflep2_branch = 0;
	if (tree->GetBranch("pflep2") != 0) {
		pflep2_branch = tree->GetBranch("pflep2");
		if (pflep2_branch) {pflep2_branch->SetAddress(&pflep2_);}
	}
	leppfjet1_branch = 0;
	if (tree->GetBranch("leppfjet1") != 0) {
		leppfjet1_branch = tree->GetBranch("leppfjet1");
		if (leppfjet1_branch) {leppfjet1_branch->SetAddress(&leppfjet1_);}
	}
	leppfjet2_branch = 0;
	if (tree->GetBranch("leppfjet2") != 0) {
		leppfjet2_branch = tree->GetBranch("leppfjet2");
		if (leppfjet2_branch) {leppfjet2_branch->SetAddress(&leppfjet2_);}
	}
	mclep1_branch = 0;
	if (tree->GetBranch("mclep1") != 0) {
		mclep1_branch = tree->GetBranch("mclep1");
		if (mclep1_branch) {mclep1_branch->SetAddress(&mclep1_);}
	}
	mclep2_branch = 0;
	if (tree->GetBranch("mclep2") != 0) {
		mclep2_branch = tree->GetBranch("mclep2");
		if (mclep2_branch) {mclep2_branch->SetAddress(&mclep2_);}
	}
	mctaud1_branch = 0;
	if (tree->GetBranch("mctaud1") != 0) {
		mctaud1_branch = tree->GetBranch("mctaud1");
		if (mctaud1_branch) {mctaud1_branch->SetAddress(&mctaud1_);}
	}
	mctaud2_branch = 0;
	if (tree->GetBranch("mctaud2") != 0) {
		mctaud2_branch = tree->GetBranch("mctaud2");
		if (mctaud2_branch) {mctaud2_branch->SetAddress(&mctaud2_);}
	}
	mctaudvis1_branch = 0;
	if (tree->GetBranch("mctaudvis1") != 0) {
		mctaudvis1_branch = tree->GetBranch("mctaudvis1");
		if (mctaudvis1_branch) {mctaudvis1_branch->SetAddress(&mctaudvis1_);}
	}
	mctaudvis2_branch = 0;
	if (tree->GetBranch("mctaudvis2") != 0) {
		mctaudvis2_branch = tree->GetBranch("mctaudvis2");
		if (mctaudvis2_branch) {mctaudvis2_branch->SetAddress(&mctaudvis2_);}
	}
	pflep_branch = 0;
	if (tree->GetBranch("pflep") != 0) {
		pflep_branch = tree->GetBranch("pflep");
		if (pflep_branch) {pflep_branch->SetAddress(&pflep_);}
	}
	pftaud_branch = 0;
	if (tree->GetBranch("pftaud") != 0) {
		pftaud_branch = tree->GetBranch("pftaud");
		if (pftaud_branch) {pftaud_branch->SetAddress(&pftaud_);}
	}
	pfcand5_branch = 0;
	if (tree->GetBranch("pfcand5") != 0) {
		pfcand5_branch = tree->GetBranch("pfcand5");
		if (pfcand5_branch) {pfcand5_branch->SetAddress(&pfcand5_);}
	}
	pfcand10_branch = 0;
	if (tree->GetBranch("pfcand10") != 0) {
		pfcand10_branch = tree->GetBranch("pfcand10");
		if (pfcand10_branch) {pfcand10_branch->SetAddress(&pfcand10_);}
	}
	pfTau15_branch = 0;
	if (tree->GetBranch("pfTau15") != 0) {
		pfTau15_branch = tree->GetBranch("pfTau15");
		if (pfTau15_branch) {pfTau15_branch->SetAddress(&pfTau15_);}
	}
	pfTau15_leadPtcand_branch = 0;
	if (tree->GetBranch("pfTau15_leadPtcand") != 0) {
		pfTau15_leadPtcand_branch = tree->GetBranch("pfTau15_leadPtcand");
		if (pfTau15_leadPtcand_branch) {pfTau15_leadPtcand_branch->SetAddress(&pfTau15_leadPtcand_);}
	}
	pfTau_branch = 0;
	if (tree->GetBranch("pfTau") != 0) {
		pfTau_branch = tree->GetBranch("pfTau");
		if (pfTau_branch) {pfTau_branch->SetAddress(&pfTau_);}
	}
	pfTau_leadPtcand_branch = 0;
	if (tree->GetBranch("pfTau_leadPtcand") != 0) {
		pfTau_leadPtcand_branch = tree->GetBranch("pfTau_leadPtcand");
		if (pfTau_leadPtcand_branch) {pfTau_leadPtcand_branch->SetAddress(&pfTau_leadPtcand_);}
	}
	pfTauLoose_branch = 0;
	if (tree->GetBranch("pfTauLoose") != 0) {
		pfTauLoose_branch = tree->GetBranch("pfTauLoose");
		if (pfTauLoose_branch) {pfTauLoose_branch->SetAddress(&pfTauLoose_);}
	}
	pfTauLoose_leadPtcand_branch = 0;
	if (tree->GetBranch("pfTauLoose_leadPtcand") != 0) {
		pfTauLoose_leadPtcand_branch = tree->GetBranch("pfTauLoose_leadPtcand");
		if (pfTauLoose_leadPtcand_branch) {pfTauLoose_leadPtcand_branch->SetAddress(&pfTauLoose_leadPtcand_);}
	}
	pfcandOS10_branch = 0;
	if (tree->GetBranch("pfcandOS10") != 0) {
		pfcandOS10_branch = tree->GetBranch("pfcandOS10");
		if (pfcandOS10_branch) {pfcandOS10_branch->SetAddress(&pfcandOS10_);}
	}
	pfcandOS10looseZ_branch = 0;
	if (tree->GetBranch("pfcandOS10looseZ") != 0) {
		pfcandOS10looseZ_branch = tree->GetBranch("pfcandOS10looseZ");
		if (pfcandOS10looseZ_branch) {pfcandOS10looseZ_branch->SetAddress(&pfcandOS10looseZ_);}
	}
	pfcand5looseZ_branch = 0;
	if (tree->GetBranch("pfcand5looseZ") != 0) {
		pfcand5looseZ_branch = tree->GetBranch("pfcand5looseZ");
		if (pfcand5looseZ_branch) {pfcand5looseZ_branch->SetAddress(&pfcand5looseZ_);}
	}
	pfcanddir10_branch = 0;
	if (tree->GetBranch("pfcanddir10") != 0) {
		pfcanddir10_branch = tree->GetBranch("pfcanddir10");
		if (pfcanddir10_branch) {pfcanddir10_branch->SetAddress(&pfcanddir10_);}
	}
	pfcandveto10_branch = 0;
	if (tree->GetBranch("pfcandveto10") != 0) {
		pfcandveto10_branch = tree->GetBranch("pfcandveto10");
		if (pfcandveto10_branch) {pfcandveto10_branch->SetAddress(&pfcandveto10_);}
	}
	pfcandvetoL10_branch = 0;
	if (tree->GetBranch("pfcandvetoL10") != 0) {
		pfcandvetoL10_branch = tree->GetBranch("pfcandvetoL10");
		if (pfcandvetoL10_branch) {pfcandvetoL10_branch->SetAddress(&pfcandvetoL10_);}
	}
	jet_branch = 0;
	if (tree->GetBranch("jet") != 0) {
		jet_branch = tree->GetBranch("jet");
		if (jet_branch) {jet_branch->SetAddress(&jet_);}
	}
	nonisoel_branch = 0;
	if (tree->GetBranch("nonisoel") != 0) {
		nonisoel_branch = tree->GetBranch("nonisoel");
		if (nonisoel_branch) {nonisoel_branch->SetAddress(&nonisoel_);}
	}
	nonisomu_branch = 0;
	if (tree->GetBranch("nonisomu") != 0) {
		nonisomu_branch = tree->GetBranch("nonisomu");
		if (nonisomu_branch) {nonisomu_branch->SetAddress(&nonisomu_);}
	}
	t_branch = 0;
	if (tree->GetBranch("t") != 0) {
		t_branch = tree->GetBranch("t");
		if (t_branch) {t_branch->SetAddress(&t_);}
	}
	tbar_branch = 0;
	if (tree->GetBranch("tbar") != 0) {
		tbar_branch = tree->GetBranch("tbar");
		if (tbar_branch) {tbar_branch->SetAddress(&tbar_);}
	}
	ttbar_branch = 0;
	if (tree->GetBranch("ttbar") != 0) {
		ttbar_branch = tree->GetBranch("ttbar");
		if (ttbar_branch) {ttbar_branch->SetAddress(&ttbar_);}
	}
	lep_t_branch = 0;
	if (tree->GetBranch("lep_t") != 0) {
		lep_t_branch = tree->GetBranch("lep_t");
		if (lep_t_branch) {lep_t_branch->SetAddress(&lep_t_);}
	}
	lep_tbar_branch = 0;
	if (tree->GetBranch("lep_tbar") != 0) {
		lep_tbar_branch = tree->GetBranch("lep_tbar");
		if (lep_tbar_branch) {lep_tbar_branch->SetAddress(&lep_tbar_);}
	}
	stop_t_branch = 0;
	if (tree->GetBranch("stop_t") != 0) {
		stop_t_branch = tree->GetBranch("stop_t");
		if (stop_t_branch) {stop_t_branch->SetAddress(&stop_t_);}
	}
	stop_tbar_branch = 0;
	if (tree->GetBranch("stop_tbar") != 0) {
		stop_tbar_branch = tree->GetBranch("stop_tbar");
		if (stop_tbar_branch) {stop_tbar_branch->SetAddress(&stop_tbar_);}
	}
	neutralino_t_branch = 0;
	if (tree->GetBranch("neutralino_t") != 0) {
		neutralino_t_branch = tree->GetBranch("neutralino_t");
		if (neutralino_t_branch) {neutralino_t_branch->SetAddress(&neutralino_t_);}
	}
	neutralino_tbar_branch = 0;
	if (tree->GetBranch("neutralino_tbar") != 0) {
		neutralino_tbar_branch = tree->GetBranch("neutralino_tbar");
		if (neutralino_tbar_branch) {neutralino_tbar_branch->SetAddress(&neutralino_tbar_);}
	}
	pfjets_branch = 0;
	if (tree->GetBranch("pfjets") != 0) {
		pfjets_branch = tree->GetBranch("pfjets");
		if (pfjets_branch) {pfjets_branch->SetAddress(&pfjets_);}
	}
	pfjets_genJet__branch = 0;
	if (tree->GetBranch("pfjets_genJet_") != 0) {
		pfjets_genJet__branch = tree->GetBranch("pfjets_genJet_");
		if (pfjets_genJet__branch) {pfjets_genJet__branch->SetAddress(&pfjets_genJet__);}
	}
	pfjets_failjetid_branch = 0;
	if (tree->GetBranch("pfjets_failjetid") != 0) {
		pfjets_failjetid_branch = tree->GetBranch("pfjets_failjetid");
		if (pfjets_failjetid_branch) {pfjets_failjetid_branch->SetAddress(&pfjets_failjetid_);}
	}
	pfjets_faillepolap_branch = 0;
	if (tree->GetBranch("pfjets_faillepolap") != 0) {
		pfjets_faillepolap_branch = tree->GetBranch("pfjets_faillepolap");
		if (pfjets_faillepolap_branch) {pfjets_faillepolap_branch->SetAddress(&pfjets_faillepolap_);}
	}
	genjets_branch = 0;
	if (tree->GetBranch("genjets") != 0) {
		genjets_branch = tree->GetBranch("genjets");
		if (genjets_branch) {genjets_branch->SetAddress(&genjets_);}
	}
	genqgs_branch = 0;
	if (tree->GetBranch("genqgs") != 0) {
		genqgs_branch = tree->GetBranch("genqgs");
		if (genqgs_branch) {genqgs_branch->SetAddress(&genqgs_);}
	}
	genbs_branch = 0;
	if (tree->GetBranch("genbs") != 0) {
		genbs_branch = tree->GetBranch("genbs");
		if (genbs_branch) {genbs_branch->SetAddress(&genbs_);}
	}
  tree->SetMakeClass(1);
	acc_2010_branch = 0;
	if (tree->GetBranch("acc_2010") != 0) {
		acc_2010_branch = tree->GetBranch("acc_2010");
		if (acc_2010_branch) {acc_2010_branch->SetAddress(&acc_2010_);}
	}
	acc_highmet_branch = 0;
	if (tree->GetBranch("acc_highmet") != 0) {
		acc_highmet_branch = tree->GetBranch("acc_highmet");
		if (acc_highmet_branch) {acc_highmet_branch->SetAddress(&acc_highmet_);}
	}
	acc_highht_branch = 0;
	if (tree->GetBranch("acc_highht") != 0) {
		acc_highht_branch = tree->GetBranch("acc_highht");
		if (acc_highht_branch) {acc_highht_branch->SetAddress(&acc_highht_);}
	}
	eldup_branch = 0;
	if (tree->GetBranch("eldup") != 0) {
		eldup_branch = tree->GetBranch("eldup");
		if (eldup_branch) {eldup_branch->SetAddress(&eldup_);}
	}
	csc_branch = 0;
	if (tree->GetBranch("csc") != 0) {
		csc_branch = tree->GetBranch("csc");
		if (csc_branch) {csc_branch->SetAddress(&csc_);}
	}
	hbhe_branch = 0;
	if (tree->GetBranch("hbhe") != 0) {
		hbhe_branch = tree->GetBranch("hbhe");
		if (hbhe_branch) {hbhe_branch->SetAddress(&hbhe_);}
	}
	hbhenew_branch = 0;
	if (tree->GetBranch("hbhenew") != 0) {
		hbhenew_branch = tree->GetBranch("hbhenew");
		if (hbhenew_branch) {hbhenew_branch->SetAddress(&hbhenew_);}
	}
	hcallaser_branch = 0;
	if (tree->GetBranch("hcallaser") != 0) {
		hcallaser_branch = tree->GetBranch("hcallaser");
		if (hcallaser_branch) {hcallaser_branch->SetAddress(&hcallaser_);}
	}
	ecaltp_branch = 0;
	if (tree->GetBranch("ecaltp") != 0) {
		ecaltp_branch = tree->GetBranch("ecaltp");
		if (ecaltp_branch) {ecaltp_branch->SetAddress(&ecaltp_);}
	}
	trkfail_branch = 0;
	if (tree->GetBranch("trkfail") != 0) {
		trkfail_branch = tree->GetBranch("trkfail");
		if (trkfail_branch) {trkfail_branch->SetAddress(&trkfail_);}
	}
	eebadsc_branch = 0;
	if (tree->GetBranch("eebadsc") != 0) {
		eebadsc_branch = tree->GetBranch("eebadsc");
		if (eebadsc_branch) {eebadsc_branch->SetAddress(&eebadsc_);}
	}
	lep1_badecallaser_branch = 0;
	if (tree->GetBranch("lep1_badecallaser") != 0) {
		lep1_badecallaser_branch = tree->GetBranch("lep1_badecallaser");
		if (lep1_badecallaser_branch) {lep1_badecallaser_branch->SetAddress(&lep1_badecallaser_);}
	}
	lep2_badecallaser_branch = 0;
	if (tree->GetBranch("lep2_badecallaser") != 0) {
		lep2_badecallaser_branch = tree->GetBranch("lep2_badecallaser");
		if (lep2_badecallaser_branch) {lep2_badecallaser_branch->SetAddress(&lep2_badecallaser_);}
	}
	isdata_branch = 0;
	if (tree->GetBranch("isdata") != 0) {
		isdata_branch = tree->GetBranch("isdata");
		if (isdata_branch) {isdata_branch->SetAddress(&isdata_);}
	}
	jetid_branch = 0;
	if (tree->GetBranch("jetid") != 0) {
		jetid_branch = tree->GetBranch("jetid");
		if (jetid_branch) {jetid_branch->SetAddress(&jetid_);}
	}
	jetid30_branch = 0;
	if (tree->GetBranch("jetid30") != 0) {
		jetid30_branch = tree->GetBranch("jetid30");
		if (jetid30_branch) {jetid30_branch->SetAddress(&jetid30_);}
	}
	json_branch = 0;
	if (tree->GetBranch("json") != 0) {
		json_branch = tree->GetBranch("json");
		if (json_branch) {json_branch->SetAddress(&json_);}
	}
	htoffset_branch = 0;
	if (tree->GetBranch("htoffset") != 0) {
		htoffset_branch = tree->GetBranch("htoffset");
		if (htoffset_branch) {htoffset_branch->SetAddress(&htoffset_);}
	}
	htuncor_branch = 0;
	if (tree->GetBranch("htuncor") != 0) {
		htuncor_branch = tree->GetBranch("htuncor");
		if (htuncor_branch) {htuncor_branch->SetAddress(&htuncor_);}
	}
	ptt_branch = 0;
	if (tree->GetBranch("ptt") != 0) {
		ptt_branch = tree->GetBranch("ptt");
		if (ptt_branch) {ptt_branch->SetAddress(&ptt_);}
	}
	pttbar_branch = 0;
	if (tree->GetBranch("pttbar") != 0) {
		pttbar_branch = tree->GetBranch("pttbar");
		if (pttbar_branch) {pttbar_branch->SetAddress(&pttbar_);}
	}
	ptttbar_branch = 0;
	if (tree->GetBranch("ptttbar") != 0) {
		ptttbar_branch = tree->GetBranch("ptttbar");
		if (ptttbar_branch) {ptttbar_branch->SetAddress(&ptttbar_);}
	}
	mttbar_branch = 0;
	if (tree->GetBranch("mttbar") != 0) {
		mttbar_branch = tree->GetBranch("mttbar");
		if (mttbar_branch) {mttbar_branch->SetAddress(&mttbar_);}
	}
	npartons_branch = 0;
	if (tree->GetBranch("npartons") != 0) {
		npartons_branch = tree->GetBranch("npartons");
		if (npartons_branch) {npartons_branch->SetAddress(&npartons_);}
	}
	nwzpartons_branch = 0;
	if (tree->GetBranch("nwzpartons") != 0) {
		nwzpartons_branch = tree->GetBranch("nwzpartons");
		if (nwzpartons_branch) {nwzpartons_branch->SetAddress(&nwzpartons_);}
	}
	hyptype_branch = 0;
	if (tree->GetBranch("hyptype") != 0) {
		hyptype_branch = tree->GetBranch("hyptype");
		if (hyptype_branch) {hyptype_branch->SetAddress(&hyptype_);}
	}
	maxpartonpt_branch = 0;
	if (tree->GetBranch("maxpartonpt") != 0) {
		maxpartonpt_branch = tree->GetBranch("maxpartonpt");
		if (maxpartonpt_branch) {maxpartonpt_branch->SetAddress(&maxpartonpt_);}
	}
	etattbar_branch = 0;
	if (tree->GetBranch("etattbar") != 0) {
		etattbar_branch = tree->GetBranch("etattbar");
		if (etattbar_branch) {etattbar_branch->SetAddress(&etattbar_);}
	}
	njetsoffset_branch = 0;
	if (tree->GetBranch("njetsoffset") != 0) {
		njetsoffset_branch = tree->GetBranch("njetsoffset");
		if (njetsoffset_branch) {njetsoffset_branch->SetAddress(&njetsoffset_);}
	}
	njetsuncor_branch = 0;
	if (tree->GetBranch("njetsuncor") != 0) {
		njetsuncor_branch = tree->GetBranch("njetsuncor");
		if (njetsuncor_branch) {njetsuncor_branch->SetAddress(&njetsuncor_);}
	}
	costhetaweight_branch = 0;
	if (tree->GetBranch("costhetaweight") != 0) {
		costhetaweight_branch = tree->GetBranch("costhetaweight");
		if (costhetaweight_branch) {costhetaweight_branch->SetAddress(&costhetaweight_);}
	}
	weight_branch = 0;
	if (tree->GetBranch("weight") != 0) {
		weight_branch = tree->GetBranch("weight");
		if (weight_branch) {weight_branch->SetAddress(&weight_);}
	}
	weightleft_branch = 0;
	if (tree->GetBranch("weightleft") != 0) {
		weightleft_branch = tree->GetBranch("weightleft");
		if (weightleft_branch) {weightleft_branch->SetAddress(&weightleft_);}
	}
	weightright_branch = 0;
	if (tree->GetBranch("weightright") != 0) {
		weightright_branch = tree->GetBranch("weightright");
		if (weightright_branch) {weightright_branch->SetAddress(&weightright_);}
	}
	mutrigweight_branch = 0;
	if (tree->GetBranch("mutrigweight") != 0) {
		mutrigweight_branch = tree->GetBranch("mutrigweight");
		if (mutrigweight_branch) {mutrigweight_branch->SetAddress(&mutrigweight_);}
	}
	mutrigweight2_branch = 0;
	if (tree->GetBranch("mutrigweight2") != 0) {
		mutrigweight2_branch = tree->GetBranch("mutrigweight2");
		if (mutrigweight2_branch) {mutrigweight2_branch->SetAddress(&mutrigweight2_);}
	}
	sltrigweight_branch = 0;
	if (tree->GetBranch("sltrigweight") != 0) {
		sltrigweight_branch = tree->GetBranch("sltrigweight");
		if (sltrigweight_branch) {sltrigweight_branch->SetAddress(&sltrigweight_);}
	}
	dltrigweight_branch = 0;
	if (tree->GetBranch("dltrigweight") != 0) {
		dltrigweight_branch = tree->GetBranch("dltrigweight");
		if (dltrigweight_branch) {dltrigweight_branch->SetAddress(&dltrigweight_);}
	}
	trgeff_branch = 0;
	if (tree->GetBranch("trgeff") != 0) {
		trgeff_branch = tree->GetBranch("trgeff");
		if (trgeff_branch) {trgeff_branch->SetAddress(&trgeff_);}
	}
	pthat_branch = 0;
	if (tree->GetBranch("pthat") != 0) {
		pthat_branch = tree->GetBranch("pthat");
		if (pthat_branch) {pthat_branch->SetAddress(&pthat_);}
	}
	qscale_branch = 0;
	if (tree->GetBranch("qscale") != 0) {
		qscale_branch = tree->GetBranch("qscale");
		if (qscale_branch) {qscale_branch->SetAddress(&qscale_);}
	}
	mgcor_branch = 0;
	if (tree->GetBranch("mgcor") != 0) {
		mgcor_branch = tree->GetBranch("mgcor");
		if (mgcor_branch) {mgcor_branch->SetAddress(&mgcor_);}
	}
	wflav_branch = 0;
	if (tree->GetBranch("wflav") != 0) {
		wflav_branch = tree->GetBranch("wflav");
		if (wflav_branch) {wflav_branch->SetAddress(&wflav_);}
	}
	ksusy_branch = 0;
	if (tree->GetBranch("ksusy") != 0) {
		ksusy_branch = tree->GetBranch("ksusy");
		if (ksusy_branch) {ksusy_branch->SetAddress(&ksusy_);}
	}
	ksusyup_branch = 0;
	if (tree->GetBranch("ksusyup") != 0) {
		ksusyup_branch = tree->GetBranch("ksusyup");
		if (ksusyup_branch) {ksusyup_branch->SetAddress(&ksusyup_);}
	}
	ksusydn_branch = 0;
	if (tree->GetBranch("ksusydn") != 0) {
		ksusydn_branch = tree->GetBranch("ksusydn");
		if (ksusydn_branch) {ksusydn_branch->SetAddress(&ksusydn_);}
	}
	xsecsusy_branch = 0;
	if (tree->GetBranch("xsecsusy") != 0) {
		xsecsusy_branch = tree->GetBranch("xsecsusy");
		if (xsecsusy_branch) {xsecsusy_branch->SetAddress(&xsecsusy_);}
	}
	xsecsusy2_branch = 0;
	if (tree->GetBranch("xsecsusy2") != 0) {
		xsecsusy2_branch = tree->GetBranch("xsecsusy2");
		if (xsecsusy2_branch) {xsecsusy2_branch->SetAddress(&xsecsusy2_);}
	}
	smeff_branch = 0;
	if (tree->GetBranch("smeff") != 0) {
		smeff_branch = tree->GetBranch("smeff");
		if (smeff_branch) {smeff_branch->SetAddress(&smeff_);}
	}
	k_branch = 0;
	if (tree->GetBranch("k") != 0) {
		k_branch = tree->GetBranch("k");
		if (k_branch) {k_branch->SetAddress(&k_);}
	}
	mllgen_branch = 0;
	if (tree->GetBranch("mllgen") != 0) {
		mllgen_branch = tree->GetBranch("mllgen");
		if (mllgen_branch) {mllgen_branch->SetAddress(&mllgen_);}
	}
	ptwgen_branch = 0;
	if (tree->GetBranch("ptwgen") != 0) {
		ptwgen_branch = tree->GetBranch("ptwgen");
		if (ptwgen_branch) {ptwgen_branch->SetAddress(&ptwgen_);}
	}
	ptzgen_branch = 0;
	if (tree->GetBranch("ptzgen") != 0) {
		ptzgen_branch = tree->GetBranch("ptzgen");
		if (ptzgen_branch) {ptzgen_branch->SetAddress(&ptzgen_);}
	}
	nlep_branch = 0;
	if (tree->GetBranch("nlep") != 0) {
		nlep_branch = tree->GetBranch("nlep");
		if (nlep_branch) {nlep_branch->SetAddress(&nlep_);}
	}
	nosel_branch = 0;
	if (tree->GetBranch("nosel") != 0) {
		nosel_branch = tree->GetBranch("nosel");
		if (nosel_branch) {nosel_branch->SetAddress(&nosel_);}
	}
	ngoodlep_branch = 0;
	if (tree->GetBranch("ngoodlep") != 0) {
		ngoodlep_branch = tree->GetBranch("ngoodlep");
		if (ngoodlep_branch) {ngoodlep_branch->SetAddress(&ngoodlep_);}
	}
	ngoodel_branch = 0;
	if (tree->GetBranch("ngoodel") != 0) {
		ngoodel_branch = tree->GetBranch("ngoodel");
		if (ngoodel_branch) {ngoodel_branch->SetAddress(&ngoodel_);}
	}
	ngoodmu_branch = 0;
	if (tree->GetBranch("ngoodmu") != 0) {
		ngoodmu_branch = tree->GetBranch("ngoodmu");
		if (ngoodmu_branch) {ngoodmu_branch->SetAddress(&ngoodmu_);}
	}
	mull_branch = 0;
	if (tree->GetBranch("mull") != 0) {
		mull_branch = tree->GetBranch("mull");
		if (mull_branch) {mull_branch->SetAddress(&mull_);}
	}
	mult_branch = 0;
	if (tree->GetBranch("mult") != 0) {
		mult_branch = tree->GetBranch("mult");
		if (mult_branch) {mult_branch->SetAddress(&mult_);}
	}
	mullgen_branch = 0;
	if (tree->GetBranch("mullgen") != 0) {
		mullgen_branch = tree->GetBranch("mullgen");
		if (mullgen_branch) {mullgen_branch->SetAddress(&mullgen_);}
	}
	multgen_branch = 0;
	if (tree->GetBranch("multgen") != 0) {
		multgen_branch = tree->GetBranch("multgen");
		if (multgen_branch) {multgen_branch->SetAddress(&multgen_);}
	}
	proc_branch = 0;
	if (tree->GetBranch("proc") != 0) {
		proc_branch = tree->GetBranch("proc");
		if (proc_branch) {proc_branch->SetAddress(&proc_);}
	}
	leptype_branch = 0;
	if (tree->GetBranch("leptype") != 0) {
		leptype_branch = tree->GetBranch("leptype");
		if (leptype_branch) {leptype_branch->SetAddress(&leptype_);}
	}
	topmass_branch = 0;
	if (tree->GetBranch("topmass") != 0) {
		topmass_branch = tree->GetBranch("topmass");
		if (topmass_branch) {topmass_branch->SetAddress(&topmass_);}
	}
	dilmass_branch = 0;
	if (tree->GetBranch("dilmass") != 0) {
		dilmass_branch = tree->GetBranch("dilmass");
		if (dilmass_branch) {dilmass_branch->SetAddress(&dilmass_);}
	}
	dilrecoil_branch = 0;
	if (tree->GetBranch("dilrecoil") != 0) {
		dilrecoil_branch = tree->GetBranch("dilrecoil");
		if (dilrecoil_branch) {dilrecoil_branch->SetAddress(&dilrecoil_);}
	}
	dilrecoilparl_branch = 0;
	if (tree->GetBranch("dilrecoilparl") != 0) {
		dilrecoilparl_branch = tree->GetBranch("dilrecoilparl");
		if (dilrecoilparl_branch) {dilrecoilparl_branch->SetAddress(&dilrecoilparl_);}
	}
	dilrecoilperp_branch = 0;
	if (tree->GetBranch("dilrecoilperp") != 0) {
		dilrecoilperp_branch = tree->GetBranch("dilrecoilperp");
		if (dilrecoilperp_branch) {dilrecoilperp_branch->SetAddress(&dilrecoilperp_);}
	}
	tcmet_branch = 0;
	if (tree->GetBranch("tcmet") != 0) {
		tcmet_branch = tree->GetBranch("tcmet");
		if (tcmet_branch) {tcmet_branch->SetAddress(&tcmet_);}
	}
	genmet_branch = 0;
	if (tree->GetBranch("genmet") != 0) {
		genmet_branch = tree->GetBranch("genmet");
		if (genmet_branch) {genmet_branch->SetAddress(&genmet_);}
	}
	gensumet_branch = 0;
	if (tree->GetBranch("gensumet") != 0) {
		gensumet_branch = tree->GetBranch("gensumet");
		if (gensumet_branch) {gensumet_branch->SetAddress(&gensumet_);}
	}
	genmetphi_branch = 0;
	if (tree->GetBranch("genmetphi") != 0) {
		genmetphi_branch = tree->GetBranch("genmetphi");
		if (genmetphi_branch) {genmetphi_branch->SetAddress(&genmetphi_);}
	}
	calomet_branch = 0;
	if (tree->GetBranch("calomet") != 0) {
		calomet_branch = tree->GetBranch("calomet");
		if (calomet_branch) {calomet_branch->SetAddress(&calomet_);}
	}
	calometphi_branch = 0;
	if (tree->GetBranch("calometphi") != 0) {
		calometphi_branch = tree->GetBranch("calometphi");
		if (calometphi_branch) {calometphi_branch->SetAddress(&calometphi_);}
	}
	trkmet_branch = 0;
	if (tree->GetBranch("trkmet") != 0) {
		trkmet_branch = tree->GetBranch("trkmet");
		if (trkmet_branch) {trkmet_branch->SetAddress(&trkmet_);}
	}
	trkmetphi_branch = 0;
	if (tree->GetBranch("trkmetphi") != 0) {
		trkmetphi_branch = tree->GetBranch("trkmetphi");
		if (trkmetphi_branch) {trkmetphi_branch->SetAddress(&trkmetphi_);}
	}
	pfmet_branch = 0;
	if (tree->GetBranch("pfmet") != 0) {
		pfmet_branch = tree->GetBranch("pfmet");
		if (pfmet_branch) {pfmet_branch->SetAddress(&pfmet_);}
	}
	pfmetveto_branch = 0;
	if (tree->GetBranch("pfmetveto") != 0) {
		pfmetveto_branch = tree->GetBranch("pfmetveto");
		if (pfmetveto_branch) {pfmetveto_branch->SetAddress(&pfmetveto_);}
	}
	pfmetsig_branch = 0;
	if (tree->GetBranch("pfmetsig") != 0) {
		pfmetsig_branch = tree->GetBranch("pfmetsig");
		if (pfmetsig_branch) {pfmetsig_branch->SetAddress(&pfmetsig_);}
	}
	pfmetphi_branch = 0;
	if (tree->GetBranch("pfmetphi") != 0) {
		pfmetphi_branch = tree->GetBranch("pfmetphi");
		if (pfmetphi_branch) {pfmetphi_branch->SetAddress(&pfmetphi_);}
	}
	pfsumet_branch = 0;
	if (tree->GetBranch("pfsumet") != 0) {
		pfsumet_branch = tree->GetBranch("pfsumet");
		if (pfsumet_branch) {pfsumet_branch->SetAddress(&pfsumet_);}
	}
	mucormet_branch = 0;
	if (tree->GetBranch("mucormet") != 0) {
		mucormet_branch = tree->GetBranch("mucormet");
		if (mucormet_branch) {mucormet_branch->SetAddress(&mucormet_);}
	}
	mucorjesmet_branch = 0;
	if (tree->GetBranch("mucorjesmet") != 0) {
		mucorjesmet_branch = tree->GetBranch("mucorjesmet");
		if (mucorjesmet_branch) {mucorjesmet_branch->SetAddress(&mucorjesmet_);}
	}
	tcmet35X_branch = 0;
	if (tree->GetBranch("tcmet35X") != 0) {
		tcmet35X_branch = tree->GetBranch("tcmet35X");
		if (tcmet35X_branch) {tcmet35X_branch->SetAddress(&tcmet35X_);}
	}
	tcmetevent_branch = 0;
	if (tree->GetBranch("tcmetevent") != 0) {
		tcmetevent_branch = tree->GetBranch("tcmetevent");
		if (tcmetevent_branch) {tcmetevent_branch->SetAddress(&tcmetevent_);}
	}
	tcmetlooper_branch = 0;
	if (tree->GetBranch("tcmetlooper") != 0) {
		tcmetlooper_branch = tree->GetBranch("tcmetlooper");
		if (tcmetlooper_branch) {tcmetlooper_branch->SetAddress(&tcmetlooper_);}
	}
	tcmetphi_branch = 0;
	if (tree->GetBranch("tcmetphi") != 0) {
		tcmetphi_branch = tree->GetBranch("tcmetphi");
		if (tcmetphi_branch) {tcmetphi_branch->SetAddress(&tcmetphi_);}
	}
	tcsumet_branch = 0;
	if (tree->GetBranch("tcsumet") != 0) {
		tcsumet_branch = tree->GetBranch("tcsumet");
		if (tcsumet_branch) {tcsumet_branch->SetAddress(&tcsumet_);}
	}
	tcmetUp_branch = 0;
	if (tree->GetBranch("tcmetUp") != 0) {
		tcmetUp_branch = tree->GetBranch("tcmetUp");
		if (tcmetUp_branch) {tcmetUp_branch->SetAddress(&tcmetUp_);}
	}
	tcmetDown_branch = 0;
	if (tree->GetBranch("tcmetDown") != 0) {
		tcmetDown_branch = tree->GetBranch("tcmetDown");
		if (tcmetDown_branch) {tcmetDown_branch->SetAddress(&tcmetDown_);}
	}
	tcmetTest_branch = 0;
	if (tree->GetBranch("tcmetTest") != 0) {
		tcmetTest_branch = tree->GetBranch("tcmetTest");
		if (tcmetTest_branch) {tcmetTest_branch->SetAddress(&tcmetTest_);}
	}
	pfmetUp_branch = 0;
	if (tree->GetBranch("pfmetUp") != 0) {
		pfmetUp_branch = tree->GetBranch("pfmetUp");
		if (pfmetUp_branch) {pfmetUp_branch->SetAddress(&pfmetUp_);}
	}
	pfmetDown_branch = 0;
	if (tree->GetBranch("pfmetDown") != 0) {
		pfmetDown_branch = tree->GetBranch("pfmetDown");
		if (pfmetDown_branch) {pfmetDown_branch->SetAddress(&pfmetDown_);}
	}
	pfmetTest_branch = 0;
	if (tree->GetBranch("pfmetTest") != 0) {
		pfmetTest_branch = tree->GetBranch("pfmetTest");
		if (pfmetTest_branch) {pfmetTest_branch->SetAddress(&pfmetTest_);}
	}
	sumjetpt_branch = 0;
	if (tree->GetBranch("sumjetpt") != 0) {
		sumjetpt_branch = tree->GetBranch("sumjetpt");
		if (sumjetpt_branch) {sumjetpt_branch->SetAddress(&sumjetpt_);}
	}
	dileta_branch = 0;
	if (tree->GetBranch("dileta") != 0) {
		dileta_branch = tree->GetBranch("dileta");
		if (dileta_branch) {dileta_branch->SetAddress(&dileta_);}
	}
	dilpt_branch = 0;
	if (tree->GetBranch("dilpt") != 0) {
		dilpt_branch = tree->GetBranch("dilpt");
		if (dilpt_branch) {dilpt_branch->SetAddress(&dilpt_);}
	}
	dildphi_branch = 0;
	if (tree->GetBranch("dildphi") != 0) {
		dildphi_branch = tree->GetBranch("dildphi");
		if (dildphi_branch) {dildphi_branch->SetAddress(&dildphi_);}
	}
	ngenjets_branch = 0;
	if (tree->GetBranch("ngenjets") != 0) {
		ngenjets_branch = tree->GetBranch("ngenjets");
		if (ngenjets_branch) {ngenjets_branch->SetAddress(&ngenjets_);}
	}
	njpt_branch = 0;
	if (tree->GetBranch("njpt") != 0) {
		njpt_branch = tree->GetBranch("njpt");
		if (njpt_branch) {njpt_branch->SetAddress(&njpt_);}
	}
	trgmu1_branch = 0;
	if (tree->GetBranch("trgmu1") != 0) {
		trgmu1_branch = tree->GetBranch("trgmu1");
		if (trgmu1_branch) {trgmu1_branch->SetAddress(&trgmu1_);}
	}
	trgmu2_branch = 0;
	if (tree->GetBranch("trgmu2") != 0) {
		trgmu2_branch = tree->GetBranch("trgmu2");
		if (trgmu2_branch) {trgmu2_branch->SetAddress(&trgmu2_);}
	}
	trgel1_branch = 0;
	if (tree->GetBranch("trgel1") != 0) {
		trgel1_branch = tree->GetBranch("trgel1");
		if (trgel1_branch) {trgel1_branch->SetAddress(&trgel1_);}
	}
	trgel2_branch = 0;
	if (tree->GetBranch("trgel2") != 0) {
		trgel2_branch = tree->GetBranch("trgel2");
		if (trgel2_branch) {trgel2_branch->SetAddress(&trgel2_);}
	}
	isomu24_branch = 0;
	if (tree->GetBranch("isomu24") != 0) {
		isomu24_branch = tree->GetBranch("isomu24");
		if (isomu24_branch) {isomu24_branch->SetAddress(&isomu24_);}
	}
	ele27wp80_branch = 0;
	if (tree->GetBranch("ele27wp80") != 0) {
		ele27wp80_branch = tree->GetBranch("ele27wp80");
		if (ele27wp80_branch) {ele27wp80_branch->SetAddress(&ele27wp80_);}
	}
	mm_branch = 0;
	if (tree->GetBranch("mm") != 0) {
		mm_branch = tree->GetBranch("mm");
		if (mm_branch) {mm_branch->SetAddress(&mm_);}
	}
	mmtk_branch = 0;
	if (tree->GetBranch("mmtk") != 0) {
		mmtk_branch = tree->GetBranch("mmtk");
		if (mmtk_branch) {mmtk_branch->SetAddress(&mmtk_);}
	}
	me_branch = 0;
	if (tree->GetBranch("me") != 0) {
		me_branch = tree->GetBranch("me");
		if (me_branch) {me_branch->SetAddress(&me_);}
	}
	em_branch = 0;
	if (tree->GetBranch("em") != 0) {
		em_branch = tree->GetBranch("em");
		if (em_branch) {em_branch->SetAddress(&em_);}
	}
	mu_branch = 0;
	if (tree->GetBranch("mu") != 0) {
		mu_branch = tree->GetBranch("mu");
		if (mu_branch) {mu_branch->SetAddress(&mu_);}
	}
	ee_branch = 0;
	if (tree->GetBranch("ee") != 0) {
		ee_branch = tree->GetBranch("ee");
		if (ee_branch) {ee_branch->SetAddress(&ee_);}
	}
	npfjets30_branch = 0;
	if (tree->GetBranch("npfjets30") != 0) {
		npfjets30_branch = tree->GetBranch("npfjets30");
		if (npfjets30_branch) {npfjets30_branch->SetAddress(&npfjets30_);}
	}
	npfjets30lepcorr_branch = 0;
	if (tree->GetBranch("npfjets30lepcorr") != 0) {
		npfjets30lepcorr_branch = tree->GetBranch("npfjets30lepcorr");
		if (npfjets30lepcorr_branch) {npfjets30lepcorr_branch->SetAddress(&npfjets30lepcorr_);}
	}
	knjets_branch = 0;
	if (tree->GetBranch("knjets") != 0) {
		knjets_branch = tree->GetBranch("knjets");
		if (knjets_branch) {knjets_branch->SetAddress(&knjets_);}
	}
	rhovor_branch = 0;
	if (tree->GetBranch("rhovor") != 0) {
		rhovor_branch = tree->GetBranch("rhovor");
		if (rhovor_branch) {rhovor_branch->SetAddress(&rhovor_);}
	}
	htpf30_branch = 0;
	if (tree->GetBranch("htpf30") != 0) {
		htpf30_branch = tree->GetBranch("htpf30");
		if (htpf30_branch) {htpf30_branch->SetAddress(&htpf30_);}
	}
	t1met10_branch = 0;
	if (tree->GetBranch("t1met10") != 0) {
		t1met10_branch = tree->GetBranch("t1met10");
		if (t1met10_branch) {t1met10_branch->SetAddress(&t1met10_);}
	}
	t1met20_branch = 0;
	if (tree->GetBranch("t1met20") != 0) {
		t1met20_branch = tree->GetBranch("t1met20");
		if (t1met20_branch) {t1met20_branch->SetAddress(&t1met20_);}
	}
	t1met30_branch = 0;
	if (tree->GetBranch("t1met30") != 0) {
		t1met30_branch = tree->GetBranch("t1met30");
		if (t1met30_branch) {t1met30_branch->SetAddress(&t1met30_);}
	}
	t1met10phi_branch = 0;
	if (tree->GetBranch("t1met10phi") != 0) {
		t1met10phi_branch = tree->GetBranch("t1met10phi");
		if (t1met10phi_branch) {t1met10phi_branch->SetAddress(&t1met10phi_);}
	}
	t1met20phi_branch = 0;
	if (tree->GetBranch("t1met20phi") != 0) {
		t1met20phi_branch = tree->GetBranch("t1met20phi");
		if (t1met20phi_branch) {t1met20phi_branch->SetAddress(&t1met20phi_);}
	}
	t1met30phi_branch = 0;
	if (tree->GetBranch("t1met30phi") != 0) {
		t1met30phi_branch = tree->GetBranch("t1met30phi");
		if (t1met30phi_branch) {t1met30phi_branch->SetAddress(&t1met30phi_);}
	}
	t1met10mt_branch = 0;
	if (tree->GetBranch("t1met10mt") != 0) {
		t1met10mt_branch = tree->GetBranch("t1met10mt");
		if (t1met10mt_branch) {t1met10mt_branch->SetAddress(&t1met10mt_);}
	}
	t1met20mt_branch = 0;
	if (tree->GetBranch("t1met20mt") != 0) {
		t1met20mt_branch = tree->GetBranch("t1met20mt");
		if (t1met20mt_branch) {t1met20mt_branch->SetAddress(&t1met20mt_);}
	}
	t1met30mt_branch = 0;
	if (tree->GetBranch("t1met30mt") != 0) {
		t1met30mt_branch = tree->GetBranch("t1met30mt");
		if (t1met30mt_branch) {t1met30mt_branch->SetAddress(&t1met30mt_);}
	}
	lepmetpt_branch = 0;
	if (tree->GetBranch("lepmetpt") != 0) {
		lepmetpt_branch = tree->GetBranch("lepmetpt");
		if (lepmetpt_branch) {lepmetpt_branch->SetAddress(&lepmetpt_);}
	}
	lept1met10pt_branch = 0;
	if (tree->GetBranch("lept1met10pt") != 0) {
		lept1met10pt_branch = tree->GetBranch("lept1met10pt");
		if (lept1met10pt_branch) {lept1met10pt_branch->SetAddress(&lept1met10pt_);}
	}
	t1met10s_branch = 0;
	if (tree->GetBranch("t1met10s") != 0) {
		t1met10s_branch = tree->GetBranch("t1met10s");
		if (t1met10s_branch) {t1met10s_branch->SetAddress(&t1met10s_);}
	}
	t1met10sphi_branch = 0;
	if (tree->GetBranch("t1met10sphi") != 0) {
		t1met10sphi_branch = tree->GetBranch("t1met10sphi");
		if (t1met10sphi_branch) {t1met10sphi_branch->SetAddress(&t1met10sphi_);}
	}
	t1met10smt_branch = 0;
	if (tree->GetBranch("t1met10smt") != 0) {
		t1met10smt_branch = tree->GetBranch("t1met10smt");
		if (t1met10smt_branch) {t1met10smt_branch->SetAddress(&t1met10smt_);}
	}
	t1metphicorr_branch = 0;
	if (tree->GetBranch("t1metphicorr") != 0) {
		t1metphicorr_branch = tree->GetBranch("t1metphicorr");
		if (t1metphicorr_branch) {t1metphicorr_branch->SetAddress(&t1metphicorr_);}
	}
	t1metphicorrup_branch = 0;
	if (tree->GetBranch("t1metphicorrup") != 0) {
		t1metphicorrup_branch = tree->GetBranch("t1metphicorrup");
		if (t1metphicorrup_branch) {t1metphicorrup_branch->SetAddress(&t1metphicorrup_);}
	}
	t1metphicorrdn_branch = 0;
	if (tree->GetBranch("t1metphicorrdn") != 0) {
		t1metphicorrdn_branch = tree->GetBranch("t1metphicorrdn");
		if (t1metphicorrdn_branch) {t1metphicorrdn_branch->SetAddress(&t1metphicorrdn_);}
	}
	t1metphicorrphi_branch = 0;
	if (tree->GetBranch("t1metphicorrphi") != 0) {
		t1metphicorrphi_branch = tree->GetBranch("t1metphicorrphi");
		if (t1metphicorrphi_branch) {t1metphicorrphi_branch->SetAddress(&t1metphicorrphi_);}
	}
	t1metphicorrphiup_branch = 0;
	if (tree->GetBranch("t1metphicorrphiup") != 0) {
		t1metphicorrphiup_branch = tree->GetBranch("t1metphicorrphiup");
		if (t1metphicorrphiup_branch) {t1metphicorrphiup_branch->SetAddress(&t1metphicorrphiup_);}
	}
	t1metphicorrphidn_branch = 0;
	if (tree->GetBranch("t1metphicorrphidn") != 0) {
		t1metphicorrphidn_branch = tree->GetBranch("t1metphicorrphidn");
		if (t1metphicorrphidn_branch) {t1metphicorrphidn_branch->SetAddress(&t1metphicorrphidn_);}
	}
	t1metphicorrlep_branch = 0;
	if (tree->GetBranch("t1metphicorrlep") != 0) {
		t1metphicorrlep_branch = tree->GetBranch("t1metphicorrlep");
		if (t1metphicorrlep_branch) {t1metphicorrlep_branch->SetAddress(&t1metphicorrlep_);}
	}
	t1metphicorrlepphi_branch = 0;
	if (tree->GetBranch("t1metphicorrlepphi") != 0) {
		t1metphicorrlepphi_branch = tree->GetBranch("t1metphicorrlepphi");
		if (t1metphicorrlepphi_branch) {t1metphicorrlepphi_branch->SetAddress(&t1metphicorrlepphi_);}
	}
	t1metphicorrmt_branch = 0;
	if (tree->GetBranch("t1metphicorrmt") != 0) {
		t1metphicorrmt_branch = tree->GetBranch("t1metphicorrmt");
		if (t1metphicorrmt_branch) {t1metphicorrmt_branch->SetAddress(&t1metphicorrmt_);}
	}
	t1metphicorrmtup_branch = 0;
	if (tree->GetBranch("t1metphicorrmtup") != 0) {
		t1metphicorrmtup_branch = tree->GetBranch("t1metphicorrmtup");
		if (t1metphicorrmtup_branch) {t1metphicorrmtup_branch->SetAddress(&t1metphicorrmtup_);}
	}
	t1metphicorrmtdn_branch = 0;
	if (tree->GetBranch("t1metphicorrmtdn") != 0) {
		t1metphicorrmtdn_branch = tree->GetBranch("t1metphicorrmtdn");
		if (t1metphicorrmtdn_branch) {t1metphicorrmtdn_branch->SetAddress(&t1metphicorrmtdn_);}
	}
	t1metphicorrlepmt_branch = 0;
	if (tree->GetBranch("t1metphicorrlepmt") != 0) {
		t1metphicorrlepmt_branch = tree->GetBranch("t1metphicorrlepmt");
		if (t1metphicorrlepmt_branch) {t1metphicorrlepmt_branch->SetAddress(&t1metphicorrlepmt_);}
	}
	t1met_off_branch = 0;
	if (tree->GetBranch("t1met_off") != 0) {
		t1met_off_branch = tree->GetBranch("t1met_off");
		if (t1met_off_branch) {t1met_off_branch->SetAddress(&t1met_off_);}
	}
	t1metphi_off_branch = 0;
	if (tree->GetBranch("t1metphi_off") != 0) {
		t1metphi_off_branch = tree->GetBranch("t1metphi_off");
		if (t1metphi_off_branch) {t1metphi_off_branch->SetAddress(&t1metphi_off_);}
	}
	t1metmt_off_branch = 0;
	if (tree->GetBranch("t1metmt_off") != 0) {
		t1metmt_off_branch = tree->GetBranch("t1metmt_off");
		if (t1metmt_off_branch) {t1metmt_off_branch->SetAddress(&t1metmt_off_);}
	}
	t1metphicorr_off_branch = 0;
	if (tree->GetBranch("t1metphicorr_off") != 0) {
		t1metphicorr_off_branch = tree->GetBranch("t1metphicorr_off");
		if (t1metphicorr_off_branch) {t1metphicorr_off_branch->SetAddress(&t1metphicorr_off_);}
	}
	t1metphicorrphi_off_branch = 0;
	if (tree->GetBranch("t1metphicorrphi_off") != 0) {
		t1metphicorrphi_off_branch = tree->GetBranch("t1metphicorrphi_off");
		if (t1metphicorrphi_off_branch) {t1metphicorrphi_off_branch->SetAddress(&t1metphicorrphi_off_);}
	}
	t1metphicorrmt_off_branch = 0;
	if (tree->GetBranch("t1metphicorrmt_off") != 0) {
		t1metphicorrmt_off_branch = tree->GetBranch("t1metphicorrmt_off");
		if (t1metphicorrmt_off_branch) {t1metphicorrmt_off_branch->SetAddress(&t1metphicorrmt_off_);}
	}
	mht15_branch = 0;
	if (tree->GetBranch("mht15") != 0) {
		mht15_branch = tree->GetBranch("mht15");
		if (mht15_branch) {mht15_branch->SetAddress(&mht15_);}
	}
	mht15phi_branch = 0;
	if (tree->GetBranch("mht15phi") != 0) {
		mht15phi_branch = tree->GetBranch("mht15phi");
		if (mht15phi_branch) {mht15phi_branch->SetAddress(&mht15phi_);}
	}
	trkmet_mht15_branch = 0;
	if (tree->GetBranch("trkmet_mht15") != 0) {
		trkmet_mht15_branch = tree->GetBranch("trkmet_mht15");
		if (trkmet_mht15_branch) {trkmet_mht15_branch->SetAddress(&trkmet_mht15_);}
	}
	trkmetphi_mht15_branch = 0;
	if (tree->GetBranch("trkmetphi_mht15") != 0) {
		trkmetphi_mht15_branch = tree->GetBranch("trkmetphi_mht15");
		if (trkmetphi_mht15_branch) {trkmetphi_mht15_branch->SetAddress(&trkmetphi_mht15_);}
	}
	mettlj15_branch = 0;
	if (tree->GetBranch("mettlj15") != 0) {
		mettlj15_branch = tree->GetBranch("mettlj15");
		if (mettlj15_branch) {mettlj15_branch->SetAddress(&mettlj15_);}
	}
	mettlj15phi_branch = 0;
	if (tree->GetBranch("mettlj15phi") != 0) {
		mettlj15phi_branch = tree->GetBranch("mettlj15phi");
		if (mettlj15phi_branch) {mettlj15phi_branch->SetAddress(&mettlj15phi_);}
	}
	mt2bmin_branch = 0;
	if (tree->GetBranch("mt2bmin") != 0) {
		mt2bmin_branch = tree->GetBranch("mt2bmin");
		if (mt2bmin_branch) {mt2bmin_branch->SetAddress(&mt2bmin_);}
	}
	mt2blmin_branch = 0;
	if (tree->GetBranch("mt2blmin") != 0) {
		mt2blmin_branch = tree->GetBranch("mt2blmin");
		if (mt2blmin_branch) {mt2blmin_branch->SetAddress(&mt2blmin_);}
	}
	mt2wmin_branch = 0;
	if (tree->GetBranch("mt2wmin") != 0) {
		mt2wmin_branch = tree->GetBranch("mt2wmin");
		if (mt2wmin_branch) {mt2wmin_branch->SetAddress(&mt2wmin_);}
	}
	chi2min_branch = 0;
	if (tree->GetBranch("chi2min") != 0) {
		chi2min_branch = tree->GetBranch("chi2min");
		if (chi2min_branch) {chi2min_branch->SetAddress(&chi2min_);}
	}
	chi2minprob_branch = 0;
	if (tree->GetBranch("chi2minprob") != 0) {
		chi2minprob_branch = tree->GetBranch("chi2minprob");
		if (chi2minprob_branch) {chi2minprob_branch->SetAddress(&chi2minprob_);}
	}
	nbtagsssv_branch = 0;
	if (tree->GetBranch("nbtagsssv") != 0) {
		nbtagsssv_branch = tree->GetBranch("nbtagsssv");
		if (nbtagsssv_branch) {nbtagsssv_branch->SetAddress(&nbtagsssv_);}
	}
	nbtagstcl_branch = 0;
	if (tree->GetBranch("nbtagstcl") != 0) {
		nbtagstcl_branch = tree->GetBranch("nbtagstcl");
		if (nbtagstcl_branch) {nbtagstcl_branch->SetAddress(&nbtagstcl_);}
	}
	nbtagstcm_branch = 0;
	if (tree->GetBranch("nbtagstcm") != 0) {
		nbtagstcm_branch = tree->GetBranch("nbtagstcm");
		if (nbtagstcm_branch) {nbtagstcm_branch->SetAddress(&nbtagstcm_);}
	}
	nbtagscsvl_branch = 0;
	if (tree->GetBranch("nbtagscsvl") != 0) {
		nbtagscsvl_branch = tree->GetBranch("nbtagscsvl");
		if (nbtagscsvl_branch) {nbtagscsvl_branch->SetAddress(&nbtagscsvl_);}
	}
	nbtagscsvm_branch = 0;
	if (tree->GetBranch("nbtagscsvm") != 0) {
		nbtagscsvm_branch = tree->GetBranch("nbtagscsvm");
		if (nbtagscsvm_branch) {nbtagscsvm_branch->SetAddress(&nbtagscsvm_);}
	}
	nbtagscsvt_branch = 0;
	if (tree->GetBranch("nbtagscsvt") != 0) {
		nbtagscsvt_branch = tree->GetBranch("nbtagscsvt");
		if (nbtagscsvt_branch) {nbtagscsvt_branch->SetAddress(&nbtagscsvt_);}
	}
	nbtagsssvcorr_branch = 0;
	if (tree->GetBranch("nbtagsssvcorr") != 0) {
		nbtagsssvcorr_branch = tree->GetBranch("nbtagsssvcorr");
		if (nbtagsssvcorr_branch) {nbtagsssvcorr_branch->SetAddress(&nbtagsssvcorr_);}
	}
	nbtagstclcorr_branch = 0;
	if (tree->GetBranch("nbtagstclcorr") != 0) {
		nbtagstclcorr_branch = tree->GetBranch("nbtagstclcorr");
		if (nbtagstclcorr_branch) {nbtagstclcorr_branch->SetAddress(&nbtagstclcorr_);}
	}
	nbtagstcmcorr_branch = 0;
	if (tree->GetBranch("nbtagstcmcorr") != 0) {
		nbtagstcmcorr_branch = tree->GetBranch("nbtagstcmcorr");
		if (nbtagstcmcorr_branch) {nbtagstcmcorr_branch->SetAddress(&nbtagstcmcorr_);}
	}
	nbtagscsvlcorr_branch = 0;
	if (tree->GetBranch("nbtagscsvlcorr") != 0) {
		nbtagscsvlcorr_branch = tree->GetBranch("nbtagscsvlcorr");
		if (nbtagscsvlcorr_branch) {nbtagscsvlcorr_branch->SetAddress(&nbtagscsvlcorr_);}
	}
	nbtagscsvmcorr_branch = 0;
	if (tree->GetBranch("nbtagscsvmcorr") != 0) {
		nbtagscsvmcorr_branch = tree->GetBranch("nbtagscsvmcorr");
		if (nbtagscsvmcorr_branch) {nbtagscsvmcorr_branch->SetAddress(&nbtagscsvmcorr_);}
	}
	nbtagscsvtcott_branch = 0;
	if (tree->GetBranch("nbtagscsvtcott") != 0) {
		nbtagscsvtcott_branch = tree->GetBranch("nbtagscsvtcott");
		if (nbtagscsvtcott_branch) {nbtagscsvtcott_branch->SetAddress(&nbtagscsvtcott_);}
	}
	njetsUp_branch = 0;
	if (tree->GetBranch("njetsUp") != 0) {
		njetsUp_branch = tree->GetBranch("njetsUp");
		if (njetsUp_branch) {njetsUp_branch->SetAddress(&njetsUp_);}
	}
	njetsDown_branch = 0;
	if (tree->GetBranch("njetsDown") != 0) {
		njetsDown_branch = tree->GetBranch("njetsDown");
		if (njetsDown_branch) {njetsDown_branch->SetAddress(&njetsDown_);}
	}
	htUp_branch = 0;
	if (tree->GetBranch("htUp") != 0) {
		htUp_branch = tree->GetBranch("htUp");
		if (htUp_branch) {htUp_branch->SetAddress(&htUp_);}
	}
	htDown_branch = 0;
	if (tree->GetBranch("htDown") != 0) {
		htDown_branch = tree->GetBranch("htDown");
		if (htDown_branch) {htDown_branch->SetAddress(&htDown_);}
	}
	ntruepu_branch = 0;
	if (tree->GetBranch("ntruepu") != 0) {
		ntruepu_branch = tree->GetBranch("ntruepu");
		if (ntruepu_branch) {ntruepu_branch->SetAddress(&ntruepu_);}
	}
	npu_branch = 0;
	if (tree->GetBranch("npu") != 0) {
		npu_branch = tree->GetBranch("npu");
		if (npu_branch) {npu_branch->SetAddress(&npu_);}
	}
	npuMinusOne_branch = 0;
	if (tree->GetBranch("npuMinusOne") != 0) {
		npuMinusOne_branch = tree->GetBranch("npuMinusOne");
		if (npuMinusOne_branch) {npuMinusOne_branch->SetAddress(&npuMinusOne_);}
	}
	npuPlusOne_branch = 0;
	if (tree->GetBranch("npuPlusOne") != 0) {
		npuPlusOne_branch = tree->GetBranch("npuPlusOne");
		if (npuPlusOne_branch) {npuPlusOne_branch->SetAddress(&npuPlusOne_);}
	}
	nvtx_branch = 0;
	if (tree->GetBranch("nvtx") != 0) {
		nvtx_branch = tree->GetBranch("nvtx");
		if (nvtx_branch) {nvtx_branch->SetAddress(&nvtx_);}
	}
	indexfirstGoodVertex__branch = 0;
	if (tree->GetBranch("indexfirstGoodVertex_") != 0) {
		indexfirstGoodVertex__branch = tree->GetBranch("indexfirstGoodVertex_");
		if (indexfirstGoodVertex__branch) {indexfirstGoodVertex__branch->SetAddress(&indexfirstGoodVertex__);}
	}
	nvtxweight_branch = 0;
	if (tree->GetBranch("nvtxweight") != 0) {
		nvtxweight_branch = tree->GetBranch("nvtxweight");
		if (nvtxweight_branch) {nvtxweight_branch->SetAddress(&nvtxweight_);}
	}
	n3dvtxweight_branch = 0;
	if (tree->GetBranch("n3dvtxweight") != 0) {
		n3dvtxweight_branch = tree->GetBranch("n3dvtxweight");
		if (n3dvtxweight_branch) {n3dvtxweight_branch->SetAddress(&n3dvtxweight_);}
	}
	pdfid1_branch = 0;
	if (tree->GetBranch("pdfid1") != 0) {
		pdfid1_branch = tree->GetBranch("pdfid1");
		if (pdfid1_branch) {pdfid1_branch->SetAddress(&pdfid1_);}
	}
	pdfid2_branch = 0;
	if (tree->GetBranch("pdfid2") != 0) {
		pdfid2_branch = tree->GetBranch("pdfid2");
		if (pdfid2_branch) {pdfid2_branch->SetAddress(&pdfid2_);}
	}
	pdfx1_branch = 0;
	if (tree->GetBranch("pdfx1") != 0) {
		pdfx1_branch = tree->GetBranch("pdfx1");
		if (pdfx1_branch) {pdfx1_branch->SetAddress(&pdfx1_);}
	}
	pdfx2_branch = 0;
	if (tree->GetBranch("pdfx2") != 0) {
		pdfx2_branch = tree->GetBranch("pdfx2");
		if (pdfx2_branch) {pdfx2_branch->SetAddress(&pdfx2_);}
	}
	pdfQ_branch = 0;
	if (tree->GetBranch("pdfQ") != 0) {
		pdfQ_branch = tree->GetBranch("pdfQ");
		if (pdfQ_branch) {pdfQ_branch->SetAddress(&pdfQ_);}
	}
	vecjetpt_branch = 0;
	if (tree->GetBranch("vecjetpt") != 0) {
		vecjetpt_branch = tree->GetBranch("vecjetpt");
		if (vecjetpt_branch) {vecjetpt_branch->SetAddress(&vecjetpt_);}
	}
	pass_branch = 0;
	if (tree->GetBranch("pass") != 0) {
		pass_branch = tree->GetBranch("pass");
		if (pass_branch) {pass_branch->SetAddress(&pass_);}
	}
	passz_branch = 0;
	if (tree->GetBranch("passz") != 0) {
		passz_branch = tree->GetBranch("passz");
		if (passz_branch) {passz_branch->SetAddress(&passz_);}
	}
	m0_branch = 0;
	if (tree->GetBranch("m0") != 0) {
		m0_branch = tree->GetBranch("m0");
		if (m0_branch) {m0_branch->SetAddress(&m0_);}
	}
	mg_branch = 0;
	if (tree->GetBranch("mg") != 0) {
		mg_branch = tree->GetBranch("mg");
		if (mg_branch) {mg_branch->SetAddress(&mg_);}
	}
	ml_branch = 0;
	if (tree->GetBranch("ml") != 0) {
		ml_branch = tree->GetBranch("ml");
		if (ml_branch) {ml_branch->SetAddress(&ml_);}
	}
	x_branch = 0;
	if (tree->GetBranch("x") != 0) {
		x_branch = tree->GetBranch("x");
		if (x_branch) {x_branch->SetAddress(&x_);}
	}
	m12_branch = 0;
	if (tree->GetBranch("m12") != 0) {
		m12_branch = tree->GetBranch("m12");
		if (m12_branch) {m12_branch->SetAddress(&m12_);}
	}
	lep1chi2ndf_branch = 0;
	if (tree->GetBranch("lep1chi2ndf") != 0) {
		lep1chi2ndf_branch = tree->GetBranch("lep1chi2ndf");
		if (lep1chi2ndf_branch) {lep1chi2ndf_branch->SetAddress(&lep1chi2ndf_);}
	}
	lep2chi2ndf_branch = 0;
	if (tree->GetBranch("lep2chi2ndf") != 0) {
		lep2chi2ndf_branch = tree->GetBranch("lep2chi2ndf");
		if (lep2chi2ndf_branch) {lep2chi2ndf_branch->SetAddress(&lep2chi2ndf_);}
	}
	lep1dpt_branch = 0;
	if (tree->GetBranch("lep1dpt") != 0) {
		lep1dpt_branch = tree->GetBranch("lep1dpt");
		if (lep1dpt_branch) {lep1dpt_branch->SetAddress(&lep1dpt_);}
	}
	lep2dpt_branch = 0;
	if (tree->GetBranch("lep2dpt") != 0) {
		lep2dpt_branch = tree->GetBranch("lep2dpt");
		if (lep2dpt_branch) {lep2dpt_branch->SetAddress(&lep2dpt_);}
	}
	id1_branch = 0;
	if (tree->GetBranch("id1") != 0) {
		id1_branch = tree->GetBranch("id1");
		if (id1_branch) {id1_branch->SetAddress(&id1_);}
	}
	id2_branch = 0;
	if (tree->GetBranch("id2") != 0) {
		id2_branch = tree->GetBranch("id2");
		if (id2_branch) {id2_branch->SetAddress(&id2_);}
	}
	leptype1_branch = 0;
	if (tree->GetBranch("leptype1") != 0) {
		leptype1_branch = tree->GetBranch("leptype1");
		if (leptype1_branch) {leptype1_branch->SetAddress(&leptype1_);}
	}
	leptype2_branch = 0;
	if (tree->GetBranch("leptype2") != 0) {
		leptype2_branch = tree->GetBranch("leptype2");
		if (leptype2_branch) {leptype2_branch->SetAddress(&leptype2_);}
	}
	w1_branch = 0;
	if (tree->GetBranch("w1") != 0) {
		w1_branch = tree->GetBranch("w1");
		if (w1_branch) {w1_branch->SetAddress(&w1_);}
	}
	w2_branch = 0;
	if (tree->GetBranch("w2") != 0) {
		w2_branch = tree->GetBranch("w2");
		if (w2_branch) {w2_branch->SetAddress(&w2_);}
	}
	iso1_branch = 0;
	if (tree->GetBranch("iso1") != 0) {
		iso1_branch = tree->GetBranch("iso1");
		if (iso1_branch) {iso1_branch->SetAddress(&iso1_);}
	}
	isont1_branch = 0;
	if (tree->GetBranch("isont1") != 0) {
		isont1_branch = tree->GetBranch("isont1");
		if (isont1_branch) {isont1_branch->SetAddress(&isont1_);}
	}
	isopfold1_branch = 0;
	if (tree->GetBranch("isopfold1") != 0) {
		isopfold1_branch = tree->GetBranch("isopfold1");
		if (isopfold1_branch) {isopfold1_branch->SetAddress(&isopfold1_);}
	}
	isopf1_branch = 0;
	if (tree->GetBranch("isopf1") != 0) {
		isopf1_branch = tree->GetBranch("isopf1");
		if (isopf1_branch) {isopf1_branch->SetAddress(&isopf1_);}
	}
	etasc1_branch = 0;
	if (tree->GetBranch("etasc1") != 0) {
		etasc1_branch = tree->GetBranch("etasc1");
		if (etasc1_branch) {etasc1_branch->SetAddress(&etasc1_);}
	}
	etasc2_branch = 0;
	if (tree->GetBranch("etasc2") != 0) {
		etasc2_branch = tree->GetBranch("etasc2");
		if (etasc2_branch) {etasc2_branch->SetAddress(&etasc2_);}
	}
	eoverpin_branch = 0;
	if (tree->GetBranch("eoverpin") != 0) {
		eoverpin_branch = tree->GetBranch("eoverpin");
		if (eoverpin_branch) {eoverpin_branch->SetAddress(&eoverpin_);}
	}
	eoverpout_branch = 0;
	if (tree->GetBranch("eoverpout") != 0) {
		eoverpout_branch = tree->GetBranch("eoverpout");
		if (eoverpout_branch) {eoverpout_branch->SetAddress(&eoverpout_);}
	}
	dEtaIn_branch = 0;
	if (tree->GetBranch("dEtaIn") != 0) {
		dEtaIn_branch = tree->GetBranch("dEtaIn");
		if (dEtaIn_branch) {dEtaIn_branch->SetAddress(&dEtaIn_);}
	}
	dPhiIn_branch = 0;
	if (tree->GetBranch("dPhiIn") != 0) {
		dPhiIn_branch = tree->GetBranch("dPhiIn");
		if (dPhiIn_branch) {dPhiIn_branch->SetAddress(&dPhiIn_);}
	}
	sigmaIEtaIEta_branch = 0;
	if (tree->GetBranch("sigmaIEtaIEta") != 0) {
		sigmaIEtaIEta_branch = tree->GetBranch("sigmaIEtaIEta");
		if (sigmaIEtaIEta_branch) {sigmaIEtaIEta_branch->SetAddress(&sigmaIEtaIEta_);}
	}
	hOverE_branch = 0;
	if (tree->GetBranch("hOverE") != 0) {
		hOverE_branch = tree->GetBranch("hOverE");
		if (hOverE_branch) {hOverE_branch->SetAddress(&hOverE_);}
	}
	ooemoop_branch = 0;
	if (tree->GetBranch("ooemoop") != 0) {
		ooemoop_branch = tree->GetBranch("ooemoop");
		if (ooemoop_branch) {ooemoop_branch->SetAddress(&ooemoop_);}
	}
	d0vtx_branch = 0;
	if (tree->GetBranch("d0vtx") != 0) {
		d0vtx_branch = tree->GetBranch("d0vtx");
		if (d0vtx_branch) {d0vtx_branch->SetAddress(&d0vtx_);}
	}
	dzvtx_branch = 0;
	if (tree->GetBranch("dzvtx") != 0) {
		dzvtx_branch = tree->GetBranch("dzvtx");
		if (dzvtx_branch) {dzvtx_branch->SetAddress(&dzvtx_);}
	}
	expinnerlayers_branch = 0;
	if (tree->GetBranch("expinnerlayers") != 0) {
		expinnerlayers_branch = tree->GetBranch("expinnerlayers");
		if (expinnerlayers_branch) {expinnerlayers_branch->SetAddress(&expinnerlayers_);}
	}
	fbrem_branch = 0;
	if (tree->GetBranch("fbrem") != 0) {
		fbrem_branch = tree->GetBranch("fbrem");
		if (fbrem_branch) {fbrem_branch->SetAddress(&fbrem_);}
	}
	pfisoch_branch = 0;
	if (tree->GetBranch("pfisoch") != 0) {
		pfisoch_branch = tree->GetBranch("pfisoch");
		if (pfisoch_branch) {pfisoch_branch->SetAddress(&pfisoch_);}
	}
	pfisoem_branch = 0;
	if (tree->GetBranch("pfisoem") != 0) {
		pfisoem_branch = tree->GetBranch("pfisoem");
		if (pfisoem_branch) {pfisoem_branch->SetAddress(&pfisoem_);}
	}
	pfisonh_branch = 0;
	if (tree->GetBranch("pfisonh") != 0) {
		pfisonh_branch = tree->GetBranch("pfisonh");
		if (pfisonh_branch) {pfisonh_branch->SetAddress(&pfisonh_);}
	}
	eSC_branch = 0;
	if (tree->GetBranch("eSC") != 0) {
		eSC_branch = tree->GetBranch("eSC");
		if (eSC_branch) {eSC_branch->SetAddress(&eSC_);}
	}
	phiSC_branch = 0;
	if (tree->GetBranch("phiSC") != 0) {
		phiSC_branch = tree->GetBranch("phiSC");
		if (phiSC_branch) {phiSC_branch->SetAddress(&phiSC_);}
	}
	eSCRaw_branch = 0;
	if (tree->GetBranch("eSCRaw") != 0) {
		eSCRaw_branch = tree->GetBranch("eSCRaw");
		if (eSCRaw_branch) {eSCRaw_branch->SetAddress(&eSCRaw_);}
	}
	eSCPresh_branch = 0;
	if (tree->GetBranch("eSCPresh") != 0) {
		eSCPresh_branch = tree->GetBranch("eSCPresh");
		if (eSCPresh_branch) {eSCPresh_branch->SetAddress(&eSCPresh_);}
	}
	lep1_scslasercormean_branch = 0;
	if (tree->GetBranch("lep1_scslasercormean") != 0) {
		lep1_scslasercormean_branch = tree->GetBranch("lep1_scslasercormean");
		if (lep1_scslasercormean_branch) {lep1_scslasercormean_branch->SetAddress(&lep1_scslasercormean_);}
	}
	lep1_scslasercormax_branch = 0;
	if (tree->GetBranch("lep1_scslasercormax") != 0) {
		lep1_scslasercormax_branch = tree->GetBranch("lep1_scslasercormax");
		if (lep1_scslasercormax_branch) {lep1_scslasercormax_branch->SetAddress(&lep1_scslasercormax_);}
	}
	eoverpin2_branch = 0;
	if (tree->GetBranch("eoverpin2") != 0) {
		eoverpin2_branch = tree->GetBranch("eoverpin2");
		if (eoverpin2_branch) {eoverpin2_branch->SetAddress(&eoverpin2_);}
	}
	eoverpout2_branch = 0;
	if (tree->GetBranch("eoverpout2") != 0) {
		eoverpout2_branch = tree->GetBranch("eoverpout2");
		if (eoverpout2_branch) {eoverpout2_branch->SetAddress(&eoverpout2_);}
	}
	dEtaIn2_branch = 0;
	if (tree->GetBranch("dEtaIn2") != 0) {
		dEtaIn2_branch = tree->GetBranch("dEtaIn2");
		if (dEtaIn2_branch) {dEtaIn2_branch->SetAddress(&dEtaIn2_);}
	}
	dPhiIn2_branch = 0;
	if (tree->GetBranch("dPhiIn2") != 0) {
		dPhiIn2_branch = tree->GetBranch("dPhiIn2");
		if (dPhiIn2_branch) {dPhiIn2_branch->SetAddress(&dPhiIn2_);}
	}
	sigmaIEtaIEta2_branch = 0;
	if (tree->GetBranch("sigmaIEtaIEta2") != 0) {
		sigmaIEtaIEta2_branch = tree->GetBranch("sigmaIEtaIEta2");
		if (sigmaIEtaIEta2_branch) {sigmaIEtaIEta2_branch->SetAddress(&sigmaIEtaIEta2_);}
	}
	hOverE2_branch = 0;
	if (tree->GetBranch("hOverE2") != 0) {
		hOverE2_branch = tree->GetBranch("hOverE2");
		if (hOverE2_branch) {hOverE2_branch->SetAddress(&hOverE2_);}
	}
	ooemoop2_branch = 0;
	if (tree->GetBranch("ooemoop2") != 0) {
		ooemoop2_branch = tree->GetBranch("ooemoop2");
		if (ooemoop2_branch) {ooemoop2_branch->SetAddress(&ooemoop2_);}
	}
	d0vtx2_branch = 0;
	if (tree->GetBranch("d0vtx2") != 0) {
		d0vtx2_branch = tree->GetBranch("d0vtx2");
		if (d0vtx2_branch) {d0vtx2_branch->SetAddress(&d0vtx2_);}
	}
	dzvtx2_branch = 0;
	if (tree->GetBranch("dzvtx2") != 0) {
		dzvtx2_branch = tree->GetBranch("dzvtx2");
		if (dzvtx2_branch) {dzvtx2_branch->SetAddress(&dzvtx2_);}
	}
	expinnerlayers2_branch = 0;
	if (tree->GetBranch("expinnerlayers2") != 0) {
		expinnerlayers2_branch = tree->GetBranch("expinnerlayers2");
		if (expinnerlayers2_branch) {expinnerlayers2_branch->SetAddress(&expinnerlayers2_);}
	}
	fbrem2_branch = 0;
	if (tree->GetBranch("fbrem2") != 0) {
		fbrem2_branch = tree->GetBranch("fbrem2");
		if (fbrem2_branch) {fbrem2_branch->SetAddress(&fbrem2_);}
	}
	pfisoch2_branch = 0;
	if (tree->GetBranch("pfisoch2") != 0) {
		pfisoch2_branch = tree->GetBranch("pfisoch2");
		if (pfisoch2_branch) {pfisoch2_branch->SetAddress(&pfisoch2_);}
	}
	pfisoem2_branch = 0;
	if (tree->GetBranch("pfisoem2") != 0) {
		pfisoem2_branch = tree->GetBranch("pfisoem2");
		if (pfisoem2_branch) {pfisoem2_branch->SetAddress(&pfisoem2_);}
	}
	pfisonh2_branch = 0;
	if (tree->GetBranch("pfisonh2") != 0) {
		pfisonh2_branch = tree->GetBranch("pfisonh2");
		if (pfisonh2_branch) {pfisonh2_branch->SetAddress(&pfisonh2_);}
	}
	eSC2_branch = 0;
	if (tree->GetBranch("eSC2") != 0) {
		eSC2_branch = tree->GetBranch("eSC2");
		if (eSC2_branch) {eSC2_branch->SetAddress(&eSC2_);}
	}
	phiSC2_branch = 0;
	if (tree->GetBranch("phiSC2") != 0) {
		phiSC2_branch = tree->GetBranch("phiSC2");
		if (phiSC2_branch) {phiSC2_branch->SetAddress(&phiSC2_);}
	}
	eSCRaw2_branch = 0;
	if (tree->GetBranch("eSCRaw2") != 0) {
		eSCRaw2_branch = tree->GetBranch("eSCRaw2");
		if (eSCRaw2_branch) {eSCRaw2_branch->SetAddress(&eSCRaw2_);}
	}
	eSCPresh2_branch = 0;
	if (tree->GetBranch("eSCPresh2") != 0) {
		eSCPresh2_branch = tree->GetBranch("eSCPresh2");
		if (eSCPresh2_branch) {eSCPresh2_branch->SetAddress(&eSCPresh2_);}
	}
	lep2_scslasercormean_branch = 0;
	if (tree->GetBranch("lep2_scslasercormean") != 0) {
		lep2_scslasercormean_branch = tree->GetBranch("lep2_scslasercormean");
		if (lep2_scslasercormean_branch) {lep2_scslasercormean_branch->SetAddress(&lep2_scslasercormean_);}
	}
	lep2_scslasercormax_branch = 0;
	if (tree->GetBranch("lep2_scslasercormax") != 0) {
		lep2_scslasercormax_branch = tree->GetBranch("lep2_scslasercormax");
		if (lep2_scslasercormax_branch) {lep2_scslasercormax_branch->SetAddress(&lep2_scslasercormax_);}
	}
	scslasercormax_branch = 0;
	if (tree->GetBranch("scslasercormax") != 0) {
		scslasercormax_branch = tree->GetBranch("scslasercormax");
		if (scslasercormax_branch) {scslasercormax_branch->SetAddress(&scslasercormax_);}
	}
	scslasercormax_pt_branch = 0;
	if (tree->GetBranch("scslasercormax_pt") != 0) {
		scslasercormax_pt_branch = tree->GetBranch("scslasercormax_pt");
		if (scslasercormax_pt_branch) {scslasercormax_pt_branch->SetAddress(&scslasercormax_pt_);}
	}
	scslasercormax_eta_branch = 0;
	if (tree->GetBranch("scslasercormax_eta") != 0) {
		scslasercormax_eta_branch = tree->GetBranch("scslasercormax_eta");
		if (scslasercormax_eta_branch) {scslasercormax_eta_branch->SetAddress(&scslasercormax_eta_);}
	}
	iso2_branch = 0;
	if (tree->GetBranch("iso2") != 0) {
		iso2_branch = tree->GetBranch("iso2");
		if (iso2_branch) {iso2_branch->SetAddress(&iso2_);}
	}
	ecalveto1_branch = 0;
	if (tree->GetBranch("ecalveto1") != 0) {
		ecalveto1_branch = tree->GetBranch("ecalveto1");
		if (ecalveto1_branch) {ecalveto1_branch->SetAddress(&ecalveto1_);}
	}
	ecalveto2_branch = 0;
	if (tree->GetBranch("ecalveto2") != 0) {
		ecalveto2_branch = tree->GetBranch("ecalveto2");
		if (ecalveto2_branch) {ecalveto2_branch->SetAddress(&ecalveto2_);}
	}
	hcalveto1_branch = 0;
	if (tree->GetBranch("hcalveto1") != 0) {
		hcalveto1_branch = tree->GetBranch("hcalveto1");
		if (hcalveto1_branch) {hcalveto1_branch->SetAddress(&hcalveto1_);}
	}
	hcalveto2_branch = 0;
	if (tree->GetBranch("hcalveto2") != 0) {
		hcalveto2_branch = tree->GetBranch("hcalveto2");
		if (hcalveto2_branch) {hcalveto2_branch->SetAddress(&hcalveto2_);}
	}
	isont2_branch = 0;
	if (tree->GetBranch("isont2") != 0) {
		isont2_branch = tree->GetBranch("isont2");
		if (isont2_branch) {isont2_branch->SetAddress(&isont2_);}
	}
	isopf2_branch = 0;
	if (tree->GetBranch("isopf2") != 0) {
		isopf2_branch = tree->GetBranch("isopf2");
		if (isopf2_branch) {isopf2_branch->SetAddress(&isopf2_);}
	}
	ptl1_branch = 0;
	if (tree->GetBranch("ptl1") != 0) {
		ptl1_branch = tree->GetBranch("ptl1");
		if (ptl1_branch) {ptl1_branch->SetAddress(&ptl1_);}
	}
	ptl2_branch = 0;
	if (tree->GetBranch("ptl2") != 0) {
		ptl2_branch = tree->GetBranch("ptl2");
		if (ptl2_branch) {ptl2_branch->SetAddress(&ptl2_);}
	}
	etal1_branch = 0;
	if (tree->GetBranch("etal1") != 0) {
		etal1_branch = tree->GetBranch("etal1");
		if (etal1_branch) {etal1_branch->SetAddress(&etal1_);}
	}
	etal2_branch = 0;
	if (tree->GetBranch("etal2") != 0) {
		etal2_branch = tree->GetBranch("etal2");
		if (etal2_branch) {etal2_branch->SetAddress(&etal2_);}
	}
	phil1_branch = 0;
	if (tree->GetBranch("phil1") != 0) {
		phil1_branch = tree->GetBranch("phil1");
		if (phil1_branch) {phil1_branch->SetAddress(&phil1_);}
	}
	phil2_branch = 0;
	if (tree->GetBranch("phil2") != 0) {
		phil2_branch = tree->GetBranch("phil2");
		if (phil2_branch) {phil2_branch->SetAddress(&phil2_);}
	}
	meff_branch = 0;
	if (tree->GetBranch("meff") != 0) {
		meff_branch = tree->GetBranch("meff");
		if (meff_branch) {meff_branch->SetAddress(&meff_);}
	}
	mt_branch = 0;
	if (tree->GetBranch("mt") != 0) {
		mt_branch = tree->GetBranch("mt");
		if (mt_branch) {mt_branch->SetAddress(&mt_);}
	}
	run_branch = 0;
	if (tree->GetBranch("run") != 0) {
		run_branch = tree->GetBranch("run");
		if (run_branch) {run_branch->SetAddress(&run_);}
	}
	lumi_branch = 0;
	if (tree->GetBranch("lumi") != 0) {
		lumi_branch = tree->GetBranch("lumi");
		if (lumi_branch) {lumi_branch->SetAddress(&lumi_);}
	}
	event_branch = 0;
	if (tree->GetBranch("event") != 0) {
		event_branch = tree->GetBranch("event");
		if (event_branch) {event_branch->SetAddress(&event_);}
	}
	y_branch = 0;
	if (tree->GetBranch("y") != 0) {
		y_branch = tree->GetBranch("y");
		if (y_branch) {y_branch->SetAddress(&y_);}
	}
	ht_branch = 0;
	if (tree->GetBranch("ht") != 0) {
		ht_branch = tree->GetBranch("ht");
		if (ht_branch) {ht_branch->SetAddress(&ht_);}
	}
	htgen_branch = 0;
	if (tree->GetBranch("htgen") != 0) {
		htgen_branch = tree->GetBranch("htgen");
		if (htgen_branch) {htgen_branch->SetAddress(&htgen_);}
	}
	htjpt_branch = 0;
	if (tree->GetBranch("htjpt") != 0) {
		htjpt_branch = tree->GetBranch("htjpt");
		if (htjpt_branch) {htjpt_branch->SetAddress(&htjpt_);}
	}
	nels_branch = 0;
	if (tree->GetBranch("nels") != 0) {
		nels_branch = tree->GetBranch("nels");
		if (nels_branch) {nels_branch->SetAddress(&nels_);}
	}
	nmus_branch = 0;
	if (tree->GetBranch("nmus") != 0) {
		nmus_branch = tree->GetBranch("nmus");
		if (nmus_branch) {nmus_branch->SetAddress(&nmus_);}
	}
	ntaus_branch = 0;
	if (tree->GetBranch("ntaus") != 0) {
		ntaus_branch = tree->GetBranch("ntaus");
		if (ntaus_branch) {ntaus_branch->SetAddress(&ntaus_);}
	}
	nleps_branch = 0;
	if (tree->GetBranch("nleps") != 0) {
		nleps_branch = tree->GetBranch("nleps");
		if (nleps_branch) {nleps_branch->SetAddress(&nleps_);}
	}
	nbs_branch = 0;
	if (tree->GetBranch("nbs") != 0) {
		nbs_branch = tree->GetBranch("nbs");
		if (nbs_branch) {nbs_branch->SetAddress(&nbs_);}
	}
	dphijm_branch = 0;
	if (tree->GetBranch("dphijm") != 0) {
		dphijm_branch = tree->GetBranch("dphijm");
		if (dphijm_branch) {dphijm_branch->SetAddress(&dphijm_);}
	}
	ptjetraw_branch = 0;
	if (tree->GetBranch("ptjetraw") != 0) {
		ptjetraw_branch = tree->GetBranch("ptjetraw");
		if (ptjetraw_branch) {ptjetraw_branch->SetAddress(&ptjetraw_);}
	}
	ptjet23_branch = 0;
	if (tree->GetBranch("ptjet23") != 0) {
		ptjet23_branch = tree->GetBranch("ptjet23");
		if (ptjet23_branch) {ptjet23_branch->SetAddress(&ptjet23_);}
	}
	ptjetF23_branch = 0;
	if (tree->GetBranch("ptjetF23") != 0) {
		ptjetF23_branch = tree->GetBranch("ptjetF23");
		if (ptjetF23_branch) {ptjetF23_branch->SetAddress(&ptjetF23_);}
	}
	ptjetO23_branch = 0;
	if (tree->GetBranch("ptjetO23") != 0) {
		ptjetO23_branch = tree->GetBranch("ptjetO23");
		if (ptjetO23_branch) {ptjetO23_branch->SetAddress(&ptjetO23_);}
	}
	mcid1_branch = 0;
	if (tree->GetBranch("mcid1") != 0) {
		mcid1_branch = tree->GetBranch("mcid1");
		if (mcid1_branch) {mcid1_branch->SetAddress(&mcid1_);}
	}
	mcdr1_branch = 0;
	if (tree->GetBranch("mcdr1") != 0) {
		mcdr1_branch = tree->GetBranch("mcdr1");
		if (mcdr1_branch) {mcdr1_branch->SetAddress(&mcdr1_);}
	}
	mcdecay1_branch = 0;
	if (tree->GetBranch("mcdecay1") != 0) {
		mcdecay1_branch = tree->GetBranch("mcdecay1");
		if (mcdecay1_branch) {mcdecay1_branch->SetAddress(&mcdecay1_);}
	}
	mcndec1_branch = 0;
	if (tree->GetBranch("mcndec1") != 0) {
		mcndec1_branch = tree->GetBranch("mcndec1");
		if (mcndec1_branch) {mcndec1_branch->SetAddress(&mcndec1_);}
	}
	mcndec2_branch = 0;
	if (tree->GetBranch("mcndec2") != 0) {
		mcndec2_branch = tree->GetBranch("mcndec2");
		if (mcndec2_branch) {mcndec2_branch->SetAddress(&mcndec2_);}
	}
	mcndeckls1_branch = 0;
	if (tree->GetBranch("mcndeckls1") != 0) {
		mcndeckls1_branch = tree->GetBranch("mcndeckls1");
		if (mcndeckls1_branch) {mcndeckls1_branch->SetAddress(&mcndeckls1_);}
	}
	mcndeckls2_branch = 0;
	if (tree->GetBranch("mcndeckls2") != 0) {
		mcndeckls2_branch = tree->GetBranch("mcndeckls2");
		if (mcndeckls2_branch) {mcndeckls2_branch->SetAddress(&mcndeckls2_);}
	}
	mcndecem1_branch = 0;
	if (tree->GetBranch("mcndecem1") != 0) {
		mcndecem1_branch = tree->GetBranch("mcndecem1");
		if (mcndecem1_branch) {mcndecem1_branch->SetAddress(&mcndecem1_);}
	}
	mcndecem2_branch = 0;
	if (tree->GetBranch("mcndecem2") != 0) {
		mcndecem2_branch = tree->GetBranch("mcndecem2");
		if (mcndecem2_branch) {mcndecem2_branch->SetAddress(&mcndecem2_);}
	}
	mcid2_branch = 0;
	if (tree->GetBranch("mcid2") != 0) {
		mcid2_branch = tree->GetBranch("mcid2");
		if (mcid2_branch) {mcid2_branch->SetAddress(&mcid2_);}
	}
	mcdr2_branch = 0;
	if (tree->GetBranch("mcdr2") != 0) {
		mcdr2_branch = tree->GetBranch("mcdr2");
		if (mcdr2_branch) {mcdr2_branch->SetAddress(&mcdr2_);}
	}
	mcdecay2_branch = 0;
	if (tree->GetBranch("mcdecay2") != 0) {
		mcdecay2_branch = tree->GetBranch("mcdecay2");
		if (mcdecay2_branch) {mcdecay2_branch->SetAddress(&mcdecay2_);}
	}
	mctaudpt1_branch = 0;
	if (tree->GetBranch("mctaudpt1") != 0) {
		mctaudpt1_branch = tree->GetBranch("mctaudpt1");
		if (mctaudpt1_branch) {mctaudpt1_branch->SetAddress(&mctaudpt1_);}
	}
	mctaudpt2_branch = 0;
	if (tree->GetBranch("mctaudpt2") != 0) {
		mctaudpt2_branch = tree->GetBranch("mctaudpt2");
		if (mctaudpt2_branch) {mctaudpt2_branch->SetAddress(&mctaudpt2_);}
	}
	mctaudid1_branch = 0;
	if (tree->GetBranch("mctaudid1") != 0) {
		mctaudid1_branch = tree->GetBranch("mctaudid1");
		if (mctaudid1_branch) {mctaudid1_branch->SetAddress(&mctaudid1_);}
	}
	mctaudid2_branch = 0;
	if (tree->GetBranch("mctaudid2") != 0) {
		mctaudid2_branch = tree->GetBranch("mctaudid2");
		if (mctaudid2_branch) {mctaudid2_branch->SetAddress(&mctaudid2_);}
	}
	mlepid_branch = 0;
	if (tree->GetBranch("mlepid") != 0) {
		mlepid_branch = tree->GetBranch("mlepid");
		if (mlepid_branch) {mlepid_branch->SetAddress(&mlepid_);}
	}
	mleppassid_branch = 0;
	if (tree->GetBranch("mleppassid") != 0) {
		mleppassid_branch = tree->GetBranch("mleppassid");
		if (mleppassid_branch) {mleppassid_branch->SetAddress(&mleppassid_);}
	}
	mleppassiso_branch = 0;
	if (tree->GetBranch("mleppassiso") != 0) {
		mleppassiso_branch = tree->GetBranch("mleppassiso");
		if (mleppassiso_branch) {mleppassiso_branch->SetAddress(&mleppassiso_);}
	}
	mlepiso_branch = 0;
	if (tree->GetBranch("mlepiso") != 0) {
		mlepiso_branch = tree->GetBranch("mlepiso");
		if (mlepiso_branch) {mlepiso_branch->SetAddress(&mlepiso_);}
	}
	mlepdr_branch = 0;
	if (tree->GetBranch("mlepdr") != 0) {
		mlepdr_branch = tree->GetBranch("mlepdr");
		if (mlepdr_branch) {mlepdr_branch->SetAddress(&mlepdr_);}
	}
	pflepiso_branch = 0;
	if (tree->GetBranch("pflepiso") != 0) {
		pflepiso_branch = tree->GetBranch("pflepiso");
		if (pflepiso_branch) {pflepiso_branch->SetAddress(&pflepiso_);}
	}
	pflepdr_branch = 0;
	if (tree->GetBranch("pflepdr") != 0) {
		pflepdr_branch = tree->GetBranch("pflepdr");
		if (pflepdr_branch) {pflepdr_branch->SetAddress(&pflepdr_);}
	}
	pfleppt_branch = 0;
	if (tree->GetBranch("pfleppt") != 0) {
		pfleppt_branch = tree->GetBranch("pfleppt");
		if (pfleppt_branch) {pfleppt_branch->SetAddress(&pfleppt_);}
	}
	pflepmindrj_branch = 0;
	if (tree->GetBranch("pflepmindrj") != 0) {
		pflepmindrj_branch = tree->GetBranch("pflepmindrj");
		if (pflepmindrj_branch) {pflepmindrj_branch->SetAddress(&pflepmindrj_);}
	}
	pftaudiso_branch = 0;
	if (tree->GetBranch("pftaudiso") != 0) {
		pftaudiso_branch = tree->GetBranch("pftaudiso");
		if (pftaudiso_branch) {pftaudiso_branch->SetAddress(&pftaudiso_);}
	}
	pftauddr_branch = 0;
	if (tree->GetBranch("pftauddr") != 0) {
		pftauddr_branch = tree->GetBranch("pftauddr");
		if (pftauddr_branch) {pftauddr_branch->SetAddress(&pftauddr_);}
	}
	pftaudpt_branch = 0;
	if (tree->GetBranch("pftaudpt") != 0) {
		pftaudpt_branch = tree->GetBranch("pftaudpt");
		if (pftaudpt_branch) {pftaudpt_branch->SetAddress(&pftaudpt_);}
	}
	pftaudmindrj_branch = 0;
	if (tree->GetBranch("pftaudmindrj") != 0) {
		pftaudmindrj_branch = tree->GetBranch("pftaudmindrj");
		if (pftaudmindrj_branch) {pftaudmindrj_branch->SetAddress(&pftaudmindrj_);}
	}
	pfcandid5_branch = 0;
	if (tree->GetBranch("pfcandid5") != 0) {
		pfcandid5_branch = tree->GetBranch("pfcandid5");
		if (pfcandid5_branch) {pfcandid5_branch->SetAddress(&pfcandid5_);}
	}
	pfcandiso5_branch = 0;
	if (tree->GetBranch("pfcandiso5") != 0) {
		pfcandiso5_branch = tree->GetBranch("pfcandiso5");
		if (pfcandiso5_branch) {pfcandiso5_branch->SetAddress(&pfcandiso5_);}
	}
	pfcandpt5_branch = 0;
	if (tree->GetBranch("pfcandpt5") != 0) {
		pfcandpt5_branch = tree->GetBranch("pfcandpt5");
		if (pfcandpt5_branch) {pfcandpt5_branch->SetAddress(&pfcandpt5_);}
	}
	pfcanddz5_branch = 0;
	if (tree->GetBranch("pfcanddz5") != 0) {
		pfcanddz5_branch = tree->GetBranch("pfcanddz5");
		if (pfcanddz5_branch) {pfcanddz5_branch->SetAddress(&pfcanddz5_);}
	}
	pfcandmindrj5_branch = 0;
	if (tree->GetBranch("pfcandmindrj5") != 0) {
		pfcandmindrj5_branch = tree->GetBranch("pfcandmindrj5");
		if (pfcandmindrj5_branch) {pfcandmindrj5_branch->SetAddress(&pfcandmindrj5_);}
	}
	pfcandid10_branch = 0;
	if (tree->GetBranch("pfcandid10") != 0) {
		pfcandid10_branch = tree->GetBranch("pfcandid10");
		if (pfcandid10_branch) {pfcandid10_branch->SetAddress(&pfcandid10_);}
	}
	pfcandiso10_branch = 0;
	if (tree->GetBranch("pfcandiso10") != 0) {
		pfcandiso10_branch = tree->GetBranch("pfcandiso10");
		if (pfcandiso10_branch) {pfcandiso10_branch->SetAddress(&pfcandiso10_);}
	}
	pfcandpt10_branch = 0;
	if (tree->GetBranch("pfcandpt10") != 0) {
		pfcandpt10_branch = tree->GetBranch("pfcandpt10");
		if (pfcandpt10_branch) {pfcandpt10_branch->SetAddress(&pfcandpt10_);}
	}
	pfcanddz10_branch = 0;
	if (tree->GetBranch("pfcanddz10") != 0) {
		pfcanddz10_branch = tree->GetBranch("pfcanddz10");
		if (pfcanddz10_branch) {pfcanddz10_branch->SetAddress(&pfcanddz10_);}
	}
	pfcandmindrj10_branch = 0;
	if (tree->GetBranch("pfcandmindrj10") != 0) {
		pfcandmindrj10_branch = tree->GetBranch("pfcandmindrj10");
		if (pfcandmindrj10_branch) {pfcandmindrj10_branch->SetAddress(&pfcandmindrj10_);}
	}
	pfcandidOS10_branch = 0;
	if (tree->GetBranch("pfcandidOS10") != 0) {
		pfcandidOS10_branch = tree->GetBranch("pfcandidOS10");
		if (pfcandidOS10_branch) {pfcandidOS10_branch->SetAddress(&pfcandidOS10_);}
	}
	pfcandisoOS10_branch = 0;
	if (tree->GetBranch("pfcandisoOS10") != 0) {
		pfcandisoOS10_branch = tree->GetBranch("pfcandisoOS10");
		if (pfcandisoOS10_branch) {pfcandisoOS10_branch->SetAddress(&pfcandisoOS10_);}
	}
	pfcandptOS10_branch = 0;
	if (tree->GetBranch("pfcandptOS10") != 0) {
		pfcandptOS10_branch = tree->GetBranch("pfcandptOS10");
		if (pfcandptOS10_branch) {pfcandptOS10_branch->SetAddress(&pfcandptOS10_);}
	}
	pfcanddzOS10_branch = 0;
	if (tree->GetBranch("pfcanddzOS10") != 0) {
		pfcanddzOS10_branch = tree->GetBranch("pfcanddzOS10");
		if (pfcanddzOS10_branch) {pfcanddzOS10_branch->SetAddress(&pfcanddzOS10_);}
	}
	pfcandid5looseZ_branch = 0;
	if (tree->GetBranch("pfcandid5looseZ") != 0) {
		pfcandid5looseZ_branch = tree->GetBranch("pfcandid5looseZ");
		if (pfcandid5looseZ_branch) {pfcandid5looseZ_branch->SetAddress(&pfcandid5looseZ_);}
	}
	pfcandiso5looseZ_branch = 0;
	if (tree->GetBranch("pfcandiso5looseZ") != 0) {
		pfcandiso5looseZ_branch = tree->GetBranch("pfcandiso5looseZ");
		if (pfcandiso5looseZ_branch) {pfcandiso5looseZ_branch->SetAddress(&pfcandiso5looseZ_);}
	}
	pfcandpt5looseZ_branch = 0;
	if (tree->GetBranch("pfcandpt5looseZ") != 0) {
		pfcandpt5looseZ_branch = tree->GetBranch("pfcandpt5looseZ");
		if (pfcandpt5looseZ_branch) {pfcandpt5looseZ_branch->SetAddress(&pfcandpt5looseZ_);}
	}
	pfcanddz5looseZ_branch = 0;
	if (tree->GetBranch("pfcanddz5looseZ") != 0) {
		pfcanddz5looseZ_branch = tree->GetBranch("pfcanddz5looseZ");
		if (pfcanddz5looseZ_branch) {pfcanddz5looseZ_branch->SetAddress(&pfcanddz5looseZ_);}
	}
	pfcandidOS10looseZ_branch = 0;
	if (tree->GetBranch("pfcandidOS10looseZ") != 0) {
		pfcandidOS10looseZ_branch = tree->GetBranch("pfcandidOS10looseZ");
		if (pfcandidOS10looseZ_branch) {pfcandidOS10looseZ_branch->SetAddress(&pfcandidOS10looseZ_);}
	}
	pfcandisoOS10looseZ_branch = 0;
	if (tree->GetBranch("pfcandisoOS10looseZ") != 0) {
		pfcandisoOS10looseZ_branch = tree->GetBranch("pfcandisoOS10looseZ");
		if (pfcandisoOS10looseZ_branch) {pfcandisoOS10looseZ_branch->SetAddress(&pfcandisoOS10looseZ_);}
	}
	pfcandptOS10looseZ_branch = 0;
	if (tree->GetBranch("pfcandptOS10looseZ") != 0) {
		pfcandptOS10looseZ_branch = tree->GetBranch("pfcandptOS10looseZ");
		if (pfcandptOS10looseZ_branch) {pfcandptOS10looseZ_branch->SetAddress(&pfcandptOS10looseZ_);}
	}
	pfcanddzOS10looseZ_branch = 0;
	if (tree->GetBranch("pfcanddzOS10looseZ") != 0) {
		pfcanddzOS10looseZ_branch = tree->GetBranch("pfcanddzOS10looseZ");
		if (pfcanddzOS10looseZ_branch) {pfcanddzOS10looseZ_branch->SetAddress(&pfcanddzOS10looseZ_);}
	}
	pfcanddirid10_branch = 0;
	if (tree->GetBranch("pfcanddirid10") != 0) {
		pfcanddirid10_branch = tree->GetBranch("pfcanddirid10");
		if (pfcanddirid10_branch) {pfcanddirid10_branch->SetAddress(&pfcanddirid10_);}
	}
	pfcanddiriso10_branch = 0;
	if (tree->GetBranch("pfcanddiriso10") != 0) {
		pfcanddiriso10_branch = tree->GetBranch("pfcanddiriso10");
		if (pfcanddiriso10_branch) {pfcanddiriso10_branch->SetAddress(&pfcanddiriso10_);}
	}
	pfcanddirpt10_branch = 0;
	if (tree->GetBranch("pfcanddirpt10") != 0) {
		pfcanddirpt10_branch = tree->GetBranch("pfcanddirpt10");
		if (pfcanddirpt10_branch) {pfcanddirpt10_branch->SetAddress(&pfcanddirpt10_);}
	}
	pfcanddirmindrj10_branch = 0;
	if (tree->GetBranch("pfcanddirmindrj10") != 0) {
		pfcanddirmindrj10_branch = tree->GetBranch("pfcanddirmindrj10");
		if (pfcanddirmindrj10_branch) {pfcanddirmindrj10_branch->SetAddress(&pfcanddirmindrj10_);}
	}
	pfcandvetoid10_branch = 0;
	if (tree->GetBranch("pfcandvetoid10") != 0) {
		pfcandvetoid10_branch = tree->GetBranch("pfcandvetoid10");
		if (pfcandvetoid10_branch) {pfcandvetoid10_branch->SetAddress(&pfcandvetoid10_);}
	}
	pfcandvetoiso10_branch = 0;
	if (tree->GetBranch("pfcandvetoiso10") != 0) {
		pfcandvetoiso10_branch = tree->GetBranch("pfcandvetoiso10");
		if (pfcandvetoiso10_branch) {pfcandvetoiso10_branch->SetAddress(&pfcandvetoiso10_);}
	}
	pfcandvetopt10_branch = 0;
	if (tree->GetBranch("pfcandvetopt10") != 0) {
		pfcandvetopt10_branch = tree->GetBranch("pfcandvetopt10");
		if (pfcandvetopt10_branch) {pfcandvetopt10_branch->SetAddress(&pfcandvetopt10_);}
	}
	pfcandvetomindrj10_branch = 0;
	if (tree->GetBranch("pfcandvetomindrj10") != 0) {
		pfcandvetomindrj10_branch = tree->GetBranch("pfcandvetomindrj10");
		if (pfcandvetomindrj10_branch) {pfcandvetomindrj10_branch->SetAddress(&pfcandvetomindrj10_);}
	}
	pfcandvetoLid10_branch = 0;
	if (tree->GetBranch("pfcandvetoLid10") != 0) {
		pfcandvetoLid10_branch = tree->GetBranch("pfcandvetoLid10");
		if (pfcandvetoLid10_branch) {pfcandvetoLid10_branch->SetAddress(&pfcandvetoLid10_);}
	}
	pfcandvetoLiso10_branch = 0;
	if (tree->GetBranch("pfcandvetoLiso10") != 0) {
		pfcandvetoLiso10_branch = tree->GetBranch("pfcandvetoLiso10");
		if (pfcandvetoLiso10_branch) {pfcandvetoLiso10_branch->SetAddress(&pfcandvetoLiso10_);}
	}
	pfcandvetoLpt10_branch = 0;
	if (tree->GetBranch("pfcandvetoLpt10") != 0) {
		pfcandvetoLpt10_branch = tree->GetBranch("pfcandvetoLpt10");
		if (pfcandvetoLpt10_branch) {pfcandvetoLpt10_branch->SetAddress(&pfcandvetoLpt10_);}
	}
	pfcandvetoLmindrj10_branch = 0;
	if (tree->GetBranch("pfcandvetoLmindrj10") != 0) {
		pfcandvetoLmindrj10_branch = tree->GetBranch("pfcandvetoLmindrj10");
		if (pfcandvetoLmindrj10_branch) {pfcandvetoLmindrj10_branch->SetAddress(&pfcandvetoLmindrj10_);}
	}
	emjet10_branch = 0;
	if (tree->GetBranch("emjet10") != 0) {
		emjet10_branch = tree->GetBranch("emjet10");
		if (emjet10_branch) {emjet10_branch->SetAddress(&emjet10_);}
	}
	mjj_branch = 0;
	if (tree->GetBranch("mjj") != 0) {
		mjj_branch = tree->GetBranch("mjj");
		if (mjj_branch) {mjj_branch->SetAddress(&mjj_);}
	}
	emjet20_branch = 0;
	if (tree->GetBranch("emjet20") != 0) {
		emjet20_branch = tree->GetBranch("emjet20");
		if (emjet20_branch) {emjet20_branch->SetAddress(&emjet20_);}
	}
	trkpt5_branch = 0;
	if (tree->GetBranch("trkpt5") != 0) {
		trkpt5_branch = tree->GetBranch("trkpt5");
		if (trkpt5_branch) {trkpt5_branch->SetAddress(&trkpt5_);}
	}
	trkpt10_branch = 0;
	if (tree->GetBranch("trkpt10") != 0) {
		trkpt10_branch = tree->GetBranch("trkpt10");
		if (trkpt10_branch) {trkpt10_branch->SetAddress(&trkpt10_);}
	}
	mleptrk5_branch = 0;
	if (tree->GetBranch("mleptrk5") != 0) {
		mleptrk5_branch = tree->GetBranch("mleptrk5");
		if (mleptrk5_branch) {mleptrk5_branch->SetAddress(&mleptrk5_);}
	}
	mleptrk10_branch = 0;
	if (tree->GetBranch("mleptrk10") != 0) {
		mleptrk10_branch = tree->GetBranch("mleptrk10");
		if (mleptrk10_branch) {mleptrk10_branch->SetAddress(&mleptrk10_);}
	}
	trkreliso5_branch = 0;
	if (tree->GetBranch("trkreliso5") != 0) {
		trkreliso5_branch = tree->GetBranch("trkreliso5");
		if (trkreliso5_branch) {trkreliso5_branch->SetAddress(&trkreliso5_);}
	}
	trkreliso10_branch = 0;
	if (tree->GetBranch("trkreliso10") != 0) {
		trkreliso10_branch = tree->GetBranch("trkreliso10");
		if (trkreliso10_branch) {trkreliso10_branch->SetAddress(&trkreliso10_);}
	}
	trkpt5loose_branch = 0;
	if (tree->GetBranch("trkpt5loose") != 0) {
		trkpt5loose_branch = tree->GetBranch("trkpt5loose");
		if (trkpt5loose_branch) {trkpt5loose_branch->SetAddress(&trkpt5loose_);}
	}
	trkpt10loose_branch = 0;
	if (tree->GetBranch("trkpt10loose") != 0) {
		trkpt10loose_branch = tree->GetBranch("trkpt10loose");
		if (trkpt10loose_branch) {trkpt10loose_branch->SetAddress(&trkpt10loose_);}
	}
	trkreliso5loose_branch = 0;
	if (tree->GetBranch("trkreliso5loose") != 0) {
		trkreliso5loose_branch = tree->GetBranch("trkreliso5loose");
		if (trkreliso5loose_branch) {trkreliso5loose_branch->SetAddress(&trkreliso5loose_);}
	}
	trkreliso10loose_branch = 0;
	if (tree->GetBranch("trkreliso10loose") != 0) {
		trkreliso10loose_branch = tree->GetBranch("trkreliso10loose");
		if (trkreliso10loose_branch) {trkreliso10loose_branch->SetAddress(&trkreliso10loose_);}
	}
	trkpt10pt0p1_branch = 0;
	if (tree->GetBranch("trkpt10pt0p1") != 0) {
		trkpt10pt0p1_branch = tree->GetBranch("trkpt10pt0p1");
		if (trkpt10pt0p1_branch) {trkpt10pt0p1_branch->SetAddress(&trkpt10pt0p1_);}
	}
	trkpt10pt0p2_branch = 0;
	if (tree->GetBranch("trkpt10pt0p2") != 0) {
		trkpt10pt0p2_branch = tree->GetBranch("trkpt10pt0p2");
		if (trkpt10pt0p2_branch) {trkpt10pt0p2_branch->SetAddress(&trkpt10pt0p2_);}
	}
	trkpt10pt0p3_branch = 0;
	if (tree->GetBranch("trkpt10pt0p3") != 0) {
		trkpt10pt0p3_branch = tree->GetBranch("trkpt10pt0p3");
		if (trkpt10pt0p3_branch) {trkpt10pt0p3_branch->SetAddress(&trkpt10pt0p3_);}
	}
	trkpt10pt0p4_branch = 0;
	if (tree->GetBranch("trkpt10pt0p4") != 0) {
		trkpt10pt0p4_branch = tree->GetBranch("trkpt10pt0p4");
		if (trkpt10pt0p4_branch) {trkpt10pt0p4_branch->SetAddress(&trkpt10pt0p4_);}
	}
	trkpt10pt0p5_branch = 0;
	if (tree->GetBranch("trkpt10pt0p5") != 0) {
		trkpt10pt0p5_branch = tree->GetBranch("trkpt10pt0p5");
		if (trkpt10pt0p5_branch) {trkpt10pt0p5_branch->SetAddress(&trkpt10pt0p5_);}
	}
	trkpt10pt0p6_branch = 0;
	if (tree->GetBranch("trkpt10pt0p6") != 0) {
		trkpt10pt0p6_branch = tree->GetBranch("trkpt10pt0p6");
		if (trkpt10pt0p6_branch) {trkpt10pt0p6_branch->SetAddress(&trkpt10pt0p6_);}
	}
	trkpt10pt0p7_branch = 0;
	if (tree->GetBranch("trkpt10pt0p7") != 0) {
		trkpt10pt0p7_branch = tree->GetBranch("trkpt10pt0p7");
		if (trkpt10pt0p7_branch) {trkpt10pt0p7_branch->SetAddress(&trkpt10pt0p7_);}
	}
	trkpt10pt0p8_branch = 0;
	if (tree->GetBranch("trkpt10pt0p8") != 0) {
		trkpt10pt0p8_branch = tree->GetBranch("trkpt10pt0p8");
		if (trkpt10pt0p8_branch) {trkpt10pt0p8_branch->SetAddress(&trkpt10pt0p8_);}
	}
	trkpt10pt0p9_branch = 0;
	if (tree->GetBranch("trkpt10pt0p9") != 0) {
		trkpt10pt0p9_branch = tree->GetBranch("trkpt10pt0p9");
		if (trkpt10pt0p9_branch) {trkpt10pt0p9_branch->SetAddress(&trkpt10pt0p9_);}
	}
	trkpt10pt1p0_branch = 0;
	if (tree->GetBranch("trkpt10pt1p0") != 0) {
		trkpt10pt1p0_branch = tree->GetBranch("trkpt10pt1p0");
		if (trkpt10pt1p0_branch) {trkpt10pt1p0_branch->SetAddress(&trkpt10pt1p0_);}
	}
	trkreliso10pt0p1_branch = 0;
	if (tree->GetBranch("trkreliso10pt0p1") != 0) {
		trkreliso10pt0p1_branch = tree->GetBranch("trkreliso10pt0p1");
		if (trkreliso10pt0p1_branch) {trkreliso10pt0p1_branch->SetAddress(&trkreliso10pt0p1_);}
	}
	trkreliso10pt0p2_branch = 0;
	if (tree->GetBranch("trkreliso10pt0p2") != 0) {
		trkreliso10pt0p2_branch = tree->GetBranch("trkreliso10pt0p2");
		if (trkreliso10pt0p2_branch) {trkreliso10pt0p2_branch->SetAddress(&trkreliso10pt0p2_);}
	}
	trkreliso10pt0p3_branch = 0;
	if (tree->GetBranch("trkreliso10pt0p3") != 0) {
		trkreliso10pt0p3_branch = tree->GetBranch("trkreliso10pt0p3");
		if (trkreliso10pt0p3_branch) {trkreliso10pt0p3_branch->SetAddress(&trkreliso10pt0p3_);}
	}
	trkreliso10pt0p4_branch = 0;
	if (tree->GetBranch("trkreliso10pt0p4") != 0) {
		trkreliso10pt0p4_branch = tree->GetBranch("trkreliso10pt0p4");
		if (trkreliso10pt0p4_branch) {trkreliso10pt0p4_branch->SetAddress(&trkreliso10pt0p4_);}
	}
	trkreliso10pt0p5_branch = 0;
	if (tree->GetBranch("trkreliso10pt0p5") != 0) {
		trkreliso10pt0p5_branch = tree->GetBranch("trkreliso10pt0p5");
		if (trkreliso10pt0p5_branch) {trkreliso10pt0p5_branch->SetAddress(&trkreliso10pt0p5_);}
	}
	trkreliso10pt0p6_branch = 0;
	if (tree->GetBranch("trkreliso10pt0p6") != 0) {
		trkreliso10pt0p6_branch = tree->GetBranch("trkreliso10pt0p6");
		if (trkreliso10pt0p6_branch) {trkreliso10pt0p6_branch->SetAddress(&trkreliso10pt0p6_);}
	}
	trkreliso10pt0p7_branch = 0;
	if (tree->GetBranch("trkreliso10pt0p7") != 0) {
		trkreliso10pt0p7_branch = tree->GetBranch("trkreliso10pt0p7");
		if (trkreliso10pt0p7_branch) {trkreliso10pt0p7_branch->SetAddress(&trkreliso10pt0p7_);}
	}
	trkreliso10pt0p8_branch = 0;
	if (tree->GetBranch("trkreliso10pt0p8") != 0) {
		trkreliso10pt0p8_branch = tree->GetBranch("trkreliso10pt0p8");
		if (trkreliso10pt0p8_branch) {trkreliso10pt0p8_branch->SetAddress(&trkreliso10pt0p8_);}
	}
	trkreliso10pt0p9_branch = 0;
	if (tree->GetBranch("trkreliso10pt0p9") != 0) {
		trkreliso10pt0p9_branch = tree->GetBranch("trkreliso10pt0p9");
		if (trkreliso10pt0p9_branch) {trkreliso10pt0p9_branch->SetAddress(&trkreliso10pt0p9_);}
	}
	trkreliso10pt1p0_branch = 0;
	if (tree->GetBranch("trkreliso10pt1p0") != 0) {
		trkreliso10pt1p0_branch = tree->GetBranch("trkreliso10pt1p0");
		if (trkreliso10pt1p0_branch) {trkreliso10pt1p0_branch->SetAddress(&trkreliso10pt1p0_);}
	}
	pfcandpt10pt0p1_branch = 0;
	if (tree->GetBranch("pfcandpt10pt0p1") != 0) {
		pfcandpt10pt0p1_branch = tree->GetBranch("pfcandpt10pt0p1");
		if (pfcandpt10pt0p1_branch) {pfcandpt10pt0p1_branch->SetAddress(&pfcandpt10pt0p1_);}
	}
	pfcandpt10pt0p2_branch = 0;
	if (tree->GetBranch("pfcandpt10pt0p2") != 0) {
		pfcandpt10pt0p2_branch = tree->GetBranch("pfcandpt10pt0p2");
		if (pfcandpt10pt0p2_branch) {pfcandpt10pt0p2_branch->SetAddress(&pfcandpt10pt0p2_);}
	}
	pfcandpt10pt0p3_branch = 0;
	if (tree->GetBranch("pfcandpt10pt0p3") != 0) {
		pfcandpt10pt0p3_branch = tree->GetBranch("pfcandpt10pt0p3");
		if (pfcandpt10pt0p3_branch) {pfcandpt10pt0p3_branch->SetAddress(&pfcandpt10pt0p3_);}
	}
	pfcandpt10pt0p4_branch = 0;
	if (tree->GetBranch("pfcandpt10pt0p4") != 0) {
		pfcandpt10pt0p4_branch = tree->GetBranch("pfcandpt10pt0p4");
		if (pfcandpt10pt0p4_branch) {pfcandpt10pt0p4_branch->SetAddress(&pfcandpt10pt0p4_);}
	}
	pfcandpt10pt0p5_branch = 0;
	if (tree->GetBranch("pfcandpt10pt0p5") != 0) {
		pfcandpt10pt0p5_branch = tree->GetBranch("pfcandpt10pt0p5");
		if (pfcandpt10pt0p5_branch) {pfcandpt10pt0p5_branch->SetAddress(&pfcandpt10pt0p5_);}
	}
	pfcandpt10pt0p6_branch = 0;
	if (tree->GetBranch("pfcandpt10pt0p6") != 0) {
		pfcandpt10pt0p6_branch = tree->GetBranch("pfcandpt10pt0p6");
		if (pfcandpt10pt0p6_branch) {pfcandpt10pt0p6_branch->SetAddress(&pfcandpt10pt0p6_);}
	}
	pfcandpt10pt0p7_branch = 0;
	if (tree->GetBranch("pfcandpt10pt0p7") != 0) {
		pfcandpt10pt0p7_branch = tree->GetBranch("pfcandpt10pt0p7");
		if (pfcandpt10pt0p7_branch) {pfcandpt10pt0p7_branch->SetAddress(&pfcandpt10pt0p7_);}
	}
	pfcandpt10pt0p8_branch = 0;
	if (tree->GetBranch("pfcandpt10pt0p8") != 0) {
		pfcandpt10pt0p8_branch = tree->GetBranch("pfcandpt10pt0p8");
		if (pfcandpt10pt0p8_branch) {pfcandpt10pt0p8_branch->SetAddress(&pfcandpt10pt0p8_);}
	}
	pfcandpt10pt0p9_branch = 0;
	if (tree->GetBranch("pfcandpt10pt0p9") != 0) {
		pfcandpt10pt0p9_branch = tree->GetBranch("pfcandpt10pt0p9");
		if (pfcandpt10pt0p9_branch) {pfcandpt10pt0p9_branch->SetAddress(&pfcandpt10pt0p9_);}
	}
	pfcandpt10pt1p0_branch = 0;
	if (tree->GetBranch("pfcandpt10pt1p0") != 0) {
		pfcandpt10pt1p0_branch = tree->GetBranch("pfcandpt10pt1p0");
		if (pfcandpt10pt1p0_branch) {pfcandpt10pt1p0_branch->SetAddress(&pfcandpt10pt1p0_);}
	}
	pfcandiso10pt0p1_branch = 0;
	if (tree->GetBranch("pfcandiso10pt0p1") != 0) {
		pfcandiso10pt0p1_branch = tree->GetBranch("pfcandiso10pt0p1");
		if (pfcandiso10pt0p1_branch) {pfcandiso10pt0p1_branch->SetAddress(&pfcandiso10pt0p1_);}
	}
	pfcandiso10pt0p2_branch = 0;
	if (tree->GetBranch("pfcandiso10pt0p2") != 0) {
		pfcandiso10pt0p2_branch = tree->GetBranch("pfcandiso10pt0p2");
		if (pfcandiso10pt0p2_branch) {pfcandiso10pt0p2_branch->SetAddress(&pfcandiso10pt0p2_);}
	}
	pfcandiso10pt0p3_branch = 0;
	if (tree->GetBranch("pfcandiso10pt0p3") != 0) {
		pfcandiso10pt0p3_branch = tree->GetBranch("pfcandiso10pt0p3");
		if (pfcandiso10pt0p3_branch) {pfcandiso10pt0p3_branch->SetAddress(&pfcandiso10pt0p3_);}
	}
	pfcandiso10pt0p4_branch = 0;
	if (tree->GetBranch("pfcandiso10pt0p4") != 0) {
		pfcandiso10pt0p4_branch = tree->GetBranch("pfcandiso10pt0p4");
		if (pfcandiso10pt0p4_branch) {pfcandiso10pt0p4_branch->SetAddress(&pfcandiso10pt0p4_);}
	}
	pfcandiso10pt0p5_branch = 0;
	if (tree->GetBranch("pfcandiso10pt0p5") != 0) {
		pfcandiso10pt0p5_branch = tree->GetBranch("pfcandiso10pt0p5");
		if (pfcandiso10pt0p5_branch) {pfcandiso10pt0p5_branch->SetAddress(&pfcandiso10pt0p5_);}
	}
	pfcandiso10pt0p6_branch = 0;
	if (tree->GetBranch("pfcandiso10pt0p6") != 0) {
		pfcandiso10pt0p6_branch = tree->GetBranch("pfcandiso10pt0p6");
		if (pfcandiso10pt0p6_branch) {pfcandiso10pt0p6_branch->SetAddress(&pfcandiso10pt0p6_);}
	}
	pfcandiso10pt0p7_branch = 0;
	if (tree->GetBranch("pfcandiso10pt0p7") != 0) {
		pfcandiso10pt0p7_branch = tree->GetBranch("pfcandiso10pt0p7");
		if (pfcandiso10pt0p7_branch) {pfcandiso10pt0p7_branch->SetAddress(&pfcandiso10pt0p7_);}
	}
	pfcandiso10pt0p8_branch = 0;
	if (tree->GetBranch("pfcandiso10pt0p8") != 0) {
		pfcandiso10pt0p8_branch = tree->GetBranch("pfcandiso10pt0p8");
		if (pfcandiso10pt0p8_branch) {pfcandiso10pt0p8_branch->SetAddress(&pfcandiso10pt0p8_);}
	}
	pfcandiso10pt0p9_branch = 0;
	if (tree->GetBranch("pfcandiso10pt0p9") != 0) {
		pfcandiso10pt0p9_branch = tree->GetBranch("pfcandiso10pt0p9");
		if (pfcandiso10pt0p9_branch) {pfcandiso10pt0p9_branch->SetAddress(&pfcandiso10pt0p9_);}
	}
	pfcandiso10pt1p0_branch = 0;
	if (tree->GetBranch("pfcandiso10pt1p0") != 0) {
		pfcandiso10pt1p0_branch = tree->GetBranch("pfcandiso10pt1p0");
		if (pfcandiso10pt1p0_branch) {pfcandiso10pt1p0_branch->SetAddress(&pfcandiso10pt1p0_);}
	}
	mbb_branch = 0;
	if (tree->GetBranch("mbb") != 0) {
		mbb_branch = tree->GetBranch("mbb");
		if (mbb_branch) {mbb_branch->SetAddress(&mbb_);}
	}
	lep1pfjetdr_branch = 0;
	if (tree->GetBranch("lep1pfjetdr") != 0) {
		lep1pfjetdr_branch = tree->GetBranch("lep1pfjetdr");
		if (lep1pfjetdr_branch) {lep1pfjetdr_branch->SetAddress(&lep1pfjetdr_);}
	}
	lep2pfjetdr_branch = 0;
	if (tree->GetBranch("lep2pfjetdr") != 0) {
		lep2pfjetdr_branch = tree->GetBranch("lep2pfjetdr");
		if (lep2pfjetdr_branch) {lep2pfjetdr_branch->SetAddress(&lep2pfjetdr_);}
	}
	mcmln_branch = 0;
	if (tree->GetBranch("mcmln") != 0) {
		mcmln_branch = tree->GetBranch("mcmln");
		if (mcmln_branch) {mcmln_branch->SetAddress(&mcmln_);}
	}
	mcmtln_branch = 0;
	if (tree->GetBranch("mcmtln") != 0) {
		mcmtln_branch = tree->GetBranch("mcmtln");
		if (mcmtln_branch) {mcmtln_branch->SetAddress(&mcmtln_);}
	}
	pfTau15_leadPtcandID_branch = 0;
	if (tree->GetBranch("pfTau15_leadPtcandID") != 0) {
		pfTau15_leadPtcandID_branch = tree->GetBranch("pfTau15_leadPtcandID");
		if (pfTau15_leadPtcandID_branch) {pfTau15_leadPtcandID_branch->SetAddress(&pfTau15_leadPtcandID_);}
	}
	pfTau_leadPtcandID_branch = 0;
	if (tree->GetBranch("pfTau_leadPtcandID") != 0) {
		pfTau_leadPtcandID_branch = tree->GetBranch("pfTau_leadPtcandID");
		if (pfTau_leadPtcandID_branch) {pfTau_leadPtcandID_branch->SetAddress(&pfTau_leadPtcandID_);}
	}
	pfTauLoose_leadPtcandID_branch = 0;
	if (tree->GetBranch("pfTauLoose_leadPtcandID") != 0) {
		pfTauLoose_leadPtcandID_branch = tree->GetBranch("pfTauLoose_leadPtcandID");
		if (pfTauLoose_leadPtcandID_branch) {pfTauLoose_leadPtcandID_branch->SetAddress(&pfTauLoose_leadPtcandID_);}
	}
	lep_t_id_branch = 0;
	if (tree->GetBranch("lep_t_id") != 0) {
		lep_t_id_branch = tree->GetBranch("lep_t_id");
		if (lep_t_id_branch) {lep_t_id_branch->SetAddress(&lep_t_id_);}
	}
	lep_tbar_id_branch = 0;
	if (tree->GetBranch("lep_tbar_id") != 0) {
		lep_tbar_id_branch = tree->GetBranch("lep_tbar_id");
		if (lep_tbar_id_branch) {lep_tbar_id_branch->SetAddress(&lep_tbar_id_);}
	}
	pfjets_csv_branch = 0;
	if (tree->GetBranch("pfjets_csv") != 0) {
		pfjets_csv_branch = tree->GetBranch("pfjets_csv");
		if (pfjets_csv_branch) {pfjets_csv_branch->SetAddress(&pfjets_csv_);}
	}
	pfjets_chEfrac_branch = 0;
	if (tree->GetBranch("pfjets_chEfrac") != 0) {
		pfjets_chEfrac_branch = tree->GetBranch("pfjets_chEfrac");
		if (pfjets_chEfrac_branch) {pfjets_chEfrac_branch->SetAddress(&pfjets_chEfrac_);}
	}
	pfjets_chm_branch = 0;
	if (tree->GetBranch("pfjets_chm") != 0) {
		pfjets_chm_branch = tree->GetBranch("pfjets_chm");
		if (pfjets_chm_branch) {pfjets_chm_branch->SetAddress(&pfjets_chm_);}
	}
	pfjets_neu_branch = 0;
	if (tree->GetBranch("pfjets_neu") != 0) {
		pfjets_neu_branch = tree->GetBranch("pfjets_neu");
		if (pfjets_neu_branch) {pfjets_neu_branch->SetAddress(&pfjets_neu_);}
	}
	pfjets_l1corr_branch = 0;
	if (tree->GetBranch("pfjets_l1corr") != 0) {
		pfjets_l1corr_branch = tree->GetBranch("pfjets_l1corr");
		if (pfjets_l1corr_branch) {pfjets_l1corr_branch->SetAddress(&pfjets_l1corr_);}
	}
	pfjets_corr_branch = 0;
	if (tree->GetBranch("pfjets_corr") != 0) {
		pfjets_corr_branch = tree->GetBranch("pfjets_corr");
		if (pfjets_corr_branch) {pfjets_corr_branch->SetAddress(&pfjets_corr_);}
	}
	pfjets_mc3_branch = 0;
	if (tree->GetBranch("pfjets_mc3") != 0) {
		pfjets_mc3_branch = tree->GetBranch("pfjets_mc3");
		if (pfjets_mc3_branch) {pfjets_mc3_branch->SetAddress(&pfjets_mc3_);}
	}
	pfjets_mcflavorAlgo_branch = 0;
	if (tree->GetBranch("pfjets_mcflavorAlgo") != 0) {
		pfjets_mcflavorAlgo_branch = tree->GetBranch("pfjets_mcflavorAlgo");
		if (pfjets_mcflavorAlgo_branch) {pfjets_mcflavorAlgo_branch->SetAddress(&pfjets_mcflavorAlgo_);}
	}
	pfjets_mcflavorPhys_branch = 0;
	if (tree->GetBranch("pfjets_mcflavorPhys") != 0) {
		pfjets_mcflavorPhys_branch = tree->GetBranch("pfjets_mcflavorPhys");
		if (pfjets_mcflavorPhys_branch) {pfjets_mcflavorPhys_branch->SetAddress(&pfjets_mcflavorPhys_);}
	}
	pfjets_uncertainty_branch = 0;
	if (tree->GetBranch("pfjets_uncertainty") != 0) {
		pfjets_uncertainty_branch = tree->GetBranch("pfjets_uncertainty");
		if (pfjets_uncertainty_branch) {pfjets_uncertainty_branch->SetAddress(&pfjets_uncertainty_);}
	}
	pfjets_flav_branch = 0;
	if (tree->GetBranch("pfjets_flav") != 0) {
		pfjets_flav_branch = tree->GetBranch("pfjets_flav");
		if (pfjets_flav_branch) {pfjets_flav_branch->SetAddress(&pfjets_flav_);}
	}
	pfjets_lrm_branch = 0;
	if (tree->GetBranch("pfjets_lrm") != 0) {
		pfjets_lrm_branch = tree->GetBranch("pfjets_lrm");
		if (pfjets_lrm_branch) {pfjets_lrm_branch->SetAddress(&pfjets_lrm_);}
	}
	pfjets_lrm2_branch = 0;
	if (tree->GetBranch("pfjets_lrm2") != 0) {
		pfjets_lrm2_branch = tree->GetBranch("pfjets_lrm2");
		if (pfjets_lrm2_branch) {pfjets_lrm2_branch->SetAddress(&pfjets_lrm2_);}
	}
	pfjets_qgtag_branch = 0;
	if (tree->GetBranch("pfjets_qgtag") != 0) {
		pfjets_qgtag_branch = tree->GetBranch("pfjets_qgtag");
		if (pfjets_qgtag_branch) {pfjets_qgtag_branch->SetAddress(&pfjets_qgtag_);}
	}
	pfjets_genJetDr_branch = 0;
	if (tree->GetBranch("pfjets_genJetDr") != 0) {
		pfjets_genJetDr_branch = tree->GetBranch("pfjets_genJetDr");
		if (pfjets_genJetDr_branch) {pfjets_genJetDr_branch->SetAddress(&pfjets_genJetDr_);}
	}
	pfjets_sigma_branch = 0;
	if (tree->GetBranch("pfjets_sigma") != 0) {
		pfjets_sigma_branch = tree->GetBranch("pfjets_sigma");
		if (pfjets_sigma_branch) {pfjets_sigma_branch->SetAddress(&pfjets_sigma_);}
	}
	pfjets_lepjet_branch = 0;
	if (tree->GetBranch("pfjets_lepjet") != 0) {
		pfjets_lepjet_branch = tree->GetBranch("pfjets_lepjet");
		if (pfjets_lepjet_branch) {pfjets_lepjet_branch->SetAddress(&pfjets_lepjet_);}
	}
	pfjets_tobtecmult_branch = 0;
	if (tree->GetBranch("pfjets_tobtecmult") != 0) {
		pfjets_tobtecmult_branch = tree->GetBranch("pfjets_tobtecmult");
		if (pfjets_tobtecmult_branch) {pfjets_tobtecmult_branch->SetAddress(&pfjets_tobtecmult_);}
	}
	pfjets_tobtecfrac_branch = 0;
	if (tree->GetBranch("pfjets_tobtecfrac") != 0) {
		pfjets_tobtecfrac_branch = tree->GetBranch("pfjets_tobtecfrac");
		if (pfjets_tobtecfrac_branch) {pfjets_tobtecfrac_branch->SetAddress(&pfjets_tobtecfrac_);}
	}
	pfjets_beta_branch = 0;
	if (tree->GetBranch("pfjets_beta") != 0) {
		pfjets_beta_branch = tree->GetBranch("pfjets_beta");
		if (pfjets_beta_branch) {pfjets_beta_branch->SetAddress(&pfjets_beta_);}
	}
	pfjets_beta2_branch = 0;
	if (tree->GetBranch("pfjets_beta2") != 0) {
		pfjets_beta2_branch = tree->GetBranch("pfjets_beta2");
		if (pfjets_beta2_branch) {pfjets_beta2_branch->SetAddress(&pfjets_beta2_);}
	}
	pfjets_beta_0p1_branch = 0;
	if (tree->GetBranch("pfjets_beta_0p1") != 0) {
		pfjets_beta_0p1_branch = tree->GetBranch("pfjets_beta_0p1");
		if (pfjets_beta_0p1_branch) {pfjets_beta_0p1_branch->SetAddress(&pfjets_beta_0p1_);}
	}
	pfjets_beta_0p2_branch = 0;
	if (tree->GetBranch("pfjets_beta_0p2") != 0) {
		pfjets_beta_0p2_branch = tree->GetBranch("pfjets_beta_0p2");
		if (pfjets_beta_0p2_branch) {pfjets_beta_0p2_branch->SetAddress(&pfjets_beta_0p2_);}
	}
	pfjets_beta2_0p1_branch = 0;
	if (tree->GetBranch("pfjets_beta2_0p1") != 0) {
		pfjets_beta2_0p1_branch = tree->GetBranch("pfjets_beta2_0p1");
		if (pfjets_beta2_0p1_branch) {pfjets_beta2_0p1_branch->SetAddress(&pfjets_beta2_0p1_);}
	}
	pfjets_beta2_0p5_branch = 0;
	if (tree->GetBranch("pfjets_beta2_0p5") != 0) {
		pfjets_beta2_0p5_branch = tree->GetBranch("pfjets_beta2_0p5");
		if (pfjets_beta2_0p5_branch) {pfjets_beta2_0p5_branch->SetAddress(&pfjets_beta2_0p5_);}
	}
	pfjets_mvaPUid_branch = 0;
	if (tree->GetBranch("pfjets_mvaPUid") != 0) {
		pfjets_mvaPUid_branch = tree->GetBranch("pfjets_mvaPUid");
		if (pfjets_mvaPUid_branch) {pfjets_mvaPUid_branch->SetAddress(&pfjets_mvaPUid_);}
	}
	pfjets_mva5xPUid_branch = 0;
	if (tree->GetBranch("pfjets_mva5xPUid") != 0) {
		pfjets_mva5xPUid_branch = tree->GetBranch("pfjets_mva5xPUid");
		if (pfjets_mva5xPUid_branch) {pfjets_mva5xPUid_branch->SetAddress(&pfjets_mva5xPUid_);}
	}
	pfjets_mvaBeta_branch = 0;
	if (tree->GetBranch("pfjets_mvaBeta") != 0) {
		pfjets_mvaBeta_branch = tree->GetBranch("pfjets_mvaBeta");
		if (pfjets_mvaBeta_branch) {pfjets_mvaBeta_branch->SetAddress(&pfjets_mvaBeta_);}
	}
	genps_pdgId_branch = 0;
	if (tree->GetBranch("genps_pdgId") != 0) {
		genps_pdgId_branch = tree->GetBranch("genps_pdgId");
		if (genps_pdgId_branch) {genps_pdgId_branch->SetAddress(&genps_pdgId_);}
	}
	genps_firstMother_branch = 0;
	if (tree->GetBranch("genps_firstMother") != 0) {
		genps_firstMother_branch = tree->GetBranch("genps_firstMother");
		if (genps_firstMother_branch) {genps_firstMother_branch->SetAddress(&genps_firstMother_);}
	}
	genps_energy_branch = 0;
	if (tree->GetBranch("genps_energy") != 0) {
		genps_energy_branch = tree->GetBranch("genps_energy");
		if (genps_energy_branch) {genps_energy_branch->SetAddress(&genps_energy_);}
	}
	genps_pt_branch = 0;
	if (tree->GetBranch("genps_pt") != 0) {
		genps_pt_branch = tree->GetBranch("genps_pt");
		if (genps_pt_branch) {genps_pt_branch->SetAddress(&genps_pt_);}
	}
	genps_eta_branch = 0;
	if (tree->GetBranch("genps_eta") != 0) {
		genps_eta_branch = tree->GetBranch("genps_eta");
		if (genps_eta_branch) {genps_eta_branch->SetAddress(&genps_eta_);}
	}
	genps_phi_branch = 0;
	if (tree->GetBranch("genps_phi") != 0) {
		genps_phi_branch = tree->GetBranch("genps_phi");
		if (genps_phi_branch) {genps_phi_branch->SetAddress(&genps_phi_);}
	}
	genps_mass_branch = 0;
	if (tree->GetBranch("genps_mass") != 0) {
		genps_mass_branch = tree->GetBranch("genps_mass");
		if (genps_mass_branch) {genps_mass_branch->SetAddress(&genps_mass_);}
	}
  tree->SetMakeClass(0);
}
void GetEntry(unsigned int idx) 
	// this only marks branches as not loaded, saving a lot of time
	{
		index = idx;
		acc_2010_isLoaded = false;
		acc_highmet_isLoaded = false;
		acc_highht_isLoaded = false;
		eldup_isLoaded = false;
		csc_isLoaded = false;
		hbhe_isLoaded = false;
		hbhenew_isLoaded = false;
		hcallaser_isLoaded = false;
		ecaltp_isLoaded = false;
		trkfail_isLoaded = false;
		eebadsc_isLoaded = false;
		lep1_badecallaser_isLoaded = false;
		lep2_badecallaser_isLoaded = false;
		isdata_isLoaded = false;
		jetid_isLoaded = false;
		jetid30_isLoaded = false;
		json_isLoaded = false;
		htoffset_isLoaded = false;
		htuncor_isLoaded = false;
		ptt_isLoaded = false;
		pttbar_isLoaded = false;
		ptttbar_isLoaded = false;
		mttbar_isLoaded = false;
		npartons_isLoaded = false;
		nwzpartons_isLoaded = false;
		hyptype_isLoaded = false;
		maxpartonpt_isLoaded = false;
		etattbar_isLoaded = false;
		njetsoffset_isLoaded = false;
		njetsuncor_isLoaded = false;
		costhetaweight_isLoaded = false;
		weight_isLoaded = false;
		weightleft_isLoaded = false;
		weightright_isLoaded = false;
		mutrigweight_isLoaded = false;
		mutrigweight2_isLoaded = false;
		sltrigweight_isLoaded = false;
		dltrigweight_isLoaded = false;
		trgeff_isLoaded = false;
		pthat_isLoaded = false;
		qscale_isLoaded = false;
		mgcor_isLoaded = false;
		wflav_isLoaded = false;
		ksusy_isLoaded = false;
		ksusyup_isLoaded = false;
		ksusydn_isLoaded = false;
		xsecsusy_isLoaded = false;
		xsecsusy2_isLoaded = false;
		smeff_isLoaded = false;
		k_isLoaded = false;
		mllgen_isLoaded = false;
		ptwgen_isLoaded = false;
		ptzgen_isLoaded = false;
		nlep_isLoaded = false;
		nosel_isLoaded = false;
		ngoodlep_isLoaded = false;
		ngoodel_isLoaded = false;
		ngoodmu_isLoaded = false;
		mull_isLoaded = false;
		mult_isLoaded = false;
		mullgen_isLoaded = false;
		multgen_isLoaded = false;
		proc_isLoaded = false;
		leptype_isLoaded = false;
		topmass_isLoaded = false;
		dilmass_isLoaded = false;
		dilrecoil_isLoaded = false;
		dilrecoilparl_isLoaded = false;
		dilrecoilperp_isLoaded = false;
		tcmet_isLoaded = false;
		genmet_isLoaded = false;
		gensumet_isLoaded = false;
		genmetphi_isLoaded = false;
		calomet_isLoaded = false;
		calometphi_isLoaded = false;
		trkmet_isLoaded = false;
		trkmetphi_isLoaded = false;
		pfmet_isLoaded = false;
		pfmetveto_isLoaded = false;
		pfmetsig_isLoaded = false;
		pfmetphi_isLoaded = false;
		pfsumet_isLoaded = false;
		mucormet_isLoaded = false;
		mucorjesmet_isLoaded = false;
		tcmet35X_isLoaded = false;
		tcmetevent_isLoaded = false;
		tcmetlooper_isLoaded = false;
		tcmetphi_isLoaded = false;
		tcsumet_isLoaded = false;
		tcmetUp_isLoaded = false;
		tcmetDown_isLoaded = false;
		tcmetTest_isLoaded = false;
		pfmetUp_isLoaded = false;
		pfmetDown_isLoaded = false;
		pfmetTest_isLoaded = false;
		sumjetpt_isLoaded = false;
		dileta_isLoaded = false;
		dilpt_isLoaded = false;
		dildphi_isLoaded = false;
		ngenjets_isLoaded = false;
		njpt_isLoaded = false;
		trgmu1_isLoaded = false;
		trgmu2_isLoaded = false;
		trgel1_isLoaded = false;
		trgel2_isLoaded = false;
		isomu24_isLoaded = false;
		ele27wp80_isLoaded = false;
		mm_isLoaded = false;
		mmtk_isLoaded = false;
		me_isLoaded = false;
		em_isLoaded = false;
		mu_isLoaded = false;
		ee_isLoaded = false;
		npfjets30_isLoaded = false;
		npfjets30lepcorr_isLoaded = false;
		knjets_isLoaded = false;
		rhovor_isLoaded = false;
		htpf30_isLoaded = false;
		t1met10_isLoaded = false;
		t1met20_isLoaded = false;
		t1met30_isLoaded = false;
		t1met10phi_isLoaded = false;
		t1met20phi_isLoaded = false;
		t1met30phi_isLoaded = false;
		t1met10mt_isLoaded = false;
		t1met20mt_isLoaded = false;
		t1met30mt_isLoaded = false;
		lepmetpt_isLoaded = false;
		lept1met10pt_isLoaded = false;
		t1met10s_isLoaded = false;
		t1met10sphi_isLoaded = false;
		t1met10smt_isLoaded = false;
		t1metphicorr_isLoaded = false;
		t1metphicorrup_isLoaded = false;
		t1metphicorrdn_isLoaded = false;
		t1metphicorrphi_isLoaded = false;
		t1metphicorrphiup_isLoaded = false;
		t1metphicorrphidn_isLoaded = false;
		t1metphicorrlep_isLoaded = false;
		t1metphicorrlepphi_isLoaded = false;
		t1metphicorrmt_isLoaded = false;
		t1metphicorrmtup_isLoaded = false;
		t1metphicorrmtdn_isLoaded = false;
		t1metphicorrlepmt_isLoaded = false;
		t1met_off_isLoaded = false;
		t1metphi_off_isLoaded = false;
		t1metmt_off_isLoaded = false;
		t1metphicorr_off_isLoaded = false;
		t1metphicorrphi_off_isLoaded = false;
		t1metphicorrmt_off_isLoaded = false;
		mht15_isLoaded = false;
		mht15phi_isLoaded = false;
		trkmet_mht15_isLoaded = false;
		trkmetphi_mht15_isLoaded = false;
		mettlj15_isLoaded = false;
		mettlj15phi_isLoaded = false;
		mt2bmin_isLoaded = false;
		mt2blmin_isLoaded = false;
		mt2wmin_isLoaded = false;
		chi2min_isLoaded = false;
		chi2minprob_isLoaded = false;
		nbtagsssv_isLoaded = false;
		nbtagstcl_isLoaded = false;
		nbtagstcm_isLoaded = false;
		nbtagscsvl_isLoaded = false;
		nbtagscsvm_isLoaded = false;
		nbtagscsvt_isLoaded = false;
		nbtagsssvcorr_isLoaded = false;
		nbtagstclcorr_isLoaded = false;
		nbtagstcmcorr_isLoaded = false;
		nbtagscsvlcorr_isLoaded = false;
		nbtagscsvmcorr_isLoaded = false;
		nbtagscsvtcott_isLoaded = false;
		njetsUp_isLoaded = false;
		njetsDown_isLoaded = false;
		htUp_isLoaded = false;
		htDown_isLoaded = false;
		ntruepu_isLoaded = false;
		npu_isLoaded = false;
		npuMinusOne_isLoaded = false;
		npuPlusOne_isLoaded = false;
		nvtx_isLoaded = false;
		indexfirstGoodVertex__isLoaded = false;
		nvtxweight_isLoaded = false;
		n3dvtxweight_isLoaded = false;
		pdfid1_isLoaded = false;
		pdfid2_isLoaded = false;
		pdfx1_isLoaded = false;
		pdfx2_isLoaded = false;
		pdfQ_isLoaded = false;
		vecjetpt_isLoaded = false;
		pass_isLoaded = false;
		passz_isLoaded = false;
		m0_isLoaded = false;
		mg_isLoaded = false;
		ml_isLoaded = false;
		x_isLoaded = false;
		m12_isLoaded = false;
		lep1chi2ndf_isLoaded = false;
		lep2chi2ndf_isLoaded = false;
		lep1dpt_isLoaded = false;
		lep2dpt_isLoaded = false;
		id1_isLoaded = false;
		id2_isLoaded = false;
		leptype1_isLoaded = false;
		leptype2_isLoaded = false;
		w1_isLoaded = false;
		w2_isLoaded = false;
		iso1_isLoaded = false;
		isont1_isLoaded = false;
		isopfold1_isLoaded = false;
		isopf1_isLoaded = false;
		etasc1_isLoaded = false;
		etasc2_isLoaded = false;
		eoverpin_isLoaded = false;
		eoverpout_isLoaded = false;
		dEtaIn_isLoaded = false;
		dPhiIn_isLoaded = false;
		sigmaIEtaIEta_isLoaded = false;
		hOverE_isLoaded = false;
		ooemoop_isLoaded = false;
		d0vtx_isLoaded = false;
		dzvtx_isLoaded = false;
		expinnerlayers_isLoaded = false;
		fbrem_isLoaded = false;
		pfisoch_isLoaded = false;
		pfisoem_isLoaded = false;
		pfisonh_isLoaded = false;
		eSC_isLoaded = false;
		phiSC_isLoaded = false;
		eSCRaw_isLoaded = false;
		eSCPresh_isLoaded = false;
		lep1_scslasercormean_isLoaded = false;
		lep1_scslasercormax_isLoaded = false;
		eoverpin2_isLoaded = false;
		eoverpout2_isLoaded = false;
		dEtaIn2_isLoaded = false;
		dPhiIn2_isLoaded = false;
		sigmaIEtaIEta2_isLoaded = false;
		hOverE2_isLoaded = false;
		ooemoop2_isLoaded = false;
		d0vtx2_isLoaded = false;
		dzvtx2_isLoaded = false;
		expinnerlayers2_isLoaded = false;
		fbrem2_isLoaded = false;
		pfisoch2_isLoaded = false;
		pfisoem2_isLoaded = false;
		pfisonh2_isLoaded = false;
		eSC2_isLoaded = false;
		phiSC2_isLoaded = false;
		eSCRaw2_isLoaded = false;
		eSCPresh2_isLoaded = false;
		lep2_scslasercormean_isLoaded = false;
		lep2_scslasercormax_isLoaded = false;
		scslasercormax_isLoaded = false;
		scslasercormax_pt_isLoaded = false;
		scslasercormax_eta_isLoaded = false;
		iso2_isLoaded = false;
		ecalveto1_isLoaded = false;
		ecalveto2_isLoaded = false;
		hcalveto1_isLoaded = false;
		hcalveto2_isLoaded = false;
		isont2_isLoaded = false;
		isopf2_isLoaded = false;
		ptl1_isLoaded = false;
		ptl2_isLoaded = false;
		etal1_isLoaded = false;
		etal2_isLoaded = false;
		phil1_isLoaded = false;
		phil2_isLoaded = false;
		meff_isLoaded = false;
		mt_isLoaded = false;
		run_isLoaded = false;
		lumi_isLoaded = false;
		event_isLoaded = false;
		y_isLoaded = false;
		ht_isLoaded = false;
		htgen_isLoaded = false;
		htjpt_isLoaded = false;
		nels_isLoaded = false;
		nmus_isLoaded = false;
		ntaus_isLoaded = false;
		nleps_isLoaded = false;
		nbs_isLoaded = false;
		dphijm_isLoaded = false;
		ptjetraw_isLoaded = false;
		ptjet23_isLoaded = false;
		ptjetF23_isLoaded = false;
		ptjetO23_isLoaded = false;
		mcid1_isLoaded = false;
		mcdr1_isLoaded = false;
		mcdecay1_isLoaded = false;
		mcndec1_isLoaded = false;
		mcndec2_isLoaded = false;
		mcndeckls1_isLoaded = false;
		mcndeckls2_isLoaded = false;
		mcndecem1_isLoaded = false;
		mcndecem2_isLoaded = false;
		mcid2_isLoaded = false;
		mcdr2_isLoaded = false;
		mcdecay2_isLoaded = false;
		mctaudpt1_isLoaded = false;
		mctaudpt2_isLoaded = false;
		mctaudid1_isLoaded = false;
		mctaudid2_isLoaded = false;
		mlepid_isLoaded = false;
		mleppassid_isLoaded = false;
		mleppassiso_isLoaded = false;
		mlepiso_isLoaded = false;
		mlepdr_isLoaded = false;
		pflepiso_isLoaded = false;
		pflepdr_isLoaded = false;
		pfleppt_isLoaded = false;
		pflepmindrj_isLoaded = false;
		pftaudiso_isLoaded = false;
		pftauddr_isLoaded = false;
		pftaudpt_isLoaded = false;
		pftaudmindrj_isLoaded = false;
		pfcandid5_isLoaded = false;
		pfcandiso5_isLoaded = false;
		pfcandpt5_isLoaded = false;
		pfcanddz5_isLoaded = false;
		pfcandmindrj5_isLoaded = false;
		pfcandid10_isLoaded = false;
		pfcandiso10_isLoaded = false;
		pfcandpt10_isLoaded = false;
		pfcanddz10_isLoaded = false;
		pfcandmindrj10_isLoaded = false;
		pfcandidOS10_isLoaded = false;
		pfcandisoOS10_isLoaded = false;
		pfcandptOS10_isLoaded = false;
		pfcanddzOS10_isLoaded = false;
		pfcandid5looseZ_isLoaded = false;
		pfcandiso5looseZ_isLoaded = false;
		pfcandpt5looseZ_isLoaded = false;
		pfcanddz5looseZ_isLoaded = false;
		pfcandidOS10looseZ_isLoaded = false;
		pfcandisoOS10looseZ_isLoaded = false;
		pfcandptOS10looseZ_isLoaded = false;
		pfcanddzOS10looseZ_isLoaded = false;
		pfcanddirid10_isLoaded = false;
		pfcanddiriso10_isLoaded = false;
		pfcanddirpt10_isLoaded = false;
		pfcanddirmindrj10_isLoaded = false;
		pfcandvetoid10_isLoaded = false;
		pfcandvetoiso10_isLoaded = false;
		pfcandvetopt10_isLoaded = false;
		pfcandvetomindrj10_isLoaded = false;
		pfcandvetoLid10_isLoaded = false;
		pfcandvetoLiso10_isLoaded = false;
		pfcandvetoLpt10_isLoaded = false;
		pfcandvetoLmindrj10_isLoaded = false;
		emjet10_isLoaded = false;
		mjj_isLoaded = false;
		emjet20_isLoaded = false;
		trkpt5_isLoaded = false;
		trkpt10_isLoaded = false;
		mleptrk5_isLoaded = false;
		mleptrk10_isLoaded = false;
		trkreliso5_isLoaded = false;
		trkreliso10_isLoaded = false;
		trkpt5loose_isLoaded = false;
		trkpt10loose_isLoaded = false;
		trkreliso5loose_isLoaded = false;
		trkreliso10loose_isLoaded = false;
		trkpt10pt0p1_isLoaded = false;
		trkpt10pt0p2_isLoaded = false;
		trkpt10pt0p3_isLoaded = false;
		trkpt10pt0p4_isLoaded = false;
		trkpt10pt0p5_isLoaded = false;
		trkpt10pt0p6_isLoaded = false;
		trkpt10pt0p7_isLoaded = false;
		trkpt10pt0p8_isLoaded = false;
		trkpt10pt0p9_isLoaded = false;
		trkpt10pt1p0_isLoaded = false;
		trkreliso10pt0p1_isLoaded = false;
		trkreliso10pt0p2_isLoaded = false;
		trkreliso10pt0p3_isLoaded = false;
		trkreliso10pt0p4_isLoaded = false;
		trkreliso10pt0p5_isLoaded = false;
		trkreliso10pt0p6_isLoaded = false;
		trkreliso10pt0p7_isLoaded = false;
		trkreliso10pt0p8_isLoaded = false;
		trkreliso10pt0p9_isLoaded = false;
		trkreliso10pt1p0_isLoaded = false;
		pfcandpt10pt0p1_isLoaded = false;
		pfcandpt10pt0p2_isLoaded = false;
		pfcandpt10pt0p3_isLoaded = false;
		pfcandpt10pt0p4_isLoaded = false;
		pfcandpt10pt0p5_isLoaded = false;
		pfcandpt10pt0p6_isLoaded = false;
		pfcandpt10pt0p7_isLoaded = false;
		pfcandpt10pt0p8_isLoaded = false;
		pfcandpt10pt0p9_isLoaded = false;
		pfcandpt10pt1p0_isLoaded = false;
		pfcandiso10pt0p1_isLoaded = false;
		pfcandiso10pt0p2_isLoaded = false;
		pfcandiso10pt0p3_isLoaded = false;
		pfcandiso10pt0p4_isLoaded = false;
		pfcandiso10pt0p5_isLoaded = false;
		pfcandiso10pt0p6_isLoaded = false;
		pfcandiso10pt0p7_isLoaded = false;
		pfcandiso10pt0p8_isLoaded = false;
		pfcandiso10pt0p9_isLoaded = false;
		pfcandiso10pt1p0_isLoaded = false;
		mbb_isLoaded = false;
		lep1pfjetdr_isLoaded = false;
		lep2pfjetdr_isLoaded = false;
		mclep_isLoaded = false;
		mcnu_isLoaded = false;
		mcmln_isLoaded = false;
		mcmtln_isLoaded = false;
		mlep_isLoaded = false;
		lep1_isLoaded = false;
		lep2_isLoaded = false;
		trklep1_isLoaded = false;
		trklep2_isLoaded = false;
		gfitlep1_isLoaded = false;
		gfitlep2_isLoaded = false;
		lepp_isLoaded = false;
		lepm_isLoaded = false;
		pflep1_isLoaded = false;
		pflep2_isLoaded = false;
		leppfjet1_isLoaded = false;
		leppfjet2_isLoaded = false;
		mclep1_isLoaded = false;
		mclep2_isLoaded = false;
		mctaud1_isLoaded = false;
		mctaud2_isLoaded = false;
		mctaudvis1_isLoaded = false;
		mctaudvis2_isLoaded = false;
		pflep_isLoaded = false;
		pftaud_isLoaded = false;
		pfcand5_isLoaded = false;
		pfcand10_isLoaded = false;
		pfTau15_leadPtcandID_isLoaded = false;
		pfTau15_isLoaded = false;
		pfTau15_leadPtcand_isLoaded = false;
		pfTau_leadPtcandID_isLoaded = false;
		pfTau_isLoaded = false;
		pfTau_leadPtcand_isLoaded = false;
		pfTauLoose_leadPtcandID_isLoaded = false;
		pfTauLoose_isLoaded = false;
		pfTauLoose_leadPtcand_isLoaded = false;
		pfcandOS10_isLoaded = false;
		pfcandOS10looseZ_isLoaded = false;
		pfcand5looseZ_isLoaded = false;
		pfcanddir10_isLoaded = false;
		pfcandveto10_isLoaded = false;
		pfcandvetoL10_isLoaded = false;
		jet_isLoaded = false;
		nonisoel_isLoaded = false;
		nonisomu_isLoaded = false;
		t_isLoaded = false;
		tbar_isLoaded = false;
		ttbar_isLoaded = false;
		lep_t_isLoaded = false;
		lep_tbar_isLoaded = false;
		stop_t_isLoaded = false;
		stop_tbar_isLoaded = false;
		neutralino_t_isLoaded = false;
		neutralino_tbar_isLoaded = false;
		lep_t_id_isLoaded = false;
		lep_tbar_id_isLoaded = false;
		pfjets_isLoaded = false;
		pfjets_genJet__isLoaded = false;
		pfjets_failjetid_isLoaded = false;
		pfjets_faillepolap_isLoaded = false;
		pfjets_csv_isLoaded = false;
		pfjets_chEfrac_isLoaded = false;
		pfjets_chm_isLoaded = false;
		pfjets_neu_isLoaded = false;
		pfjets_l1corr_isLoaded = false;
		pfjets_corr_isLoaded = false;
		pfjets_mc3_isLoaded = false;
		pfjets_mcflavorAlgo_isLoaded = false;
		pfjets_mcflavorPhys_isLoaded = false;
		pfjets_uncertainty_isLoaded = false;
		pfjets_flav_isLoaded = false;
		pfjets_lrm_isLoaded = false;
		pfjets_lrm2_isLoaded = false;
		pfjets_qgtag_isLoaded = false;
		pfjets_genJetDr_isLoaded = false;
		pfjets_sigma_isLoaded = false;
		pfjets_lepjet_isLoaded = false;
		pfjets_tobtecmult_isLoaded = false;
		pfjets_tobtecfrac_isLoaded = false;
		pfjets_beta_isLoaded = false;
		pfjets_beta2_isLoaded = false;
		pfjets_beta_0p1_isLoaded = false;
		pfjets_beta_0p2_isLoaded = false;
		pfjets_beta2_0p1_isLoaded = false;
		pfjets_beta2_0p5_isLoaded = false;
		pfjets_mvaPUid_isLoaded = false;
		pfjets_mva5xPUid_isLoaded = false;
		pfjets_mvaBeta_isLoaded = false;
		genjets_isLoaded = false;
		genqgs_isLoaded = false;
		genbs_isLoaded = false;
		genps_pdgId_isLoaded = false;
		genps_firstMother_isLoaded = false;
		genps_energy_isLoaded = false;
		genps_pt_isLoaded = false;
		genps_eta_isLoaded = false;
		genps_phi_isLoaded = false;
		genps_mass_isLoaded = false;
	}

void LoadAllBranches() 
	// load all branches
{
	if (acc_2010_branch != 0) acc_2010();
	if (acc_highmet_branch != 0) acc_highmet();
	if (acc_highht_branch != 0) acc_highht();
	if (eldup_branch != 0) eldup();
	if (csc_branch != 0) csc();
	if (hbhe_branch != 0) hbhe();
	if (hbhenew_branch != 0) hbhenew();
	if (hcallaser_branch != 0) hcallaser();
	if (ecaltp_branch != 0) ecaltp();
	if (trkfail_branch != 0) trkfail();
	if (eebadsc_branch != 0) eebadsc();
	if (lep1_badecallaser_branch != 0) lep1_badecallaser();
	if (lep2_badecallaser_branch != 0) lep2_badecallaser();
	if (isdata_branch != 0) isdata();
	if (jetid_branch != 0) jetid();
	if (jetid30_branch != 0) jetid30();
	if (json_branch != 0) json();
	if (htoffset_branch != 0) htoffset();
	if (htuncor_branch != 0) htuncor();
	if (ptt_branch != 0) ptt();
	if (pttbar_branch != 0) pttbar();
	if (ptttbar_branch != 0) ptttbar();
	if (mttbar_branch != 0) mttbar();
	if (npartons_branch != 0) npartons();
	if (nwzpartons_branch != 0) nwzpartons();
	if (hyptype_branch != 0) hyptype();
	if (maxpartonpt_branch != 0) maxpartonpt();
	if (etattbar_branch != 0) etattbar();
	if (njetsoffset_branch != 0) njetsoffset();
	if (njetsuncor_branch != 0) njetsuncor();
	if (costhetaweight_branch != 0) costhetaweight();
	if (weight_branch != 0) weight();
	if (weightleft_branch != 0) weightleft();
	if (weightright_branch != 0) weightright();
	if (mutrigweight_branch != 0) mutrigweight();
	if (mutrigweight2_branch != 0) mutrigweight2();
	if (sltrigweight_branch != 0) sltrigweight();
	if (dltrigweight_branch != 0) dltrigweight();
	if (trgeff_branch != 0) trgeff();
	if (pthat_branch != 0) pthat();
	if (qscale_branch != 0) qscale();
	if (mgcor_branch != 0) mgcor();
	if (wflav_branch != 0) wflav();
	if (ksusy_branch != 0) ksusy();
	if (ksusyup_branch != 0) ksusyup();
	if (ksusydn_branch != 0) ksusydn();
	if (xsecsusy_branch != 0) xsecsusy();
	if (xsecsusy2_branch != 0) xsecsusy2();
	if (smeff_branch != 0) smeff();
	if (k_branch != 0) k();
	if (mllgen_branch != 0) mllgen();
	if (ptwgen_branch != 0) ptwgen();
	if (ptzgen_branch != 0) ptzgen();
	if (nlep_branch != 0) nlep();
	if (nosel_branch != 0) nosel();
	if (ngoodlep_branch != 0) ngoodlep();
	if (ngoodel_branch != 0) ngoodel();
	if (ngoodmu_branch != 0) ngoodmu();
	if (mull_branch != 0) mull();
	if (mult_branch != 0) mult();
	if (mullgen_branch != 0) mullgen();
	if (multgen_branch != 0) multgen();
	if (proc_branch != 0) proc();
	if (leptype_branch != 0) leptype();
	if (topmass_branch != 0) topmass();
	if (dilmass_branch != 0) dilmass();
	if (dilrecoil_branch != 0) dilrecoil();
	if (dilrecoilparl_branch != 0) dilrecoilparl();
	if (dilrecoilperp_branch != 0) dilrecoilperp();
	if (tcmet_branch != 0) tcmet();
	if (genmet_branch != 0) genmet();
	if (gensumet_branch != 0) gensumet();
	if (genmetphi_branch != 0) genmetphi();
	if (calomet_branch != 0) calomet();
	if (calometphi_branch != 0) calometphi();
	if (trkmet_branch != 0) trkmet();
	if (trkmetphi_branch != 0) trkmetphi();
	if (pfmet_branch != 0) pfmet();
	if (pfmetveto_branch != 0) pfmetveto();
	if (pfmetsig_branch != 0) pfmetsig();
	if (pfmetphi_branch != 0) pfmetphi();
	if (pfsumet_branch != 0) pfsumet();
	if (mucormet_branch != 0) mucormet();
	if (mucorjesmet_branch != 0) mucorjesmet();
	if (tcmet35X_branch != 0) tcmet35X();
	if (tcmetevent_branch != 0) tcmetevent();
	if (tcmetlooper_branch != 0) tcmetlooper();
	if (tcmetphi_branch != 0) tcmetphi();
	if (tcsumet_branch != 0) tcsumet();
	if (tcmetUp_branch != 0) tcmetUp();
	if (tcmetDown_branch != 0) tcmetDown();
	if (tcmetTest_branch != 0) tcmetTest();
	if (pfmetUp_branch != 0) pfmetUp();
	if (pfmetDown_branch != 0) pfmetDown();
	if (pfmetTest_branch != 0) pfmetTest();
	if (sumjetpt_branch != 0) sumjetpt();
	if (dileta_branch != 0) dileta();
	if (dilpt_branch != 0) dilpt();
	if (dildphi_branch != 0) dildphi();
	if (ngenjets_branch != 0) ngenjets();
	if (njpt_branch != 0) njpt();
	if (trgmu1_branch != 0) trgmu1();
	if (trgmu2_branch != 0) trgmu2();
	if (trgel1_branch != 0) trgel1();
	if (trgel2_branch != 0) trgel2();
	if (isomu24_branch != 0) isomu24();
	if (ele27wp80_branch != 0) ele27wp80();
	if (mm_branch != 0) mm();
	if (mmtk_branch != 0) mmtk();
	if (me_branch != 0) me();
	if (em_branch != 0) em();
	if (mu_branch != 0) mu();
	if (ee_branch != 0) ee();
	if (npfjets30_branch != 0) npfjets30();
	if (npfjets30lepcorr_branch != 0) npfjets30lepcorr();
	if (knjets_branch != 0) knjets();
	if (rhovor_branch != 0) rhovor();
	if (htpf30_branch != 0) htpf30();
	if (t1met10_branch != 0) t1met10();
	if (t1met20_branch != 0) t1met20();
	if (t1met30_branch != 0) t1met30();
	if (t1met10phi_branch != 0) t1met10phi();
	if (t1met20phi_branch != 0) t1met20phi();
	if (t1met30phi_branch != 0) t1met30phi();
	if (t1met10mt_branch != 0) t1met10mt();
	if (t1met20mt_branch != 0) t1met20mt();
	if (t1met30mt_branch != 0) t1met30mt();
	if (lepmetpt_branch != 0) lepmetpt();
	if (lept1met10pt_branch != 0) lept1met10pt();
	if (t1met10s_branch != 0) t1met10s();
	if (t1met10sphi_branch != 0) t1met10sphi();
	if (t1met10smt_branch != 0) t1met10smt();
	if (t1metphicorr_branch != 0) t1metphicorr();
	if (t1metphicorrup_branch != 0) t1metphicorrup();
	if (t1metphicorrdn_branch != 0) t1metphicorrdn();
	if (t1metphicorrphi_branch != 0) t1metphicorrphi();
	if (t1metphicorrphiup_branch != 0) t1metphicorrphiup();
	if (t1metphicorrphidn_branch != 0) t1metphicorrphidn();
	if (t1metphicorrlep_branch != 0) t1metphicorrlep();
	if (t1metphicorrlepphi_branch != 0) t1metphicorrlepphi();
	if (t1metphicorrmt_branch != 0) t1metphicorrmt();
	if (t1metphicorrmtup_branch != 0) t1metphicorrmtup();
	if (t1metphicorrmtdn_branch != 0) t1metphicorrmtdn();
	if (t1metphicorrlepmt_branch != 0) t1metphicorrlepmt();
	if (t1met_off_branch != 0) t1met_off();
	if (t1metphi_off_branch != 0) t1metphi_off();
	if (t1metmt_off_branch != 0) t1metmt_off();
	if (t1metphicorr_off_branch != 0) t1metphicorr_off();
	if (t1metphicorrphi_off_branch != 0) t1metphicorrphi_off();
	if (t1metphicorrmt_off_branch != 0) t1metphicorrmt_off();
	if (mht15_branch != 0) mht15();
	if (mht15phi_branch != 0) mht15phi();
	if (trkmet_mht15_branch != 0) trkmet_mht15();
	if (trkmetphi_mht15_branch != 0) trkmetphi_mht15();
	if (mettlj15_branch != 0) mettlj15();
	if (mettlj15phi_branch != 0) mettlj15phi();
	if (mt2bmin_branch != 0) mt2bmin();
	if (mt2blmin_branch != 0) mt2blmin();
	if (mt2wmin_branch != 0) mt2wmin();
	if (chi2min_branch != 0) chi2min();
	if (chi2minprob_branch != 0) chi2minprob();
	if (nbtagsssv_branch != 0) nbtagsssv();
	if (nbtagstcl_branch != 0) nbtagstcl();
	if (nbtagstcm_branch != 0) nbtagstcm();
	if (nbtagscsvl_branch != 0) nbtagscsvl();
	if (nbtagscsvm_branch != 0) nbtagscsvm();
	if (nbtagscsvt_branch != 0) nbtagscsvt();
	if (nbtagsssvcorr_branch != 0) nbtagsssvcorr();
	if (nbtagstclcorr_branch != 0) nbtagstclcorr();
	if (nbtagstcmcorr_branch != 0) nbtagstcmcorr();
	if (nbtagscsvlcorr_branch != 0) nbtagscsvlcorr();
	if (nbtagscsvmcorr_branch != 0) nbtagscsvmcorr();
	if (nbtagscsvtcott_branch != 0) nbtagscsvtcott();
	if (njetsUp_branch != 0) njetsUp();
	if (njetsDown_branch != 0) njetsDown();
	if (htUp_branch != 0) htUp();
	if (htDown_branch != 0) htDown();
	if (ntruepu_branch != 0) ntruepu();
	if (npu_branch != 0) npu();
	if (npuMinusOne_branch != 0) npuMinusOne();
	if (npuPlusOne_branch != 0) npuPlusOne();
	if (nvtx_branch != 0) nvtx();
	if (indexfirstGoodVertex__branch != 0) indexfirstGoodVertex_();
	if (nvtxweight_branch != 0) nvtxweight();
	if (n3dvtxweight_branch != 0) n3dvtxweight();
	if (pdfid1_branch != 0) pdfid1();
	if (pdfid2_branch != 0) pdfid2();
	if (pdfx1_branch != 0) pdfx1();
	if (pdfx2_branch != 0) pdfx2();
	if (pdfQ_branch != 0) pdfQ();
	if (vecjetpt_branch != 0) vecjetpt();
	if (pass_branch != 0) pass();
	if (passz_branch != 0) passz();
	if (m0_branch != 0) m0();
	if (mg_branch != 0) mg();
	if (ml_branch != 0) ml();
	if (x_branch != 0) x();
	if (m12_branch != 0) m12();
	if (lep1chi2ndf_branch != 0) lep1chi2ndf();
	if (lep2chi2ndf_branch != 0) lep2chi2ndf();
	if (lep1dpt_branch != 0) lep1dpt();
	if (lep2dpt_branch != 0) lep2dpt();
	if (id1_branch != 0) id1();
	if (id2_branch != 0) id2();
	if (leptype1_branch != 0) leptype1();
	if (leptype2_branch != 0) leptype2();
	if (w1_branch != 0) w1();
	if (w2_branch != 0) w2();
	if (iso1_branch != 0) iso1();
	if (isont1_branch != 0) isont1();
	if (isopfold1_branch != 0) isopfold1();
	if (isopf1_branch != 0) isopf1();
	if (etasc1_branch != 0) etasc1();
	if (etasc2_branch != 0) etasc2();
	if (eoverpin_branch != 0) eoverpin();
	if (eoverpout_branch != 0) eoverpout();
	if (dEtaIn_branch != 0) dEtaIn();
	if (dPhiIn_branch != 0) dPhiIn();
	if (sigmaIEtaIEta_branch != 0) sigmaIEtaIEta();
	if (hOverE_branch != 0) hOverE();
	if (ooemoop_branch != 0) ooemoop();
	if (d0vtx_branch != 0) d0vtx();
	if (dzvtx_branch != 0) dzvtx();
	if (expinnerlayers_branch != 0) expinnerlayers();
	if (fbrem_branch != 0) fbrem();
	if (pfisoch_branch != 0) pfisoch();
	if (pfisoem_branch != 0) pfisoem();
	if (pfisonh_branch != 0) pfisonh();
	if (eSC_branch != 0) eSC();
	if (phiSC_branch != 0) phiSC();
	if (eSCRaw_branch != 0) eSCRaw();
	if (eSCPresh_branch != 0) eSCPresh();
	if (lep1_scslasercormean_branch != 0) lep1_scslasercormean();
	if (lep1_scslasercormax_branch != 0) lep1_scslasercormax();
	if (eoverpin2_branch != 0) eoverpin2();
	if (eoverpout2_branch != 0) eoverpout2();
	if (dEtaIn2_branch != 0) dEtaIn2();
	if (dPhiIn2_branch != 0) dPhiIn2();
	if (sigmaIEtaIEta2_branch != 0) sigmaIEtaIEta2();
	if (hOverE2_branch != 0) hOverE2();
	if (ooemoop2_branch != 0) ooemoop2();
	if (d0vtx2_branch != 0) d0vtx2();
	if (dzvtx2_branch != 0) dzvtx2();
	if (expinnerlayers2_branch != 0) expinnerlayers2();
	if (fbrem2_branch != 0) fbrem2();
	if (pfisoch2_branch != 0) pfisoch2();
	if (pfisoem2_branch != 0) pfisoem2();
	if (pfisonh2_branch != 0) pfisonh2();
	if (eSC2_branch != 0) eSC2();
	if (phiSC2_branch != 0) phiSC2();
	if (eSCRaw2_branch != 0) eSCRaw2();
	if (eSCPresh2_branch != 0) eSCPresh2();
	if (lep2_scslasercormean_branch != 0) lep2_scslasercormean();
	if (lep2_scslasercormax_branch != 0) lep2_scslasercormax();
	if (scslasercormax_branch != 0) scslasercormax();
	if (scslasercormax_pt_branch != 0) scslasercormax_pt();
	if (scslasercormax_eta_branch != 0) scslasercormax_eta();
	if (iso2_branch != 0) iso2();
	if (ecalveto1_branch != 0) ecalveto1();
	if (ecalveto2_branch != 0) ecalveto2();
	if (hcalveto1_branch != 0) hcalveto1();
	if (hcalveto2_branch != 0) hcalveto2();
	if (isont2_branch != 0) isont2();
	if (isopf2_branch != 0) isopf2();
	if (ptl1_branch != 0) ptl1();
	if (ptl2_branch != 0) ptl2();
	if (etal1_branch != 0) etal1();
	if (etal2_branch != 0) etal2();
	if (phil1_branch != 0) phil1();
	if (phil2_branch != 0) phil2();
	if (meff_branch != 0) meff();
	if (mt_branch != 0) mt();
	if (run_branch != 0) run();
	if (lumi_branch != 0) lumi();
	if (event_branch != 0) event();
	if (y_branch != 0) y();
	if (ht_branch != 0) ht();
	if (htgen_branch != 0) htgen();
	if (htjpt_branch != 0) htjpt();
	if (nels_branch != 0) nels();
	if (nmus_branch != 0) nmus();
	if (ntaus_branch != 0) ntaus();
	if (nleps_branch != 0) nleps();
	if (nbs_branch != 0) nbs();
	if (dphijm_branch != 0) dphijm();
	if (ptjetraw_branch != 0) ptjetraw();
	if (ptjet23_branch != 0) ptjet23();
	if (ptjetF23_branch != 0) ptjetF23();
	if (ptjetO23_branch != 0) ptjetO23();
	if (mcid1_branch != 0) mcid1();
	if (mcdr1_branch != 0) mcdr1();
	if (mcdecay1_branch != 0) mcdecay1();
	if (mcndec1_branch != 0) mcndec1();
	if (mcndec2_branch != 0) mcndec2();
	if (mcndeckls1_branch != 0) mcndeckls1();
	if (mcndeckls2_branch != 0) mcndeckls2();
	if (mcndecem1_branch != 0) mcndecem1();
	if (mcndecem2_branch != 0) mcndecem2();
	if (mcid2_branch != 0) mcid2();
	if (mcdr2_branch != 0) mcdr2();
	if (mcdecay2_branch != 0) mcdecay2();
	if (mctaudpt1_branch != 0) mctaudpt1();
	if (mctaudpt2_branch != 0) mctaudpt2();
	if (mctaudid1_branch != 0) mctaudid1();
	if (mctaudid2_branch != 0) mctaudid2();
	if (mlepid_branch != 0) mlepid();
	if (mleppassid_branch != 0) mleppassid();
	if (mleppassiso_branch != 0) mleppassiso();
	if (mlepiso_branch != 0) mlepiso();
	if (mlepdr_branch != 0) mlepdr();
	if (pflepiso_branch != 0) pflepiso();
	if (pflepdr_branch != 0) pflepdr();
	if (pfleppt_branch != 0) pfleppt();
	if (pflepmindrj_branch != 0) pflepmindrj();
	if (pftaudiso_branch != 0) pftaudiso();
	if (pftauddr_branch != 0) pftauddr();
	if (pftaudpt_branch != 0) pftaudpt();
	if (pftaudmindrj_branch != 0) pftaudmindrj();
	if (pfcandid5_branch != 0) pfcandid5();
	if (pfcandiso5_branch != 0) pfcandiso5();
	if (pfcandpt5_branch != 0) pfcandpt5();
	if (pfcanddz5_branch != 0) pfcanddz5();
	if (pfcandmindrj5_branch != 0) pfcandmindrj5();
	if (pfcandid10_branch != 0) pfcandid10();
	if (pfcandiso10_branch != 0) pfcandiso10();
	if (pfcandpt10_branch != 0) pfcandpt10();
	if (pfcanddz10_branch != 0) pfcanddz10();
	if (pfcandmindrj10_branch != 0) pfcandmindrj10();
	if (pfcandidOS10_branch != 0) pfcandidOS10();
	if (pfcandisoOS10_branch != 0) pfcandisoOS10();
	if (pfcandptOS10_branch != 0) pfcandptOS10();
	if (pfcanddzOS10_branch != 0) pfcanddzOS10();
	if (pfcandid5looseZ_branch != 0) pfcandid5looseZ();
	if (pfcandiso5looseZ_branch != 0) pfcandiso5looseZ();
	if (pfcandpt5looseZ_branch != 0) pfcandpt5looseZ();
	if (pfcanddz5looseZ_branch != 0) pfcanddz5looseZ();
	if (pfcandidOS10looseZ_branch != 0) pfcandidOS10looseZ();
	if (pfcandisoOS10looseZ_branch != 0) pfcandisoOS10looseZ();
	if (pfcandptOS10looseZ_branch != 0) pfcandptOS10looseZ();
	if (pfcanddzOS10looseZ_branch != 0) pfcanddzOS10looseZ();
	if (pfcanddirid10_branch != 0) pfcanddirid10();
	if (pfcanddiriso10_branch != 0) pfcanddiriso10();
	if (pfcanddirpt10_branch != 0) pfcanddirpt10();
	if (pfcanddirmindrj10_branch != 0) pfcanddirmindrj10();
	if (pfcandvetoid10_branch != 0) pfcandvetoid10();
	if (pfcandvetoiso10_branch != 0) pfcandvetoiso10();
	if (pfcandvetopt10_branch != 0) pfcandvetopt10();
	if (pfcandvetomindrj10_branch != 0) pfcandvetomindrj10();
	if (pfcandvetoLid10_branch != 0) pfcandvetoLid10();
	if (pfcandvetoLiso10_branch != 0) pfcandvetoLiso10();
	if (pfcandvetoLpt10_branch != 0) pfcandvetoLpt10();
	if (pfcandvetoLmindrj10_branch != 0) pfcandvetoLmindrj10();
	if (emjet10_branch != 0) emjet10();
	if (mjj_branch != 0) mjj();
	if (emjet20_branch != 0) emjet20();
	if (trkpt5_branch != 0) trkpt5();
	if (trkpt10_branch != 0) trkpt10();
	if (mleptrk5_branch != 0) mleptrk5();
	if (mleptrk10_branch != 0) mleptrk10();
	if (trkreliso5_branch != 0) trkreliso5();
	if (trkreliso10_branch != 0) trkreliso10();
	if (trkpt5loose_branch != 0) trkpt5loose();
	if (trkpt10loose_branch != 0) trkpt10loose();
	if (trkreliso5loose_branch != 0) trkreliso5loose();
	if (trkreliso10loose_branch != 0) trkreliso10loose();
	if (trkpt10pt0p1_branch != 0) trkpt10pt0p1();
	if (trkpt10pt0p2_branch != 0) trkpt10pt0p2();
	if (trkpt10pt0p3_branch != 0) trkpt10pt0p3();
	if (trkpt10pt0p4_branch != 0) trkpt10pt0p4();
	if (trkpt10pt0p5_branch != 0) trkpt10pt0p5();
	if (trkpt10pt0p6_branch != 0) trkpt10pt0p6();
	if (trkpt10pt0p7_branch != 0) trkpt10pt0p7();
	if (trkpt10pt0p8_branch != 0) trkpt10pt0p8();
	if (trkpt10pt0p9_branch != 0) trkpt10pt0p9();
	if (trkpt10pt1p0_branch != 0) trkpt10pt1p0();
	if (trkreliso10pt0p1_branch != 0) trkreliso10pt0p1();
	if (trkreliso10pt0p2_branch != 0) trkreliso10pt0p2();
	if (trkreliso10pt0p3_branch != 0) trkreliso10pt0p3();
	if (trkreliso10pt0p4_branch != 0) trkreliso10pt0p4();
	if (trkreliso10pt0p5_branch != 0) trkreliso10pt0p5();
	if (trkreliso10pt0p6_branch != 0) trkreliso10pt0p6();
	if (trkreliso10pt0p7_branch != 0) trkreliso10pt0p7();
	if (trkreliso10pt0p8_branch != 0) trkreliso10pt0p8();
	if (trkreliso10pt0p9_branch != 0) trkreliso10pt0p9();
	if (trkreliso10pt1p0_branch != 0) trkreliso10pt1p0();
	if (pfcandpt10pt0p1_branch != 0) pfcandpt10pt0p1();
	if (pfcandpt10pt0p2_branch != 0) pfcandpt10pt0p2();
	if (pfcandpt10pt0p3_branch != 0) pfcandpt10pt0p3();
	if (pfcandpt10pt0p4_branch != 0) pfcandpt10pt0p4();
	if (pfcandpt10pt0p5_branch != 0) pfcandpt10pt0p5();
	if (pfcandpt10pt0p6_branch != 0) pfcandpt10pt0p6();
	if (pfcandpt10pt0p7_branch != 0) pfcandpt10pt0p7();
	if (pfcandpt10pt0p8_branch != 0) pfcandpt10pt0p8();
	if (pfcandpt10pt0p9_branch != 0) pfcandpt10pt0p9();
	if (pfcandpt10pt1p0_branch != 0) pfcandpt10pt1p0();
	if (pfcandiso10pt0p1_branch != 0) pfcandiso10pt0p1();
	if (pfcandiso10pt0p2_branch != 0) pfcandiso10pt0p2();
	if (pfcandiso10pt0p3_branch != 0) pfcandiso10pt0p3();
	if (pfcandiso10pt0p4_branch != 0) pfcandiso10pt0p4();
	if (pfcandiso10pt0p5_branch != 0) pfcandiso10pt0p5();
	if (pfcandiso10pt0p6_branch != 0) pfcandiso10pt0p6();
	if (pfcandiso10pt0p7_branch != 0) pfcandiso10pt0p7();
	if (pfcandiso10pt0p8_branch != 0) pfcandiso10pt0p8();
	if (pfcandiso10pt0p9_branch != 0) pfcandiso10pt0p9();
	if (pfcandiso10pt1p0_branch != 0) pfcandiso10pt1p0();
	if (mbb_branch != 0) mbb();
	if (lep1pfjetdr_branch != 0) lep1pfjetdr();
	if (lep2pfjetdr_branch != 0) lep2pfjetdr();
	if (mclep_branch != 0) mclep();
	if (mcnu_branch != 0) mcnu();
	if (mcmln_branch != 0) mcmln();
	if (mcmtln_branch != 0) mcmtln();
	if (mlep_branch != 0) mlep();
	if (lep1_branch != 0) lep1();
	if (lep2_branch != 0) lep2();
	if (trklep1_branch != 0) trklep1();
	if (trklep2_branch != 0) trklep2();
	if (gfitlep1_branch != 0) gfitlep1();
	if (gfitlep2_branch != 0) gfitlep2();
	if (lepp_branch != 0) lepp();
	if (lepm_branch != 0) lepm();
	if (pflep1_branch != 0) pflep1();
	if (pflep2_branch != 0) pflep2();
	if (leppfjet1_branch != 0) leppfjet1();
	if (leppfjet2_branch != 0) leppfjet2();
	if (mclep1_branch != 0) mclep1();
	if (mclep2_branch != 0) mclep2();
	if (mctaud1_branch != 0) mctaud1();
	if (mctaud2_branch != 0) mctaud2();
	if (mctaudvis1_branch != 0) mctaudvis1();
	if (mctaudvis2_branch != 0) mctaudvis2();
	if (pflep_branch != 0) pflep();
	if (pftaud_branch != 0) pftaud();
	if (pfcand5_branch != 0) pfcand5();
	if (pfcand10_branch != 0) pfcand10();
	if (pfTau15_leadPtcandID_branch != 0) pfTau15_leadPtcandID();
	if (pfTau15_branch != 0) pfTau15();
	if (pfTau15_leadPtcand_branch != 0) pfTau15_leadPtcand();
	if (pfTau_leadPtcandID_branch != 0) pfTau_leadPtcandID();
	if (pfTau_branch != 0) pfTau();
	if (pfTau_leadPtcand_branch != 0) pfTau_leadPtcand();
	if (pfTauLoose_leadPtcandID_branch != 0) pfTauLoose_leadPtcandID();
	if (pfTauLoose_branch != 0) pfTauLoose();
	if (pfTauLoose_leadPtcand_branch != 0) pfTauLoose_leadPtcand();
	if (pfcandOS10_branch != 0) pfcandOS10();
	if (pfcandOS10looseZ_branch != 0) pfcandOS10looseZ();
	if (pfcand5looseZ_branch != 0) pfcand5looseZ();
	if (pfcanddir10_branch != 0) pfcanddir10();
	if (pfcandveto10_branch != 0) pfcandveto10();
	if (pfcandvetoL10_branch != 0) pfcandvetoL10();
	if (jet_branch != 0) jet();
	if (nonisoel_branch != 0) nonisoel();
	if (nonisomu_branch != 0) nonisomu();
	if (t_branch != 0) t();
	if (tbar_branch != 0) tbar();
	if (ttbar_branch != 0) ttbar();
	if (lep_t_branch != 0) lep_t();
	if (lep_tbar_branch != 0) lep_tbar();
	if (stop_t_branch != 0) stop_t();
	if (stop_tbar_branch != 0) stop_tbar();
	if (neutralino_t_branch != 0) neutralino_t();
	if (neutralino_tbar_branch != 0) neutralino_tbar();
	if (lep_t_id_branch != 0) lep_t_id();
	if (lep_tbar_id_branch != 0) lep_tbar_id();
	if (pfjets_branch != 0) pfjets();
	if (pfjets_genJet__branch != 0) pfjets_genJet_();
	if (pfjets_failjetid_branch != 0) pfjets_failjetid();
	if (pfjets_faillepolap_branch != 0) pfjets_faillepolap();
	if (pfjets_csv_branch != 0) pfjets_csv();
	if (pfjets_chEfrac_branch != 0) pfjets_chEfrac();
	if (pfjets_chm_branch != 0) pfjets_chm();
	if (pfjets_neu_branch != 0) pfjets_neu();
	if (pfjets_l1corr_branch != 0) pfjets_l1corr();
	if (pfjets_corr_branch != 0) pfjets_corr();
	if (pfjets_mc3_branch != 0) pfjets_mc3();
	if (pfjets_mcflavorAlgo_branch != 0) pfjets_mcflavorAlgo();
	if (pfjets_mcflavorPhys_branch != 0) pfjets_mcflavorPhys();
	if (pfjets_uncertainty_branch != 0) pfjets_uncertainty();
	if (pfjets_flav_branch != 0) pfjets_flav();
	if (pfjets_lrm_branch != 0) pfjets_lrm();
	if (pfjets_lrm2_branch != 0) pfjets_lrm2();
	if (pfjets_qgtag_branch != 0) pfjets_qgtag();
	if (pfjets_genJetDr_branch != 0) pfjets_genJetDr();
	if (pfjets_sigma_branch != 0) pfjets_sigma();
	if (pfjets_lepjet_branch != 0) pfjets_lepjet();
	if (pfjets_tobtecmult_branch != 0) pfjets_tobtecmult();
	if (pfjets_tobtecfrac_branch != 0) pfjets_tobtecfrac();
	if (pfjets_beta_branch != 0) pfjets_beta();
	if (pfjets_beta2_branch != 0) pfjets_beta2();
	if (pfjets_beta_0p1_branch != 0) pfjets_beta_0p1();
	if (pfjets_beta_0p2_branch != 0) pfjets_beta_0p2();
	if (pfjets_beta2_0p1_branch != 0) pfjets_beta2_0p1();
	if (pfjets_beta2_0p5_branch != 0) pfjets_beta2_0p5();
	if (pfjets_mvaPUid_branch != 0) pfjets_mvaPUid();
	if (pfjets_mva5xPUid_branch != 0) pfjets_mva5xPUid();
	if (pfjets_mvaBeta_branch != 0) pfjets_mvaBeta();
	if (genjets_branch != 0) genjets();
	if (genqgs_branch != 0) genqgs();
	if (genbs_branch != 0) genbs();
	if (genps_pdgId_branch != 0) genps_pdgId();
	if (genps_firstMother_branch != 0) genps_firstMother();
	if (genps_energy_branch != 0) genps_energy();
	if (genps_pt_branch != 0) genps_pt();
	if (genps_eta_branch != 0) genps_eta();
	if (genps_phi_branch != 0) genps_phi();
	if (genps_mass_branch != 0) genps_mass();
}

	int &acc_2010()
	{
		if (not acc_2010_isLoaded) {
			if (acc_2010_branch != 0) {
				acc_2010_branch->GetEntry(index);
			} else { 
				printf("branch acc_2010_branch does not exist!\n");
				exit(1);
			}
			acc_2010_isLoaded = true;
		}
		return acc_2010_;
	}
	int &acc_highmet()
	{
		if (not acc_highmet_isLoaded) {
			if (acc_highmet_branch != 0) {
				acc_highmet_branch->GetEntry(index);
			} else { 
				printf("branch acc_highmet_branch does not exist!\n");
				exit(1);
			}
			acc_highmet_isLoaded = true;
		}
		return acc_highmet_;
	}
	int &acc_highht()
	{
		if (not acc_highht_isLoaded) {
			if (acc_highht_branch != 0) {
				acc_highht_branch->GetEntry(index);
			} else { 
				printf("branch acc_highht_branch does not exist!\n");
				exit(1);
			}
			acc_highht_isLoaded = true;
		}
		return acc_highht_;
	}
	int &eldup()
	{
		if (not eldup_isLoaded) {
			if (eldup_branch != 0) {
				eldup_branch->GetEntry(index);
			} else { 
				printf("branch eldup_branch does not exist!\n");
				exit(1);
			}
			eldup_isLoaded = true;
		}
		return eldup_;
	}
	int &csc()
	{
		if (not csc_isLoaded) {
			if (csc_branch != 0) {
				csc_branch->GetEntry(index);
			} else { 
				printf("branch csc_branch does not exist!\n");
				exit(1);
			}
			csc_isLoaded = true;
		}
		return csc_;
	}
	int &hbhe()
	{
		if (not hbhe_isLoaded) {
			if (hbhe_branch != 0) {
				hbhe_branch->GetEntry(index);
			} else { 
				printf("branch hbhe_branch does not exist!\n");
				exit(1);
			}
			hbhe_isLoaded = true;
		}
		return hbhe_;
	}
	int &hbhenew()
	{
		if (not hbhenew_isLoaded) {
			if (hbhenew_branch != 0) {
				hbhenew_branch->GetEntry(index);
			} else { 
				printf("branch hbhenew_branch does not exist!\n");
				exit(1);
			}
			hbhenew_isLoaded = true;
		}
		return hbhenew_;
	}
	int &hcallaser()
	{
		if (not hcallaser_isLoaded) {
			if (hcallaser_branch != 0) {
				hcallaser_branch->GetEntry(index);
			} else { 
				printf("branch hcallaser_branch does not exist!\n");
				exit(1);
			}
			hcallaser_isLoaded = true;
		}
		return hcallaser_;
	}
	int &ecaltp()
	{
		if (not ecaltp_isLoaded) {
			if (ecaltp_branch != 0) {
				ecaltp_branch->GetEntry(index);
			} else { 
				printf("branch ecaltp_branch does not exist!\n");
				exit(1);
			}
			ecaltp_isLoaded = true;
		}
		return ecaltp_;
	}
	int &trkfail()
	{
		if (not trkfail_isLoaded) {
			if (trkfail_branch != 0) {
				trkfail_branch->GetEntry(index);
			} else { 
				printf("branch trkfail_branch does not exist!\n");
				exit(1);
			}
			trkfail_isLoaded = true;
		}
		return trkfail_;
	}
	int &eebadsc()
	{
		if (not eebadsc_isLoaded) {
			if (eebadsc_branch != 0) {
				eebadsc_branch->GetEntry(index);
			} else { 
				printf("branch eebadsc_branch does not exist!\n");
				exit(1);
			}
			eebadsc_isLoaded = true;
		}
		return eebadsc_;
	}
	int &lep1_badecallaser()
	{
		if (not lep1_badecallaser_isLoaded) {
			if (lep1_badecallaser_branch != 0) {
				lep1_badecallaser_branch->GetEntry(index);
			} else { 
				printf("branch lep1_badecallaser_branch does not exist!\n");
				exit(1);
			}
			lep1_badecallaser_isLoaded = true;
		}
		return lep1_badecallaser_;
	}
	int &lep2_badecallaser()
	{
		if (not lep2_badecallaser_isLoaded) {
			if (lep2_badecallaser_branch != 0) {
				lep2_badecallaser_branch->GetEntry(index);
			} else { 
				printf("branch lep2_badecallaser_branch does not exist!\n");
				exit(1);
			}
			lep2_badecallaser_isLoaded = true;
		}
		return lep2_badecallaser_;
	}
	int &isdata()
	{
		if (not isdata_isLoaded) {
			if (isdata_branch != 0) {
				isdata_branch->GetEntry(index);
			} else { 
				printf("branch isdata_branch does not exist!\n");
				exit(1);
			}
			isdata_isLoaded = true;
		}
		return isdata_;
	}
	int &jetid()
	{
		if (not jetid_isLoaded) {
			if (jetid_branch != 0) {
				jetid_branch->GetEntry(index);
			} else { 
				printf("branch jetid_branch does not exist!\n");
				exit(1);
			}
			jetid_isLoaded = true;
		}
		return jetid_;
	}
	int &jetid30()
	{
		if (not jetid30_isLoaded) {
			if (jetid30_branch != 0) {
				jetid30_branch->GetEntry(index);
			} else { 
				printf("branch jetid30_branch does not exist!\n");
				exit(1);
			}
			jetid30_isLoaded = true;
		}
		return jetid30_;
	}
	int &json()
	{
		if (not json_isLoaded) {
			if (json_branch != 0) {
				json_branch->GetEntry(index);
			} else { 
				printf("branch json_branch does not exist!\n");
				exit(1);
			}
			json_isLoaded = true;
		}
		return json_;
	}
	float &htoffset()
	{
		if (not htoffset_isLoaded) {
			if (htoffset_branch != 0) {
				htoffset_branch->GetEntry(index);
			} else { 
				printf("branch htoffset_branch does not exist!\n");
				exit(1);
			}
			htoffset_isLoaded = true;
		}
		return htoffset_;
	}
	float &htuncor()
	{
		if (not htuncor_isLoaded) {
			if (htuncor_branch != 0) {
				htuncor_branch->GetEntry(index);
			} else { 
				printf("branch htuncor_branch does not exist!\n");
				exit(1);
			}
			htuncor_isLoaded = true;
		}
		return htuncor_;
	}
	float &ptt()
	{
		if (not ptt_isLoaded) {
			if (ptt_branch != 0) {
				ptt_branch->GetEntry(index);
			} else { 
				printf("branch ptt_branch does not exist!\n");
				exit(1);
			}
			ptt_isLoaded = true;
		}
		return ptt_;
	}
	float &pttbar()
	{
		if (not pttbar_isLoaded) {
			if (pttbar_branch != 0) {
				pttbar_branch->GetEntry(index);
			} else { 
				printf("branch pttbar_branch does not exist!\n");
				exit(1);
			}
			pttbar_isLoaded = true;
		}
		return pttbar_;
	}
	float &ptttbar()
	{
		if (not ptttbar_isLoaded) {
			if (ptttbar_branch != 0) {
				ptttbar_branch->GetEntry(index);
			} else { 
				printf("branch ptttbar_branch does not exist!\n");
				exit(1);
			}
			ptttbar_isLoaded = true;
		}
		return ptttbar_;
	}
	float &mttbar()
	{
		if (not mttbar_isLoaded) {
			if (mttbar_branch != 0) {
				mttbar_branch->GetEntry(index);
			} else { 
				printf("branch mttbar_branch does not exist!\n");
				exit(1);
			}
			mttbar_isLoaded = true;
		}
		return mttbar_;
	}
	int &npartons()
	{
		if (not npartons_isLoaded) {
			if (npartons_branch != 0) {
				npartons_branch->GetEntry(index);
			} else { 
				printf("branch npartons_branch does not exist!\n");
				exit(1);
			}
			npartons_isLoaded = true;
		}
		return npartons_;
	}
	int &nwzpartons()
	{
		if (not nwzpartons_isLoaded) {
			if (nwzpartons_branch != 0) {
				nwzpartons_branch->GetEntry(index);
			} else { 
				printf("branch nwzpartons_branch does not exist!\n");
				exit(1);
			}
			nwzpartons_isLoaded = true;
		}
		return nwzpartons_;
	}
	int &hyptype()
	{
		if (not hyptype_isLoaded) {
			if (hyptype_branch != 0) {
				hyptype_branch->GetEntry(index);
			} else { 
				printf("branch hyptype_branch does not exist!\n");
				exit(1);
			}
			hyptype_isLoaded = true;
		}
		return hyptype_;
	}
	float &maxpartonpt()
	{
		if (not maxpartonpt_isLoaded) {
			if (maxpartonpt_branch != 0) {
				maxpartonpt_branch->GetEntry(index);
			} else { 
				printf("branch maxpartonpt_branch does not exist!\n");
				exit(1);
			}
			maxpartonpt_isLoaded = true;
		}
		return maxpartonpt_;
	}
	float &etattbar()
	{
		if (not etattbar_isLoaded) {
			if (etattbar_branch != 0) {
				etattbar_branch->GetEntry(index);
			} else { 
				printf("branch etattbar_branch does not exist!\n");
				exit(1);
			}
			etattbar_isLoaded = true;
		}
		return etattbar_;
	}
	int &njetsoffset()
	{
		if (not njetsoffset_isLoaded) {
			if (njetsoffset_branch != 0) {
				njetsoffset_branch->GetEntry(index);
			} else { 
				printf("branch njetsoffset_branch does not exist!\n");
				exit(1);
			}
			njetsoffset_isLoaded = true;
		}
		return njetsoffset_;
	}
	int &njetsuncor()
	{
		if (not njetsuncor_isLoaded) {
			if (njetsuncor_branch != 0) {
				njetsuncor_branch->GetEntry(index);
			} else { 
				printf("branch njetsuncor_branch does not exist!\n");
				exit(1);
			}
			njetsuncor_isLoaded = true;
		}
		return njetsuncor_;
	}
	float &costhetaweight()
	{
		if (not costhetaweight_isLoaded) {
			if (costhetaweight_branch != 0) {
				costhetaweight_branch->GetEntry(index);
			} else { 
				printf("branch costhetaweight_branch does not exist!\n");
				exit(1);
			}
			costhetaweight_isLoaded = true;
		}
		return costhetaweight_;
	}
	float &weight()
	{
		if (not weight_isLoaded) {
			if (weight_branch != 0) {
				weight_branch->GetEntry(index);
			} else { 
				printf("branch weight_branch does not exist!\n");
				exit(1);
			}
			weight_isLoaded = true;
		}
		return weight_;
	}
	float &weightleft()
	{
		if (not weightleft_isLoaded) {
			if (weightleft_branch != 0) {
				weightleft_branch->GetEntry(index);
			} else { 
				printf("branch weightleft_branch does not exist!\n");
				exit(1);
			}
			weightleft_isLoaded = true;
		}
		return weightleft_;
	}
	float &weightright()
	{
		if (not weightright_isLoaded) {
			if (weightright_branch != 0) {
				weightright_branch->GetEntry(index);
			} else { 
				printf("branch weightright_branch does not exist!\n");
				exit(1);
			}
			weightright_isLoaded = true;
		}
		return weightright_;
	}
	float &mutrigweight()
	{
		if (not mutrigweight_isLoaded) {
			if (mutrigweight_branch != 0) {
				mutrigweight_branch->GetEntry(index);
			} else { 
				printf("branch mutrigweight_branch does not exist!\n");
				exit(1);
			}
			mutrigweight_isLoaded = true;
		}
		return mutrigweight_;
	}
	float &mutrigweight2()
	{
		if (not mutrigweight2_isLoaded) {
			if (mutrigweight2_branch != 0) {
				mutrigweight2_branch->GetEntry(index);
			} else { 
				printf("branch mutrigweight2_branch does not exist!\n");
				exit(1);
			}
			mutrigweight2_isLoaded = true;
		}
		return mutrigweight2_;
	}
	float &sltrigweight()
	{
		if (not sltrigweight_isLoaded) {
			if (sltrigweight_branch != 0) {
				sltrigweight_branch->GetEntry(index);
			} else { 
				printf("branch sltrigweight_branch does not exist!\n");
				exit(1);
			}
			sltrigweight_isLoaded = true;
		}
		return sltrigweight_;
	}
	float &dltrigweight()
	{
		if (not dltrigweight_isLoaded) {
			if (dltrigweight_branch != 0) {
				dltrigweight_branch->GetEntry(index);
			} else { 
				printf("branch dltrigweight_branch does not exist!\n");
				exit(1);
			}
			dltrigweight_isLoaded = true;
		}
		return dltrigweight_;
	}
	float &trgeff()
	{
		if (not trgeff_isLoaded) {
			if (trgeff_branch != 0) {
				trgeff_branch->GetEntry(index);
			} else { 
				printf("branch trgeff_branch does not exist!\n");
				exit(1);
			}
			trgeff_isLoaded = true;
		}
		return trgeff_;
	}
	float &pthat()
	{
		if (not pthat_isLoaded) {
			if (pthat_branch != 0) {
				pthat_branch->GetEntry(index);
			} else { 
				printf("branch pthat_branch does not exist!\n");
				exit(1);
			}
			pthat_isLoaded = true;
		}
		return pthat_;
	}
	float &qscale()
	{
		if (not qscale_isLoaded) {
			if (qscale_branch != 0) {
				qscale_branch->GetEntry(index);
			} else { 
				printf("branch qscale_branch does not exist!\n");
				exit(1);
			}
			qscale_isLoaded = true;
		}
		return qscale_;
	}
	float &mgcor()
	{
		if (not mgcor_isLoaded) {
			if (mgcor_branch != 0) {
				mgcor_branch->GetEntry(index);
			} else { 
				printf("branch mgcor_branch does not exist!\n");
				exit(1);
			}
			mgcor_isLoaded = true;
		}
		return mgcor_;
	}
	int &wflav()
	{
		if (not wflav_isLoaded) {
			if (wflav_branch != 0) {
				wflav_branch->GetEntry(index);
			} else { 
				printf("branch wflav_branch does not exist!\n");
				exit(1);
			}
			wflav_isLoaded = true;
		}
		return wflav_;
	}
	float &ksusy()
	{
		if (not ksusy_isLoaded) {
			if (ksusy_branch != 0) {
				ksusy_branch->GetEntry(index);
			} else { 
				printf("branch ksusy_branch does not exist!\n");
				exit(1);
			}
			ksusy_isLoaded = true;
		}
		return ksusy_;
	}
	float &ksusyup()
	{
		if (not ksusyup_isLoaded) {
			if (ksusyup_branch != 0) {
				ksusyup_branch->GetEntry(index);
			} else { 
				printf("branch ksusyup_branch does not exist!\n");
				exit(1);
			}
			ksusyup_isLoaded = true;
		}
		return ksusyup_;
	}
	float &ksusydn()
	{
		if (not ksusydn_isLoaded) {
			if (ksusydn_branch != 0) {
				ksusydn_branch->GetEntry(index);
			} else { 
				printf("branch ksusydn_branch does not exist!\n");
				exit(1);
			}
			ksusydn_isLoaded = true;
		}
		return ksusydn_;
	}
	float &xsecsusy()
	{
		if (not xsecsusy_isLoaded) {
			if (xsecsusy_branch != 0) {
				xsecsusy_branch->GetEntry(index);
			} else { 
				printf("branch xsecsusy_branch does not exist!\n");
				exit(1);
			}
			xsecsusy_isLoaded = true;
		}
		return xsecsusy_;
	}
	float &xsecsusy2()
	{
		if (not xsecsusy2_isLoaded) {
			if (xsecsusy2_branch != 0) {
				xsecsusy2_branch->GetEntry(index);
			} else { 
				printf("branch xsecsusy2_branch does not exist!\n");
				exit(1);
			}
			xsecsusy2_isLoaded = true;
		}
		return xsecsusy2_;
	}
	float &smeff()
	{
		if (not smeff_isLoaded) {
			if (smeff_branch != 0) {
				smeff_branch->GetEntry(index);
			} else { 
				printf("branch smeff_branch does not exist!\n");
				exit(1);
			}
			smeff_isLoaded = true;
		}
		return smeff_;
	}
	float &k()
	{
		if (not k_isLoaded) {
			if (k_branch != 0) {
				k_branch->GetEntry(index);
			} else { 
				printf("branch k_branch does not exist!\n");
				exit(1);
			}
			k_isLoaded = true;
		}
		return k_;
	}
	float &mllgen()
	{
		if (not mllgen_isLoaded) {
			if (mllgen_branch != 0) {
				mllgen_branch->GetEntry(index);
			} else { 
				printf("branch mllgen_branch does not exist!\n");
				exit(1);
			}
			mllgen_isLoaded = true;
		}
		return mllgen_;
	}
	float &ptwgen()
	{
		if (not ptwgen_isLoaded) {
			if (ptwgen_branch != 0) {
				ptwgen_branch->GetEntry(index);
			} else { 
				printf("branch ptwgen_branch does not exist!\n");
				exit(1);
			}
			ptwgen_isLoaded = true;
		}
		return ptwgen_;
	}
	float &ptzgen()
	{
		if (not ptzgen_isLoaded) {
			if (ptzgen_branch != 0) {
				ptzgen_branch->GetEntry(index);
			} else { 
				printf("branch ptzgen_branch does not exist!\n");
				exit(1);
			}
			ptzgen_isLoaded = true;
		}
		return ptzgen_;
	}
	int &nlep()
	{
		if (not nlep_isLoaded) {
			if (nlep_branch != 0) {
				nlep_branch->GetEntry(index);
			} else { 
				printf("branch nlep_branch does not exist!\n");
				exit(1);
			}
			nlep_isLoaded = true;
		}
		return nlep_;
	}
	int &nosel()
	{
		if (not nosel_isLoaded) {
			if (nosel_branch != 0) {
				nosel_branch->GetEntry(index);
			} else { 
				printf("branch nosel_branch does not exist!\n");
				exit(1);
			}
			nosel_isLoaded = true;
		}
		return nosel_;
	}
	int &ngoodlep()
	{
		if (not ngoodlep_isLoaded) {
			if (ngoodlep_branch != 0) {
				ngoodlep_branch->GetEntry(index);
			} else { 
				printf("branch ngoodlep_branch does not exist!\n");
				exit(1);
			}
			ngoodlep_isLoaded = true;
		}
		return ngoodlep_;
	}
	int &ngoodel()
	{
		if (not ngoodel_isLoaded) {
			if (ngoodel_branch != 0) {
				ngoodel_branch->GetEntry(index);
			} else { 
				printf("branch ngoodel_branch does not exist!\n");
				exit(1);
			}
			ngoodel_isLoaded = true;
		}
		return ngoodel_;
	}
	int &ngoodmu()
	{
		if (not ngoodmu_isLoaded) {
			if (ngoodmu_branch != 0) {
				ngoodmu_branch->GetEntry(index);
			} else { 
				printf("branch ngoodmu_branch does not exist!\n");
				exit(1);
			}
			ngoodmu_isLoaded = true;
		}
		return ngoodmu_;
	}
	int &mull()
	{
		if (not mull_isLoaded) {
			if (mull_branch != 0) {
				mull_branch->GetEntry(index);
			} else { 
				printf("branch mull_branch does not exist!\n");
				exit(1);
			}
			mull_isLoaded = true;
		}
		return mull_;
	}
	int &mult()
	{
		if (not mult_isLoaded) {
			if (mult_branch != 0) {
				mult_branch->GetEntry(index);
			} else { 
				printf("branch mult_branch does not exist!\n");
				exit(1);
			}
			mult_isLoaded = true;
		}
		return mult_;
	}
	int &mullgen()
	{
		if (not mullgen_isLoaded) {
			if (mullgen_branch != 0) {
				mullgen_branch->GetEntry(index);
			} else { 
				printf("branch mullgen_branch does not exist!\n");
				exit(1);
			}
			mullgen_isLoaded = true;
		}
		return mullgen_;
	}
	int &multgen()
	{
		if (not multgen_isLoaded) {
			if (multgen_branch != 0) {
				multgen_branch->GetEntry(index);
			} else { 
				printf("branch multgen_branch does not exist!\n");
				exit(1);
			}
			multgen_isLoaded = true;
		}
		return multgen_;
	}
	int &proc()
	{
		if (not proc_isLoaded) {
			if (proc_branch != 0) {
				proc_branch->GetEntry(index);
			} else { 
				printf("branch proc_branch does not exist!\n");
				exit(1);
			}
			proc_isLoaded = true;
		}
		return proc_;
	}
	int &leptype()
	{
		if (not leptype_isLoaded) {
			if (leptype_branch != 0) {
				leptype_branch->GetEntry(index);
			} else { 
				printf("branch leptype_branch does not exist!\n");
				exit(1);
			}
			leptype_isLoaded = true;
		}
		return leptype_;
	}
	float &topmass()
	{
		if (not topmass_isLoaded) {
			if (topmass_branch != 0) {
				topmass_branch->GetEntry(index);
			} else { 
				printf("branch topmass_branch does not exist!\n");
				exit(1);
			}
			topmass_isLoaded = true;
		}
		return topmass_;
	}
	float &dilmass()
	{
		if (not dilmass_isLoaded) {
			if (dilmass_branch != 0) {
				dilmass_branch->GetEntry(index);
			} else { 
				printf("branch dilmass_branch does not exist!\n");
				exit(1);
			}
			dilmass_isLoaded = true;
		}
		return dilmass_;
	}
	float &dilrecoil()
	{
		if (not dilrecoil_isLoaded) {
			if (dilrecoil_branch != 0) {
				dilrecoil_branch->GetEntry(index);
			} else { 
				printf("branch dilrecoil_branch does not exist!\n");
				exit(1);
			}
			dilrecoil_isLoaded = true;
		}
		return dilrecoil_;
	}
	float &dilrecoilparl()
	{
		if (not dilrecoilparl_isLoaded) {
			if (dilrecoilparl_branch != 0) {
				dilrecoilparl_branch->GetEntry(index);
			} else { 
				printf("branch dilrecoilparl_branch does not exist!\n");
				exit(1);
			}
			dilrecoilparl_isLoaded = true;
		}
		return dilrecoilparl_;
	}
	float &dilrecoilperp()
	{
		if (not dilrecoilperp_isLoaded) {
			if (dilrecoilperp_branch != 0) {
				dilrecoilperp_branch->GetEntry(index);
			} else { 
				printf("branch dilrecoilperp_branch does not exist!\n");
				exit(1);
			}
			dilrecoilperp_isLoaded = true;
		}
		return dilrecoilperp_;
	}
	float &tcmet()
	{
		if (not tcmet_isLoaded) {
			if (tcmet_branch != 0) {
				tcmet_branch->GetEntry(index);
			} else { 
				printf("branch tcmet_branch does not exist!\n");
				exit(1);
			}
			tcmet_isLoaded = true;
		}
		return tcmet_;
	}
	float &genmet()
	{
		if (not genmet_isLoaded) {
			if (genmet_branch != 0) {
				genmet_branch->GetEntry(index);
			} else { 
				printf("branch genmet_branch does not exist!\n");
				exit(1);
			}
			genmet_isLoaded = true;
		}
		return genmet_;
	}
	float &gensumet()
	{
		if (not gensumet_isLoaded) {
			if (gensumet_branch != 0) {
				gensumet_branch->GetEntry(index);
			} else { 
				printf("branch gensumet_branch does not exist!\n");
				exit(1);
			}
			gensumet_isLoaded = true;
		}
		return gensumet_;
	}
	float &genmetphi()
	{
		if (not genmetphi_isLoaded) {
			if (genmetphi_branch != 0) {
				genmetphi_branch->GetEntry(index);
			} else { 
				printf("branch genmetphi_branch does not exist!\n");
				exit(1);
			}
			genmetphi_isLoaded = true;
		}
		return genmetphi_;
	}
	float &calomet()
	{
		if (not calomet_isLoaded) {
			if (calomet_branch != 0) {
				calomet_branch->GetEntry(index);
			} else { 
				printf("branch calomet_branch does not exist!\n");
				exit(1);
			}
			calomet_isLoaded = true;
		}
		return calomet_;
	}
	float &calometphi()
	{
		if (not calometphi_isLoaded) {
			if (calometphi_branch != 0) {
				calometphi_branch->GetEntry(index);
			} else { 
				printf("branch calometphi_branch does not exist!\n");
				exit(1);
			}
			calometphi_isLoaded = true;
		}
		return calometphi_;
	}
	float &trkmet()
	{
		if (not trkmet_isLoaded) {
			if (trkmet_branch != 0) {
				trkmet_branch->GetEntry(index);
			} else { 
				printf("branch trkmet_branch does not exist!\n");
				exit(1);
			}
			trkmet_isLoaded = true;
		}
		return trkmet_;
	}
	float &trkmetphi()
	{
		if (not trkmetphi_isLoaded) {
			if (trkmetphi_branch != 0) {
				trkmetphi_branch->GetEntry(index);
			} else { 
				printf("branch trkmetphi_branch does not exist!\n");
				exit(1);
			}
			trkmetphi_isLoaded = true;
		}
		return trkmetphi_;
	}
	float &pfmet()
	{
		if (not pfmet_isLoaded) {
			if (pfmet_branch != 0) {
				pfmet_branch->GetEntry(index);
			} else { 
				printf("branch pfmet_branch does not exist!\n");
				exit(1);
			}
			pfmet_isLoaded = true;
		}
		return pfmet_;
	}
	float &pfmetveto()
	{
		if (not pfmetveto_isLoaded) {
			if (pfmetveto_branch != 0) {
				pfmetveto_branch->GetEntry(index);
			} else { 
				printf("branch pfmetveto_branch does not exist!\n");
				exit(1);
			}
			pfmetveto_isLoaded = true;
		}
		return pfmetveto_;
	}
	float &pfmetsig()
	{
		if (not pfmetsig_isLoaded) {
			if (pfmetsig_branch != 0) {
				pfmetsig_branch->GetEntry(index);
			} else { 
				printf("branch pfmetsig_branch does not exist!\n");
				exit(1);
			}
			pfmetsig_isLoaded = true;
		}
		return pfmetsig_;
	}
	float &pfmetphi()
	{
		if (not pfmetphi_isLoaded) {
			if (pfmetphi_branch != 0) {
				pfmetphi_branch->GetEntry(index);
			} else { 
				printf("branch pfmetphi_branch does not exist!\n");
				exit(1);
			}
			pfmetphi_isLoaded = true;
		}
		return pfmetphi_;
	}
	float &pfsumet()
	{
		if (not pfsumet_isLoaded) {
			if (pfsumet_branch != 0) {
				pfsumet_branch->GetEntry(index);
			} else { 
				printf("branch pfsumet_branch does not exist!\n");
				exit(1);
			}
			pfsumet_isLoaded = true;
		}
		return pfsumet_;
	}
	float &mucormet()
	{
		if (not mucormet_isLoaded) {
			if (mucormet_branch != 0) {
				mucormet_branch->GetEntry(index);
			} else { 
				printf("branch mucormet_branch does not exist!\n");
				exit(1);
			}
			mucormet_isLoaded = true;
		}
		return mucormet_;
	}
	float &mucorjesmet()
	{
		if (not mucorjesmet_isLoaded) {
			if (mucorjesmet_branch != 0) {
				mucorjesmet_branch->GetEntry(index);
			} else { 
				printf("branch mucorjesmet_branch does not exist!\n");
				exit(1);
			}
			mucorjesmet_isLoaded = true;
		}
		return mucorjesmet_;
	}
	float &tcmet35X()
	{
		if (not tcmet35X_isLoaded) {
			if (tcmet35X_branch != 0) {
				tcmet35X_branch->GetEntry(index);
			} else { 
				printf("branch tcmet35X_branch does not exist!\n");
				exit(1);
			}
			tcmet35X_isLoaded = true;
		}
		return tcmet35X_;
	}
	float &tcmetevent()
	{
		if (not tcmetevent_isLoaded) {
			if (tcmetevent_branch != 0) {
				tcmetevent_branch->GetEntry(index);
			} else { 
				printf("branch tcmetevent_branch does not exist!\n");
				exit(1);
			}
			tcmetevent_isLoaded = true;
		}
		return tcmetevent_;
	}
	float &tcmetlooper()
	{
		if (not tcmetlooper_isLoaded) {
			if (tcmetlooper_branch != 0) {
				tcmetlooper_branch->GetEntry(index);
			} else { 
				printf("branch tcmetlooper_branch does not exist!\n");
				exit(1);
			}
			tcmetlooper_isLoaded = true;
		}
		return tcmetlooper_;
	}
	float &tcmetphi()
	{
		if (not tcmetphi_isLoaded) {
			if (tcmetphi_branch != 0) {
				tcmetphi_branch->GetEntry(index);
			} else { 
				printf("branch tcmetphi_branch does not exist!\n");
				exit(1);
			}
			tcmetphi_isLoaded = true;
		}
		return tcmetphi_;
	}
	float &tcsumet()
	{
		if (not tcsumet_isLoaded) {
			if (tcsumet_branch != 0) {
				tcsumet_branch->GetEntry(index);
			} else { 
				printf("branch tcsumet_branch does not exist!\n");
				exit(1);
			}
			tcsumet_isLoaded = true;
		}
		return tcsumet_;
	}
	float &tcmetUp()
	{
		if (not tcmetUp_isLoaded) {
			if (tcmetUp_branch != 0) {
				tcmetUp_branch->GetEntry(index);
			} else { 
				printf("branch tcmetUp_branch does not exist!\n");
				exit(1);
			}
			tcmetUp_isLoaded = true;
		}
		return tcmetUp_;
	}
	float &tcmetDown()
	{
		if (not tcmetDown_isLoaded) {
			if (tcmetDown_branch != 0) {
				tcmetDown_branch->GetEntry(index);
			} else { 
				printf("branch tcmetDown_branch does not exist!\n");
				exit(1);
			}
			tcmetDown_isLoaded = true;
		}
		return tcmetDown_;
	}
	float &tcmetTest()
	{
		if (not tcmetTest_isLoaded) {
			if (tcmetTest_branch != 0) {
				tcmetTest_branch->GetEntry(index);
			} else { 
				printf("branch tcmetTest_branch does not exist!\n");
				exit(1);
			}
			tcmetTest_isLoaded = true;
		}
		return tcmetTest_;
	}
	float &pfmetUp()
	{
		if (not pfmetUp_isLoaded) {
			if (pfmetUp_branch != 0) {
				pfmetUp_branch->GetEntry(index);
			} else { 
				printf("branch pfmetUp_branch does not exist!\n");
				exit(1);
			}
			pfmetUp_isLoaded = true;
		}
		return pfmetUp_;
	}
	float &pfmetDown()
	{
		if (not pfmetDown_isLoaded) {
			if (pfmetDown_branch != 0) {
				pfmetDown_branch->GetEntry(index);
			} else { 
				printf("branch pfmetDown_branch does not exist!\n");
				exit(1);
			}
			pfmetDown_isLoaded = true;
		}
		return pfmetDown_;
	}
	float &pfmetTest()
	{
		if (not pfmetTest_isLoaded) {
			if (pfmetTest_branch != 0) {
				pfmetTest_branch->GetEntry(index);
			} else { 
				printf("branch pfmetTest_branch does not exist!\n");
				exit(1);
			}
			pfmetTest_isLoaded = true;
		}
		return pfmetTest_;
	}
	float &sumjetpt()
	{
		if (not sumjetpt_isLoaded) {
			if (sumjetpt_branch != 0) {
				sumjetpt_branch->GetEntry(index);
			} else { 
				printf("branch sumjetpt_branch does not exist!\n");
				exit(1);
			}
			sumjetpt_isLoaded = true;
		}
		return sumjetpt_;
	}
	float &dileta()
	{
		if (not dileta_isLoaded) {
			if (dileta_branch != 0) {
				dileta_branch->GetEntry(index);
			} else { 
				printf("branch dileta_branch does not exist!\n");
				exit(1);
			}
			dileta_isLoaded = true;
		}
		return dileta_;
	}
	float &dilpt()
	{
		if (not dilpt_isLoaded) {
			if (dilpt_branch != 0) {
				dilpt_branch->GetEntry(index);
			} else { 
				printf("branch dilpt_branch does not exist!\n");
				exit(1);
			}
			dilpt_isLoaded = true;
		}
		return dilpt_;
	}
	float &dildphi()
	{
		if (not dildphi_isLoaded) {
			if (dildphi_branch != 0) {
				dildphi_branch->GetEntry(index);
			} else { 
				printf("branch dildphi_branch does not exist!\n");
				exit(1);
			}
			dildphi_isLoaded = true;
		}
		return dildphi_;
	}
	int &ngenjets()
	{
		if (not ngenjets_isLoaded) {
			if (ngenjets_branch != 0) {
				ngenjets_branch->GetEntry(index);
			} else { 
				printf("branch ngenjets_branch does not exist!\n");
				exit(1);
			}
			ngenjets_isLoaded = true;
		}
		return ngenjets_;
	}
	int &njpt()
	{
		if (not njpt_isLoaded) {
			if (njpt_branch != 0) {
				njpt_branch->GetEntry(index);
			} else { 
				printf("branch njpt_branch does not exist!\n");
				exit(1);
			}
			njpt_isLoaded = true;
		}
		return njpt_;
	}
	int &trgmu1()
	{
		if (not trgmu1_isLoaded) {
			if (trgmu1_branch != 0) {
				trgmu1_branch->GetEntry(index);
			} else { 
				printf("branch trgmu1_branch does not exist!\n");
				exit(1);
			}
			trgmu1_isLoaded = true;
		}
		return trgmu1_;
	}
	int &trgmu2()
	{
		if (not trgmu2_isLoaded) {
			if (trgmu2_branch != 0) {
				trgmu2_branch->GetEntry(index);
			} else { 
				printf("branch trgmu2_branch does not exist!\n");
				exit(1);
			}
			trgmu2_isLoaded = true;
		}
		return trgmu2_;
	}
	int &trgel1()
	{
		if (not trgel1_isLoaded) {
			if (trgel1_branch != 0) {
				trgel1_branch->GetEntry(index);
			} else { 
				printf("branch trgel1_branch does not exist!\n");
				exit(1);
			}
			trgel1_isLoaded = true;
		}
		return trgel1_;
	}
	int &trgel2()
	{
		if (not trgel2_isLoaded) {
			if (trgel2_branch != 0) {
				trgel2_branch->GetEntry(index);
			} else { 
				printf("branch trgel2_branch does not exist!\n");
				exit(1);
			}
			trgel2_isLoaded = true;
		}
		return trgel2_;
	}
	int &isomu24()
	{
		if (not isomu24_isLoaded) {
			if (isomu24_branch != 0) {
				isomu24_branch->GetEntry(index);
			} else { 
				printf("branch isomu24_branch does not exist!\n");
				exit(1);
			}
			isomu24_isLoaded = true;
		}
		return isomu24_;
	}
	int &ele27wp80()
	{
		if (not ele27wp80_isLoaded) {
			if (ele27wp80_branch != 0) {
				ele27wp80_branch->GetEntry(index);
			} else { 
				printf("branch ele27wp80_branch does not exist!\n");
				exit(1);
			}
			ele27wp80_isLoaded = true;
		}
		return ele27wp80_;
	}
	int &mm()
	{
		if (not mm_isLoaded) {
			if (mm_branch != 0) {
				mm_branch->GetEntry(index);
			} else { 
				printf("branch mm_branch does not exist!\n");
				exit(1);
			}
			mm_isLoaded = true;
		}
		return mm_;
	}
	int &mmtk()
	{
		if (not mmtk_isLoaded) {
			if (mmtk_branch != 0) {
				mmtk_branch->GetEntry(index);
			} else { 
				printf("branch mmtk_branch does not exist!\n");
				exit(1);
			}
			mmtk_isLoaded = true;
		}
		return mmtk_;
	}
	int &me()
	{
		if (not me_isLoaded) {
			if (me_branch != 0) {
				me_branch->GetEntry(index);
			} else { 
				printf("branch me_branch does not exist!\n");
				exit(1);
			}
			me_isLoaded = true;
		}
		return me_;
	}
	int &em()
	{
		if (not em_isLoaded) {
			if (em_branch != 0) {
				em_branch->GetEntry(index);
			} else { 
				printf("branch em_branch does not exist!\n");
				exit(1);
			}
			em_isLoaded = true;
		}
		return em_;
	}
	int &mu()
	{
		if (not mu_isLoaded) {
			if (mu_branch != 0) {
				mu_branch->GetEntry(index);
			} else { 
				printf("branch mu_branch does not exist!\n");
				exit(1);
			}
			mu_isLoaded = true;
		}
		return mu_;
	}
	int &ee()
	{
		if (not ee_isLoaded) {
			if (ee_branch != 0) {
				ee_branch->GetEntry(index);
			} else { 
				printf("branch ee_branch does not exist!\n");
				exit(1);
			}
			ee_isLoaded = true;
		}
		return ee_;
	}
	int &npfjets30()
	{
		if (not npfjets30_isLoaded) {
			if (npfjets30_branch != 0) {
				npfjets30_branch->GetEntry(index);
			} else { 
				printf("branch npfjets30_branch does not exist!\n");
				exit(1);
			}
			npfjets30_isLoaded = true;
		}
		return npfjets30_;
	}
	int &npfjets30lepcorr()
	{
		if (not npfjets30lepcorr_isLoaded) {
			if (npfjets30lepcorr_branch != 0) {
				npfjets30lepcorr_branch->GetEntry(index);
			} else { 
				printf("branch npfjets30lepcorr_branch does not exist!\n");
				exit(1);
			}
			npfjets30lepcorr_isLoaded = true;
		}
		return npfjets30lepcorr_;
	}
	float &knjets()
	{
		if (not knjets_isLoaded) {
			if (knjets_branch != 0) {
				knjets_branch->GetEntry(index);
			} else { 
				printf("branch knjets_branch does not exist!\n");
				exit(1);
			}
			knjets_isLoaded = true;
		}
		return knjets_;
	}
	float &rhovor()
	{
		if (not rhovor_isLoaded) {
			if (rhovor_branch != 0) {
				rhovor_branch->GetEntry(index);
			} else { 
				printf("branch rhovor_branch does not exist!\n");
				exit(1);
			}
			rhovor_isLoaded = true;
		}
		return rhovor_;
	}
	float &htpf30()
	{
		if (not htpf30_isLoaded) {
			if (htpf30_branch != 0) {
				htpf30_branch->GetEntry(index);
			} else { 
				printf("branch htpf30_branch does not exist!\n");
				exit(1);
			}
			htpf30_isLoaded = true;
		}
		return htpf30_;
	}
	float &t1met10()
	{
		if (not t1met10_isLoaded) {
			if (t1met10_branch != 0) {
				t1met10_branch->GetEntry(index);
			} else { 
				printf("branch t1met10_branch does not exist!\n");
				exit(1);
			}
			t1met10_isLoaded = true;
		}
		return t1met10_;
	}
	float &t1met20()
	{
		if (not t1met20_isLoaded) {
			if (t1met20_branch != 0) {
				t1met20_branch->GetEntry(index);
			} else { 
				printf("branch t1met20_branch does not exist!\n");
				exit(1);
			}
			t1met20_isLoaded = true;
		}
		return t1met20_;
	}
	float &t1met30()
	{
		if (not t1met30_isLoaded) {
			if (t1met30_branch != 0) {
				t1met30_branch->GetEntry(index);
			} else { 
				printf("branch t1met30_branch does not exist!\n");
				exit(1);
			}
			t1met30_isLoaded = true;
		}
		return t1met30_;
	}
	float &t1met10phi()
	{
		if (not t1met10phi_isLoaded) {
			if (t1met10phi_branch != 0) {
				t1met10phi_branch->GetEntry(index);
			} else { 
				printf("branch t1met10phi_branch does not exist!\n");
				exit(1);
			}
			t1met10phi_isLoaded = true;
		}
		return t1met10phi_;
	}
	float &t1met20phi()
	{
		if (not t1met20phi_isLoaded) {
			if (t1met20phi_branch != 0) {
				t1met20phi_branch->GetEntry(index);
			} else { 
				printf("branch t1met20phi_branch does not exist!\n");
				exit(1);
			}
			t1met20phi_isLoaded = true;
		}
		return t1met20phi_;
	}
	float &t1met30phi()
	{
		if (not t1met30phi_isLoaded) {
			if (t1met30phi_branch != 0) {
				t1met30phi_branch->GetEntry(index);
			} else { 
				printf("branch t1met30phi_branch does not exist!\n");
				exit(1);
			}
			t1met30phi_isLoaded = true;
		}
		return t1met30phi_;
	}
	float &t1met10mt()
	{
		if (not t1met10mt_isLoaded) {
			if (t1met10mt_branch != 0) {
				t1met10mt_branch->GetEntry(index);
			} else { 
				printf("branch t1met10mt_branch does not exist!\n");
				exit(1);
			}
			t1met10mt_isLoaded = true;
		}
		return t1met10mt_;
	}
	float &t1met20mt()
	{
		if (not t1met20mt_isLoaded) {
			if (t1met20mt_branch != 0) {
				t1met20mt_branch->GetEntry(index);
			} else { 
				printf("branch t1met20mt_branch does not exist!\n");
				exit(1);
			}
			t1met20mt_isLoaded = true;
		}
		return t1met20mt_;
	}
	float &t1met30mt()
	{
		if (not t1met30mt_isLoaded) {
			if (t1met30mt_branch != 0) {
				t1met30mt_branch->GetEntry(index);
			} else { 
				printf("branch t1met30mt_branch does not exist!\n");
				exit(1);
			}
			t1met30mt_isLoaded = true;
		}
		return t1met30mt_;
	}
	float &lepmetpt()
	{
		if (not lepmetpt_isLoaded) {
			if (lepmetpt_branch != 0) {
				lepmetpt_branch->GetEntry(index);
			} else { 
				printf("branch lepmetpt_branch does not exist!\n");
				exit(1);
			}
			lepmetpt_isLoaded = true;
		}
		return lepmetpt_;
	}
	float &lept1met10pt()
	{
		if (not lept1met10pt_isLoaded) {
			if (lept1met10pt_branch != 0) {
				lept1met10pt_branch->GetEntry(index);
			} else { 
				printf("branch lept1met10pt_branch does not exist!\n");
				exit(1);
			}
			lept1met10pt_isLoaded = true;
		}
		return lept1met10pt_;
	}
	float &t1met10s()
	{
		if (not t1met10s_isLoaded) {
			if (t1met10s_branch != 0) {
				t1met10s_branch->GetEntry(index);
			} else { 
				printf("branch t1met10s_branch does not exist!\n");
				exit(1);
			}
			t1met10s_isLoaded = true;
		}
		return t1met10s_;
	}
	float &t1met10sphi()
	{
		if (not t1met10sphi_isLoaded) {
			if (t1met10sphi_branch != 0) {
				t1met10sphi_branch->GetEntry(index);
			} else { 
				printf("branch t1met10sphi_branch does not exist!\n");
				exit(1);
			}
			t1met10sphi_isLoaded = true;
		}
		return t1met10sphi_;
	}
	float &t1met10smt()
	{
		if (not t1met10smt_isLoaded) {
			if (t1met10smt_branch != 0) {
				t1met10smt_branch->GetEntry(index);
			} else { 
				printf("branch t1met10smt_branch does not exist!\n");
				exit(1);
			}
			t1met10smt_isLoaded = true;
		}
		return t1met10smt_;
	}
	float &t1metphicorr()
	{
		if (not t1metphicorr_isLoaded) {
			if (t1metphicorr_branch != 0) {
				t1metphicorr_branch->GetEntry(index);
			} else { 
				printf("branch t1metphicorr_branch does not exist!\n");
				exit(1);
			}
			t1metphicorr_isLoaded = true;
		}
		return t1metphicorr_;
	}
	float &t1metphicorrup()
	{
		if (not t1metphicorrup_isLoaded) {
			if (t1metphicorrup_branch != 0) {
				t1metphicorrup_branch->GetEntry(index);
			} else { 
				printf("branch t1metphicorrup_branch does not exist!\n");
				exit(1);
			}
			t1metphicorrup_isLoaded = true;
		}
		return t1metphicorrup_;
	}
	float &t1metphicorrdn()
	{
		if (not t1metphicorrdn_isLoaded) {
			if (t1metphicorrdn_branch != 0) {
				t1metphicorrdn_branch->GetEntry(index);
			} else { 
				printf("branch t1metphicorrdn_branch does not exist!\n");
				exit(1);
			}
			t1metphicorrdn_isLoaded = true;
		}
		return t1metphicorrdn_;
	}
	float &t1metphicorrphi()
	{
		if (not t1metphicorrphi_isLoaded) {
			if (t1metphicorrphi_branch != 0) {
				t1metphicorrphi_branch->GetEntry(index);
			} else { 
				printf("branch t1metphicorrphi_branch does not exist!\n");
				exit(1);
			}
			t1metphicorrphi_isLoaded = true;
		}
		return t1metphicorrphi_;
	}
	float &t1metphicorrphiup()
	{
		if (not t1metphicorrphiup_isLoaded) {
			if (t1metphicorrphiup_branch != 0) {
				t1metphicorrphiup_branch->GetEntry(index);
			} else { 
				printf("branch t1metphicorrphiup_branch does not exist!\n");
				exit(1);
			}
			t1metphicorrphiup_isLoaded = true;
		}
		return t1metphicorrphiup_;
	}
	float &t1metphicorrphidn()
	{
		if (not t1metphicorrphidn_isLoaded) {
			if (t1metphicorrphidn_branch != 0) {
				t1metphicorrphidn_branch->GetEntry(index);
			} else { 
				printf("branch t1metphicorrphidn_branch does not exist!\n");
				exit(1);
			}
			t1metphicorrphidn_isLoaded = true;
		}
		return t1metphicorrphidn_;
	}
	float &t1metphicorrlep()
	{
		if (not t1metphicorrlep_isLoaded) {
			if (t1metphicorrlep_branch != 0) {
				t1metphicorrlep_branch->GetEntry(index);
			} else { 
				printf("branch t1metphicorrlep_branch does not exist!\n");
				exit(1);
			}
			t1metphicorrlep_isLoaded = true;
		}
		return t1metphicorrlep_;
	}
	float &t1metphicorrlepphi()
	{
		if (not t1metphicorrlepphi_isLoaded) {
			if (t1metphicorrlepphi_branch != 0) {
				t1metphicorrlepphi_branch->GetEntry(index);
			} else { 
				printf("branch t1metphicorrlepphi_branch does not exist!\n");
				exit(1);
			}
			t1metphicorrlepphi_isLoaded = true;
		}
		return t1metphicorrlepphi_;
	}
	float &t1metphicorrmt()
	{
		if (not t1metphicorrmt_isLoaded) {
			if (t1metphicorrmt_branch != 0) {
				t1metphicorrmt_branch->GetEntry(index);
			} else { 
				printf("branch t1metphicorrmt_branch does not exist!\n");
				exit(1);
			}
			t1metphicorrmt_isLoaded = true;
		}
		return t1metphicorrmt_;
	}
	float &t1metphicorrmtup()
	{
		if (not t1metphicorrmtup_isLoaded) {
			if (t1metphicorrmtup_branch != 0) {
				t1metphicorrmtup_branch->GetEntry(index);
			} else { 
				printf("branch t1metphicorrmtup_branch does not exist!\n");
				exit(1);
			}
			t1metphicorrmtup_isLoaded = true;
		}
		return t1metphicorrmtup_;
	}
	float &t1metphicorrmtdn()
	{
		if (not t1metphicorrmtdn_isLoaded) {
			if (t1metphicorrmtdn_branch != 0) {
				t1metphicorrmtdn_branch->GetEntry(index);
			} else { 
				printf("branch t1metphicorrmtdn_branch does not exist!\n");
				exit(1);
			}
			t1metphicorrmtdn_isLoaded = true;
		}
		return t1metphicorrmtdn_;
	}
	float &t1metphicorrlepmt()
	{
		if (not t1metphicorrlepmt_isLoaded) {
			if (t1metphicorrlepmt_branch != 0) {
				t1metphicorrlepmt_branch->GetEntry(index);
			} else { 
				printf("branch t1metphicorrlepmt_branch does not exist!\n");
				exit(1);
			}
			t1metphicorrlepmt_isLoaded = true;
		}
		return t1metphicorrlepmt_;
	}
	float &t1met_off()
	{
		if (not t1met_off_isLoaded) {
			if (t1met_off_branch != 0) {
				t1met_off_branch->GetEntry(index);
			} else { 
				printf("branch t1met_off_branch does not exist!\n");
				exit(1);
			}
			t1met_off_isLoaded = true;
		}
		return t1met_off_;
	}
	float &t1metphi_off()
	{
		if (not t1metphi_off_isLoaded) {
			if (t1metphi_off_branch != 0) {
				t1metphi_off_branch->GetEntry(index);
			} else { 
				printf("branch t1metphi_off_branch does not exist!\n");
				exit(1);
			}
			t1metphi_off_isLoaded = true;
		}
		return t1metphi_off_;
	}
	float &t1metmt_off()
	{
		if (not t1metmt_off_isLoaded) {
			if (t1metmt_off_branch != 0) {
				t1metmt_off_branch->GetEntry(index);
			} else { 
				printf("branch t1metmt_off_branch does not exist!\n");
				exit(1);
			}
			t1metmt_off_isLoaded = true;
		}
		return t1metmt_off_;
	}
	float &t1metphicorr_off()
	{
		if (not t1metphicorr_off_isLoaded) {
			if (t1metphicorr_off_branch != 0) {
				t1metphicorr_off_branch->GetEntry(index);
			} else { 
				printf("branch t1metphicorr_off_branch does not exist!\n");
				exit(1);
			}
			t1metphicorr_off_isLoaded = true;
		}
		return t1metphicorr_off_;
	}
	float &t1metphicorrphi_off()
	{
		if (not t1metphicorrphi_off_isLoaded) {
			if (t1metphicorrphi_off_branch != 0) {
				t1metphicorrphi_off_branch->GetEntry(index);
			} else { 
				printf("branch t1metphicorrphi_off_branch does not exist!\n");
				exit(1);
			}
			t1metphicorrphi_off_isLoaded = true;
		}
		return t1metphicorrphi_off_;
	}
	float &t1metphicorrmt_off()
	{
		if (not t1metphicorrmt_off_isLoaded) {
			if (t1metphicorrmt_off_branch != 0) {
				t1metphicorrmt_off_branch->GetEntry(index);
			} else { 
				printf("branch t1metphicorrmt_off_branch does not exist!\n");
				exit(1);
			}
			t1metphicorrmt_off_isLoaded = true;
		}
		return t1metphicorrmt_off_;
	}
	float &mht15()
	{
		if (not mht15_isLoaded) {
			if (mht15_branch != 0) {
				mht15_branch->GetEntry(index);
			} else { 
				printf("branch mht15_branch does not exist!\n");
				exit(1);
			}
			mht15_isLoaded = true;
		}
		return mht15_;
	}
	float &mht15phi()
	{
		if (not mht15phi_isLoaded) {
			if (mht15phi_branch != 0) {
				mht15phi_branch->GetEntry(index);
			} else { 
				printf("branch mht15phi_branch does not exist!\n");
				exit(1);
			}
			mht15phi_isLoaded = true;
		}
		return mht15phi_;
	}
	float &trkmet_mht15()
	{
		if (not trkmet_mht15_isLoaded) {
			if (trkmet_mht15_branch != 0) {
				trkmet_mht15_branch->GetEntry(index);
			} else { 
				printf("branch trkmet_mht15_branch does not exist!\n");
				exit(1);
			}
			trkmet_mht15_isLoaded = true;
		}
		return trkmet_mht15_;
	}
	float &trkmetphi_mht15()
	{
		if (not trkmetphi_mht15_isLoaded) {
			if (trkmetphi_mht15_branch != 0) {
				trkmetphi_mht15_branch->GetEntry(index);
			} else { 
				printf("branch trkmetphi_mht15_branch does not exist!\n");
				exit(1);
			}
			trkmetphi_mht15_isLoaded = true;
		}
		return trkmetphi_mht15_;
	}
	float &mettlj15()
	{
		if (not mettlj15_isLoaded) {
			if (mettlj15_branch != 0) {
				mettlj15_branch->GetEntry(index);
			} else { 
				printf("branch mettlj15_branch does not exist!\n");
				exit(1);
			}
			mettlj15_isLoaded = true;
		}
		return mettlj15_;
	}
	float &mettlj15phi()
	{
		if (not mettlj15phi_isLoaded) {
			if (mettlj15phi_branch != 0) {
				mettlj15phi_branch->GetEntry(index);
			} else { 
				printf("branch mettlj15phi_branch does not exist!\n");
				exit(1);
			}
			mettlj15phi_isLoaded = true;
		}
		return mettlj15phi_;
	}
	float &mt2bmin()
	{
		if (not mt2bmin_isLoaded) {
			if (mt2bmin_branch != 0) {
				mt2bmin_branch->GetEntry(index);
			} else { 
				printf("branch mt2bmin_branch does not exist!\n");
				exit(1);
			}
			mt2bmin_isLoaded = true;
		}
		return mt2bmin_;
	}
	float &mt2blmin()
	{
		if (not mt2blmin_isLoaded) {
			if (mt2blmin_branch != 0) {
				mt2blmin_branch->GetEntry(index);
			} else { 
				printf("branch mt2blmin_branch does not exist!\n");
				exit(1);
			}
			mt2blmin_isLoaded = true;
		}
		return mt2blmin_;
	}
	float &mt2wmin()
	{
		if (not mt2wmin_isLoaded) {
			if (mt2wmin_branch != 0) {
				mt2wmin_branch->GetEntry(index);
			} else { 
				printf("branch mt2wmin_branch does not exist!\n");
				exit(1);
			}
			mt2wmin_isLoaded = true;
		}
		return mt2wmin_;
	}
	float &chi2min()
	{
		if (not chi2min_isLoaded) {
			if (chi2min_branch != 0) {
				chi2min_branch->GetEntry(index);
			} else { 
				printf("branch chi2min_branch does not exist!\n");
				exit(1);
			}
			chi2min_isLoaded = true;
		}
		return chi2min_;
	}
	float &chi2minprob()
	{
		if (not chi2minprob_isLoaded) {
			if (chi2minprob_branch != 0) {
				chi2minprob_branch->GetEntry(index);
			} else { 
				printf("branch chi2minprob_branch does not exist!\n");
				exit(1);
			}
			chi2minprob_isLoaded = true;
		}
		return chi2minprob_;
	}
	int &nbtagsssv()
	{
		if (not nbtagsssv_isLoaded) {
			if (nbtagsssv_branch != 0) {
				nbtagsssv_branch->GetEntry(index);
			} else { 
				printf("branch nbtagsssv_branch does not exist!\n");
				exit(1);
			}
			nbtagsssv_isLoaded = true;
		}
		return nbtagsssv_;
	}
	int &nbtagstcl()
	{
		if (not nbtagstcl_isLoaded) {
			if (nbtagstcl_branch != 0) {
				nbtagstcl_branch->GetEntry(index);
			} else { 
				printf("branch nbtagstcl_branch does not exist!\n");
				exit(1);
			}
			nbtagstcl_isLoaded = true;
		}
		return nbtagstcl_;
	}
	int &nbtagstcm()
	{
		if (not nbtagstcm_isLoaded) {
			if (nbtagstcm_branch != 0) {
				nbtagstcm_branch->GetEntry(index);
			} else { 
				printf("branch nbtagstcm_branch does not exist!\n");
				exit(1);
			}
			nbtagstcm_isLoaded = true;
		}
		return nbtagstcm_;
	}
	int &nbtagscsvl()
	{
		if (not nbtagscsvl_isLoaded) {
			if (nbtagscsvl_branch != 0) {
				nbtagscsvl_branch->GetEntry(index);
			} else { 
				printf("branch nbtagscsvl_branch does not exist!\n");
				exit(1);
			}
			nbtagscsvl_isLoaded = true;
		}
		return nbtagscsvl_;
	}
	int &nbtagscsvm()
	{
		if (not nbtagscsvm_isLoaded) {
			if (nbtagscsvm_branch != 0) {
				nbtagscsvm_branch->GetEntry(index);
			} else { 
				printf("branch nbtagscsvm_branch does not exist!\n");
				exit(1);
			}
			nbtagscsvm_isLoaded = true;
		}
		return nbtagscsvm_;
	}
	int &nbtagscsvt()
	{
		if (not nbtagscsvt_isLoaded) {
			if (nbtagscsvt_branch != 0) {
				nbtagscsvt_branch->GetEntry(index);
			} else { 
				printf("branch nbtagscsvt_branch does not exist!\n");
				exit(1);
			}
			nbtagscsvt_isLoaded = true;
		}
		return nbtagscsvt_;
	}
	int &nbtagsssvcorr()
	{
		if (not nbtagsssvcorr_isLoaded) {
			if (nbtagsssvcorr_branch != 0) {
				nbtagsssvcorr_branch->GetEntry(index);
			} else { 
				printf("branch nbtagsssvcorr_branch does not exist!\n");
				exit(1);
			}
			nbtagsssvcorr_isLoaded = true;
		}
		return nbtagsssvcorr_;
	}
	int &nbtagstclcorr()
	{
		if (not nbtagstclcorr_isLoaded) {
			if (nbtagstclcorr_branch != 0) {
				nbtagstclcorr_branch->GetEntry(index);
			} else { 
				printf("branch nbtagstclcorr_branch does not exist!\n");
				exit(1);
			}
			nbtagstclcorr_isLoaded = true;
		}
		return nbtagstclcorr_;
	}
	int &nbtagstcmcorr()
	{
		if (not nbtagstcmcorr_isLoaded) {
			if (nbtagstcmcorr_branch != 0) {
				nbtagstcmcorr_branch->GetEntry(index);
			} else { 
				printf("branch nbtagstcmcorr_branch does not exist!\n");
				exit(1);
			}
			nbtagstcmcorr_isLoaded = true;
		}
		return nbtagstcmcorr_;
	}
	int &nbtagscsvlcorr()
	{
		if (not nbtagscsvlcorr_isLoaded) {
			if (nbtagscsvlcorr_branch != 0) {
				nbtagscsvlcorr_branch->GetEntry(index);
			} else { 
				printf("branch nbtagscsvlcorr_branch does not exist!\n");
				exit(1);
			}
			nbtagscsvlcorr_isLoaded = true;
		}
		return nbtagscsvlcorr_;
	}
	int &nbtagscsvmcorr()
	{
		if (not nbtagscsvmcorr_isLoaded) {
			if (nbtagscsvmcorr_branch != 0) {
				nbtagscsvmcorr_branch->GetEntry(index);
			} else { 
				printf("branch nbtagscsvmcorr_branch does not exist!\n");
				exit(1);
			}
			nbtagscsvmcorr_isLoaded = true;
		}
		return nbtagscsvmcorr_;
	}
	int &nbtagscsvtcott()
	{
		if (not nbtagscsvtcott_isLoaded) {
			if (nbtagscsvtcott_branch != 0) {
				nbtagscsvtcott_branch->GetEntry(index);
			} else { 
				printf("branch nbtagscsvtcott_branch does not exist!\n");
				exit(1);
			}
			nbtagscsvtcott_isLoaded = true;
		}
		return nbtagscsvtcott_;
	}
	int &njetsUp()
	{
		if (not njetsUp_isLoaded) {
			if (njetsUp_branch != 0) {
				njetsUp_branch->GetEntry(index);
			} else { 
				printf("branch njetsUp_branch does not exist!\n");
				exit(1);
			}
			njetsUp_isLoaded = true;
		}
		return njetsUp_;
	}
	int &njetsDown()
	{
		if (not njetsDown_isLoaded) {
			if (njetsDown_branch != 0) {
				njetsDown_branch->GetEntry(index);
			} else { 
				printf("branch njetsDown_branch does not exist!\n");
				exit(1);
			}
			njetsDown_isLoaded = true;
		}
		return njetsDown_;
	}
	float &htUp()
	{
		if (not htUp_isLoaded) {
			if (htUp_branch != 0) {
				htUp_branch->GetEntry(index);
			} else { 
				printf("branch htUp_branch does not exist!\n");
				exit(1);
			}
			htUp_isLoaded = true;
		}
		return htUp_;
	}
	float &htDown()
	{
		if (not htDown_isLoaded) {
			if (htDown_branch != 0) {
				htDown_branch->GetEntry(index);
			} else { 
				printf("branch htDown_branch does not exist!\n");
				exit(1);
			}
			htDown_isLoaded = true;
		}
		return htDown_;
	}
	int &ntruepu()
	{
		if (not ntruepu_isLoaded) {
			if (ntruepu_branch != 0) {
				ntruepu_branch->GetEntry(index);
			} else { 
				printf("branch ntruepu_branch does not exist!\n");
				exit(1);
			}
			ntruepu_isLoaded = true;
		}
		return ntruepu_;
	}
	int &npu()
	{
		if (not npu_isLoaded) {
			if (npu_branch != 0) {
				npu_branch->GetEntry(index);
			} else { 
				printf("branch npu_branch does not exist!\n");
				exit(1);
			}
			npu_isLoaded = true;
		}
		return npu_;
	}
	int &npuMinusOne()
	{
		if (not npuMinusOne_isLoaded) {
			if (npuMinusOne_branch != 0) {
				npuMinusOne_branch->GetEntry(index);
			} else { 
				printf("branch npuMinusOne_branch does not exist!\n");
				exit(1);
			}
			npuMinusOne_isLoaded = true;
		}
		return npuMinusOne_;
	}
	int &npuPlusOne()
	{
		if (not npuPlusOne_isLoaded) {
			if (npuPlusOne_branch != 0) {
				npuPlusOne_branch->GetEntry(index);
			} else { 
				printf("branch npuPlusOne_branch does not exist!\n");
				exit(1);
			}
			npuPlusOne_isLoaded = true;
		}
		return npuPlusOne_;
	}
	int &nvtx()
	{
		if (not nvtx_isLoaded) {
			if (nvtx_branch != 0) {
				nvtx_branch->GetEntry(index);
			} else { 
				printf("branch nvtx_branch does not exist!\n");
				exit(1);
			}
			nvtx_isLoaded = true;
		}
		return nvtx_;
	}
	int &indexfirstGoodVertex_()
	{
		if (not indexfirstGoodVertex__isLoaded) {
			if (indexfirstGoodVertex__branch != 0) {
				indexfirstGoodVertex__branch->GetEntry(index);
			} else { 
				printf("branch indexfirstGoodVertex__branch does not exist!\n");
				exit(1);
			}
			indexfirstGoodVertex__isLoaded = true;
		}
		return indexfirstGoodVertex__;
	}
	float &nvtxweight()
	{
		if (not nvtxweight_isLoaded) {
			if (nvtxweight_branch != 0) {
				nvtxweight_branch->GetEntry(index);
			} else { 
				printf("branch nvtxweight_branch does not exist!\n");
				exit(1);
			}
			nvtxweight_isLoaded = true;
		}
		return nvtxweight_;
	}
	float &n3dvtxweight()
	{
		if (not n3dvtxweight_isLoaded) {
			if (n3dvtxweight_branch != 0) {
				n3dvtxweight_branch->GetEntry(index);
			} else { 
				printf("branch n3dvtxweight_branch does not exist!\n");
				exit(1);
			}
			n3dvtxweight_isLoaded = true;
		}
		return n3dvtxweight_;
	}
	int &pdfid1()
	{
		if (not pdfid1_isLoaded) {
			if (pdfid1_branch != 0) {
				pdfid1_branch->GetEntry(index);
			} else { 
				printf("branch pdfid1_branch does not exist!\n");
				exit(1);
			}
			pdfid1_isLoaded = true;
		}
		return pdfid1_;
	}
	int &pdfid2()
	{
		if (not pdfid2_isLoaded) {
			if (pdfid2_branch != 0) {
				pdfid2_branch->GetEntry(index);
			} else { 
				printf("branch pdfid2_branch does not exist!\n");
				exit(1);
			}
			pdfid2_isLoaded = true;
		}
		return pdfid2_;
	}
	float &pdfx1()
	{
		if (not pdfx1_isLoaded) {
			if (pdfx1_branch != 0) {
				pdfx1_branch->GetEntry(index);
			} else { 
				printf("branch pdfx1_branch does not exist!\n");
				exit(1);
			}
			pdfx1_isLoaded = true;
		}
		return pdfx1_;
	}
	float &pdfx2()
	{
		if (not pdfx2_isLoaded) {
			if (pdfx2_branch != 0) {
				pdfx2_branch->GetEntry(index);
			} else { 
				printf("branch pdfx2_branch does not exist!\n");
				exit(1);
			}
			pdfx2_isLoaded = true;
		}
		return pdfx2_;
	}
	float &pdfQ()
	{
		if (not pdfQ_isLoaded) {
			if (pdfQ_branch != 0) {
				pdfQ_branch->GetEntry(index);
			} else { 
				printf("branch pdfQ_branch does not exist!\n");
				exit(1);
			}
			pdfQ_isLoaded = true;
		}
		return pdfQ_;
	}
	float &vecjetpt()
	{
		if (not vecjetpt_isLoaded) {
			if (vecjetpt_branch != 0) {
				vecjetpt_branch->GetEntry(index);
			} else { 
				printf("branch vecjetpt_branch does not exist!\n");
				exit(1);
			}
			vecjetpt_isLoaded = true;
		}
		return vecjetpt_;
	}
	int &pass()
	{
		if (not pass_isLoaded) {
			if (pass_branch != 0) {
				pass_branch->GetEntry(index);
			} else { 
				printf("branch pass_branch does not exist!\n");
				exit(1);
			}
			pass_isLoaded = true;
		}
		return pass_;
	}
	int &passz()
	{
		if (not passz_isLoaded) {
			if (passz_branch != 0) {
				passz_branch->GetEntry(index);
			} else { 
				printf("branch passz_branch does not exist!\n");
				exit(1);
			}
			passz_isLoaded = true;
		}
		return passz_;
	}
	float &m0()
	{
		if (not m0_isLoaded) {
			if (m0_branch != 0) {
				m0_branch->GetEntry(index);
			} else { 
				printf("branch m0_branch does not exist!\n");
				exit(1);
			}
			m0_isLoaded = true;
		}
		return m0_;
	}
	float &mg()
	{
		if (not mg_isLoaded) {
			if (mg_branch != 0) {
				mg_branch->GetEntry(index);
			} else { 
				printf("branch mg_branch does not exist!\n");
				exit(1);
			}
			mg_isLoaded = true;
		}
		return mg_;
	}
	float &ml()
	{
		if (not ml_isLoaded) {
			if (ml_branch != 0) {
				ml_branch->GetEntry(index);
			} else { 
				printf("branch ml_branch does not exist!\n");
				exit(1);
			}
			ml_isLoaded = true;
		}
		return ml_;
	}
	float &x()
	{
		if (not x_isLoaded) {
			if (x_branch != 0) {
				x_branch->GetEntry(index);
			} else { 
				printf("branch x_branch does not exist!\n");
				exit(1);
			}
			x_isLoaded = true;
		}
		return x_;
	}
	float &m12()
	{
		if (not m12_isLoaded) {
			if (m12_branch != 0) {
				m12_branch->GetEntry(index);
			} else { 
				printf("branch m12_branch does not exist!\n");
				exit(1);
			}
			m12_isLoaded = true;
		}
		return m12_;
	}
	float &lep1chi2ndf()
	{
		if (not lep1chi2ndf_isLoaded) {
			if (lep1chi2ndf_branch != 0) {
				lep1chi2ndf_branch->GetEntry(index);
			} else { 
				printf("branch lep1chi2ndf_branch does not exist!\n");
				exit(1);
			}
			lep1chi2ndf_isLoaded = true;
		}
		return lep1chi2ndf_;
	}
	float &lep2chi2ndf()
	{
		if (not lep2chi2ndf_isLoaded) {
			if (lep2chi2ndf_branch != 0) {
				lep2chi2ndf_branch->GetEntry(index);
			} else { 
				printf("branch lep2chi2ndf_branch does not exist!\n");
				exit(1);
			}
			lep2chi2ndf_isLoaded = true;
		}
		return lep2chi2ndf_;
	}
	float &lep1dpt()
	{
		if (not lep1dpt_isLoaded) {
			if (lep1dpt_branch != 0) {
				lep1dpt_branch->GetEntry(index);
			} else { 
				printf("branch lep1dpt_branch does not exist!\n");
				exit(1);
			}
			lep1dpt_isLoaded = true;
		}
		return lep1dpt_;
	}
	float &lep2dpt()
	{
		if (not lep2dpt_isLoaded) {
			if (lep2dpt_branch != 0) {
				lep2dpt_branch->GetEntry(index);
			} else { 
				printf("branch lep2dpt_branch does not exist!\n");
				exit(1);
			}
			lep2dpt_isLoaded = true;
		}
		return lep2dpt_;
	}
	int &id1()
	{
		if (not id1_isLoaded) {
			if (id1_branch != 0) {
				id1_branch->GetEntry(index);
			} else { 
				printf("branch id1_branch does not exist!\n");
				exit(1);
			}
			id1_isLoaded = true;
		}
		return id1_;
	}
	int &id2()
	{
		if (not id2_isLoaded) {
			if (id2_branch != 0) {
				id2_branch->GetEntry(index);
			} else { 
				printf("branch id2_branch does not exist!\n");
				exit(1);
			}
			id2_isLoaded = true;
		}
		return id2_;
	}
	int &leptype1()
	{
		if (not leptype1_isLoaded) {
			if (leptype1_branch != 0) {
				leptype1_branch->GetEntry(index);
			} else { 
				printf("branch leptype1_branch does not exist!\n");
				exit(1);
			}
			leptype1_isLoaded = true;
		}
		return leptype1_;
	}
	int &leptype2()
	{
		if (not leptype2_isLoaded) {
			if (leptype2_branch != 0) {
				leptype2_branch->GetEntry(index);
			} else { 
				printf("branch leptype2_branch does not exist!\n");
				exit(1);
			}
			leptype2_isLoaded = true;
		}
		return leptype2_;
	}
	int &w1()
	{
		if (not w1_isLoaded) {
			if (w1_branch != 0) {
				w1_branch->GetEntry(index);
			} else { 
				printf("branch w1_branch does not exist!\n");
				exit(1);
			}
			w1_isLoaded = true;
		}
		return w1_;
	}
	int &w2()
	{
		if (not w2_isLoaded) {
			if (w2_branch != 0) {
				w2_branch->GetEntry(index);
			} else { 
				printf("branch w2_branch does not exist!\n");
				exit(1);
			}
			w2_isLoaded = true;
		}
		return w2_;
	}
	float &iso1()
	{
		if (not iso1_isLoaded) {
			if (iso1_branch != 0) {
				iso1_branch->GetEntry(index);
			} else { 
				printf("branch iso1_branch does not exist!\n");
				exit(1);
			}
			iso1_isLoaded = true;
		}
		return iso1_;
	}
	float &isont1()
	{
		if (not isont1_isLoaded) {
			if (isont1_branch != 0) {
				isont1_branch->GetEntry(index);
			} else { 
				printf("branch isont1_branch does not exist!\n");
				exit(1);
			}
			isont1_isLoaded = true;
		}
		return isont1_;
	}
	float &isopfold1()
	{
		if (not isopfold1_isLoaded) {
			if (isopfold1_branch != 0) {
				isopfold1_branch->GetEntry(index);
			} else { 
				printf("branch isopfold1_branch does not exist!\n");
				exit(1);
			}
			isopfold1_isLoaded = true;
		}
		return isopfold1_;
	}
	float &isopf1()
	{
		if (not isopf1_isLoaded) {
			if (isopf1_branch != 0) {
				isopf1_branch->GetEntry(index);
			} else { 
				printf("branch isopf1_branch does not exist!\n");
				exit(1);
			}
			isopf1_isLoaded = true;
		}
		return isopf1_;
	}
	float &etasc1()
	{
		if (not etasc1_isLoaded) {
			if (etasc1_branch != 0) {
				etasc1_branch->GetEntry(index);
			} else { 
				printf("branch etasc1_branch does not exist!\n");
				exit(1);
			}
			etasc1_isLoaded = true;
		}
		return etasc1_;
	}
	float &etasc2()
	{
		if (not etasc2_isLoaded) {
			if (etasc2_branch != 0) {
				etasc2_branch->GetEntry(index);
			} else { 
				printf("branch etasc2_branch does not exist!\n");
				exit(1);
			}
			etasc2_isLoaded = true;
		}
		return etasc2_;
	}
	float &eoverpin()
	{
		if (not eoverpin_isLoaded) {
			if (eoverpin_branch != 0) {
				eoverpin_branch->GetEntry(index);
			} else { 
				printf("branch eoverpin_branch does not exist!\n");
				exit(1);
			}
			eoverpin_isLoaded = true;
		}
		return eoverpin_;
	}
	float &eoverpout()
	{
		if (not eoverpout_isLoaded) {
			if (eoverpout_branch != 0) {
				eoverpout_branch->GetEntry(index);
			} else { 
				printf("branch eoverpout_branch does not exist!\n");
				exit(1);
			}
			eoverpout_isLoaded = true;
		}
		return eoverpout_;
	}
	float &dEtaIn()
	{
		if (not dEtaIn_isLoaded) {
			if (dEtaIn_branch != 0) {
				dEtaIn_branch->GetEntry(index);
			} else { 
				printf("branch dEtaIn_branch does not exist!\n");
				exit(1);
			}
			dEtaIn_isLoaded = true;
		}
		return dEtaIn_;
	}
	float &dPhiIn()
	{
		if (not dPhiIn_isLoaded) {
			if (dPhiIn_branch != 0) {
				dPhiIn_branch->GetEntry(index);
			} else { 
				printf("branch dPhiIn_branch does not exist!\n");
				exit(1);
			}
			dPhiIn_isLoaded = true;
		}
		return dPhiIn_;
	}
	float &sigmaIEtaIEta()
	{
		if (not sigmaIEtaIEta_isLoaded) {
			if (sigmaIEtaIEta_branch != 0) {
				sigmaIEtaIEta_branch->GetEntry(index);
			} else { 
				printf("branch sigmaIEtaIEta_branch does not exist!\n");
				exit(1);
			}
			sigmaIEtaIEta_isLoaded = true;
		}
		return sigmaIEtaIEta_;
	}
	float &hOverE()
	{
		if (not hOverE_isLoaded) {
			if (hOverE_branch != 0) {
				hOverE_branch->GetEntry(index);
			} else { 
				printf("branch hOverE_branch does not exist!\n");
				exit(1);
			}
			hOverE_isLoaded = true;
		}
		return hOverE_;
	}
	float &ooemoop()
	{
		if (not ooemoop_isLoaded) {
			if (ooemoop_branch != 0) {
				ooemoop_branch->GetEntry(index);
			} else { 
				printf("branch ooemoop_branch does not exist!\n");
				exit(1);
			}
			ooemoop_isLoaded = true;
		}
		return ooemoop_;
	}
	float &d0vtx()
	{
		if (not d0vtx_isLoaded) {
			if (d0vtx_branch != 0) {
				d0vtx_branch->GetEntry(index);
			} else { 
				printf("branch d0vtx_branch does not exist!\n");
				exit(1);
			}
			d0vtx_isLoaded = true;
		}
		return d0vtx_;
	}
	float &dzvtx()
	{
		if (not dzvtx_isLoaded) {
			if (dzvtx_branch != 0) {
				dzvtx_branch->GetEntry(index);
			} else { 
				printf("branch dzvtx_branch does not exist!\n");
				exit(1);
			}
			dzvtx_isLoaded = true;
		}
		return dzvtx_;
	}
	float &expinnerlayers()
	{
		if (not expinnerlayers_isLoaded) {
			if (expinnerlayers_branch != 0) {
				expinnerlayers_branch->GetEntry(index);
			} else { 
				printf("branch expinnerlayers_branch does not exist!\n");
				exit(1);
			}
			expinnerlayers_isLoaded = true;
		}
		return expinnerlayers_;
	}
	float &fbrem()
	{
		if (not fbrem_isLoaded) {
			if (fbrem_branch != 0) {
				fbrem_branch->GetEntry(index);
			} else { 
				printf("branch fbrem_branch does not exist!\n");
				exit(1);
			}
			fbrem_isLoaded = true;
		}
		return fbrem_;
	}
	float &pfisoch()
	{
		if (not pfisoch_isLoaded) {
			if (pfisoch_branch != 0) {
				pfisoch_branch->GetEntry(index);
			} else { 
				printf("branch pfisoch_branch does not exist!\n");
				exit(1);
			}
			pfisoch_isLoaded = true;
		}
		return pfisoch_;
	}
	float &pfisoem()
	{
		if (not pfisoem_isLoaded) {
			if (pfisoem_branch != 0) {
				pfisoem_branch->GetEntry(index);
			} else { 
				printf("branch pfisoem_branch does not exist!\n");
				exit(1);
			}
			pfisoem_isLoaded = true;
		}
		return pfisoem_;
	}
	float &pfisonh()
	{
		if (not pfisonh_isLoaded) {
			if (pfisonh_branch != 0) {
				pfisonh_branch->GetEntry(index);
			} else { 
				printf("branch pfisonh_branch does not exist!\n");
				exit(1);
			}
			pfisonh_isLoaded = true;
		}
		return pfisonh_;
	}
	float &eSC()
	{
		if (not eSC_isLoaded) {
			if (eSC_branch != 0) {
				eSC_branch->GetEntry(index);
			} else { 
				printf("branch eSC_branch does not exist!\n");
				exit(1);
			}
			eSC_isLoaded = true;
		}
		return eSC_;
	}
	float &phiSC()
	{
		if (not phiSC_isLoaded) {
			if (phiSC_branch != 0) {
				phiSC_branch->GetEntry(index);
			} else { 
				printf("branch phiSC_branch does not exist!\n");
				exit(1);
			}
			phiSC_isLoaded = true;
		}
		return phiSC_;
	}
	float &eSCRaw()
	{
		if (not eSCRaw_isLoaded) {
			if (eSCRaw_branch != 0) {
				eSCRaw_branch->GetEntry(index);
			} else { 
				printf("branch eSCRaw_branch does not exist!\n");
				exit(1);
			}
			eSCRaw_isLoaded = true;
		}
		return eSCRaw_;
	}
	float &eSCPresh()
	{
		if (not eSCPresh_isLoaded) {
			if (eSCPresh_branch != 0) {
				eSCPresh_branch->GetEntry(index);
			} else { 
				printf("branch eSCPresh_branch does not exist!\n");
				exit(1);
			}
			eSCPresh_isLoaded = true;
		}
		return eSCPresh_;
	}
	float &lep1_scslasercormean()
	{
		if (not lep1_scslasercormean_isLoaded) {
			if (lep1_scslasercormean_branch != 0) {
				lep1_scslasercormean_branch->GetEntry(index);
			} else { 
				printf("branch lep1_scslasercormean_branch does not exist!\n");
				exit(1);
			}
			lep1_scslasercormean_isLoaded = true;
		}
		return lep1_scslasercormean_;
	}
	float &lep1_scslasercormax()
	{
		if (not lep1_scslasercormax_isLoaded) {
			if (lep1_scslasercormax_branch != 0) {
				lep1_scslasercormax_branch->GetEntry(index);
			} else { 
				printf("branch lep1_scslasercormax_branch does not exist!\n");
				exit(1);
			}
			lep1_scslasercormax_isLoaded = true;
		}
		return lep1_scslasercormax_;
	}
	float &eoverpin2()
	{
		if (not eoverpin2_isLoaded) {
			if (eoverpin2_branch != 0) {
				eoverpin2_branch->GetEntry(index);
			} else { 
				printf("branch eoverpin2_branch does not exist!\n");
				exit(1);
			}
			eoverpin2_isLoaded = true;
		}
		return eoverpin2_;
	}
	float &eoverpout2()
	{
		if (not eoverpout2_isLoaded) {
			if (eoverpout2_branch != 0) {
				eoverpout2_branch->GetEntry(index);
			} else { 
				printf("branch eoverpout2_branch does not exist!\n");
				exit(1);
			}
			eoverpout2_isLoaded = true;
		}
		return eoverpout2_;
	}
	float &dEtaIn2()
	{
		if (not dEtaIn2_isLoaded) {
			if (dEtaIn2_branch != 0) {
				dEtaIn2_branch->GetEntry(index);
			} else { 
				printf("branch dEtaIn2_branch does not exist!\n");
				exit(1);
			}
			dEtaIn2_isLoaded = true;
		}
		return dEtaIn2_;
	}
	float &dPhiIn2()
	{
		if (not dPhiIn2_isLoaded) {
			if (dPhiIn2_branch != 0) {
				dPhiIn2_branch->GetEntry(index);
			} else { 
				printf("branch dPhiIn2_branch does not exist!\n");
				exit(1);
			}
			dPhiIn2_isLoaded = true;
		}
		return dPhiIn2_;
	}
	float &sigmaIEtaIEta2()
	{
		if (not sigmaIEtaIEta2_isLoaded) {
			if (sigmaIEtaIEta2_branch != 0) {
				sigmaIEtaIEta2_branch->GetEntry(index);
			} else { 
				printf("branch sigmaIEtaIEta2_branch does not exist!\n");
				exit(1);
			}
			sigmaIEtaIEta2_isLoaded = true;
		}
		return sigmaIEtaIEta2_;
	}
	float &hOverE2()
	{
		if (not hOverE2_isLoaded) {
			if (hOverE2_branch != 0) {
				hOverE2_branch->GetEntry(index);
			} else { 
				printf("branch hOverE2_branch does not exist!\n");
				exit(1);
			}
			hOverE2_isLoaded = true;
		}
		return hOverE2_;
	}
	float &ooemoop2()
	{
		if (not ooemoop2_isLoaded) {
			if (ooemoop2_branch != 0) {
				ooemoop2_branch->GetEntry(index);
			} else { 
				printf("branch ooemoop2_branch does not exist!\n");
				exit(1);
			}
			ooemoop2_isLoaded = true;
		}
		return ooemoop2_;
	}
	float &d0vtx2()
	{
		if (not d0vtx2_isLoaded) {
			if (d0vtx2_branch != 0) {
				d0vtx2_branch->GetEntry(index);
			} else { 
				printf("branch d0vtx2_branch does not exist!\n");
				exit(1);
			}
			d0vtx2_isLoaded = true;
		}
		return d0vtx2_;
	}
	float &dzvtx2()
	{
		if (not dzvtx2_isLoaded) {
			if (dzvtx2_branch != 0) {
				dzvtx2_branch->GetEntry(index);
			} else { 
				printf("branch dzvtx2_branch does not exist!\n");
				exit(1);
			}
			dzvtx2_isLoaded = true;
		}
		return dzvtx2_;
	}
	float &expinnerlayers2()
	{
		if (not expinnerlayers2_isLoaded) {
			if (expinnerlayers2_branch != 0) {
				expinnerlayers2_branch->GetEntry(index);
			} else { 
				printf("branch expinnerlayers2_branch does not exist!\n");
				exit(1);
			}
			expinnerlayers2_isLoaded = true;
		}
		return expinnerlayers2_;
	}
	float &fbrem2()
	{
		if (not fbrem2_isLoaded) {
			if (fbrem2_branch != 0) {
				fbrem2_branch->GetEntry(index);
			} else { 
				printf("branch fbrem2_branch does not exist!\n");
				exit(1);
			}
			fbrem2_isLoaded = true;
		}
		return fbrem2_;
	}
	float &pfisoch2()
	{
		if (not pfisoch2_isLoaded) {
			if (pfisoch2_branch != 0) {
				pfisoch2_branch->GetEntry(index);
			} else { 
				printf("branch pfisoch2_branch does not exist!\n");
				exit(1);
			}
			pfisoch2_isLoaded = true;
		}
		return pfisoch2_;
	}
	float &pfisoem2()
	{
		if (not pfisoem2_isLoaded) {
			if (pfisoem2_branch != 0) {
				pfisoem2_branch->GetEntry(index);
			} else { 
				printf("branch pfisoem2_branch does not exist!\n");
				exit(1);
			}
			pfisoem2_isLoaded = true;
		}
		return pfisoem2_;
	}
	float &pfisonh2()
	{
		if (not pfisonh2_isLoaded) {
			if (pfisonh2_branch != 0) {
				pfisonh2_branch->GetEntry(index);
			} else { 
				printf("branch pfisonh2_branch does not exist!\n");
				exit(1);
			}
			pfisonh2_isLoaded = true;
		}
		return pfisonh2_;
	}
	float &eSC2()
	{
		if (not eSC2_isLoaded) {
			if (eSC2_branch != 0) {
				eSC2_branch->GetEntry(index);
			} else { 
				printf("branch eSC2_branch does not exist!\n");
				exit(1);
			}
			eSC2_isLoaded = true;
		}
		return eSC2_;
	}
	float &phiSC2()
	{
		if (not phiSC2_isLoaded) {
			if (phiSC2_branch != 0) {
				phiSC2_branch->GetEntry(index);
			} else { 
				printf("branch phiSC2_branch does not exist!\n");
				exit(1);
			}
			phiSC2_isLoaded = true;
		}
		return phiSC2_;
	}
	float &eSCRaw2()
	{
		if (not eSCRaw2_isLoaded) {
			if (eSCRaw2_branch != 0) {
				eSCRaw2_branch->GetEntry(index);
			} else { 
				printf("branch eSCRaw2_branch does not exist!\n");
				exit(1);
			}
			eSCRaw2_isLoaded = true;
		}
		return eSCRaw2_;
	}
	float &eSCPresh2()
	{
		if (not eSCPresh2_isLoaded) {
			if (eSCPresh2_branch != 0) {
				eSCPresh2_branch->GetEntry(index);
			} else { 
				printf("branch eSCPresh2_branch does not exist!\n");
				exit(1);
			}
			eSCPresh2_isLoaded = true;
		}
		return eSCPresh2_;
	}
	float &lep2_scslasercormean()
	{
		if (not lep2_scslasercormean_isLoaded) {
			if (lep2_scslasercormean_branch != 0) {
				lep2_scslasercormean_branch->GetEntry(index);
			} else { 
				printf("branch lep2_scslasercormean_branch does not exist!\n");
				exit(1);
			}
			lep2_scslasercormean_isLoaded = true;
		}
		return lep2_scslasercormean_;
	}
	float &lep2_scslasercormax()
	{
		if (not lep2_scslasercormax_isLoaded) {
			if (lep2_scslasercormax_branch != 0) {
				lep2_scslasercormax_branch->GetEntry(index);
			} else { 
				printf("branch lep2_scslasercormax_branch does not exist!\n");
				exit(1);
			}
			lep2_scslasercormax_isLoaded = true;
		}
		return lep2_scslasercormax_;
	}
	float &scslasercormax()
	{
		if (not scslasercormax_isLoaded) {
			if (scslasercormax_branch != 0) {
				scslasercormax_branch->GetEntry(index);
			} else { 
				printf("branch scslasercormax_branch does not exist!\n");
				exit(1);
			}
			scslasercormax_isLoaded = true;
		}
		return scslasercormax_;
	}
	float &scslasercormax_pt()
	{
		if (not scslasercormax_pt_isLoaded) {
			if (scslasercormax_pt_branch != 0) {
				scslasercormax_pt_branch->GetEntry(index);
			} else { 
				printf("branch scslasercormax_pt_branch does not exist!\n");
				exit(1);
			}
			scslasercormax_pt_isLoaded = true;
		}
		return scslasercormax_pt_;
	}
	float &scslasercormax_eta()
	{
		if (not scslasercormax_eta_isLoaded) {
			if (scslasercormax_eta_branch != 0) {
				scslasercormax_eta_branch->GetEntry(index);
			} else { 
				printf("branch scslasercormax_eta_branch does not exist!\n");
				exit(1);
			}
			scslasercormax_eta_isLoaded = true;
		}
		return scslasercormax_eta_;
	}
	float &iso2()
	{
		if (not iso2_isLoaded) {
			if (iso2_branch != 0) {
				iso2_branch->GetEntry(index);
			} else { 
				printf("branch iso2_branch does not exist!\n");
				exit(1);
			}
			iso2_isLoaded = true;
		}
		return iso2_;
	}
	float &ecalveto1()
	{
		if (not ecalveto1_isLoaded) {
			if (ecalveto1_branch != 0) {
				ecalveto1_branch->GetEntry(index);
			} else { 
				printf("branch ecalveto1_branch does not exist!\n");
				exit(1);
			}
			ecalveto1_isLoaded = true;
		}
		return ecalveto1_;
	}
	float &ecalveto2()
	{
		if (not ecalveto2_isLoaded) {
			if (ecalveto2_branch != 0) {
				ecalveto2_branch->GetEntry(index);
			} else { 
				printf("branch ecalveto2_branch does not exist!\n");
				exit(1);
			}
			ecalveto2_isLoaded = true;
		}
		return ecalveto2_;
	}
	float &hcalveto1()
	{
		if (not hcalveto1_isLoaded) {
			if (hcalveto1_branch != 0) {
				hcalveto1_branch->GetEntry(index);
			} else { 
				printf("branch hcalveto1_branch does not exist!\n");
				exit(1);
			}
			hcalveto1_isLoaded = true;
		}
		return hcalveto1_;
	}
	float &hcalveto2()
	{
		if (not hcalveto2_isLoaded) {
			if (hcalveto2_branch != 0) {
				hcalveto2_branch->GetEntry(index);
			} else { 
				printf("branch hcalveto2_branch does not exist!\n");
				exit(1);
			}
			hcalveto2_isLoaded = true;
		}
		return hcalveto2_;
	}
	float &isont2()
	{
		if (not isont2_isLoaded) {
			if (isont2_branch != 0) {
				isont2_branch->GetEntry(index);
			} else { 
				printf("branch isont2_branch does not exist!\n");
				exit(1);
			}
			isont2_isLoaded = true;
		}
		return isont2_;
	}
	float &isopf2()
	{
		if (not isopf2_isLoaded) {
			if (isopf2_branch != 0) {
				isopf2_branch->GetEntry(index);
			} else { 
				printf("branch isopf2_branch does not exist!\n");
				exit(1);
			}
			isopf2_isLoaded = true;
		}
		return isopf2_;
	}
	float &ptl1()
	{
		if (not ptl1_isLoaded) {
			if (ptl1_branch != 0) {
				ptl1_branch->GetEntry(index);
			} else { 
				printf("branch ptl1_branch does not exist!\n");
				exit(1);
			}
			ptl1_isLoaded = true;
		}
		return ptl1_;
	}
	float &ptl2()
	{
		if (not ptl2_isLoaded) {
			if (ptl2_branch != 0) {
				ptl2_branch->GetEntry(index);
			} else { 
				printf("branch ptl2_branch does not exist!\n");
				exit(1);
			}
			ptl2_isLoaded = true;
		}
		return ptl2_;
	}
	float &etal1()
	{
		if (not etal1_isLoaded) {
			if (etal1_branch != 0) {
				etal1_branch->GetEntry(index);
			} else { 
				printf("branch etal1_branch does not exist!\n");
				exit(1);
			}
			etal1_isLoaded = true;
		}
		return etal1_;
	}
	float &etal2()
	{
		if (not etal2_isLoaded) {
			if (etal2_branch != 0) {
				etal2_branch->GetEntry(index);
			} else { 
				printf("branch etal2_branch does not exist!\n");
				exit(1);
			}
			etal2_isLoaded = true;
		}
		return etal2_;
	}
	float &phil1()
	{
		if (not phil1_isLoaded) {
			if (phil1_branch != 0) {
				phil1_branch->GetEntry(index);
			} else { 
				printf("branch phil1_branch does not exist!\n");
				exit(1);
			}
			phil1_isLoaded = true;
		}
		return phil1_;
	}
	float &phil2()
	{
		if (not phil2_isLoaded) {
			if (phil2_branch != 0) {
				phil2_branch->GetEntry(index);
			} else { 
				printf("branch phil2_branch does not exist!\n");
				exit(1);
			}
			phil2_isLoaded = true;
		}
		return phil2_;
	}
	float &meff()
	{
		if (not meff_isLoaded) {
			if (meff_branch != 0) {
				meff_branch->GetEntry(index);
			} else { 
				printf("branch meff_branch does not exist!\n");
				exit(1);
			}
			meff_isLoaded = true;
		}
		return meff_;
	}
	float &mt()
	{
		if (not mt_isLoaded) {
			if (mt_branch != 0) {
				mt_branch->GetEntry(index);
			} else { 
				printf("branch mt_branch does not exist!\n");
				exit(1);
			}
			mt_isLoaded = true;
		}
		return mt_;
	}
	unsigned int &run()
	{
		if (not run_isLoaded) {
			if (run_branch != 0) {
				run_branch->GetEntry(index);
			} else { 
				printf("branch run_branch does not exist!\n");
				exit(1);
			}
			run_isLoaded = true;
		}
		return run_;
	}
	unsigned int &lumi()
	{
		if (not lumi_isLoaded) {
			if (lumi_branch != 0) {
				lumi_branch->GetEntry(index);
			} else { 
				printf("branch lumi_branch does not exist!\n");
				exit(1);
			}
			lumi_isLoaded = true;
		}
		return lumi_;
	}
	unsigned int &event()
	{
		if (not event_isLoaded) {
			if (event_branch != 0) {
				event_branch->GetEntry(index);
			} else { 
				printf("branch event_branch does not exist!\n");
				exit(1);
			}
			event_isLoaded = true;
		}
		return event_;
	}
	float &y()
	{
		if (not y_isLoaded) {
			if (y_branch != 0) {
				y_branch->GetEntry(index);
			} else { 
				printf("branch y_branch does not exist!\n");
				exit(1);
			}
			y_isLoaded = true;
		}
		return y_;
	}
	float &ht()
	{
		if (not ht_isLoaded) {
			if (ht_branch != 0) {
				ht_branch->GetEntry(index);
			} else { 
				printf("branch ht_branch does not exist!\n");
				exit(1);
			}
			ht_isLoaded = true;
		}
		return ht_;
	}
	float &htgen()
	{
		if (not htgen_isLoaded) {
			if (htgen_branch != 0) {
				htgen_branch->GetEntry(index);
			} else { 
				printf("branch htgen_branch does not exist!\n");
				exit(1);
			}
			htgen_isLoaded = true;
		}
		return htgen_;
	}
	float &htjpt()
	{
		if (not htjpt_isLoaded) {
			if (htjpt_branch != 0) {
				htjpt_branch->GetEntry(index);
			} else { 
				printf("branch htjpt_branch does not exist!\n");
				exit(1);
			}
			htjpt_isLoaded = true;
		}
		return htjpt_;
	}
	int &nels()
	{
		if (not nels_isLoaded) {
			if (nels_branch != 0) {
				nels_branch->GetEntry(index);
			} else { 
				printf("branch nels_branch does not exist!\n");
				exit(1);
			}
			nels_isLoaded = true;
		}
		return nels_;
	}
	int &nmus()
	{
		if (not nmus_isLoaded) {
			if (nmus_branch != 0) {
				nmus_branch->GetEntry(index);
			} else { 
				printf("branch nmus_branch does not exist!\n");
				exit(1);
			}
			nmus_isLoaded = true;
		}
		return nmus_;
	}
	int &ntaus()
	{
		if (not ntaus_isLoaded) {
			if (ntaus_branch != 0) {
				ntaus_branch->GetEntry(index);
			} else { 
				printf("branch ntaus_branch does not exist!\n");
				exit(1);
			}
			ntaus_isLoaded = true;
		}
		return ntaus_;
	}
	int &nleps()
	{
		if (not nleps_isLoaded) {
			if (nleps_branch != 0) {
				nleps_branch->GetEntry(index);
			} else { 
				printf("branch nleps_branch does not exist!\n");
				exit(1);
			}
			nleps_isLoaded = true;
		}
		return nleps_;
	}
	int &nbs()
	{
		if (not nbs_isLoaded) {
			if (nbs_branch != 0) {
				nbs_branch->GetEntry(index);
			} else { 
				printf("branch nbs_branch does not exist!\n");
				exit(1);
			}
			nbs_isLoaded = true;
		}
		return nbs_;
	}
	float &dphijm()
	{
		if (not dphijm_isLoaded) {
			if (dphijm_branch != 0) {
				dphijm_branch->GetEntry(index);
			} else { 
				printf("branch dphijm_branch does not exist!\n");
				exit(1);
			}
			dphijm_isLoaded = true;
		}
		return dphijm_;
	}
	float &ptjetraw()
	{
		if (not ptjetraw_isLoaded) {
			if (ptjetraw_branch != 0) {
				ptjetraw_branch->GetEntry(index);
			} else { 
				printf("branch ptjetraw_branch does not exist!\n");
				exit(1);
			}
			ptjetraw_isLoaded = true;
		}
		return ptjetraw_;
	}
	float &ptjet23()
	{
		if (not ptjet23_isLoaded) {
			if (ptjet23_branch != 0) {
				ptjet23_branch->GetEntry(index);
			} else { 
				printf("branch ptjet23_branch does not exist!\n");
				exit(1);
			}
			ptjet23_isLoaded = true;
		}
		return ptjet23_;
	}
	float &ptjetF23()
	{
		if (not ptjetF23_isLoaded) {
			if (ptjetF23_branch != 0) {
				ptjetF23_branch->GetEntry(index);
			} else { 
				printf("branch ptjetF23_branch does not exist!\n");
				exit(1);
			}
			ptjetF23_isLoaded = true;
		}
		return ptjetF23_;
	}
	float &ptjetO23()
	{
		if (not ptjetO23_isLoaded) {
			if (ptjetO23_branch != 0) {
				ptjetO23_branch->GetEntry(index);
			} else { 
				printf("branch ptjetO23_branch does not exist!\n");
				exit(1);
			}
			ptjetO23_isLoaded = true;
		}
		return ptjetO23_;
	}
	int &mcid1()
	{
		if (not mcid1_isLoaded) {
			if (mcid1_branch != 0) {
				mcid1_branch->GetEntry(index);
			} else { 
				printf("branch mcid1_branch does not exist!\n");
				exit(1);
			}
			mcid1_isLoaded = true;
		}
		return mcid1_;
	}
	float &mcdr1()
	{
		if (not mcdr1_isLoaded) {
			if (mcdr1_branch != 0) {
				mcdr1_branch->GetEntry(index);
			} else { 
				printf("branch mcdr1_branch does not exist!\n");
				exit(1);
			}
			mcdr1_isLoaded = true;
		}
		return mcdr1_;
	}
	int &mcdecay1()
	{
		if (not mcdecay1_isLoaded) {
			if (mcdecay1_branch != 0) {
				mcdecay1_branch->GetEntry(index);
			} else { 
				printf("branch mcdecay1_branch does not exist!\n");
				exit(1);
			}
			mcdecay1_isLoaded = true;
		}
		return mcdecay1_;
	}
	int &mcndec1()
	{
		if (not mcndec1_isLoaded) {
			if (mcndec1_branch != 0) {
				mcndec1_branch->GetEntry(index);
			} else { 
				printf("branch mcndec1_branch does not exist!\n");
				exit(1);
			}
			mcndec1_isLoaded = true;
		}
		return mcndec1_;
	}
	int &mcndec2()
	{
		if (not mcndec2_isLoaded) {
			if (mcndec2_branch != 0) {
				mcndec2_branch->GetEntry(index);
			} else { 
				printf("branch mcndec2_branch does not exist!\n");
				exit(1);
			}
			mcndec2_isLoaded = true;
		}
		return mcndec2_;
	}
	int &mcndeckls1()
	{
		if (not mcndeckls1_isLoaded) {
			if (mcndeckls1_branch != 0) {
				mcndeckls1_branch->GetEntry(index);
			} else { 
				printf("branch mcndeckls1_branch does not exist!\n");
				exit(1);
			}
			mcndeckls1_isLoaded = true;
		}
		return mcndeckls1_;
	}
	int &mcndeckls2()
	{
		if (not mcndeckls2_isLoaded) {
			if (mcndeckls2_branch != 0) {
				mcndeckls2_branch->GetEntry(index);
			} else { 
				printf("branch mcndeckls2_branch does not exist!\n");
				exit(1);
			}
			mcndeckls2_isLoaded = true;
		}
		return mcndeckls2_;
	}
	int &mcndecem1()
	{
		if (not mcndecem1_isLoaded) {
			if (mcndecem1_branch != 0) {
				mcndecem1_branch->GetEntry(index);
			} else { 
				printf("branch mcndecem1_branch does not exist!\n");
				exit(1);
			}
			mcndecem1_isLoaded = true;
		}
		return mcndecem1_;
	}
	int &mcndecem2()
	{
		if (not mcndecem2_isLoaded) {
			if (mcndecem2_branch != 0) {
				mcndecem2_branch->GetEntry(index);
			} else { 
				printf("branch mcndecem2_branch does not exist!\n");
				exit(1);
			}
			mcndecem2_isLoaded = true;
		}
		return mcndecem2_;
	}
	int &mcid2()
	{
		if (not mcid2_isLoaded) {
			if (mcid2_branch != 0) {
				mcid2_branch->GetEntry(index);
			} else { 
				printf("branch mcid2_branch does not exist!\n");
				exit(1);
			}
			mcid2_isLoaded = true;
		}
		return mcid2_;
	}
	float &mcdr2()
	{
		if (not mcdr2_isLoaded) {
			if (mcdr2_branch != 0) {
				mcdr2_branch->GetEntry(index);
			} else { 
				printf("branch mcdr2_branch does not exist!\n");
				exit(1);
			}
			mcdr2_isLoaded = true;
		}
		return mcdr2_;
	}
	int &mcdecay2()
	{
		if (not mcdecay2_isLoaded) {
			if (mcdecay2_branch != 0) {
				mcdecay2_branch->GetEntry(index);
			} else { 
				printf("branch mcdecay2_branch does not exist!\n");
				exit(1);
			}
			mcdecay2_isLoaded = true;
		}
		return mcdecay2_;
	}
	float &mctaudpt1()
	{
		if (not mctaudpt1_isLoaded) {
			if (mctaudpt1_branch != 0) {
				mctaudpt1_branch->GetEntry(index);
			} else { 
				printf("branch mctaudpt1_branch does not exist!\n");
				exit(1);
			}
			mctaudpt1_isLoaded = true;
		}
		return mctaudpt1_;
	}
	float &mctaudpt2()
	{
		if (not mctaudpt2_isLoaded) {
			if (mctaudpt2_branch != 0) {
				mctaudpt2_branch->GetEntry(index);
			} else { 
				printf("branch mctaudpt2_branch does not exist!\n");
				exit(1);
			}
			mctaudpt2_isLoaded = true;
		}
		return mctaudpt2_;
	}
	int &mctaudid1()
	{
		if (not mctaudid1_isLoaded) {
			if (mctaudid1_branch != 0) {
				mctaudid1_branch->GetEntry(index);
			} else { 
				printf("branch mctaudid1_branch does not exist!\n");
				exit(1);
			}
			mctaudid1_isLoaded = true;
		}
		return mctaudid1_;
	}
	int &mctaudid2()
	{
		if (not mctaudid2_isLoaded) {
			if (mctaudid2_branch != 0) {
				mctaudid2_branch->GetEntry(index);
			} else { 
				printf("branch mctaudid2_branch does not exist!\n");
				exit(1);
			}
			mctaudid2_isLoaded = true;
		}
		return mctaudid2_;
	}
	int &mlepid()
	{
		if (not mlepid_isLoaded) {
			if (mlepid_branch != 0) {
				mlepid_branch->GetEntry(index);
			} else { 
				printf("branch mlepid_branch does not exist!\n");
				exit(1);
			}
			mlepid_isLoaded = true;
		}
		return mlepid_;
	}
	int &mleppassid()
	{
		if (not mleppassid_isLoaded) {
			if (mleppassid_branch != 0) {
				mleppassid_branch->GetEntry(index);
			} else { 
				printf("branch mleppassid_branch does not exist!\n");
				exit(1);
			}
			mleppassid_isLoaded = true;
		}
		return mleppassid_;
	}
	int &mleppassiso()
	{
		if (not mleppassiso_isLoaded) {
			if (mleppassiso_branch != 0) {
				mleppassiso_branch->GetEntry(index);
			} else { 
				printf("branch mleppassiso_branch does not exist!\n");
				exit(1);
			}
			mleppassiso_isLoaded = true;
		}
		return mleppassiso_;
	}
	float &mlepiso()
	{
		if (not mlepiso_isLoaded) {
			if (mlepiso_branch != 0) {
				mlepiso_branch->GetEntry(index);
			} else { 
				printf("branch mlepiso_branch does not exist!\n");
				exit(1);
			}
			mlepiso_isLoaded = true;
		}
		return mlepiso_;
	}
	float &mlepdr()
	{
		if (not mlepdr_isLoaded) {
			if (mlepdr_branch != 0) {
				mlepdr_branch->GetEntry(index);
			} else { 
				printf("branch mlepdr_branch does not exist!\n");
				exit(1);
			}
			mlepdr_isLoaded = true;
		}
		return mlepdr_;
	}
	float &pflepiso()
	{
		if (not pflepiso_isLoaded) {
			if (pflepiso_branch != 0) {
				pflepiso_branch->GetEntry(index);
			} else { 
				printf("branch pflepiso_branch does not exist!\n");
				exit(1);
			}
			pflepiso_isLoaded = true;
		}
		return pflepiso_;
	}
	float &pflepdr()
	{
		if (not pflepdr_isLoaded) {
			if (pflepdr_branch != 0) {
				pflepdr_branch->GetEntry(index);
			} else { 
				printf("branch pflepdr_branch does not exist!\n");
				exit(1);
			}
			pflepdr_isLoaded = true;
		}
		return pflepdr_;
	}
	float &pfleppt()
	{
		if (not pfleppt_isLoaded) {
			if (pfleppt_branch != 0) {
				pfleppt_branch->GetEntry(index);
			} else { 
				printf("branch pfleppt_branch does not exist!\n");
				exit(1);
			}
			pfleppt_isLoaded = true;
		}
		return pfleppt_;
	}
	float &pflepmindrj()
	{
		if (not pflepmindrj_isLoaded) {
			if (pflepmindrj_branch != 0) {
				pflepmindrj_branch->GetEntry(index);
			} else { 
				printf("branch pflepmindrj_branch does not exist!\n");
				exit(1);
			}
			pflepmindrj_isLoaded = true;
		}
		return pflepmindrj_;
	}
	float &pftaudiso()
	{
		if (not pftaudiso_isLoaded) {
			if (pftaudiso_branch != 0) {
				pftaudiso_branch->GetEntry(index);
			} else { 
				printf("branch pftaudiso_branch does not exist!\n");
				exit(1);
			}
			pftaudiso_isLoaded = true;
		}
		return pftaudiso_;
	}
	float &pftauddr()
	{
		if (not pftauddr_isLoaded) {
			if (pftauddr_branch != 0) {
				pftauddr_branch->GetEntry(index);
			} else { 
				printf("branch pftauddr_branch does not exist!\n");
				exit(1);
			}
			pftauddr_isLoaded = true;
		}
		return pftauddr_;
	}
	float &pftaudpt()
	{
		if (not pftaudpt_isLoaded) {
			if (pftaudpt_branch != 0) {
				pftaudpt_branch->GetEntry(index);
			} else { 
				printf("branch pftaudpt_branch does not exist!\n");
				exit(1);
			}
			pftaudpt_isLoaded = true;
		}
		return pftaudpt_;
	}
	float &pftaudmindrj()
	{
		if (not pftaudmindrj_isLoaded) {
			if (pftaudmindrj_branch != 0) {
				pftaudmindrj_branch->GetEntry(index);
			} else { 
				printf("branch pftaudmindrj_branch does not exist!\n");
				exit(1);
			}
			pftaudmindrj_isLoaded = true;
		}
		return pftaudmindrj_;
	}
	int &pfcandid5()
	{
		if (not pfcandid5_isLoaded) {
			if (pfcandid5_branch != 0) {
				pfcandid5_branch->GetEntry(index);
			} else { 
				printf("branch pfcandid5_branch does not exist!\n");
				exit(1);
			}
			pfcandid5_isLoaded = true;
		}
		return pfcandid5_;
	}
	float &pfcandiso5()
	{
		if (not pfcandiso5_isLoaded) {
			if (pfcandiso5_branch != 0) {
				pfcandiso5_branch->GetEntry(index);
			} else { 
				printf("branch pfcandiso5_branch does not exist!\n");
				exit(1);
			}
			pfcandiso5_isLoaded = true;
		}
		return pfcandiso5_;
	}
	float &pfcandpt5()
	{
		if (not pfcandpt5_isLoaded) {
			if (pfcandpt5_branch != 0) {
				pfcandpt5_branch->GetEntry(index);
			} else { 
				printf("branch pfcandpt5_branch does not exist!\n");
				exit(1);
			}
			pfcandpt5_isLoaded = true;
		}
		return pfcandpt5_;
	}
	float &pfcanddz5()
	{
		if (not pfcanddz5_isLoaded) {
			if (pfcanddz5_branch != 0) {
				pfcanddz5_branch->GetEntry(index);
			} else { 
				printf("branch pfcanddz5_branch does not exist!\n");
				exit(1);
			}
			pfcanddz5_isLoaded = true;
		}
		return pfcanddz5_;
	}
	float &pfcandmindrj5()
	{
		if (not pfcandmindrj5_isLoaded) {
			if (pfcandmindrj5_branch != 0) {
				pfcandmindrj5_branch->GetEntry(index);
			} else { 
				printf("branch pfcandmindrj5_branch does not exist!\n");
				exit(1);
			}
			pfcandmindrj5_isLoaded = true;
		}
		return pfcandmindrj5_;
	}
	int &pfcandid10()
	{
		if (not pfcandid10_isLoaded) {
			if (pfcandid10_branch != 0) {
				pfcandid10_branch->GetEntry(index);
			} else { 
				printf("branch pfcandid10_branch does not exist!\n");
				exit(1);
			}
			pfcandid10_isLoaded = true;
		}
		return pfcandid10_;
	}
	float &pfcandiso10()
	{
		if (not pfcandiso10_isLoaded) {
			if (pfcandiso10_branch != 0) {
				pfcandiso10_branch->GetEntry(index);
			} else { 
				printf("branch pfcandiso10_branch does not exist!\n");
				exit(1);
			}
			pfcandiso10_isLoaded = true;
		}
		return pfcandiso10_;
	}
	float &pfcandpt10()
	{
		if (not pfcandpt10_isLoaded) {
			if (pfcandpt10_branch != 0) {
				pfcandpt10_branch->GetEntry(index);
			} else { 
				printf("branch pfcandpt10_branch does not exist!\n");
				exit(1);
			}
			pfcandpt10_isLoaded = true;
		}
		return pfcandpt10_;
	}
	float &pfcanddz10()
	{
		if (not pfcanddz10_isLoaded) {
			if (pfcanddz10_branch != 0) {
				pfcanddz10_branch->GetEntry(index);
			} else { 
				printf("branch pfcanddz10_branch does not exist!\n");
				exit(1);
			}
			pfcanddz10_isLoaded = true;
		}
		return pfcanddz10_;
	}
	float &pfcandmindrj10()
	{
		if (not pfcandmindrj10_isLoaded) {
			if (pfcandmindrj10_branch != 0) {
				pfcandmindrj10_branch->GetEntry(index);
			} else { 
				printf("branch pfcandmindrj10_branch does not exist!\n");
				exit(1);
			}
			pfcandmindrj10_isLoaded = true;
		}
		return pfcandmindrj10_;
	}
	int &pfcandidOS10()
	{
		if (not pfcandidOS10_isLoaded) {
			if (pfcandidOS10_branch != 0) {
				pfcandidOS10_branch->GetEntry(index);
			} else { 
				printf("branch pfcandidOS10_branch does not exist!\n");
				exit(1);
			}
			pfcandidOS10_isLoaded = true;
		}
		return pfcandidOS10_;
	}
	float &pfcandisoOS10()
	{
		if (not pfcandisoOS10_isLoaded) {
			if (pfcandisoOS10_branch != 0) {
				pfcandisoOS10_branch->GetEntry(index);
			} else { 
				printf("branch pfcandisoOS10_branch does not exist!\n");
				exit(1);
			}
			pfcandisoOS10_isLoaded = true;
		}
		return pfcandisoOS10_;
	}
	float &pfcandptOS10()
	{
		if (not pfcandptOS10_isLoaded) {
			if (pfcandptOS10_branch != 0) {
				pfcandptOS10_branch->GetEntry(index);
			} else { 
				printf("branch pfcandptOS10_branch does not exist!\n");
				exit(1);
			}
			pfcandptOS10_isLoaded = true;
		}
		return pfcandptOS10_;
	}
	float &pfcanddzOS10()
	{
		if (not pfcanddzOS10_isLoaded) {
			if (pfcanddzOS10_branch != 0) {
				pfcanddzOS10_branch->GetEntry(index);
			} else { 
				printf("branch pfcanddzOS10_branch does not exist!\n");
				exit(1);
			}
			pfcanddzOS10_isLoaded = true;
		}
		return pfcanddzOS10_;
	}
	int &pfcandid5looseZ()
	{
		if (not pfcandid5looseZ_isLoaded) {
			if (pfcandid5looseZ_branch != 0) {
				pfcandid5looseZ_branch->GetEntry(index);
			} else { 
				printf("branch pfcandid5looseZ_branch does not exist!\n");
				exit(1);
			}
			pfcandid5looseZ_isLoaded = true;
		}
		return pfcandid5looseZ_;
	}
	float &pfcandiso5looseZ()
	{
		if (not pfcandiso5looseZ_isLoaded) {
			if (pfcandiso5looseZ_branch != 0) {
				pfcandiso5looseZ_branch->GetEntry(index);
			} else { 
				printf("branch pfcandiso5looseZ_branch does not exist!\n");
				exit(1);
			}
			pfcandiso5looseZ_isLoaded = true;
		}
		return pfcandiso5looseZ_;
	}
	float &pfcandpt5looseZ()
	{
		if (not pfcandpt5looseZ_isLoaded) {
			if (pfcandpt5looseZ_branch != 0) {
				pfcandpt5looseZ_branch->GetEntry(index);
			} else { 
				printf("branch pfcandpt5looseZ_branch does not exist!\n");
				exit(1);
			}
			pfcandpt5looseZ_isLoaded = true;
		}
		return pfcandpt5looseZ_;
	}
	float &pfcanddz5looseZ()
	{
		if (not pfcanddz5looseZ_isLoaded) {
			if (pfcanddz5looseZ_branch != 0) {
				pfcanddz5looseZ_branch->GetEntry(index);
			} else { 
				printf("branch pfcanddz5looseZ_branch does not exist!\n");
				exit(1);
			}
			pfcanddz5looseZ_isLoaded = true;
		}
		return pfcanddz5looseZ_;
	}
	int &pfcandidOS10looseZ()
	{
		if (not pfcandidOS10looseZ_isLoaded) {
			if (pfcandidOS10looseZ_branch != 0) {
				pfcandidOS10looseZ_branch->GetEntry(index);
			} else { 
				printf("branch pfcandidOS10looseZ_branch does not exist!\n");
				exit(1);
			}
			pfcandidOS10looseZ_isLoaded = true;
		}
		return pfcandidOS10looseZ_;
	}
	float &pfcandisoOS10looseZ()
	{
		if (not pfcandisoOS10looseZ_isLoaded) {
			if (pfcandisoOS10looseZ_branch != 0) {
				pfcandisoOS10looseZ_branch->GetEntry(index);
			} else { 
				printf("branch pfcandisoOS10looseZ_branch does not exist!\n");
				exit(1);
			}
			pfcandisoOS10looseZ_isLoaded = true;
		}
		return pfcandisoOS10looseZ_;
	}
	float &pfcandptOS10looseZ()
	{
		if (not pfcandptOS10looseZ_isLoaded) {
			if (pfcandptOS10looseZ_branch != 0) {
				pfcandptOS10looseZ_branch->GetEntry(index);
			} else { 
				printf("branch pfcandptOS10looseZ_branch does not exist!\n");
				exit(1);
			}
			pfcandptOS10looseZ_isLoaded = true;
		}
		return pfcandptOS10looseZ_;
	}
	float &pfcanddzOS10looseZ()
	{
		if (not pfcanddzOS10looseZ_isLoaded) {
			if (pfcanddzOS10looseZ_branch != 0) {
				pfcanddzOS10looseZ_branch->GetEntry(index);
			} else { 
				printf("branch pfcanddzOS10looseZ_branch does not exist!\n");
				exit(1);
			}
			pfcanddzOS10looseZ_isLoaded = true;
		}
		return pfcanddzOS10looseZ_;
	}
	int &pfcanddirid10()
	{
		if (not pfcanddirid10_isLoaded) {
			if (pfcanddirid10_branch != 0) {
				pfcanddirid10_branch->GetEntry(index);
			} else { 
				printf("branch pfcanddirid10_branch does not exist!\n");
				exit(1);
			}
			pfcanddirid10_isLoaded = true;
		}
		return pfcanddirid10_;
	}
	float &pfcanddiriso10()
	{
		if (not pfcanddiriso10_isLoaded) {
			if (pfcanddiriso10_branch != 0) {
				pfcanddiriso10_branch->GetEntry(index);
			} else { 
				printf("branch pfcanddiriso10_branch does not exist!\n");
				exit(1);
			}
			pfcanddiriso10_isLoaded = true;
		}
		return pfcanddiriso10_;
	}
	float &pfcanddirpt10()
	{
		if (not pfcanddirpt10_isLoaded) {
			if (pfcanddirpt10_branch != 0) {
				pfcanddirpt10_branch->GetEntry(index);
			} else { 
				printf("branch pfcanddirpt10_branch does not exist!\n");
				exit(1);
			}
			pfcanddirpt10_isLoaded = true;
		}
		return pfcanddirpt10_;
	}
	float &pfcanddirmindrj10()
	{
		if (not pfcanddirmindrj10_isLoaded) {
			if (pfcanddirmindrj10_branch != 0) {
				pfcanddirmindrj10_branch->GetEntry(index);
			} else { 
				printf("branch pfcanddirmindrj10_branch does not exist!\n");
				exit(1);
			}
			pfcanddirmindrj10_isLoaded = true;
		}
		return pfcanddirmindrj10_;
	}
	int &pfcandvetoid10()
	{
		if (not pfcandvetoid10_isLoaded) {
			if (pfcandvetoid10_branch != 0) {
				pfcandvetoid10_branch->GetEntry(index);
			} else { 
				printf("branch pfcandvetoid10_branch does not exist!\n");
				exit(1);
			}
			pfcandvetoid10_isLoaded = true;
		}
		return pfcandvetoid10_;
	}
	float &pfcandvetoiso10()
	{
		if (not pfcandvetoiso10_isLoaded) {
			if (pfcandvetoiso10_branch != 0) {
				pfcandvetoiso10_branch->GetEntry(index);
			} else { 
				printf("branch pfcandvetoiso10_branch does not exist!\n");
				exit(1);
			}
			pfcandvetoiso10_isLoaded = true;
		}
		return pfcandvetoiso10_;
	}
	float &pfcandvetopt10()
	{
		if (not pfcandvetopt10_isLoaded) {
			if (pfcandvetopt10_branch != 0) {
				pfcandvetopt10_branch->GetEntry(index);
			} else { 
				printf("branch pfcandvetopt10_branch does not exist!\n");
				exit(1);
			}
			pfcandvetopt10_isLoaded = true;
		}
		return pfcandvetopt10_;
	}
	float &pfcandvetomindrj10()
	{
		if (not pfcandvetomindrj10_isLoaded) {
			if (pfcandvetomindrj10_branch != 0) {
				pfcandvetomindrj10_branch->GetEntry(index);
			} else { 
				printf("branch pfcandvetomindrj10_branch does not exist!\n");
				exit(1);
			}
			pfcandvetomindrj10_isLoaded = true;
		}
		return pfcandvetomindrj10_;
	}
	int &pfcandvetoLid10()
	{
		if (not pfcandvetoLid10_isLoaded) {
			if (pfcandvetoLid10_branch != 0) {
				pfcandvetoLid10_branch->GetEntry(index);
			} else { 
				printf("branch pfcandvetoLid10_branch does not exist!\n");
				exit(1);
			}
			pfcandvetoLid10_isLoaded = true;
		}
		return pfcandvetoLid10_;
	}
	float &pfcandvetoLiso10()
	{
		if (not pfcandvetoLiso10_isLoaded) {
			if (pfcandvetoLiso10_branch != 0) {
				pfcandvetoLiso10_branch->GetEntry(index);
			} else { 
				printf("branch pfcandvetoLiso10_branch does not exist!\n");
				exit(1);
			}
			pfcandvetoLiso10_isLoaded = true;
		}
		return pfcandvetoLiso10_;
	}
	float &pfcandvetoLpt10()
	{
		if (not pfcandvetoLpt10_isLoaded) {
			if (pfcandvetoLpt10_branch != 0) {
				pfcandvetoLpt10_branch->GetEntry(index);
			} else { 
				printf("branch pfcandvetoLpt10_branch does not exist!\n");
				exit(1);
			}
			pfcandvetoLpt10_isLoaded = true;
		}
		return pfcandvetoLpt10_;
	}
	float &pfcandvetoLmindrj10()
	{
		if (not pfcandvetoLmindrj10_isLoaded) {
			if (pfcandvetoLmindrj10_branch != 0) {
				pfcandvetoLmindrj10_branch->GetEntry(index);
			} else { 
				printf("branch pfcandvetoLmindrj10_branch does not exist!\n");
				exit(1);
			}
			pfcandvetoLmindrj10_isLoaded = true;
		}
		return pfcandvetoLmindrj10_;
	}
	float &emjet10()
	{
		if (not emjet10_isLoaded) {
			if (emjet10_branch != 0) {
				emjet10_branch->GetEntry(index);
			} else { 
				printf("branch emjet10_branch does not exist!\n");
				exit(1);
			}
			emjet10_isLoaded = true;
		}
		return emjet10_;
	}
	float &mjj()
	{
		if (not mjj_isLoaded) {
			if (mjj_branch != 0) {
				mjj_branch->GetEntry(index);
			} else { 
				printf("branch mjj_branch does not exist!\n");
				exit(1);
			}
			mjj_isLoaded = true;
		}
		return mjj_;
	}
	float &emjet20()
	{
		if (not emjet20_isLoaded) {
			if (emjet20_branch != 0) {
				emjet20_branch->GetEntry(index);
			} else { 
				printf("branch emjet20_branch does not exist!\n");
				exit(1);
			}
			emjet20_isLoaded = true;
		}
		return emjet20_;
	}
	float &trkpt5()
	{
		if (not trkpt5_isLoaded) {
			if (trkpt5_branch != 0) {
				trkpt5_branch->GetEntry(index);
			} else { 
				printf("branch trkpt5_branch does not exist!\n");
				exit(1);
			}
			trkpt5_isLoaded = true;
		}
		return trkpt5_;
	}
	float &trkpt10()
	{
		if (not trkpt10_isLoaded) {
			if (trkpt10_branch != 0) {
				trkpt10_branch->GetEntry(index);
			} else { 
				printf("branch trkpt10_branch does not exist!\n");
				exit(1);
			}
			trkpt10_isLoaded = true;
		}
		return trkpt10_;
	}
	float &mleptrk5()
	{
		if (not mleptrk5_isLoaded) {
			if (mleptrk5_branch != 0) {
				mleptrk5_branch->GetEntry(index);
			} else { 
				printf("branch mleptrk5_branch does not exist!\n");
				exit(1);
			}
			mleptrk5_isLoaded = true;
		}
		return mleptrk5_;
	}
	float &mleptrk10()
	{
		if (not mleptrk10_isLoaded) {
			if (mleptrk10_branch != 0) {
				mleptrk10_branch->GetEntry(index);
			} else { 
				printf("branch mleptrk10_branch does not exist!\n");
				exit(1);
			}
			mleptrk10_isLoaded = true;
		}
		return mleptrk10_;
	}
	float &trkreliso5()
	{
		if (not trkreliso5_isLoaded) {
			if (trkreliso5_branch != 0) {
				trkreliso5_branch->GetEntry(index);
			} else { 
				printf("branch trkreliso5_branch does not exist!\n");
				exit(1);
			}
			trkreliso5_isLoaded = true;
		}
		return trkreliso5_;
	}
	float &trkreliso10()
	{
		if (not trkreliso10_isLoaded) {
			if (trkreliso10_branch != 0) {
				trkreliso10_branch->GetEntry(index);
			} else { 
				printf("branch trkreliso10_branch does not exist!\n");
				exit(1);
			}
			trkreliso10_isLoaded = true;
		}
		return trkreliso10_;
	}
	float &trkpt5loose()
	{
		if (not trkpt5loose_isLoaded) {
			if (trkpt5loose_branch != 0) {
				trkpt5loose_branch->GetEntry(index);
			} else { 
				printf("branch trkpt5loose_branch does not exist!\n");
				exit(1);
			}
			trkpt5loose_isLoaded = true;
		}
		return trkpt5loose_;
	}
	float &trkpt10loose()
	{
		if (not trkpt10loose_isLoaded) {
			if (trkpt10loose_branch != 0) {
				trkpt10loose_branch->GetEntry(index);
			} else { 
				printf("branch trkpt10loose_branch does not exist!\n");
				exit(1);
			}
			trkpt10loose_isLoaded = true;
		}
		return trkpt10loose_;
	}
	float &trkreliso5loose()
	{
		if (not trkreliso5loose_isLoaded) {
			if (trkreliso5loose_branch != 0) {
				trkreliso5loose_branch->GetEntry(index);
			} else { 
				printf("branch trkreliso5loose_branch does not exist!\n");
				exit(1);
			}
			trkreliso5loose_isLoaded = true;
		}
		return trkreliso5loose_;
	}
	float &trkreliso10loose()
	{
		if (not trkreliso10loose_isLoaded) {
			if (trkreliso10loose_branch != 0) {
				trkreliso10loose_branch->GetEntry(index);
			} else { 
				printf("branch trkreliso10loose_branch does not exist!\n");
				exit(1);
			}
			trkreliso10loose_isLoaded = true;
		}
		return trkreliso10loose_;
	}
	float &trkpt10pt0p1()
	{
		if (not trkpt10pt0p1_isLoaded) {
			if (trkpt10pt0p1_branch != 0) {
				trkpt10pt0p1_branch->GetEntry(index);
			} else { 
				printf("branch trkpt10pt0p1_branch does not exist!\n");
				exit(1);
			}
			trkpt10pt0p1_isLoaded = true;
		}
		return trkpt10pt0p1_;
	}
	float &trkpt10pt0p2()
	{
		if (not trkpt10pt0p2_isLoaded) {
			if (trkpt10pt0p2_branch != 0) {
				trkpt10pt0p2_branch->GetEntry(index);
			} else { 
				printf("branch trkpt10pt0p2_branch does not exist!\n");
				exit(1);
			}
			trkpt10pt0p2_isLoaded = true;
		}
		return trkpt10pt0p2_;
	}
	float &trkpt10pt0p3()
	{
		if (not trkpt10pt0p3_isLoaded) {
			if (trkpt10pt0p3_branch != 0) {
				trkpt10pt0p3_branch->GetEntry(index);
			} else { 
				printf("branch trkpt10pt0p3_branch does not exist!\n");
				exit(1);
			}
			trkpt10pt0p3_isLoaded = true;
		}
		return trkpt10pt0p3_;
	}
	float &trkpt10pt0p4()
	{
		if (not trkpt10pt0p4_isLoaded) {
			if (trkpt10pt0p4_branch != 0) {
				trkpt10pt0p4_branch->GetEntry(index);
			} else { 
				printf("branch trkpt10pt0p4_branch does not exist!\n");
				exit(1);
			}
			trkpt10pt0p4_isLoaded = true;
		}
		return trkpt10pt0p4_;
	}
	float &trkpt10pt0p5()
	{
		if (not trkpt10pt0p5_isLoaded) {
			if (trkpt10pt0p5_branch != 0) {
				trkpt10pt0p5_branch->GetEntry(index);
			} else { 
				printf("branch trkpt10pt0p5_branch does not exist!\n");
				exit(1);
			}
			trkpt10pt0p5_isLoaded = true;
		}
		return trkpt10pt0p5_;
	}
	float &trkpt10pt0p6()
	{
		if (not trkpt10pt0p6_isLoaded) {
			if (trkpt10pt0p6_branch != 0) {
				trkpt10pt0p6_branch->GetEntry(index);
			} else { 
				printf("branch trkpt10pt0p6_branch does not exist!\n");
				exit(1);
			}
			trkpt10pt0p6_isLoaded = true;
		}
		return trkpt10pt0p6_;
	}
	float &trkpt10pt0p7()
	{
		if (not trkpt10pt0p7_isLoaded) {
			if (trkpt10pt0p7_branch != 0) {
				trkpt10pt0p7_branch->GetEntry(index);
			} else { 
				printf("branch trkpt10pt0p7_branch does not exist!\n");
				exit(1);
			}
			trkpt10pt0p7_isLoaded = true;
		}
		return trkpt10pt0p7_;
	}
	float &trkpt10pt0p8()
	{
		if (not trkpt10pt0p8_isLoaded) {
			if (trkpt10pt0p8_branch != 0) {
				trkpt10pt0p8_branch->GetEntry(index);
			} else { 
				printf("branch trkpt10pt0p8_branch does not exist!\n");
				exit(1);
			}
			trkpt10pt0p8_isLoaded = true;
		}
		return trkpt10pt0p8_;
	}
	float &trkpt10pt0p9()
	{
		if (not trkpt10pt0p9_isLoaded) {
			if (trkpt10pt0p9_branch != 0) {
				trkpt10pt0p9_branch->GetEntry(index);
			} else { 
				printf("branch trkpt10pt0p9_branch does not exist!\n");
				exit(1);
			}
			trkpt10pt0p9_isLoaded = true;
		}
		return trkpt10pt0p9_;
	}
	float &trkpt10pt1p0()
	{
		if (not trkpt10pt1p0_isLoaded) {
			if (trkpt10pt1p0_branch != 0) {
				trkpt10pt1p0_branch->GetEntry(index);
			} else { 
				printf("branch trkpt10pt1p0_branch does not exist!\n");
				exit(1);
			}
			trkpt10pt1p0_isLoaded = true;
		}
		return trkpt10pt1p0_;
	}
	float &trkreliso10pt0p1()
	{
		if (not trkreliso10pt0p1_isLoaded) {
			if (trkreliso10pt0p1_branch != 0) {
				trkreliso10pt0p1_branch->GetEntry(index);
			} else { 
				printf("branch trkreliso10pt0p1_branch does not exist!\n");
				exit(1);
			}
			trkreliso10pt0p1_isLoaded = true;
		}
		return trkreliso10pt0p1_;
	}
	float &trkreliso10pt0p2()
	{
		if (not trkreliso10pt0p2_isLoaded) {
			if (trkreliso10pt0p2_branch != 0) {
				trkreliso10pt0p2_branch->GetEntry(index);
			} else { 
				printf("branch trkreliso10pt0p2_branch does not exist!\n");
				exit(1);
			}
			trkreliso10pt0p2_isLoaded = true;
		}
		return trkreliso10pt0p2_;
	}
	float &trkreliso10pt0p3()
	{
		if (not trkreliso10pt0p3_isLoaded) {
			if (trkreliso10pt0p3_branch != 0) {
				trkreliso10pt0p3_branch->GetEntry(index);
			} else { 
				printf("branch trkreliso10pt0p3_branch does not exist!\n");
				exit(1);
			}
			trkreliso10pt0p3_isLoaded = true;
		}
		return trkreliso10pt0p3_;
	}
	float &trkreliso10pt0p4()
	{
		if (not trkreliso10pt0p4_isLoaded) {
			if (trkreliso10pt0p4_branch != 0) {
				trkreliso10pt0p4_branch->GetEntry(index);
			} else { 
				printf("branch trkreliso10pt0p4_branch does not exist!\n");
				exit(1);
			}
			trkreliso10pt0p4_isLoaded = true;
		}
		return trkreliso10pt0p4_;
	}
	float &trkreliso10pt0p5()
	{
		if (not trkreliso10pt0p5_isLoaded) {
			if (trkreliso10pt0p5_branch != 0) {
				trkreliso10pt0p5_branch->GetEntry(index);
			} else { 
				printf("branch trkreliso10pt0p5_branch does not exist!\n");
				exit(1);
			}
			trkreliso10pt0p5_isLoaded = true;
		}
		return trkreliso10pt0p5_;
	}
	float &trkreliso10pt0p6()
	{
		if (not trkreliso10pt0p6_isLoaded) {
			if (trkreliso10pt0p6_branch != 0) {
				trkreliso10pt0p6_branch->GetEntry(index);
			} else { 
				printf("branch trkreliso10pt0p6_branch does not exist!\n");
				exit(1);
			}
			trkreliso10pt0p6_isLoaded = true;
		}
		return trkreliso10pt0p6_;
	}
	float &trkreliso10pt0p7()
	{
		if (not trkreliso10pt0p7_isLoaded) {
			if (trkreliso10pt0p7_branch != 0) {
				trkreliso10pt0p7_branch->GetEntry(index);
			} else { 
				printf("branch trkreliso10pt0p7_branch does not exist!\n");
				exit(1);
			}
			trkreliso10pt0p7_isLoaded = true;
		}
		return trkreliso10pt0p7_;
	}
	float &trkreliso10pt0p8()
	{
		if (not trkreliso10pt0p8_isLoaded) {
			if (trkreliso10pt0p8_branch != 0) {
				trkreliso10pt0p8_branch->GetEntry(index);
			} else { 
				printf("branch trkreliso10pt0p8_branch does not exist!\n");
				exit(1);
			}
			trkreliso10pt0p8_isLoaded = true;
		}
		return trkreliso10pt0p8_;
	}
	float &trkreliso10pt0p9()
	{
		if (not trkreliso10pt0p9_isLoaded) {
			if (trkreliso10pt0p9_branch != 0) {
				trkreliso10pt0p9_branch->GetEntry(index);
			} else { 
				printf("branch trkreliso10pt0p9_branch does not exist!\n");
				exit(1);
			}
			trkreliso10pt0p9_isLoaded = true;
		}
		return trkreliso10pt0p9_;
	}
	float &trkreliso10pt1p0()
	{
		if (not trkreliso10pt1p0_isLoaded) {
			if (trkreliso10pt1p0_branch != 0) {
				trkreliso10pt1p0_branch->GetEntry(index);
			} else { 
				printf("branch trkreliso10pt1p0_branch does not exist!\n");
				exit(1);
			}
			trkreliso10pt1p0_isLoaded = true;
		}
		return trkreliso10pt1p0_;
	}
	float &pfcandpt10pt0p1()
	{
		if (not pfcandpt10pt0p1_isLoaded) {
			if (pfcandpt10pt0p1_branch != 0) {
				pfcandpt10pt0p1_branch->GetEntry(index);
			} else { 
				printf("branch pfcandpt10pt0p1_branch does not exist!\n");
				exit(1);
			}
			pfcandpt10pt0p1_isLoaded = true;
		}
		return pfcandpt10pt0p1_;
	}
	float &pfcandpt10pt0p2()
	{
		if (not pfcandpt10pt0p2_isLoaded) {
			if (pfcandpt10pt0p2_branch != 0) {
				pfcandpt10pt0p2_branch->GetEntry(index);
			} else { 
				printf("branch pfcandpt10pt0p2_branch does not exist!\n");
				exit(1);
			}
			pfcandpt10pt0p2_isLoaded = true;
		}
		return pfcandpt10pt0p2_;
	}
	float &pfcandpt10pt0p3()
	{
		if (not pfcandpt10pt0p3_isLoaded) {
			if (pfcandpt10pt0p3_branch != 0) {
				pfcandpt10pt0p3_branch->GetEntry(index);
			} else { 
				printf("branch pfcandpt10pt0p3_branch does not exist!\n");
				exit(1);
			}
			pfcandpt10pt0p3_isLoaded = true;
		}
		return pfcandpt10pt0p3_;
	}
	float &pfcandpt10pt0p4()
	{
		if (not pfcandpt10pt0p4_isLoaded) {
			if (pfcandpt10pt0p4_branch != 0) {
				pfcandpt10pt0p4_branch->GetEntry(index);
			} else { 
				printf("branch pfcandpt10pt0p4_branch does not exist!\n");
				exit(1);
			}
			pfcandpt10pt0p4_isLoaded = true;
		}
		return pfcandpt10pt0p4_;
	}
	float &pfcandpt10pt0p5()
	{
		if (not pfcandpt10pt0p5_isLoaded) {
			if (pfcandpt10pt0p5_branch != 0) {
				pfcandpt10pt0p5_branch->GetEntry(index);
			} else { 
				printf("branch pfcandpt10pt0p5_branch does not exist!\n");
				exit(1);
			}
			pfcandpt10pt0p5_isLoaded = true;
		}
		return pfcandpt10pt0p5_;
	}
	float &pfcandpt10pt0p6()
	{
		if (not pfcandpt10pt0p6_isLoaded) {
			if (pfcandpt10pt0p6_branch != 0) {
				pfcandpt10pt0p6_branch->GetEntry(index);
			} else { 
				printf("branch pfcandpt10pt0p6_branch does not exist!\n");
				exit(1);
			}
			pfcandpt10pt0p6_isLoaded = true;
		}
		return pfcandpt10pt0p6_;
	}
	float &pfcandpt10pt0p7()
	{
		if (not pfcandpt10pt0p7_isLoaded) {
			if (pfcandpt10pt0p7_branch != 0) {
				pfcandpt10pt0p7_branch->GetEntry(index);
			} else { 
				printf("branch pfcandpt10pt0p7_branch does not exist!\n");
				exit(1);
			}
			pfcandpt10pt0p7_isLoaded = true;
		}
		return pfcandpt10pt0p7_;
	}
	float &pfcandpt10pt0p8()
	{
		if (not pfcandpt10pt0p8_isLoaded) {
			if (pfcandpt10pt0p8_branch != 0) {
				pfcandpt10pt0p8_branch->GetEntry(index);
			} else { 
				printf("branch pfcandpt10pt0p8_branch does not exist!\n");
				exit(1);
			}
			pfcandpt10pt0p8_isLoaded = true;
		}
		return pfcandpt10pt0p8_;
	}
	float &pfcandpt10pt0p9()
	{
		if (not pfcandpt10pt0p9_isLoaded) {
			if (pfcandpt10pt0p9_branch != 0) {
				pfcandpt10pt0p9_branch->GetEntry(index);
			} else { 
				printf("branch pfcandpt10pt0p9_branch does not exist!\n");
				exit(1);
			}
			pfcandpt10pt0p9_isLoaded = true;
		}
		return pfcandpt10pt0p9_;
	}
	float &pfcandpt10pt1p0()
	{
		if (not pfcandpt10pt1p0_isLoaded) {
			if (pfcandpt10pt1p0_branch != 0) {
				pfcandpt10pt1p0_branch->GetEntry(index);
			} else { 
				printf("branch pfcandpt10pt1p0_branch does not exist!\n");
				exit(1);
			}
			pfcandpt10pt1p0_isLoaded = true;
		}
		return pfcandpt10pt1p0_;
	}
	float &pfcandiso10pt0p1()
	{
		if (not pfcandiso10pt0p1_isLoaded) {
			if (pfcandiso10pt0p1_branch != 0) {
				pfcandiso10pt0p1_branch->GetEntry(index);
			} else { 
				printf("branch pfcandiso10pt0p1_branch does not exist!\n");
				exit(1);
			}
			pfcandiso10pt0p1_isLoaded = true;
		}
		return pfcandiso10pt0p1_;
	}
	float &pfcandiso10pt0p2()
	{
		if (not pfcandiso10pt0p2_isLoaded) {
			if (pfcandiso10pt0p2_branch != 0) {
				pfcandiso10pt0p2_branch->GetEntry(index);
			} else { 
				printf("branch pfcandiso10pt0p2_branch does not exist!\n");
				exit(1);
			}
			pfcandiso10pt0p2_isLoaded = true;
		}
		return pfcandiso10pt0p2_;
	}
	float &pfcandiso10pt0p3()
	{
		if (not pfcandiso10pt0p3_isLoaded) {
			if (pfcandiso10pt0p3_branch != 0) {
				pfcandiso10pt0p3_branch->GetEntry(index);
			} else { 
				printf("branch pfcandiso10pt0p3_branch does not exist!\n");
				exit(1);
			}
			pfcandiso10pt0p3_isLoaded = true;
		}
		return pfcandiso10pt0p3_;
	}
	float &pfcandiso10pt0p4()
	{
		if (not pfcandiso10pt0p4_isLoaded) {
			if (pfcandiso10pt0p4_branch != 0) {
				pfcandiso10pt0p4_branch->GetEntry(index);
			} else { 
				printf("branch pfcandiso10pt0p4_branch does not exist!\n");
				exit(1);
			}
			pfcandiso10pt0p4_isLoaded = true;
		}
		return pfcandiso10pt0p4_;
	}
	float &pfcandiso10pt0p5()
	{
		if (not pfcandiso10pt0p5_isLoaded) {
			if (pfcandiso10pt0p5_branch != 0) {
				pfcandiso10pt0p5_branch->GetEntry(index);
			} else { 
				printf("branch pfcandiso10pt0p5_branch does not exist!\n");
				exit(1);
			}
			pfcandiso10pt0p5_isLoaded = true;
		}
		return pfcandiso10pt0p5_;
	}
	float &pfcandiso10pt0p6()
	{
		if (not pfcandiso10pt0p6_isLoaded) {
			if (pfcandiso10pt0p6_branch != 0) {
				pfcandiso10pt0p6_branch->GetEntry(index);
			} else { 
				printf("branch pfcandiso10pt0p6_branch does not exist!\n");
				exit(1);
			}
			pfcandiso10pt0p6_isLoaded = true;
		}
		return pfcandiso10pt0p6_;
	}
	float &pfcandiso10pt0p7()
	{
		if (not pfcandiso10pt0p7_isLoaded) {
			if (pfcandiso10pt0p7_branch != 0) {
				pfcandiso10pt0p7_branch->GetEntry(index);
			} else { 
				printf("branch pfcandiso10pt0p7_branch does not exist!\n");
				exit(1);
			}
			pfcandiso10pt0p7_isLoaded = true;
		}
		return pfcandiso10pt0p7_;
	}
	float &pfcandiso10pt0p8()
	{
		if (not pfcandiso10pt0p8_isLoaded) {
			if (pfcandiso10pt0p8_branch != 0) {
				pfcandiso10pt0p8_branch->GetEntry(index);
			} else { 
				printf("branch pfcandiso10pt0p8_branch does not exist!\n");
				exit(1);
			}
			pfcandiso10pt0p8_isLoaded = true;
		}
		return pfcandiso10pt0p8_;
	}
	float &pfcandiso10pt0p9()
	{
		if (not pfcandiso10pt0p9_isLoaded) {
			if (pfcandiso10pt0p9_branch != 0) {
				pfcandiso10pt0p9_branch->GetEntry(index);
			} else { 
				printf("branch pfcandiso10pt0p9_branch does not exist!\n");
				exit(1);
			}
			pfcandiso10pt0p9_isLoaded = true;
		}
		return pfcandiso10pt0p9_;
	}
	float &pfcandiso10pt1p0()
	{
		if (not pfcandiso10pt1p0_isLoaded) {
			if (pfcandiso10pt1p0_branch != 0) {
				pfcandiso10pt1p0_branch->GetEntry(index);
			} else { 
				printf("branch pfcandiso10pt1p0_branch does not exist!\n");
				exit(1);
			}
			pfcandiso10pt1p0_isLoaded = true;
		}
		return pfcandiso10pt1p0_;
	}
	float &mbb()
	{
		if (not mbb_isLoaded) {
			if (mbb_branch != 0) {
				mbb_branch->GetEntry(index);
			} else { 
				printf("branch mbb_branch does not exist!\n");
				exit(1);
			}
			mbb_isLoaded = true;
		}
		return mbb_;
	}
	float &lep1pfjetdr()
	{
		if (not lep1pfjetdr_isLoaded) {
			if (lep1pfjetdr_branch != 0) {
				lep1pfjetdr_branch->GetEntry(index);
			} else { 
				printf("branch lep1pfjetdr_branch does not exist!\n");
				exit(1);
			}
			lep1pfjetdr_isLoaded = true;
		}
		return lep1pfjetdr_;
	}
	float &lep2pfjetdr()
	{
		if (not lep2pfjetdr_isLoaded) {
			if (lep2pfjetdr_branch != 0) {
				lep2pfjetdr_branch->GetEntry(index);
			} else { 
				printf("branch lep2pfjetdr_branch does not exist!\n");
				exit(1);
			}
			lep2pfjetdr_isLoaded = true;
		}
		return lep2pfjetdr_;
	}
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &mclep()
	{
		if (not mclep_isLoaded) {
			if (mclep_branch != 0) {
				mclep_branch->GetEntry(index);
			} else { 
				printf("branch mclep_branch does not exist!\n");
				exit(1);
			}
			mclep_isLoaded = true;
		}
		return *mclep_;
	}
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &mcnu()
	{
		if (not mcnu_isLoaded) {
			if (mcnu_branch != 0) {
				mcnu_branch->GetEntry(index);
			} else { 
				printf("branch mcnu_branch does not exist!\n");
				exit(1);
			}
			mcnu_isLoaded = true;
		}
		return *mcnu_;
	}
	float &mcmln()
	{
		if (not mcmln_isLoaded) {
			if (mcmln_branch != 0) {
				mcmln_branch->GetEntry(index);
			} else { 
				printf("branch mcmln_branch does not exist!\n");
				exit(1);
			}
			mcmln_isLoaded = true;
		}
		return mcmln_;
	}
	float &mcmtln()
	{
		if (not mcmtln_isLoaded) {
			if (mcmtln_branch != 0) {
				mcmtln_branch->GetEntry(index);
			} else { 
				printf("branch mcmtln_branch does not exist!\n");
				exit(1);
			}
			mcmtln_isLoaded = true;
		}
		return mcmtln_;
	}
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &mlep()
	{
		if (not mlep_isLoaded) {
			if (mlep_branch != 0) {
				mlep_branch->GetEntry(index);
			} else { 
				printf("branch mlep_branch does not exist!\n");
				exit(1);
			}
			mlep_isLoaded = true;
		}
		return *mlep_;
	}
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &lep1()
	{
		if (not lep1_isLoaded) {
			if (lep1_branch != 0) {
				lep1_branch->GetEntry(index);
			} else { 
				printf("branch lep1_branch does not exist!\n");
				exit(1);
			}
			lep1_isLoaded = true;
		}
		return *lep1_;
	}
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &lep2()
	{
		if (not lep2_isLoaded) {
			if (lep2_branch != 0) {
				lep2_branch->GetEntry(index);
			} else { 
				printf("branch lep2_branch does not exist!\n");
				exit(1);
			}
			lep2_isLoaded = true;
		}
		return *lep2_;
	}
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &trklep1()
	{
		if (not trklep1_isLoaded) {
			if (trklep1_branch != 0) {
				trklep1_branch->GetEntry(index);
			} else { 
				printf("branch trklep1_branch does not exist!\n");
				exit(1);
			}
			trklep1_isLoaded = true;
		}
		return *trklep1_;
	}
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &trklep2()
	{
		if (not trklep2_isLoaded) {
			if (trklep2_branch != 0) {
				trklep2_branch->GetEntry(index);
			} else { 
				printf("branch trklep2_branch does not exist!\n");
				exit(1);
			}
			trklep2_isLoaded = true;
		}
		return *trklep2_;
	}
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &gfitlep1()
	{
		if (not gfitlep1_isLoaded) {
			if (gfitlep1_branch != 0) {
				gfitlep1_branch->GetEntry(index);
			} else { 
				printf("branch gfitlep1_branch does not exist!\n");
				exit(1);
			}
			gfitlep1_isLoaded = true;
		}
		return *gfitlep1_;
	}
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &gfitlep2()
	{
		if (not gfitlep2_isLoaded) {
			if (gfitlep2_branch != 0) {
				gfitlep2_branch->GetEntry(index);
			} else { 
				printf("branch gfitlep2_branch does not exist!\n");
				exit(1);
			}
			gfitlep2_isLoaded = true;
		}
		return *gfitlep2_;
	}
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &lepp()
	{
		if (not lepp_isLoaded) {
			if (lepp_branch != 0) {
				lepp_branch->GetEntry(index);
			} else { 
				printf("branch lepp_branch does not exist!\n");
				exit(1);
			}
			lepp_isLoaded = true;
		}
		return *lepp_;
	}
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &lepm()
	{
		if (not lepm_isLoaded) {
			if (lepm_branch != 0) {
				lepm_branch->GetEntry(index);
			} else { 
				printf("branch lepm_branch does not exist!\n");
				exit(1);
			}
			lepm_isLoaded = true;
		}
		return *lepm_;
	}
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &pflep1()
	{
		if (not pflep1_isLoaded) {
			if (pflep1_branch != 0) {
				pflep1_branch->GetEntry(index);
			} else { 
				printf("branch pflep1_branch does not exist!\n");
				exit(1);
			}
			pflep1_isLoaded = true;
		}
		return *pflep1_;
	}
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &pflep2()
	{
		if (not pflep2_isLoaded) {
			if (pflep2_branch != 0) {
				pflep2_branch->GetEntry(index);
			} else { 
				printf("branch pflep2_branch does not exist!\n");
				exit(1);
			}
			pflep2_isLoaded = true;
		}
		return *pflep2_;
	}
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &leppfjet1()
	{
		if (not leppfjet1_isLoaded) {
			if (leppfjet1_branch != 0) {
				leppfjet1_branch->GetEntry(index);
			} else { 
				printf("branch leppfjet1_branch does not exist!\n");
				exit(1);
			}
			leppfjet1_isLoaded = true;
		}
		return *leppfjet1_;
	}
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &leppfjet2()
	{
		if (not leppfjet2_isLoaded) {
			if (leppfjet2_branch != 0) {
				leppfjet2_branch->GetEntry(index);
			} else { 
				printf("branch leppfjet2_branch does not exist!\n");
				exit(1);
			}
			leppfjet2_isLoaded = true;
		}
		return *leppfjet2_;
	}
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &mclep1()
	{
		if (not mclep1_isLoaded) {
			if (mclep1_branch != 0) {
				mclep1_branch->GetEntry(index);
			} else { 
				printf("branch mclep1_branch does not exist!\n");
				exit(1);
			}
			mclep1_isLoaded = true;
		}
		return *mclep1_;
	}
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &mclep2()
	{
		if (not mclep2_isLoaded) {
			if (mclep2_branch != 0) {
				mclep2_branch->GetEntry(index);
			} else { 
				printf("branch mclep2_branch does not exist!\n");
				exit(1);
			}
			mclep2_isLoaded = true;
		}
		return *mclep2_;
	}
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &mctaud1()
	{
		if (not mctaud1_isLoaded) {
			if (mctaud1_branch != 0) {
				mctaud1_branch->GetEntry(index);
			} else { 
				printf("branch mctaud1_branch does not exist!\n");
				exit(1);
			}
			mctaud1_isLoaded = true;
		}
		return *mctaud1_;
	}
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &mctaud2()
	{
		if (not mctaud2_isLoaded) {
			if (mctaud2_branch != 0) {
				mctaud2_branch->GetEntry(index);
			} else { 
				printf("branch mctaud2_branch does not exist!\n");
				exit(1);
			}
			mctaud2_isLoaded = true;
		}
		return *mctaud2_;
	}
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &mctaudvis1()
	{
		if (not mctaudvis1_isLoaded) {
			if (mctaudvis1_branch != 0) {
				mctaudvis1_branch->GetEntry(index);
			} else { 
				printf("branch mctaudvis1_branch does not exist!\n");
				exit(1);
			}
			mctaudvis1_isLoaded = true;
		}
		return *mctaudvis1_;
	}
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &mctaudvis2()
	{
		if (not mctaudvis2_isLoaded) {
			if (mctaudvis2_branch != 0) {
				mctaudvis2_branch->GetEntry(index);
			} else { 
				printf("branch mctaudvis2_branch does not exist!\n");
				exit(1);
			}
			mctaudvis2_isLoaded = true;
		}
		return *mctaudvis2_;
	}
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &pflep()
	{
		if (not pflep_isLoaded) {
			if (pflep_branch != 0) {
				pflep_branch->GetEntry(index);
			} else { 
				printf("branch pflep_branch does not exist!\n");
				exit(1);
			}
			pflep_isLoaded = true;
		}
		return *pflep_;
	}
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &pftaud()
	{
		if (not pftaud_isLoaded) {
			if (pftaud_branch != 0) {
				pftaud_branch->GetEntry(index);
			} else { 
				printf("branch pftaud_branch does not exist!\n");
				exit(1);
			}
			pftaud_isLoaded = true;
		}
		return *pftaud_;
	}
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &pfcand5()
	{
		if (not pfcand5_isLoaded) {
			if (pfcand5_branch != 0) {
				pfcand5_branch->GetEntry(index);
			} else { 
				printf("branch pfcand5_branch does not exist!\n");
				exit(1);
			}
			pfcand5_isLoaded = true;
		}
		return *pfcand5_;
	}
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &pfcand10()
	{
		if (not pfcand10_isLoaded) {
			if (pfcand10_branch != 0) {
				pfcand10_branch->GetEntry(index);
			} else { 
				printf("branch pfcand10_branch does not exist!\n");
				exit(1);
			}
			pfcand10_isLoaded = true;
		}
		return *pfcand10_;
	}
	int &pfTau15_leadPtcandID()
	{
		if (not pfTau15_leadPtcandID_isLoaded) {
			if (pfTau15_leadPtcandID_branch != 0) {
				pfTau15_leadPtcandID_branch->GetEntry(index);
			} else { 
				printf("branch pfTau15_leadPtcandID_branch does not exist!\n");
				exit(1);
			}
			pfTau15_leadPtcandID_isLoaded = true;
		}
		return pfTau15_leadPtcandID_;
	}
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &pfTau15()
	{
		if (not pfTau15_isLoaded) {
			if (pfTau15_branch != 0) {
				pfTau15_branch->GetEntry(index);
			} else { 
				printf("branch pfTau15_branch does not exist!\n");
				exit(1);
			}
			pfTau15_isLoaded = true;
		}
		return *pfTau15_;
	}
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &pfTau15_leadPtcand()
	{
		if (not pfTau15_leadPtcand_isLoaded) {
			if (pfTau15_leadPtcand_branch != 0) {
				pfTau15_leadPtcand_branch->GetEntry(index);
			} else { 
				printf("branch pfTau15_leadPtcand_branch does not exist!\n");
				exit(1);
			}
			pfTau15_leadPtcand_isLoaded = true;
		}
		return *pfTau15_leadPtcand_;
	}
	int &pfTau_leadPtcandID()
	{
		if (not pfTau_leadPtcandID_isLoaded) {
			if (pfTau_leadPtcandID_branch != 0) {
				pfTau_leadPtcandID_branch->GetEntry(index);
			} else { 
				printf("branch pfTau_leadPtcandID_branch does not exist!\n");
				exit(1);
			}
			pfTau_leadPtcandID_isLoaded = true;
		}
		return pfTau_leadPtcandID_;
	}
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &pfTau()
	{
		if (not pfTau_isLoaded) {
			if (pfTau_branch != 0) {
				pfTau_branch->GetEntry(index);
			} else { 
				printf("branch pfTau_branch does not exist!\n");
				exit(1);
			}
			pfTau_isLoaded = true;
		}
		return *pfTau_;
	}
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &pfTau_leadPtcand()
	{
		if (not pfTau_leadPtcand_isLoaded) {
			if (pfTau_leadPtcand_branch != 0) {
				pfTau_leadPtcand_branch->GetEntry(index);
			} else { 
				printf("branch pfTau_leadPtcand_branch does not exist!\n");
				exit(1);
			}
			pfTau_leadPtcand_isLoaded = true;
		}
		return *pfTau_leadPtcand_;
	}
	int &pfTauLoose_leadPtcandID()
	{
		if (not pfTauLoose_leadPtcandID_isLoaded) {
			if (pfTauLoose_leadPtcandID_branch != 0) {
				pfTauLoose_leadPtcandID_branch->GetEntry(index);
			} else { 
				printf("branch pfTauLoose_leadPtcandID_branch does not exist!\n");
				exit(1);
			}
			pfTauLoose_leadPtcandID_isLoaded = true;
		}
		return pfTauLoose_leadPtcandID_;
	}
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &pfTauLoose()
	{
		if (not pfTauLoose_isLoaded) {
			if (pfTauLoose_branch != 0) {
				pfTauLoose_branch->GetEntry(index);
			} else { 
				printf("branch pfTauLoose_branch does not exist!\n");
				exit(1);
			}
			pfTauLoose_isLoaded = true;
		}
		return *pfTauLoose_;
	}
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &pfTauLoose_leadPtcand()
	{
		if (not pfTauLoose_leadPtcand_isLoaded) {
			if (pfTauLoose_leadPtcand_branch != 0) {
				pfTauLoose_leadPtcand_branch->GetEntry(index);
			} else { 
				printf("branch pfTauLoose_leadPtcand_branch does not exist!\n");
				exit(1);
			}
			pfTauLoose_leadPtcand_isLoaded = true;
		}
		return *pfTauLoose_leadPtcand_;
	}
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &pfcandOS10()
	{
		if (not pfcandOS10_isLoaded) {
			if (pfcandOS10_branch != 0) {
				pfcandOS10_branch->GetEntry(index);
			} else { 
				printf("branch pfcandOS10_branch does not exist!\n");
				exit(1);
			}
			pfcandOS10_isLoaded = true;
		}
		return *pfcandOS10_;
	}
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &pfcandOS10looseZ()
	{
		if (not pfcandOS10looseZ_isLoaded) {
			if (pfcandOS10looseZ_branch != 0) {
				pfcandOS10looseZ_branch->GetEntry(index);
			} else { 
				printf("branch pfcandOS10looseZ_branch does not exist!\n");
				exit(1);
			}
			pfcandOS10looseZ_isLoaded = true;
		}
		return *pfcandOS10looseZ_;
	}
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &pfcand5looseZ()
	{
		if (not pfcand5looseZ_isLoaded) {
			if (pfcand5looseZ_branch != 0) {
				pfcand5looseZ_branch->GetEntry(index);
			} else { 
				printf("branch pfcand5looseZ_branch does not exist!\n");
				exit(1);
			}
			pfcand5looseZ_isLoaded = true;
		}
		return *pfcand5looseZ_;
	}
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &pfcanddir10()
	{
		if (not pfcanddir10_isLoaded) {
			if (pfcanddir10_branch != 0) {
				pfcanddir10_branch->GetEntry(index);
			} else { 
				printf("branch pfcanddir10_branch does not exist!\n");
				exit(1);
			}
			pfcanddir10_isLoaded = true;
		}
		return *pfcanddir10_;
	}
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &pfcandveto10()
	{
		if (not pfcandveto10_isLoaded) {
			if (pfcandveto10_branch != 0) {
				pfcandveto10_branch->GetEntry(index);
			} else { 
				printf("branch pfcandveto10_branch does not exist!\n");
				exit(1);
			}
			pfcandveto10_isLoaded = true;
		}
		return *pfcandveto10_;
	}
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &pfcandvetoL10()
	{
		if (not pfcandvetoL10_isLoaded) {
			if (pfcandvetoL10_branch != 0) {
				pfcandvetoL10_branch->GetEntry(index);
			} else { 
				printf("branch pfcandvetoL10_branch does not exist!\n");
				exit(1);
			}
			pfcandvetoL10_isLoaded = true;
		}
		return *pfcandvetoL10_;
	}
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &jet()
	{
		if (not jet_isLoaded) {
			if (jet_branch != 0) {
				jet_branch->GetEntry(index);
			} else { 
				printf("branch jet_branch does not exist!\n");
				exit(1);
			}
			jet_isLoaded = true;
		}
		return *jet_;
	}
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &nonisoel()
	{
		if (not nonisoel_isLoaded) {
			if (nonisoel_branch != 0) {
				nonisoel_branch->GetEntry(index);
			} else { 
				printf("branch nonisoel_branch does not exist!\n");
				exit(1);
			}
			nonisoel_isLoaded = true;
		}
		return *nonisoel_;
	}
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &nonisomu()
	{
		if (not nonisomu_isLoaded) {
			if (nonisomu_branch != 0) {
				nonisomu_branch->GetEntry(index);
			} else { 
				printf("branch nonisomu_branch does not exist!\n");
				exit(1);
			}
			nonisomu_isLoaded = true;
		}
		return *nonisomu_;
	}
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &t()
	{
		if (not t_isLoaded) {
			if (t_branch != 0) {
				t_branch->GetEntry(index);
			} else { 
				printf("branch t_branch does not exist!\n");
				exit(1);
			}
			t_isLoaded = true;
		}
		return *t_;
	}
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &tbar()
	{
		if (not tbar_isLoaded) {
			if (tbar_branch != 0) {
				tbar_branch->GetEntry(index);
			} else { 
				printf("branch tbar_branch does not exist!\n");
				exit(1);
			}
			tbar_isLoaded = true;
		}
		return *tbar_;
	}
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &ttbar()
	{
		if (not ttbar_isLoaded) {
			if (ttbar_branch != 0) {
				ttbar_branch->GetEntry(index);
			} else { 
				printf("branch ttbar_branch does not exist!\n");
				exit(1);
			}
			ttbar_isLoaded = true;
		}
		return *ttbar_;
	}
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &lep_t()
	{
		if (not lep_t_isLoaded) {
			if (lep_t_branch != 0) {
				lep_t_branch->GetEntry(index);
			} else { 
				printf("branch lep_t_branch does not exist!\n");
				exit(1);
			}
			lep_t_isLoaded = true;
		}
		return *lep_t_;
	}
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &lep_tbar()
	{
		if (not lep_tbar_isLoaded) {
			if (lep_tbar_branch != 0) {
				lep_tbar_branch->GetEntry(index);
			} else { 
				printf("branch lep_tbar_branch does not exist!\n");
				exit(1);
			}
			lep_tbar_isLoaded = true;
		}
		return *lep_tbar_;
	}
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &stop_t()
	{
		if (not stop_t_isLoaded) {
			if (stop_t_branch != 0) {
				stop_t_branch->GetEntry(index);
			} else { 
				printf("branch stop_t_branch does not exist!\n");
				exit(1);
			}
			stop_t_isLoaded = true;
		}
		return *stop_t_;
	}
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &stop_tbar()
	{
		if (not stop_tbar_isLoaded) {
			if (stop_tbar_branch != 0) {
				stop_tbar_branch->GetEntry(index);
			} else { 
				printf("branch stop_tbar_branch does not exist!\n");
				exit(1);
			}
			stop_tbar_isLoaded = true;
		}
		return *stop_tbar_;
	}
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &neutralino_t()
	{
		if (not neutralino_t_isLoaded) {
			if (neutralino_t_branch != 0) {
				neutralino_t_branch->GetEntry(index);
			} else { 
				printf("branch neutralino_t_branch does not exist!\n");
				exit(1);
			}
			neutralino_t_isLoaded = true;
		}
		return *neutralino_t_;
	}
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &neutralino_tbar()
	{
		if (not neutralino_tbar_isLoaded) {
			if (neutralino_tbar_branch != 0) {
				neutralino_tbar_branch->GetEntry(index);
			} else { 
				printf("branch neutralino_tbar_branch does not exist!\n");
				exit(1);
			}
			neutralino_tbar_isLoaded = true;
		}
		return *neutralino_tbar_;
	}
	int &lep_t_id()
	{
		if (not lep_t_id_isLoaded) {
			if (lep_t_id_branch != 0) {
				lep_t_id_branch->GetEntry(index);
			} else { 
				printf("branch lep_t_id_branch does not exist!\n");
				exit(1);
			}
			lep_t_id_isLoaded = true;
		}
		return lep_t_id_;
	}
	int &lep_tbar_id()
	{
		if (not lep_tbar_id_isLoaded) {
			if (lep_tbar_id_branch != 0) {
				lep_tbar_id_branch->GetEntry(index);
			} else { 
				printf("branch lep_tbar_id_branch does not exist!\n");
				exit(1);
			}
			lep_tbar_id_isLoaded = true;
		}
		return lep_tbar_id_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &pfjets()
	{
		if (not pfjets_isLoaded) {
			if (pfjets_branch != 0) {
				pfjets_branch->GetEntry(index);
			} else { 
				printf("branch pfjets_branch does not exist!\n");
				exit(1);
			}
			pfjets_isLoaded = true;
		}
		return *pfjets_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &pfjets_genJet_()
	{
		if (not pfjets_genJet__isLoaded) {
			if (pfjets_genJet__branch != 0) {
				pfjets_genJet__branch->GetEntry(index);
			} else { 
				printf("branch pfjets_genJet__branch does not exist!\n");
				exit(1);
			}
			pfjets_genJet__isLoaded = true;
		}
		return *pfjets_genJet__;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &pfjets_failjetid()
	{
		if (not pfjets_failjetid_isLoaded) {
			if (pfjets_failjetid_branch != 0) {
				pfjets_failjetid_branch->GetEntry(index);
			} else { 
				printf("branch pfjets_failjetid_branch does not exist!\n");
				exit(1);
			}
			pfjets_failjetid_isLoaded = true;
		}
		return *pfjets_failjetid_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &pfjets_faillepolap()
	{
		if (not pfjets_faillepolap_isLoaded) {
			if (pfjets_faillepolap_branch != 0) {
				pfjets_faillepolap_branch->GetEntry(index);
			} else { 
				printf("branch pfjets_faillepolap_branch does not exist!\n");
				exit(1);
			}
			pfjets_faillepolap_isLoaded = true;
		}
		return *pfjets_faillepolap_;
	}
	vector<float> &pfjets_csv()
	{
		if (not pfjets_csv_isLoaded) {
			if (pfjets_csv_branch != 0) {
				pfjets_csv_branch->GetEntry(index);
			} else { 
				printf("branch pfjets_csv_branch does not exist!\n");
				exit(1);
			}
			pfjets_csv_isLoaded = true;
		}
		return *pfjets_csv_;
	}
	vector<float> &pfjets_chEfrac()
	{
		if (not pfjets_chEfrac_isLoaded) {
			if (pfjets_chEfrac_branch != 0) {
				pfjets_chEfrac_branch->GetEntry(index);
			} else { 
				printf("branch pfjets_chEfrac_branch does not exist!\n");
				exit(1);
			}
			pfjets_chEfrac_isLoaded = true;
		}
		return *pfjets_chEfrac_;
	}
	vector<float> &pfjets_chm()
	{
		if (not pfjets_chm_isLoaded) {
			if (pfjets_chm_branch != 0) {
				pfjets_chm_branch->GetEntry(index);
			} else { 
				printf("branch pfjets_chm_branch does not exist!\n");
				exit(1);
			}
			pfjets_chm_isLoaded = true;
		}
		return *pfjets_chm_;
	}
	vector<float> &pfjets_neu()
	{
		if (not pfjets_neu_isLoaded) {
			if (pfjets_neu_branch != 0) {
				pfjets_neu_branch->GetEntry(index);
			} else { 
				printf("branch pfjets_neu_branch does not exist!\n");
				exit(1);
			}
			pfjets_neu_isLoaded = true;
		}
		return *pfjets_neu_;
	}
	vector<float> &pfjets_l1corr()
	{
		if (not pfjets_l1corr_isLoaded) {
			if (pfjets_l1corr_branch != 0) {
				pfjets_l1corr_branch->GetEntry(index);
			} else { 
				printf("branch pfjets_l1corr_branch does not exist!\n");
				exit(1);
			}
			pfjets_l1corr_isLoaded = true;
		}
		return *pfjets_l1corr_;
	}
	vector<float> &pfjets_corr()
	{
		if (not pfjets_corr_isLoaded) {
			if (pfjets_corr_branch != 0) {
				pfjets_corr_branch->GetEntry(index);
			} else { 
				printf("branch pfjets_corr_branch does not exist!\n");
				exit(1);
			}
			pfjets_corr_isLoaded = true;
		}
		return *pfjets_corr_;
	}
	vector<int> &pfjets_mc3()
	{
		if (not pfjets_mc3_isLoaded) {
			if (pfjets_mc3_branch != 0) {
				pfjets_mc3_branch->GetEntry(index);
			} else { 
				printf("branch pfjets_mc3_branch does not exist!\n");
				exit(1);
			}
			pfjets_mc3_isLoaded = true;
		}
		return *pfjets_mc3_;
	}
	vector<int> &pfjets_mcflavorAlgo()
	{
		if (not pfjets_mcflavorAlgo_isLoaded) {
			if (pfjets_mcflavorAlgo_branch != 0) {
				pfjets_mcflavorAlgo_branch->GetEntry(index);
			} else { 
				printf("branch pfjets_mcflavorAlgo_branch does not exist!\n");
				exit(1);
			}
			pfjets_mcflavorAlgo_isLoaded = true;
		}
		return *pfjets_mcflavorAlgo_;
	}
	vector<int> &pfjets_mcflavorPhys()
	{
		if (not pfjets_mcflavorPhys_isLoaded) {
			if (pfjets_mcflavorPhys_branch != 0) {
				pfjets_mcflavorPhys_branch->GetEntry(index);
			} else { 
				printf("branch pfjets_mcflavorPhys_branch does not exist!\n");
				exit(1);
			}
			pfjets_mcflavorPhys_isLoaded = true;
		}
		return *pfjets_mcflavorPhys_;
	}
	vector<float> &pfjets_uncertainty()
	{
		if (not pfjets_uncertainty_isLoaded) {
			if (pfjets_uncertainty_branch != 0) {
				pfjets_uncertainty_branch->GetEntry(index);
			} else { 
				printf("branch pfjets_uncertainty_branch does not exist!\n");
				exit(1);
			}
			pfjets_uncertainty_isLoaded = true;
		}
		return *pfjets_uncertainty_;
	}
	vector<int> &pfjets_flav()
	{
		if (not pfjets_flav_isLoaded) {
			if (pfjets_flav_branch != 0) {
				pfjets_flav_branch->GetEntry(index);
			} else { 
				printf("branch pfjets_flav_branch does not exist!\n");
				exit(1);
			}
			pfjets_flav_isLoaded = true;
		}
		return *pfjets_flav_;
	}
	vector<float> &pfjets_lrm()
	{
		if (not pfjets_lrm_isLoaded) {
			if (pfjets_lrm_branch != 0) {
				pfjets_lrm_branch->GetEntry(index);
			} else { 
				printf("branch pfjets_lrm_branch does not exist!\n");
				exit(1);
			}
			pfjets_lrm_isLoaded = true;
		}
		return *pfjets_lrm_;
	}
	vector<float> &pfjets_lrm2()
	{
		if (not pfjets_lrm2_isLoaded) {
			if (pfjets_lrm2_branch != 0) {
				pfjets_lrm2_branch->GetEntry(index);
			} else { 
				printf("branch pfjets_lrm2_branch does not exist!\n");
				exit(1);
			}
			pfjets_lrm2_isLoaded = true;
		}
		return *pfjets_lrm2_;
	}
	vector<float> &pfjets_qgtag()
	{
		if (not pfjets_qgtag_isLoaded) {
			if (pfjets_qgtag_branch != 0) {
				pfjets_qgtag_branch->GetEntry(index);
			} else { 
				printf("branch pfjets_qgtag_branch does not exist!\n");
				exit(1);
			}
			pfjets_qgtag_isLoaded = true;
		}
		return *pfjets_qgtag_;
	}
	vector<float> &pfjets_genJetDr()
	{
		if (not pfjets_genJetDr_isLoaded) {
			if (pfjets_genJetDr_branch != 0) {
				pfjets_genJetDr_branch->GetEntry(index);
			} else { 
				printf("branch pfjets_genJetDr_branch does not exist!\n");
				exit(1);
			}
			pfjets_genJetDr_isLoaded = true;
		}
		return *pfjets_genJetDr_;
	}
	vector<float> &pfjets_sigma()
	{
		if (not pfjets_sigma_isLoaded) {
			if (pfjets_sigma_branch != 0) {
				pfjets_sigma_branch->GetEntry(index);
			} else { 
				printf("branch pfjets_sigma_branch does not exist!\n");
				exit(1);
			}
			pfjets_sigma_isLoaded = true;
		}
		return *pfjets_sigma_;
	}
	vector<int> &pfjets_lepjet()
	{
		if (not pfjets_lepjet_isLoaded) {
			if (pfjets_lepjet_branch != 0) {
				pfjets_lepjet_branch->GetEntry(index);
			} else { 
				printf("branch pfjets_lepjet_branch does not exist!\n");
				exit(1);
			}
			pfjets_lepjet_isLoaded = true;
		}
		return *pfjets_lepjet_;
	}
	vector<float> &pfjets_tobtecmult()
	{
		if (not pfjets_tobtecmult_isLoaded) {
			if (pfjets_tobtecmult_branch != 0) {
				pfjets_tobtecmult_branch->GetEntry(index);
			} else { 
				printf("branch pfjets_tobtecmult_branch does not exist!\n");
				exit(1);
			}
			pfjets_tobtecmult_isLoaded = true;
		}
		return *pfjets_tobtecmult_;
	}
	vector<float> &pfjets_tobtecfrac()
	{
		if (not pfjets_tobtecfrac_isLoaded) {
			if (pfjets_tobtecfrac_branch != 0) {
				pfjets_tobtecfrac_branch->GetEntry(index);
			} else { 
				printf("branch pfjets_tobtecfrac_branch does not exist!\n");
				exit(1);
			}
			pfjets_tobtecfrac_isLoaded = true;
		}
		return *pfjets_tobtecfrac_;
	}
	vector<float> &pfjets_beta()
	{
		if (not pfjets_beta_isLoaded) {
			if (pfjets_beta_branch != 0) {
				pfjets_beta_branch->GetEntry(index);
			} else { 
				printf("branch pfjets_beta_branch does not exist!\n");
				exit(1);
			}
			pfjets_beta_isLoaded = true;
		}
		return *pfjets_beta_;
	}
	vector<float> &pfjets_beta2()
	{
		if (not pfjets_beta2_isLoaded) {
			if (pfjets_beta2_branch != 0) {
				pfjets_beta2_branch->GetEntry(index);
			} else { 
				printf("branch pfjets_beta2_branch does not exist!\n");
				exit(1);
			}
			pfjets_beta2_isLoaded = true;
		}
		return *pfjets_beta2_;
	}
	vector<float> &pfjets_beta_0p1()
	{
		if (not pfjets_beta_0p1_isLoaded) {
			if (pfjets_beta_0p1_branch != 0) {
				pfjets_beta_0p1_branch->GetEntry(index);
			} else { 
				printf("branch pfjets_beta_0p1_branch does not exist!\n");
				exit(1);
			}
			pfjets_beta_0p1_isLoaded = true;
		}
		return *pfjets_beta_0p1_;
	}
	vector<float> &pfjets_beta_0p2()
	{
		if (not pfjets_beta_0p2_isLoaded) {
			if (pfjets_beta_0p2_branch != 0) {
				pfjets_beta_0p2_branch->GetEntry(index);
			} else { 
				printf("branch pfjets_beta_0p2_branch does not exist!\n");
				exit(1);
			}
			pfjets_beta_0p2_isLoaded = true;
		}
		return *pfjets_beta_0p2_;
	}
	vector<float> &pfjets_beta2_0p1()
	{
		if (not pfjets_beta2_0p1_isLoaded) {
			if (pfjets_beta2_0p1_branch != 0) {
				pfjets_beta2_0p1_branch->GetEntry(index);
			} else { 
				printf("branch pfjets_beta2_0p1_branch does not exist!\n");
				exit(1);
			}
			pfjets_beta2_0p1_isLoaded = true;
		}
		return *pfjets_beta2_0p1_;
	}
	vector<float> &pfjets_beta2_0p5()
	{
		if (not pfjets_beta2_0p5_isLoaded) {
			if (pfjets_beta2_0p5_branch != 0) {
				pfjets_beta2_0p5_branch->GetEntry(index);
			} else { 
				printf("branch pfjets_beta2_0p5_branch does not exist!\n");
				exit(1);
			}
			pfjets_beta2_0p5_isLoaded = true;
		}
		return *pfjets_beta2_0p5_;
	}
	vector<float> &pfjets_mvaPUid()
	{
		if (not pfjets_mvaPUid_isLoaded) {
			if (pfjets_mvaPUid_branch != 0) {
				pfjets_mvaPUid_branch->GetEntry(index);
			} else { 
				printf("branch pfjets_mvaPUid_branch does not exist!\n");
				exit(1);
			}
			pfjets_mvaPUid_isLoaded = true;
		}
		return *pfjets_mvaPUid_;
	}
	vector<float> &pfjets_mva5xPUid()
	{
		if (not pfjets_mva5xPUid_isLoaded) {
			if (pfjets_mva5xPUid_branch != 0) {
				pfjets_mva5xPUid_branch->GetEntry(index);
			} else { 
				printf("branch pfjets_mva5xPUid_branch does not exist!\n");
				exit(1);
			}
			pfjets_mva5xPUid_isLoaded = true;
		}
		return *pfjets_mva5xPUid_;
	}
	vector<float> &pfjets_mvaBeta()
	{
		if (not pfjets_mvaBeta_isLoaded) {
			if (pfjets_mvaBeta_branch != 0) {
				pfjets_mvaBeta_branch->GetEntry(index);
			} else { 
				printf("branch pfjets_mvaBeta_branch does not exist!\n");
				exit(1);
			}
			pfjets_mvaBeta_isLoaded = true;
		}
		return *pfjets_mvaBeta_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &genjets()
	{
		if (not genjets_isLoaded) {
			if (genjets_branch != 0) {
				genjets_branch->GetEntry(index);
			} else { 
				printf("branch genjets_branch does not exist!\n");
				exit(1);
			}
			genjets_isLoaded = true;
		}
		return *genjets_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &genqgs()
	{
		if (not genqgs_isLoaded) {
			if (genqgs_branch != 0) {
				genqgs_branch->GetEntry(index);
			} else { 
				printf("branch genqgs_branch does not exist!\n");
				exit(1);
			}
			genqgs_isLoaded = true;
		}
		return *genqgs_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &genbs()
	{
		if (not genbs_isLoaded) {
			if (genbs_branch != 0) {
				genbs_branch->GetEntry(index);
			} else { 
				printf("branch genbs_branch does not exist!\n");
				exit(1);
			}
			genbs_isLoaded = true;
		}
		return *genbs_;
	}
	vector<int> &genps_pdgId()
	{
		if (not genps_pdgId_isLoaded) {
			if (genps_pdgId_branch != 0) {
				genps_pdgId_branch->GetEntry(index);
			} else { 
				printf("branch genps_pdgId_branch does not exist!\n");
				exit(1);
			}
			genps_pdgId_isLoaded = true;
		}
		return *genps_pdgId_;
	}
	vector<int> &genps_firstMother()
	{
		if (not genps_firstMother_isLoaded) {
			if (genps_firstMother_branch != 0) {
				genps_firstMother_branch->GetEntry(index);
			} else { 
				printf("branch genps_firstMother_branch does not exist!\n");
				exit(1);
			}
			genps_firstMother_isLoaded = true;
		}
		return *genps_firstMother_;
	}
	vector<float> &genps_energy()
	{
		if (not genps_energy_isLoaded) {
			if (genps_energy_branch != 0) {
				genps_energy_branch->GetEntry(index);
			} else { 
				printf("branch genps_energy_branch does not exist!\n");
				exit(1);
			}
			genps_energy_isLoaded = true;
		}
		return *genps_energy_;
	}
	vector<float> &genps_pt()
	{
		if (not genps_pt_isLoaded) {
			if (genps_pt_branch != 0) {
				genps_pt_branch->GetEntry(index);
			} else { 
				printf("branch genps_pt_branch does not exist!\n");
				exit(1);
			}
			genps_pt_isLoaded = true;
		}
		return *genps_pt_;
	}
	vector<float> &genps_eta()
	{
		if (not genps_eta_isLoaded) {
			if (genps_eta_branch != 0) {
				genps_eta_branch->GetEntry(index);
			} else { 
				printf("branch genps_eta_branch does not exist!\n");
				exit(1);
			}
			genps_eta_isLoaded = true;
		}
		return *genps_eta_;
	}
	vector<float> &genps_phi()
	{
		if (not genps_phi_isLoaded) {
			if (genps_phi_branch != 0) {
				genps_phi_branch->GetEntry(index);
			} else { 
				printf("branch genps_phi_branch does not exist!\n");
				exit(1);
			}
			genps_phi_isLoaded = true;
		}
		return *genps_phi_;
	}
	vector<float> &genps_mass()
	{
		if (not genps_mass_isLoaded) {
			if (genps_mass_branch != 0) {
				genps_mass_branch->GetEntry(index);
			} else { 
				printf("branch genps_mass_branch does not exist!\n");
				exit(1);
			}
			genps_mass_isLoaded = true;
		}
		return *genps_mass_;
	}

  static void progress( int nEventsTotal, int nEventsChain ){
    int period = 1000;
    if(nEventsTotal%1000 == 0) {
      // xterm magic from L. Vacavant and A. Cerri
      if (isatty(1)) {
        if( ( nEventsChain - nEventsTotal ) > period ){
          float frac = (float)nEventsTotal/(nEventsChain*0.01);
          printf("\015\033[32m ---> \033[1m\033[31m%4.1f%%"
               "\033[0m\033[32m <---\033[0m\015", frac);
          fflush(stdout);
        }
        else {
          printf("\015\033[32m ---> \033[1m\033[31m%4.1f%%"
                 "\033[0m\033[32m <---\033[0m\015", 100.);
          cout << endl;
        }
      }
    }
  }
  
};

#ifndef __CINT__
extern STOPT stopt;
#endif

namespace Stop {
	int &acc_2010();
	int &acc_highmet();
	int &acc_highht();
	int &eldup();
	int &csc();
	int &hbhe();
	int &hbhenew();
	int &hcallaser();
	int &ecaltp();
	int &trkfail();
	int &eebadsc();
	int &lep1_badecallaser();
	int &lep2_badecallaser();
	int &isdata();
	int &jetid();
	int &jetid30();
	int &json();
	float &htoffset();
	float &htuncor();
	float &ptt();
	float &pttbar();
	float &ptttbar();
	float &mttbar();
	int &npartons();
	int &nwzpartons();
	int &hyptype();
	float &maxpartonpt();
	float &etattbar();
	int &njetsoffset();
	int &njetsuncor();
	float &costhetaweight();
	float &weight();
	float &weightleft();
	float &weightright();
	float &mutrigweight();
	float &mutrigweight2();
	float &sltrigweight();
	float &dltrigweight();
	float &trgeff();
	float &pthat();
	float &qscale();
	float &mgcor();
	int &wflav();
	float &ksusy();
	float &ksusyup();
	float &ksusydn();
	float &xsecsusy();
	float &xsecsusy2();
	float &smeff();
	float &k();
	float &mllgen();
	float &ptwgen();
	float &ptzgen();
	int &nlep();
	int &nosel();
	int &ngoodlep();
	int &ngoodel();
	int &ngoodmu();
	int &mull();
	int &mult();
	int &mullgen();
	int &multgen();
	int &proc();
	int &leptype();
	float &topmass();
	float &dilmass();
	float &dilrecoil();
	float &dilrecoilparl();
	float &dilrecoilperp();
	float &tcmet();
	float &genmet();
	float &gensumet();
	float &genmetphi();
	float &calomet();
	float &calometphi();
	float &trkmet();
	float &trkmetphi();
	float &pfmet();
	float &pfmetveto();
	float &pfmetsig();
	float &pfmetphi();
	float &pfsumet();
	float &mucormet();
	float &mucorjesmet();
	float &tcmet35X();
	float &tcmetevent();
	float &tcmetlooper();
	float &tcmetphi();
	float &tcsumet();
	float &tcmetUp();
	float &tcmetDown();
	float &tcmetTest();
	float &pfmetUp();
	float &pfmetDown();
	float &pfmetTest();
	float &sumjetpt();
	float &dileta();
	float &dilpt();
	float &dildphi();
	int &ngenjets();
	int &njpt();
	int &trgmu1();
	int &trgmu2();
	int &trgel1();
	int &trgel2();
	int &isomu24();
	int &ele27wp80();
	int &mm();
	int &mmtk();
	int &me();
	int &em();
	int &mu();
	int &ee();
	int &npfjets30();
	int &npfjets30lepcorr();
	float &knjets();
	float &rhovor();
	float &htpf30();
	float &t1met10();
	float &t1met20();
	float &t1met30();
	float &t1met10phi();
	float &t1met20phi();
	float &t1met30phi();
	float &t1met10mt();
	float &t1met20mt();
	float &t1met30mt();
	float &lepmetpt();
	float &lept1met10pt();
	float &t1met10s();
	float &t1met10sphi();
	float &t1met10smt();
	float &t1metphicorr();
	float &t1metphicorrup();
	float &t1metphicorrdn();
	float &t1metphicorrphi();
	float &t1metphicorrphiup();
	float &t1metphicorrphidn();
	float &t1metphicorrlep();
	float &t1metphicorrlepphi();
	float &t1metphicorrmt();
	float &t1metphicorrmtup();
	float &t1metphicorrmtdn();
	float &t1metphicorrlepmt();
	float &t1met_off();
	float &t1metphi_off();
	float &t1metmt_off();
	float &t1metphicorr_off();
	float &t1metphicorrphi_off();
	float &t1metphicorrmt_off();
	float &mht15();
	float &mht15phi();
	float &trkmet_mht15();
	float &trkmetphi_mht15();
	float &mettlj15();
	float &mettlj15phi();
	float &mt2bmin();
	float &mt2blmin();
	float &mt2wmin();
	float &chi2min();
	float &chi2minprob();
	int &nbtagsssv();
	int &nbtagstcl();
	int &nbtagstcm();
	int &nbtagscsvl();
	int &nbtagscsvm();
	int &nbtagscsvt();
	int &nbtagsssvcorr();
	int &nbtagstclcorr();
	int &nbtagstcmcorr();
	int &nbtagscsvlcorr();
	int &nbtagscsvmcorr();
	int &nbtagscsvtcott();
	int &njetsUp();
	int &njetsDown();
	float &htUp();
	float &htDown();
	int &ntruepu();
	int &npu();
	int &npuMinusOne();
	int &npuPlusOne();
	int &nvtx();
	int &indexfirstGoodVertex_();
	float &nvtxweight();
	float &n3dvtxweight();
	int &pdfid1();
	int &pdfid2();
	float &pdfx1();
	float &pdfx2();
	float &pdfQ();
	float &vecjetpt();
	int &pass();
	int &passz();
	float &m0();
	float &mg();
	float &ml();
	float &x();
	float &m12();
	float &lep1chi2ndf();
	float &lep2chi2ndf();
	float &lep1dpt();
	float &lep2dpt();
	int &id1();
	int &id2();
	int &leptype1();
	int &leptype2();
	int &w1();
	int &w2();
	float &iso1();
	float &isont1();
	float &isopfold1();
	float &isopf1();
	float &etasc1();
	float &etasc2();
	float &eoverpin();
	float &eoverpout();
	float &dEtaIn();
	float &dPhiIn();
	float &sigmaIEtaIEta();
	float &hOverE();
	float &ooemoop();
	float &d0vtx();
	float &dzvtx();
	float &expinnerlayers();
	float &fbrem();
	float &pfisoch();
	float &pfisoem();
	float &pfisonh();
	float &eSC();
	float &phiSC();
	float &eSCRaw();
	float &eSCPresh();
	float &lep1_scslasercormean();
	float &lep1_scslasercormax();
	float &eoverpin2();
	float &eoverpout2();
	float &dEtaIn2();
	float &dPhiIn2();
	float &sigmaIEtaIEta2();
	float &hOverE2();
	float &ooemoop2();
	float &d0vtx2();
	float &dzvtx2();
	float &expinnerlayers2();
	float &fbrem2();
	float &pfisoch2();
	float &pfisoem2();
	float &pfisonh2();
	float &eSC2();
	float &phiSC2();
	float &eSCRaw2();
	float &eSCPresh2();
	float &lep2_scslasercormean();
	float &lep2_scslasercormax();
	float &scslasercormax();
	float &scslasercormax_pt();
	float &scslasercormax_eta();
	float &iso2();
	float &ecalveto1();
	float &ecalveto2();
	float &hcalveto1();
	float &hcalveto2();
	float &isont2();
	float &isopf2();
	float &ptl1();
	float &ptl2();
	float &etal1();
	float &etal2();
	float &phil1();
	float &phil2();
	float &meff();
	float &mt();
	unsigned int &run();
	unsigned int &lumi();
	unsigned int &event();
	float &y();
	float &ht();
	float &htgen();
	float &htjpt();
	int &nels();
	int &nmus();
	int &ntaus();
	int &nleps();
	int &nbs();
	float &dphijm();
	float &ptjetraw();
	float &ptjet23();
	float &ptjetF23();
	float &ptjetO23();
	int &mcid1();
	float &mcdr1();
	int &mcdecay1();
	int &mcndec1();
	int &mcndec2();
	int &mcndeckls1();
	int &mcndeckls2();
	int &mcndecem1();
	int &mcndecem2();
	int &mcid2();
	float &mcdr2();
	int &mcdecay2();
	float &mctaudpt1();
	float &mctaudpt2();
	int &mctaudid1();
	int &mctaudid2();
	int &mlepid();
	int &mleppassid();
	int &mleppassiso();
	float &mlepiso();
	float &mlepdr();
	float &pflepiso();
	float &pflepdr();
	float &pfleppt();
	float &pflepmindrj();
	float &pftaudiso();
	float &pftauddr();
	float &pftaudpt();
	float &pftaudmindrj();
	int &pfcandid5();
	float &pfcandiso5();
	float &pfcandpt5();
	float &pfcanddz5();
	float &pfcandmindrj5();
	int &pfcandid10();
	float &pfcandiso10();
	float &pfcandpt10();
	float &pfcanddz10();
	float &pfcandmindrj10();
	int &pfcandidOS10();
	float &pfcandisoOS10();
	float &pfcandptOS10();
	float &pfcanddzOS10();
	int &pfcandid5looseZ();
	float &pfcandiso5looseZ();
	float &pfcandpt5looseZ();
	float &pfcanddz5looseZ();
	int &pfcandidOS10looseZ();
	float &pfcandisoOS10looseZ();
	float &pfcandptOS10looseZ();
	float &pfcanddzOS10looseZ();
	int &pfcanddirid10();
	float &pfcanddiriso10();
	float &pfcanddirpt10();
	float &pfcanddirmindrj10();
	int &pfcandvetoid10();
	float &pfcandvetoiso10();
	float &pfcandvetopt10();
	float &pfcandvetomindrj10();
	int &pfcandvetoLid10();
	float &pfcandvetoLiso10();
	float &pfcandvetoLpt10();
	float &pfcandvetoLmindrj10();
	float &emjet10();
	float &mjj();
	float &emjet20();
	float &trkpt5();
	float &trkpt10();
	float &mleptrk5();
	float &mleptrk10();
	float &trkreliso5();
	float &trkreliso10();
	float &trkpt5loose();
	float &trkpt10loose();
	float &trkreliso5loose();
	float &trkreliso10loose();
	float &trkpt10pt0p1();
	float &trkpt10pt0p2();
	float &trkpt10pt0p3();
	float &trkpt10pt0p4();
	float &trkpt10pt0p5();
	float &trkpt10pt0p6();
	float &trkpt10pt0p7();
	float &trkpt10pt0p8();
	float &trkpt10pt0p9();
	float &trkpt10pt1p0();
	float &trkreliso10pt0p1();
	float &trkreliso10pt0p2();
	float &trkreliso10pt0p3();
	float &trkreliso10pt0p4();
	float &trkreliso10pt0p5();
	float &trkreliso10pt0p6();
	float &trkreliso10pt0p7();
	float &trkreliso10pt0p8();
	float &trkreliso10pt0p9();
	float &trkreliso10pt1p0();
	float &pfcandpt10pt0p1();
	float &pfcandpt10pt0p2();
	float &pfcandpt10pt0p3();
	float &pfcandpt10pt0p4();
	float &pfcandpt10pt0p5();
	float &pfcandpt10pt0p6();
	float &pfcandpt10pt0p7();
	float &pfcandpt10pt0p8();
	float &pfcandpt10pt0p9();
	float &pfcandpt10pt1p0();
	float &pfcandiso10pt0p1();
	float &pfcandiso10pt0p2();
	float &pfcandiso10pt0p3();
	float &pfcandiso10pt0p4();
	float &pfcandiso10pt0p5();
	float &pfcandiso10pt0p6();
	float &pfcandiso10pt0p7();
	float &pfcandiso10pt0p8();
	float &pfcandiso10pt0p9();
	float &pfcandiso10pt1p0();
	float &mbb();
	float &lep1pfjetdr();
	float &lep2pfjetdr();
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &mclep();
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &mcnu();
	float &mcmln();
	float &mcmtln();
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &mlep();
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &lep1();
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &lep2();
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &trklep1();
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &trklep2();
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &gfitlep1();
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &gfitlep2();
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &lepp();
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &lepm();
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &pflep1();
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &pflep2();
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &leppfjet1();
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &leppfjet2();
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &mclep1();
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &mclep2();
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &mctaud1();
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &mctaud2();
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &mctaudvis1();
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &mctaudvis2();
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &pflep();
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &pftaud();
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &pfcand5();
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &pfcand10();
	int &pfTau15_leadPtcandID();
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &pfTau15();
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &pfTau15_leadPtcand();
	int &pfTau_leadPtcandID();
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &pfTau();
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &pfTau_leadPtcand();
	int &pfTauLoose_leadPtcandID();
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &pfTauLoose();
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &pfTauLoose_leadPtcand();
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &pfcandOS10();
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &pfcandOS10looseZ();
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &pfcand5looseZ();
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &pfcanddir10();
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &pfcandveto10();
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &pfcandvetoL10();
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &jet();
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &nonisoel();
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &nonisomu();
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &t();
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &tbar();
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &ttbar();
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &lep_t();
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &lep_tbar();
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &stop_t();
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &stop_tbar();
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &neutralino_t();
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &neutralino_tbar();
	int &lep_t_id();
	int &lep_tbar_id();
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &pfjets();
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &pfjets_genJet_();
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &pfjets_failjetid();
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &pfjets_faillepolap();
	vector<float> &pfjets_csv();
	vector<float> &pfjets_chEfrac();
	vector<float> &pfjets_chm();
	vector<float> &pfjets_neu();
	vector<float> &pfjets_l1corr();
	vector<float> &pfjets_corr();
	vector<int> &pfjets_mc3();
	vector<int> &pfjets_mcflavorAlgo();
	vector<int> &pfjets_mcflavorPhys();
	vector<float> &pfjets_uncertainty();
	vector<int> &pfjets_flav();
	vector<float> &pfjets_lrm();
	vector<float> &pfjets_lrm2();
	vector<float> &pfjets_qgtag();
	vector<float> &pfjets_genJetDr();
	vector<float> &pfjets_sigma();
	vector<int> &pfjets_lepjet();
	vector<float> &pfjets_tobtecmult();
	vector<float> &pfjets_tobtecfrac();
	vector<float> &pfjets_beta();
	vector<float> &pfjets_beta2();
	vector<float> &pfjets_beta_0p1();
	vector<float> &pfjets_beta_0p2();
	vector<float> &pfjets_beta2_0p1();
	vector<float> &pfjets_beta2_0p5();
	vector<float> &pfjets_mvaPUid();
	vector<float> &pfjets_mva5xPUid();
	vector<float> &pfjets_mvaBeta();
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &genjets();
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &genqgs();
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &genbs();
	vector<int> &genps_pdgId();
	vector<int> &genps_firstMother();
	vector<float> &genps_energy();
	vector<float> &genps_pt();
	vector<float> &genps_eta();
	vector<float> &genps_phi();
	vector<float> &genps_mass();
}
#endif
