#ifndef PROTONINEL_H
#define PROTONINEL_H

#include "TFile.h"
#include "SliceParams.h"
#include "HadAna.h"
#include "Unfold.h"
#include "BetheBloch.h"
#include "TRandom3.h"

class anavar;
class HadAna;

class ProtonInel {

 public:

  ProtonInel();

  HadAna hadana;
  BetheBloch bb;
  BetheBloch bb_mu;
  TRandom3 *grdm = new TRandom3(1); // fixed seed

  int reco_sliceID;
  int reco_end_sliceID; // differ with reco_sliceID if reco_beam_calo_endZ > 220
  int true_sliceID;
  int reco_ini_sliceID;
  int true_ini_sliceID;
  double ff_energy_reco;
  double ff_energy_true;
  double ini_energy_reco;
  double ini_energy_true;
  double int_energy_reco;
  double int_energy_true;
  double beam_inst_P;

  bool isTestSample;

  bool selectCosmics;

  TH1D *reco_incE[p::nthinslices];
  TH1D *true_incE[p::nthinslices];
  TH1D *reco_AngCorr;
  TH1D *true_AngCorr;
  
  TH3D *h_true3Dsliceid_pion_all;
  TH1D *h_truesliceid_pion_all;
  TH1D *h_trueinisliceid_pion_all;
  TH1D *h_truesliceid_pion_uf;
  //TH1D *h_trueinisliceid_pion_uf;
  TH1D *h_truesliceid_pion_cuts;
  TH1D *h_trueinisliceid_pion_cuts;
  TH1D *h_truesliceid_pioninelastic_all;
  TH1D *h_truesliceid_pioninelastic_uf;
  TH1D *h_truesliceid_pioninelastic_cuts;
  TH1D *h_recosliceid_allevts_cuts;
  TH1D *h_recoinisliceid_allevts_cuts;
  TH1D *h_recosliceid_pion_cuts;
  TH1D *h_recoinisliceid_pion_cuts;
  TH1D *h_recosliceid_pioninelastic_cuts;

  double true_interactions[p::true_nbins];
  double true_incidents[p::true_nbins];

  // energy used to calculate slice ID
  TH1D *hbeam_inst_P[p::nCuts][p::nIntTypes+1];
  TH1D *hini_recoE[p::nCuts][p::nIntTypes+1];
  TH1D *hint_recoE[p::nCuts][p::nIntTypes+1];
  TH1D *hini_trueE[p::nCuts][p::nIntTypes+1];
  TH1D *hint_trueE[p::nCuts][p::nIntTypes+1];
  
  TH1D *hreco_beam_type[p::nCuts][p::nIntTypes+1];
  TH1D *hreco_reconstructable_beam_event[p::nCuts][p::nIntTypes+1];
  
  TH1D *htrue_beam_endZ[p::nCuts][p::nIntTypes+1];
  TH1D *hreco_beam_endZ[p::nCuts][p::nIntTypes+1];
  TH1D *hreco_true_beam_endZ[p::nCuts][p::nIntTypes+1];
  TH2D *hreco_vs_true_beam_endZ[p::nCuts][p::nIntTypes+1];
  TH2D *hreco_true_vs_true_beam_endZ[p::nCuts][p::nIntTypes+1];

  // after SCE correction
  TH1D *htrue_beam_endZ_SCE[p::nCuts][p::nIntTypes+1];
  TH1D *hreco_beam_endZ_SCE[p::nCuts][p::nIntTypes+1];
  TH1D *hreco_true_beam_endZ_SCE[p::nCuts][p::nIntTypes+1];
  TH2D *hreco_vs_true_beam_endZ_SCE[p::nCuts][p::nIntTypes+1];
  TH2D *hreco_true_vs_true_beam_endZ_SCE[p::nCuts][p::nIntTypes+1];

  TH1D *htrue_sliceID[p::nCuts][p::nIntTypes+1];
  TH1D *hreco_sliceID[p::nCuts][p::nIntTypes+1];
  //TH1D *htrue_inisliceID[p::nCuts][p::nIntTypes+1];
  TH1D *hreco_incsliceID[p::nCuts][p::nIntTypes+1];
  TH1D *hreco_inisliceID[p::nCuts][p::nIntTypes+1];
  TH3D *hreco_3DsliceID[p::nCuts][p::nIntTypes+1];
  TH1D *hreco_true_sliceID[p::nCuts][p::nIntTypes+1];
  TH2D *hreco_vs_true_sliceID[p::nCuts][p::nIntTypes+1];
  TH2D *hreco_true_vs_true_sliceID[p::nCuts][p::nIntTypes+1];

  TH1D *hmediandEdx[p::nCuts][p::nIntTypes+1];
  TH1D *hdaughter_michel_score[p::nCuts][p::nIntTypes+1];
  TH1D *henergy_calorimetry_SCE[p::nCuts][p::nIntTypes+1];
  TH1D *hdEdx_SCE[p::nCuts][p::nIntTypes+1];
  TH1D *hdaughter_michel_scoreMu[p::nCuts][p::nIntTypes+1];
  TH1D *hdaughter_michel_scorePi[p::nCuts][p::nIntTypes+1];
  TH1D *hdaughter_michel_score2Mu[p::nCuts][p::nIntTypes+1];
  TH1D *htrackscore[p::nCuts][p::nIntTypes+1];
  TH1D *hemscore[p::nCuts][p::nIntTypes+1];
  TH1D *hdEdx_5cm[p::nCuts][p::nIntTypes+1];

  TH1D *hmediandEdx_bkg[p::nCuts][p::nIntTypes+1];
  TH1D *hChi2_proton_bkg[p::nCuts][p::nIntTypes+1];
  TH1D *hdaughter_michel_score_bkg[p::nCuts][p::nIntTypes+1];
  TH1D *hcostheta_bkg[p::nCuts][p::nIntTypes+1];
  TH1D *hmediandEdxSlice[p::reco_nbins][p::nCuts][p::nIntTypes+1];
  TH1D *hChi2_protonSlice[p::reco_nbins][p::nCuts][p::nIntTypes+1];
  TH1D *hdaughter_michel_scoreSlice[p::reco_nbins][p::nCuts][p::nIntTypes+1];
  TH1D *hcosthetaSlice[p::reco_nbins][p::nCuts][p::nIntTypes+1];

  TH1D *hdeltax[p::nCuts][p::nIntTypes+1];
  TH1D *hdeltay[p::nCuts][p::nIntTypes+1];
  TH1D *hdeltaz[p::nCuts][p::nIntTypes+1];
  TH1D *hcostheta[p::nCuts][p::nIntTypes+1];
  TH1D *hdeltax_inst[p::nCuts][p::nIntTypes+1];
  TH1D *hdeltay_inst[p::nCuts][p::nIntTypes+1];
  TH2D *hxy_inst[p::nCuts][p::nIntTypes+1];

  TH1D *hreco_beam_true_byE_matched[p::nCuts][p::nIntTypes+1];
  TH1D *hreco_trklen[p::nCuts][p::nIntTypes+1];
  TH1D *htrue_trklen[p::nCuts][p::nIntTypes+1];
  TH1D *hdiff_trklen[p::nCuts][p::nIntTypes+1];
  TH2D *hreco_vs_true_trklen[p::nCuts][p::nIntTypes+1];
  TH1D *hbeam_score[p::nCuts][p::nIntTypes+1];
  TH2D *beam_score_vs_hreco_trklen[p::nCuts][p::nIntTypes+1];

  TH1D *hreco_beam_startX_SCE[p::nCuts][p::nIntTypes+1];
  TH1D *hreco_beam_startY_SCE[p::nCuts][p::nIntTypes+1];
  TH1D *hreco_beam_startZ_SCE[p::nCuts][p::nIntTypes+1];
  TH1D *hreco_beam_dcosX_SCE[p::nCuts][p::nIntTypes+1];
  TH1D *hreco_beam_dcosY_SCE[p::nCuts][p::nIntTypes+1];
  TH1D *hreco_beam_dcosZ_SCE[p::nCuts][p::nIntTypes+1];
  TH1D *hreco_beam_angleX_SCE[p::nCuts][p::nIntTypes+1];
  TH1D *hreco_beam_angleY_SCE[p::nCuts][p::nIntTypes+1];
  TH1D *hreco_beam_angleZ_SCE[p::nCuts][p::nIntTypes+1];

  TH2D *hreco_beam_startXY_SCE[p::nCuts][p::nIntTypes+1];
  
  TH1D *htrklen_csda_proton[p::nCuts][p::nIntTypes+1];
  TH1D *hChi2_proton[p::nCuts][p::nIntTypes+1];
  TH1D *h_diff_reco_true_Eint[p::nCuts][p::nIntTypes+1];
  TH2D *h_diff_reco_true_vs_true_Eint[p::nCuts][p::nIntTypes+1];
  TProfile *pf_diff_reco_true_vs_true_Eint[p::nCuts][p::nIntTypes+1];
  TH2D *h_stopPID_trklenCSDA[p::nCuts][p::nIntTypes+1];
  
  TH2D *hreco_iniE_trklen[p::nCuts][p::nIntTypes+1];
  
  TH1D *h_beam_inst_KE;
  TH1D *h_true_ffKE;
  TH1D *h_upstream_Eloss;
  TH1D *h_upstream_Eloss_mu;
  TH2D *h_upstream_Eloss_vs_true_Eff;
  TH2D *h_upstream_Eloss_vs_Einst;
  TH2D *h_diff_startKE_vs_Einst;
  TH2D *h_true_upstream_Eloss;
  TH1D *h_diff_Eint;
  TH2D *h_diff_Eint_vs_true_Eint;
  TH1D *h_trklen_noint;
  TH1D *h_trklen_noint_mu;
  TH1D *h_test1, *h_test2, *h_test3;
  TH1D *h_diff_Ebeam, *h_diff_Eloss, *h_diff_Edepo;
  

  std::string fOutputFileName;
  TFile *outputFile;
  void SetOutputFileName(std::string name){fOutputFileName = name;};
  void BookHistograms();
  void FillHistograms(int cut, const anavar & evt, double weight=1.);
  void FillSliceHist(const anavar & evt, int constraint_type, double weight=1., int cut=6);
  void SaveHistograms();

  void ProcessEvent(const anavar & evt, Unfold & uf, double weight, double g4rw, double bkgw, double rdm_Eshift, double rdm_Eresol);
  void CalcXS(const Unfold & uf);

  void Run(anavar & evt, Unfold & uf, Long64_t nentries=-1, bool random=false, bool savetree=false);

  void SetSelectCosmics(bool sc);

};

#endif
