#include "ThinSlice.h"
#include "HadAna.h"
#include "anavar.h"
#include "TGraphErrors.h"
#include "TVector3.h"
#include "TRandom3.h"
#include "RooUnfoldBayes.h"
#include "RooUnfoldSvd.h"
#include "util.h"
#include "BetheBloch.h"

#include <iostream>

ThinSlice::ThinSlice(){
  hadana.InitPi();
  selectCosmics = false;
  bb.SetPdgCode(211);
}

void ThinSlice::BookHistograms(){

  outputFile = TFile::Open(fOutputFileName.c_str(), "recreate");
  
  for (int i = 0; i<pi::nthinslices; ++i){ // energy distribution in each thin slice
    reco_incE[i] = new TH1D(Form("reco_incE_%d",i),Form("Reco incident energy, %.1f < z < %.1f (cm)",i*pi::thinslicewidth, (i+1)*pi::thinslicewidth), pi::nbinse, 0, 1200.);
    true_incE[i] = new TH1D(Form("true_incE_%d",i),Form("True incident energy, %.1f < z < %.1f (cm)",i*pi::thinslicewidth, (i+1)*pi::thinslicewidth), pi::nbinse, 0, 1200.);
    reco_incE[i]->Sumw2();
    true_incE[i]->Sumw2();
  }

  reco_AngCorr = new TH1D("reco_AngCorr","Reco angle correction", 100, 0, 1.);
  true_AngCorr = new TH1D("true_AngCorr","true angle correction", 100, 0, 1.);
  reco_AngCorr->Sumw2();
  true_AngCorr->Sumw2();

  h_truesliceid_pion_all = new TH1D("h_truesliceid_pion_all","h_truesliceid_pion_all;True SliceID", pi::nthinslices + 2, -1, pi::nthinslices + 1);
  h_trueinisliceid_pion_all = new TH1D("h_trueinisliceid_pion_all","h_trueinisliceid_pion_all;True SliceID", pi::nthinslices + 2, -1, pi::nthinslices + 1);
  h_truesliceid_pion_cuts = new TH1D("h_truesliceid_pion_cuts","h_truesliceid_pion_cuts;True SliceID", pi::nthinslices + 2, -1, pi::nthinslices + 1);
  h_trueinisliceid_pion_cuts = new TH1D("h_trueinisliceid_pion_cuts","h_trueinisliceid_pion_cuts;True SliceID", pi::nthinslices + 2, -1, pi::nthinslices + 1);
  h_truesliceid_pioninelastic_all = new TH1D("h_truesliceid_pioninelastic_all","h_truesliceid_pioninelastic_all;True SliceID", pi::nthinslices + 2, -1, pi::nthinslices + 1);
  h_truesliceid_pioninelastic_cuts = new TH1D("h_truesliceid_pioninelastic_cuts","h_truesliceid_pioninelastic_cuts;True SliceID", pi::nthinslices + 2, -1, pi::nthinslices + 1);
  h_recosliceid_allevts_cuts = new TH1D("h_recosliceid_allevts_cuts","h_recosliceid_allevts_cuts;Reco SliceID", pi::nthinslices + 2, -1, pi::nthinslices + 1);
  h_recoinisliceid_allevts_cuts = new TH1D("h_recoinisliceid_allevts_cuts","h_recoinisliceid_allevts_cuts;Reco SliceID", pi::nthinslices + 2, -1, pi::nthinslices + 1);
  h_recosliceid_pion_cuts = new TH1D("h_recosliceid_pion_cuts","h_recosliceid_pion_cuts;Reco SliceID", pi::nthinslices + 2, -1, pi::nthinslices + 1);
  h_recoinisliceid_pion_cuts = new TH1D("h_recoinisliceid_pion_cuts","h_recoinisliceid_pion_cuts;Reco SliceID", pi::nthinslices + 2, -1, pi::nthinslices + 1);
  h_recosliceid_pioninelastic_cuts = new TH1D("h_recosliceid_pioninelastic_cuts","h_recosliceid_pioninelastic_cuts;Reco SliceID", pi::nthinslices + 2, -1, pi::nthinslices + 1);

  h_truesliceid_pion_all->Sumw2();
  h_trueinisliceid_pion_all->Sumw2();
  h_truesliceid_pion_cuts->Sumw2();
  h_trueinisliceid_pion_cuts->Sumw2();
  h_truesliceid_pioninelastic_all->Sumw2();
  h_truesliceid_pioninelastic_cuts->Sumw2();
  h_recosliceid_allevts_cuts->Sumw2();
  h_recoinisliceid_allevts_cuts->Sumw2();
  h_recosliceid_pion_cuts->Sumw2();
  h_recoinisliceid_pion_cuts->Sumw2();
  h_recosliceid_pioninelastic_cuts->Sumw2();

  for (int i = 0; i < pi::nCuts; ++i){
    for (int j = 0; j < pi::nIntTypes+1; ++j){
      hreco_beam_type[i][j] = new TH1D(Form("hreco_beam_type_%d_%d",i,j),Form("hreco_beam_type, %s, %s;hreco_beam_type", pi::cutName[i], pi::intTypeName[j]), 21, -1, 20);
      hreco_beam_type[i][j]->Sumw2();
      hreco_reconstructable_beam_event[i][j] = new TH1D(Form("hreco_reconstructable_beam_event_%d_%d",i,j),Form("hreco_reconstructable_beam_event, %s, %s;hreco_reconstructable_beam_event", pi::cutName[i], pi::intTypeName[j]), 21, -1, 20);
      hreco_reconstructable_beam_event[i][j]->Sumw2();
      
      htrue_beam_endZ[i][j] = new TH1D(Form("htrue_beam_endZ_%d_%d",i,j),Form("true_beam_endZ, %s, %s;true_beam_endZ (cm)", pi::cutName[i], pi::intTypeName[j]), 70, -100, 600);
      htrue_beam_endZ[i][j]->Sumw2();
      hreco_beam_endZ[i][j] = new TH1D(Form("hreco_beam_endZ_%d_%d",i,j),Form("reco_beam_endZ, %s, %s;reco_beam_endZ (cm)", pi::cutName[i], pi::intTypeName[j]), 70, -100, 600);
      hreco_beam_endZ[i][j]->Sumw2();
      hreco_true_beam_endZ[i][j] = new TH1D(Form("hreco_true_beam_endZ_%d_%d",i,j), Form("reco_true_beam_endZ, %s, %s;reco_beam_endZ - true_beam_endZ (cm)", pi::cutName[i], pi::intTypeName[j]), 100, -100, 100);
      hreco_true_beam_endZ[i][j]->Sumw2();
      hreco_vs_true_beam_endZ[i][j]= new TH2D(Form("hreco_vs_true_beam_endZ_%d_%d",i,j), Form("%s, %s;true_beam_endZ (cm);reco_beam_endZ (cm)", pi::cutName[i], pi::intTypeName[j]), 70, -100, 600, 70, -100, 600);
      hreco_true_vs_true_beam_endZ[i][j]= new TH2D(Form("hreco_true_vs_true_beam_endZ_%d_%d",i,j), Form("%s, %s;true_beam_endZ (cm);reco - true_beam_endZ (cm)", pi::cutName[i], pi::intTypeName[j]), 70, -100, 600, 100, -100, 100);

      htrue_beam_endZ_SCE[i][j] = new TH1D(Form("htrue_beam_endZ_SCE_%d_%d",i,j),Form("true_beam_endZ_SCE, %s, %s;true_beam_endZ_SCE (cm)", pi::cutName[i], pi::intTypeName[j]), 70, -100, 600);
      htrue_beam_endZ_SCE[i][j]->Sumw2();
      hreco_beam_endZ_SCE[i][j] = new TH1D(Form("hreco_beam_endZ_SCE_%d_%d",i,j),Form("reco_beam_endZ_SCE, %s, %s;reco_beam_endZ_SCE (cm)", pi::cutName[i], pi::intTypeName[j]), 70, -100, 600);
      hreco_beam_endZ_SCE[i][j]->Sumw2();
      hreco_true_beam_endZ_SCE[i][j] = new TH1D(Form("hreco_true_beam_endZ_SCE_%d_%d",i,j), Form("reco_true_beam_endZ_SCE, %s, %s;reco_beam_endZ_SCE - true_beam_endZ_SCE (cm)", pi::cutName[i], pi::intTypeName[j]), 100, -100, 100);
      hreco_true_beam_endZ_SCE[i][j]->Sumw2();
      hreco_vs_true_beam_endZ_SCE[i][j]= new TH2D(Form("hreco_vs_true_beam_endZ_SCE_%d_%d",i,j), Form("%s, %s;true_beam_endZ_SCE (cm);reco_beam_endZ_SCE (cm)", pi::cutName[i], pi::intTypeName[j]), 70, -100, 600, 70, -100, 600);
      hreco_true_vs_true_beam_endZ_SCE[i][j]= new TH2D(Form("hreco_true_vs_true_beam_endZ_SCE_%d_%d",i,j), Form("%s, %s;true_beam_endZ_SCE (cm);reco - true_beam_endZ_SCE (cm)", pi::cutName[i], pi::intTypeName[j]), 70, -100, 600, 100, -100, 100);

      htrue_sliceID[i][j] = new TH1D(Form("htrue_sliceID_%d_%d",i,j),Form("true_sliceID, %s, %s;true_sliceID (cm)", pi::cutName[i], pi::intTypeName[j]), pi::nthinslices+2, -1, pi::nthinslices+1);
      htrue_sliceID[i][j]->Sumw2();
      hreco_sliceID[i][j] = new TH1D(Form("hreco_sliceID_%d_%d",i,j),Form("reco_sliceID, %s, %s;reco_sliceID", pi::cutName[i], pi::intTypeName[j]), pi::nthinslices+2, -1, pi::nthinslices+1);
      hreco_sliceID[i][j]->Sumw2();
      hreco_incsliceID[i][j] = new TH1D(Form("hreco_incsliceID_%d_%d",i,j),Form("reco_incsliceID, %s, %s;reco_incsliceID", pi::cutName[i], pi::intTypeName[j]), pi::nthinslices+2, -1, pi::nthinslices+1);
      hreco_incsliceID[i][j]->Sumw2();
      hreco_inisliceID[i][j] = new TH1D(Form("hreco_inisliceID_%d_%d",i,j),Form("reco_inisliceID, %s, %s;reco_inisliceID", pi::cutName[i], pi::intTypeName[j]), pi::nthinslices+2, -1, pi::nthinslices+1);
      hreco_inisliceID[i][j]->Sumw2();
      hreco_true_sliceID[i][j] = new TH1D(Form("hreco_true_sliceID_%d_%d",i,j), Form("reco_true_sliceID, %s, %s;reco_sliceID - true_sliceID", pi::cutName[i], pi::intTypeName[j]), 20, -10, 10);
      hreco_true_sliceID[i][j]->Sumw2();
      hreco_vs_true_sliceID[i][j]= new TH2D(Form("hreco_vs_true_sliceID_%d_%d",i,j), Form("%s, %s;true_sliceID;reco_sliceID", pi::cutName[i], pi::intTypeName[j]), pi::nthinslices+2, -1, pi::nthinslices+1, pi::nthinslices+2, -1, pi::nthinslices+1);
      hreco_true_vs_true_sliceID[i][j]= new TH2D(Form("hreco_true_vs_true_sliceID_%d_%d",i,j), Form("%s, %s;true_sliceID;reco_sliceID - true_sliceID", pi::cutName[i], pi::intTypeName[j]), pi::nthinslices+2, -1, pi::nthinslices+1, 20, -10, 10);

      hmediandEdx[i][j] = new TH1D(Form("hmediandEdx_%d_%d",i,j), Form("mediandEdx, %s, %s;Median dE/dx (MeV/cm)", pi::cutName[i], pi::intTypeName[j]), 100, 0, 5);
      hmediandEdx[i][j]->Sumw2();

      hdaughter_michel_score[i][j] = new TH1D(Form("hdaughter_michel_score_%d_%d",i,j), Form("daughter_michel_score, %s, %s;Michel score", pi::cutName[i], pi::intTypeName[j]), 110, -0.1, 1);
      hdaughter_michel_score[i][j]->Sumw2();
      
      henergy_calorimetry_SCE[i][j] = new TH1D(Form("henergy_calorimetry_SCE_%d_%d",i,j), Form("Energy_calorimetry_SCE_corrected, %s, %s;Energy (MeV)", pi::cutName[i], pi::intTypeName[j]), 100, 0, 500);
      henergy_calorimetry_SCE[i][j]->Sumw2();
      hdEdx_SCE[i][j] = new TH1D(Form("hdEdx_SCE_%d_%d",i,j), Form("dEdx_SCE_corrected, %s, %s; dE/dx (MeV/cm)", pi::cutName[i], pi::intTypeName[j]), 100, 1.7, 4.2);
      hdEdx_SCE[i][j]->Sumw2();

      hdaughter_michel_scoreMu[i][j] = new TH1D(Form("hdaughter_michel_scoreMu_%d_%d",i,j), Form("daughter_michel_scoreMu, %s, %s;Michel score", pi::cutName[i], pi::intTypeName[j]), 10, 0, 1);
      hdaughter_michel_scoreMu[i][j]->Sumw2();

      hdaughter_michel_score2Mu[i][j] = new TH1D(Form("hdaughter_michel_score2Mu_%d_%d",i,j), Form("daughter_michel_score2Mu, %s, %s;Michel score", pi::cutName[i], pi::intTypeName[j]), 10, 0, 1);
      hdaughter_michel_score2Mu[i][j]->Sumw2();

      hdaughter_michel_scorePi[i][j] = new TH1D(Form("hdaughter_michel_scorePi_%d_%d",i,j), Form("daughter_michel_scorePi, %s, %s;Michel score", pi::cutName[i], pi::intTypeName[j]), 10, 0, 1);
      hdaughter_michel_scorePi[i][j]->Sumw2();

      hmediandEdx_bkg[i][j] = new TH1D(Form("hmediandEdx_bkg_%d_%d",i,j), Form("mediandEdx_bkg, %s, %s;Median dE/dx (MeV/cm)", pi::cutName[i], pi::intTypeName[j]), 100, 0, 8);
      hmediandEdx_bkg[i][j]->Sumw2();
      hChi2_proton_bkg[i][j] = new TH1D(Form("hChi2_proton_bkg_%d_%d",i,j), Form("Chi2_proton_bkg, %s, %s;Chi2/Ndof", pi::cutName[i], pi::intTypeName[j]), 100, 0, 100);
      hChi2_proton_bkg[i][j]->Sumw2();
      hdaughter_michel_score_bkg[i][j] = new TH1D(Form("hdaughter_michel_score_bkg_%d_%d",i,j), Form("daughter_michel_score_bkg, %s, %s;Michel score", pi::cutName[i], pi::intTypeName[j]), 100, 0, 1);
      hdaughter_michel_score_bkg[i][j]->Sumw2();
      hcostheta_bkg[i][j] = new TH1D(Form("hcostheta_bkg_%d_%d",i,j), Form("costheta_bkg, %s, %s;cos#theta", pi::cutName[i], pi::intTypeName[j]), 75, 0.85, 1);
      hcostheta_bkg[i][j]->Sumw2();
      for (int k = 0; k<pi::nthinslices; ++k){
        hmediandEdxSlice[k][i][j] = new TH1D(Form("hmediandEdxSlice_%d_%d_%d",k,i,j), Form("mediandEdx, %s, %s, sliceID = %d;Median dE/dx (MeV/cm)", pi::cutName[i], pi::intTypeName[j], k), 14, 1, 8);
        hmediandEdxSlice[k][i][j]->Sumw2();
        
        hChi2_protonSlice[k][i][j] = new TH1D(Form("hChi2_protonSlice_%d_%d_%d",k,i,j), Form("Chi2_proton, %s, %s, sliceID = %d;Chi2_proton/Ndf", pi::cutName[i], pi::intTypeName[j], k), 10, 0, 100);
        hChi2_protonSlice[k][i][j]->Sumw2();

        hdaughter_michel_scoreSlice[k][i][j] = new TH1D(Form("hdaughter_michel_scoreSlice_%d_%d_%d",k,i,j), Form("daughter_michel_score, %s, %s, sliceID = %d;Michel score", pi::cutName[i], pi::intTypeName[j], k), 10, 0, 1);
        hdaughter_michel_scoreSlice[k][i][j]->Sumw2();
        
        hcosthetaSlice[k][i][j] = new TH1D(Form("hcosthetaSlice_%d_%d_%d",k,i,j), Form("costheta, %s, %s, sliceID = %d;Cos(theta)", pi::cutName[i], pi::intTypeName[j], k), 15, 0.85, 1);
        hcosthetaSlice[k][i][j]->Sumw2();
      }        

      htrackscore[i][j] = new TH1D(Form("htrackscore_%d_%d",i,j), Form("trackscore, %s, %s;Track score", pi::cutName[i], pi::intTypeName[j]), 110, -0.1, 1);
      htrackscore[i][j]->Sumw2();

      hemscore[i][j] = new TH1D(Form("hemscore_%d_%d",i,j), Form("emscore, %s, %s;Em score", pi::cutName[i], pi::intTypeName[j]), 50, 0, 1);
      hemscore[i][j]->Sumw2();

      hdEdx_5cm[i][j] = new TH1D(Form("hdEdx_5cm_%d_%d",i,j), Form("dEdx_5cm, %s, %s;dE/dx_5cm (MeV/cm)", pi::cutName[i], pi::intTypeName[j]), 100, 0, 5);
      hdEdx_5cm[i][j]->Sumw2();

      hdeltax[i][j] = new TH1D(Form("hdeltax_%d_%d",i,j), Form("deltax, %s, %s;#Deltax/#sigma_{x}", pi::cutName[i], pi::intTypeName[j]), 100, -10, 10);
      hdeltax[i][j]->Sumw2();

      hdeltay[i][j] = new TH1D(Form("hdeltay_%d_%d",i,j), Form("deltay, %s, %s;#Deltay/#sigma_{y}", pi::cutName[i], pi::intTypeName[j]), 100, -10, 10);
      hdeltay[i][j]->Sumw2();
      
      hdeltax_inst[i][j] = new TH1D(Form("hdeltax_inst_%d_%d",i,j), Form("deltax_inst, %s, %s;#Deltax_inst/#sigma_{x}", pi::cutName[i], pi::intTypeName[j]), 100, -45, -15);
      hdeltax_inst[i][j]->Sumw2();
      hdeltay_inst[i][j] = new TH1D(Form("hdeltay_inst_%d_%d",i,j), Form("deltay_inst, %s, %s;#Deltay_inst/#sigma_{y}", pi::cutName[i], pi::intTypeName[j]), 100, 405, 440);
      hdeltay_inst[i][j]->Sumw2();
      hxy_inst[i][j] = new TH2D(Form("hxy_inst_%d_%d",i,j), Form("%s, %s;x (cm);y (cm)", pi::cutName[i], pi::intTypeName[j]), 100, -45, -15, 100, 405, 440);

      hdeltaz[i][j] = new TH1D(Form("hdeltaz_%d_%d",i,j), Form("deltaz, %s, %s;#Deltaz/#sigma_{z}", pi::cutName[i], pi::intTypeName[j]), 100, -10, 10);
      hdeltaz[i][j]->Sumw2();

      hcostheta[i][j] = new TH1D(Form("hcostheta_%d_%d",i,j), Form("costheta, %s, %s;cos#theta", pi::cutName[i], pi::intTypeName[j]), 100, 0.9, 1);
      hcostheta[i][j]->Sumw2();

      hreco_beam_true_byE_matched[i][j] = new TH1D(Form("hreco_beam_true_byE_matched_%d_%d",i,j), Form("reco_beam_true_byE_matched, %s, %s;Truth matched", pi::cutName[i], pi::intTypeName[j]), 2, 0, 2);
      hreco_beam_true_byE_matched[i][j]->Sumw2();
      //const double xbins[25] = {-10.,0.,10.,20.,30.,40.,50.,60.,70.,80.,90.,100.,110.,120.,130.,140.,150.,160.,170.,180.,190.,200.,210.,220.,600}; // a user-defined binning
      //hreco_trklen[i][j] = new TH1D(Form("hreco_trklen_%d_%d",i,j), Form("reco_trklen, %s, %s;Track length (cm)", pi::cutName[i], pi::intTypeName[j]), 24, xbins);
      hreco_trklen[i][j] = new TH1D(Form("hreco_trklen_%d_%d",i,j), Form("reco_trklen, %s, %s;Track length (cm)", pi::cutName[i], pi::intTypeName[j]), 61, -10, 600);
      hreco_trklen[i][j]->Sumw2();
      htrue_trklen[i][j] = new TH1D(Form("htrue_trklen_%d_%d",i,j), Form("true_trklen, %s, %s;Track length (cm)", pi::cutName[i], pi::intTypeName[j]), 61, -10, 600);
      htrue_trklen[i][j]->Sumw2();
      hdiff_trklen[i][j] = new TH1D(Form("hdiff_trklen_%d_%d",i,j), Form("diff_trklen, %s, %s;Track length (cm)", pi::cutName[i], pi::intTypeName[j]), 60, -600, 600);
      hdiff_trklen[i][j]->Sumw2();
      hreco_vs_true_trklen[i][j]= new TH2D(Form("hreco_vs_true_trklen_%d_%d",i,j), Form("%s, %s;true_trklen (cm);reco_trklen (cm)", pi::cutName[i], pi::intTypeName[j]), 61, -10, 600, 61, -10, 600);
      hbeam_score[i][j] = new TH1D(Form("hbeam_score_%d_%d",i,j), Form("Beam_score, %s, %s;Beam score", pi::cutName[i], pi::intTypeName[j]), 140, -0.23, 0.12);
      hbeam_score[i][j]->Sumw2();
      beam_score_vs_hreco_trklen[i][j]= new TH2D(Form("beam_score_vs_hreco_trklen_%d_%d",i,j), Form("%s, %s;reco trklen (cm);beam_score", pi::cutName[i], pi::intTypeName[j]), 50, -100, 400, 140, -0.23, 0.12);

      hreco_beam_startX_SCE[i][j] = new TH1D(Form("hreco_beam_startX_SCE_%d_%d",i,j), Form("reco_beam_startX_SCE, %s, %s; reco_beam_startX_SCE (cm)", pi::cutName[i], pi::intTypeName[j]), 100, -80, 20);
      hreco_beam_startX_SCE[i][j]->Sumw2();

      hreco_beam_startY_SCE[i][j] = new TH1D(Form("hreco_beam_startY_SCE_%d_%d",i,j), Form("reco_beam_startY_SCE, %s, %s; reco_beam_startY_SCE (cm)", pi::cutName[i], pi::intTypeName[j]), 100, 350, 500);
      hreco_beam_startY_SCE[i][j]->Sumw2();

      hreco_beam_startZ_SCE[i][j] = new TH1D(Form("hreco_beam_startZ_SCE_%d_%d",i,j), Form("reco_beam_startZ_SCE, %s, %s; reco_beam_startZ_SCE (cm)", pi::cutName[i], pi::intTypeName[j]), 100, -5, 10);
      hreco_beam_startZ_SCE[i][j]->Sumw2();

      hreco_beam_dcosX_SCE[i][j] = new TH1D(Form("hreco_beam_dcosX_SCE_%d_%d",i,j), Form("hreco_beam_dcosX_SCE, %s, %s; reco_beam_dcosX_SCE", pi::cutName[i], pi::intTypeName[j]), 100, -1, 1);
      hreco_beam_dcosX_SCE[i][j]->Sumw2();

      hreco_beam_dcosY_SCE[i][j] = new TH1D(Form("hreco_beam_dcosY_SCE_%d_%d",i,j), Form("hreco_beam_dcosY_SCE, %s, %s; reco_beam_dcosY_SCE", pi::cutName[i], pi::intTypeName[j]), 100, -1, 1);
      hreco_beam_dcosY_SCE[i][j]->Sumw2();

      hreco_beam_dcosZ_SCE[i][j] = new TH1D(Form("hreco_beam_dcosZ_SCE_%d_%d",i,j), Form("hreco_beam_dcosZ_SCE, %s, %s; reco_beam_dcosZ_SCE", pi::cutName[i], pi::intTypeName[j]), 100, -1, 1);
      hreco_beam_dcosZ_SCE[i][j]->Sumw2();

      hreco_beam_angleX_SCE[i][j] = new TH1D(Form("hreco_beam_angleX_SCE_%d_%d",i,j), Form("hreco_beam_angleX_SCE, %s, %s; #theta_{x} (deg)", pi::cutName[i], pi::intTypeName[j]), 180, 0, 180);
      hreco_beam_angleX_SCE[i][j]->Sumw2();

      hreco_beam_angleY_SCE[i][j] = new TH1D(Form("hreco_beam_angleY_SCE_%d_%d",i,j), Form("hreco_beam_angleY_SCE, %s, %s; #theta_{y} (deg)", pi::cutName[i], pi::intTypeName[j]), 180, 0, 180);
      hreco_beam_angleY_SCE[i][j]->Sumw2();

      hreco_beam_angleZ_SCE[i][j] = new TH1D(Form("hreco_beam_angleZ_SCE_%d_%d",i,j), Form("hreco_beam_angleZ_SCE, %s, %s; #theta_{z} (deg)", pi::cutName[i], pi::intTypeName[j]), 90, 0, 90);
      hreco_beam_angleZ_SCE[i][j]->Sumw2();

      hreco_beam_startXY_SCE[i][j] = new TH2D(Form("hreco_beam_startXY_SCE_%d_%d",i,j), Form("reco_beam_startXY_SCE, %s, %s;reco_beam_startX_SCE (cm);reco_beam_startY_SCE (cm)", pi::cutName[i], pi::intTypeName[j]), 1000, -360, 360, 1000, 0, 700);
      
      htrklen_csda_proton[i][j] = new TH1D(Form("htrklen_csda_proton_%d_%d",i,j), Form("trklen_csda_proton, %s, %s;Track length / CSDA", pi::cutName[i], pi::intTypeName[j]), 61, -0.1, 6);
      htrklen_csda_proton[i][j]->Sumw2();
      hChi2_proton[i][j] = new TH1D(Form("hChi2_proton_%d_%d",i,j), Form("Chi2_proton, %s, %s;Chi2/Ndof", pi::cutName[i], pi::intTypeName[j]), 101, -1, 100);
      hChi2_proton[i][j]->Sumw2();
      h_diff_reco_true_Eint[i][j] = new TH1D(Form("h_diff_reco_true_Eint_%d_%d",i,j), Form("h_diff_reco_true_Eint, %s, %s;MeV", pi::cutName[i], pi::intTypeName[j]), 100, -300, 300);
      h_diff_reco_true_Eint[i][j]->Sumw2();
      h_diff_reco_true_vs_true_Eint[i][j] = new TH2D(Form("h_diff_reco_true_vs_true_Eint_%d_%d",i,j), Form("h_diff_reco_true_vs_true_Eint, %s, %s;MeV", pi::cutName[i], pi::intTypeName[j]), 100, 0, 1000, 100, -300, 300);
      pf_diff_reco_true_vs_true_Eint[i][j] = new TProfile(Form("pf_diff_reco_true_vs_true_Eint_%d_%d",i,j), Form("pf_diff_reco_true_vs_true_Eint, %s, %s;MeV", pi::cutName[i], pi::intTypeName[j]), 100, 0, 1000);
    }
  }
  h_beam_inst_KE = new TH1D("h_beam_inst_KE","h_beam_inst_KE;MeV", 100, 0, 2000);
  h_beam_inst_KE->Sumw2();
  h_true_ffKE = new TH1D("h_true_ffKE","h_true_ffKE;MeV", 100, 0, 2000);
  h_true_ffKE->Sumw2();
  h_upstream_Eloss = new TH1D("h_upstream_Eloss","h_upstream_Eloss;MeV", 100, -100, 100);
  h_upstream_Eloss->Sumw2();
  h_upstream_Eloss_vs_true_Eff = new TH2D("h_upstream_Eloss_vs_true_Eff","h_upstream_Eloss_vs_true_Eff;MeV;MeV", 100, 600, 1100, 100, -100, 100);
  h_upstream_Eloss_vs_Einst = new TH2D("h_upstream_Eloss_vs_Einst","h_upstream_Eloss_vs_Einst;MeV;MeV", 100, 600, 1100, 100, -100, 100);
  h_diff_startKE_vs_Einst = new TH2D("h_diff_startKE_vs_Einst","h_diff_startKE_vs_Einst;MeV;MeV", 100, 600, 1100, 100, -100, 100);
  h_true_upstream_Eloss = new TH2D("h_true_upstream_Eloss","h_true_upstream_Eloss;MeV;MeV", 100, 600, 1100, 100, -100, 100);
  h_diff_Eint = new TH1D("h_diff_Eint","h_diff_Eint;MeV", 100, -100, 100);
  h_diff_Eint->Sumw2();
  h_diff_Eint_vs_true_Eint = new TH2D("h_diff_Eint_vs_true_Eint","h_diff_Eint_vs_true_Eint;MeV;MeV", 100, 0, 1000, 100, -100, 100);

   for (int i = 0; i<pi::nthinslices; ++i){
     true_interactions[i] = 0;
     true_incidents[i] = 0;
   }
   
   //response_SliceID_Pion = new RooUnfoldResponse(pi::nthinslices+2, -1, pi::nthinslices+1, "response_SliceID_Pion");
   //response_SliceID_PionInEl = new RooUnfoldResponse(pi::nthinslices+2, -1, pi::nthinslices+1, "response_SliceID_PionInEl");

}

void ThinSlice::ProcessEvent(const anavar & evt, Unfold & uf, double g4rw, double bkgw){
  //hadana.ProcessEvent(evt);
  reco_sliceID = -99;
  reco_end_sliceID = -99;
  true_sliceID = -99;
  reco_ini_sliceID = -99;
  true_ini_sliceID = -99;
  int_energy_reco = 999999.;
  int_energy_true = 999999.;

  isTestSample = (hadana.pitype == pi::kData); // fake data
  //if (evt.MC && evt.event%2 == 0) isTestSample = false;

  if (hadana.fAllTrackCheck) {} // removed for not in use
  else {//not using all track reconstruction
    double pimass = 139.57;
    double beam_inst_KE = sqrt(pow(evt.beam_inst_P*1000,2)+pow(pimass,2)) - pimass;
    if (evt.MC){
      if (hadana.true_ffKE < 999999) { // entered TPC
        if (evt.true_beam_PDG == 211) {
          h_beam_inst_KE->Fill(beam_inst_KE);
          h_true_ffKE->Fill(hadana.true_ffKE);
          h_upstream_Eloss->Fill(beam_inst_KE - hadana.true_ffKE);
          h_upstream_Eloss_vs_true_Eff->Fill(hadana.true_ffKE, beam_inst_KE - hadana.true_ffKE);
          h_upstream_Eloss_vs_Einst->Fill(beam_inst_KE, beam_inst_KE - hadana.true_ffKE);
          double true_beam_startKE = sqrt(pow(evt.true_beam_startP*1000,2)+pow(pimass,2)) - pimass;
          h_diff_startKE_vs_Einst->Fill(beam_inst_KE, beam_inst_KE - true_beam_startKE);
          h_true_upstream_Eloss->Fill(hadana.true_ffKE, true_beam_startKE - hadana.true_ffKE);
        }
        
        int traj_max = evt.true_beam_traj_Z->size()-1;
        if ((*evt.true_beam_traj_KE)[traj_max] != 0) {
          int_energy_true = (*evt.true_beam_traj_KE)[traj_max];
        }
        else {
          int temp = traj_max-1;
          while ((*evt.true_beam_traj_KE)[temp] == 0) temp--;
          //int_energy_true = bb.KEAtLength((*evt.true_beam_traj_KE)[temp], (hadana.true_trklen_accum)[traj_max]-(hadana.true_trklen_accum)[temp]);
          int_energy_true = (*evt.true_beam_traj_KE)[temp] - 2.1*((hadana.true_trklen_accum)[traj_max]-(hadana.true_trklen_accum)[temp]); // 2.1 MeV/cm
          //cout<<"int_energy_true"<<(*evt.true_beam_traj_KE)[temp]<<"\t"<<sqrt(pow(evt.true_beam_endP*1000,2)+pow(139.57,2)) - 139.57<<endl; // almost the same
        }
        double int_energy_true_trklen = bb.KEAtLength(hadana.true_ffKE, hadana.true_trklen);
        h_diff_Eint->Fill(int_energy_true_trklen - int_energy_true);
        h_diff_Eint_vs_true_Eint->Fill(int_energy_true, int_energy_true_trklen - int_energy_true);
      }
      // true initial sliceID
      true_ini_sliceID = int(ceil( (pi::plim - hadana.true_ffKE)/pi::Eslicewidth )); // ignore incomplete slices
      //if (true_ini_sliceID <= -99) true_ini_sliceID = -99;
      if (true_ini_sliceID < 0) true_ini_sliceID = -1; // both physical and unphysical underflow
      if (true_ini_sliceID >= pi::nthinslices) true_ini_sliceID = pi::nthinslices; // overflow (Eff<pi::Eslicewidth)
      // true interaction sliceID
      true_sliceID = int(floor( (pi::plim-int_energy_true)/pi::Eslicewidth ));
      //if (true_sliceID <= -99) true_sliceID = -99;
      if (true_sliceID < 0) true_sliceID = -1; // unphysical underflow
      if (true_sliceID >= pi::nthinslices) true_sliceID = pi::nthinslices; // overflow (int_energy_true <= 0)
      // ignore incomplete slices
      /*if (true_sliceID < true_ini_sliceID) {
        true_ini_sliceID = -1;
        true_sliceID = -1;
      } // if true_sliceID==-1, this event should not be used when calculating true XS (but should it be used in unfolding???)*/

      if (evt.true_beam_PDG == 211){
        int starti = true_ini_sliceID;
        if (starti == -1) starti = 0;
        for (int i = starti; i<=true_sliceID; ++i){
          if (i<pi::nthinslices) ++true_incidents[i]; // count incident events
        }
        
        if ((*evt.true_beam_endProcess) == "pi+Inelastic"){
          if (true_sliceID >= starti){
            ++true_interactions[true_sliceID]; // count interaction events
          }
          // Reco info
          if (!(evt.reco_beam_calo_wire->empty()) && evt.reco_beam_true_byE_matched){ // truth matched. so it must be the true track?
            std::vector<std::vector<double>> vincE(pi::nthinslices);
            for (size_t i = 0; i<evt.reco_beam_calo_wire->size(); ++i){
              //int this_sliceID = int((hadana.reco_trklen_accum)[i]/pi::thinslicewidth);
              int this_sliceID = int((*evt.reco_beam_incidentEnergies)[i]/pi::Eslicewidth);
              if (this_sliceID>=pi::nthinslices) continue;
              if (this_sliceID<0) continue;
              double this_incE = (*evt.reco_beam_incidentEnergies)[i];
              vincE[this_sliceID].push_back(this_incE);
            }
            for (size_t i = 0; i<vincE.size(); ++i){
              if (!vincE[i].empty()){
                double sum_incE = 0;
                for (size_t j = 0; j<vincE[i].size(); ++j){
                  sum_incE += vincE[i][j];
                }
                reco_incE[i]->Fill(sum_incE/vincE[i].size());
              }
            }
            TVector3 pt0(evt.reco_beam_calo_startX,
                         evt.reco_beam_calo_startY,
                         evt.reco_beam_calo_startZ);
            TVector3 pt1(evt.reco_beam_calo_endX,
                         evt.reco_beam_calo_endY,
                         evt.reco_beam_calo_endZ);
            TVector3 dir = pt1 - pt0; // direction of the track is determine by the start/end point
            dir = dir.Unit();
            reco_AngCorr->Fill(dir.Z()); // projection to Z of the direction of the track
          }

          // True info
          if (!(evt.true_beam_traj_Z->empty())){
            std::vector<std::vector<double>> vincE(pi::nthinslices);
            for (size_t i = 0; i<evt.true_beam_traj_Z->size()-1; ++i){//last point always has KE = 0
              //int this_sliceID = int((hadana.true_trklen_accum)[i]/pi::thinslicewidth);
              int this_sliceID = int((*evt.true_beam_traj_KE)[i]/pi::Eslicewidth);
              double this_incE = (*evt.true_beam_traj_KE)[i];
              if (this_sliceID>=pi::nthinslices) continue;
              if (this_sliceID<0) continue;
              vincE[this_sliceID].push_back(this_incE);
            }
            for (size_t i = 0; i<vincE.size(); ++i){
              if (!vincE[i].empty()){
                double sum_incE = 0;
                for (size_t j = 0; j<vincE[i].size(); ++j){
                  sum_incE += vincE[i][j];
                }
                true_incE[i]->Fill(sum_incE/vincE[i].size());
              }
            }
            TVector3 pt0(evt.true_beam_startX,
                         evt.true_beam_startY,
                         evt.true_beam_startZ);
            TVector3 pt1(evt.true_beam_endX,
                         evt.true_beam_endY,
                         evt.true_beam_endZ);
            TVector3 dir = pt1 - pt0;
            dir = dir.Unit();
            true_AngCorr->Fill(dir.Z());
          }
        }
      }
    }

    if (!evt.reco_beam_calo_wire->empty()){
      double inc_energy_reco = beam_inst_KE - 12.74;
      /*if (beam_inst_KE < 800) inc_energy_reco = beam_inst_KE - 0.95; // 0.9465 \pm 0.3051
      else if (beam_inst_KE < 850) inc_energy_reco = beam_inst_KE - 7.12; // 7.119 \pm 0.210
      else if (beam_inst_KE < 900) inc_energy_reco = beam_inst_KE - 11.87; // 11.87 \pm 0.22
      else if (beam_inst_KE < 950) inc_energy_reco = beam_inst_KE - 17.31; // 17.31 \pm 0.27
      else inc_energy_reco = beam_inst_KE - 29.28; // 29.28 \pm 0.37*/

      int_energy_reco = bb.KEAtLength(inc_energy_reco, hadana.reco_trklen);
      // reco initial sliceID
      reco_ini_sliceID = int(ceil( (pi::plim - inc_energy_reco)/pi::Eslicewidth ));
      if (reco_ini_sliceID < 0) reco_ini_sliceID = -1;
      if (reco_ini_sliceID >= pi::nthinslices) reco_ini_sliceID = pi::nthinslices;
      // reco interaction sliceID
      reco_sliceID = int(floor( (pi::plim-int_energy_reco)/pi::Eslicewidth ));
      if (reco_sliceID < 0) reco_sliceID = -1;
      //if (hadana.reco_trklen < 0) reco_sliceID = -1;
      if (reco_sliceID >= pi::nthinslices) { // overflow (int_energy_reco <= 0)
        reco_sliceID = pi::nthinslices;
        //cout<<"reco_sliceID >= pi::nthinslices"<<int_energy_reco<<endl;
        //cout<<inc_energy_reco<<"\t"<<hadana.reco_trklen<<endl;
      }
      // ignore incomplete slices
      /*if (reco_sliceID < reco_ini_sliceID) {
        reco_ini_sliceID = -1;
        reco_sliceID = -1;
      } // if reco_sliceID==-1, this event should not be used when calculating reco XS*/
      reco_end_sliceID = reco_sliceID;
      /*if (evt.reco_beam_calo_endZ > 220) { // APA3 cut
        int idx = evt.reco_beam_calo_Z->size()-1;
        for (; (*evt.reco_beam_calo_Z)[idx]>220; --idx) {}
        double energy_reco = bb.KEAtLength(inc_energy_reco, hadana.reco_trklen_accum[idx]);
        if (idx >= 0) {
          reco_end_sliceID = int(floor( (pi::plim-energy_reco)/pi::Eslicewidth )) - 1;
          reco_sliceID = pi::nthinslices;
        }
        else { // evt.reco_beam_calo_startZ > 220
          reco_ini_sliceID = -1;
          reco_end_sliceID = -1;
          reco_sliceID = -1;
        }
        //reco_ini_sliceID = -1;
        //reco_end_sliceID = -1;
        //reco_sliceID = -1;
      }*/
    }
    //for(int i=0;i<evt.reco_beam_incidentEnergies->size();i++) cout<<(*evt.reco_beam_incidentEnergies)[i]<<"\t";
    //cout<<endl<<int_energy_reco<<endl; // very different with the last point of evt.reco_beam_incidentEnergies
  }

  double weight = g4rw * bkgw;
  if (evt.MC){
    if (evt.true_beam_PDG == 211){ // true pion beam incident event
      if (isTestSample){ // fake data
        h_truesliceid_pion_all->Fill(true_sliceID, g4rw);
        h_trueinisliceid_pion_all->Fill(true_ini_sliceID, g4rw);
      }
      else{
        //uf.eff_den_Inc->Fill(true_sliceID);
      }
      if (hadana.PassPiCuts(evt) && evt.reco_beam_true_byE_matched){ // the beam pion passed full selections (reco_beam_true_byE_matched is used to veto secondary particles. Only in MC)
        if (isTestSample){
          h_recosliceid_pion_cuts->Fill(reco_end_sliceID, g4rw);
          h_truesliceid_pion_cuts->Fill(true_sliceID, g4rw);
          h_recoinisliceid_pion_cuts->Fill(reco_ini_sliceID, g4rw);
          h_trueinisliceid_pion_cuts->Fill(true_ini_sliceID, g4rw);
        }
        else{
          //uf.eff_num_Inc->Fill(true_sliceID);
          //uf.pur_num_Inc->Fill(reco_sliceID);
          uf.response_SliceID_Inc.Fill(reco_end_sliceID, true_sliceID, g4rw);
          uf.response_SliceID_Ini.Fill(reco_ini_sliceID, true_ini_sliceID, g4rw);
        }
      }
      else { // this beam pion event is not selected
        if (!isTestSample) {
          uf.response_SliceID_Inc.Miss(true_sliceID, weight); // missed event need bkgw
          uf.response_SliceID_Ini.Miss(true_ini_sliceID, weight);
        }
      }
      
      if ((*evt.true_beam_endProcess) == "pi+Inelastic"){ // true pion beam interaction event (exclude elastics)
        if (isTestSample){
          h_truesliceid_pioninelastic_all->Fill(true_sliceID, g4rw);
        }
        else{
          //uf.eff_den_Int->Fill(true_sliceID);
        }
        if (hadana.PassPiCuts(evt) && evt.reco_beam_true_byE_matched){
          if (isTestSample){
            h_recosliceid_pioninelastic_cuts->Fill(reco_sliceID, g4rw);
            h_truesliceid_pioninelastic_cuts->Fill(true_sliceID, g4rw);
          }
          else{
            //uf.eff_num_Int->Fill(true_sliceID);
            //uf.pur_num_Int->Fill(reco_sliceID);
            uf.response_SliceID_Int.Fill(reco_sliceID, true_sliceID, g4rw);
          }
        }
        else{
          if (!isTestSample) uf.response_SliceID_Int.Miss(true_sliceID, weight);
        }
      }
    }
    if (hadana.PassPiCuts(evt)){ // the event passed full selections
      if (isTestSample){
        //h_recosliceid_allevts_cuts->Fill(reco_sliceID, weight);
        //h_recoinisliceid_allevts_cuts->Fill(reco_ini_sliceID, weight);
      }
      else {
        //uf.pur_den->Fill(reco_sliceID);
      }
    }
  }
}

void ThinSlice::FillHistograms(int cut, const anavar & evt, double weight){
  if (hadana.fAllTrackCheck) {} // removed for not in use
  else{
    if (cut>=0 && cut < pi::nCuts){
      FillHistVec1D(hreco_beam_type[cut], evt.reco_beam_type, hadana.pitype, weight);
      FillHistVec1D(htrklen_csda_proton[cut], hadana.trklen_csda_proton, hadana.pitype, weight);
      FillHistVec1D(hChi2_proton[cut], hadana.chi2_proton, hadana.pitype, weight);
      FillHistVec1D(hreco_reconstructable_beam_event[cut], evt.reco_reconstructable_beam_event, hadana.pitype, weight);
      FillHistVec1D(h_diff_reco_true_Eint[cut], int_energy_reco - int_energy_true, hadana.pitype);
      FillHistVec2D(h_diff_reco_true_vs_true_Eint[cut], int_energy_true, int_energy_reco - int_energy_true, hadana.pitype);
      FillProfVec(pf_diff_reco_true_vs_true_Eint[cut], int_energy_true, int_energy_reco - int_energy_true, hadana.pitype);
      
      FillHistVec1D(htrue_beam_endZ[cut], evt.true_beam_endZ_SCE, hadana.pitype, weight);
      FillHistVec1D(htrue_beam_endZ_SCE[cut], evt.true_beam_endZ, hadana.pitype, weight); // it seems SCE is reversed? and I didn't find true_beam_endZ_SCE on wiki?
      FillHistVec1D(htrue_sliceID[cut], true_sliceID, hadana.pitype, weight);
      //    if (!evt.reco_beam_calo_wire->empty()){
      FillHistVec1D(hreco_beam_endZ[cut], evt.reco_beam_endZ, hadana.pitype, weight);
      FillHistVec1D(hreco_true_beam_endZ[cut], evt.reco_beam_endZ - evt.true_beam_endZ_SCE, hadana.pitype, weight);
      FillHistVec2D(hreco_vs_true_beam_endZ[cut], evt.true_beam_endZ_SCE, evt.reco_beam_endZ, hadana.pitype, weight);
      FillHistVec2D(hreco_true_vs_true_beam_endZ[cut], evt.true_beam_endZ_SCE, evt.reco_beam_endZ - evt.true_beam_endZ_SCE, hadana.pitype, weight);
      
      FillHistVec1D(hreco_beam_endZ_SCE[cut], evt.reco_beam_calo_endZ, hadana.pitype, weight);
      FillHistVec1D(hreco_true_beam_endZ_SCE[cut], evt.reco_beam_calo_endZ - evt.true_beam_endZ, hadana.pitype, weight);
      FillHistVec2D(hreco_vs_true_beam_endZ_SCE[cut], evt.true_beam_endZ, evt.reco_beam_calo_endZ, hadana.pitype, weight);
      FillHistVec2D(hreco_true_vs_true_beam_endZ_SCE[cut], evt.true_beam_endZ, evt.reco_beam_calo_endZ - evt.true_beam_endZ, hadana.pitype, weight);
      
      FillHistVec1D(hreco_sliceID[cut], reco_sliceID, hadana.pitype, weight);
      FillHistVec1D(hreco_incsliceID[cut], reco_end_sliceID, hadana.pitype, weight);
      FillHistVec1D(hreco_inisliceID[cut], reco_ini_sliceID, hadana.pitype, weight);
      FillHistVec1D(hreco_true_sliceID[cut], reco_sliceID - true_sliceID, hadana.pitype, weight);
      FillHistVec2D(hreco_vs_true_sliceID[cut], true_sliceID, reco_sliceID, hadana.pitype, weight);
      FillHistVec2D(hreco_true_vs_true_sliceID[cut], true_sliceID, reco_sliceID - true_sliceID, hadana.pitype, weight);
      
      // below are variables not provided by evt directly (calculated in hadana)
      FillHistVec1D(hmediandEdx[cut], hadana.median_dEdx, hadana.pitype, weight);
      FillHistVec1D(hdaughter_michel_score[cut], hadana.daughter_michel_score, hadana.pitype, weight);
      FillHistVec1D(henergy_calorimetry_SCE[cut], hadana.energy_calorimetry_SCE, hadana.pitype, weight);
      FillHistVec1D(hdEdx_SCE[cut], hadana.energy_calorimetry_SCE/hadana.reco_trklen, hadana.pitype, weight);
      if (evt.reco_beam_calo_endZ>300 && hadana.median_dEdx<2.4){ // likely to be a cosmic muon?
        if (hadana.daughter_michel_score>=0){
          FillHistVec1D(hdaughter_michel_scoreMu[cut], hadana.daughter_michel_score, hadana.pitype, weight);
          //if (!evt.MC) cout<<evt.run<<" "<<evt.subrun<<" "<<evt.event<<" "<<hadana.daughter_michel_score<<" "<<evt.reco_beam_calo_wire->back()<<" "<<evt.reco_beam_calo_tick->back()<<" "<<evt.reco_beam_calo_wire->front()<<" "<<evt.reco_beam_calo_tick->front()<<endl;
        }
        int nhits = 0;
        double michelscore = 0;
        for (size_t i = 0; i<evt.reco_daughter_PFP_michelScore_collection->size(); ++i){
          nhits += (*evt.reco_daughter_PFP_nHits_collection)[i];
          michelscore += (*evt.reco_daughter_PFP_michelScore_collection)[0] * (*evt.reco_daughter_PFP_nHits_collection)[i];
        }
        if (nhits && michelscore>=0){
          michelscore/=nhits;
          FillHistVec1D(hdaughter_michel_score2Mu[cut], michelscore, hadana.pitype, weight); // what's PFP and what's difference between hdaughter_michel_scoreMu and hdaughter_michel_score2Mu?
        }
  //      if (hadana.pitype == kMuon && hadana.daughter_michel_score < 0.01){
  //        cout<<evt.run<<" "<<evt.subrun<<" "<<evt.event<<endl;
  //      }
      }
      if (evt.reco_beam_calo_endZ<100 && hadana.median_dEdx<2.4){
        if (hadana.daughter_michel_score>=0){
          //if (!evt.MC) cout<<evt.run<<" "<<evt.subrun<<" "<<evt.event<<" "<<hadana.daughter_michel_score<<" "<<evt.reco_beam_calo_wire->back()<<" "<<evt.reco_beam_calo_tick->back()<<" "<<evt.reco_beam_calo_wire->front()<<" "<<evt.reco_beam_calo_tick->front()<<endl;
          FillHistVec1D(hdaughter_michel_scorePi[cut], hadana.daughter_michel_score, hadana.pitype, weight);
        }
      }
      /*if (reco_sliceID>=0 && reco_sliceID<pi::nthinslices){
        FillHistVec1D(hmediandEdxSlice[reco_sliceID][cut], hadana.median_dEdx, hadana.pitype);
        FillHistVec1D(hChi2_protonSlice[reco_sliceID][cut], hadana.chi2_proton, hadana.pitype);
        FillHistVec1D(hdaughter_michel_scoreSlice[reco_sliceID][cut], hadana.daughter_michel_score, hadana.pitype);
        FillHistVec1D(hcosthetaSlice[reco_sliceID][cut], hadana.beam_costh, hadana.pitype);
      }*/

      FillHistVec1D(htrackscore[cut], evt.reco_beam_PFP_trackScore_collection, hadana.pitype, weight);
      FillHistVec1D(hemscore[cut], evt.reco_beam_PFP_emScore_collection, hadana.pitype, weight);
  //    if (cut == kAPA3 && evt.reco_beam_PFP_emScore_collection > 0.9){
  //      cout<<evt.run<<" "<<evt.event<<" "<<evt.reco_beam_PFP_emScore_collection<<" "<<evt.reco_beam_calo_wire->front()<<" "<<evt.reco_beam_calo_tick->back()<<endl;
  //    }
      FillHistVec1D(hdEdx_5cm[cut], hadana.dEdx_5cm, hadana.pitype, weight);

      FillHistVec1D(hdeltax[cut], hadana.beam_dx, hadana.pitype, weight);
      FillHistVec1D(hdeltay[cut], hadana.beam_dy, hadana.pitype, weight);
      FillHistVec1D(hdeltaz[cut], hadana.beam_dz, hadana.pitype, weight);
      FillHistVec1D(hcostheta[cut], hadana.beam_costh, hadana.pitype, weight);
      FillHistVec1D(hdeltax_inst[cut], evt.beam_inst_X, hadana.pitype, weight);
      FillHistVec1D(hdeltay_inst[cut], evt.beam_inst_Y, hadana.pitype, weight);
      FillHistVec2D(hxy_inst[cut], evt.beam_inst_X, evt.beam_inst_Y, hadana.pitype, weight);

      FillHistVec1D(hreco_beam_true_byE_matched[cut], evt.reco_beam_true_byE_matched, hadana.pitype, weight);
      FillHistVec1D(hreco_trklen[cut], hadana.reco_trklen, hadana.pitype, weight);
      FillHistVec1D(htrue_trklen[cut], hadana.true_trklen, hadana.pitype, weight);
      FillHistVec1D(hdiff_trklen[cut], hadana.reco_trklen - hadana.true_trklen, hadana.pitype, weight);
      FillHistVec2D(hreco_vs_true_trklen[cut], hadana.true_trklen, hadana.reco_trklen, hadana.pitype, weight);
      FillHistVec1D(hbeam_score[cut], hadana.beam_score, hadana.pitype, weight);
      FillHistVec2D(beam_score_vs_hreco_trklen[cut], hadana.reco_trklen, hadana.beam_score, hadana.pitype, weight);
      
      //$$$temp
     /* if ( hadana.true_trklen>20 && evt.reco_beam_alt_len>20){
        int printout = kFALSE;
        if ( hadana.true_trklen>250 && abs(evt.reco_beam_alt_len-230)<5 ){
          cout<<"$$$$$ red bar ("<<cut<<" "<<pi::cutName[cut]<<",\t"<<pi::intTypeName[hadana.pitype]<<")\n";
          printout = kTRUE;
        }
        if ( hadana.true_trklen<200 && abs(evt.reco_beam_alt_len-230)<5 ){
          cout<<"$$$$$ blue bar ("<<cut<<" "<<pi::cutName[cut]<<",\t"<<pi::intTypeName[hadana.pitype]<<")\n";
          printout = kTRUE;
        }
        if ( abs(hadana.true_trklen-evt.reco_beam_alt_len-230)<5 ){
          cout<<"$$$$$ green bar ("<<cut<<" "<<pi::cutName[cut]<<",\t"<<pi::intTypeName[hadana.pitype]<<")\n";
          printout = kTRUE;
        }
        if (printout == kTRUE){
          cout<<"Run: "<<evt.run<<";\t"
          <<"SubRun: "<<evt.subrun<<";\t"
          <<"Event: "<<evt.event<<endl;
          cout<<"True trklen: "<<hadana.true_trklen<<";\t"
          <<"Reco trklen: "<<evt.reco_beam_alt_len<<endl;
          cout<<"Start: ("<<evt.reco_beam_calo_startX<<", "<<evt.reco_beam_calo_startY<<", "<<evt.reco_beam_calo_startZ<<");\n";
          cout<<"End: ("<<evt.reco_beam_calo_endX<<", "<<evt.reco_beam_calo_endY<<", "<<evt.reco_beam_calo_endZ<<")\n";
        }
      }*/
      /*if (!evt.MC){
        if (hadana.beam_dy>4 && hadana.beam_dy<6){
          cout<<"$$$$$ delta Y [4,6] ("<<cut<<" "<<pi::cutName[cut]<<",\t"<<pi::intTypeName[hadana.pitype]<<")\n";
          cout<<"Run: "<<evt.run<<";\t"
          <<"SubRun: "<<evt.subrun<<";\t"
          <<"Event: "<<evt.event<<endl;
          cout<<"Start: ("<<evt.reco_beam_calo_startX<<", "<<evt.reco_beam_calo_startY<<", "<<evt.reco_beam_calo_startZ<<");\t trklen: "<<hadana.reco_trklen<<endl;
          cout<<"End: ("<<evt.reco_beam_calo_endX<<", "<<evt.reco_beam_calo_endY<<", "<<evt.reco_beam_calo_endZ<<")\n";
          
        }
      }*/

      FillHistVec1D(hreco_beam_startX_SCE[cut], evt.reco_beam_calo_startX, hadana.pitype, weight);
      FillHistVec1D(hreco_beam_startY_SCE[cut], evt.reco_beam_calo_startY, hadana.pitype, weight);
      FillHistVec1D(hreco_beam_startZ_SCE[cut], evt.reco_beam_calo_startZ, hadana.pitype, weight);

      if (!evt.reco_beam_calo_wire->empty()){
        TVector3 pt0(evt.reco_beam_calo_startX,
                     evt.reco_beam_calo_startY,
                     evt.reco_beam_calo_startZ);
        TVector3 pt1(evt.reco_beam_calo_endX,
                     evt.reco_beam_calo_endY,
                     evt.reco_beam_calo_endZ);
        TVector3 dir = pt1 - pt0;
        dir = dir.Unit();
        FillHistVec1D(hreco_beam_dcosX_SCE[cut], dir.X(), hadana.pitype, weight);
        FillHistVec1D(hreco_beam_dcosY_SCE[cut], dir.Y(), hadana.pitype, weight);
        FillHistVec1D(hreco_beam_dcosZ_SCE[cut], dir.Z(), hadana.pitype, weight);
        FillHistVec1D(hreco_beam_angleX_SCE[cut], acos(dir.X())*180/TMath::Pi(), hadana.pitype, weight);
        FillHistVec1D(hreco_beam_angleY_SCE[cut], acos(dir.Y())*180/TMath::Pi(), hadana.pitype, weight);
        FillHistVec1D(hreco_beam_angleZ_SCE[cut], acos(dir.Z())*180/TMath::Pi(), hadana.pitype, weight);
      }

      FillHistVec2D(hreco_beam_startXY_SCE[cut], evt.reco_beam_calo_startX, evt.reco_beam_calo_startY, hadana.pitype, weight);

    }
  }
}

void ThinSlice::FillSliceHist(const anavar & evt, int constraint_type, double weight, int cut){
  if (hadana.fAllTrackCheck) {
    cout<<"Warning: AllTrackCheck hasn't been fully implemented!!!"<<endl; // Do we need AllTrackCheck in main branch?
  }
  else {
    if (constraint_type == 1) { // muon
      FillHistVec1D(hdaughter_michel_score_bkg[cut], hadana.daughter_michel_score, hadana.pitype, weight, false, false);
    }
    else if (constraint_type == 2) { // proton
      FillHistVec1D(hmediandEdx_bkg[cut], hadana.median_dEdx, hadana.pitype, weight, false, true);
      FillHistVec1D(hChi2_proton_bkg[cut], hadana.chi2_proton, hadana.pitype, weight, false, true);
    }
    else if (constraint_type == 3) { // secondary pion
      FillHistVec1D(hcostheta_bkg[cut], hadana.beam_costh, hadana.pitype, weight, true, false);
    }
    // in each slice
    if (reco_sliceID>=0 && reco_sliceID<pi::nthinslices){
      if (constraint_type == 1) { // muon
        FillHistVec1D(hdaughter_michel_scoreSlice[reco_sliceID][cut], hadana.daughter_michel_score, hadana.pitype, weight, false, false);
      }
      else if (constraint_type == 2) { // proton
        FillHistVec1D(hmediandEdxSlice[reco_sliceID][cut], hadana.median_dEdx, hadana.pitype, weight, false, true);
        FillHistVec1D(hChi2_protonSlice[reco_sliceID][cut], hadana.chi2_proton, hadana.pitype, weight, false, true);
      }
      else if (constraint_type == 3) { // secondary pion
        FillHistVec1D(hcosthetaSlice[reco_sliceID][cut], hadana.beam_costh, hadana.pitype, weight, true, false);
      }
    }
  }
}

void ThinSlice::SaveHistograms(){
  outputFile->cd();
  outputFile->Write();
  h_truesliceid_pion_uf->Write("h_truesliceid_pion_uf");
  h_truesliceid_pioninelastic_uf->Write("h_truesliceid_pioninelastic_uf");
  //h_trueinisliceid_pion_uf->Write("h_trueinisliceid_pion_uf");
  //response_SliceID_Pion->Write("response_SliceID_Pion");
  //response_SliceID_PionInEl->Write("response_SliceID_PionInEl");
}

void ThinSlice::CalcXS(const Unfold & uf){

  double slcid[pi::nthinslices] = {0};
  double Eslice[pi::nthinslices] = {0};
  double Einterval[pi::nthinslices] = {0};
  double dEdx[pi::nthinslices] = {0};
  double avg_trueincE[pi::nthinslices] = {0};
  double avg_recoincE[pi::nthinslices] = {0};
  double err_trueincE[pi::nthinslices] = {0};
  double err_recoincE[pi::nthinslices] = {0};
  double reco_trueincE[pi::nthinslices] = {0};
  double err_reco_trueincE[pi::nthinslices] = {0};
  double truexs[pi::nthinslices] = {0};
  double err_truexs[pi::nthinslices] = {0};
  double true_cosangle = 1.;

  double NA=6.02214076e23;
  double MAr=39.95; //gmol
  double Density = 1.4; // 1.396 g/cm^3
  for (int i = 0; i<pi::nthinslices; ++i){
    
    slcid[i] = i;
    Eslice[i] = pi::plim - (i+0.5)*pi::Eslicewidth;
    Einterval[i] = pi::Eslicewidth/2;
    dEdx[i] = bb.meandEdx(Eslice[i]); // MeV/cm
    //cout<<dEdx[i]<<"\t";
    avg_trueincE[i] = true_incE[i]->GetMean();
    err_trueincE[i] = true_incE[i]->GetMeanError();
    avg_recoincE[i] = reco_incE[i]->GetMean();
    err_recoincE[i] = reco_incE[i]->GetMeanError();
    reco_trueincE[i] = avg_recoincE[i] - avg_trueincE[i];
    err_reco_trueincE[i] = sqrt(pow(err_trueincE[i],2)+pow(err_recoincE[i],2)); // is it proper to simply use root_sum_square, since the two seem not independent?
    //std::cout<<i<<" "<<avg_trueincE[i]<<std::endl;
    if (true_incidents[i] && true_interactions[i]){
      //true_cosangle = true_AngCorr->GetMean(); // no need to include angle correction
      truexs[i] = dEdx[i]*MAr/(Density*NA*pi::Eslicewidth/true_cosangle)*log(true_incidents[i]/(true_incidents[i]-true_interactions[i]))*1e27;
      err_truexs[i] = dEdx[i]*MAr/(Density*NA*pi::Eslicewidth/true_cosangle)*1e27*sqrt(true_interactions[i]+pow(true_interactions[i],2)/true_incidents[i])/true_incidents[i];
    }
  }

  TGraphErrors *gr_trueincE = new TGraphErrors(pi::nthinslices, &(slcid[0]), &(avg_trueincE[0]), 0, &(err_trueincE[0]));
  TGraphErrors *gr_recoincE = new TGraphErrors(pi::nthinslices, &(slcid[0]), &(avg_recoincE[0]), 0, &(err_recoincE[0]));
  TGraphErrors *gr_reco_trueincE = new TGraphErrors(pi::nthinslices, &(slcid[0]), &(reco_trueincE[0]), 0, &(err_reco_trueincE[0]));

  gr_trueincE->Write("gr_trueincE");
  gr_recoincE->Write("gr_recoincE");
  gr_reco_trueincE->Write("gr_reco_trueincE");

  TGraphErrors *gr_truexs = new TGraphErrors(pi::nthinslices, &(Eslice[0]), &(truexs[0]), &(Einterval[0]), &(err_truexs[0]));
  
  gr_truexs->Write("gr_truexs");

  TH1D *hinc = (TH1D*)h_recosliceid_allevts_cuts->Clone("hinc");
  TH1D *hint = (TH1D*)h_recosliceid_allevts_cuts->Clone("hint");
  hinc->Multiply(uf.pur_Inc);
  hint->Multiply(uf.pur_Int);

//  RooUnfoldBayes   unfold_Inc (&uf.response_SliceID_Inc, uf.pur_num_Inc, 4);
//  RooUnfoldBayes   unfold_Int (&uf.response_SliceID_Int, uf.pur_num_Int, 4);

  RooUnfoldBayes   unfold_Inc (&uf.response_SliceID_Inc, hinc, 4);
  RooUnfoldBayes   unfold_Int (&uf.response_SliceID_Int, hint, 4);
  //RooUnfoldBayes   unfold_Ini (&uf.response_SliceID_Ini, hini, 4);

//  RooUnfoldSvd     unfold_Inc (&uf.response_SliceID_Inc, hinc, 20);   // OR
//  RooUnfoldSvd     unfold_Int (&uf.response_SliceID_Int, hint, 20);   // OR

//  RooUnfoldSvd     unfold_Inc (&uf.response_SliceID_Inc, uf.pur_num_Inc, 20);   // OR
//  RooUnfoldSvd     unfold_Int (&uf.response_SliceID_Int, uf.pur_num_Int, 20);   // OR

  h_truesliceid_pion_uf = (TH1D*) unfold_Inc.Hreco();
  h_truesliceid_pioninelastic_uf = (TH1D*) unfold_Int.Hreco();
  //h_trueinisliceid_pion_uf = (TH1D*) unfold_Ini.Hreco();
}

void ThinSlice::Run(anavar & evt, Unfold & uf, Long64_t nentries, bool random){

  BookHistograms();

  //Long64_t nentries = evt.fChain->GetEntriesFast();
  if (nentries == -1) nentries = evt.fChain->GetEntries();
  
  Long64_t nbytes = 0, nb = 0;
  TRandom3 *r3 = new TRandom3(0);
  for (Long64_t num=0; num<nentries; num++) {
    if (num%10000==0) std::cout<<num<<"/"<<nentries<<std::endl;
    Long64_t jentry = num;
    if (random) jentry = TMath::FloorNint(r3->Rndm()*302141);
    //cout<<jentry<<endl;
    Long64_t ientry = evt.LoadTree(jentry);
    if (ientry < 0) break;
    nb = evt.fChain->GetEntry(jentry);   nbytes += nb;
    // if (Cut(ientry) < 0) continue;
    //std::cout<<evt.run<<" "<<evt.event<<" "<<evt.MC<<" "<<evt.reco_beam_true_byE_matched<<" "<<evt.true_beam_PDG<<" "<<(*evt.true_beam_endProcess)<<std::endl;
    //std::cout<<GetParType(ana)<<std::endl;
    if (selectCosmics){
      if (!hadana.isCosmics(evt)) continue;
    }
    else{
      if (!hadana.isSelectedPart(evt)) continue;
    }
    
    hadana.ProcessEvent(evt); // hadana.pitype is assigned in this step
    
    double weight = CalWeight(evt, hadana.pitype); // muon reweight; momentum reweight (to reconcile real data and MC)
    double g4rw = 1; // geant4reweight (to fake data for test; it turns out it should also be applied to true MC to make unfolding reliable)
    double bkgw = 1; // bkg fraction variation (to fake data for test)
    
    if (evt.MC) {
      double weiarr_fd[20] = {
        1,1,1,1,1,1,1,1,1,1,
        1,1,1,1,1,1,1,1,1,1,
      };
      double weiarr_mc[20] = {
        1,1,1,1,1,1,1,1,1,1,
        1,1,1,1,1,1,1,1,1,1,
      };
      /*double weiarr_fd[20] = {
        1.80, 1.80, 1.80, 1.80, 1.70,
        1.60, 1.50, 1.40, 1.30, 1.20,
        1.10, 1.00, 0.90, 0.84, 0.78,
        0.80, 0.83, 0.85, 0.86, 0.90,
      };
      double weiarr_mc[20] = {
        1.80, 1.80, 1.80, 1.80, 1.70,
        1.60, 1.50, 1.40, 1.30, 1.20,
        1.10, 1.00, 0.90, 0.84, 0.78,
        0.80, 0.83, 0.85, 0.86, 0.90,
      };*/
      if (hadana.pitype == 0) { // fake data
        g4rw = CalG4RW(evt, weiarr_fd);
        bkgw = CalBkgW(evt, 1., 1., 1.);
      }
      else { // true MC
        g4rw = CalG4RW(evt, weiarr_mc);
        bkgw = CalBkgW(evt);
      }
    }
    ProcessEvent(evt, uf, g4rw, bkgw);
    weight *= g4rw;
    weight *= bkgw;
    // can change order of cuts
    FillHistograms(pi::kNocut, evt, weight);
    if (hadana.PassPandoraSliceCut(evt)){
      FillHistograms(pi::kPandoraSlice, evt, weight);
      if (hadana.PassCaloSizeCut(evt)){
        FillHistograms(pi::kCaloSize, evt, weight);
        if (hadana.PassBeamQualityCut(evt)){
          FillHistograms(pi::kBeamQuality, evt, weight);
          if (hadana.PassProtonCut()){
            FillHistograms(pi::kProtonCut, evt, weight);
            if (hadana.PassMichelScoreCut()){
              FillHistograms(pi::kMichelScore, evt, weight);
              if (hadana.PassAPA3Cut(evt)){
                FillHistograms(pi::kAPA3, evt, weight);
              }
            }
          }
        }
        // for background constraints
        if (hadana.PassAPA3Cut(evt)){
          if (hadana.PassBeamQualityCut(evt) && hadana.PassProtonCut()) { // to constrain muon
            FillSliceHist(evt, 1, weight);
          }
          if (hadana.PassBeamQualityCut(evt) && hadana.PassMichelScoreCut()) { // to constrain proton
            FillSliceHist(evt, 2, weight);
          }
          if (hadana.PassBeamQualityCut(evt, false) && hadana.PassProtonCut() && hadana.PassMichelScoreCut()) { // to constrain secondary pion
            FillSliceHist(evt, 3, weight);
          }
        }
      }
    }
    FillSliceHist(evt, 1, weight, 0);
    FillSliceHist(evt, 2, weight, 0);
    FillSliceHist(evt, 3, weight, 0);
  }
  
  uf.SaveHistograms();
  CalcXS(uf);
  SaveHistograms();
}

void ThinSlice::SetSelectCosmics(bool sc){

  selectCosmics = sc;

}
