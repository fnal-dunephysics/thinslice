#include "anavar.h"
#include "HadAna.h"
#include "TMath.h"
#include "TVector3.h"
#include <iostream>

HadAna::HadAna(){
  map_BB[13] = new BetheBloch(13);
  map_BB[211] = new BetheBloch(211);
  map_BB[321] = new BetheBloch(321);
  map_BB[2212] = new BetheBloch(2212);

  if (fProtonCSDACheck) {
    TFile *file_mom2csda = TFile::Open("root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/dune/archive/sam_managed_users/tjyang/data/3/a/7/a/51075a8b-cc10-4d2d-b4c6-8e4c4d2817d7-proton_mom_csda_converter.root");
    csda_range_vs_mom_sm = (TGraph *)file_mom2csda->Get("csda_range_vs_mom_sm");
  }
}

void HadAna::InitPi(){

  AddTruePDG(-13);
  AddTruePDG(13);
  AddTruePDG(211);

  SetPandoraSlicePDG(13);

  SetBeamQualityCuts();
}

void HadAna::InitP(){

  AddTruePDG(2212);

  SetPandoraSlicePDG(13);

  SetBeamQualityCuts();
}


void HadAna::AddTruePDG(int pdg){
  truepdglist.push_back(pdg);
};

bool HadAna::isSelectedPart(const anavar& evt) const{
  if (evt.reco_reconstructable_beam_event == 0) return false; // remove empty events first
  if (evt.MC){
    for (size_t i = 0; i<truepdglist.size(); ++i){
      if (evt.true_beam_PDG == truepdglist[i]) return true; // truth matched
    }
    return false;
  }
  else{ // real data
    if (evt.beam_inst_trigger == 8) return false; // is cosmics
    if (evt.beam_inst_nMomenta != 1 || evt.beam_inst_nTracks != 1) return false;
    for (size_t i = 0; i<truepdglist.size(); ++i){
      for (size_t j = 0; j<evt.beam_inst_PDG_candidates->size(); ++j){
        if ((*evt.beam_inst_PDG_candidates)[j] == truepdglist[i]) return true; // what's this used for?
      }
    }
    return false;
  }
}

bool HadAna::isCosmics(const anavar& evt) const{
  if (evt.MC) return false;
  else{
    if (evt.beam_inst_trigger == 8) return true;
    else return false;
  }
  return false;
}

void HadAna::SetPandoraSlicePDG(int pdg){
  pandora_slice_pdg = pdg;
};

void HadAna::SetBeamQualityCuts(double dx_min, double dx_max,
                                double dy_min, double dy_max,
                                double dz_min, double dz_max,
                                double dxy_min, double dxy_max,
                                double costh_min, double costh_max){
  beamcut_dx_min = dx_min; beamcut_dx_max = dx_max;
  beamcut_dy_min = dy_min; beamcut_dy_max = dy_max;
  beamcut_dz_min = dz_min; beamcut_dz_max = dz_max;
  beamcut_dxy_min = dxy_min; beamcut_dxy_max = dxy_max;
  beamcut_costh_min = costh_min; beamcut_costh_max = costh_max;
}

int HadAna::GetPiParType(const anavar& evt){

  if (!evt.MC){
    return pi::kData;
  }
  else if (evt.event%2){ // divide half of MC as fake data
    return pi::kData;
  }
  else if (!evt.reco_beam_true_byE_matched){ // the true beam track is not selected
    if (evt.reco_beam_true_byE_origin == 2) {
      return pi::kMIDcosmic;
    }
    else if (std::abs(evt.reco_beam_true_byE_PDG) == 211){ // the selected track is a pion (but not true beam pion, so it is a secondary pion)
      return pi::kMIDpi;
    }
    else if (evt.reco_beam_true_byE_PDG == 2212){
      return pi::kMIDp;
    }
    else if (std::abs(evt.reco_beam_true_byE_PDG) == 13){
      return pi::kMIDmu;
    }
    else if (std::abs(evt.reco_beam_true_byE_PDG) == 11 ||
             evt.reco_beam_true_byE_PDG == 22){
      return pi::kMIDeg;
    }
    else {
      //cout<<evt.reco_beam_true_byE_PDG<<endl;
      return pi::kMIDother;
    }
  }
  else if (evt.true_beam_PDG == -13){
    return pi::kMuon;
  }
  else if (evt.true_beam_PDG == 211){
    if ((*evt.true_beam_endProcess) == "pi+Inelastic"){
      return pi::kPiInel;
    }
    else return pi::kPiElas;
  }
  
  return pi::kMIDother;
}

int HadAna::GetPParType(const anavar& evt){

  if (!evt.MC){
    return p::kData;
  }
  else if (evt.event%2){
    return p::kData;
  }
  else if (!evt.reco_beam_true_byE_matched){
    if (evt.reco_beam_true_byE_origin == 2) {
      return p::kMIDcosmic;
    }
    else if (std::abs(evt.reco_beam_true_byE_PDG) == 211){
      return p::kMIDpi;
    }
    else if (evt.reco_beam_true_byE_PDG == 2212){
      return p::kMIDp;
    }
    else if (std::abs(evt.reco_beam_true_byE_PDG) == 13){
      return p::kMIDmu;
    }
    else if (std::abs(evt.reco_beam_true_byE_PDG) == 11 ||
             evt.reco_beam_true_byE_PDG == 22){
      return p::kMIDeg;
    }
    else {
      //std::cout<<reco_beam_true_byE_PDG<<std::endl;
      return p::kMIDother;
    }
  }
  else if (evt.true_beam_PDG == 2212){
    if ((*evt.true_beam_endProcess) == "protonInelastic"){
      return p::kPInel;
    }
    else return p::kPElas;
  }
  
  return p::kMIDother;
}

bool HadAna::PassPandoraSliceCut(const anavar& evt) const{ // whether recognized by Pandora correctly

  if (fAllTrackCheck) return evt.reco_beam_type > 0;
  else return (evt.reco_beam_type == pandora_slice_pdg);
}

bool HadAna::PassBeamQualityCut(bool has_angle_cut) const{ // cut on beam entrance location and beam angle

  if (beamcut_dx_min<beamcut_dx_max){
    if (beam_dx<beamcut_dx_min)
      return false;
    if (beam_dx>beamcut_dx_max)
      return false;
  }

  if (beamcut_dy_min<beamcut_dy_max){
    if (beam_dy<beamcut_dy_min)
      return false;
    if (beam_dy>beamcut_dy_max)
      return false;
  }

  if (beamcut_dz_min<beamcut_dz_max){
    if (beam_dz<beamcut_dz_min)
      return false;
    if (beam_dz>beamcut_dz_max)
      return false;
  }

  if (beamcut_dxy_min<beamcut_dxy_max){
    if (beam_dxy<beamcut_dxy_min)
      return false;
    if (beam_dxy>beamcut_dxy_max)
      return false;
  }

  if (has_angle_cut && beamcut_costh_min<beamcut_costh_max){
    if (beam_costh<beamcut_costh_min)
      return false;
    if (beam_costh>beamcut_costh_max)
      return false;
  }

  return true;
}

bool HadAna::PassBeamXYCut(const anavar& evt) const{

  //ref: https://indico.fnal.gov/event/55048/contributions/244291/attachments/156262/203827/20220616_hy.pdf
  if (!evt.MC){//data
    if ((pow(((evt.beam_inst_X-meanX_data)/(1.5*rmsX_data)),2)+pow(((evt.beam_inst_Y-meanY_data)/(1.5*rmsY_data)),2))<=1.) return true;
    else return false;
  }
  else{//mc
    if ((pow(((evt.beam_inst_X-meanX_mc)/(1.5*rmsX_mc)),2)+pow(((evt.beam_inst_Y-meanY_mc)/(1.5*rmsY_mc)),2))<=1.) return true;
    else return false;
  }
}

bool HadAna::PassAPA3Cut(const anavar& evt) const{ // only use track in the first TPC
  return true;
  double cutAPA3_Z = 220.;
  
  if (fAllTrackCheck) return evt.reco_beam_calo_endZ_allTrack < cutAPA3_Z;
  else return evt.reco_beam_calo_endZ < cutAPA3_Z;
}

bool HadAna::PassCaloSizeCut(const anavar& evt) const{ // Require hits information in collection plane
  
  if (fAllTrackCheck) return !(evt.reco_beam_calo_wire_allTrack->empty());
  else return !(evt.reco_beam_calo_wire->empty());
}

bool HadAna::PassMichelScoreCut() const{ // further veto muon tracks according to Michel score
  
  return daughter_michel_score < 0.55;
}

bool HadAna::PassProtonCut() const{ // to remove proton background

  return chi2_proton > 80; //median_dEdx < 2.4;
}

bool HadAna::PassPiCuts(const anavar& evt) const{
  return PassPandoraSliceCut(evt)&&
    PassCaloSizeCut(evt)&&
    PassBeamQualityCut()&&
    PassAPA3Cut(evt)&&
    PassMichelScoreCut()&&
    PassProtonCut();
}

bool HadAna::PassPCuts(const anavar& evt) const{
  return PassPandoraSliceCut(evt)&&
    PassCaloSizeCut(evt)&&
    PassBeamQualityCut()&&
    PassBeamXYCut(evt);
}

double HadAna::Get_true_ffKE(const anavar& evt, double KE_in_TPC, double length_to_ff){
  double this_dEdx = map_BB[abs(evt.true_beam_PDG)] -> meandEdx(KE_in_TPC);
  return KE_in_TPC + this_dEdx * length_to_ff;
}

double HadAna::Fit_dEdx_Residual_Length(const anavar& evt, const vector<double> & dEdx, const vector<double> & ResRange, int PID, bool save_graph){

  int N_max = 20; // == Maximum number of hits used for the Bethe-Bloch fitting

  // == PID input : mass hypothesis, valid only for muons, charged pions, and protons
  int abs_PID = abs(PID);
  if(!(abs(PID) == 13 || PID == 2212 || abs(PID) == 211)){
    //cout << "[HadAna::Fit_dEdx_Residual_Length] Not a valid PID!" << endl;
    return -9999.;
  }

  double best_additional_res_length = -0.1;
  double best_chi2 = 99999.;
  double min_additional_res_length = 0.; // == [cm]
  double max_additional_res_length = 60.; // == [cm]
  double res_length_step = 0.1; // == [cm]
  int res_length_trial = (max_additional_res_length - min_additional_res_length) / res_length_step;
  int this_N_calo = dEdx.size();
  if(this_N_calo <= 15){
    //cout << "[HadAna::Fit_dEdx_Residual_Length] Too small number of hits!" << endl;
    return -9999.; // == Too small number of hits
  }
  int this_N_hits = TMath::Min(this_N_calo, N_max); // == Use how many hits
  int i_bestfit = -1;
  vector<double> chi2_vector;
  vector<double> additional_res_legnth_vector;
  for(int i = 0; i < res_length_trial; i++){
    double this_additional_res_length = min_additional_res_length + (i + 0.) * res_length_step;
    double this_chi2 = 0.;
    for(int j = 5; j < this_N_hits - 5; j++){ // == Do not use first and last 5 hits
      int this_index = this_N_calo - 1 - j;
      double this_res_length = ResRange.at(this_index) - ResRange.at(this_N_calo - this_N_hits) + this_additional_res_length;
      //double this_KE = ResLength_to_KE_BB(this_res_length, this_mass);
      double this_KE = map_BB[abs_PID]->KEFromRangeSpline(this_res_length);
      //double dEdx_theory = dEdx_Bethe_Bloch(this_KE, this_mass);
      double dEdx_theory = map_BB[abs_PID]->meandEdx(this_KE);
      double dEdx_measured = dEdx.at(this_index);
      if(dEdx_measured < 0.5 || dEdx_measured > 5.0) continue; // == Truncate, it should be modified to consider protons

      // == Gaussian approx.
      //double dEdx_theory_err = dEdx_theory * 0.02;
      this_chi2 += pow(dEdx_measured - dEdx_theory, 2);
    }
    this_chi2 = this_chi2 / (this_N_hits + 0.); // == chi2 / n.d.f
    if(this_chi2 < best_chi2){
      best_chi2 = this_chi2;
      best_additional_res_length = this_additional_res_length;
      i_bestfit = i;
    }
    if(save_graph){
      // == Save vectors for graphes
      chi2_vector.push_back(this_chi2);
      additional_res_legnth_vector.push_back(this_additional_res_length);
    }
  }

  if(save_graph){
    // == Vectors for graphes
    vector<double> range_original;
    vector<double> range_bestfit;
    vector<double> range_reco;
    vector<double> dEdx_ordered;
    for(int i = 5; i < this_N_hits - 5; i++){
      int this_index = this_N_calo - 1 - i;
      range_original.push_back(ResRange.at(this_index) - ResRange.at(this_N_calo - this_N_hits));
      range_bestfit.push_back(ResRange.at(this_index) - ResRange.at(this_N_calo - this_N_hits) + best_additional_res_length);
      range_reco.push_back(ResRange.at(this_index));
      dEdx_ordered.push_back(dEdx.at(this_index));
    }
    TGraph *dEdx_gr = new TGraph(this_N_hits - 10, &range_original[0], &dEdx_ordered[0]);
    dEdx_gr -> SetName(Form("dEdx_Run%d_Evt%d_Nhit%d", evt.run, evt.event, this_N_hits));
    dEdx_gr -> Write();
    delete dEdx_gr;

    TGraph *dEdx_bestfit_gr = new TGraph(this_N_hits - 10,&range_bestfit[0], &dEdx_ordered[0]);
    dEdx_bestfit_gr -> SetName(Form("dEdx_bestfit_Run%d_Evt%d_Nhit%d", evt.run, evt.event, this_N_hits));
    dEdx_bestfit_gr -> Write();
    delete dEdx_bestfit_gr;

    TGraph *dEdx_reco_gr = new TGraph(this_N_hits - 10,&range_reco[0], &dEdx_ordered[0]);
    dEdx_reco_gr -> SetName(Form("dEdx_reco_Run%d_Evt%d_Nhit%d", evt.run, evt.event, this_N_hits));
    dEdx_reco_gr -> Write();
    delete dEdx_reco_gr;

    TGraph *chi2_gr = new TGraph(additional_res_legnth_vector.size(), &additional_res_legnth_vector[0], &chi2_vector[0]);
    chi2_gr -> SetName(Form("Chi2_Run%d_Evt%d_Nhit%d", evt.run, evt.event, this_N_hits));
    chi2_gr -> Write();
    chi2_vector.clear();
    additional_res_legnth_vector.clear();
    delete chi2_gr;
  }

  double original_res_length = ResRange.at(this_N_calo - 1) - ResRange.at(this_N_calo - this_N_hits); // == [cm]
  double best_total_res_length = best_additional_res_length + original_res_length;
  double best_KE = map_BB[abs_PID]->KEFromRangeSpline(best_total_res_length);
  double best_mom = map_BB[abs_PID]->KEtoMomentum(best_KE);

  // == Define fitting failed cases
  if(i_bestfit == res_length_trial - 1){
    //cout << "[HadAna::Fit_Beam_Hit_dEdx_Residual_Length] Fit failed : no mimumum" << endl;
    return -9999.;
  }
  else if(best_chi2 > 99990.){
    //cout << "[HadAna::Fit_Beam_Hit_dEdx_Residual_Length] Fit failed : best_chi2 > 99990." << endl;
    return -9999.;
  }
  else if(best_chi2 < 1.0e-11){
    //cout << "[HadAna::Fit_Beam_Hit_dEdx_Residual_Length] Fit failed : best_chi2 < 1.0e-11" << endl;
    return -9999.;
  }

  return best_mom;
}

void HadAna::ProcessEvent(const anavar& evt){

  pitype = GetPiParType(evt);
  ptype = GetPParType(evt);

  median_dEdx = -1;
  chi2_proton = -1;
  daughter_michel_score = -999;

  if (fAllTrackCheck) {} // removed for not in use
  else {
    if (!evt.reco_beam_calo_wire->empty()){
      chi2_proton = evt.reco_beam_Chi2_proton/evt.reco_beam_Chi2_ndof;
      median_dEdx = TMath::Median(evt.reco_beam_calibrated_dEdX_SCE->size(), &(*evt.reco_beam_calibrated_dEdX_SCE)[0]);//TMath::Median(evt.reco_beam_calibrated_dEdX_SCE->size(), &(*evt.reco_beam_calibrated_dEdX_SCE)[0]);
  //    daughter_michel_score = 0;
  //    int nhits = 0;
  //    for (size_t i = 0; i<reco_daughter_PFP_michelScore_collection->size(); ++i){
  //      nhits += (*reco_daughter_PFP_nHits_collection)[i];
  //      daughter_michel_score += (*reco_daughter_PFP_michelScore_collection)[0] * (*reco_daughter_PFP_nHits_collection)[i];
  //    }
  //    if (nhits) daughter_michel_score/=nhits;
  //    else daughter_michel_score = -999;
      if (evt.reco_beam_vertex_nHits) daughter_michel_score = evt.reco_beam_vertex_michel_score_weight_by_charge;//evt.reco_beam_vertex_michel_score/evt.reco_beam_vertex_nHits;
    }

    beam_dx = -999;
    beam_dy = -999;
    beam_dz = -999;
    beam_dxy = -999;
    beam_costh = -999;

    if (!evt.reco_beam_calo_wire->empty()){

      /*TVector3 pt0(evt.reco_beam_calo_startX,
                   evt.reco_beam_calo_startY,
                   evt.reco_beam_calo_startZ);
      TVector3 pt1(evt.reco_beam_calo_endX,
                   evt.reco_beam_calo_endY,
                   evt.reco_beam_calo_endZ);
      TVector3 dir = pt1 - pt0;
      dir = dir.Unit();*/
      //forced track info
      TVector3 pt0(evt.reco_beam_calo_startX,
                   evt.reco_beam_calo_startY,
                   evt.reco_beam_calo_startZ);
      TVector3 pt1(evt.reco_beam_calo_endX,
                   evt.reco_beam_calo_endY,
                   evt.reco_beam_calo_endZ);
      TVector3 dir = pt1 - pt0;
      dir = dir.Unit();

      if (evt.MC){
        TVector3 beamdir(cos(beam_angleX_mc*TMath::Pi()/180),
                         cos(beam_angleY_mc*TMath::Pi()/180),
                         cos(beam_angleZ_mc*TMath::Pi()/180));
        beamdir = beamdir.Unit();
        /*beam_dx = (evt.reco_beam_calo_startX - beam_startX_mc)/beam_startX_rms_mc;
        beam_dy = (evt.reco_beam_calo_startY - beam_startY_mc)/beam_startY_rms_mc;
        beam_dz = (evt.reco_beam_calo_startZ - beam_startZ_mc)/beam_startZ_rms_mc;
        beam_dxy = sqrt(pow(beam_dx,2) + pow(beam_dy,2));
        beam_costh = dir.Dot(beamdir);*/
        //forced track info
        beam_dx = (evt.reco_beam_calo_startX - beam_startX_mc)/beam_startX_rms_mc;
        beam_dy = (evt.reco_beam_calo_startY - beam_startY_mc)/beam_startY_rms_mc;
        beam_dz = (evt.reco_beam_calo_startZ - beam_startZ_mc)/beam_startZ_rms_mc;
        beam_dxy = sqrt(pow(beam_dx,2) + pow(beam_dy,2));
        beam_costh = dir.Dot(beamdir);
      }
      else{
        TVector3 beamdir(cos(beam_angleX_data*TMath::Pi()/180),
                         cos(beam_angleY_data*TMath::Pi()/180),
                         cos(beam_angleZ_data*TMath::Pi()/180));
        beamdir = beamdir.Unit();
        /*beam_dx = (evt.reco_beam_calo_startX - beam_startX_data)/beam_startX_rms_data;
        beam_dy = (evt.reco_beam_calo_startY - beam_startY_data)/beam_startY_rms_data;
        beam_dz = (evt.reco_beam_calo_startZ - beam_startZ_data)/beam_startZ_rms_data;
        beam_dxy = sqrt(pow(beam_dx,2) + pow(beam_dy,2));
        beam_costh = dir.Dot(beamdir);*/
        //forced track info
        beam_dx = (evt.reco_beam_calo_startX - beam_startX_data)/beam_startX_rms_data;
        beam_dy = (evt.reco_beam_calo_startY - beam_startY_data)/beam_startY_rms_data;
        beam_dz = (evt.reco_beam_calo_startZ - beam_startZ_data)/beam_startZ_rms_data;
        beam_dxy = sqrt(pow(beam_dx,2) + pow(beam_dy,2));
        beam_costh = dir.Dot(beamdir);
      }
    }

    dEdx_5cm = -1;
    /*
    if (!reco_beam_calibrated_dEdX_SCE->empty()){
      dEdx_5cm = 0;
      int nhits = 0;
      for (int i = 0; i<reco_beam_calibrated_dEdX_SCE->size(); ++i){
        std::cout<<i<<" "<<reco_beam_resRange->back()-(*reco_beam_resRange)[i]<<" "<<(*reco_beam_calibrated_dEdX_SCE)[i]<<endl;
        if (std::abs(reco_beam_resRange->back()-(*reco_beam_resRange)[i])<5){
          dEdx_5cm += (*reco_beam_calibrated_dEdX_SCE)[i];
          ++nhits;
        }
      }
      if (nhits) dEdx_5cm/=nhits;
      else dEdx_5cm = -1;
    }
    */

    //if (event == 78467) cout<<reco_beam_calibrated_dEdX_SCE->size()<<endl;
    //cout<<reco_beam_calibrated_dEdX_SCE->size()<<endl;
    if (!evt.reco_beam_calibrated_dEdX_SCE->empty()){ // what's this used for?
      //dEdx_5cm = 0;
      //int nhits = 0;
      std::vector<double> vdEdx;
      for (int i = 0; i<evt.reco_beam_calibrated_dEdX_SCE->size(); ++i){
        //std::cout<<i<<" "<<reco_beam_resRange->back()-(*reco_beam_resRange)[i]<<" "<<(*reco_beam_calibrated_dEdX_SCE)[i]<<endl;
        //if (event == 78467) cout<<(*reco_beam_resRange)[i]<<" "<<(*reco_beam_calo_Z)[i]<<" "<<(*reco_beam_calibrated_dEdX_SCE)[i]<<endl;
        if ((*evt.reco_beam_resRange_SCE)[i]<5){
          vdEdx.push_back((*evt.reco_beam_calibrated_dEdX_SCE)[i]);
          //dEdx_5cm += (*reco_beam_calibrated_dEdX_SCE)[i];
          //++nhits;
        }
      }
      //if (nhits) dEdx_5cm/=nhits;
      //else dEdx_5cm = -1;
      if (!vdEdx.empty()){
        dEdx_5cm = TMath::Median(vdEdx.size(), &vdEdx[0]);
      }
    }
  //  if (!MC && reco_beam_PFP_trackScore_collection>=0 && reco_beam_PFP_trackScore_collection<0.01){
  //    cout<<run<<" "<<event<<endl;
  //  }

    // calculate true track length
    int start_idx = -1;
    true_trklen_accum.reserve(evt.true_beam_traj_Z->size()); // initialize true_trklen_accum
    for (int i=0; i<evt.true_beam_traj_Z->size(); i++){
      if ((*evt.true_beam_traj_Z)[i] >= 0){
        start_idx = i-1; // the trajectory point before entering the TPC
        if (start_idx < 0) start_idx = -1;
        break;
      }
      true_trklen_accum[i] = 0.; // initialize true_trklen_accum
    }
    true_trklen = -1999; // initialize
    if (start_idx >= 0){
      for (int i=start_idx+1; i<evt.true_beam_traj_Z->size(); i++){
        if (i == start_idx+1) {
          true_trklen = sqrt( pow( (*evt.true_beam_traj_X)[i]-(*evt.true_beam_traj_X)[i-1], 2)
                              + pow( (*evt.true_beam_traj_Y)[i]-(*evt.true_beam_traj_Y)[i-1], 2)
                              + pow( (*evt.true_beam_traj_Z)[i]-(*evt.true_beam_traj_Z)[i-1], 2)
                              ) * (*evt.true_beam_traj_Z)[i]/((*evt.true_beam_traj_Z)[i]-(*evt.true_beam_traj_Z)[i-1]);
        }
        else{
          true_trklen += sqrt( pow( (*evt.true_beam_traj_X)[i]-(*evt.true_beam_traj_X)[i-1], 2)
                              + pow( (*evt.true_beam_traj_Y)[i]-(*evt.true_beam_traj_Y)[i-1], 2)
                              + pow( (*evt.true_beam_traj_Z)[i]-(*evt.true_beam_traj_Z)[i-1], 2)
                              );
        }
        true_trklen_accum[i] = true_trklen;
      }
    }

    // calculate reco track length
    reco_trklen_accum.reserve(evt.reco_beam_calo_Z->size());
    reco_trklen = -1999;
    for (int i=1; i<evt.reco_beam_calo_Z->size(); i++){
      if (i == 1) reco_trklen = 0;
      reco_trklen += sqrt( pow( (*evt.reco_beam_calo_X)[i]-(*evt.reco_beam_calo_X)[i-1], 2)
                          + pow( (*evt.reco_beam_calo_Y)[i]-(*evt.reco_beam_calo_Y)[i-1], 2)
                          + pow( (*evt.reco_beam_calo_Z)[i]-(*evt.reco_beam_calo_Z)[i-1], 2)
                          );
      reco_trklen_accum[i] = reco_trklen;
    }
    // front-face energy
    true_ffKE = 999999.;
    if (start_idx >= 0){
      true_ffKE = Get_true_ffKE(evt, (*evt.true_beam_traj_KE)[start_idx+1], (true_trklen_accum)[start_idx+1]);
    }
  }
  
  energy_calorimetry_SCE = 0; //MeV
  for (int i=0; i<evt.reco_beam_calibrated_dEdX_SCE->size(); i++){
    energy_calorimetry_SCE += (*evt.reco_beam_calibrated_dEdX_SCE)[i]*(*evt.reco_beam_TrkPitch_SCE)[i];
  }
  //cout<<evt.beam_particle_scores->size()<<"\t"<<(*evt.beam_particle_scores)[0]<<endl;
  if (evt.beam_particle_scores->size())
    beam_score = (*evt.beam_particle_scores)[0];
  else
    beam_score = -999.;
  //cout<<"$$$"<<evt.reco_beam_alt_len<<"\t"<<reco_trklen<<endl;//the two are the same
  // reco_trklen = evt.reco_beam_alt_len; // they should be the same
  if (fProtonCSDACheck)
    trklen_csda_proton = reco_trklen / csda_range_vs_mom_sm->Eval(evt.beam_inst_P);
  else
    trklen_csda_proton = -999;
}
