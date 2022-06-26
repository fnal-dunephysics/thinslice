#ifndef UNFOLD_H
#define UNFOLD_H

#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"
#include "TH2D.h"
#include "TH3D.h"

class Unfold {
  
 public:

  Unfold(TH3D* hist_reco, TH3D* hist_true);

  RooUnfoldResponse response_SliceID_Int;  //Interaction
  RooUnfoldResponse response_SliceID_Inc;  //Incident
  RooUnfoldResponse response_SliceID_Ini;  //Initial
  RooUnfoldResponse response_SliceID_3D;  //3D Initial vs Incident vs Interacting

  TH1D *eff_num_Int; //Interaction efficiency numerator
  TH1D *eff_den_Int; //Interaction efficiency denominator

  TH1D *eff_num_Inc; //Incident efficiency numerator
  TH1D *eff_den_Inc; //Incident efficiency denominator

  TH1D *pur_num_Int; //Interaction purity numerator
  TH1D *pur_num_Inc; //Incident purity numerator
  TH1D *pur_den;     //Interaction/incident purity denominator

  TH1D *eff_Int;
  TH1D *eff_Inc;
  TH1D *pur_Int;
  TH1D *pur_Inc;

  void SaveHistograms();

};

#endif
