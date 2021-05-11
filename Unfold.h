#ifndef UNFOLD_H
#define UNFOLD_H

#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"

class Unfold {
  
 public:

  Unfold(int nb, double xlo, double xhi);

  RooUnfoldResponse response_SliceID_Int;  //Interaction
  RooUnfoldResponse response_SliceID_Inc;  //Incident

  TH1D *eff_num_Int; //Interaction efficiency numerator
  TH1D *eff_den_Int; //Interaction efficiency denominator

  TH1D *eff_num_Inc; //Incident efficiency numerator
  TH1D *eff_den_Inc; //Incident efficiency denominator

  TH1D *pur_num_Int; //Interaction purity numerator
  TH1D *pur_num_Inc; //Incident purity numerator
  TH1D *pur_den;     //Interaction/incident purity denominator


};

#endif
