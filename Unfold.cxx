#include "Unfold.h"

Unfold::Unfold(int nb, double xlo, double xhi)
  : response_SliceID_Int(nb, xlo, xhi)
  , response_SliceID_Inc(nb, xlo, xhi)
{

  response_SliceID_Int.UseOverflow(false);//do not use the -999 points
  response_SliceID_Inc.UseOverflow(false);

  eff_num_Int = new TH1D("eff_num_Int", "eff_num_Int", nb, xlo, xhi);
  eff_den_Int = new TH1D("eff_den_Int", "eff_den_Int", nb, xlo, xhi);
  eff_num_Inc = new TH1D("eff_num_Inc", "eff_num_Inc", nb, xlo, xhi);
  eff_den_Inc = new TH1D("eff_den_Inc", "eff_den_Inc", nb, xlo, xhi);
  pur_num_Int = new TH1D("pur_num_Int", "pur_num_Int", nb, xlo, xhi);
  pur_num_Inc = new TH1D("pur_num_Inc", "pur_num_Inc", nb, xlo, xhi);
  pur_den     = new TH1D("pur_den",     "pur_den",     nb, xlo, xhi);

  eff_num_Int->Sumw2();//Create structure to store sum of squares of weights (for error estimation)
  eff_den_Int->Sumw2();
  eff_num_Inc->Sumw2();
  eff_den_Inc->Sumw2();
  pur_num_Int->Sumw2();
  pur_num_Inc->Sumw2();
  pur_den->Sumw2();

}  

void Unfold::SaveHistograms(){

  eff_num_Int->Write("eff_num_Int");//Write this object to the current directory
  eff_den_Int->Write("eff_den_Int");
  eff_num_Inc->Write("eff_num_Inc");
  eff_den_Inc->Write("eff_den_Inc");
  pur_num_Int->Write("pur_num_Int");
  pur_num_Inc->Write("pur_num_Inc");
  pur_den->Write("pur_den");

  eff_Int = (TH1D*)eff_num_Int->Clone("eff_Int");
  eff_Int->Divide(eff_den_Int);
  eff_Int->Write("eff_Int");

  eff_Inc = (TH1D*)eff_num_Inc->Clone("eff_Inc");
  eff_Inc->Divide(eff_den_Inc);
  eff_Inc->Write("eff_Inc");

  pur_Int = (TH1D*)pur_num_Int->Clone("pur_Int");
  pur_Int->Divide(pur_den);
  pur_Int->Write("pur_Int");

  pur_Inc = (TH1D*)pur_num_Inc->Clone("pur_Inc");
  pur_Inc->Divide(pur_den);
  pur_Inc->Write("pur_Inc");

  TH2D *hint = (TH2D*)response_SliceID_Int.Hresponse();
  hint->SetTitle("Interactions;Reco Slice ID;True Slice ID");
  hint->Write("response_SliceID_Int");
  TH2D *hinc = (TH2D*)response_SliceID_Inc.Hresponse();
  hinc->SetTitle("Incidents; Reco Slice ID; True Slice ID");
  hinc->Write("response_SliceID_Inc");
}
