#include "Unfold.h"
#include "SliceParams.h"

Unfold::Unfold(TH2D* hist_reco, TH2D* hist_true)
  : response_SliceID_Int(hist_reco->ProjectionX(), hist_true->ProjectionX())
  , response_SliceID_Inc(hist_reco->ProjectionX(), hist_true->ProjectionX())
  , response_SliceID_Ini(hist_reco->ProjectionX(), hist_true->ProjectionX())
  , response_SliceID_2D(hist_reco, hist_true)
{

  response_SliceID_Int.UseOverflow(false);
  response_SliceID_Inc.UseOverflow(false);
  response_SliceID_Ini.UseOverflow(false);
  response_SliceID_2D.UseOverflow(false);

  eff_num_Int = new TH1D("eff_num_Int", "eff_num_Int", pi::true_nbins, pi::true_bins);
  eff_den_Int = new TH1D("eff_den_Int", "eff_den_Int", pi::true_nbins, pi::true_bins);
  eff_num_Inc = new TH1D("eff_num_Inc", "eff_num_Inc", pi::true_nbins, pi::true_bins);
  eff_den_Inc = new TH1D("eff_den_Inc", "eff_den_Inc", pi::true_nbins, pi::true_bins);
  pur_num_Int = new TH1D("pur_num_Int", "pur_num_Int", pi::reco_nbins, pi::reco_bins);
  pur_num_Inc = new TH1D("pur_num_Inc", "pur_num_Inc", pi::reco_nbins, pi::reco_bins);
  pur_den     = new TH1D("pur_den",     "pur_den",     pi::reco_nbins, pi::reco_bins);

  eff_num_Int->Sumw2();
  eff_den_Int->Sumw2();
  eff_num_Inc->Sumw2();
  eff_den_Inc->Sumw2();
  pur_num_Int->Sumw2();
  pur_num_Inc->Sumw2();
  pur_den->Sumw2();

}  

void Unfold::SaveHistograms(){

  eff_num_Int->Write("eff_num_Int");
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
  hint->Write("hresponse_SliceID_Int");
  TH2D *hinc = (TH2D*)response_SliceID_Inc.Hresponse();
  hinc->SetTitle("Incidents; Reco Slice ID; True Slice ID");
  hinc->Write("hresponse_SliceID_Inc");
  TH2D *hini = (TH2D*)response_SliceID_Ini.Hresponse();
  hini->SetTitle("Initial; Reco Slice ID; True Slice ID");
  hini->Write("hresponse_SliceID_Ini");
  TH2D *h2D = (TH2D*)response_SliceID_2D.Hresponse();
  h2D->SetTitle("2D; Reco Slice ID; True Slice ID");
  h2D->Write("hresponse_SliceID_2D");

  response_SliceID_Int.Write("response_SliceID_Int");
  response_SliceID_Inc.Write("response_SliceID_Inc");
  response_SliceID_Ini.Write("response_SliceID_Ini");
  response_SliceID_2D.Write("response_SliceID_2D");
}
