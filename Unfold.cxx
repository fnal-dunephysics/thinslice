#include "Unfold.h"
#include "SliceParams.h"

Unfold::Unfold(TH3D* hist_reco, TH3D* hist_true)
  : response_SliceID_Int(hist_reco->Project3D("x"), hist_true->Project3D("x"))
  , response_SliceID_Inc(hist_reco->Project3D("x"), hist_true->Project3D("x"))
  , response_SliceID_Ini(hist_reco->Project3D("x"), hist_true->Project3D("x"))
  , response_SliceID_3D(hist_reco, hist_true)
  , response_SliceID_1D(53,0,53, 94,0,94) // edit here for 1D unfolding
  , response_SliceID_1D_noeff(53,0,53, 94,0,94)
{

  response_SliceID_Int.UseOverflow(false);
  response_SliceID_Inc.UseOverflow(false);
  response_SliceID_Ini.UseOverflow(false);
  response_SliceID_3D.UseOverflow(false);
  response_SliceID_1D.UseOverflow(false);
  response_SliceID_1D_noeff.UseOverflow(false);

  eff_num_Int = new TH1D("eff_num_Int", "eff_num_Int", pi::true_nbins, pi::true_bins);
  eff_den_Int = new TH1D("eff_den_Int", "eff_den_Int", pi::true_nbins, pi::true_bins);
  pur_num_Int = new TH1D("pur_num_Int", "pur_num_Int", pi::reco_nbins, pi::reco_bins);
  pur_den_Int = new TH1D("pur_den_Int", "pur_den_Int", pi::reco_nbins, pi::reco_bins);
  eff_num_Inc = new TH1D("eff_num_Inc", "eff_num_Inc", pi::true_nbins, pi::true_bins);
  eff_den_Inc = new TH1D("eff_den_Inc", "eff_den_Inc", pi::true_nbins, pi::true_bins);
  pur_num_Inc = new TH1D("pur_num_Inc", "pur_num_Inc", pi::reco_nbins, pi::reco_bins);
  pur_den_Inc = new TH1D("pur_den_Inc", "pur_den_Inc", pi::reco_nbins, pi::reco_bins);
  eff_num_Ini = new TH1D("eff_num_Ini", "eff_num_Ini", pi::true_nbins, pi::true_bins);
  eff_den_Ini = new TH1D("eff_den_Ini", "eff_den_Ini", pi::true_nbins, pi::true_bins);
  pur_num_Ini = new TH1D("pur_num_Ini", "pur_num_Ini", pi::reco_nbins, pi::reco_bins);
  pur_den_Ini = new TH1D("pur_den_Ini", "pur_den_Ini", pi::reco_nbins, pi::reco_bins);
  

  eff_num_Int->Sumw2();
  eff_den_Int->Sumw2();
  pur_num_Int->Sumw2();
  pur_den_Int->Sumw2();
  eff_num_Inc->Sumw2();
  eff_den_Inc->Sumw2();
  pur_num_Inc->Sumw2();
  pur_den_Inc->Sumw2();
  eff_num_Ini->Sumw2();
  eff_den_Ini->Sumw2();
  pur_num_Ini->Sumw2();
  pur_den_Ini->Sumw2();

}  

void Unfold::SaveHistograms(){

  eff_num_Int->Write("eff_num_Int");
  eff_den_Int->Write("eff_den_Int");
  pur_num_Int->Write("pur_num_Int");
  pur_den_Int->Write("pur_den_Int");
  eff_num_Inc->Write("eff_num_Inc");
  eff_den_Inc->Write("eff_den_Inc");
  pur_num_Inc->Write("pur_num_Inc");
  pur_den_Inc->Write("pur_den_Inc");
  eff_num_Ini->Write("eff_num_Ini");
  eff_den_Ini->Write("eff_den_Ini");
  pur_num_Ini->Write("pur_num_Ini");
  pur_den_Ini->Write("pur_den_Ini");

  eff_Int = (TH1D*)eff_num_Int->Clone("eff_Int");
  eff_Int->Divide(eff_den_Int);
  eff_Int->SetTitle("eff_Int");
  eff_Int->Write("eff_Int");

  eff_Inc = (TH1D*)eff_num_Inc->Clone("eff_Inc");
  eff_Inc->Divide(eff_den_Inc);
  eff_Inc->SetTitle("eff_Inc");
  eff_Inc->Write("eff_Inc");
  
  eff_Ini = (TH1D*)eff_num_Ini->Clone("eff_Ini");
  eff_Ini->Divide(eff_den_Ini);
  eff_Ini->SetTitle("eff_Ini");
  eff_Ini->Write("eff_Ini");

  pur_Int = (TH1D*)pur_num_Int->Clone("pur_Int");
  pur_Int->Divide(pur_den_Int);
  pur_Int->SetTitle("pur_Int");
  pur_Int->Write("pur_Int");

  pur_Inc = (TH1D*)pur_num_Inc->Clone("pur_Inc");
  pur_Inc->Divide(pur_den_Inc);
  pur_Inc->SetTitle("pur_Inc");
  pur_Inc->Write("pur_Inc");
  
  pur_Ini = (TH1D*)pur_num_Ini->Clone("pur_Ini");
  pur_Ini->Divide(pur_den_Ini);
  pur_Ini->SetTitle("pur_Ini");
  pur_Ini->Write("pur_Ini");

  TH2D *hint = (TH2D*)response_SliceID_Int.Hresponse();
  hint->SetTitle("Interactions;Reco Slice ID;True Slice ID");
  hint->Write("hresponse_SliceID_Int");
  TH2D *hinc = (TH2D*)response_SliceID_Inc.Hresponse();
  hinc->SetTitle("Incidents; Reco Slice ID; True Slice ID");
  hinc->Write("hresponse_SliceID_Inc");
  TH2D *hini = (TH2D*)response_SliceID_Ini.Hresponse();
  hini->SetTitle("Initial; Reco Slice ID; True Slice ID");
  hini->Write("hresponse_SliceID_Ini");
  TH2D *h3D = (TH2D*)response_SliceID_3D.Hresponse();
  h3D->SetTitle("3D; Reco Slice ID; True Slice ID");
  h3D->Write("hresponse_SliceID_3D");
  TH2D *h1D = (TH2D*)response_SliceID_1D.Hresponse();
  h1D->SetTitle("1D; Reco Slice ID; True Slice ID");
  h1D->Write("hresponse_SliceID_1D");
  TH2D *h1D_noeff = (TH2D*)response_SliceID_1D_noeff.Hresponse();
  h1D_noeff->SetTitle("1D; Reco Slice ID; True Slice ID");
  h1D_noeff->Write("hresponse_SliceID_1D_noeff");

  response_SliceID_Int.Write("response_SliceID_Int");
  response_SliceID_Inc.Write("response_SliceID_Inc");
  response_SliceID_Ini.Write("response_SliceID_Ini");
  response_SliceID_3D.Write("response_SliceID_3D");
  response_SliceID_1D.Write("response_SliceID_1D");
  response_SliceID_1D_noeff.Write("response_SliceID_1D_noeff");
}
