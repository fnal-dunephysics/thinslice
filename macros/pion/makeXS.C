#include "../../SliceParams.h"

void makeXS(){

  gStyle->SetOptStat(0);

  TFile *file = TFile::Open("../../build/mcprod4a.root");

  TH1D *eff_num_Int = (TH1D*)file->Get("eff_num_Int");
  TH1D *eff_den_Int = (TH1D*)file->Get("eff_den_Int");
  TH1D *eff_num_Inc = (TH1D*)file->Get("eff_num_Inc");
  TH1D *eff_den_Inc = (TH1D*)file->Get("eff_den_Inc");
  TH1D *pur_num_Int = (TH1D*)file->Get("pur_num_Int");
  TH1D *pur_num_Inc = (TH1D*)file->Get("pur_num_Inc");
  TH1D *pur_den = (TH1D*)file->Get("pur_den");

  TH1D *reco_AngCorr = (TH1D*)file->Get("reco_AngCorr");
  TH1D *true_AngCorr = (TH1D*)file->Get("true_AngCorr");

  TH1D *h_truesliceid_pion_uf = (TH1D*)file->Get("h_truesliceid_pion_uf");
  TH1D *h_truesliceid_pioninelastic_uf = (TH1D*)file->Get("h_truesliceid_pioninelastic_uf");

  TH1D *eff_Int = (TH1D*)eff_num_Int->Clone("eff_Int");
  eff_Int->Divide(eff_den_Int);

  TH1D *eff_Inc = (TH1D*)eff_num_Inc->Clone("eff_Inc");
  eff_Inc->Divide(eff_den_Inc);

  TH1D *pur_Int = (TH1D*)pur_num_Int->Clone("pur_Int");
  pur_Int->Divide(pur_den);

  TH1D *pur_Inc = (TH1D*)pur_num_Inc->Clone("pur_Inc");
  pur_Inc->Divide(pur_den);

  TH1D *h_truesliceid_pion_all = (TH1D*)file->Get("h_truesliceid_pion_all");
  TH1D *h_truesliceid_pion_cuts = (TH1D*)file->Get("h_truesliceid_pion_cuts");
  TH1D *h_truesliceid_pioninelastic_all = (TH1D*)file->Get("h_truesliceid_pioninelastic_all");
  TH1D *h_truesliceid_pioninelastic_cuts = (TH1D*)file->Get("h_truesliceid_pioninelastic_cuts");
  TH1D *h_recosliceid_allevts_cuts = (TH1D*)file->Get("h_recosliceid_allevts_cuts");
  TH1D *h_recosliceid_pion_cuts = (TH1D*)file->Get("h_recosliceid_pion_cuts");
  TH1D *h_recosliceid_pioninelastic_cuts = (TH1D*)file->Get("h_recosliceid_pioninelastic_cuts");

  TH1D *hinc = (TH1D*)h_recosliceid_allevts_cuts->Clone("hinc");
  TH1D *hint = (TH1D*)h_recosliceid_allevts_cuts->Clone("hint");

  hinc->Multiply(pur_Inc);
  /*TCanvas *c1 = new TCanvas("c1","c1");
  h_recosliceid_allevts_cuts->SetLineColor(3);
  h_recosliceid_allevts_cuts->SetMarkerColor(3);
  h_recosliceid_allevts_cuts->SetTitle("All Pions;Reco SliceID;Events");
  h_recosliceid_allevts_cuts->DrawCopy();
  hinc->DrawCopy("same");
  h_recosliceid_pion_cuts->SetLineColor(2);
  h_recosliceid_pion_cuts->Draw("same hist");
  TLegend *leg1 = new TLegend(0.5,0.6,0.8,0.9);
  leg1->SetFillStyle(0);
  leg1->AddEntry(h_recosliceid_allevts_cuts,"Selected","ple");
  leg1->AddEntry(hinc, "Selected #times purity","ple");
  leg1->AddEntry(h_recosliceid_pion_cuts,"Selected true pions","l");
  leg1->Draw();

  TCanvas *c2 = new TCanvas("c2","c2");
  h_truesliceid_pion_uf->SetLineColor(4);
  h_truesliceid_pion_uf->SetMarkerColor(4);
  h_truesliceid_pion_uf->SetTitle("All Pions; True SliceID; Events");
  h_truesliceid_pion_uf->Draw();
  //hinc->SetLineColor(3);
  //hinc->SetMarkerColor(3);
  hinc->DrawCopy("same");
//  eff_den_Inc->SetLineColor(2);
//  eff_den_Inc->SetMarkerColor(2);
//  eff_den_Inc->Draw("same hist");
  h_truesliceid_pion_all->SetLineColor(2);
  h_truesliceid_pion_all->SetMarkerColor(2);
  h_truesliceid_pion_all->Draw("same hist");
  TLegend *leg2 = new TLegend(0.5,0.6,0.8,0.9);
  leg2->SetFillStyle(0);
  leg2->AddEntry(hinc, "Selected #times purity","ple");
  leg2->AddEntry(h_truesliceid_pion_uf,"Unfolded pions","ple");
  leg2->AddEntry(h_truesliceid_pion_all,"True pions","l");
  leg2->Draw();

  hint->Multiply(pur_Int);
  TCanvas *c3 = new TCanvas("c3","c3");
  h_recosliceid_allevts_cuts->SetLineColor(3);
  h_recosliceid_allevts_cuts->SetMarkerColor(3);
  h_recosliceid_allevts_cuts->SetTitle("Pion Inelastic Scatterings;Reco SliceID;Events");
  h_recosliceid_allevts_cuts->Draw();
  hint->DrawCopy("same");
  h_recosliceid_pioninelastic_cuts->SetLineColor(2);
  h_recosliceid_pioninelastic_cuts->Draw("same hist");
  TLegend *leg3 = new TLegend(0.5,0.6,0.8,0.9);
  leg3->SetFillStyle(0);
  leg3->AddEntry(h_recosliceid_allevts_cuts,"Selected","ple");
  leg3->AddEntry(hint, "Selected #times purity","ple");
  leg3->AddEntry(h_recosliceid_pioninelastic_cuts,"Selected true pions","l");
  leg3->Draw();

  TCanvas *c4 = new TCanvas("c4","c4");
  h_truesliceid_pioninelastic_uf->SetLineColor(4);
  h_truesliceid_pioninelastic_uf->SetMarkerColor(4);
  h_truesliceid_pioninelastic_uf->SetTitle("Pion Inelastic Scatterings; True SliceID; Events");
  h_truesliceid_pioninelastic_uf->Draw();
  //hint->SetLineColor(3);
  //hint->SetMarkerColor(3);
  hint->DrawCopy("same");
  h_truesliceid_pioninelastic_all->SetLineColor(2);
  h_truesliceid_pioninelastic_all->SetMarkerColor(2);
  h_truesliceid_pioninelastic_all->Draw("same hist");
  TLegend *leg4 = new TLegend(0.5,0.6,0.8,0.9);
  leg4->SetFillStyle(0);
  leg4->AddEntry(hint, "Selected #times purity","ple");
  leg4->AddEntry(h_truesliceid_pioninelastic_uf,"Unfolded pions","ple");
  leg4->AddEntry(h_truesliceid_pioninelastic_all,"True pions","l");
  leg4->Draw();*/

  /*double Ninc[pi::nthinslices] = {0};
  double Nint[pi::nthinslices] = {0};
  double err_inc[pi::nthinslices] = {0};
  double err_int[pi::nthinslices] = {0};
  double SliceID[pi::nthinslices] = {0};

  for (int i = 0; i<pi::nthinslices; ++i){
    SliceID[i] = i;
    Nint[i] = h_truesliceid_pioninelastic_uf->GetBinContent(i+2);
    err_int[i] = h_truesliceid_pioninelastic_uf->GetBinError(i+2);
    for (int j = i; j<=pi::nthinslices; ++j){
      Ninc[i] += h_truesliceid_pion_uf->GetBinContent(j+2);
      err_inc[i] += pow(h_truesliceid_pion_uf->GetBinError(j+2),2);
      //cout<<i<<" "<<j<<" "<<h_truesliceid_pion_uf->GetBinContent(j+2)<<" "<<Ninc[i]<<endl;
    }
    err_inc[i] = sqrt(err_inc[i]);
  }

  TGraphErrors *gr_inc = new TGraphErrors(pi::nthinslices, SliceID, Ninc, 0, err_inc);
  TGraphErrors *gr_int = new TGraphErrors(pi::nthinslices, SliceID, Nint, 0, err_int);
  TCanvas *c6 = new TCanvas("c6","c6");
  gr_inc->SetTitle("");
  gr_inc->SetLineWidth(2);
  gr_inc->SetLineColor(4);
  gr_inc->SetMarkerColor(4);
  gr_inc->GetXaxis()->SetTitle("Slice ID");
  gr_inc->GetYaxis()->SetTitle("N_{Inc}");
  gr_inc->GetXaxis()->SetRangeUser(-1, 23);
  gr_inc->Draw("ape");

  TCanvas *c7 = new TCanvas("c7","c7");
  gr_int->SetTitle("");
  gr_int->SetLineWidth(2);
  gr_int->SetLineColor(4);
  gr_int->SetMarkerColor(4);
  gr_int->GetXaxis()->SetTitle("Slice ID");
  gr_int->GetYaxis()->SetTitle("N_{Int}");
  gr_int->GetXaxis()->SetRangeUser(-1, 23);
  gr_int->Draw("ape");

  double NA=6.02214076e23;
  double MAr=39.95; //gmol
  double Density = 1.4; // g/cm^3

  double xs[pi::nthinslices] = {0};
  double err_xs[pi::nthinslices] = {0};
  double incE[pi::nthinslices] = {0};

  TGraphErrors *gr_trueincE = (TGraphErrors*)file->Get("gr_trueincE");*/
  TGraphErrors *gr_truexs = (TGraphErrors*)file->Get("gr_truexs");
  /*for (int i = 0; i<pi::nthinslices; ++i){
    xs[i] = MAr/(Density*NA*pi::thinslicewidth)*log(Ninc[i]/(Ninc[i]-Nint[i]))*1e27;
    //err_xs[i] = MAr/(Density*NA*pi::thinslicewidth)*1e27*sqrt(N_int[i]+pow(N_int[i],2)/N_inc[i])/N_incidents[i];
    err_xs[i] = MAr/(Density*NA*pi::thinslicewidth)*1e27*sqrt(pow(Nint[i]*err_inc[i]/Ninc[i]/(Ninc[i]-Nint[i]),2)+pow(err_int[i]/(Ninc[i]-Nint[i]),2));
    incE[i] = gr_trueincE->GetPointY(i);
    //std::cout<<i<<" "<<Ninc[i]<<" "<<Nint[i]<<" "<<xs[i]<<" "<<incE[i]<<std::endl;
  }*/

  TFile f2("../../files/exclusive_xsec.root");
  TGraph *total_inel_KE = (TGraph*)f2.Get("total_inel_KE");
//  TGraph *abs_KE = (TGraph*)f2.Get("abs_KE");
//  TGraph *cex_KE = (TGraph*)f2.Get("cex_KE");

  //TGraphErrors *gr_recoxs = new TGraphErrors(pi::nthinslices, incE, xs, 0, err_xs);
  TCanvas *c5 = new TCanvas("c5", "c5", 1200, 500);
  gr_truexs->SetTitle("Pion Inelastic Cross Section");
  gr_truexs->GetXaxis()->SetTitle("Pion Kinetic Energy (MeV)");
  gr_truexs->GetXaxis()->SetRangeUser(10, 1000);
  gr_truexs->GetYaxis()->SetTitle("#sigma_{inelastic} (mb)");
  gr_truexs->GetYaxis()->SetRangeUser(0, 1000);
  gr_truexs->SetLineWidth(2);
  gr_truexs->Draw("ape");
  gr_truexs->SetMarkerColor(3);
  gr_truexs->SetLineColor(3);
  //gr_truexs->Draw("pe");
  total_inel_KE->SetLineColor(2);
  total_inel_KE->Draw("c");
  TLegend *leg5 = new TLegend(0.3,0.65,0.8,0.9);
  leg5->SetFillStyle(0);
  //leg5->AddEntry(gr_recoxs, "MC with reconstruction", "pe");
  leg5->AddEntry(gr_truexs, "MC truth", "pe");
  leg5->AddEntry(total_inel_KE, "Geant4 (theory prediction)", "l");
  leg5->Draw();

  /*c1->Print("plots/xs_sliceidinc_reco.pdf");
  c2->Print("plots/xs_sliceidinc_true.pdf");
  c3->Print("plots/xs_sliceidint_reco.pdf");
  c4->Print("plots/xs_sliceidint_true.pdf");*/
  c5->Print("plots/xs_pi+inel.pdf");
  //c6->Print("plots/xs_Ninc.pdf");
  //c7->Print("plots/xs_Nint.pdf");

  /*c1->Print("plots/xs_sliceidinc_reco.png");
  c2->Print("plots/xs_sliceidinc_true.png");
  c3->Print("plots/xs_sliceidint_reco.png");
  c4->Print("plots/xs_sliceidint_true.png");*/
  c5->Print("plots/xs_pi+inel.png");
  //c6->Print("plots/xs_Ninc.png");
  //c7->Print("plots/xs_Nint.png");

}
  
