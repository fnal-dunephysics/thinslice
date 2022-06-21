#include "../../SliceParams.h"

void plotXS(){

  gStyle->SetOptStat(0);

  TFile *file = TFile::Open("../../build/XSMC.root");
  
  TH1D *h_sel_data = (TH1D*)file->Get("hdata");
  TH1D *h_sel_sig_int = (TH1D*)file->Get("hsignal");
  TH1D *hval_sel_sig_int = (TH1D*)file->Get("hval_signal_reco");
  TH1D *h_sel_sig_int_uf = (TH1D*)file->Get("hsignal_uf");
  TH1D *hval_sel_sig_int_uf = (TH1D*)file->Get("hval_trueint");
  TH1D *h_sel_data_inc = (TH1D*)file->Get("hdata_inc");
  TH1D *h_sel_sig_inc = (TH1D*)file->Get("hsiginc");
  TH1D *hval_sel_sig_inc = (TH1D*)file->Get("hval_siginc_reco");
  TH1D *h_sel_sig_inc_uf = (TH1D*)file->Get("hsiginc_uf");
  TH1D *hval_sel_sig_inc_uf = (TH1D*)file->Get("hval_trueinc");
  TH1D *h_sel_data_ini = (TH1D*)file->Get("hdata_ini");
  TH1D *h_sel_sig_ini = (TH1D*)file->Get("hsigini");
  TH1D *hval_sel_sig_ini = (TH1D*)file->Get("hval_sigini_reco");
  TH1D *h_sel_sig_ini_uf = (TH1D*)file->Get("hsigini_uf");
  TH1D *hval_sel_sig_ini_uf = (TH1D*)file->Get("hval_trueini");
  TGraphErrors *gr_inc = (TGraphErrors*)file->Get("gr_inc");
  TGraphErrors *gr_int = (TGraphErrors*)file->Get("gr_int");
  TGraphErrors *gr_ini = (TGraphErrors*)file->Get("gr_ini");
  TGraphErrors *gr_ina = (TGraphErrors*)file->Get("gr_ina");
  TGraphErrors *gr_inc_t = (TGraphErrors*)file->Get("gr_inc_t");
  TGraphErrors *gr_int_t = (TGraphErrors*)file->Get("gr_int_t");
  TGraphErrors *gr_ini_t = (TGraphErrors*)file->Get("gr_ini_t");
  TGraphErrors *gr_ina_t = (TGraphErrors*)file->Get("gr_ina_t");
  
  /*for (int i=0; i<pi::nthinslices; ++i){
    cout<<gr_inc->GetPointY(i)<<",";
    //cout<<hval_sel_sig_int_uf->GetBinContent(i+2)<<",";
  }*/
  
  // incident (all pion)
  TCanvas *c1 = new TCanvas("c1","c1");
  h_sel_data_inc->SetTitle("All Pions;Reco SliceID;Events");
  h_sel_data_inc->DrawCopy();
  h_sel_sig_inc->SetLineColor(3);
  h_sel_sig_inc->SetMarkerColor(3);
  h_sel_sig_inc->DrawCopy("same");
  hval_sel_sig_inc->SetLineColor(2);
  //hval_sel_sig_inc->Scale(h_sel_sig_inc->Integral()/hval_sel_sig_inc->Integral());
  hval_sel_sig_inc->Draw("same hist");
  TLegend *leg1 = new TLegend(0.5,0.6,0.8,0.9);
  leg1->SetFillStyle(0);
  leg1->AddEntry(h_sel_data_inc,"Selected data","ple");
  leg1->AddEntry(h_sel_sig_inc, "After bkg subtraction","ple");
  leg1->AddEntry(hval_sel_sig_inc,"Selected true pions","l");
  leg1->Draw();

  TCanvas *c2 = new TCanvas("c2","c2");
  h_sel_sig_inc_uf->SetLineColor(4);
  h_sel_sig_inc_uf->SetMarkerColor(4);
  h_sel_sig_inc_uf->SetTitle("All Pions; True SliceID; Events");
  h_sel_sig_inc_uf->Draw();
  h_sel_sig_inc->DrawCopy("same");
  hval_sel_sig_inc_uf->SetLineColor(2);
  hval_sel_sig_inc_uf->SetMarkerColor(2);
  //hval_sel_sig_inc_uf->Scale(h_sel_sig_inc_uf->Integral()/hval_sel_sig_inc_uf->Integral());
  hval_sel_sig_inc_uf->Draw("same hist");
  TLegend *leg2 = new TLegend(0.5,0.6,0.8,0.9);
  leg2->SetFillStyle(0);
  leg2->AddEntry(h_sel_sig_inc, "Selected signal","ple");
  leg2->AddEntry(h_sel_sig_inc_uf,"After unfolding","ple");
  leg2->AddEntry(hval_sel_sig_inc_uf,"True pions","l");
  leg2->Draw();

  // interaction (pion inelastic)
  TCanvas *c3 = new TCanvas("c3","c3"); // compared to c1, only miss a few pion elastic events
  h_sel_data->SetTitle("Pion Inelastic Scatterings;Reco SliceID;Events");
  h_sel_data->Draw();
  h_sel_sig_int->SetLineColor(3);
  h_sel_sig_int->SetMarkerColor(3);
  h_sel_sig_int->DrawCopy("same");
  hval_sel_sig_int->SetLineColor(2);
  //hval_sel_sig_int->Scale(h_sel_sig_int->Integral()/hval_sel_sig_int->Integral());
  hval_sel_sig_int->Draw("same hist");
  TLegend *leg3 = new TLegend(0.5,0.6,0.8,0.9);
  leg3->SetFillStyle(0);
  leg3->AddEntry(h_sel_data,"Selected data","ple");
  leg3->AddEntry(h_sel_sig_int, "After bkg subtraction","ple");
  leg3->AddEntry(hval_sel_sig_int,"Selected true pions","l");
  leg3->Draw();

  TCanvas *c4 = new TCanvas("c4","c4");
  h_sel_sig_int_uf->SetLineColor(4);
  h_sel_sig_int_uf->SetMarkerColor(4);
  h_sel_sig_int_uf->SetTitle("Pion Inelastic Scatterings; True SliceID; Events");
  h_sel_sig_int_uf->Draw();
  h_sel_sig_int->DrawCopy("same");
  hval_sel_sig_int_uf->SetLineColor(2);
  hval_sel_sig_int_uf->SetMarkerColor(2);
  //hval_sel_sig_int_uf->Scale(h_sel_sig_int_uf->Integral()/hval_sel_sig_int_uf->Integral());
  hval_sel_sig_int_uf->Draw("same hist");
  TLegend *leg4 = new TLegend(0.5,0.6,0.8,0.9);
  leg4->SetFillStyle(0);
  leg4->AddEntry(h_sel_sig_int, "Selected signal","ple");
  leg4->AddEntry(h_sel_sig_int_uf,"After unfolding","ple");
  leg4->AddEntry(hval_sel_sig_int_uf,"True pions","l");
  leg4->Draw();
  
  // start
  TCanvas *c5 = new TCanvas("c5","c5");
  h_sel_data_ini->SetTitle("All Pions;Reco SliceID;Events");
  h_sel_data_ini->Draw();
  h_sel_sig_ini->SetLineColor(3);
  h_sel_sig_ini->SetMarkerColor(3);
  h_sel_sig_ini->DrawCopy("same");
  hval_sel_sig_ini->SetLineColor(2);
  //hval_sel_sig_ini->Scale(h_sel_sig_ini->Integral()/hval_sel_sig_ini->Integral());
  hval_sel_sig_ini->Draw("same hist");
  TLegend *leg5 = new TLegend(0.5,0.6,0.8,0.9);
  leg5->SetFillStyle(0);
  leg5->AddEntry(h_sel_data_ini,"Selected data","ple");
  leg5->AddEntry(h_sel_sig_ini, "After bkg subtraction","ple");
  leg5->AddEntry(hval_sel_sig_ini,"Selected true pions","l");
  leg5->Draw();

  TCanvas *c6 = new TCanvas("c6","c6");
  h_sel_sig_ini_uf->SetLineColor(4);
  h_sel_sig_ini_uf->SetMarkerColor(4);
  h_sel_sig_ini_uf->SetTitle("All Pions; True SliceID; Events");
  h_sel_sig_ini_uf->SetMinimum(0);
  h_sel_sig_ini_uf->Draw();
  h_sel_sig_ini->DrawCopy("same");
  hval_sel_sig_ini_uf->SetLineColor(2);
  hval_sel_sig_ini_uf->SetMarkerColor(2);
  //hval_sel_sig_ini_uf->Scale(h_sel_sig_ini_uf->Integral()/hval_sel_sig_ini_uf->Integral());
  hval_sel_sig_ini_uf->Draw("same hist");
  TLegend *leg6 = new TLegend(0.5,0.6,0.8,0.9);
  leg6->SetFillStyle(0);
  leg6->AddEntry(h_sel_sig_ini, "Selected signal","ple");
  leg6->AddEntry(h_sel_sig_ini_uf,"After unfolding","ple");
  leg6->AddEntry(hval_sel_sig_ini_uf,"True pions","l");
  leg6->Draw();
  /*for (int i=1; i<=22; ++i) {
    cout<<h_sel_sig_ini_uf->GetBinContent(i)<<"\t"<<h_sel_sig_ini_uf->GetBinError(i)<<endl;
  }*/

  // incident histogram
  TCanvas *c7 = new TCanvas("c7","c7");
  gr_inc->SetTitle("");
  gr_inc->SetLineWidth(2);
  gr_inc->GetXaxis()->SetTitle("Slice ID");
  gr_inc->GetYaxis()->SetTitle("N_{Inc}");
  gr_inc->GetXaxis()->SetRangeUser(0, 21);
  gr_inc->SetMinimum(0);
  gr_inc->Draw("ape");
  gr_inc_t->SetLineWidth(2);
  gr_inc_t->SetLineColor(3);
  gr_inc_t->SetMarkerColor(3);
  gr_inc_t->Draw("pe");

  // interaction histogram
  TCanvas *c8 = new TCanvas("c8","c8");
  gr_int->SetTitle("");
  gr_int->SetLineWidth(2);
  gr_int->GetXaxis()->SetTitle("Slice ID");
  gr_int->GetYaxis()->SetTitle("N_{Int}");
  gr_int->GetXaxis()->SetRangeUser(0, 21);
  gr_int->SetMinimum(0);
  gr_int->Draw("ape");
  gr_int_t->SetLineWidth(2);
  gr_int_t->SetLineColor(3);
  gr_int_t->SetMarkerColor(3);
  gr_int_t->Draw("pe");
  
  // initial histogram
  TCanvas *c9 = new TCanvas("c9","c9");
  gr_ini->SetTitle("");
  gr_ini->SetLineWidth(2);
  gr_ini->GetXaxis()->SetTitle("Slice ID");
  gr_ini->GetYaxis()->SetTitle("N_{Ini}");
  gr_ini->GetXaxis()->SetRangeUser(0, 21);
  gr_ini->SetMinimum(0);
  gr_ini->Draw("ape");
  gr_ini_t->SetLineWidth(2);
  gr_ini_t->SetLineColor(3);
  gr_ini_t->SetMarkerColor(3);
  gr_ini_t->Draw("pe");

  // cross-section
  /*double xs[pi::nthinslices] = {0};
  double err_xs[pi::nthinslices] = {0};
  double incE[pi::nthinslices] = {0};
  double err_incE[pi::nthinslices] = {0};
  TGraph *xs_new = (TGraph*)file->Get("gr_recoxs");
  TGraphErrors *gr_trueincE = (TGraphErrors*)file->Get("gr_trueincE");
  TGraphErrors *gr_truexs = (TGraphErrors*)file->Get("gr_truexs");
  for (int i = 0; i<pi::nthinslices; ++i){
    xs[i] = xs_new->GetPointY(i);
    err_xs[i] = xs_new->GetErrorY(i);
    incE[i] = gr_trueincE->GetPointY(i);
    err_incE[i] = gr_trueincE->GetErrorY(i);
  }
  TGraphErrors *gr_recoxs = new TGraphErrors(pi::nthinslices, incE, xs, 0, err_xs);*/
  
  TFile f2("../../files/exclusive_xsec.root");
  TGraph *total_inel_KE = (TGraph*)f2.Get("total_inel_KE");
  TGraphErrors *gr_truexs = (TGraphErrors*)file->Get("gr_truexs");
  TGraphErrors *gr_recoxs = (TGraphErrors*)file->Get("gr_recoxs");

  TCanvas *c10 = new TCanvas("c10", "c10", 1200, 500);
  gr_recoxs->SetTitle("Pion Inelastic Cross Section");
  gr_recoxs->GetXaxis()->SetTitle("Pion Kinetic Energy (MeV)");
  gr_recoxs->GetXaxis()->SetRangeUser(0, 1000);
  gr_recoxs->GetYaxis()->SetTitle("#sigma_{inelastic} (mb)");
  gr_recoxs->GetYaxis()->SetRangeUser(0, 1000);
  gr_recoxs->SetLineWidth(2);
  gr_recoxs->Draw("ape");
  gr_truexs->SetMarkerColor(3);
  gr_truexs->SetLineColor(3);
  gr_truexs->Draw("pe");
  total_inel_KE->SetLineColor(2);
  
  /*for (int i=0;i<total_inel_KE->GetN();i++) {
    if (total_inel_KE->GetX()[i] <= 476.44931) {
      total_inel_KE->GetY()[i] *= 1;
    }
    else total_inel_KE->GetY()[i] *= 1;
  }*/
  double chi2 = 0;
  int nbins = 0;
  cout<<"KE\tData XS\t\tData XS_err\t\tMC XS\t\tMC XS_err\t\tChi2"<<endl;
  for (int i=1; i<19; ++i) {
    double KE = 975 - 50*i;
    double xs_curve = total_inel_KE->Eval(KE);
    double xs_MC = gr_truexs->GetPointY(i);
    double xserr_MC = gr_truexs->GetErrorY(i);
    double xs_data = gr_recoxs->GetPointY(i);
    double xserr_data = gr_recoxs->GetErrorY(i);
    double c2 = 0;
    if (true) {
      c2 = pow( (xs_data-xs_MC)/sqrt(pow(xserr_data,2)+pow(xserr_MC,2)) , 2);
      chi2 += c2;
      ++nbins;
    }
    cout<<KE<<"\t"<<xs_data<<"\t\t"<<xserr_data<<"\t\t\t"<<xs_MC<<"\t\t"<<xserr_MC<<"\t\t\t"<<c2<<endl;
  }
  cout<<"Chi2/Ndf = "<<chi2/nbins<<endl;
  
  total_inel_KE->Draw("L");
  TLegend *leg10 = new TLegend(0.5,0.15,0.85,0.4);
  leg10->SetFillStyle(0);
  leg10->AddEntry(gr_recoxs, "MC with reconstruction", "pe");
  leg10->AddEntry(gr_truexs, "MC truth", "pe");
  //leg10->AddEntry(total_inel_KE, "Geant4 (theory prediction)", "l");
  leg10->Draw();

  gSystem->Exec("rm -rf plots");
  gSystem->Exec("mkdir plots");
  c1->Print("plots/xs_sliceidinc_reco.pdf");
  c2->Print("plots/xs_sliceidinc_true.pdf");
  c3->Print("plots/xs_sliceidint_reco.pdf");
  c4->Print("plots/xs_sliceidint_true.pdf");
  c5->Print("plots/xs_sliceidini_reco.pdf");
  c6->Print("plots/xs_sliceidini_true.pdf");
  c7->Print("plots/xs_Ninc.pdf");
  c8->Print("plots/xs_Nint.pdf");
  c9->Print("plots/xs_Nini.pdf");
  c10->Print("plots/xs_pi+inel.pdf");

  c1->Print("plots/xs_sliceidinc_reco.png");
  c2->Print("plots/xs_sliceidinc_true.png");
  c3->Print("plots/xs_sliceidint_reco.png");
  c4->Print("plots/xs_sliceidint_true.png");
  c5->Print("plots/xs_sliceidini_reco.png");
  c6->Print("plots/xs_sliceidini_true.png");
  c7->Print("plots/xs_Ninc.png");
  c8->Print("plots/xs_Nint.png");
  c9->Print("plots/xs_Nini.png");
  c10->Print("plots/xs_pi+inel.png");
}
  
