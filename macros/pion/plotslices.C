#include "../../SliceParams.h"
void plotslices(){


  gStyle->SetOptStat(0);
  gStyle->SetPadLeftMargin(0.05);
  gStyle->SetPadRightMargin(0.05);


  TCanvas *c1 = new TCanvas("c1","c1",800,500);

  TH1D *h = new TH1D("h",";z (cm)", (pi::nthinslices + 2), -10, (pi::nthinslices+1)*10);
  for (int i = 1; i<=h->GetNbinsX(); ++i){
    h->SetBinContent(i,1);
  }
  h->SetFillColor(38);
  h->SetBarWidth(.9);
  h->Draw("b");
  h->GetYaxis()->SetRangeUser(0, 1.);
  h->GetYaxis()->SetLabelSize(0);
  h->GetYaxis()->SetTickLength(0);
  h->GetXaxis()->CenterTitle();
  gPad->RedrawAxis();
  
  TLatex l;
  l.SetTextAngle(90);
  l.SetTextColor(0);
  for (int i = 1; i<=h->GetNbinsX(); ++i){
    l.DrawLatex(h->GetBinCenter(i)+2, 0.1, Form("Slice %d",i-2));
  }

  TLine l1(-10, 0.8, 105, 0.7);
  l1.SetLineStyle(2);
  l1.SetLineWidth(4);
  l1.SetLineColor(kOrange);
  l1.Draw();

  TLine l2(105, 0.7, 155, 0.9);
  l2.SetLineStyle(2);
  l2.SetLineWidth(4);
  l2.SetLineColor(kOrange);
  l2.Draw();

  TLine l3(105, 0.7, 125, 0.5);
  l3.SetLineStyle(2);
  l3.SetLineWidth(4);
  l3.SetLineColor(kOrange);
  l3.Draw();

  TMarker m(105, 0.7, 29);
  m.SetMarkerSize(2);
  m.SetMarkerColor(kOrange);
  m.Draw();

  c1->Print("plots/slices.pdf");
}
