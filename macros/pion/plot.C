#include "../../EventType.h"
#include "../../EventSelection.h"
#include "../../SliceParams.h"
#include "TH1D.h"
#include "TH2D.h"
#include "THStack.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TError.h"
#include <iostream>

using namespace std;

int nc = 0;
double totaldata = 0;
double totalmc = 0;
double totalcosmics = 0;

int colors[pi::nIntTypes] = {2, 3, 800, 7, 33, 9, 46, 28, 41};

TFile *fmc;
TFile *fdata;
TFile *fcosmics;

void PrintEvents(string name){

  TH1D *hdata[pi::nCuts];
  TH1D *hmc[pi::nCuts][pi::nIntTypes+1];
  for (int i = 0; i < pi::nCuts; ++i){
    for (int j = 0; j < pi::nIntTypes+1; ++j){
      hmc[i][j] = (TH1D*)fmc->Get(Form("%s_%d_%d",name.c_str(),i,j));
    }
    hdata[i] = (TH1D*)fdata->Get(Form("%s_%d_%d",name.c_str(),i,0));
  }
  
  for (int i = 0; i<pi::nCuts; ++i){
    cout<<"=========="<<pi::cutName[i]<<"=========="<<endl;
    double totalmc = 0;
    for (int j =0; j<pi::nIntTypes+1; ++j){
      if (j==0){
        cout<<pi::intTypeName[j]<<" "<<hdata[i]->Integral()<<endl;
      }
      else{
        cout<<pi::intTypeName[j]<<" "<<hmc[i][j]->Integral()<<endl;
        totalmc += hmc[i][j]->Integral();
      }
    }
    cout<<"Total MC: "<<totalmc<<endl;
  }
}

void plot1d(string name, int cut, string xtitle, string ytitle){

  static bool first = true;

  TH1D *hdata[pi::nCuts];
  TH1D *hcosmics[pi::nCuts];
  TH1D *hmc[pi::nCuts][pi::nIntTypes+1];
  for (int i = 0; i < pi::nCuts; ++i){
    for (int j = 0; j < pi::nIntTypes+1; ++j){
      hmc[i][j] = (TH1D*)fmc->Get(Form("%s_%d_%d",name.c_str(),i,j));
      if (i==0 && j!=0 && first){
        totalmc += hmc[i][j]->Integral();
      }
    }
    hdata[i] = (TH1D*)fdata->Get(Form("%s_%d_%d",name.c_str(),i,0));
    hcosmics[i] = (TH1D*)fcosmics->Get(Form("%s_%d_%d",name.c_str(),i,0));
    if (i==0 && first){
      totaldata = hdata[i]->Integral();
      totalcosmics = hcosmics[i]->Integral();
    }
  }
  //std::cout<<"$$$$$"<<totaldata<<"\t"<<totalmc<<"\t"<<totalcosmics<<std::endl;
  first = false;
  
  TH1D *htotmc;
  THStack *hs = new THStack(Form("hs_%d",nc),"");
  for (int i = 0; i<pi::nIntTypes; ++i){
    hmc[cut][i+1]->SetFillColor(colors[i]);
    hmc[cut][i+1]->SetLineWidth(0);
    hmc[cut][i+1]->Scale(totaldata/totalmc*100150./105575.);
    hs->Add(hmc[cut][i+1]);
    if (i==0) htotmc = (TH1D*)hmc[cut][1]->Clone(Form("htotmc_%d",nc));
    else htotmc->Add(hmc[cut][i+1]);
  }

  if (!htotmc->Integral() || !hdata[cut]->Integral()) return;
  hcosmics[cut]->Scale(totaldata/totalcosmics);

  TCanvas *can = new TCanvas(Form("can_%d",nc), Form("can_%d",nc), 800, 800);
  TPad *pad1 = new TPad(Form("pad1_%d",nc), Form("pad1_%d",nc), 0, 0.2, 1, 1.);
  pad1->SetBottomMargin(0.03);
//  pad1->SetGridx();
//  pad1->SetGridy();
  pad1->Draw();
  pad1->cd();
  

  double max = TMath::Max(hs->GetMaximum(), hdata[cut]->GetMaximum());
  max = TMath::Max(max, hcosmics[cut]->GetMaximum());
  hs->SetMaximum(8000);
  hs->Draw("hist");
  hs->SetTitle(hdata[cut]->GetTitle());
  hs->GetXaxis()->SetTitle(xtitle.c_str());
  hs->GetXaxis()->SetTitleSize(0);
  hs->GetXaxis()->SetLabelSize(0);
  hs->GetYaxis()->SetTitle(ytitle.c_str());
  hs->GetYaxis()->SetTitleSize(0.04);
  hs->GetYaxis()->SetTitleOffset(1.35);
  hs->GetYaxis()->SetLabelSize(0.04);
  hdata[cut]->SetMarkerStyle(20);
  hdata[cut]->Draw("same");
  hcosmics[cut]->SetLineStyle(2);
  hcosmics[cut]->Draw("hist same");
  TLegend *leg = new TLegend(0.15, 0.75, 0.88, 0.9);
//  if (hmc[cut][0]->GetMaximumBin() < hmc[cut][0]->GetNbinsX()/2){
//    leg = new TLegend(0.55,0.7,0.9,0.9);
//  }
//  else{
//    leg = new TLegend(0.15,0.7,0.5,0.9);
//  }
  leg->SetFillStyle(0);
  leg->SetNColumns(3);
  leg->AddEntry(hdata[cut],Form("%s %.0f",pi::intTypeName[0],hdata[cut]->Integral()),"ple");
  leg->AddEntry((TObject*)0,Form("TotalMC %.0f",htotmc->Integral()),"");
  for (int i = 0; i<pi::nIntTypes; ++i){
    leg->AddEntry(hmc[cut][i+1], Form("%s %.0f",pi::intTypeName[i+1],hmc[cut][i+1]->Integral()), "f");
  }
  leg->AddEntry(hcosmics[cut], Form("%s %.0f","Cosmics",hcosmics[cut]->Integral()),"l");
  leg->Draw();

  can->cd();
  TPad *pad2 = new TPad(Form("pad2_%d",nc), Form("pad2_%d",nc), 0, 0, 1, 0.2);
  pad2->SetTopMargin(0.02);
  pad2->SetBottomMargin(0.3);
  pad2->SetGridx();
  pad2->SetGridy();
  pad2->Draw();
  pad2->cd();
  
  TH1D *hratio = (TH1D*)hdata[cut]->Clone(Form("hratio_%d",nc));
  hratio->SetTitle("");
  hratio->Divide(htotmc);
  hratio->GetYaxis()->SetTitle("Data/MC");
  hratio->GetXaxis()->SetLabelSize(0.15);
  hratio->GetXaxis()->SetTitleSize(0.15);
  hratio->GetXaxis()->SetTitleOffset(1.);
  hratio->GetYaxis()->SetLabelSize(0.15);
  hratio->GetYaxis()->SetTitleSize(0.15);
  hratio->GetYaxis()->SetTitleOffset(.3);
//  if (hratio->GetMaximum()>5){
//    hratio->GetYaxis()->SetRangeUser(0,hratio->GetMaximum()*1.1);
//  }
//  else{
  hratio->GetYaxis()->SetRangeUser(0,2);
    //  }
  hratio->GetYaxis()->SetNdivisions(505);
  hratio->Draw();
  
  /*TH1D *hpull = (TH1D*)hdata[cut]->Clone(Form("hpull_%d",nc));
  hpull->SetTitle("");
  for (int ii=1; ii<hdata[cut]->GetNbinsX(); ii++) {
    double dataerr = hdata[cut]->GetBinError(ii);
    if (dataerr == 0) dataerr = 1;
    double pull = (hdata[cut]->GetBinContent(ii)-htotmc->GetBinContent(ii)) / dataerr;
    hpull->SetBinContent(ii, pull); // why there is error even I set it as constant
  }
  hpull->GetYaxis()->SetTitle("Pull");
  hpull->GetXaxis()->SetLabelSize(0.15);
  hpull->GetXaxis()->SetTitleSize(0.15);
  hpull->GetXaxis()->SetTitleOffset(1.);
  hpull->GetYaxis()->SetLabelSize(0.15);
  hpull->GetYaxis()->SetTitleSize(0.15);
  hpull->GetYaxis()->SetTitleOffset(.3);
  double upp = 5;
  double low = -5;
  //if (upp < hpull->GetMaximum()) upp = hpull->GetMaximum();
  //if (low > hpull->GetMinimum()) low = hpull->GetMinimum();
  hpull->GetYaxis()->SetRangeUser(low, upp);
  hpull->GetYaxis()->SetNdivisions(505);
  hpull->Draw();*/

  can->Print(Form("../plots/can_%s_%d_%s.png", name.c_str(), cut, pi::cutName[cut]));
  can->Print(Form("../plots/can_%s_%d_%s.pdf", name.c_str(), cut, pi::cutName[cut]));
  hs->SetMaximum(8*max);
  hs->SetMinimum(1);
  pad1->SetLogy();
  can->Print(Form("../plots/canlog_%s_%d_%s.png", name.c_str(), cut, pi::cutName[cut]));
  can->Print(Form("../plots/canlog_%s_%d_%s.pdf", name.c_str(), cut, pi::cutName[cut]));

  ++nc;
}

void plot1dslice(string name, int cut, string xtitle, string ytitle){

  TH1D *hdata[pi::nthinslices][pi::nCuts];
  TH1D *hmc[pi::nthinslices][pi::nCuts][pi::nIntTypes+1];
  for (int i = 0; i < pi::nCuts; ++i){
    for (int j = 0; j < pi::nIntTypes+1; ++j){
      for (int k = 0; k < pi::nthinslices; ++k){
        hmc[k][i][j] = (TH1D*)fmc->Get(Form("%s_%d_%d_%d",name.c_str(),k,i,j));
        if (j==0) hdata[k][i] = (TH1D*)fdata->Get(Form("%s_%d_%d_%d",name.c_str(),k,i,0));
      }
    }
  }
  //std::cout<<totaldata<<" "<<totalmc<<std::endl;

  for (int j = 0; j<pi::nthinslices; ++j){
    TCanvas *can = new TCanvas(Form("can_%d",nc), Form("can_%d",nc), 800, 800);
    TPad *pad1 = new TPad(Form("pad1_%d",nc), Form("pad1_%d",nc), 0, 0.2, 1, 1.);
    pad1->SetBottomMargin(0.03);
    pad1->SetGridx();
    pad1->SetGridy();
    pad1->Draw();
    pad1->cd();
    
    TH1D *htotmc;
    THStack *hs = new THStack(Form("hs_%d",nc),"");
    for (int i = 0; i<pi::nIntTypes; ++i){
      hmc[j][cut][i+1]->SetFillColor(colors[i]);
      hmc[j][cut][i+1]->SetLineWidth(0);
      hmc[j][cut][i+1]->Scale(totaldata/totalmc);
      hs->Add(hmc[j][cut][i+1]);
      if (i==0) htotmc = (TH1D*)hmc[j][cut][1]->Clone(Form("htotmc_%d",nc));
      else htotmc->Add(hmc[j][cut][i+1]);
    }
    
    double max = TMath::Max(hs->GetMaximum(), hdata[j][cut]->GetMaximum());
    hs->SetMaximum(1.2*max);
    hs->Draw("hist");
    hs->SetTitle(hdata[j][cut]->GetTitle());
    hs->GetXaxis()->SetTitle(xtitle.c_str());
    hs->GetXaxis()->SetTitleSize(0);
    hs->GetXaxis()->SetLabelSize(0);
    hs->GetYaxis()->SetTitle(ytitle.c_str());
    hs->GetYaxis()->SetTitleSize(0.04);
    hs->GetYaxis()->SetTitleOffset(1.8);
    hs->GetYaxis()->SetLabelSize(0.04);
    hdata[j][cut]->SetMarkerStyle(20);
    hdata[j][cut]->Draw("same");
    TLegend *leg = new TLegend(0.15, 0.75, 0.88, 0.9);
    //  if (hmc[j][cut][0]->GetMaximumBin() < hmc[j][cut][0]->GetNbinsX()/2){
    //    leg = new TLegend(0.55,0.7,0.9,0.9);
    //  }
    //  else{
    //    leg = new TLegend(0.15,0.7,0.5,0.9);
    //  }
    //leg->SetFillStyle(0);
    leg->SetNColumns(3);
    leg->AddEntry(hdata[j][cut],Form("%s %.0f",pi::intTypeName[0],hdata[j][cut]->Integral()),"ple");
    leg->AddEntry((TObject*)0,Form("TotalMC %.0f",htotmc->Integral()),"");
    for (int i = 0; i<pi::nIntTypes; ++i){
      leg->AddEntry(hmc[j][cut][i+1], Form("%s %.0f",pi::intTypeName[i+1],hmc[j][cut][i+1]->Integral()), "f");
    }
    leg->Draw();
    
    can->cd();
    TPad *pad2 = new TPad(Form("pad2_%d",nc), Form("pad2_%d",nc), 0, 0, 1, 0.2);
    pad2->SetTopMargin(0.02);
    pad2->SetBottomMargin(0.25);
    pad2->SetGridx();
    pad2->SetGridy();
    pad2->Draw();
    pad2->cd();
    
    TH1D *hratio = (TH1D*)hdata[j][cut]->Clone(Form("hratio_%d",nc));
    hratio->SetTitle("");
    hratio->Divide(htotmc);
    hratio->GetYaxis()->SetTitle("Data/MC");
    hratio->GetXaxis()->SetLabelSize(0.12);
    hratio->GetXaxis()->SetTitleSize(0.12);
    hratio->GetXaxis()->SetTitleOffset(1.);
    hratio->GetYaxis()->SetLabelSize(0.12);
    hratio->GetYaxis()->SetTitleSize(0.12);
    hratio->GetYaxis()->SetTitleOffset(.5);
    hratio->GetYaxis()->SetRangeUser(0,5);
    hratio->GetYaxis()->SetNdivisions(505);
    hratio->Draw();
    
    /*TH1D *hpull = (TH1D*)hdata[j][cut]->Clone(Form("hpull_%d",nc));
    hpull->SetTitle("");
    for (int ii=1; ii<hdata[j][cut]->GetNbinsX(); ii++) {
      double dataerr = hdata[j][cut]->GetBinError(ii);
      if (dataerr == 0) dataerr = 1;
      double pull = (hdata[j][cut]->GetBinContent(ii)-htotmc->GetBinContent(ii)) / dataerr;
      hpull->SetBinContent(ii, pull); // why there is error even I set it as constant
    }
    hpull->GetYaxis()->SetTitle("Pull");
    hpull->GetXaxis()->SetLabelSize(0.12);
    hpull->GetXaxis()->SetTitleSize(0.12);
    hpull->GetXaxis()->SetTitleOffset(1.);
    hpull->GetYaxis()->SetLabelSize(0.12);
    hpull->GetYaxis()->SetTitleSize(0.12);
    hpull->GetYaxis()->SetTitleOffset(.3);
    double upp = 5;
    double low = -5;
    //if (upp < hpull->GetMaximum()) upp = hpull->GetMaximum();
    //if (low > hpull->GetMinimum()) low = hpull->GetMinimum();
    hpull->GetYaxis()->SetRangeUser(low, upp);
    hpull->GetYaxis()->SetNdivisions(505);
    hpull->Draw();*/
    
    can->Print(Form("../plots/can_%s_%d_%s_%02d.png", name.c_str(), cut, pi::cutName[cut], j));
    can->Print(Form("../plots/can_%s_%d_%s_%02d.pdf", name.c_str(), cut, pi::cutName[cut], j));
    hs->SetMaximum(8*max);
    hs->SetMinimum(1);
    pad1->SetLogy();
    can->Print(Form("../plots/canlog_%s_%d_%s_%02d.png", name.c_str(), cut, pi::cutName[cut], j));
    can->Print(Form("../plots/canlog_%s_%d_%s_%02d.pdf", name.c_str(), cut, pi::cutName[cut], j));
    
    ++nc;
    
  }
}


void plot2d(string name, int cut){

  TH2D *h[pi::nCuts][pi::nIntTypes+1];
  for (int i = 0; i < pi::nCuts; ++i){
    for (int j = 0; j < pi::nIntTypes+1; ++j){
      h[i][j] = (TH2D*)fmc->Get(Form("%s_%d_%d",name.c_str(),i,j));
    }
  }
  for (int i = 1; i<=4; ++i){
    TCanvas *can = new TCanvas(Form("can_%d",nc), Form("can_%d",nc));
    h[cut][i]->Draw("colz");
    can->Print(Form("../plots/can_%s_%s_%s.png",name.c_str(), pi::cutName[cut], pi::intTypeName[i]));
    can->Print(Form("../plots/can_%s_%s_%s.pdf",name.c_str(), pi::cutName[cut], pi::intTypeName[i]));
    ++nc;
  }
}
void plot(){
  gSystem->Exec("rm -rf ../plots");
  gSystem->Exec("mkdir ../plots");

  gStyle->SetOptStat(0);
  gErrorIgnoreLevel = kWarning;

  fmc = TFile::Open("../../build/mcprod4a_0721_76sigma.root");
  //fdata = TFile::Open("../install/bin/mcprod4a.root");
  fdata = TFile::Open("../../build/data_0721_nominal.root");
  fcosmics = TFile::Open("../../build/cosmics.root");

  for (int i = 0; i<pi::nCuts; ++i){
    plot1d("hmediandEdx", i, "Median dE/dx (MeV/cm)", "Events");
    plot1d("hChi2_proton", i, "Chi2_p/Ndof", "Events");
    plot1d("hdaughter_michel_score", i, "Daughter Michel Score", "Events");
    plot1d("hdaughter_michel_scoreMu", i, "Muon Daughter Michel Score", "Events");
    plot1d("hdaughter_michel_score2Mu", i, "Muon Daughter Michel Score2", "Events");
    plot1d("hdaughter_michel_scorePi", i, "Pion Daughter Michel Score", "Events");
    plot1d("htrackscore", i, "Track Score", "Events");
    plot1d("hemscore", i, "EM Score", "Events");
    plot1d("hdEdx_5cm", i, "dEdx_5cm (MeV/cm)", "Events");
    plot1d("hreco_beam_endZ_SCE", i, "Reco track end Z (cm)", "Events");
    plot1d("htrue_beam_endZ_SCE", i, "True pion/muon end Z (cm)", "Events");
    plot1d("htrue_sliceID", i, "True Slice ID", "Events");
    plot1d("hreco_sliceID", i, "Reco Slice ID", "Events");
    plot1d("hreco_true_beam_endZ_SCE", i, "Reco - True End Z SCE (cm)", "Events");
    plot2d("hreco_vs_true_beam_endZ_SCE",i);
    plot1d("hdeltax", i, "#Deltax/#sigma_{x}", "Events");
    plot1d("hdeltay", i, "#Deltay/#sigma_{y}", "Events");
    plot1d("hdeltaz", i, "#Deltaz/#sigma_{z}", "Events");
    plot1d("hcostheta", i, "cos#theta", "Events");
    plot1d("hini_recoE", i, "Reco initial energy (MeV)", "Events");
    plot1d("hint_recoE", i, "Reco interacting energy (MeV)", "Events");
    plot1d("hreco_trklen", i, "Reco track length (cm)", "Events");
    plot1d("hreco_beam_startX_SCE", i, "Reco track start X (cm)", "Events");
    plot1d("hreco_beam_startY_SCE", i, "Reco track start Y (cm)", "Events");
    plot1d("hreco_beam_startZ_SCE", i, "Reco track start Z (cm)", "Events");

    plot1d("hreco_beam_dcosX_SCE", i, "dcosX", "Events");
    plot1d("hreco_beam_dcosY_SCE", i, "dcosY", "Events");
    plot1d("hreco_beam_dcosZ_SCE", i, "dcosZ", "Events");

    plot1d("hreco_beam_angleX_SCE", i, "angleX", "Events");
    plot1d("hreco_beam_angleY_SCE", i, "angleY", "Events");
    plot1d("hreco_beam_angleZ_SCE", i, "angleZ", "Events");
  }

  plot1d("hmediandEdx_bkg", 0, "Median dE/dx (MeV/cm)", "Events");
  plot1d("hmediandEdx_bkg", 6, "Median dE/dx (MeV/cm)", "Events");
  plot1d("hChi2_proton_bkg", 0, "Chi2_p/Ndf", "Events");
  plot1d("hChi2_proton_bkg", 6, "Chi2_p/Ndf", "Events");
  plot1d("hdaughter_michel_score_bkg", 0, "Michel score", "Events");
  plot1d("hdaughter_michel_score_bkg", 6, "Michel score", "Events");
  plot1d("hcostheta_bkg", 0, "Cos(theta)", "Events");
  plot1d("hcostheta_bkg", 6, "Cos(theta)", "Events");
  /*plot1dslice("hmediandEdxSlice", 6, "Median dE/dx (MeV/cm)", "Events");
  plot1dslice("hChi2_protonSlice", 6, "Chi2_p/Ndf", "Events");
  plot1dslice("hdaughter_michel_scoreSlice", 6, "Michel score", "Events");
  plot1dslice("hcosthetaSlice", 6, "Cos(theta)", "Events");
*/
  PrintEvents("hreco_beam_endZ_SCE");

}
