#include "SliceParams.h"
#include "EventType.h"
#include "EventSelection.h"
#include "TemplateFitter.h"
#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TGraphErrors.h"
#include "TVectorD.h"
#include "TSystem.h"
#include "json/json.h"
#include <fstream>
#include <iostream>
#include <string>

static void show_usage(std::string name)
{
  std::cerr << "Usage: " << name <<" <option(s)>\n"
            << "Options:\n"
            << "\t-h,--help\t\tShow this help message.\n"
            << "\t-c config.json\t\tSpecify configuration file."
            << std::endl;
}

void save_results(vector<double> vslice, vector<double> vcorr, vector<double> vcorrerr, const char* particle, string outfile, double par=1, double parerr=0){
  TCanvas *c1 = new TCanvas("c1","c1");
  c1->SetGrid();
  TGraphErrors *gr_corr = new TGraphErrors(vslice.size(), &vslice[0], &vcorr[0], 0, &vcorrerr[0]);
  gr_corr->SetMarkerStyle(106);
  gr_corr->SetTitle(Form("%s bkg constraint (overall fit %.2f #pm  %.2f)", particle, par, parerr));
  gr_corr->GetXaxis()->SetTitle("Slice");
  gr_corr->GetXaxis()->SetLimits(-1, pi::reco_nbins);
  gr_corr->GetYaxis()->SetTitle("Scale factor");
  gr_corr->GetYaxis()->SetRangeUser(0, 5);
  gr_corr->Draw("AP");
  TLine *line = new TLine(-1, par, pi::reco_nbins, par);
  line->SetLineColor(kRed);
  line->Draw("same");
  TLine *line_upp = new TLine(-1, par+parerr, pi::reco_nbins, par+parerr);
  line_upp->SetLineColor(kRed);
  line_upp->SetLineStyle(2);
  line_upp->Draw("same");
  TLine *line_low = new TLine(-1, par-parerr, pi::reco_nbins, par-parerr);
  line_low->SetLineColor(kRed);
  line_low->SetLineStyle(2);
  line_low->Draw("same");
  c1->Print(Form("%s_%s.png", outfile.substr(0,outfile.find(".root")).c_str(), particle));
  c1->Print(Form("%s_%s.pdf", outfile.substr(0,outfile.find(".root")).c_str(), particle));
  gr_corr->Write(Form("gr_corr_%s", particle));
  TVectorD sf(2);
  sf[0] = par;
  sf[1] = parerr;
  sf.Write(Form("sf_%s", particle));
}

double bkgFit_mu(TFile *fmc, TFile *fdata, string outfile){
  const char varname[50] = "daughter_michel_score";
  const char particle[10] = "mu";
  cout<<"##### Constrain muon bkg using "<<varname<<endl;
  
  double totaldata = 0;
  double totalmc = 0;
  //double ndata[pi::reco_nbins] = {0};
  //double nmc[pi::reco_nbins] = {0};
  TH1D *hvarSlice[pi::reco_nbins][pi::nCuts][pi::nIntTypes+1];
  TH1D *hvar[pi::nCuts][pi::nIntTypes+1];
  //TH1D *hsliceID[pi::nCuts][pi::nIntTypes+1];
 
  for (int i = 0; i < pi::nCuts; ++i){
    for (int j = 0; j < pi::nIntTypes+1; ++j){
      for (int k = 0; k < pi::reco_nbins; ++k){
        if (j==0){ // data
          hvarSlice[k][i][j] = (TH1D*)fdata->Get(Form("h%sSlice_%d_%d_%d",varname,k,i,j));
          if (i==6) totaldata += hvarSlice[k][i][j]->Integral(); // normalization should be done before any cut or after all cuts?
        }
        else{ // MC
          hvarSlice[k][i][j] = (TH1D*)fmc->Get(Form("h%sSlice_%d_%d_%d",varname,k,i,j));
          if (i==6) totalmc += hvarSlice[k][i][j]->Integral();
        }
        if (k==0){
          if (j==0){ // data
            hvar[i][j] = (TH1D*)fdata->Get(Form("h%s_bkg_%d_%d",varname,i,j));
            //hsliceID[i][j] = (TH1D*)fdata->Get(Form("hreco_sliceID_%d_%d",i,j));
          }
          else{ // MC
            hvar[i][j] = (TH1D*)fmc->Get(Form("h%s_bkg_%d_%d",varname,i,j));
            //hsliceID[i][j] = (TH1D*)fmc->Get(Form("hreco_sliceID_%d_%d",i,j));
          }
        }
        if (i==6){
          if (j==0){// data
            //ndata[k] += hsliceID[i][j]->GetBinContent(hsliceID[i][j]->FindBin(k+0.5));
          }
          else {//MC
            //nmc[k] += hsliceID[i][j]->GetBinContent(hsliceID[i][j]->FindBin(k+0.5));
          }
        }
      }
    }
  }
  std::cout<<"#Data "<<totaldata<<"; #MC "<<totalmc<<std::endl;
  TemplateFitter fitter;

  std::vector<double> vslice;
  std::vector<double> vcorr;
  std::vector<double> vcorrerr;

  // muon constraint begins
  for (int i = 0; i<pi::reco_nbins; ++i){
    std::cout<<"##### Slice "<<i<<std::endl;
    TH1D *h0 = hvarSlice[i][6][pi::kData];
    TH1D *h1 = hvarSlice[i][6][pi::kPiInel];
    h1->Add(hvarSlice[i][6][pi::kPiElas]);
    h1->Add(hvarSlice[i][6][pi::kMIDp]);
    h1->Add(hvarSlice[i][6][pi::kMIDcosmic]);
    h1->Add(hvarSlice[i][6][pi::kMIDpi]);
    h1->Add(hvarSlice[i][6][pi::kMIDeg]);
    h1->Add(hvarSlice[i][6][pi::kMIDother]);
    // components to be rescaled
    TH1D *h2 = hvarSlice[i][6][pi::kMuon];
    h2->Add(hvarSlice[i][6][pi::kMIDmu]);
    
    h1->Scale(totaldata/totalmc);
    h2->Scale(totaldata/totalmc);
//    h1->Scale(ndata[i]/nmc[i]);
//    h2->Scale(ndata[i]/nmc[i]);

    fitter.SetHistograms(h0, h1, h2);
    fitter.SetFitRange(7, 9);
    fitter.Fit();
    if (fitter.GetFitStatus()){
      vslice.push_back(i);
      vcorr.push_back(fitter.GetPar());
      vcorrerr.push_back(fitter.GetParError());
    }
  }
  cout<<"##### Global fit #####"<<endl;
  TH1D *h0 = hvar[6][pi::kData];
  TH1D *h1 = hvar[6][pi::kPiInel];
  h1->Add(hvar[6][pi::kPiElas]);
  h1->Add(hvar[6][pi::kMIDp]);
  h1->Add(hvar[6][pi::kMIDcosmic]);
  h1->Add(hvar[6][pi::kMIDpi]);
  h1->Add(hvar[6][pi::kMIDeg]);
  h1->Add(hvar[6][pi::kMIDother]);
  // components to be rescaled
  TH1D *h2 = hvar[6][pi::kMuon];
  h2->Add(hvar[6][pi::kMIDmu]);
  
  h1->Scale(totaldata/totalmc);
  h2->Scale(totaldata/totalmc);
  fitter.SetHistograms(h0, h1, h2);
  fitter.SetFitRange(h0->FindBin(0.6), h0->FindBin(0.9));
  fitter.Fit();
  double par = fitter.GetPar();
  double parerr = fitter.GetParError();
  std::cout<<par<<" "<<parerr<<std::endl;
  
  save_results(vslice, vcorr, vcorrerr, particle, outfile, par, parerr);
  return fitter.GetPar();
}

double bkgFit_p(TFile *fmc, TFile *fdata, string outfile, double muscale = 1.){
  const char varname[50] = "Chi2_proton"; // "Chi2_proton"/"mediandEdx"
  const char particle[10] = "p";
  cout<<"##### Constrain proton bkg using "<<varname<<endl;
  
  double totaldata = 0;
  double totalmc = 0;
  TH1D *hvarSlice[pi::reco_nbins][pi::nCuts][pi::nIntTypes+1];
  TH1D *hvar[pi::nCuts][pi::nIntTypes+1];

  for (int i = 0; i < pi::nCuts; ++i){
    for (int j = 0; j < pi::nIntTypes+1; ++j){
      for (int k = 0; k < pi::reco_nbins; ++k){
        if (j==0){ // data
          hvarSlice[k][i][j] = (TH1D*)fdata->Get(Form("h%sSlice_%d_%d_%d",varname,k,i,j));
          if (i==6) totaldata += hvarSlice[k][i][j]->Integral(); // normalization should be done before any cut or after all cuts?
        }
        else{ // MC
          hvarSlice[k][i][j] = (TH1D*)fmc->Get(Form("h%sSlice_%d_%d_%d",varname,k,i,j));
          if (i==6) totalmc += hvarSlice[k][i][j]->Integral();
        }
        if (k==0){
          if (j==0){ // data
            hvar[i][j] = (TH1D*)fdata->Get(Form("h%s_bkg_%d_%d",varname,i,j));
          }
          else{ // MC
            hvar[i][j] = (TH1D*)fmc->Get(Form("h%s_bkg_%d_%d",varname,i,j));
          }
        }
      }
    }
  }
  std::cout<<"#Data "<<totaldata<<"; #MC "<<totalmc<<std::endl;
  TemplateFitter fitter;

  std::vector<double> vslice;
  std::vector<double> vcorr;
  std::vector<double> vcorrerr;

  // proton constraint begins
  TH1D *htemp;
  for (int i = 0; i<pi::reco_nbins; ++i){
    std::cout<<"##### Slice "<<i<<std::endl;
    TH1D *h0 = hvarSlice[i][6][pi::kData];
    TH1D *h1 = hvarSlice[i][6][pi::kPiInel];
    h1->Add(hvarSlice[i][6][pi::kPiElas]);
    htemp = (TH1D*)hvarSlice[i][6][pi::kMuon]->Clone();
    htemp->Add(hvarSlice[i][6][pi::kMIDmu]);
    htemp->Scale(muscale);
    h1->Add(htemp);
    h1->Add(hvarSlice[i][6][pi::kMIDcosmic]);
    h1->Add(hvarSlice[i][6][pi::kMIDpi]);
    h1->Add(hvarSlice[i][6][pi::kMIDeg]);
    h1->Add(hvarSlice[i][6][pi::kMIDother]);
    // components to be rescaled
    TH1D *h2 = hvarSlice[i][6][pi::kMIDp];
    
    h1->Scale(totaldata/totalmc);
    h2->Scale(totaldata/totalmc);
    fitter.SetHistograms(h0, h1, h2);
    fitter.SetFitRange(3, 7);
    fitter.Fit();
    if (fitter.GetFitStatus()){
      vslice.push_back(i);
      vcorr.push_back(fitter.GetPar());
      vcorrerr.push_back(fitter.GetParError());
    }
  }
  cout<<"##### Global fit #####"<<endl;
  TH1D *h0 = hvar[6][pi::kData];
  TH1D *h1 = hvar[6][pi::kPiInel];
  h1->Add(hvar[6][pi::kPiElas]);
  htemp = (TH1D*)hvar[6][pi::kMuon]->Clone();
  htemp->Add(hvar[6][pi::kMIDmu]);
  htemp->Scale(muscale);
  h1->Add(htemp);
  h1->Add(hvar[6][pi::kMIDcosmic]);
  h1->Add(hvar[6][pi::kMIDpi]);
  h1->Add(hvar[6][pi::kMIDeg]);
  h1->Add(hvar[6][pi::kMIDother]);
  
  TH1D *h2 = hvar[6][pi::kMIDp];
  h1->Scale(totaldata/totalmc);
  h2->Scale(totaldata/totalmc);
  fitter.SetHistograms(h0, h1, h2);
  fitter.SetFitRange(h0->FindBin(20), h0->FindBin(70));
  fitter.Fit();
  double par = fitter.GetPar();
  double parerr = fitter.GetParError();
  std::cout<<par<<" "<<parerr<<std::endl;
  
  save_results(vslice, vcorr, vcorrerr, particle, outfile, par, parerr);
  return fitter.GetPar();
}

double bkgFit_spi(TFile *fmc, TFile *fdata, string outfile, double muscale = 1., double pscale = 1.){
  const char varname[50] = "costheta";
  const char particle[10] = "spi";
  cout<<"##### Constrain secondary pion bkg using "<<varname<<endl;
  
  double totaldata = 0;
  double totalmc = 0;
  TH1D *hvarSlice[pi::reco_nbins][pi::nCuts][pi::nIntTypes+1];
  TH1D *hvar[pi::nCuts][pi::nIntTypes+1];

  for (int i = 0; i < pi::nCuts; ++i){
    for (int j = 0; j < pi::nIntTypes+1; ++j){
      for (int k = 0; k < pi::reco_nbins; ++k){
        if (j==0){ // data
          hvarSlice[k][i][j] = (TH1D*)fdata->Get(Form("h%sSlice_%d_%d_%d",varname,k,i,j));
          if (i==6) totaldata += hvarSlice[k][i][j]->Integral(); // normalization should be done before any cut or after all cuts?
        }
        else{ // MC
          hvarSlice[k][i][j] = (TH1D*)fmc->Get(Form("h%sSlice_%d_%d_%d",varname,k,i,j));
          if (i==6) totalmc += hvarSlice[k][i][j]->Integral();
        }
        if (k==0){
          if (j==0){ // data
            hvar[i][j] = (TH1D*)fdata->Get(Form("h%s_bkg_%d_%d",varname,i,j));
          }
          else{ // MC
            hvar[i][j] = (TH1D*)fmc->Get(Form("h%s_bkg_%d_%d",varname,i,j));
          }
        }
      }
    }
  }
  std::cout<<"#Data "<<totaldata<<"; #MC "<<totalmc<<std::endl;
  TemplateFitter fitter;

  std::vector<double> vslice;
  std::vector<double> vcorr;
  std::vector<double> vcorrerr;

  // secondary pion constraint begins
  TH1D *htemp;
  for (int i = 0; i<pi::reco_nbins; ++i){
    std::cout<<"##### Slice "<<i<<std::endl;
    TH1D *h0 = hvarSlice[i][6][pi::kData];
    TH1D *h1 = hvarSlice[i][6][pi::kPiInel];
    h1->Add(hvarSlice[i][6][pi::kPiElas]);
    htemp = (TH1D*)hvarSlice[i][6][pi::kMuon]->Clone();
    htemp->Add(hvarSlice[i][6][pi::kMIDmu]);
    htemp->Scale(muscale);
    h1->Add(htemp);
    htemp = (TH1D*)hvarSlice[i][6][pi::kMIDp]->Clone();
    htemp->Scale(pscale);
    h1->Add(htemp);
    h1->Add(hvarSlice[i][6][pi::kMIDcosmic]);
    h1->Add(hvarSlice[i][6][pi::kMIDeg]);
    h1->Add(hvarSlice[i][6][pi::kMIDother]);
    // components to be rescaled
    TH1D *h2 = hvarSlice[i][6][pi::kMIDpi];
    
    h1->Scale(totaldata/totalmc);
    h2->Scale(totaldata/totalmc);
    fitter.SetHistograms(h0, h1, h2);
    fitter.SetFitRange(6, 10); // don't use underflow bin
    fitter.Fit();
    if (fitter.GetFitStatus()){
      vslice.push_back(i);
      vcorr.push_back(fitter.GetPar());
      vcorrerr.push_back(fitter.GetParError());
    }
  }
  cout<<"##### Global fit #####"<<endl;
  TH1D *h0 = hvar[6][pi::kData];
  TH1D *h1 = hvar[6][pi::kPiInel];
  h1->Add(hvar[6][pi::kPiElas]);
  htemp = (TH1D*)hvar[6][pi::kMuon]->Clone();
  htemp->Add(hvar[6][pi::kMIDmu]);
  htemp->Scale(muscale);
  h1->Add(htemp);
  htemp = (TH1D*)hvar[6][pi::kMIDp]->Clone();
  htemp->Scale(pscale);
  h1->Add(htemp);
  h1->Add(hvar[6][pi::kMIDcosmic]);
  h1->Add(hvar[6][pi::kMIDeg]);
  h1->Add(hvar[6][pi::kMIDother]);
  
  TH1D *h2 = hvar[6][pi::kMIDpi];
  h1->Scale(totaldata/totalmc);
  h2->Scale(totaldata/totalmc);
  fitter.SetHistograms(h0, h1, h2);
  fitter.SetFitRange(h0->FindBin(0.9), h0->FindBin(0.95));
  fitter.Fit();
  double par = fitter.GetPar();
  double parerr = fitter.GetParError();
  std::cout<<par<<" "<<parerr<<std::endl;
  
  save_results(vslice, vcorr, vcorrerr, particle, outfile, par, parerr);
  return fitter.GetPar();
}


int main(int argc, char** argv){

  //bool fitfakedata = false;

  bool found_config = false;

  string config_file;

  //for (int i = 1; i < argc; ++i) {
  for (int iArg = 1; iArg < argc; iArg++) {
    if (!strcasecmp(argv[iArg],"-c")) {
     config_file = argv[++iArg];
     found_config = true;
    }
    if (!strcasecmp(argv[iArg],"-h")) {
      show_usage(argv[0]);
      return 1;
    }
  }

  if (!found_config){
    show_usage(argv[0]);
    return 1;
  }

  Json::Value root;
  ifstream file(config_file);
  file >> root;
  cout<<root<<endl;
  
  TFile *fdata = TFile::Open(root["datafile"].asString().c_str());
  TFile *fmc = TFile::Open(root["mcfile"].asString().c_str());
  TFile *fout = TFile::Open(root["outfile"].asString().c_str(), "recreate");

  double muscale = bkgFit_mu(fmc, fdata, root["outfile"].asString());
  double pscale = bkgFit_p(fmc, fdata, root["outfile"].asString(), muscale);
  double spiscale = bkgFit_spi(fmc, fdata, root["outfile"].asString(), muscale, pscale);

  fout->Close();
  return 0;
}
