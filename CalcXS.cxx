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
#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"
#include "TH2D.h"
#include "BetheBloch.h"

static void show_usage(std::string name)
{
  std::cerr << "Usage: " << name <<" <option(s)>\n"
            << "Options:\n"
            << "\t-h,--help\t\tShow this help message.\n"
            << "\t-c config.json\t\tSpecify configuration file."
            << std::endl;
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
  TFile *fmc   = TFile::Open(root["mcfile"].asString().c_str());
  TFile *fbkg  = TFile::Open(root["bkgfile"].asString().c_str());
  TFile *fout  = TFile::Open(root["outfile"].asString().c_str(), "recreate");

  // background constraint
  //TVectorD sf(2); sf[0] = 1.; sf[1] = 0.;
  TVectorD *sf_mu = (TVectorD*)fbkg->Get("sf_mu");
  TVectorD *sf_p = (TVectorD*)fbkg->Get("sf_p");
  TVectorD *sf_spi = (TVectorD*)fbkg->Get("sf_spi");

  cout<<"Muon scaling factor: "<<(*sf_mu)[0]<<"+-"<<(*sf_mu)[1]<<endl;
  cout<<"Proton scaling factor: "<<(*sf_p)[0]<<"+-"<<(*sf_p)[1]<<endl;
  cout<<"Pion scaling factor: "<<(*sf_spi)[0]<<"+-"<<(*sf_spi)[1]<<endl;

  TH1D *hsliceID[pi::nIntTypes+1];
  TH1D *hdata = new TH1D("hdata","Data;Slice ID;Events",pi::nthinslices+2,-1,pi::nthinslices+1); // h_recosliceid_allevts_cuts (hreco_sliceID_6_0)
  TH1D *hproton = new TH1D("hproton","Proton background;Slice ID;Events",pi::nthinslices+2,-1,pi::nthinslices+1);
  TH1D *hmu = new TH1D("hmu","Muon background;Slice ID;Events",pi::nthinslices+2,-1,pi::nthinslices+1);
  TH1D *hspi = new TH1D("hspi","Secondary pion background;Slice ID;Events",pi::nthinslices+2,-1,pi::nthinslices+1);
  TH1D *hpiel = new TH1D("hpiel","Pion elastic;Slice ID;Events",pi::nthinslices+2,-1,pi::nthinslices+1);
  TH1D *hother = new TH1D("hother","Other backgrounds;Slice ID;Events",pi::nthinslices+2,-1,pi::nthinslices+1);
  hdata->Sumw2();
  hproton->Sumw2();
  hmu->Sumw2();
  hspi->Sumw2();
  hpiel->Sumw2();
  hother->Sumw2();
  
  TH1D *hincsliceID[pi::nIntTypes+1];
  TH1D *hdata_inc = new TH1D("hdata_inc","Data_inc;Slice ID;Events",pi::nthinslices+2,-1,pi::nthinslices+1);
  TH1D *hproton_inc = new TH1D("hproton_inc","Proton_inc background;Slice ID;Events",pi::nthinslices+2,-1,pi::nthinslices+1);
  TH1D *hmu_inc = new TH1D("hmu_inc","Muon_inc background;Slice ID;Events",pi::nthinslices+2,-1,pi::nthinslices+1);
  TH1D *hspi_inc = new TH1D("hspi_inc","Secondary pion_inc background;Slice ID;Events",pi::nthinslices+2,-1,pi::nthinslices+1);
  //TH1D *hpiel_inc = new TH1D("hpiel_inc","Pion elastic_inc;Slice ID;Events",pi::nthinslices+2,-1,pi::nthinslices+1);
  TH1D *hother_inc = new TH1D("hother_inc","Other_inc backgrounds;Slice ID;Events",pi::nthinslices+2,-1,pi::nthinslices+1);
  hdata_inc->Sumw2();
  hproton_inc->Sumw2();
  hmu_inc->Sumw2();
  hspi_inc->Sumw2();
  //hpiel_inc->Sumw2();
  hother_inc->Sumw2();
  
  TH1D *hinisliceID[pi::nIntTypes+1];
  TH1D *hdata_ini = new TH1D("hdata_ini","Data_ini;Slice ID;Events",pi::nthinslices+2,-1,pi::nthinslices+1);
  TH1D *hproton_ini = new TH1D("hproton_ini","Proton_ini background;Slice ID;Events",pi::nthinslices+2,-1,pi::nthinslices+1);
  TH1D *hmu_ini = new TH1D("hmu_ini","Muon_ini background;Slice ID;Events",pi::nthinslices+2,-1,pi::nthinslices+1);
  TH1D *hspi_ini = new TH1D("hspi_ini","Secondary pion_ini background;Slice ID;Events",pi::nthinslices+2,-1,pi::nthinslices+1);
  //TH1D *hpiel_ini = new TH1D("hpiel_ini","Pion elastic_ini;Slice ID;Events",pi::nthinslices+2,-1,pi::nthinslices+1);
  TH1D *hother_ini = new TH1D("hother_ini","Other_ini backgrounds;Slice ID;Events",pi::nthinslices+2,-1,pi::nthinslices+1);
  hdata_ini->Sumw2();
  hproton_ini->Sumw2();
  hmu_ini->Sumw2();
  hspi_ini->Sumw2();
  //hpiel_ini->Sumw2();
  hother_ini->Sumw2();

  // first loop to get total number of events
  for (int i = 0; i < pi::nIntTypes+1; ++i){
    if (i==0){
      hsliceID[i] = (TH1D*)fdata->Get(Form("hreco_sliceID_%d_%d",pi::nCuts-1,i));
      hincsliceID[i] = (TH1D*)fdata->Get(Form("hreco_incsliceID_%d_%d",pi::nCuts-1,i));
      hinisliceID[i] = (TH1D*)fdata->Get(Form("hreco_inisliceID_%d_%d",pi::nCuts-1,i));
    }
    else {
      hsliceID[i] = (TH1D*)fmc->Get(Form("hreco_sliceID_%d_%d",pi::nCuts-1,i));
      hincsliceID[i] = (TH1D*)fmc->Get(Form("hreco_incsliceID_%d_%d",pi::nCuts-1,i));
      hinisliceID[i] = (TH1D*)fmc->Get(Form("hreco_inisliceID_%d_%d",pi::nCuts-1,i));
    }
  }
  TH1D *hmc = new TH1D("hmc","MC;Slice ID;Events",hsliceID[0]->GetNbinsX(),-1,hsliceID[0]->GetNbinsX()-1);
  hmc->Sumw2();
  TH1D *hmc_inc = new TH1D("hmc_inc","MC_inc;Slice ID;Events",hsliceID[0]->GetNbinsX(),-1,hsliceID[0]->GetNbinsX()-1);
  hmc_inc->Sumw2();
  TH1D *hmc_ini = new TH1D("hmc_ini","MC_ini;Slice ID;Events",hsliceID[0]->GetNbinsX(),-1,hsliceID[0]->GetNbinsX()-1);
  hmc_ini->Sumw2();
  for (int j = 0; j < pi::nthinslices+2; ++j){
    int bin = hsliceID[0]->FindBin(j-0.5);
    hdata->SetBinContent(j+1, hsliceID[0]->GetBinContent(bin));
    hdata->SetBinError(j+1, hsliceID[0]->GetBinError(bin));
    double nmc = 0;
    for (int i = 1; i <= pi::nIntTypes; ++i){
      nmc += hsliceID[i]->GetBinContent(bin);
    }
    hmc->SetBinContent(j+1, nmc);
    hmc->SetBinError(j+1, sqrt(nmc));
    
    hdata_inc->SetBinContent(j+1, hincsliceID[0]->GetBinContent(bin));
    hdata_inc->SetBinError(j+1, hincsliceID[0]->GetBinError(bin));
    nmc = 0;
    for (int i = 1; i <= pi::nIntTypes; ++i){
      nmc += hincsliceID[i]->GetBinContent(bin);
    }
    hmc_inc->SetBinContent(j+1, nmc);
    hmc_inc->SetBinError(j+1, sqrt(nmc));
    
    hdata_ini->SetBinContent(j+1, hinisliceID[0]->GetBinContent(bin));
    hdata_ini->SetBinError(j+1, hinisliceID[0]->GetBinError(bin));
    nmc = 0;
    for (int i = 1; i <= pi::nIntTypes; ++i){
      nmc += hinisliceID[i]->GetBinContent(bin);
    }
    hmc_ini->SetBinContent(j+1, nmc);
    hmc_ini->SetBinError(j+1, sqrt(nmc));
  }

  // second loop to fill histograms
  for (int i = 0; i < pi::nIntTypes+1; ++i){
    if (i==0){
      //hsliceID[i] = (TH1D*)fdata->Get(Form("hreco_sliceID_%d_%d",pi::nCuts-1,i));
    }
    else {
      //hsliceID[i] = (TH1D*)fmc->Get(Form("hreco_sliceID_%d_%d",pi::nCuts-1,i));
      //hsliceID[i]->Scale(totaldata/totalmc);
      hsliceID[i]->Multiply(hsliceID[0]);
      hsliceID[i]->Divide(hmc);
      hincsliceID[i]->Multiply(hincsliceID[0]);
      hincsliceID[i]->Divide(hmc_inc);
      hinisliceID[i]->Multiply(hinisliceID[0]);
      hinisliceID[i]->Divide(hmc_ini);
    }
    for (int j = 0; j < pi::nthinslices+2; ++j){
      int bin = hsliceID[i]->FindBin(j-0.5);
      if (i!=0){
        double binc;
        double bine;
        if (i == pi::kMuon || i == pi::kMIDmu){
          binc = hmu->GetBinContent(bin);
          bine = hmu->GetBinError(bin);
          binc += hsliceID[i]->GetBinContent(bin)*(*sf_mu)[0];
          bine = sqrt(pow(bine,2)
                      + pow(hsliceID[i]->GetBinError(bin)*(*sf_mu)[0],2)
                      + pow(hsliceID[i]->GetBinContent(bin)*(*sf_mu)[1],2));
          hmu->SetBinContent(j+1, binc);
          hmu->SetBinError(j+1, bine);
          
          binc = hmu_inc->GetBinContent(bin);
          bine = hmu_inc->GetBinError(bin);
          binc += hincsliceID[i]->GetBinContent(bin)*(*sf_mu)[0];
          bine = sqrt(pow(bine,2)
                      + pow(hincsliceID[i]->GetBinError(bin)*(*sf_mu)[0],2)
                      + pow(hincsliceID[i]->GetBinContent(bin)*(*sf_mu)[1],2));
          hmu_inc->SetBinContent(j+1, binc);
          hmu_inc->SetBinError(j+1, bine);
          
          binc = hmu_ini->GetBinContent(bin);
          bine = hmu_ini->GetBinError(bin);
          binc += hinisliceID[i]->GetBinContent(bin)*(*sf_mu)[0];
          bine = sqrt(pow(bine,2)
                      + pow(hinisliceID[i]->GetBinError(bin)*(*sf_mu)[0],2)
                      + pow(hinisliceID[i]->GetBinContent(bin)*(*sf_mu)[1],2));
          hmu_ini->SetBinContent(j+1, binc);
          hmu_ini->SetBinError(j+1, bine);
        }
        else if (i == pi::kMIDp){
          binc = hproton->GetBinContent(bin);
          bine = hproton->GetBinError(bin);
          binc += hsliceID[i]->GetBinContent(bin)*(*sf_p)[0];
          bine = sqrt(pow(bine,2)
                      + pow(hsliceID[i]->GetBinError(bin)*(*sf_p)[0],2)
                      + pow(hsliceID[i]->GetBinContent(bin)*(*sf_p)[1],2));
          hproton->SetBinContent(j+1, binc);
          hproton->SetBinError(j+1, bine);
          
          binc = hproton_inc->GetBinContent(bin);
          bine = hproton_inc->GetBinError(bin);
          binc += hincsliceID[i]->GetBinContent(bin)*(*sf_p)[0];
          bine = sqrt(pow(bine,2)
                      + pow(hincsliceID[i]->GetBinError(bin)*(*sf_p)[0],2)
                      + pow(hincsliceID[i]->GetBinContent(bin)*(*sf_p)[1],2));
          hproton_inc->SetBinContent(j+1, binc);
          hproton_inc->SetBinError(j+1, bine);
          
          binc = hproton_ini->GetBinContent(bin);
          bine = hproton_ini->GetBinError(bin);
          binc += hinisliceID[i]->GetBinContent(bin)*(*sf_p)[0];
          bine = sqrt(pow(bine,2)
                      + pow(hinisliceID[i]->GetBinError(bin)*(*sf_p)[0],2)
                      + pow(hinisliceID[i]->GetBinContent(bin)*(*sf_p)[1],2));
          hproton_ini->SetBinContent(j+1, binc);
          hproton_ini->SetBinError(j+1, bine);
        }
        else if (i == pi::kMIDpi){
          binc = hspi->GetBinContent(bin);
          bine = hspi->GetBinError(bin);
          binc += hsliceID[i]->GetBinContent(bin)*(*sf_spi)[0];
          bine = sqrt(pow(bine,2)
                      + pow(hsliceID[i]->GetBinError(bin)*(*sf_spi)[0],2)
                      + pow(hsliceID[i]->GetBinContent(bin)*(*sf_spi)[1],2));
          hspi->SetBinContent(j+1, binc);
          hspi->SetBinError(j+1, bine);
          
          binc = hspi_inc->GetBinContent(bin);
          bine = hspi_inc->GetBinError(bin);
          binc += hincsliceID[i]->GetBinContent(bin)*(*sf_spi)[0];
          bine = sqrt(pow(bine,2)
                      + pow(hincsliceID[i]->GetBinError(bin)*(*sf_spi)[0],2)
                      + pow(hincsliceID[i]->GetBinContent(bin)*(*sf_spi)[1],2));
          hspi_inc->SetBinContent(j+1, binc);
          hspi_inc->SetBinError(j+1, bine);
          
          binc = hspi_ini->GetBinContent(bin);
          bine = hspi_ini->GetBinError(bin);
          binc += hinisliceID[i]->GetBinContent(bin)*(*sf_spi)[0];
          bine = sqrt(pow(bine,2)
                      + pow(hinisliceID[i]->GetBinError(bin)*(*sf_spi)[0],2)
                      + pow(hinisliceID[i]->GetBinContent(bin)*(*sf_spi)[1],2));
          hspi_ini->SetBinContent(j+1, binc);
          hspi_ini->SetBinError(j+1, bine);
        }
        else if (i == pi::kPiElas){
          binc = hpiel->GetBinContent(bin);
          bine = hpiel->GetBinError(bin);
          binc += hsliceID[i]->GetBinContent(bin);
          bine = sqrt(pow(bine,2)
                      + pow(hsliceID[i]->GetBinError(bin),2));
          hpiel->SetBinContent(j+1, binc);
          hpiel->SetBinError(j+1, bine);
          
          /*binc = hpiel_ini->GetBinContent(bin);
          bine = hpiel_ini->GetBinError(bin);
          binc += hinisliceID[i]->GetBinContent(bin);
          bine = sqrt(pow(bine,2)
                      + pow(hinisliceID[i]->GetBinError(bin),2));
          hpiel_ini->SetBinContent(j+1, binc);
          hpiel_ini->SetBinError(j+1, bine);*/
        }
        else if (i != pi::kPiInel){
          binc = hother->GetBinContent(bin);
          bine = hother->GetBinError(bin);
          binc += hsliceID[i]->GetBinContent(bin);
          bine = sqrt(pow(bine,2)
                      + pow(hsliceID[i]->GetBinError(bin),2));
          hother->SetBinContent(j+1, binc);
          hother->SetBinError(j+1, bine);
          
          binc = hother_inc->GetBinContent(bin);
          bine = hother_inc->GetBinError(bin);
          binc += hincsliceID[i]->GetBinContent(bin);
          bine = sqrt(pow(bine,2)
                      + pow(hincsliceID[i]->GetBinError(bin),2));
          hother_inc->SetBinContent(j+1, binc);
          hother_inc->SetBinError(j+1, bine);
          
          binc = hother_ini->GetBinContent(bin);
          bine = hother_ini->GetBinError(bin);
          binc += hinisliceID[i]->GetBinContent(bin);
          bine = sqrt(pow(bine,2)
                      + pow(hinisliceID[i]->GetBinError(bin),2));
          hother_ini->SetBinContent(j+1, binc);
          hother_ini->SetBinError(j+1, bine);
        }
      }
    }
  }

  TH1D *hsignal = (TH1D*)hdata->Clone("hsignal");
  hsignal->SetTitle("Pion interaction signal;Slice ID;Events");
  hsignal->Add(hmu,-1);
  hsignal->Add(hproton,-1);
  hsignal->Add(hspi,-1);
  hsignal->Add(hpiel,-1);
  hsignal->Add(hother,-1);
  
  TH1D *hsiginc = (TH1D*)hdata_inc->Clone("hsiginc");
  hsiginc->SetTitle("All pion incident;Slice ID;Events");
  hsiginc->Add(hmu_inc,-1);
  hsiginc->Add(hproton_inc,-1);
  hsiginc->Add(hspi_inc,-1);
  hsiginc->Add(hother_inc,-1);
  
  TH1D *hsigini = (TH1D*)hdata_ini->Clone("hsigini");
  hsigini->SetTitle("All pion initial;Slice ID;Events");
  hsigini->Add(hmu_ini,-1);
  hsigini->Add(hproton_ini,-1);
  hsigini->Add(hspi_ini,-1);
  hsigini->Add(hother_ini,-1);
  
  // unfolding
  RooUnfoldResponse *response_SliceID_Inc = (RooUnfoldResponse*)fmc->Get("response_SliceID_Inc");
  RooUnfoldBayes unfold_Inc (response_SliceID_Inc, hsiginc, 20);
  RooUnfoldResponse *response_SliceID_Int = (RooUnfoldResponse*)fmc->Get("response_SliceID_Int");
  RooUnfoldBayes unfold_Int (response_SliceID_Int, hsignal, 20);
  RooUnfoldResponse *response_SliceID_Ini = (RooUnfoldResponse*)fmc->Get("response_SliceID_Ini");
  RooUnfoldBayes unfold_Ini (response_SliceID_Ini, hsigini, 40);
  
  TH1D *hsiginc_uf;
  TH1D *hsignal_uf;
  TH1D *hsigini_uf;
  hsiginc_uf = (TH1D*)unfold_Inc.Hreco();
  hsiginc_uf->SetNameTitle("hsiginc_uf", "Unfolded incident signal;Slice ID;Events");
  hsignal_uf = (TH1D*)unfold_Int.Hreco();
  hsignal_uf->SetNameTitle("hsignal_uf", "Unfolded interaction signal;Slice ID;Events");
  hsigini_uf = (TH1D*)unfold_Ini.Hreco();
  hsigini_uf->SetNameTitle("hsigini_uf", "Unfolded initial signal;Slice ID;Events");
  
  // get Ninc, Nint, Nini
  double Ninc[pi::nthinslices] = {0};
  double Nint[pi::nthinslices] = {0};
  double Nini[pi::nthinslices] = {0};
  double Nina[pi::nthinslices] = {0};
  double err_inc[pi::nthinslices] = {0};
  double err_int[pi::nthinslices] = {0};
  double err_ini[pi::nthinslices] = {0};
  double err_ina[pi::nthinslices] = {0};
  double SliceID[pi::nthinslices] = {0};

  for (int i = 0; i<pi::nthinslices; ++i){
    SliceID[i] = i+1;
    Nint[i] = hsignal_uf->GetBinContent(i+2);
    err_int[i] = hsignal_uf->GetBinError(i+2);
    Nini[i] = hsigini_uf->GetBinContent(i+2);
    err_ini[i] = hsigini_uf->GetBinError(i+2);
    Nina[i] = hsiginc_uf->GetBinContent(i+2);
    err_ina[i] = hsiginc_uf->GetBinError(i+2);
    for (int j = 0; j<=i; ++j){
      Ninc[i] += hsigini_uf->GetBinContent(j+2);
      err_inc[i] += pow(hsigini_uf->GetBinError(j+2),2);
    }
    for (int j = 0; j<=i-1; ++j){
      Ninc[i] -= hsiginc_uf->GetBinContent(j+2);
      err_inc[i] += pow(hsiginc_uf->GetBinError(j+2),2);
    }
    err_inc[i] = sqrt(err_inc[i]);
  }
  TGraphErrors *gr_inc = new TGraphErrors(pi::nthinslices, SliceID, Ninc, 0, err_inc);
  gr_inc->SetNameTitle("gr_inc", "Incident number;Slice ID;Events");
  gr_inc->Write();
  TGraphErrors *gr_int = new TGraphErrors(pi::nthinslices, SliceID, Nint, 0, err_int);
  gr_int->SetNameTitle("gr_int", "Interaction number;Slice ID;Events");
  gr_int->Write();
  TGraphErrors *gr_ini = new TGraphErrors(pi::nthinslices, SliceID, Nini, 0, err_ini);
  gr_ini->SetNameTitle("gr_ini", "Initial number;Slice ID;Events");
  gr_ini->Write();
  TGraphErrors *gr_ina = new TGraphErrors(pi::nthinslices, SliceID, Nina, 0, err_ina);
  gr_ina->SetNameTitle("gr_ina", "Interaction_allpion number;Slice ID;Events");
  gr_ina->Write();
  /*TGraphErrors *gr_trueincE = (TGraphErrors*)fmc->Get("gr_trueincE");
  gr_trueincE->SetNameTitle("gr_trueincE", "True incident energy;Slice ID;Energy (MeV)");
  gr_trueincE->Write();
  TGraphErrors *gr_recoincE = (TGraphErrors*)fdata->Get("gr_recoincE");
  gr_recoincE->SetNameTitle("gr_recoincE", "Reco incident energy;Slice ID;Energy (MeV)");
  gr_recoincE->Write();*/
  
  // Calculate cross-section
  BetheBloch bb(211);
  double NA=6.02214076e23;
  double MAr=39.95; //gmol
  double Density = 1.4; // g/cm^3
  double xs[pi::nthinslices] = {0};
  double err_xs[pi::nthinslices] = {0};
  double KE[pi::nthinslices] = {0};
  double err_KE[pi::nthinslices] = {0};
  double dEdx[pi::nthinslices] = {0};
  for (int i = 0; i<pi::nthinslices; ++i){
    KE[i] = pi::plim - (i+0.5)*pi::Eslicewidth;
    err_KE[i] = pi::Eslicewidth/2;
    dEdx[i] = bb.meandEdx(KE[i]); // MeV/cm
    xs[i] = dEdx[i]*MAr/(Density*NA*pi::Eslicewidth)*log(Ninc[i]/(Ninc[i]-Nint[i]))*1e27;
    err_xs[i] = dEdx[i]*MAr/(Density*NA*pi::Eslicewidth)*sqrt(pow(Nint[i]*err_inc[i]/Ninc[i]/(Ninc[i]-Nint[i]),2)+pow(err_int[i]/(Ninc[i]-Nint[i]),2))*1e27;
  }
  TGraphErrors *gr_recoxs = new TGraphErrors(pi::nthinslices, KE, xs, err_KE, err_xs);
  gr_recoxs->SetNameTitle("gr_recoxs", "Reco cross-section;Energy (MeV); Cross-section (mb)");
  gr_recoxs->Write();
  /*TFile f2("../files/exclusive_xsec.root");
  TGraph *total_inel_KE = (TGraph*)f2.Get("total_inel_KE");
  TCanvas *c1 = new TCanvas("c1", "c1");
  gr_recoxs->SetTitle("Pion Inelastic Cross Section");
  gr_recoxs->GetXaxis()->SetTitle("Pion Kinetic Energy (MeV)");
  gr_recoxs->GetXaxis()->SetRangeUser(360, 900);
  gr_recoxs->GetYaxis()->SetTitle("#sigma_{inelastic} (mb)");
  gr_recoxs->GetYaxis()->SetRangeUser(400, 900);
  gr_recoxs->SetLineWidth(2);
  gr_recoxs->Draw("ape");
  total_inel_KE->SetLineColor(2);
  total_inel_KE->Draw("c");
  c1->Print("recoxs.png");*/
  
  // test sample validation
  TH1D *hval_siginc_reco = (TH1D*)fmc->Get("h_recosliceid_pion_cuts");
  hval_siginc_reco->Write("hval_siginc_reco");
  TH1D *hval_signal_reco = (TH1D*)fmc->Get("h_recosliceid_pioninelastic_cuts");
  hval_signal_reco->Write("hval_signal_reco");
  TH1D *hval_sigini_reco = (TH1D*)fmc->Get("h_recoinisliceid_pion_cuts");
  hval_sigini_reco->Write("hval_sigini_reco");
  TH1D *hval_trueinc = (TH1D*)fmc->Get("h_truesliceid_pion_all");
  hval_trueinc->Write("hval_trueinc");
  TH1D *hval_trueint = (TH1D*)fmc->Get("h_truesliceid_pioninelastic_all");
  hval_trueint->Write("hval_trueint");
  TH1D *hval_trueini = (TH1D*)fmc->Get("h_trueinisliceid_pion_all");
  hval_trueini->Write("hval_trueini");
  double Ninc_t[pi::nthinslices] = {0};
  double Nint_t[pi::nthinslices] = {0};
  double Nini_t[pi::nthinslices] = {0};
  double Nina_t[pi::nthinslices] = {0};
  double err_inc_t[pi::nthinslices] = {0};
  double err_int_t[pi::nthinslices] = {0};
  double err_ini_t[pi::nthinslices] = {0};
  double err_ina_t[pi::nthinslices] = {0};
  for (int i = 0; i<pi::nthinslices; ++i){
    Nint_t[i] = hval_trueint->GetBinContent(i+2);
    err_int_t[i] = hval_trueint->GetBinError(i+2);
    Nini_t[i] = hval_trueini->GetBinContent(i+2);
    err_ini_t[i] = hval_trueini->GetBinError(i+2);
    Nina_t[i] = hval_trueinc->GetBinContent(i+2);
    err_ina_t[i] = hval_trueinc->GetBinError(i+2);
    for (int j = 0; j<=i; ++j){
      Ninc_t[i] += hval_trueini->GetBinContent(j+2);
      err_inc_t[i] += pow(hval_trueini->GetBinError(j+2),2);
    }
    for (int j = 0; j<=i-1; ++j){
      Ninc_t[i] -= hval_trueinc->GetBinContent(j+2);
      err_inc_t[i] += pow(hval_trueinc->GetBinError(j+2),2);
    }
    err_inc_t[i] = sqrt(err_inc_t[i]);
  }
  TGraphErrors *gr_inc_t = new TGraphErrors(pi::nthinslices, SliceID, Ninc_t, 0, err_inc_t);
  gr_inc_t->SetNameTitle("gr_inc_t", "Incident number;Slice ID;Events");
  gr_inc_t->Write();
  TGraphErrors *gr_int_t = new TGraphErrors(pi::nthinslices, SliceID, Nint_t, 0, err_int_t);
  gr_int_t->SetNameTitle("gr_int_t", "Interaction number;Slice ID;Events");
  gr_int_t->Write();
  TGraphErrors *gr_ini_t = new TGraphErrors(pi::nthinslices, SliceID, Nini_t, 0, err_ini_t);
  gr_ini_t->SetNameTitle("gr_ini_t", "Initial number;Slice ID;Events");
  gr_ini_t->Write();
  TGraphErrors *gr_ina_t = new TGraphErrors(pi::nthinslices, SliceID, Nina_t, 0, err_ina_t);
  gr_ina_t->SetNameTitle("gr_ina_t", "Interaction_allpion number;Slice ID;Events");
  gr_ina_t->Write();
  double xs_t[pi::nthinslices] = {0};
  double err_xs_t[pi::nthinslices] = {0};
  for (int i = 0; i<pi::nthinslices; ++i){
    xs_t[i] = dEdx[i]*MAr/(Density*NA*pi::Eslicewidth)*log(Ninc_t[i]/(Ninc_t[i]-Nint_t[i]))*1e27;
    err_xs_t[i] = dEdx[i]*MAr/(Density*NA*pi::Eslicewidth)*sqrt(pow(Nint_t[i]*err_inc_t[i]/Ninc_t[i]/(Ninc_t[i]-Nint_t[i]),2)+pow(err_int_t[i]/(Ninc_t[i]-Nint_t[i]),2))*1e27;
  }
  TGraphErrors *gr_truexs = new TGraphErrors(pi::nthinslices, KE, xs_t, err_KE, err_xs_t);
  gr_truexs->SetNameTitle("gr_truexs", "Reco cross-section;Energy (MeV); Cross-section (mb)");
  gr_truexs->Write();
  TGraphErrors *gr_truexs_allMC = (TGraphErrors*)fmc->Get("gr_truexs");
  gr_truexs_allMC->Write("gr_truexs_allMC");
  
  fout->Write();
  fout->Close();
  return 0;
}
