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
#include "TStyle.h"
#include "json/json.h"
#include <fstream>
#include <iostream>
#include <string>
#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"
#include "RooUnfoldSvd.h"
#include "RooUnfoldBinByBin.h"
#include "TH2D.h"
#include "TH3D.h"
#include "BetheBloch.h"
#include "TRandom3.h"

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

  //////////  get backgrounds scale factors  //////////
  TVectorD sf(2); sf[0] = 1.; sf[1] = 0.;
  TVectorD *sf_mu = (TVectorD*)fbkg->Get("sf_mu");
  TVectorD *sf_p = (TVectorD*)fbkg->Get("sf_p");
  TVectorD *sf_spi = (TVectorD*)fbkg->Get("sf_spi");
  /*TRandom3 *r3 = new TRandom3(0);
  (*sf_mu)[0] = r3->Gaus((*sf_mu)[0], (*sf_mu)[1]);
  (*sf_p)[0] = r3->Gaus((*sf_p)[0], (*sf_p)[1]);
  (*sf_spi)[0] = r3->Gaus((*sf_spi)[0], (*sf_spi)[1]);
  cout<<"$$$rdm_bkgsc<<endl;*/
  cout<<"Muon scaling factor: "<<(*sf_mu)[0]<<"+-"<<(*sf_mu)[1]<<endl;
  cout<<"Proton scaling factor: "<<(*sf_p)[0]<<"+-"<<(*sf_p)[1]<<endl;
  cout<<"Pion scaling factor: "<<(*sf_spi)[0]<<"+-"<<(*sf_spi)[1]<<endl;
  
  (*sf_mu)[1] = 0; // the uncertainty should be treated at the end (otherwise we should add this nuisance parameter into the total covvariance matrix in the error propagation)
  (*sf_p)[1] = 0;
  (*sf_spi)[1] = 0;
  ///END  get backgrounds scale factors  //////////

  
  //////////  define data histograms  //////////
  TH1D *hsliceID[pi::nIntTypes+1];
  TH1D *hdata = new TH1D("hdata","Data;Slice ID;Events",pi::reco_nbins,pi::reco_bins); // h_recosliceid_allevts_cuts (hreco_sliceID_6_0)
  TH1D *hproton = new TH1D("hproton","Proton background;Slice ID;Events",pi::reco_nbins,pi::reco_bins);
  TH1D *hmu = new TH1D("hmu","Muon background;Slice ID;Events",pi::reco_nbins,pi::reco_bins);
  TH1D *hspi = new TH1D("hspi","Secondary pion background;Slice ID;Events",pi::reco_nbins,pi::reco_bins);
  TH1D *hpiel = new TH1D("hpiel","Pion elastic;Slice ID;Events",pi::reco_nbins,pi::reco_bins);
  TH1D *hother = new TH1D("hother","Other backgrounds;Slice ID;Events",pi::reco_nbins,pi::reco_bins);
  hdata->Sumw2();
  hproton->Sumw2();
  hmu->Sumw2();
  hspi->Sumw2();
  hpiel->Sumw2();
  hother->Sumw2();
  
  TH1D *hincsliceID[pi::nIntTypes+1];
  TH1D *hdata_inc = new TH1D("hdata_inc","Data_inc;Slice ID;Events",pi::reco_nbins,pi::reco_bins);
  TH1D *hproton_inc = new TH1D("hproton_inc","Proton_inc background;Slice ID;Events",pi::reco_nbins,pi::reco_bins);
  TH1D *hmu_inc = new TH1D("hmu_inc","Muon_inc background;Slice ID;Events",pi::reco_nbins,pi::reco_bins);
  TH1D *hspi_inc = new TH1D("hspi_inc","Secondary pion_inc background;Slice ID;Events",pi::reco_nbins,pi::reco_bins);
  //TH1D *hpiel_inc = new TH1D("hpiel_inc","Pion elastic_inc;Slice ID;Events",pi::reco_nbins,pi::reco_bins);
  TH1D *hother_inc = new TH1D("hother_inc","Other_inc backgrounds;Slice ID;Events",pi::reco_nbins,pi::reco_bins);
  hdata_inc->Sumw2();
  hproton_inc->Sumw2();
  hmu_inc->Sumw2();
  hspi_inc->Sumw2();
  //hpiel_inc->Sumw2();
  hother_inc->Sumw2();
  
  TH1D *hinisliceID[pi::nIntTypes+1];
  TH1D *hdata_ini = new TH1D("hdata_ini","Data_ini;Slice ID;Events",pi::reco_nbins,pi::reco_bins);
  TH1D *hproton_ini = new TH1D("hproton_ini","Proton_ini background;Slice ID;Events",pi::reco_nbins,pi::reco_bins);
  TH1D *hmu_ini = new TH1D("hmu_ini","Muon_ini background;Slice ID;Events",pi::reco_nbins,pi::reco_bins);
  TH1D *hspi_ini = new TH1D("hspi_ini","Secondary pion_ini background;Slice ID;Events",pi::reco_nbins,pi::reco_bins);
  //TH1D *hpiel_ini = new TH1D("hpiel_ini","Pion elastic_ini;Slice ID;Events",pi::reco_nbins,pi::reco_bins);
  TH1D *hother_ini = new TH1D("hother_ini","Other_ini backgrounds;Slice ID;Events",pi::reco_nbins,pi::reco_bins);
  hdata_ini->Sumw2();
  hproton_ini->Sumw2();
  hmu_ini->Sumw2();
  hspi_ini->Sumw2();
  //hpiel_ini->Sumw2();
  hother_ini->Sumw2();
  
  TH3D *h3DsliceID[pi::nIntTypes+1];
  TH3D *hdata_3D = new TH3D("hdata_3D","Data_3D;Slice ID;Events",pi::reco_nbins,pi::reco_bins,pi::reco_nbins,pi::reco_bins,pi::reco_nbins,pi::reco_bins);
  TH3D *hproton_3D = new TH3D("hproton_3D","Proton_3D background;Slice ID;Events",pi::reco_nbins,pi::reco_bins,pi::reco_nbins,pi::reco_bins,pi::reco_nbins,pi::reco_bins);
  TH3D *hmu_3D = new TH3D("hmu_3D","Muon_3D background;Slice ID;Events",pi::reco_nbins,pi::reco_bins,pi::reco_nbins,pi::reco_bins,pi::reco_nbins,pi::reco_bins);
  TH3D *hspi_3D = new TH3D("hspi_3D","Secondary pion_3D background;Slice ID;Events",pi::reco_nbins,pi::reco_bins,pi::reco_nbins,pi::reco_bins,pi::reco_nbins,pi::reco_bins);
  TH3D *hpiel_3D = new TH3D("hpiel_3D","Pion elastic_3D;Slice ID;Events",pi::reco_nbins,pi::reco_bins,pi::reco_nbins,pi::reco_bins,pi::reco_nbins,pi::reco_bins);
  TH3D *hother_3D = new TH3D("hother_3D","Other_3D backgrounds;Slice ID;Events",pi::reco_nbins,pi::reco_bins,pi::reco_nbins,pi::reco_bins,pi::reco_nbins,pi::reco_bins);

  
  /// first loop to get total number of events
  for (int i = 0; i < pi::nIntTypes+1; ++i){
    if (i==0){
      hsliceID[i] = (TH1D*)fdata->Get(Form("hreco_sliceID_%d_%d",pi::nCuts-1,i));
      hincsliceID[i] = (TH1D*)fdata->Get(Form("hreco_incsliceID_%d_%d",pi::nCuts-1,i));
      hinisliceID[i] = (TH1D*)fdata->Get(Form("hreco_inisliceID_%d_%d",pi::nCuts-1,i));
      h3DsliceID[i] = (TH3D*)fdata->Get(Form("hreco_3DsliceID_%d_%d",pi::nCuts-1,i));
    }
    else {
      hsliceID[i] = (TH1D*)fmc->Get(Form("hreco_sliceID_%d_%d",pi::nCuts-1,i));
      hincsliceID[i] = (TH1D*)fmc->Get(Form("hreco_incsliceID_%d_%d",pi::nCuts-1,i));
      hinisliceID[i] = (TH1D*)fmc->Get(Form("hreco_inisliceID_%d_%d",pi::nCuts-1,i));
      h3DsliceID[i] = (TH3D*)fmc->Get(Form("hreco_3DsliceID_%d_%d",pi::nCuts-1,i));
    }
  }
  TH1D *hmc = new TH1D("hmc","MC;Slice ID;Events",pi::reco_nbins,pi::reco_bins);
  hmc->Sumw2();
  TH1D *hmc_inc = new TH1D("hmc_inc","MC_inc;Slice ID;Events",pi::reco_nbins,pi::reco_bins);
  hmc_inc->Sumw2();
  TH1D *hmc_ini = new TH1D("hmc_ini","MC_ini;Slice ID;Events",pi::reco_nbins,pi::reco_bins);
  hmc_ini->Sumw2();
  TH3D *hmc_3D = new TH3D("hmc_3D","MC_3D;Slice ID;Events",pi::reco_nbins,pi::reco_bins,pi::reco_nbins,pi::reco_bins,pi::reco_nbins,pi::reco_bins);
  for (int j = 0; j < pi::reco_nbins; ++j){
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
    
    for (int k = 0; k < pi::reco_nbins; ++k){
      int bink = hsliceID[0]->FindBin(k-0.5);
      for (int l = 0; l < pi::reco_nbins; ++l){
        int binl = hsliceID[0]->FindBin(l-0.5);
        hdata_3D->SetBinContent(j+1, k+1, l+1, h3DsliceID[0]->GetBinContent(bin, bink, binl));
        hdata_3D->SetBinError(j+1, k+1, l+1, h3DsliceID[0]->GetBinError(bin, bink, binl));
        nmc = 0;
        for (int i = 1; i <= pi::nIntTypes; ++i){
          nmc += h3DsliceID[i]->GetBinContent(bin, bink, binl);
        }
        hmc_3D->SetBinContent(j+1, k+1, l+1, nmc);
        hmc_3D->SetBinError(j+1, k+1, l+1, sqrt(nmc));
      }
    }
  }
  ///END  define data histograms  //////////


  //////////  second loop to fill histograms  //////////
  for (int i = 0; i < pi::nIntTypes+1; ++i){
    if (i==0){
    }
    else {
      // try which scale method to use?
      /*hsliceID[i]->Multiply(hsliceID[0]);
      hsliceID[i]->Divide(hmc);
      hincsliceID[i]->Multiply(hincsliceID[0]);
      hincsliceID[i]->Divide(hmc_inc);
      hinisliceID[i]->Multiply(hinisliceID[0]);
      hinisliceID[i]->Divide(hmc_ini);
      h3DsliceID[i]->Multiply(h3DsliceID[0]);
      h3DsliceID[i]->Divide(hmc_3D);*/
      /*hsliceID[i]->Scale(hdata->Integral()/hmc->Integral());
      hincsliceID[i]->Scale(hdata_inc->Integral()/hmc_inc->Integral());
      hinisliceID[i]->Scale(hdata_ini->Integral()/hmc_ini->Integral());
      h3DsliceID[i]->Scale(hdata_3D->Integral()/hmc_3D->Integral());*/
    }
    for (int j = 0; j < pi::reco_nbins; ++j){
      int bin = hsliceID[i]->FindBin(j-0.5);
      if (i!=0){
        double binc;
        double bine;
        if (i == pi::kMuon || i == pi::kMIDmu){
          binc = hmu->GetBinContent(bin); // only possibly nonzero here because there are two pi::types
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
          for (int k = 0; k < pi::reco_nbins; ++k){
            int bink = hsliceID[i]->FindBin(k-0.5);
            for (int l = 0; l < pi::reco_nbins; ++l){
              int binl = hsliceID[i]->FindBin(l-0.5);
              binc = hmu_3D->GetBinContent(bin, bink, binl);
              bine = hmu_3D->GetBinError(bin, bink, binl);
              binc += h3DsliceID[i]->GetBinContent(bin, bink, binl)*(*sf_mu)[0];
              bine = sqrt(pow(bine,2)
                          + pow(h3DsliceID[i]->GetBinError(bin, bink, binl)*(*sf_mu)[0],2)
                          + pow(h3DsliceID[i]->GetBinContent(bin, bink, binl)*(*sf_mu)[1],2));
              hmu_3D->SetBinContent(j+1, k+1, l+1, binc);
              hmu_3D->SetBinError(j+1, k+1, l+1, bine);
            }
          }
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
          
          for (int k = 0; k < pi::reco_nbins; ++k){
            int bink = hsliceID[i]->FindBin(k-0.5);
            for (int l = 0; l < pi::reco_nbins; ++l){
              int binl = hsliceID[i]->FindBin(l-0.5);
              binc = hproton_3D->GetBinContent(bin, bink, binl);
              bine = hproton_3D->GetBinError(bin, bink, binl);
              binc += h3DsliceID[i]->GetBinContent(bin, bink, binl)*(*sf_p)[0];
              bine = sqrt(pow(bine,2)
                          + pow(h3DsliceID[i]->GetBinError(bin, bink, binl)*(*sf_p)[0],2)
                          + pow(h3DsliceID[i]->GetBinContent(bin, bink, binl)*(*sf_p)[1],2));
              hproton_3D->SetBinContent(j+1, k+1, l+1, binc);
              hproton_3D->SetBinError(j+1, k+1, l+1, bine);
            }
          }
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
          
          for (int k = 0; k < pi::reco_nbins; ++k){
            int bink = hsliceID[i]->FindBin(k-0.5);
            for (int l = 0; l < pi::reco_nbins; ++l){
              int binl = hsliceID[i]->FindBin(l-0.5);
              binc = hspi_3D->GetBinContent(bin, bink, binl);
              bine = hspi_3D->GetBinError(bin, bink, binl);
              binc += h3DsliceID[i]->GetBinContent(bin, bink, binl)*(*sf_spi)[0];
              bine = sqrt(pow(bine,2)
                          + pow(h3DsliceID[i]->GetBinError(bin, bink, binl)*(*sf_spi)[0],2)
                          + pow(h3DsliceID[i]->GetBinContent(bin, bink, binl)*(*sf_spi)[1],2));
              hspi_3D->SetBinContent(j+1, k+1, l+1, binc);
              hspi_3D->SetBinError(j+1, k+1, l+1, bine);
            }
          }
        }
        else if (i == pi::kPiElas){
          if (j == 0) {
            hpiel->SetBinContent(1, 0);
            hpiel->SetBinError(1, 0);
          }
          else {
            double pielval = hsliceID[i]->GetBinContent(bin);
            double pielerr = hsliceID[i]->GetBinError(bin);
            binc = hpiel->GetBinContent(1);
            bine = hpiel->GetBinError(1);
            binc -= pielval;
            bine = sqrt(pow(bine,2) + pow(pielerr,2));
            hpiel->SetBinContent(1, binc);
            hpiel->SetBinError(1, bine);
            hpiel->SetBinContent(j+1, pielval);
            hpiel->SetBinError(j+1, pielerr);
          }
          for (int k = 0; k < pi::reco_nbins; ++k){
            int bink = hsliceID[i]->FindBin(k-0.5);
            for (int l = 0; l < pi::reco_nbins; ++l){
              if (l == 0) {
                hpiel_3D->SetBinContent(j+1, k+1, 1, 0);
                hpiel_3D->SetBinError(j+1, k+1, 1, 0);
              }
              else {
                int binl = hsliceID[i]->FindBin(l-0.5);
                double pielval = h3DsliceID[i]->GetBinContent(bin, bink, binl);
                double pielerr = h3DsliceID[i]->GetBinError(bin, bink, binl);
                binc = hpiel_3D->GetBinContent(bin, bink, 1);
                bine = hpiel_3D->GetBinError(bin, bink, 1);
                binc -= pielval;
                bine = sqrt(pow(bine,2) + pow(pielerr,2));
                hpiel_3D->SetBinContent(j+1, k+1, 1, binc);
                hpiel_3D->SetBinError(j+1, k+1, 1, bine);
                hpiel_3D->SetBinContent(j+1, k+1, l+1, pielval);
                hpiel_3D->SetBinError(j+1, k+1, l+1, pielerr);
              }
            }
          }
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
          
          for (int k = 0; k < pi::reco_nbins; ++k){
            int bink = hsliceID[i]->FindBin(k-0.5);
            for (int l = 0; l < pi::reco_nbins; ++l){
              int binl = hsliceID[i]->FindBin(l-0.5);
              binc = hother_3D->GetBinContent(bin, bink, binl);
              bine = hother_3D->GetBinError(bin, bink, binl);
              binc += h3DsliceID[i]->GetBinContent(bin, bink, binl);
              bine = sqrt(pow(bine,2)
                          + pow(h3DsliceID[i]->GetBinError(bin, bink, binl),2));
              hother_3D->SetBinContent(j+1, k+1, l+1, binc);
              hother_3D->SetBinError(j+1, k+1, l+1, bine);
            }
          }
        }
      }
    }
  }
  ///END  second loop to fill histograms  //////////
  
  
  //////////  get signal histograms by subtracting MCbkg histograms  //////////
  double norm_mc2data = hdata_3D->Integral()/hmc_3D->Integral();
  
  hmu->Scale(norm_mc2data);
  hproton->Scale(norm_mc2data);
  hspi->Scale(norm_mc2data);
  hpiel->Scale(norm_mc2data);
  hother->Scale(norm_mc2data);
  TH1D *hsignal = (TH1D*)hdata->Clone("hsignal");
  hsignal->Add(hmu,-1);
  hsignal->Add(hproton,-1);
  hsignal->Add(hspi,-1);
  hsignal->Add(hpiel,-1);
  hsignal->Add(hother,-1);
  
  hmu_inc->Scale(norm_mc2data);
  hproton_inc->Scale(norm_mc2data);
  hspi_inc->Scale(norm_mc2data);
  hother_inc->Scale(norm_mc2data);
  TH1D *hsiginc = (TH1D*)hdata_inc->Clone("hsiginc");
  hsiginc->Add(hmu_inc,-1);
  hsiginc->Add(hproton_inc,-1);
  hsiginc->Add(hspi_inc,-1);
  hsiginc->Add(hother_inc,-1);
  
  hmu_ini->Scale(norm_mc2data);
  hproton_ini->Scale(norm_mc2data);
  hspi_ini->Scale(norm_mc2data);
  hother_ini->Scale(norm_mc2data);
  TH1D *hsigini = (TH1D*)hdata_ini->Clone("hsigini");
  hsigini->Add(hmu_ini,-1);
  hsigini->Add(hproton_ini,-1);
  hsigini->Add(hspi_ini,-1);
  hsigini->Add(hother_ini,-1);
  
  hmu_3D->Scale(norm_mc2data);
  hproton_3D->Scale(norm_mc2data);
  hspi_3D->Scale(norm_mc2data);
  hpiel_3D->Scale(norm_mc2data);
  hother_3D->Scale(norm_mc2data);
  TH3D *hsig3D = (TH3D*)hdata_3D->Clone("hsig3D");
  hsig3D->SetTitle("All pion 3D;Slice ID;Events");
  hsig3D->Add(hmu_3D,-1);
  hsig3D->Add(hproton_3D,-1);
  hsig3D->Add(hspi_3D,-1);
  hsig3D->Add(hpiel_3D,-1);
  hsig3D->Add(hother_3D,-1);
  ///END  get signal histograms by subtracting MCbkg histograms  //////////
  
  
  const int reco_nbins3D = pi::reco_nbins3D;
  const int true_nbins3D = pi::true_nbins3D;
  double central_datasig[reco_nbins3D] = {41.0708,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1.68426,0,0,0,0,0,0,0,0,0,0,5.93295,3.24624,0,0,0,0,0,0,0,0,0,8.91404,3.60095,1.07531,0,0,0,0,0,0,0,0,3.1897,4.73233,4.31457,1.25746,0,0,0,0,0,0,0,4.07467,6.36021,2.12795,1.70263,0,0,0,0,0,0,0,2.30315,7.02867,3.25404,2.73375,0.710337,0,0,0,0,0,0,5.50353,3.94857,6.09099,2.28615,0,0,0,0,0,0,0,11.8293,18.8412,10.1564,4.81078,0,0,0,0,0,0,0,1.40144,5.77928,7.98471,4.41632,1.77615,0.583961,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,157.828,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,31.0058,761.973,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,35.1708,2168.08,577.559,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,36.5369,1954.96,2226.48,517.454,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,26.4788,1445.16,2051.63,1560.15,202.416,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,3.36103,1153.98,1697.15,1400.53,505.525,32.8944,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1.26203,865.993,1283.91,1171.79,425.373,76.928,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-2.98093,711.046,969.519,816.92,352.32,68.9363,5.78502,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1.37209,849.608,1654.87,1570.48,630.057,132.189,18.3109,5.28728,2.09576,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,13.3916,38.0636,278.587,187.588,88.2262,6.08636,10.6703,7.19152,6.38304}; // data - bkg
  TVectorD vcentral_truth(reco_nbins3D);
  for (int i=0; i<pi::reco_nbins; ++i)
    for (int j=0; j<pi::reco_nbins; ++j)
      for (int k=0; k<pi::reco_nbins; ++k) {
        int idx = i*pow(pi::reco_nbins,2) + j*pi::reco_nbins + k;
        //vcentral_truth(idx) = central_truth[idx];
        //hsig3D->SetBinContent(k+1, j+1, i+1, central_datasig[idx]);
        //hsig3D->SetBinError(k+1, j+1, i+1, hdata_3D->GetBinError(k+1, j+1, i+1)); // exclude bkg sys error
      }
  /// use Cholesky decomposition to generate sample of histogram with bin correlations
  /*FILE *fcholsqrt_Inc=fopen("../../thinslice_sys_nominal/build/toys/cholsqrt_MCXS_inc_1110.txt","r");
  if (!fcholsqrt_Inc) {
    cout<<"cholsqrt_inc_MCXS not found!"<<endl;
    return 1;
  }
  TMatrixD mcholsqrt_Inc_input(pi::reco_nbins, pi::reco_nbins);
  double vvchol;
  for(int i=0; i<pi::reco_nbins; i++) {
    for(int j=0; j<pi::reco_nbins; j++) {
      fscanf(fcholsqrt_Inc, "%lf", &vvchol);
      mcholsqrt_Inc_input(i, j) = vvchol;
    }
  }
  TRandom3 *r3 = new TRandom3(0);
  TVectorD rdmtt(pi::reco_nbins);
  for (int i=0; i<pi::reco_nbins; i++) {
    rdmtt(i) = r3->Gaus(0, 1);
  }
  TVectorD rdmss(pi::reco_nbins);
  rdmss = mcholsqrt_Inc_input * rdmtt;
  for (int i=1; i<=pi::reco_nbins; i++) {
    hsiginc->SetBinContent(i, rdmss(i-1) + hsiginc->GetBinContent(i));
  }*/
  
  
  RooUnfoldResponse *response_SliceID_3D = (RooUnfoldResponse*)fmc->Get("response_SliceID_3D");
  RooUnfoldResponse *response_SliceID_1D_eff = (RooUnfoldResponse*)fmc->Get("response_SliceID_1D"); // should be the same with 3D excluding empty bins
  RooUnfoldResponse *response_SliceID_1D = (RooUnfoldResponse*)fmc->Get("response_SliceID_1D_noeff");
  TH3D *hmeas_3D = (TH3D*)response_SliceID_3D->Hmeasured();
  TH3D *htruth_3D = (TH3D*)response_SliceID_3D->Htruth();
  
  
  //////////  get the 1D index and the hsig1D histogram for q3D unfolding  //////////
  int nmeas_3D = 0; // number of non-empty bins of measure spectrum
  int ntruth_3D = 0; // number of non-empty bins of truth spectrum
  int ntruth_3D_eff = 0; // number of non-empty bins of truth spectrum including tiny bins
  int idx_meas1D[reco_nbins3D]; // indices of non-empty bins of measure spectrum (0 if emtpy)
  int idx_truth1D[true_nbins3D]; // indices of non-empty bins of truth spectrum (0 if emtpy)
  int idx_truth1D_eff[true_nbins3D]; // indices of non-empty bins of truth spectrum (0 if emtpy)
  
  bool get_1Didx = false; /// true in the case when you have whole MC as truth MC, and get the non-empty indices; false when you already pasted the tmp idx below from whole MC and ready for fake or real data study
  if (get_1Didx) { /// when you have whole MC as truth MC, and get the non-empty indices
    cout<<"### idx_meas1D"<<endl;
    for (int i=0; i<pi::reco_nbins; ++i)
      for (int j=0; j<pi::reco_nbins; ++j)
        for (int k=0; k<pi::reco_nbins; ++k) {
          int idx = i*pow(pi::reco_nbins,2) + j*pi::reco_nbins + k;
          idx_meas1D[idx] = 0;
          if (hmeas_3D->GetBinContent(k+1, j+1, i+1)!=0) {
            ++nmeas_3D;
            idx_meas1D[idx] = nmeas_3D;
          }
          cout<<idx_meas1D[idx]<<",";
        }
    cout<<endl<<"### idx_truth1D (has eff)"<<endl;
    for (int i=0; i<pi::true_nbins; ++i)
      for (int j=0; j<pi::true_nbins; ++j)
        for (int k=0; k<pi::true_nbins; ++k) {
          int idx = i*pow(pi::true_nbins,2) + j*pi::true_nbins + k;
          idx_truth1D_eff[idx] = 0;
          if (htruth_3D->GetBinContent(k+1, j+1, i+1)!=0) {
            ++ntruth_3D_eff;
            idx_truth1D_eff[idx] = ntruth_3D_eff;
          }
          cout<<idx_truth1D_eff[idx]<<",";
        }
    cout<<endl;
    fout->Write();
    fout->Close();
    return 0;
  }
  else { /// when you already pasted the tmp idx below from whole MC and ready for fake or real data study
    int tmp_idx_meas1D[pi::reco_nbins3D]={1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,3,4,0,0,0,0,0,0,0,0,0,0,0,5,6,7,0,0,0,0,0,0,0,0,0,0,8,9,10,11,0,0,0,0,0,0,0,0,0,12,13,14,15,0,0,0,0,0,0,0,0,0,0,16,17,18,19,0,0,0,0,0,0,0,0,0,20,21,22,0,0,0,0,0,0,0,0,0,23,24,25,26,0,0,0,0,0,0,0,0,0,0,27,28,29,0,0,0,0,0,0,0,0,0,0,30,31,32,33,0,0,0,0,0,0,0,0,0,34,35,36,37,38,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,39,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,40,41,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,42,43,44,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,45,46,47,48,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,49,50,51,52,53,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,54,55,56,57,58,59,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,60,61,62,63,64,65,66,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,67,68,69,70,71,72,73,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,74,75,76,77,78,79,80,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,81,82,83,84,85,86,87,88,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,89,90,91,92,93,94,0,95,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,96,97,98,99,100,101,102,103,0,0,0};
    int tmp_idx_truth1D[pi::true_nbins3D]={1,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,3,4,0,0,0,0,0,0,0,0,0,0,0,0,5,6,0,0,0,0,0,0,0,0,0,0,0,7,8,9,0,0,0,0,0,0,0,0,0,0,10,11,12,13,0,0,0,0,0,0,0,0,0,14,15,16,17,18,0,0,0,0,0,0,0,19,20,21,22,23,24,25,0,0,0,0,0,0,0,26,27,28,29,30,31,32,0,0,0,0,0,33,34,35,36,37,38,39,40,41,0,0,0,0,0,42,43,44,45,46,47,48,49,50,0,0,0,0,51,52,53,54,55,56,57,58,0,59,0,0,60,61,62,63,64,65,66,67,68,69,70,71,0,0,0,0,0,0,0,0,0,0,0,0,0,0,72,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,73,74,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,75,76,77,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,78,79,80,81,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,82,83,84,85,86,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,87,88,89,90,91,92,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,93,94,95,96,97,98,99,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,100,101,102,103,104,105,106,107,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,108,109,110,111,112,113,114,115,116,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,117,118,119,120,121,122,123,124,125,126,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,127,128,129,130,131,132,133,134,135,136,137,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,138,139,140,141,142,143,144,145,146,147,148,149};
    for (int i=0; i<pi::reco_nbins3D; ++i) {
      idx_meas1D[i] = tmp_idx_meas1D[i];
      if (idx_meas1D[i] != 0) ++nmeas_3D;
    }
    for (int i=0; i<pi::true_nbins3D; ++i) {
      idx_truth1D_eff[i] = tmp_idx_truth1D[i];
      if (idx_truth1D_eff[i] != 0) ++ntruth_3D_eff;
    }
  }
  //cout<<"### idx_truth1D (separate eff)"<<endl;
  TH1D *htruth_1D = (TH1D*)response_SliceID_1D->Htruth();
  TH1D *htruth_1D_eff = (TH1D*)response_SliceID_1D_eff->Htruth();
  double eff1D[ntruth_3D_eff];
  for (int i=1;i<=ntruth_3D_eff;++i) {
    if (htruth_1D_eff->GetBinContent(i) != 0)
      eff1D[i-1] = htruth_1D->GetBinContent(i)/htruth_1D_eff->GetBinContent(i);
    else eff1D[i-1] = 0;
  }
  //eff1D = {0.00986798,0,0,0.0222041,0,0.0191202,0.0521255,0,0,0.110889,0.0254234,0.0181585,0.0351556,0.0782013,0.0771751,0.0399252,0,0,0,0.13343,0.06179,0.096454,0,0,0.201811,0.140998,0.101281,0.0720633,0,0,0,0,0.0878781,0.0579772,0.0946022,0.0351521,0.0766265,0.178665,0,0,0,0.134462,0.187227,0.135899,0.121898,0.126813,0,0,0,0,0.230769,0.666667,0.210014,0.5,0.59481,0.161341,0.588235,0.688751,0.598386,0.168825,0.6,0.71612,0.700204,0.597549,0.177631,0.571429,0.71505,0.722723,0.668319,0.554158,0.127931,0,0.756531,0.7524,0.695452,0.602508,0.337168,0.0308479,0.8,0.718756,0.716518,0.683866,0.642291,0.464901,0.112423,0.00976884,0,0.629916,0.726249,0.719209,0.668244,0.431872,0.13628,0.0409745,0.0634207,0,0.118215,0.225619,0.308048,0.356286,0.404898,0.12399,0.0247327,0.059936,0.109793};
  for (int i=0; i<pi::true_nbins; ++i)
    for (int j=0; j<pi::true_nbins; ++j)
      for (int k=0; k<pi::true_nbins; ++k) {
        int idx = i*pow(pi::true_nbins,2) + j*pi::true_nbins + k;
        idx_truth1D[idx] = 0;
        if (idx_truth1D_eff[idx] != 0) {
          if (eff1D[idx_truth1D_eff[idx]-1]!=0 && eff1D[idx_truth1D_eff[idx]-1]*htruth_3D->GetBinContent(k+1, j+1, i+1)>0) { // edit here for 1D unfolding (bins for 1D unfolding)
            ++ntruth_3D;
            idx_truth1D[idx] = ntruth_3D;
          }
        }
        //cout<<idx_truth1D[idx]<<",";
      }
  cout<<endl;
  /*int nmeas_dt = 0;
  int idx_meas1D_dt[pi::reco_nbins][pi::reco_nbins][pi::reco_nbins];
  for (int i=0; i<pi::reco_nbins; ++i)
    for (int j=0; j<pi::reco_nbins; ++j)
      for (int k=0; k<pi::reco_nbins; ++k) {
        if (hsig3D->GetBinContent(k+1, j+1, i+1)!=0) {
          ++nmeas_dt;
          idx_meas1D_dt[k+1][j+1][i+1] = nmeas_dt;
          cout<<hmeas_3D->GetBinContent(k+1, j+1, i+1)<<"\t"<<hsig3D->GetBinContent(k+1, j+1, i+1)<<endl;
        }
      }
  cout<<nmeas_dt<<endl; // 88
  cout<<nmeas_3D<<endl; // 86
  cout<<ntruth_3D<<endl; // 106
  return 0;*/
  TH1D *hmeas_1D = (TH1D*)response_SliceID_1D->Hmeasured();
  TH1D* hsig1D = new TH1D("hsig1D", "hsig1D", nmeas_3D, 0, nmeas_3D);
  for (int i=0; i<pi::reco_nbins; ++i)
    for (int j=0; j<pi::reco_nbins; ++j)
      for (int k=0; k<pi::reco_nbins; ++k) {
        int idx = i*pow(pi::reco_nbins,2) + j*pi::reco_nbins + k;
        if (hsig3D->GetBinContent(k+1, j+1, i+1)!=0) {
          if (idx_meas1D[idx] != 0) {
            hsig1D->SetBinContent(idx_meas1D[idx], hsig3D->GetBinContent(k+1, j+1, i+1));
            hsig1D->SetBinError(idx_meas1D[idx], hdata_3D->GetBinError(k+1, j+1, i+1)); // don't include bkg stat error
          }
          else {
            cout<<"Warning: this bin in data is emtpy in MC: (hmeas_3D, hsig3D, hdata_3D, hmc_3D)  "<<hmeas_3D->GetBinContent(k+1, j+1, i+1)<<", "<<hsig3D->GetBinContent(k+1, j+1, i+1)<<", "<<hdata_3D->GetBinContent(k+1, j+1, i+1)<<", "<<hmc_3D->GetBinContent(k+1, j+1, i+1)<<endl;
          }
        }
        else {
          if (idx_meas1D[idx] != 0) cout<<"Warning: this bin in MC is emtpy in data: (hmeas_3D, hsig3D)  "<<hmeas_3D->GetBinContent(k+1, j+1, i+1)<<", "<<hsig3D->GetBinContent(k+1, j+1, i+1)<<endl;
        }
      }
  ///END  get the 1D index and the hsig1D histogram for q3D unfolding  //////////
  
  
  //////////  get outputs of unfolding  //////////
  RooUnfoldBayes unfold_3D (response_SliceID_3D, hsig3D, 0);
  RooUnfoldResponse *response_SliceID_Inc = (RooUnfoldResponse*)fmc->Get("response_SliceID_Inc");
  RooUnfoldBayes unfold_Inc (response_SliceID_Inc, hsiginc, 10);
  RooUnfoldResponse *response_SliceID_Int = (RooUnfoldResponse*)fmc->Get("response_SliceID_Int");
  RooUnfoldBayes unfold_Int (response_SliceID_Int, hsignal, 10);
  RooUnfoldResponse *response_SliceID_Ini = (RooUnfoldResponse*)fmc->Get("response_SliceID_Ini");
  RooUnfoldBayes unfold_Ini (response_SliceID_Ini, hsigini, 10);
  
  RooUnfoldBayes unfold_1D (response_SliceID_1D, hsig1D, 4);
  /// to determine the number of iterations in unfolding using toys
  RooUnfoldResponse *response_SliceID_1D_1 = (RooUnfoldResponse*)fmc->Get("response_SliceID_1D_noeff");
  RooUnfoldBayes unfold_1D_1 (response_SliceID_1D_1, hsig1D, 1);
  RooUnfoldResponse *response_SliceID_1D_2 = (RooUnfoldResponse*)fmc->Get("response_SliceID_1D_noeff");
  RooUnfoldBayes unfold_1D_2 (response_SliceID_1D_2, hsig1D, 2);
  RooUnfoldResponse *response_SliceID_1D_3 = (RooUnfoldResponse*)fmc->Get("response_SliceID_1D_noeff");
  RooUnfoldBayes unfold_1D_3 (response_SliceID_1D_3, hsig1D, 3);
  RooUnfoldResponse *response_SliceID_1D_4 = (RooUnfoldResponse*)fmc->Get("response_SliceID_1D_noeff");
  RooUnfoldBayes unfold_1D_4 (response_SliceID_1D_4, hsig1D, 4);
  RooUnfoldResponse *response_SliceID_1D_5 = (RooUnfoldResponse*)fmc->Get("response_SliceID_1D_noeff");
  RooUnfoldBayes unfold_1D_5 (response_SliceID_1D_5, hsig1D, 5);
  RooUnfoldResponse *response_SliceID_1D_6 = (RooUnfoldResponse*)fmc->Get("response_SliceID_1D_noeff");
  RooUnfoldBayes unfold_1D_6 (response_SliceID_1D_6, hsig1D, 6);
  RooUnfoldResponse *response_SliceID_1D_7 = (RooUnfoldResponse*)fmc->Get("response_SliceID_1D_noeff");
  RooUnfoldBayes unfold_1D_7 (response_SliceID_1D_7, hsig1D, 7);
  RooUnfoldResponse *response_SliceID_1D_8 = (RooUnfoldResponse*)fmc->Get("response_SliceID_1D_noeff");
  RooUnfoldBayes unfold_1D_8 (response_SliceID_1D_8, hsig1D, 8);
  
  
  //hsigini = (TH1D*)hsig3D->Project3D("x");
  //hsiginc = (TH1D*)hsig3D->Project3D("y");
  //hsignal = (TH1D*)hsig3D->Project3D("z");
  hsigini->SetNameTitle("hsigini","All pion initial;Slice ID;Events");
  hsiginc->SetNameTitle("hsiginc","All pion incident;Slice ID;Events");
  hsignal->SetNameTitle("hsignal","Pion interaction signal;Slice ID;Events");
  
  double Ndata = hsig3D->Integral();
  double Nmc = hmeas_3D->Integral();
  double Eweight = Ndata/Nmc;

  /*TH2D *hresponse_3D = (TH2D*)response_SliceID_3D->Hresponse();
  TMatrixD mR_3D(pow(pi::reco_nbins, 3), pow(pi::reco_nbins, 3));
  for (int i=0; i<pi::reco_nbins; i++) {
    for (int j=0; j<pi::reco_nbins; j++) {
      for (int k=0; k<pi::reco_nbins; k++) {
        double factor = 1/htruth_3D->GetBinContent(k+1, j+1, i+1); // to normalize R matrix
        for (int l=0; l<pow(pi::reco_nbins, 3); l++) {
          int idx = i*pow(pi::reco_nbins,2) + j*pi::reco_nbins + k;
          double tmp = hresponse_3D->GetBinContent(l+1, idx+1);
          if (tmp != 0)
            mR_3D(l, idx) = tmp*factor;
        }
      }
    }
  }
  TVectorD vmeas_3D(pow(pi::reco_nbins, 3));
  vmeas_3D = mR_3D * vcentral_truth;
  if (false) { // toy
    ofstream myfile_sys;
    myfile_sys.open ("../../thinslice_sys_nominal/build/toys/syscov_nominal_0115.txt", ios::app); // input toys
    //3D
    for (int i=0; i<pi::reco_nbins; ++i)
      for (int j=0; j<pi::reco_nbins; ++j)
        for (int k=0; k<pi::reco_nbins; ++k) {
          int idx = i*pow(pi::reco_nbins,2) + j*pi::reco_nbins + k;
          myfile_sys<<vmeas_3D(idx) - hsig3D->GetBinContent(k+1, j+1, i+1)<<"\t";
          //cout<<hsig3D->GetBinContent(k+1, j+1, i+1)<<",";
        }
    myfile_sys<<endl<<endl;
    myfile_sys.close();
    
    /*fout->Write();
     fout->Close();
     return 0;
  }*/
  
  TH1D *hval_siginc_reco = (TH1D*)fdata->Get("h_recosliceid_pion_cuts");
  TH1D *hval_signal_reco = (TH1D*)fdata->Get("h_recosliceid_pioninelastic_cuts");
  TH1D *hval_sigini_reco = (TH1D*)fdata->Get("h_recoinisliceid_pion_cuts");
  TH1D *hval_trueinc = (TH1D*)fdata->Get("h_truesliceid_pion_all");
  TH1D *hval_trueint = (TH1D*)fdata->Get("h_truesliceid_pioninelastic_all");
  TH1D *hval_trueini = (TH1D*)fdata->Get("h_trueinisliceid_pion_all");
  TH3D *hval_true3D = (TH3D*)fdata->Get("h_true3Dsliceid_pion_all");
  double Nfakedata = hval_siginc_reco->Integral();
  cout<<"### Nsig (Ndata-Nbkg): "<<Ndata<<"; NtruthMC: "<<Nmc<<"; Nfakedata_sig (judge by truth): "<<Nfakedata<<endl;
  //hval_siginc_reco->Scale(Ndata/Nfakedata);
  //hval_signal_reco->Scale(Ndata/Nfakedata);
  //hval_sigini_reco->Scale(Ndata/Nfakedata);
  //hval_trueinc->Scale(Ndata/Nfakedata);
  //hval_trueint->Scale(Ndata/Nfakedata);
  //hval_trueini->Scale(Ndata/Nfakedata);

  // include MC stat error
  bool include_MCstat = false;
  if (include_MCstat) {
    TMatrixD mcov_3D_stat(pow(pi::reco_nbins,3), pow(pi::reco_nbins,3));
    TMatrixD mcov_3D(pow(pi::reco_nbins,3), pow(pi::reco_nbins,3));
    mcov_3D_stat = unfold_3D.GetMeasuredCov();
    TVectorD Emeas_MC = response_SliceID_3D->Emeasured();
    for(int i=0; i<pow(pi::reco_nbins,3); i++) {
      mcov_3D(i, i) = mcov_3D_stat(i, i) + pow(Emeas_MC(i)*Eweight, 2);
      //if (mcov_3D(i, i)!=0)
      //  cout<<mcov_3D_stat(i, i)<<"\t"<<mcov_3D(i, i)<<endl;
      unfold_3D.SetMeasuredCov(mcov_3D);
    }
  }
  /*if (false) { // has_cov_input
    FILE *fcov_3D=fopen("../../thinslice_sys_nominal/build/toys/cov_nominal_0115.txt","r"); // input covariance
    if (!fcov_3D) {
      cout<<"cov_3D_input not found!"<<endl;
      return 1;
    }
    TMatrixD mcov_3D_input(pow(pi::reco_nbins, 3), pow(pi::reco_nbins, 3));
    double vv;
    for(int i=0; i<pow(pi::reco_nbins, 3); i++) {
      for(int j=0; j<pow(pi::reco_nbins, 3); j++) {
        fscanf(fcov_3D, "%lf", &vv);
        mcov_3D_input(i, j) = vv;
      }
    }
    //mcov_3D = mcov_3D_input;
    unfold_3D.SetMeasuredCov(mcov_3D_input);
  }*/
  
  TH3D *hsig3D_uf;
  TH1D *hsiginc_uf;
  TH1D *hsignal_uf;
  TH1D *hsigini_uf;
  TH1D *hsig1D_uf;
  
  for (int i=0; i<true_nbins3D; ++i) {
    if (idx_truth1D_eff[i]!=0) {
      if (idx_truth1D[i]==0) {
        if (eff1D[idx_truth1D_eff[i]-1] != 0) {
          cout<<"*** Check here "<<idx_truth1D_eff[i]<<endl;
          eff1D[idx_truth1D_eff[i]-1] = 0;
        }
      }
    }
  }
  hsig1D_uf = (TH1D*)unfold_1D.Hreco();
  /// to determine the number of iterations in unfolding using toys
  bool toy_study = false;
  TH1D *hsig1D_uf_1;
  TH1D *hsig1D_uf_2;
  TH1D *hsig1D_uf_3;
  TH1D *hsig1D_uf_4;
  TH1D *hsig1D_uf_5;
  TH1D *hsig1D_uf_6;
  TH1D *hsig1D_uf_7;
  TH1D *hsig1D_uf_8;
  if (toy_study) {
    hsig1D_uf_1 = (TH1D*)unfold_1D_1.Hreco();
    hsig1D_uf_2 = (TH1D*)unfold_1D_2.Hreco();
    hsig1D_uf_3 = (TH1D*)unfold_1D_3.Hreco();
    hsig1D_uf_4 = (TH1D*)unfold_1D_4.Hreco();
    hsig1D_uf_5 = (TH1D*)unfold_1D_5.Hreco();
    hsig1D_uf_6 = (TH1D*)unfold_1D_6.Hreco();
    hsig1D_uf_7 = (TH1D*)unfold_1D_7.Hreco();
    hsig1D_uf_8 = (TH1D*)unfold_1D_8.Hreco();
  }
  /// define input and output histograms for external unfolding methods
  TH1D* h1measdata_1D = new TH1D("h1measdata_1D", "h1measdata_1D", nmeas_3D, 0, nmeas_3D);
  TH1D* h1unfdata_1D = new TH1D("h1unfdata_1D", "h1unfdata_1D", ntruth_3D_eff, 0, ntruth_3D_eff);
  TH1D* h1measMC_1D = new TH1D("h1measMC_1D", "h1measMC_1D", nmeas_3D, 0, nmeas_3D);
  TH1D* h1truthMC_1D = new TH1D("h1truthMC_1D", "h1truthMC_1D", ntruth_3D_eff, 0, ntruth_3D_eff);
  TH1D* h1seltruthMC_1D = new TH1D("h1seltruthMC_1D", "h1seltruthMC_1D", ntruth_3D_eff, 0, ntruth_3D_eff);
  for (int i=0; i<nmeas_3D; ++i) {
    h1measMC_1D->SetBinContent(i+1, response_SliceID_1D->Hmeasured()->GetBinContent(i+1));
    h1measdata_1D->SetBinContent(i+1, hsig1D->GetBinContent(i+1));
  }
  for (int i=0; i<ntruth_3D_eff; ++i) {
    h1truthMC_1D->SetBinContent(i+1, response_SliceID_1D_eff->Htruth()->GetBinContent(i+1));
    h1seltruthMC_1D->SetBinContent(i+1, response_SliceID_1D->Htruth()->GetBinContent(i+1));//->Hresponse()->ProjectionY()->GetBinContent(i+1));
    h1unfdata_1D->SetBinContent(i+1, hsig1D_uf->GetBinContent(i+1));
  }
  h1measdata_1D->Write("h1measdata_1D");
  h1unfdata_1D->Write("h1unfdata_1D");
  h1measMC_1D->Write("h1measMC_1D");
  h1truthMC_1D->Write("h1truthMC_1D");
  h1seltruthMC_1D->Write("h1seltruthMC_1D");
  TH2D* covinput_1D = new TH2D(unfold_1D.GetMeasuredCov()); // input covariance matrix (should be diagonal)
  covinput_1D->Write("covinput_1D");
  TMatrixD cov_matrix_1D = unfold_1D.Ereco();
  TH1D *heff_1D = new TH1D("heff_1D", "heff_1D", ntruth_3D_eff, 0, ntruth_3D_eff);
  for (int i=0; i<ntruth_3D_eff; ++i) {
    heff_1D->SetBinContent(i+1, eff1D[i]);
    if (eff1D[i] != 0) {
      hsig1D_uf->SetBinContent(i+1, hsig1D_uf->GetBinContent(i+1)/eff1D[i]);
      double bine = hsig1D_uf->GetBinError(i+1);
      hsig1D_uf->SetBinError(i+1, bine/eff1D[i]); // should estimate efficiency uncertainty using Clopper-Pearson (as a systematic)
    } // h1seltruthMC_1D : h1truthMC_1D == h1unfdata_1D : hsig1D_uf
    else {
      hsig1D_uf->SetBinContent(i+1, response_SliceID_1D_eff->Htruth()->GetBinContent(i+1)*Eweight);
      hsig1D_uf->SetBinError(i+1, response_SliceID_1D_eff->Htruth()->GetBinError(i+1)*Eweight);
    }
  }
  hsig1D_uf->SetNameTitle("hsig1D_uf", "Unfolded 1D signal;Slice ID;Events");
  ///END  get outputs of unfolding  //////////

  
  //////////  save toy unfolded histograms to txt for systematic evaluations  //////////
  if (toy_study) {
    for (int i=0; i<ntruth_3D_eff; ++i) {
      if (eff1D[i] != 0) {
        double bine;
        hsig1D_uf_1->SetBinContent(i+1, hsig1D_uf_1->GetBinContent(i+1)/eff1D[i]);
        bine = hsig1D_uf_1->GetBinError(i+1);
        hsig1D_uf_1->SetBinError(i+1, bine/eff1D[i]);
        
        hsig1D_uf_2->SetBinContent(i+1, hsig1D_uf_2->GetBinContent(i+1)/eff1D[i]);
        bine = hsig1D_uf_2->GetBinError(i+1);
        hsig1D_uf_2->SetBinError(i+1, bine/eff1D[i]);
        
        hsig1D_uf_3->SetBinContent(i+1, hsig1D_uf_3->GetBinContent(i+1)/eff1D[i]);
        bine = hsig1D_uf_3->GetBinError(i+1);
        hsig1D_uf_3->SetBinError(i+1, bine/eff1D[i]);
        
        hsig1D_uf_4->SetBinContent(i+1, hsig1D_uf_4->GetBinContent(i+1)/eff1D[i]);
        bine = hsig1D_uf_4->GetBinError(i+1);
        hsig1D_uf_4->SetBinError(i+1, bine/eff1D[i]);
        
        hsig1D_uf_5->SetBinContent(i+1, hsig1D_uf_5->GetBinContent(i+1)/eff1D[i]);
        bine = hsig1D_uf_5->GetBinError(i+1);
        hsig1D_uf_5->SetBinError(i+1, bine/eff1D[i]);
        
        hsig1D_uf_6->SetBinContent(i+1, hsig1D_uf_6->GetBinContent(i+1)/eff1D[i]);
        bine = hsig1D_uf_6->GetBinError(i+1);
        hsig1D_uf_6->SetBinError(i+1, bine/eff1D[i]);
        
        hsig1D_uf_7->SetBinContent(i+1, hsig1D_uf_7->GetBinContent(i+1)/eff1D[i]);
        bine = hsig1D_uf_7->GetBinError(i+1);
        hsig1D_uf_7->SetBinError(i+1, bine/eff1D[i]);
        
        hsig1D_uf_8->SetBinContent(i+1, hsig1D_uf_8->GetBinContent(i+1)/eff1D[i]);
        bine = hsig1D_uf_8->GetBinError(i+1);
        hsig1D_uf_8->SetBinError(i+1, bine/eff1D[i]);
      }
      else {
        double binc = response_SliceID_1D_eff->Htruth()->GetBinContent(i+1)*Eweight;
        double bine = response_SliceID_1D_eff->Htruth()->GetBinError(i+1)*Eweight;
        hsig1D_uf_1->SetBinContent(i+1, binc);
        hsig1D_uf_1->SetBinError(i+1, bine);
        
        hsig1D_uf_2->SetBinContent(i+1, binc);
        hsig1D_uf_2->SetBinError(i+1, bine);
        
        hsig1D_uf_3->SetBinContent(i+1, binc);
        hsig1D_uf_3->SetBinError(i+1, bine);
        
        hsig1D_uf_4->SetBinContent(i+1, binc);
        hsig1D_uf_4->SetBinError(i+1, bine);
        
        hsig1D_uf_5->SetBinContent(i+1, binc);
        hsig1D_uf_5->SetBinError(i+1, bine);
        
        hsig1D_uf_6->SetBinContent(i+1, binc);
        hsig1D_uf_6->SetBinError(i+1, bine);
        
        hsig1D_uf_7->SetBinContent(i+1, binc);
        hsig1D_uf_7->SetBinError(i+1, bine);
        
        hsig1D_uf_8->SetBinContent(i+1, binc);
        hsig1D_uf_8->SetBinError(i+1, bine);
      }
    }
    
    const char feat[20] = "MCstat_0430";
    ofstream myfile_sys_unfold_1;
    myfile_sys_unfold_1.open (Form("../../thinslice_sys_nominal/build/toys/hsig1Duf_%s_iter1.txt",feat), ios::app); //output hsig1D_uf of each toy
    for (int i=0; i<ntruth_3D_eff; ++i) {
      myfile_sys_unfold_1<<hsig1D_uf_1->GetBinContent(i+1)<<"\t";
    }
    myfile_sys_unfold_1<<endl<<endl;
    myfile_sys_unfold_1.close();
    
    ofstream myfile_sys_unfold_2;
    myfile_sys_unfold_2.open (Form("../../thinslice_sys_nominal/build/toys/hsig1Duf_%s_iter2.txt",feat), ios::app); //output hsig1D_uf of each toy
    for (int i=0; i<ntruth_3D_eff; ++i) {
      myfile_sys_unfold_2<<hsig1D_uf_2->GetBinContent(i+1)<<"\t";
    }
    myfile_sys_unfold_2<<endl<<endl;
    myfile_sys_unfold_2.close();
    
    ofstream myfile_sys_unfold_3;
    myfile_sys_unfold_3.open (Form("../../thinslice_sys_nominal/build/toys/hsig1Duf_%s_iter3.txt",feat), ios::app); //output hsig1D_uf of each toy
    for (int i=0; i<ntruth_3D_eff; ++i) {
      myfile_sys_unfold_3<<hsig1D_uf_3->GetBinContent(i+1)<<"\t";
    }
    myfile_sys_unfold_3<<endl<<endl;
    myfile_sys_unfold_3.close();
    
    ofstream myfile_sys_unfold_4;
    myfile_sys_unfold_4.open (Form("../../thinslice_sys_nominal/build/toys/hsig1Duf_%s_iter4.txt",feat), ios::app); //output hsig1D_uf of each toy
    for (int i=0; i<ntruth_3D_eff; ++i) {
      myfile_sys_unfold_4<<hsig1D_uf_4->GetBinContent(i+1)<<"\t";
    }
    myfile_sys_unfold_4<<endl<<endl;
    myfile_sys_unfold_4.close();
    
    ofstream myfile_sys_unfold_5;
    myfile_sys_unfold_5.open (Form("../../thinslice_sys_nominal/build/toys/hsig1Duf_%s_iter5.txt",feat), ios::app); //output hsig1D_uf of each toy
    for (int i=0; i<ntruth_3D_eff; ++i) {
      myfile_sys_unfold_5<<hsig1D_uf_5->GetBinContent(i+1)<<"\t";
    }
    myfile_sys_unfold_5<<endl<<endl;
    myfile_sys_unfold_5.close();
    
    ofstream myfile_sys_unfold_6;
    myfile_sys_unfold_6.open (Form("../../thinslice_sys_nominal/build/toys/hsig1Duf_%s_iter6.txt",feat), ios::app); //output hsig1D_uf of each toy
    for (int i=0; i<ntruth_3D_eff; ++i) {
      myfile_sys_unfold_6<<hsig1D_uf_6->GetBinContent(i+1)<<"\t";
    }
    myfile_sys_unfold_6<<endl<<endl;
    myfile_sys_unfold_6.close();
    
    ofstream myfile_sys_unfold_7;
    myfile_sys_unfold_7.open (Form("../../thinslice_sys_nominal/build/toys/hsig1Duf_%s_iter7.txt",feat), ios::app); //output hsig1D_uf of each toy
    for (int i=0; i<ntruth_3D_eff; ++i) {
      myfile_sys_unfold_7<<hsig1D_uf_7->GetBinContent(i+1)<<"\t";
    }
    myfile_sys_unfold_7<<endl<<endl;
    myfile_sys_unfold_7.close();
    
    ofstream myfile_sys_unfold_8;
    myfile_sys_unfold_8.open (Form("../../thinslice_sys_nominal/build/toys/hsig1Duf_%s_iter8.txt",feat), ios::app); //output hsig1D_uf of each toy
    for (int i=0; i<ntruth_3D_eff; ++i) {
      myfile_sys_unfold_8<<hsig1D_uf_8->GetBinContent(i+1)<<"\t";
    }
    myfile_sys_unfold_8<<endl<<endl;
    myfile_sys_unfold_8.close();
    
    fout->Write();
    fout->Close();
    return 0;
  }
  
  //hsig3D_uf = (TH3D*)unfold_3D.Hreco(); // 3D unfolding (can be time-consuming)
  //hsig3D_uf->SetNameTitle("hsig3D_uf", "Unfolded 3D signal;Slice ID;Events");
  hsig3D_uf = new TH3D("hsig3D_uf","Unfolded 3D signal;Slice ID;Events",pi::true_nbins,pi::true_bins,pi::true_nbins,pi::true_bins,pi::true_nbins,pi::true_bins);
  // 1D index back to 3D
  for (int i=0; i<pi::true_nbins; ++i)
    for (int j=0; j<pi::true_nbins; ++j)
      for (int k=0; k<pi::true_nbins; ++k) {
        int idx = i*pow(pi::true_nbins,2) + j*pi::true_nbins + k;
        hsig3D_uf->SetBinContent(k+1, j+1, i+1, 0);
        if (idx_truth1D_eff[idx] != 0) {
          hsig3D_uf->SetBinContent(k+1, j+1, i+1, hsig1D_uf->GetBinContent(idx_truth1D_eff[idx]));
          hsig3D_uf->SetBinError(k+1, j+1, i+1, hsig1D_uf->GetBinError(idx_truth1D_eff[idx])); // should be the same with cov_matrix_3D below
        }
      }
  /// save the 3D unfolded histogram as Unfold3DTable.cxx
  /*std::filebuf fb;
  fb.open(root["UnfoldTable"].asString().c_str(),std::ios::out);
  std::ostream os(&fb);
  unfold_3D.PrintTable(os);
  fb.close();*/
  hsiginc_uf = (TH1D*)unfold_Inc.Hreco();
  hsignal_uf = (TH1D*)unfold_Int.Hreco();
  hsigini_uf = (TH1D*)unfold_Ini.Hreco();
  cout<<"hsiginc_uf: "<<hsiginc_uf->Integral()<<endl;
  cout<<"hsigini_uf: "<<hsigini_uf->Integral()<<endl;
  //hsigini_uf->Scale(hsiginc_uf->Integral()/hsigini_uf->Integral());
  hsigini_uf = (TH1D*)hsig3D_uf->Project3D("x");
  hsiginc_uf = (TH1D*)hsig3D_uf->Project3D("y");
  hsignal_uf = (TH1D*)hsig3D_uf->Project3D("z");
  hsigini_uf->SetNameTitle("hsigini_uf","Unfolded initial signal;Slice ID;Events");
  hsiginc_uf->SetNameTitle("hsiginc_uf","Unfolded incident signal;Slice ID;Events");
  hsignal_uf->SetNameTitle("hsignal_uf","Unfolded interaction signal;Slice ID;Events");
  cout<<"hsiginc_uf: "<<hsiginc_uf->Integral()<<endl;
  cout<<"hsigini_uf: "<<hsigini_uf->Integral()<<endl;
  /*if (toy_study) {
    ofstream myfile_sys_unfold;
    myfile_sys_unfold.open ("../../thinslice_sys_nominal/build/toys/syscov_RooUnfold_nominal_0426.txt", ios::app); //output toys
    //3D
    for (int i=0; i<pi::true_nbins; ++i)
      for (int j=0; j<pi::true_nbins; ++j)
        for (int k=0; k<pi::true_nbins; ++k) {
          myfile_sys_unfold<<hsig3D_uf->GetBinContent(k+1, j+1, i+1)<<"\t";
        }
    myfile_sys_unfold<<endl<<endl;
    myfile_sys_unfold.close();

    fout->Write();
    fout->Close();
    return 0;
  }*/
  ///END  save toy unfolded histograms to txt for systematic evaluations  //////////
  
  
  /// 1D covariance matrix (statistical only; will be updated later with external outcov)
  for (int i=0; i<ntruth_3D_eff; ++i)
    for (int j=0; j<ntruth_3D_eff; ++j) {
      if (eff1D[i] != 0 && eff1D[j] != 0) { // non-zero efficiency
        cov_matrix_1D(i, j) /= (eff1D[i]*eff1D[j]);
      }
      else if (i == j) { // zero-eff bins are directly estimated using truthMC normalized to data, and no correlation with other bins considered
        cov_matrix_1D(i, j) = response_SliceID_1D_eff->Htruth()->GetBinError(i+1)*Eweight*Eweight; // simply estimate as Poisson error
      }
      else cov_matrix_1D(i, j) = 0;
    }
  TMatrixD cov_matrix_3D(pi::true_nbins3D, pi::true_nbins3D);
  //cov_matrix_3D = unfold_3D.Ereco();
  
  
  //////////  update using the external covariance matrix after unfolding (derived using toys)  //////////
  bool has_external_cov = false;
  if (has_external_cov) { // has_external_cov after unfolding
    const char feat[20] = "MCstat_0430_iter4";
    FILE *fcov_1D=fopen(Form("../../thinslice_sys_nominal/build/toys/outCov_%s.txt", feat), "r"); //output covariance
    if (!fcov_1D) {
      cout<<"*** outCov file not found!"<<endl;
      return 1;
    }
    cout<<Form("### Using external covariance matrix outCov_%s.txt", feat)<<endl;
    double vv;
    for(int i=0; i<ntruth_3D_eff; i++) {
      for(int j=0; j<ntruth_3D_eff; j++) {
        fscanf(fcov_1D, "%lf", &vv);
        cov_matrix_1D(i, j) += vv;
      }
    }
  }
  ///END  update using the external covariance matrix after unfolding (derived using toys)  //////////
  
  /// define input and output histograms for external unfolding methods (3D)
  /*TH1D* h1measdata = new TH1D("h1measdata", "h1measdata", pi::reco_nbins3D, 0, pi::reco_nbins3D);
  TH1D* h1unfdata = new TH1D("h1unfdata", "h1unfdata", pi::true_nbins3D, 0, pi::true_nbins3D);
  TH1D* h1measMC = new TH1D("h1measMC", "h1measMC", pi::reco_nbins3D, 0, pi::reco_nbins3D);
  TH1D* h1truthMC = new TH1D("h1truthMC", "h1truthMC", pi::true_nbins3D, 0, pi::true_nbins3D);
  for (int i=0; i<pi::reco_nbins; ++i)
    for (int j=0; j<pi::reco_nbins; ++j)
      for (int k=0; k<pi::reco_nbins; ++k) {
        int idx = i*pow(pi::reco_nbins,2) + j*pi::reco_nbins + k;
        h1measMC->SetBinContent(idx+1, hmeas_3D->GetBinContent(k+1, j+1, i+1));
        h1measdata->SetBinContent(idx+1, hsig3D->GetBinContent(k+1, j+1, i+1));
      }
  for (int i=0; i<pi::true_nbins; ++i)
    for (int j=0; j<pi::true_nbins; ++j)
      for (int k=0; k<pi::true_nbins; ++k) {
        int idx = i*pow(pi::true_nbins,2) + j*pi::true_nbins + k;
        h1truthMC->SetBinContent(idx+1, htruth_3D->GetBinContent(k+1, j+1, i+1));
        h1unfdata->SetBinContent(idx+1, hsig3D_uf->GetBinContent(k+1, j+1, i+1));
      }
  h1measdata->Write("h1measdata");
  h1unfdata->Write("h1unfdata");
  h1measMC->Write("h1measMC");
  h1truthMC->Write("h1truthMC");
  TMatrixD mcov_3D_meas(pi::reco_nbins3D, pi::reco_nbins3D);
  mcov_3D_meas = unfold_3D.GetMeasuredCov();
  TH2D* covinput_3D = new TH2D(mcov_3D_meas);
  covinput_3D->Write("covinput_3D");*/
  
  /// output indices for q3D unfolding
  cout<<"##### idx_meas1D[idx]; (k, j, i); MC_meas; data_meas"<<endl;
  for (int i=0; i<pi::reco_nbins; ++i)
    for (int j=0; j<pi::reco_nbins; ++j)
      for (int k=0; k<pi::reco_nbins; ++k) {
        int idx = i*pow(pi::reco_nbins,2) + j*pi::reco_nbins + k;
        if (hmeas_3D->GetBinContent(k+1, j+1, i+1)!=0) { // print out 3D <-> 1D map
          cout<<idx_meas1D[idx]<<"\t"<<k<<"\t"<<j<<"\t"<<i<<"\t\t"<<h1measMC_1D->GetBinContent(idx_meas1D[idx])<<"\t"<<h1measdata_1D->GetBinContent(idx_meas1D[idx])<<endl;
        }
      }
  cout<<endl;
  cout<<"##### idx_truth1D[idx]; (k, j, i); MC_truth; data_unf"<<endl;
  for (int i=0; i<pi::true_nbins; ++i)
    for (int j=0; j<pi::true_nbins; ++j)
      for (int k=0; k<pi::true_nbins; ++k) {
        int idx = i*pow(pi::true_nbins,2) + j*pi::true_nbins + k;
        if (htruth_3D->GetBinContent(k+1, j+1, i+1)!=0) { // print out 3D <-> 1D map
          cout<<idx_truth1D_eff[idx]<<"\t"<<k<<"\t"<<j<<"\t"<<i<<"\t\t"<<h1truthMC_1D->GetBinContent(idx_truth1D_eff[idx])<<"\t"<<hsig1D_uf->GetBinContent(idx_truth1D_eff[idx])<<endl;
        }
      }
  cout<<endl;
  

  //////////  define input and output histograms for external unfolding methods for ini, end, int  //////////
  TH1D* h1measdata_ini = new TH1D("h1measdata_ini", "h1measdata_ini", pi::reco_nbins, 0, pi::reco_nbins);
  TH1D* h1unfdata_ini = new TH1D("h1unfdata_ini", "h1unfdata_ini", pi::true_nbins, 0, pi::true_nbins);
  TH1D* h1measMC_ini = new TH1D("h1measMC_ini", "h1measMC_ini", pi::reco_nbins, 0, pi::reco_nbins);
  TH1D* h1truthMC_ini = new TH1D("h1truthMC_ini", "h1truthMC_ini", pi::true_nbins, 0, pi::true_nbins);
  TH1D* h1measdata_inc = new TH1D("h1measdata_inc", "h1measdata_inc", pi::reco_nbins, 0, pi::reco_nbins);
  TH1D* h1unfdata_inc = new TH1D("h1unfdata_inc", "h1unfdata_inc", pi::true_nbins, 0, pi::true_nbins);
  TH1D* h1measMC_inc = new TH1D("h1measMC_inc", "h1measMC_inc", pi::reco_nbins, 0, pi::reco_nbins);
  TH1D* h1truthMC_inc = new TH1D("h1truthMC_inc", "h1truthMC_inc", pi::true_nbins, 0, pi::true_nbins);
  TH1D* h1measdata_int = new TH1D("h1measdata_int", "h1measdata_int", pi::reco_nbins, 0, pi::reco_nbins);
  TH1D* h1unfdata_int = new TH1D("h1unfdata_int", "h1unfdata_int", pi::true_nbins, 0, pi::true_nbins);
  TH1D* h1measMC_int = new TH1D("h1measMC_int", "h1measMC_int", pi::reco_nbins, 0, pi::reco_nbins);
  TH1D* h1truthMC_int = new TH1D("h1truthMC_int", "h1truthMC_int", pi::true_nbins, 0, pi::true_nbins);
  for (int i=0; i<pi::reco_nbins; ++i) {
    h1measdata_ini->SetBinContent(i+1, hsigini->GetBinContent(i+1));
    h1measMC_ini->SetBinContent(i+1, response_SliceID_Ini->Hmeasured()->GetBinContent(i+1));
    h1measdata_inc->SetBinContent(i+1, hsiginc->GetBinContent(i+1));
    h1measMC_inc->SetBinContent(i+1, response_SliceID_Inc->Hmeasured()->GetBinContent(i+1));
    h1measdata_int->SetBinContent(i+1, hsignal->GetBinContent(i+1));
    h1measMC_int->SetBinContent(i+1, response_SliceID_Int->Hmeasured()->GetBinContent(i+1));
  }
  for (int i=0; i<pi::true_nbins; ++i) {
    h1unfdata_ini->SetBinContent(i+1, hsigini_uf->GetBinContent(i+1));
    h1truthMC_ini->SetBinContent(i+1, response_SliceID_Ini->Htruth()->GetBinContent(i+1));
    h1unfdata_inc->SetBinContent(i+1, hsiginc_uf->GetBinContent(i+1));
    h1truthMC_inc->SetBinContent(i+1, response_SliceID_Inc->Htruth()->GetBinContent(i+1));
    h1unfdata_int->SetBinContent(i+1, hsignal_uf->GetBinContent(i+1));
    h1truthMC_int->SetBinContent(i+1, response_SliceID_Int->Htruth()->GetBinContent(i+1));
  }
  h1measdata_ini->Write("h1measdata_ini");
  h1unfdata_ini->Write("h1unfdata_ini");
  h1measMC_ini->Write("h1measMC_ini");
  h1truthMC_ini->Write("h1truthMC_ini");
  TH2D* covinput_ini = new TH2D(unfold_Ini.GetMeasuredCov());
  covinput_ini->Write("covinput_ini");
  h1measdata_inc->Write("h1measdata_inc");
  h1unfdata_inc->Write("h1unfdata_inc");
  h1measMC_inc->Write("h1measMC_inc");
  h1truthMC_inc->Write("h1truthMC_inc");
  TH2D* covinput_inc = new TH2D(unfold_Inc.GetMeasuredCov());
  covinput_inc->Write("covinput_inc");
  h1measdata_int->Write("h1measdata_int");
  h1unfdata_int->Write("h1unfdata_int");
  h1measMC_int->Write("h1measMC_int");
  h1truthMC_int->Write("h1truthMC_int");
  TH2D* covinput_int = new TH2D(unfold_Int.GetMeasuredCov());
  covinput_int->Write("covinput_int");
  ///END  define input and output histograms for external unfolding methods for ini, end, int  //////////

  
  /// output covariance matrix for q3D variable
  TH2D *covariance_1D = new TH2D(cov_matrix_1D); // output covariance matrix
  covariance_1D->Write("covariance_1D");
  TH2D *correlation_1D = (TH2D*)covariance_1D->Clone();
  vector<double> sigma_1D;
  for (int i=1; i<=covariance_1D->GetNbinsX(); ++i) {
    sigma_1D.push_back(sqrt(covariance_1D->GetBinContent(i,i)));
  }
  for (int i=1; i<=covariance_1D->GetNbinsX(); ++i) {
    for (int j=1; j<=covariance_1D->GetNbinsY(); ++j) {
      if (covariance_1D->GetBinContent(i,j) == 0) correlation_1D->SetBinContent(i, j, 0);
      else correlation_1D->SetBinContent(i, j, covariance_1D->GetBinContent(i,j)/sigma_1D.at(i-1)/sigma_1D.at(j-1));
    }
  }
  correlation_1D->Write("correlation_1D");
  
  /// transform covariance matrix from q3D to 3D
  for (int i=0; i<pi::true_nbins3D; ++i)
    for (int j=0; j<pi::true_nbins3D; ++j) {
      cov_matrix_3D(i, j) = 0;
      if (idx_truth1D_eff[i]!=0 && idx_truth1D_eff[j]!=0) { // non-empty bin
        cov_matrix_3D(i, j) = cov_matrix_1D(idx_truth1D_eff[i]-1, idx_truth1D_eff[j]-1);
      }
    }
  TH2D *covariance_3D = new TH2D(cov_matrix_3D);
  covariance_3D->Write("covariance_3D");
  TH2D *correlation_3D = (TH2D*)covariance_3D->Clone();
  vector<double> sigma_3D;
  for (int i=1; i<=covariance_3D->GetNbinsX(); ++i) {
    sigma_3D.push_back(sqrt(covariance_3D->GetBinContent(i,i)));
  }
  for (int i=1; i<=covariance_3D->GetNbinsX(); ++i) {
    for (int j=1; j<=covariance_3D->GetNbinsY(); ++j) {
      if (covariance_3D->GetBinContent(i,j) == 0) correlation_3D->SetBinContent(i, j, 0);
      else correlation_3D->SetBinContent(i, j, covariance_3D->GetBinContent(i,j)/sigma_3D.at(i-1)/sigma_3D.at(j-1));
    }
  }
  correlation_3D->Write("correlation_3D");
  
  
  //////////  error propagation from N^3 x N^3 to 3N x 3N  //////////
  TMatrixD mhist(3*pi::true_nbins, pow(pi::true_nbins,3));
  for (int bi=0; bi<pow(pi::true_nbins,3); ++bi) {
    int bix = bi%pi::true_nbins;
    int biy = (bi/pi::true_nbins)%pi::true_nbins;
    int biz = (bi/pi::true_nbins/pi::true_nbins)%pi::true_nbins;
    mhist(bix, bi) = 1;
    mhist(pi::true_nbins+biy, bi) = 1;
    mhist(2*pi::true_nbins+biz, bi) = 1;
  }
  TMatrixD Vhist_tmp(3*pi::true_nbins, pow(pi::true_nbins,3));
  TMatrixD Vhist(3*pi::true_nbins, 3*pi::true_nbins);
  Vhist_tmp.Mult(mhist, cov_matrix_3D);
  Vhist.MultT(Vhist_tmp, mhist); // 3N*3N covariance matrix
  TH2D *Vhist_3D = new TH2D(Vhist);
  Vhist_3D->Write("Vhist_3D");
  for (int bi=0; bi<pi::true_nbins; ++bi) {
    int tbi = bi;
    hsigini_uf->SetBinError(bi+1, sqrt(Vhist(tbi,tbi)));
    tbi += pi::true_nbins;
    hsiginc_uf->SetBinError(bi+1, sqrt(Vhist(tbi,tbi)));
    tbi += pi::true_nbins;
    hsignal_uf->SetBinError(bi+1, sqrt(Vhist(tbi,tbi)));
  }
  ///END  error propagation from N^3 x N^3 to 3N x 3N  //////////
  
  
  //////////  output covariance matrix for ini, end, int  //////////
  TMatrixD cov_matrix_inc = unfold_Inc.Ereco();
  TH2D *covariance_inc = new TH2D(cov_matrix_inc);
  covariance_inc->Write("covariance_inc");
  TH2D *correlation_inc = (TH2D*)covariance_inc->Clone();
  vector<double> sigma_inc;
  for (int i=1; i<=covariance_inc->GetNbinsX(); ++i) {
    sigma_inc.push_back(sqrt(covariance_inc->GetBinContent(i,i)));
  }
  for (int i=1; i<=covariance_inc->GetNbinsX(); ++i) {
    for (int j=1; j<=covariance_inc->GetNbinsY(); ++j) {
      if (covariance_inc->GetBinContent(i,j) == 0) correlation_inc->SetBinContent(i, j, 0);
      else correlation_inc->SetBinContent(i, j, covariance_inc->GetBinContent(i,j)/sigma_inc.at(i-1)/sigma_inc.at(j-1));
    }
  }
  correlation_inc->Write("correlation_inc");
  
  TMatrixD cov_matrix_int = unfold_Int.Ereco();
  TH2D *covariance_int = new TH2D(cov_matrix_int);
  covariance_int->Write("covariance_int");
  TH2D *correlation_int = (TH2D*)covariance_int->Clone();
  vector<double> sigma_int;
  for (int i=1; i<=covariance_int->GetNbinsX(); ++i) {
    sigma_int.push_back(sqrt(covariance_int->GetBinContent(i,i)));
  }
  for (int i=1; i<=covariance_int->GetNbinsX(); ++i) {
    for (int j=1; j<=covariance_int->GetNbinsY(); ++j) {
      if (covariance_int->GetBinContent(i,j) == 0) correlation_int->SetBinContent(i, j, 0);
      else correlation_int->SetBinContent(i, j, covariance_int->GetBinContent(i,j)/sigma_int.at(i-1)/sigma_int.at(j-1));
    }
  }
  correlation_int->Write("correlation_int");
  
  TMatrixD cov_matrix_ini = unfold_Ini.Ereco();
  TH2D *covariance_ini = new TH2D(cov_matrix_ini);
  covariance_ini->Write("covariance_ini");
  TH2D *correlation_ini = (TH2D*)covariance_ini->Clone();
  vector<double> sigma_ini;
  for (int i=1; i<=covariance_ini->GetNbinsX(); ++i) {
    sigma_ini.push_back(sqrt(covariance_ini->GetBinContent(i,i)));
  }
  for (int i=1; i<=covariance_ini->GetNbinsX(); ++i) {
    for (int j=1; j<=covariance_ini->GetNbinsY(); ++j) {
      if (covariance_ini->GetBinContent(i,j) == 0) correlation_ini->SetBinContent(i, j, 0);
      else correlation_ini->SetBinContent(i, j, covariance_ini->GetBinContent(i,j)/sigma_ini.at(i-1)/sigma_ini.at(j-1));
    }
  }
  correlation_ini->Write("correlation_ini");
  /*TVectorD covar_diag = unfold_Inc.ErecoV();
  TH1D *cov_diag = new TH1D(covar_diag);
  cov_diag->Write("Mcov_diag");*/
  ///END  output covariance matrix for ini, end, int  //////////
  
  
  //////////  calculate the histograms Ninc, Nint, Nini as well as their covariance matrix  //////////
  double Ninc[pi::true_nbins-1] = {0};
  double Nint[pi::true_nbins-1] = {0};
  double Nini[pi::true_nbins-1] = {0};
  double Nend[pi::true_nbins-1] = {0};
  double err_inc[pi::true_nbins-1] = {0};
  double err_int[pi::true_nbins-1] = {0};
  double err_ini[pi::true_nbins-1] = {0};
  double err_end[pi::true_nbins-1] = {0};
  double SliceID[pi::true_nbins-1] = {0};

  for (int i = 0; i<pi::true_nbins-1; ++i){
    SliceID[i] = i+1;
    Nint[i] = hsignal_uf->GetBinContent(i+2);
    err_int[i] = hsignal_uf->GetBinError(i+2);
    Nini[i] = hsigini_uf->GetBinContent(i+2);
    err_ini[i] = hsigini_uf->GetBinError(i+2);
    Nend[i] = hsiginc_uf->GetBinContent(i+2);
    err_end[i] = hsiginc_uf->GetBinError(i+2);
    for (int j = 0; j<=i; ++j){
      Ninc[i] += hsigini_uf->GetBinContent(j+2);
      //err_inc[i] += pow(hsigini_uf->GetBinError(j+2),2);
    }
    for (int j = 0; j<=i-1; ++j){
      Ninc[i] -= hsiginc_uf->GetBinContent(j+2);
      //err_inc[i] += pow(hsiginc_uf->GetBinError(j+2),2);
    }
    /*for (int j = i; j<pi::true_nbins-1; ++j){
      Ninc[i] += hsiginc_uf->GetBinContent(j+2);
      //err_inc[i] += pow(hsiginc_uf->GetBinError(j+2),2);
    }
    for (int j = i+1; j<pi::true_nbins-1; ++j){
      Ninc[i] -= hsigini_uf->GetBinContent(j+2);
      //err_inc[i] += pow(hsigini_uf->GetBinError(j+2),2);
    }*/
    //err_inc[i] = sqrt(err_inc[i]);
    //err_inc[i] = sqrt(Ninc[i]);
  }
  /// error propagation from (Nini; Nend; Nint) to (Ninc; Nend; Nint)
  TMatrixD mNin(3*(pi::true_nbins-1), 3*pi::true_nbins);
  for (int bi=0; bi<pi::true_nbins-1; ++bi) { // Ninc
    for (int ti=bi+2; ti<pi::true_nbins; ++ti)
      mNin(bi, ti) = -1;
    for (int ti=bi+1+pi::true_nbins; ti<2*pi::true_nbins; ++ti)
      mNin(bi, ti) = 1;
  }
  /*for (int bi=0; bi<pi::true_nbins-1; ++bi) { // Ninc (equivalent)
    for (int ti=1; ti<bi+2; ++ti)
      mNin(bi, ti) = 1;
    for (int ti=1+pi::true_nbins; ti<bi+1+pi::true_nbins; ++ti)
      mNin(bi, ti) = -1;
  }*/
  for (int bi=pi::true_nbins-1; bi<2*(pi::true_nbins-1); ++bi) { // Nend
    mNin(bi, bi+2) = 1;
  }
  for (int bi=2*(pi::true_nbins-1); bi<3*(pi::true_nbins-1); ++bi) { // Nint
    mNin(bi, bi+3) = 1;
  }
  TMatrixD VNin_tmp(3*(pi::true_nbins-1), 3*pi::true_nbins);
  TMatrixD VNin(3*(pi::true_nbins-1), 3*(pi::true_nbins-1));
  VNin_tmp.Mult(mNin, Vhist);
  VNin.MultT(VNin_tmp, mNin); // 3(N-1)*3(N-1) covariance matrix
  TH2D *VNin_3D = new TH2D(VNin);
  VNin_3D->Write("VNin_3D");
  for (int bi=0; bi<pi::true_nbins-1; ++bi) {
    err_inc[bi] = sqrt(VNin(bi,bi));
  }
  
  TGraphErrors *gr_inc = new TGraphErrors(pi::true_nbins-1, SliceID, Ninc, 0, err_inc);
  gr_inc->SetNameTitle("gr_inc", "Incident number;Slice ID;Events");
  gr_inc->Write();
  TGraphErrors *gr_int = new TGraphErrors(pi::true_nbins-1, SliceID, Nint, 0, err_int);
  gr_int->SetNameTitle("gr_int", "Interaction number;Slice ID;Events");
  gr_int->Write();
  TGraphErrors *gr_ini = new TGraphErrors(pi::true_nbins-1, SliceID, Nini, 0, err_ini);
  gr_ini->SetNameTitle("gr_ini", "Initial number;Slice ID;Events");
  gr_ini->Write();
  TGraphErrors *gr_end = new TGraphErrors(pi::true_nbins-1, SliceID, Nend, 0, err_end);
  gr_end->SetNameTitle("gr_end", "Interaction_allpion number;Slice ID;Events");
  gr_end->Write();
  /*TGraphErrors *gr_trueincE = (TGraphErrors*)fmc->Get("gr_trueincE");
  gr_trueincE->SetNameTitle("gr_trueincE", "True incident energy;Slice ID;Energy (MeV)");
  gr_trueincE->Write();
  TGraphErrors *gr_recoincE = (TGraphErrors*)fdata->Get("gr_recoincE");
  gr_recoincE->SetNameTitle("gr_recoincE", "Reco incident energy;Slice ID;Energy (MeV)");
  gr_recoincE->Write();*/
  ///END  calculate the histograms Ninc, Nint, Nini as well as their covariance matrix  //////////
  
  
  //////////  calculate the cross-section as well as its covariance matrix  //////////
  BetheBloch bb(211);
  double NA=6.02214076e23;
  double MAr=39.95; //gmol
  double Density = 1.4; // g/cm^3
  double xs[pi::true_nbins-1] = {0};
  double err_xs[pi::true_nbins-1] = {0};
  double KE[pi::true_nbins-1] = {0};
  double err_KE[pi::true_nbins-1] = {0};
  double dEdx[pi::true_nbins-1] = {0};
  TMatrixD mXS(pi::true_nbins-1, 3*(pi::true_nbins-1)); // Jacobian for (Ninc; Nend; Nint) to XS
  for (int i = 0; i<pi::true_nbins-1; ++i){
    KE[i] = (pi::true_KE[i]+pi::true_KE[i+1])/2;
    err_KE[i] = (pi::true_KE[i]-pi::true_KE[i+1])/2;
    dEdx[i] = bb.meandEdx(KE[i]); // MeV/cm
    double confac = dEdx[i]*MAr/(Density*NA*2*err_KE[i])*1e27;
    xs[i] = confac * Nint[i]/Nend[i] * log(Ninc[i]/(Ninc[i]-Nend[i]));
    //err_xs[i] = dEdx[i]*MAr/(Density*NA*2*err_KE[i])*sqrt(pow(Nint[i]*err_inc[i]/Ninc[i]/(Ninc[i]-Nint[i]),2)+pow(err_int[i]/(Ninc[i]-Nint[i]),2))*1e27;
    mXS(i,i) = confac * Nint[i] / Ninc[i] / (Nend[i]-Ninc[i]); // ∂σ/∂Ninc
    mXS(i,i+pi::true_nbins-1) = confac * Nint[i]/Nend[i] * ( 1/(Ninc[i]-Nend[i]) - 1/Nend[i]*log(Ninc[i]/(Ninc[i]-Nend[i])) ); // ∂σ/∂Nend
    mXS(i,i+2*(pi::true_nbins-1)) = confac * 1 / Nend[i] * log(Ninc[i]/(Ninc[i]-Nend[i])); // ∂σ/∂Nint
  }
  TMatrixD VXS_tmp(pi::true_nbins-1, 3*(pi::true_nbins-1));
  TMatrixD VXS(pi::true_nbins-1, pi::true_nbins-1);
  VXS_tmp.Mult(mXS, VNin);
  VXS.MultT(VXS_tmp, mXS); // (N-1)*(N-1) covariance matrix (minus 1 means excluding the unphysics bin)
  TH2D *VXS_3D = new TH2D(VXS);
  VXS_3D->Write("VXS_3D");
  for (int bi=0; bi<pi::true_nbins-1; ++bi) {
    err_xs[bi] = sqrt(VXS(bi,bi));
  }
  TH2D *corr_matrix_XS = (TH2D*)VXS_3D->Clone(); // correlation matrix
  for (int i=1; i<=VXS_3D->GetNbinsX(); ++i) {
    for (int j=1; j<=VXS_3D->GetNbinsY(); ++j) {
      if (VXS_3D->GetBinContent(i,j) == 0) corr_matrix_XS->SetBinContent(i, j, 0);
      else corr_matrix_XS->SetBinContent(i, j, VXS_3D->GetBinContent(i,j)/err_xs[i-1]/err_xs[j-1]);
    }
  }
  corr_matrix_XS->GetXaxis()->SetTitle("Reco slice ID");
  corr_matrix_XS->GetYaxis()->SetTitle("Reco slice ID");
  /// BEGIN draw correlation matrices option
  /*// (put it in the interactive ROOT)
  const Int_t NRGBs = 3;
  const Int_t NCont = 200;
  Double_t stops[NRGBs] = { 0.0, 0.5, 1.0};
  Double_t red[NRGBs]   = { 0.0, 1.0, 1.0};
  Double_t green[NRGBs]   = { 0.0, 1.0, 0.0};
  Double_t blue[NRGBs]   = { 1.0, 1.0, 0.0};
  TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
  gStyle->SetNumberContours(NCont);
  const Double_t mymin = -1;
  const Double_t mymax = 1;
  Double_t levels[NCont];
  for(int i = 1; i < NCont; i++) {
    levels[i] = mymin + (mymax - mymin) / (NCont - 1) * (i);
  }
  levels[0] = -1;
  corr_matrix_XS->SetContour((sizeof(levels)/sizeof(Double_t)), levels);*/
  /// END draw correlation matrices option
  corr_matrix_XS->Write("corr_matrix_XS");

  TGraphErrors *gr_recoxs = new TGraphErrors(pi::true_nbins-1, KE, xs, err_KE, err_xs);
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
  ///END  calculate the cross-section as well as its covariance matrix  //////////

  
  //////////  calculate truth histograms and cross-sections (similar to procedures for recos above)  //////////
  hval_siginc_reco->Write("hval_siginc_reco");
  hval_signal_reco->Write("hval_signal_reco");
  hval_sigini_reco->Write("hval_sigini_reco");
  hval_trueinc->Write("hval_trueinc");
  hval_trueint->Write("hval_trueint");
  hval_trueini->Write("hval_trueini");
  hval_true3D->Write("hval_true3D");
  TMatrixD cov3D_t(pow(pi::true_nbins,3), pow(pi::true_nbins,3)); // cov3D_t is diagonal
  for (int i=0; i<pi::true_nbins; ++i)
    for (int j=0; j<pi::true_nbins; ++j)
      for (int k=0; k<pi::true_nbins; ++k) {
        int idx = i*pow(pi::true_nbins,2) + j*pi::true_nbins + k;
        cov3D_t(idx,idx) = pow(hval_true3D->GetBinError(k+1, j+1, i+1), 2);
      }
  
  double Ninc_t[pi::true_nbins-1] = {0};
  double Nint_t[pi::true_nbins-1] = {0};
  double Nini_t[pi::true_nbins-1] = {0};
  double Nend_t[pi::true_nbins-1] = {0};
  double err_inc_t[pi::true_nbins-1] = {0};
  double err_int_t[pi::true_nbins-1] = {0};
  double err_ini_t[pi::true_nbins-1] = {0};
  double err_end_t[pi::true_nbins-1] = {0};
  for (int i = 0; i<pi::true_nbins-1; ++i){
    Nint_t[i] = hval_trueint->GetBinContent(i+2);
    err_int_t[i] = hval_trueint->GetBinError(i+2);
    Nini_t[i] = hval_trueini->GetBinContent(i+2);
    err_ini_t[i] = hval_trueini->GetBinError(i+2);
    Nend_t[i] = hval_trueinc->GetBinContent(i+2);
    err_end_t[i] = hval_trueinc->GetBinError(i+2);
    for (int j = 0; j<=i; ++j){
      Ninc_t[i] += hval_trueini->GetBinContent(j+2);
      //err_inc_t[i] += pow(hval_trueini->GetBinError(j+2),2);
    }
    for (int j = 0; j<=i-1; ++j){
      Ninc_t[i] -= hval_trueinc->GetBinContent(j+2);
      //err_inc_t[i] += pow(hval_trueinc->GetBinError(j+2),2);
    }
    /*for (int j = i; j<pi::true_nbins-1; ++j){
      Ninc_t[i] += hval_trueinc->GetBinContent(j+2);
      //err_inc_t[i] += pow(hval_trueinc->GetBinError(j+2),2);
    }
    for (int j = i+1; j<pi::true_nbins-1; ++j){
      Ninc_t[i] -= hval_trueini->GetBinContent(j+2);
      //err_inc_t[i] += pow(hval_trueini->GetBinError(j+2),2);
    }*/
    //err_inc_t[i] = sqrt(err_inc_t[i]);
    //err_inc_t[i] = sqrt(Ninc_t[i]);
  }
  TMatrixD Vhist_t(3*pi::true_nbins, 3*pi::true_nbins);
  /*for (int bi=0; bi<pi::true_nbins; ++bi) { // Vhist_t is diagonal
    int tbi = bi;
    Vhist_t(tbi,tbi) = pow(hval_trueini->GetBinError(bi+1),2);
    tbi += pi::true_nbins;
    Vhist_t(tbi,tbi) = pow(hval_trueinc->GetBinError(bi+1),2);
    tbi += pi::true_nbins;
    Vhist_t(tbi,tbi) = pow(hval_trueint->GetBinError(bi+1),2);
  }*/
  TMatrixD Vhist_tmp_t(3*pi::true_nbins, pow(pi::true_nbins,3));
  Vhist_tmp_t.Mult(mhist, cov3D_t);
  Vhist_t.MultT(Vhist_tmp_t, mhist); // 3N*3N covariance matrix
  TH2D *Vhist_3D_t = new TH2D(Vhist_t);
  Vhist_3D_t->Write("Vhist_3D_t");
  
  TMatrixD VNin_tmp_t(3*(pi::true_nbins-1), 3*pi::true_nbins);
  TMatrixD VNin_t(3*(pi::true_nbins-1), 3*(pi::true_nbins-1));
  VNin_tmp_t.Mult(mNin, Vhist_t);
  VNin_t.MultT(VNin_tmp_t, mNin); // 3(N-1)*3(N-1) covariance matrix
  TH2D *VNin_3D_t = new TH2D(VNin_t);
  VNin_3D_t->Write("VNin_3D_t");
  for (int bi=0; bi<pi::true_nbins-1; ++bi) {
    err_inc_t[bi] = sqrt(VNin_t(bi,bi));
  }
  
  TGraphErrors *gr_inc_t = new TGraphErrors(pi::true_nbins-1, SliceID, Ninc_t, 0, err_inc_t);
  gr_inc_t->SetNameTitle("gr_inc_t", "Incident number;Slice ID;Events");
  gr_inc_t->Write();
  TGraphErrors *gr_int_t = new TGraphErrors(pi::true_nbins-1, SliceID, Nint_t, 0, err_int_t);
  gr_int_t->SetNameTitle("gr_int_t", "Interaction number;Slice ID;Events");
  gr_int_t->Write();
  TGraphErrors *gr_ini_t = new TGraphErrors(pi::true_nbins-1, SliceID, Nini_t, 0, err_ini_t);
  gr_ini_t->SetNameTitle("gr_ini_t", "Initial number;Slice ID;Events");
  gr_ini_t->Write();
  TGraphErrors *gr_end_t = new TGraphErrors(pi::true_nbins-1, SliceID, Nend_t, 0, err_end_t);
  gr_end_t->SetNameTitle("gr_end_t", "Interaction_allpion number;Slice ID;Events");
  gr_end_t->Write();
  double xs_t[pi::true_nbins-1] = {0};
  double err_xs_t[pi::true_nbins-1] = {0};
  TMatrixD mXS_t(pi::true_nbins-1, 3*(pi::true_nbins-1)); // Jacobian for (Ninc; Nend; Nint) to XS
  for (int i = 0; i<pi::true_nbins-1; ++i){
    double confac = dEdx[i]*MAr/(Density*NA*2*err_KE[i])*1e27;
    xs_t[i] = confac * Nint_t[i]/Nend_t[i] * log(Ninc_t[i]/(Ninc_t[i]-Nend_t[i]));
    mXS_t(i,i) = confac * Nint_t[i] / Ninc_t[i] / (Nend_t[i]-Ninc_t[i]); // ∂σ/∂Ninc
    mXS_t(i,i+pi::true_nbins-1) = confac * Nint_t[i]/Nend_t[i] * ( 1/(Ninc_t[i]-Nend_t[i]) - 1/Nend_t[i]*log(Ninc_t[i]/(Ninc_t[i]-Nend_t[i])) ); // ∂σ/∂Nend
    mXS_t(i,i+2*(pi::true_nbins-1)) = confac * 1 / Nend_t[i] * log(Ninc_t[i]/(Ninc_t[i]-Nend_t[i])); // ∂σ/∂Nint
  }
  TMatrixD VXS_tmp_t(pi::true_nbins-1, 3*(pi::true_nbins-1));
  TMatrixD VXS_t(pi::true_nbins-1, pi::true_nbins-1);
  VXS_tmp_t.Mult(mXS_t, VNin_t);
  VXS_t.MultT(VXS_tmp_t, mXS_t); // (N-1)*(N-1) covariance matrix
  TH2D *VXS_3D_t = new TH2D(VXS_t);
  VXS_3D_t->Write("VXS_3D_t");
  for (int bi=0; bi<pi::true_nbins-1; ++bi) {
    err_xs_t[bi] = sqrt(VXS_t(bi,bi));
  }
  TH2D *corr_matrix_XS_t = (TH2D*)VXS_3D_t->Clone(); // correlation matrix
  for (int i=1; i<=VXS_3D_t->GetNbinsX(); ++i) {
    for (int j=1; j<=VXS_3D_t->GetNbinsY(); ++j) {
      if (VXS_3D_t->GetBinContent(i,j) == 0) corr_matrix_XS_t->SetBinContent(i, j, 0);
      else corr_matrix_XS_t->SetBinContent(i, j, VXS_3D_t->GetBinContent(i,j)/err_xs_t[i-1]/err_xs_t[j-1]);
    }
  }
  corr_matrix_XS_t->Write("corr_matrix_XS_t");
  
  TGraphErrors *gr_truexs = new TGraphErrors(pi::true_nbins-1, KE, xs_t, err_KE, err_xs_t);
  gr_truexs->SetNameTitle("gr_truexs", "Reco cross-section;Energy (MeV); Cross-section (mb)");
  gr_truexs->Write();
  TGraphErrors *gr_truexs_allMC = (TGraphErrors*)fdata->Get("gr_truexs");
  gr_truexs_allMC->Write("gr_truexs_allMC");
  ///END  calculate truth histograms and cross-sections (similar to procedures for recos above)  //////////
  
  fout->Write();
  fout->Close();
  return 0;
}
