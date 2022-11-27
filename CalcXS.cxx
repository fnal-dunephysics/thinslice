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

  // background constraint
  TVectorD sf(2); sf[0] = 1.; sf[1] = 0.;
  TVectorD *sf_mu = (TVectorD*)fbkg->Get("sf_mu");
  TVectorD *sf_p = (TVectorD*)fbkg->Get("sf_p");
  TVectorD *sf_spi = (TVectorD*)fbkg->Get("sf_spi");
  //(*sf_mu)[1] = 0;
  //(*sf_p)[1] = 0;
  //(*sf_spi)[1] = 0;

  cout<<"Muon scaling factor: "<<(*sf_mu)[0]<<"+-"<<(*sf_mu)[1]<<endl;
  cout<<"Proton scaling factor: "<<(*sf_p)[0]<<"+-"<<(*sf_p)[1]<<endl;
  cout<<"Pion scaling factor: "<<(*sf_spi)[0]<<"+-"<<(*sf_spi)[1]<<endl;

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
  
  // first loop to get total number of events
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

  // second loop to fill histograms
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
      hsliceID[i]->Scale(hdata->Integral()/hmc->Integral());
      hincsliceID[i]->Scale(hdata_inc->Integral()/hmc_inc->Integral());
      hinisliceID[i]->Scale(hdata_ini->Integral()/hmc_ini->Integral());
      h3DsliceID[i]->Scale(hdata_3D->Integral()/hmc_3D->Integral());
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
          
          for (int k = 0; k < pi::reco_nbins; ++k){
            int bink = hsliceID[i]->FindBin(k-0.5);
            for (int l = 0; l < pi::reco_nbins; ++l){
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

  TH1D *hsignal = (TH1D*)hdata->Clone("hsignal");
  hsignal->Add(hmu,-1);
  hsignal->Add(hproton,-1);
  hsignal->Add(hspi,-1);
  hsignal->Add(hpiel,-1);
  hsignal->Add(hother,-1);
  
  TH1D *hsiginc = (TH1D*)hdata_inc->Clone("hsiginc");
  hsiginc->Add(hmu_inc,-1);
  hsiginc->Add(hproton_inc,-1);
  hsiginc->Add(hspi_inc,-1);
  hsiginc->Add(hother_inc,-1);
  
  TH1D *hsigini = (TH1D*)hdata_ini->Clone("hsigini");
  hsigini->Add(hmu_ini,-1);
  hsigini->Add(hproton_ini,-1);
  hsigini->Add(hspi_ini,-1);
  hsigini->Add(hother_ini,-1);
  
  TH3D *hsig3D = (TH3D*)hdata_3D->Clone("hsig3D");
  hsig3D->SetTitle("All pion 3D;Slice ID;Events");
  hsig3D->Add(hmu_3D,-1);
  hsig3D->Add(hproton_3D,-1);
  hsig3D->Add(hspi_3D,-1);
  hsig3D->Add(hpiel_3D,-1);
  hsig3D->Add(hother_3D,-1);
  
  const int reco_nbins3D = pow(pi::reco_nbins,3);
  //double central_truth[reco_nbins3D] = {12619.3,0,0,0,0,0,0,0,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,3,63.7984,0,0,0,0,0,0,0,0,0,0,47.2254,62.6778,0,0,0,0,0,0,0,0,0,33.9785,40.5839,58.7379,0,0,0,0,0,0,0,0,35.9773,42.7963,51.1829,28.0609,0,0,0,0,0,0,0,24.3006,31.3072,23.2588,24.8591,9.4323,0,0,0,0,0,0,9.33343,31.6074,31.4622,10.4221,6.5242,0.934514,0,0,0,0,0,22.6129,21.5267,27.603,26.2599,8.55322,7.39511,1.37333,0,0,0,3,48.1213,54.0348,39.236,30.6509,14.7971,5.97039,2.2277,3.61514,0,0,6,101.566,148.338,144.86,88.6699,36.9555,15.7319,23.4759,8.51539,4.41719,0,0,0,0,0,0,0,0,0,0,0,0,468,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,54,4168.05,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,54,3237,4324.29,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,51,2485.49,3097.55,3432.5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,30,1749.63,2333.61,2565.43,1549.35,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,21,1410.37,1751.49,1958.13,1133.56,443.661,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9,1006.67,1317.27,1401.12,875.032,312.145,207.383,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,15,722.342,968.469,1070.15,603.55,247.81,137.842,101.188,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,12,1404.98,1731.73,1862.59,1099.73,504.049,241.862,160.591,69.6013,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,18,857.556,1087.4,1196.98,720.545,303.528,150.842,121.297,50.0534,54.6485}; // MC truth
  double central_truth[reco_nbins3D] = {16558.4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,71.6412,0,0,0,0,0,0,0,0,0,0,0,70.3828,0,0,0,0,0,0,0,0,0,38.1555,0,0,0,0,0,0,0,0,0,0,40.4001,48.0574,57.4749,31.5105,0,0,0,0,0,0,0,27.2879,35.1558,26.1181,0,0,0,0,0,0,0,0,0,35.493,35.3299,11.7032,0,0,0,0,0,0,0,25.3927,24.173,30.9963,29.4881,0,0,0,0,0,0,0,54.0369,60.6773,44.0593,34.419,16.6161,6.70433,0,0,0,0,0,114.052,166.573,162.668,99.5702,41.4985,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,897.908,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,61.2597,4372.24,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,32.0605,3297.1,4167.05,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,47.7562,2692.64,3152.71,3567.88,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,25.7251,1889.8,2461.51,2809.74,1772.24,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,20.0316,1643.68,1909.4,2210.23,1293.24,612.961,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1180.83,1560.71,1663.24,952.884,323.231,244.888,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,13.5387,770.538,1094.8,1335.06,768.724,241.049,163.292,79.4288,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1146.18,1815.74,2297.27,1355.97,605,368.408,197.899,47.0622,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,707.788,555.23,1056.82,887.475,334.908,158.523,157.946,140.319,59.9514}; // data unfolded (iterated)
  double central_datasig[reco_nbins3D] = {102.188,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,7.86334,1.35575,0,0,0,0,0,0,0,0,0,4.65634,3.87494,0,0,0,0,0,0,0,0,0,3.8442,3.85163,4.45631,0,0,0,0,0,0,0,0,1.53215,4.60658,3.24738,1.15254,0,0,0,0,0,0,0,1.28249,2.20855,2.18762,2.25887,1.2099,0,0,0,0,0,0,10.4693,5.80001,4.44475,3.31591,1.42079,0,0,0,0,0,0,8.94783,18.2658,12.0404,9.09564,1.29002,0,0,0,0,0,0,2.80312,3.30549,4.36473,5.32057,2.18515,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,310.714,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,147.318,834.746,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,81.9655,2066.59,450.734,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,104.483,1828.19,1778.72,452.705,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,65.5469,1311.37,1681.16,1561.45,244.045,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,39.9205,1031.36,1335.9,1396.95,842.675,96.5756,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,5.79838,848.525,1073.88,1170.54,692.045,208.212,22.0357,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,32.1614,635.665,780.309,860.417,544.473,196.018,20.1658,2.13673,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,6.41122,783.072,1390.98,1579.53,1020.4,392.251,64.8325,9.13673,2.13673,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,12.9656,16.6283,301.131,358.069,167.541,55.1432,11.4102,8.41019,8.82038}; // data - bkg
  TVectorD vcentral_truth(reco_nbins3D);
  for (int i=0; i<pi::reco_nbins; ++i)
    for (int j=0; j<pi::reco_nbins; ++j)
      for (int k=0; k<pi::reco_nbins; ++k) {
        int idx = i*pow(pi::reco_nbins,2) + j*pi::reco_nbins + k;
        vcentral_truth(idx) = central_truth[idx];
        hsig3D->SetBinContent(k+1, j+1, i+1, central_datasig[idx]);
      }
  // use Cholesky decomposition to generate correlated toys
  /*FILE *fcholsqrt_Inc=fopen("/dune/app/users/yinrui/Wiener-SVD-Unfolding/toys/cholsqrt_MCXS_inc_1110.txt","r");
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
  // unfolding
  RooUnfoldResponse *response_SliceID_3D = (RooUnfoldResponse*)fmc->Get("response_SliceID_3D");
  RooUnfoldBayes unfold_3D (response_SliceID_3D, hsig3D, 1);
  RooUnfoldResponse *response_SliceID_Inc = (RooUnfoldResponse*)fmc->Get("response_SliceID_Inc");
  RooUnfoldBayes unfold_Inc (response_SliceID_Inc, hsiginc, 10);
  RooUnfoldResponse *response_SliceID_Int = (RooUnfoldResponse*)fmc->Get("response_SliceID_Int");
  RooUnfoldBayes unfold_Int (response_SliceID_Int, hsignal, 10);
  RooUnfoldResponse *response_SliceID_Ini = (RooUnfoldResponse*)fmc->Get("response_SliceID_Ini");
  RooUnfoldBayes unfold_Ini (response_SliceID_Ini, hsigini, 10);
  
  TH3D *hmeas_MC = (TH3D*)response_SliceID_3D->Hmeasured();
  hsigini = (TH1D*)hsig3D->Project3D("x");
  hsiginc = (TH1D*)hsig3D->Project3D("y");
  hsignal = (TH1D*)hsig3D->Project3D("z");
  hsigini->SetNameTitle("hsigini","All pion initial;Slice ID;Events");
  hsiginc->SetNameTitle("hsiginc","All pion incident;Slice ID;Events");
  hsignal->SetNameTitle("hsignal","Pion interaction signal;Slice ID;Events");
  
  double Ndata = hsig3D->Integral();
  double Nmc = hmeas_MC->Integral();
  double Eweight = Ndata/Nmc;

  TH3D *htruth_3D = (TH3D*)response_SliceID_3D->Htruth(); // to normalize R matrix
  TH2D *hresponse_3D = (TH2D*)response_SliceID_3D->Hresponse();
  TMatrixD mR_3D(pow(pi::reco_nbins, 3), pow(pi::reco_nbins, 3));
  for (int i=0; i<pi::reco_nbins; i++) {
    for (int j=0; j<pi::reco_nbins; j++) {
      for (int k=0; k<pi::reco_nbins; k++) {
        double factor = 1/htruth_3D->GetBinContent(k+1, j+1, i+1);
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
    myfile_sys.open ("/dune/app/users/yinrui/Wiener-SVD-Unfolding/toys/syscov_MCXS_1127.txt", ios::app);
    /*//Ninc
     for (int i=0; i<pi::reco_nbins; ++i) {
     myfile_sys<<vmeas_inc(i) - hsiginc->GetBinContent(i+1)<<"\t";
     }
     myfile_sys<<endl;
     //Nini
     for (int i=0; i<pi::reco_nbins; ++i)
     myfile_sys<<hmeas_MC->ProjectionX()->GetBinContent(i+1)*Eweight - hsig3D->ProjectionX()->GetBinContent(i+1)<<"\t";
     myfile_sys<<endl;
     //Nint
     for (int i=0; i<pi::reco_nbins; ++i)
     myfile_sys<<hmeas_MC->ProjectionZ()->GetBinContent(i+1)*Eweight - hsig3D->ProjectionZ()->GetBinContent(i+1)<<"\t";
     myfile_sys<<endl;*/
    //3D
    for (int i=0; i<pi::reco_nbins; ++i)
      for (int j=0; j<pi::reco_nbins; ++j)
        for (int k=0; k<pi::reco_nbins; ++k) {
          int idx = i*pow(pi::reco_nbins,2) + j*pi::reco_nbins + k;
          myfile_sys<<vmeas_3D(idx) - hsig3D->GetBinContent(k+1, j+1, i+1)<<"\t";
        }
    myfile_sys<<endl<<endl;
    myfile_sys.close();
    
    /*fout->Write();
     fout->Close();
     return 0;*/
  }
  
  TH1D *hval_siginc_reco = (TH1D*)fmc->Get("h_recosliceid_pion_cuts");
  TH1D *hval_signal_reco = (TH1D*)fmc->Get("h_recosliceid_pioninelastic_cuts");
  TH1D *hval_sigini_reco = (TH1D*)fmc->Get("h_recoinisliceid_pion_cuts");
  TH1D *hval_trueinc = (TH1D*)fmc->Get("h_truesliceid_pion_all");
  TH1D *hval_trueint = (TH1D*)fmc->Get("h_truesliceid_pioninelastic_all");
  TH1D *hval_trueini = (TH1D*)fmc->Get("h_trueinisliceid_pion_all");
  double Nfakedata = hval_siginc_reco->Integral();
  hval_siginc_reco->Scale(Ndata/Nfakedata);
  hval_signal_reco->Scale(Ndata/Nfakedata);
  hval_sigini_reco->Scale(Ndata/Nfakedata);
  hval_trueinc->Scale(Ndata/Nfakedata);
  hval_trueint->Scale(Ndata/Nfakedata);
  hval_trueini->Scale(Ndata/Nfakedata);

  /*TMatrixD mcov_3D_stat(pow(pi::reco_nbins,3), pow(pi::reco_nbins,3));
  TMatrixD mcov_3D(pow(pi::reco_nbins,3), pow(pi::reco_nbins,3));
  mcov_3D_stat = unfold_3D.GetMeasuredCov();
  TVectorD Emeas_MC = response_SliceID_3D->Emeasured();
  cout<<"Ndata_sig: "<<Ndata<<"; Nmc_sig: "<<Nmc<<"; Nfakedata_sig: "<<Nfakedata<<endl;
  for(int i=0; i<pow(pi::reco_nbins,3); i++) {
    mcov_3D(i, i) = mcov_3D_stat(i, i) + pow(Emeas_MC(i)*Eweight, 2);
    //if (mcov_3D(i, i)!=0)
    //  cout<<mcov_3D_stat(i, i)<<"\t"<<mcov_3D(i, i)<<endl;
  }*/
  if (false) { // has_cov_input
    FILE *fcov_3D=fopen("/dune/app/users/yinrui/Wiener-SVD-Unfolding/toys/cov_MCXS_1127.txt","r");
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
  }
  
  TH3D *hsig3D_uf;
  TH1D *hsiginc_uf;
  TH1D *hsignal_uf;
  TH1D *hsigini_uf;
  
  hsig3D_uf = (TH3D*)unfold_3D.Hreco();
  hsig3D_uf->SetNameTitle("hsig3D_uf", "Unfolded 3D signal;Slice ID;Events");
  std::filebuf fb;
  fb.open(root["UnfoldTable"].asString().c_str(),std::ios::out);
  std::ostream os(&fb);
  unfold_3D.PrintTable(os);
  fb.close();
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
  
  if (false) { // toy
    ofstream myfile_sys_unfold;
    myfile_sys_unfold.open ("/dune/app/users/yinrui/Wiener-SVD-Unfolding/toys/syscov_MCXS_RooUnfold_1127.txt", ios::app);
    /*//Ninc
     for (int i=0; i<pi::reco_nbins; ++i)
     myfile_sys_unfold<<hsiginc_uf->GetBinContent(i+1)<<"\t";
     myfile_sys_unfold<<endl;
     //Nini
     for (int i=0; i<pi::reco_nbins; ++i)
     myfile_sys_unfold<<hsigini_uf->GetBinContent(i+1)<<"\t";
     myfile_sys_unfold<<endl;
     //Nint
     for (int i=0; i<pi::reco_nbins; ++i)
     myfile_sys_unfold<<hsignal_uf->GetBinContent(i+1)<<"\t";
     myfile_sys_unfold<<endl;*/
    //3D
    for (int i=0; i<pi::reco_nbins; ++i)
      for (int j=0; j<pi::reco_nbins; ++j)
        for (int k=0; k<pi::reco_nbins; ++k) {
          myfile_sys_unfold<<hsig3D_uf->GetBinContent(k+1, j+1, i+1)<<"\t";
        }
    myfile_sys_unfold<<endl<<endl;
    myfile_sys_unfold.close();
    
    fout->Write();
    fout->Close();
    return 0;
  }
  
  TMatrixD cov_matrix_3D = unfold_3D.Ereco();
  if (false) { // has_cov_input after unfolding
    FILE *fcov_3D=fopen("/dune/app/users/yinrui/Wiener-SVD-Unfolding/toys/cov_MCXS_RooUnfold_1127.txt","r");
    if (!fcov_3D) {
      cout<<"cov_3D_RooUnfold not found!"<<endl;
      return 1;
    }
    double vv;
    for(int i=0; i<pow(pi::reco_nbins, 3); i++) {
      for(int j=0; j<pow(pi::reco_nbins, 3); j++) {
        fscanf(fcov_3D, "%lf", &vv);
        cov_matrix_3D(i, j) = vv;
      }
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
  
  // get Ninc, Nint, Nini
  double Ninc[pi::true_nbins-1] = {0};
  double Nint[pi::true_nbins-1] = {0};
  double Nini[pi::true_nbins-1] = {0};
  double Nina[pi::true_nbins-1] = {0};
  double err_inc[pi::true_nbins-1] = {0};
  double err_int[pi::true_nbins-1] = {0};
  double err_ini[pi::true_nbins-1] = {0};
  double err_ina[pi::true_nbins-1] = {0};
  double SliceID[pi::true_nbins-1] = {0};

  for (int i = 0; i<pi::true_nbins-1; ++i){
    SliceID[i] = i+1;
    Nint[i] = hsignal_uf->GetBinContent(i+2);
    err_int[i] = hsignal_uf->GetBinError(i+2);
    Nini[i] = hsigini_uf->GetBinContent(i+2);
    err_ini[i] = hsigini_uf->GetBinError(i+2);
    Nina[i] = hsiginc_uf->GetBinContent(i+2);
    err_ina[i] = hsiginc_uf->GetBinError(i+2);
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
  // calculate covariance matrix for Ninc and Nint
  TMatrixD mNin(2*(pi::true_nbins-1), 3*pi::true_nbins);
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
  for (int bi=pi::true_nbins-1; bi<2*(pi::true_nbins-1); ++bi) { // Nint
    mNin(bi, pi::true_nbins+2+bi) = 1;
  }
  TMatrixD VNin_tmp(2*(pi::true_nbins-1), 3*pi::true_nbins);
  TMatrixD VNin(2*(pi::true_nbins-1), 2*(pi::true_nbins-1));
  VNin_tmp.Mult(mNin, Vhist);
  VNin.MultT(VNin_tmp, mNin); // 2(N-1)*2(N-1) covariance matrix
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
  TGraphErrors *gr_ina = new TGraphErrors(pi::true_nbins-1, SliceID, Nina, 0, err_ina);
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
  double xs[pi::true_nbins-1] = {0};
  double err_xs[pi::true_nbins-1] = {0};
  double KE[pi::true_nbins-1] = {0};
  double err_KE[pi::true_nbins-1] = {0};
  double dEdx[pi::true_nbins-1] = {0};
  TMatrixD mXS(pi::true_nbins-1, 2*(pi::true_nbins-1)); // Jacobian for (Ninc, Nint) to XS
  for (int i = 0; i<pi::true_nbins-1; ++i){
    KE[i] = (pi::true_KE[i]+pi::true_KE[i+1])/2;
    err_KE[i] = (pi::true_KE[i]-pi::true_KE[i+1])/2;
    dEdx[i] = bb.meandEdx(KE[i]); // MeV/cm
    xs[i] = dEdx[i]*MAr/(Density*NA*2*err_KE[i])*log(Ninc[i]/(Ninc[i]-Nint[i]))*1e27;
    //err_xs[i] = dEdx[i]*MAr/(Density*NA*2*err_KE[i])*sqrt(pow(Nint[i]*err_inc[i]/Ninc[i]/(Ninc[i]-Nint[i]),2)+pow(err_int[i]/(Ninc[i]-Nint[i]),2))*1e27;
    mXS(i,i) = dEdx[i]*MAr/(Density*NA*2*err_KE[i])*1e27 * (-Nint[i]/Ninc[i])/(Ninc[i]-Nint[i]);
    mXS(i,i+pi::true_nbins-1) = dEdx[i]*MAr/(Density*NA*2*err_KE[i])*1e27 * 1/(Ninc[i]-Nint[i]);
  }
  TMatrixD VXS_tmp(pi::true_nbins-1, 2*(pi::true_nbins-1));
  TMatrixD VXS(pi::true_nbins-1, pi::true_nbins-1);
  VXS_tmp.Mult(mXS, VNin);
  VXS.MultT(VXS_tmp, mXS); // (N-1)*(N-1) covariance matrix
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
  
  // test sample validation
  hval_siginc_reco->Write("hval_siginc_reco");
  hval_signal_reco->Write("hval_signal_reco");
  hval_sigini_reco->Write("hval_sigini_reco");
  hval_trueinc->Write("hval_trueinc");
  hval_trueint->Write("hval_trueint");
  hval_trueini->Write("hval_trueini");
  
  double Ninc_t[pi::true_nbins-1] = {0};
  double Nint_t[pi::true_nbins-1] = {0};
  double Nini_t[pi::true_nbins-1] = {0};
  double Nina_t[pi::true_nbins-1] = {0};
  double err_inc_t[pi::true_nbins-1] = {0};
  double err_int_t[pi::true_nbins-1] = {0};
  double err_ini_t[pi::true_nbins-1] = {0};
  double err_ina_t[pi::true_nbins-1] = {0};
  for (int i = 0; i<pi::true_nbins-1; ++i){
    Nint_t[i] = hval_trueint->GetBinContent(i+2);
    err_int_t[i] = hval_trueint->GetBinError(i+2);
    Nini_t[i] = hval_trueini->GetBinContent(i+2);
    err_ini_t[i] = hval_trueini->GetBinError(i+2);
    Nina_t[i] = hval_trueinc->GetBinContent(i+2);
    err_ina_t[i] = hval_trueinc->GetBinError(i+2);
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
      err_inc_t[i] += pow(hval_trueinc->GetBinError(j+2),2);
    }
    for (int j = i+1; j<pi::true_nbins-1; ++j){
      Ninc_t[i] -= hval_trueini->GetBinContent(j+2);
      err_inc_t[i] += pow(hval_trueini->GetBinError(j+2),2);
    }*/
    //err_inc_t[i] = sqrt(err_inc_t[i]);
    //err_inc_t[i] = sqrt(Ninc_t[i]);
  }
  TMatrixD Vhist_t(3*pi::true_nbins, 3*pi::true_nbins); // Vhist_t is diagonal
  for (int bi=0; bi<pi::true_nbins; ++bi) {
    int tbi = bi;
    Vhist_t(tbi,tbi) = pow(hval_trueini->GetBinError(bi+1),2);
    tbi += pi::true_nbins;
    Vhist_t(tbi,tbi) = pow(hval_trueinc->GetBinError(bi+1),2);
    tbi += pi::true_nbins;
    Vhist_t(tbi,tbi) = pow(hval_trueint->GetBinError(bi+1),2);
  }
  TMatrixD VNin_tmp_t(2*(pi::true_nbins-1), 3*pi::true_nbins);
  TMatrixD VNin_t(2*(pi::true_nbins-1), 2*(pi::true_nbins-1));
  VNin_tmp_t.Mult(mNin, Vhist_t);
  VNin_t.MultT(VNin_tmp_t, mNin); // 2(N-1)*2(N-1) covariance matrix
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
  TGraphErrors *gr_ina_t = new TGraphErrors(pi::true_nbins-1, SliceID, Nina_t, 0, err_ina_t);
  gr_ina_t->SetNameTitle("gr_ina_t", "Interaction_allpion number;Slice ID;Events");
  gr_ina_t->Write();
  double xs_t[pi::true_nbins-1] = {0};
  double err_xs_t[pi::true_nbins-1] = {0};
  TMatrixD mXS_t(pi::true_nbins-1, 2*(pi::true_nbins-1)); // Jacobian for (Ninc, Nint) to XS
  for (int i = 0; i<pi::true_nbins-1; ++i){
    xs_t[i] = dEdx[i]*MAr/(Density*NA*2*err_KE[i])*log(Ninc_t[i]/(Ninc_t[i]-Nint_t[i]))*1e27;
    //err_xs_t[i] = dEdx[i]*MAr/(Density*NA*2*err_KE[i])*sqrt(pow(Nint_t[i]*err_inc_t[i]/Ninc_t[i]/(Ninc_t[i]-Nint_t[i]),2)+pow(err_int_t[i]/(Ninc_t[i]-Nint_t[i]),2))*1e27;
    mXS_t(i,i) = dEdx[i]*MAr/(Density*NA*2*err_KE[i])*1e27 * (-Nint_t[i]/Ninc_t[i])/(Ninc_t[i]-Nint_t[i]);
    mXS_t(i,i+pi::true_nbins-1) = dEdx[i]*MAr/(Density*NA*2*err_KE[i])*1e27 * 1/(Ninc_t[i]-Nint_t[i]);
  }
  TMatrixD VXS_tmp_t(pi::true_nbins-1, 2*(pi::true_nbins-1));
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
  TGraphErrors *gr_truexs_allMC = (TGraphErrors*)fmc->Get("gr_truexs");
  gr_truexs_allMC->Write("gr_truexs_allMC");
  
  fout->Write();
  fout->Close();
  return 0;
}
