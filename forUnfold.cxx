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

int main(){
  TFile *mcfile = TFile::Open("mcprod4a.root");
  RooUnfoldResponse *response_SliceID = (RooUnfoldResponse*)mcfile->Get("response_SliceID_Inc");
  TH1D *sig = (TH1D*)response_SliceID->Htruth(); // sig
  TH1D *measure = (TH1D*)response_SliceID->Hmeasured(); // measure
  TH2D *hresponse = (TH2D*)response_SliceID->Hresponse();
  for (int i=1; i<=pi::reco_nbins; i++) {
    double factor = 1/sig->GetBinContent(i);
    for (int j=1; j<=pi::reco_nbins; j++) {
      double tmp = hresponse->GetBinContent(j, i);
      hresponse->SetBinContent(j, i, tmp*factor);
    }
  } // hresponse
  
  TFile *xsfile = TFile::Open("XS.root");
  TH1D *Rmeasure = (TH1D*)xsfile->Get("hsiginc"); // Rmeasure
  
  TMatrixD mcov(pi::reco_nbins, pi::reco_nbins);
  TH2D *hcov = new TH2D(mcov);
  /*for(int i=1; i<=pi::reco_nbins; i++) {
    hcov->SetBinContent(i, i, pow(Rmeasure->GetBinError(i), 2));
  } // hcov */
  FILE *covfile=fopen("/dune/app/users/yinrui/Wiener-SVD-Unfolding/toys/cov_MCXS_inc_1110.txt","r");
  if (!covfile) {
    cout<<"covfile not found!"<<endl;
    return 1;
  }
  double vv;
  for(int i=1; i<=pi::reco_nbins; i++) {
    for(int j=1; j<=pi::reco_nbins; j++) {
      fscanf(covfile, "%lf", &vv);
      hcov->SetBinContent(i, j, vv);
    }
  } // hcov */
  
  TFile *fout  = TFile::Open("forUnfold.root", "recreate");
  sig->Write("sig");
  Rmeasure->Write("Rmeasure");
  hresponse->Write("hresponse");
  hcov->Write("hcov");
  measure->Write("measure");
  fout->Write();
  fout->Close();
  return 0;
}

