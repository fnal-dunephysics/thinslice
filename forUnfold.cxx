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

int main(){
  TFile *mcfile = TFile::Open("mcprod4a.root");
  RooUnfoldResponse *response_SliceID = (RooUnfoldResponse*)mcfile->Get("response_SliceID_3D");
  TH3D *sig3D = (TH3D*)response_SliceID->Htruth();
  int nn = pi::reco_nbins;
  TH1D *sig = new TH1D("sig", "sig", pow(nn,3), 0, pow(nn,3)); // sig
  TH3D *measure3D = (TH3D*)response_SliceID->Hmeasured();
  TH1D *measure = new TH1D("measure", "measure", pow(nn,3), 0, pow(nn,3)); // measure
  
  TFile *xsfile = TFile::Open("XS.root");
  TH3D *Rmeasure3D = (TH3D*)xsfile->Get("hsig3D");
  TH1D *Rmeasure = new TH1D("Rmeasure", "Rmeasure", pow(nn,3), 0, pow(nn,3)); // Rmeasure
  TRandom3 *r3 = new TRandom3(0);
  int idx = 1;
  for (int i=1; i<=nn; ++i)
    for (int j=1; j<=nn; ++j)
      for (int k=1; k<=nn; ++k) {
        double binc = Rmeasure3D->GetBinContent(k, j, i);
        double bine = Rmeasure3D->GetBinError(k, j, i);
        //Rmeasure->SetBinContent(idx, binc);
        //Rmeasure->SetBinError(idx, bine);
        binc = sig3D->GetBinContent(k, j, i);
        bine = sig3D->GetBinError(k, j, i);
        sig->SetBinContent(idx, binc);
        sig->SetBinError(idx, bine);
        binc = measure3D->GetBinContent(k, j, i);
        bine = measure3D->GetBinError(k, j, i);
        measure->SetBinContent(idx, binc);
        measure->SetBinError(idx, bine);
        Rmeasure->SetBinContent(idx, r3->Gaus(binc, 2*bine));
        Rmeasure->SetBinError(idx, 2*bine);
        ++idx;
      }
  
  TH2D *hresponse = (TH2D*)response_SliceID->Hresponse();
  //TH2D *hresponse = new TH2D("hresponse", "hresponse", pow(nn,3), 0, pow(nn,3), pow(nn,3), 0, pow(nn,3));
  for (int i=1; i<=pow(nn,3); i++) {
    double factor = 1/sig->GetBinContent(i);
    for (int j=1; j<=pow(nn,3); j++) {
      double tmp = hresponse->GetBinContent(j, i);
      if (tmp != 0)
        hresponse->SetBinContent(j, i, tmp*factor);
    }
  } // hresponse
  
  TMatrixD mcov(pow(nn,3), pow(nn,3));
  TH2D *hcov = new TH2D(mcov);
  for(int i=1; i<=pow(nn,3); i++) {
    hcov->SetBinContent(i, i, pow(Rmeasure->GetBinError(i), 2));
  } // hcov */
  /*FILE *covfile=fopen("/dune/app/users/yinrui/Wiener-SVD-Unfolding/toys/cov_MCXS_3D_1110.txt","r");
  if (!covfile) {
    cout<<"covfile not found!"<<endl;
    return 1;
  }
  double vv;
  for(int i=1; i<=pow(nn,3); i++) {
    for(int j=1; j<=pow(nn,3); j++) {
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

