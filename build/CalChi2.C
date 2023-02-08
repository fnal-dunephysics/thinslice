#include "../SliceParams.h"

void CalChi2(){
  
  TFile *file_N = TFile::Open("XS.root");
  /*TGraphErrors *gxs_N = (TGraphErrors*)file_N->Get("gr_recoxs");
  TH2D *hcovxs_N = (TH2D*)file_N->Get("VXS_3D");
  
  //TFile *file_n = TFile::Open("XS_2iter.root");
  //TGraphErrors *gxs_n = (TGraphErrors*)file_n->Get("gr_recoxs");
  TGraphErrors *gxs_n = (TGraphErrors*)file_N->Get("gr_truexs");
  
  double nbins = pi::reco_nbins - 1;*/
  
  /*TH2D *hcovxs_N = (TH2D*)file_N->Get("covinput_1D");
  TH1D *hxs_N = (TH1D*)file_N->Get("h1measdata_1D");
  TH1D *hxs_n = (TH1D*)file_N->Get("h1measMC_1D");*/
  TH2D *hcovxs_N = (TH2D*)file_N->Get("covariance_1D");
  TH1D *hxs_N = (TH1D*)file_N->Get("h1unfdata_1D");
  TH1D *hxs_n = (TH1D*)file_N->Get("h1truthMC_1D");
  int nbins = hxs_N->GetNbinsX();
  int rank = 0;
  vector<int> idx;
  for (int i=1; i<=nbins; ++i) {
    if (hcovxs_N->GetBinContent(i, i) != 0) {
      idx.push_back(i);
      ++rank;
    }
    else cout<<"$ "<<i<<endl;
  }
  cout<<"$$ nbins "<<nbins<<"; rank "<<rank<<endl;
  
  TVectorD xs_N(rank);
  TMatrixD Mcovxs_N(rank, rank);
  TVectorD xs_n(rank);
  for (int i=0; i<rank; ++i) {
    int idxi = idx.at(i);
    //xs_N(i) = gxs_N->GetPointY(i);
    //xs_n(i) = gxs_n->GetPointY(i);
    xs_N(i) = hxs_N->GetBinContent(idxi);
    xs_n(i) = hxs_n->GetBinContent(idxi);
    for (int j=0; j<rank; ++j) {
      int idxj = idx.at(j);
      Mcovxs_N(i, j) = hcovxs_N->GetBinContent(idxi, idxj);
    }
  }
  TMatrixD Minvxs_N(Mcovxs_N);
  Minvxs_N = Minvxs_N.Invert();
  TVectorD chixs(rank);
  chixs = xs_n - xs_N;
  double chi2 = (Minvxs_N * chixs) * chixs;
  //cout<<"Chi2/nbins = "<<chi2<<"/"<<nbins<<" = "<<chi2/nbins<<endl;
  cout<<"Chi2 = "<<chi2<<endl;
}
