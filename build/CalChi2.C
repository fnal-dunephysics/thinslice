#include "../SliceParams.h"
#include <iostream>

void CalChi2(){
  
  TFile *file_N = TFile::Open("XSMC.root");
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
  TH1D *hxs_n = (TH1D*)file_N->Get("h1seltruthMC_1D");
  int nbins = hxs_N->GetNbinsX();
  int rank = 0;
  vector<int> idx;
  for (int i=1; i<=nbins; ++i) {
    if (hxs_n->GetBinContent(i) > 0) {
      idx.push_back(i);
      ++rank;
    }
    else cout<<"$ "<<i<<"\t"<<hxs_N->GetBinContent(i)<<"\t"<<hxs_n->GetBinContent(i)<<endl;
  }
  cout<<"$$ nbins "<<nbins<<"; rank "<<rank<<endl;
  
  int first=0;
  //rank=5;
  TVectorD xs_N(rank);
  TMatrixD Mcovxs_N(rank, rank);
  TVectorD xs_n(rank);
  for (int i=0; i<rank; ++i) {
    int idxi = idx.at(i+first);
    //xs_N(i) = gxs_N->GetPointY(i);
    //xs_n(i) = gxs_n->GetPointY(i);
    xs_N(i) = hxs_N->GetBinContent(idxi);
    xs_n(i) = hxs_n->GetBinContent(idxi);
    for (int j=0; j<rank; ++j) {
      int idxj = idx.at(j+first);
      Mcovxs_N(i, j) = hcovxs_N->GetBinContent(idxi, idxj);
    }
    //cout<<xs_N(i)<<"\t"<<xs_n(i)<<"\t\t"<<sqrt(Mcovxs_N(i, i))<<endl;
  }
  TMatrixD Minvxs_N(Mcovxs_N);
  Minvxs_N.Invert();
  cout<<"Determinant: "<<Mcovxs_N.Determinant()<<"\t"<<Minvxs_N.Determinant()<<endl;
  TVectorD diffxs(rank);
  diffxs = xs_N - xs_n;
  double chi2 = (Minvxs_N * diffxs) * diffxs;
  //cout<<"Chi2/nbins = "<<chi2<<"/"<<nbins<<" = "<<chi2/nbins<<endl;
  cout<<"Chi2 = "<<chi2<<endl;
  
  TFile *fout = TFile::Open("testCalChi2.root", "recreate");
  TH2D *tcov = new TH2D(Mcovxs_N);
  tcov->Write("tcov");
  TH2D *tinv = new TH2D(Minvxs_N);
  tinv->Write("tinv");
  TH1D *txsN = new TH1D(xs_N);
  txsN->Write("txsN");
  TH1D *txsn = new TH1D(xs_n);
  txsn->Write("txsn");
  TH1D *tdxs = new TH1D(diffxs);
  tdxs->Write("tdxs");
  fout->Write();
  fout->Close();
  
  ofstream txtcov;
  txtcov.open ("txtcov.txt", ios::trunc);
  for (int i=0; i<rank; ++i){
    for (int j=0; j<rank; ++j){
      txtcov<<Mcovxs_N(i,j)<<"\t";
    }
    txtcov<<endl;
  }
  txtcov.close();
  
  for (int i=0; i<rank; ++i) cout<<diffxs(i)<<",";
  cout<<endl;
  
  return 0;
}
