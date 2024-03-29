
const int typenum = 9;
const int boundbin = 22; // frontier of first TPC bin
const int uppbound = 50; // frontier of all interested TPC bin
const double xbins[boundbin+2] = {0.,10.,20.,30.,40.,50.,60.,70.,80.,90.,100.,110.,120.,130.,140.,150.,160.,170.,180.,190.,200.,210.,220.,10.*uppbound};

double calchi2(double muweight, TH1D *data_dist, TH1D *mc_dist_sep[typenum], TH1D *hdata, TH1D *hmc, TH1D *hmc0, bool save_plot=false){
  const double typeweight[typenum] = {1,1, muweight, 1,1,1,1,1,1};
  //double mc_inte = mc_dist->Integral(0, uppbound);
  double data_inte = data_dist->Integral(0, uppbound);
  double mc_inte_sep;
  double mc_inte_sep0;
  for (int j=0; j<typenum; ++j){
    mc_inte_sep += mc_dist_sep[j]->Integral(0, uppbound)*typeweight[j];
    mc_inte_sep0 += mc_dist_sep[j]->Integral(0, uppbound);
  }
  //cout<<mc_inte<<"\t"<<mc_inte_sep<<"\t"<<data_inte<<endl;
  double mcnorm = data_inte/mc_inte_sep;
  double mcnorm0 = data_inte/mc_inte_sep0;
  double fom = 0;
  double chi = 0;
  double mcbin;
  double mcbin0;
  for (int i=1; i<=boundbin; ++i){
    //mcbin = mc_dist->GetBinContent(i+1) * mcnorm;
    mcbin = 0;
    mcbin0 = 0;
    for (int j=0; j<typenum; ++j){
      mcbin += mc_dist_sep[j]->GetBinContent(i+1)*typeweight[j];
      mcbin0 += mc_dist_sep[j]->GetBinContent(i+1);
    }
    mcbin *= mcnorm;
    mcbin0 *= mcnorm0;
    
    chi = (data_dist->GetBinContent(i+1) - mcbin)/sqrt(data_dist->GetBinContent(i+1) + mcbin);
    fom += chi*chi;
    //cout<<i<<endl;
    //cout<<mcbin<<endl;
    //cout<<data_dist->GetBinError(i)<<endl;
    if (save_plot) {
      hdata->SetBinContent(i, data_dist->GetBinContent(i+1)/10);
      hdata->SetBinError(i, data_dist->GetBinError(i+1)/10);
      hmc->SetBinContent(i, mcbin/10);
      hmc0->SetBinContent(i, mcbin0/10);
    }
  }
  double over_data = data_dist->Integral(23, uppbound);
  //mcbin = mc_dist->Integral(23, uppbound) * mcnorm;
  mcbin = 0;
  mcbin0 = 0;
  for (int j=0; j<typenum; ++j){
    mcbin += mc_dist_sep[j]->Integral(23, uppbound)*typeweight[j];
    mcbin0 += mc_dist_sep[j]->Integral(23, uppbound);
  }
  mcbin *= mcnorm;
  mcbin0 *= mcnorm0;
  
  chi = (over_data - mcbin)/sqrt(over_data + mcbin);
  fom += chi*chi;
  fom /= (boundbin+1);
  if (save_plot) {
    cout<<"\nWeight = "<<muweight<<"\t FOM = "<<fom<<endl;
    
    hdata->SetBinContent(boundbin+1, over_data/280);
    hdata->SetBinError(boundbin+1, sqrt(over_data)/280);
    hmc->SetBinContent(boundbin+1, mcbin/280);
    hmc0->SetBinContent(boundbin+1, mcbin0/280);
    
    TCanvas* c1 = new TCanvas("c1","c1");
    hmc->SetLineColor(kRed);
    hmc0->SetLineColor(kBlue);
    hdata->GetYaxis()->SetRangeUser(0, 1.1*hmc0->GetBinContent(hmc0->GetMaximumBin()));
    hdata->Draw();
    hdata->GetXaxis()->SetTitle("Reco_track_length [cm]");
    hdata->GetYaxis()->SetTitle("Event/cm");
    hmc->Draw("same");
    hmc0->Draw("same");
    hdata->SetTitle("Muon reweighting");
    TLegend *leg1 = new TLegend(0.5,0.6,0.8,0.9);
    leg1->SetFillStyle(0);
    leg1->AddEntry(hdata,"Data","ple");
    leg1->AddEntry(hmc0, "MC default","l");
    leg1->AddEntry(hmc,"MC reweighted","l");
    leg1->Draw();
    c1->Print("muon_reweight.png");
  }
  return fom;
}

void muon_reweight(){
  TFile *mcfile = TFile::Open("../../build/mcprod4a.root");
  TFile *datafile = TFile::Open("../../build/mcprod4a.root"); // data.root
  
  int cut = 4;
  TH1D *data_dist = (TH1D*)datafile->Get(Form("hreco_trklen_%d_0", cut));
  //TH1D *mc_dist = (TH1D*)mcfile->Get(Form("hreco_trklen_%d_0", cut));
  TH1D *mc_dist_sep[typenum];
  for (int j=0; j<typenum; ++j){
    mc_dist_sep[j] = (TH1D*)mcfile->Get(Form("hreco_trklen_%d_%d", cut, j+1));
  }

  TH1D *hdata = new TH1D("hdata", "hdata", boundbin+1, xbins);
  TH1D *hmc = new TH1D("hmc", "hmc", boundbin+1, xbins);
  TH1D *hmc0 = new TH1D("hmc0", "hmc0", boundbin+1, xbins);
  
  vector<double> fom_list;
  double selw = 1.;
  double minfom = 999.;
  for (double w = 0.9; w <= 1.1; w+=0.01) {
    double fom = calchi2(w, data_dist, mc_dist_sep, hdata, hmc, hmc0);
    fom_list.push_back(fom);
    //cout<<fom<<", "; // python list format for plotting
    cout<<"weight = "<<w<<"\t fom = "<<fom<<endl;
    if (fom<minfom) {
      minfom = fom;
      selw = w;
    }
  }
  calchi2(selw, data_dist, mc_dist_sep, hdata, hmc, hmc0, true);
}
