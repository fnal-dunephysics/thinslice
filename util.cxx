#include "util.h"
#include "TVector3.h"
#include "TMath.h"
#include "TStyle.h"
#include "TGaxis.h"
#include <iostream>

using namespace std;

double CalWeight(const anavar & evt, const int &partype, double rdm_radius, double rdm_angle){
  double weight = 1.;
  //return weight;
  
  double mufrac = 1;//.71;
  
  double mom_mu0 = 0.9940;
  double mom_sigma0 = 0.0545;
  double mom_mu = 0.98939;
  double mom_sigma = 0.06394;
  if (rdm_radius != 0) {
    double oval_cx = 1.01853;
    double oval_cy = 0.07245;
    double oval_a = 0.00557;
    double oval_b = 0.00309;
    double oval_phi = 0.776;
    mom_mu = oval_cx + rdm_radius * (oval_a*TMath::Cos(oval_phi)*TMath::Cos(rdm_angle) - oval_b*TMath::Sin(oval_phi)*TMath::Sin(rdm_angle));
    mom_sigma = oval_cy + rdm_radius * (oval_a*TMath::Sin(oval_phi)*TMath::Cos(rdm_angle) + oval_b*TMath::Cos(oval_phi)*TMath::Sin(rdm_angle));
  }
  double wlimit = 3;
  
  if (evt.MC) {
    // momentum reweight (outlier weights are set to 1e-5)
    double deno = exp(-pow((evt.true_beam_startP-mom_mu0)/mom_sigma0,2)/2);
    //if (deno < wlimit) deno = wlimit;
    double numo = exp(-pow((evt.true_beam_startP-mom_mu)/mom_sigma,2)/2);
    //if (numo < wlimit) numo = wlimit;
    weight *= numo;
    weight /= deno;
    if (weight>wlimit) weight=wlimit;
    if (weight<1./wlimit) weight=1./wlimit;
    if (evt.true_beam_PDG == -13) weight = 1;
    // muon reweight
    if (evt.true_beam_PDG == -13 && evt.reco_beam_true_byE_matched) { // kMuon
      weight *= mufrac;
    }
  }

  return weight;
}

double CalBkgW(const anavar & evt, double mu_weight, double p_weight, double spi_weight){
  double weight = 1.;

  if (evt.MC) {
    // if (fakedata) { // it seems true MC does not need bkgw
    if (evt.reco_beam_true_byE_matched) {
      if (evt.true_beam_PDG == -13) { // muon bkg
        weight = mu_weight;
      }
    }
    else {
      if (evt.reco_beam_true_byE_PDG == -13) { // secondary muon bkg
        weight = mu_weight;
      }
      else if (evt.reco_beam_true_byE_PDG == 2212) { // secondary proton bkg
        weight = p_weight;
      }
      else if (evt.reco_beam_true_byE_PDG == 211) { // secondary pion bkg
        weight = spi_weight;
      }
    }
  }
  return weight;
}
double CalG4RW_ori(const anavar & evt, const double w1, const double w2){
  vector<double> tot_inel_0_600 = (*evt.g4rw_full_grid_piplus_coeffs)[7];
  vector<double> tot_inel_600_2000 = (*evt.g4rw_full_grid_piplus_coeffs)[8];
  double w_tot_inel_0_600 = w1; // KE 0 ~ 476.44931 MeV
  double w_tot_inel_600_2000 = w2; // KE 476.44931 ~ 1865.2940 MeV
  
  double g4rw = 1;
  if (evt.true_beam_PDG == 211 && (tot_inel_0_600.size() > 0) && (tot_inel_600_2000.size() > 0) && (w_tot_inel_0_600 != 1 || w_tot_inel_600_2000 != 1)){
    double g4rw_tot_inel_0_600 = 0;
    for (size_t i = 0; i < tot_inel_0_600.size(); ++i) {
      g4rw_tot_inel_0_600 += tot_inel_0_600[i] * pow(w_tot_inel_0_600, i);
    }
    double g4rw_tot_inel_600_2000 = 0;
    for (size_t i = 0; i < tot_inel_600_2000.size(); ++i) {
      g4rw_tot_inel_600_2000 += tot_inel_600_2000[i] * pow(w_tot_inel_600_2000, i);
    }
    g4rw *= g4rw_tot_inel_0_600;
    g4rw *= g4rw_tot_inel_600_2000;
  }
  return g4rw;
}
double CalG4RW(const anavar & evt, double weight[]){
  vector<vector<double>> tot_inel;

  double g4rw = 1;
  if (evt.true_beam_PDG == 211) {
    for (int i=0; i<20; ++i) {
      double g4rw_tot_inel = 0;
      tot_inel.push_back((*evt.g4rw_full_grid_piplus_coeffs)[i]);
      if (tot_inel[i].size() > 0) {
        double sum = 0;
        for (size_t j = 0; j < tot_inel[i].size(); ++j) {
          sum += tot_inel[i][j];
          g4rw_tot_inel += tot_inel[i][j] * pow(weight[i], j);
        }

        if (abs(sum-1)>0.1) {
          cout<<"$$$$$"<<i<<"\t"<<sum<<endl;
          cout<<(*evt.g4rw_full_grid_piplus_weights)[i][9]<<endl;
          cout<<evt.run<<"\t"<<evt.subrun<<"\t"<<evt.event<<endl;
          g4rw_tot_inel = 1;
        }
        g4rw *= g4rw_tot_inel;
      }
    }
  }
  
  return g4rw;
}

void FillHistVec1D(TH1D *hist[pi::nIntTypes+1], const double &value, const int &partype, double weight, bool fill_underflow, bool fill_overflow){
  //hist[0]->Fill(value);
  if (partype>=0 && partype < pi::nIntTypes+1){
    if (value < hist[partype]->GetXaxis()->GetXmin()){ // underflow values
      if (fill_underflow) hist[partype]->Fill(hist[partype]->GetXaxis()->GetXmin(), weight);
    }
    else if (value < hist[partype]->GetXaxis()->GetXmax()){
      hist[partype]->Fill(value, weight);
    }
    else{ // overflow values
      if (fill_overflow) hist[partype]->Fill(hist[partype]->GetXaxis()->GetXmax()-0.000001, weight);
    }
  }
}

void FillHistVec2D(TH2D *hist[pi::nIntTypes+1], const double &value1, const double &value2, const int &partype, double weight){
  //hist[0]->Fill(value1, value2);
  if (partype>=0 && partype < pi::nIntTypes+1){
    hist[partype]->Fill(value1, value2, weight);
  }
}
void FillHistVec3D(TH3D *hist[pi::nIntTypes+1], const double &value1, const double &value2, const double &value3, const int &partype, double weight){
  //hist[0]->Fill(value1, value2);
  if (partype>=0 && partype < pi::nIntTypes+1){
    hist[partype]->Fill(value1, value2, value3, weight);
  }
}
void FillProfVec(TProfile *profile[pi::nIntTypes+1], const double &value1, const double &value2, const int &partype){
  //hist[0]->Fill(value1, value2);
  if (partype>=0 && partype < pi::nIntTypes+1){
    profile[partype]->Fill(value1, value2);
  }
}

void FillHist1D(TH1D *hist, const double &value, const double &wei){
  if (value<hist->GetXaxis()->GetXmin()){
    hist->Fill(hist->GetXaxis()->GetXmin(), wei);
  }
  else if (value<hist->GetXaxis()->GetXmax()){
    hist->Fill(value, wei);
  }
  else{
    hist->Fill(hist->GetXaxis()->GetXmax()-0.000001, wei);
  }
}

double GetPionKE(double length){
  if (length>0){
    // https://arxiv.org/abs/1306.1712
    return (8/(1-0.37))*pow(length,1-0.37);
  }
  else{
    cout<<"Invalid length "<<length<<", returning -1"<<endl;
    return -1;
  }
}

double GetTheta(double x, double y, double z){
  TVector3 v(x,y,z);
  if (!v.Mag()){
    cout<<"Empty vector, returning -1"<<endl;
    return -1;
  }
  else{
    return v.Theta()*180/TMath::Pi();
  }
}

double GetPhi(double x, double y, double z){
  TVector3 v(x,y,z);
  if (!v.Mag()){
    cout<<"Empty vector, returning -1"<<endl;
    return -1;
  }
  else{
    return v.Phi()*180/TMath::Pi();
  }
}

double GetThetaxz(double x, double y, double z){
  TVector3 v(x,y,z);
  if (!v.Mag()){
    cout<<"Empty vector, returning -1"<<endl;
    return -1;
  }
  else{
    return atan2(v.X(), v.Z())*180/TMath::Pi();
  }
}

double GetThetayz(double x, double y, double z){
  TVector3 v(x,y,z);
  if (!v.Mag()){
    cout<<"Empty vector, returning -1"<<endl;
    return -1;
  }
  else{
    return atan2(v.Y(), v.Z())*180/TMath::Pi();
  }
}

void SetProtoDUNEStyle(){

  // This is the protoDUNE style file

  TStyle* protoDUNEStyle = new  TStyle("protoDUNEStyle", "ProtoDUNE Style");

  // Colors

  //set the background color to white
  protoDUNEStyle->SetFillColor(10);
  protoDUNEStyle->SetFrameFillColor(10);
  protoDUNEStyle->SetCanvasColor(10);
  protoDUNEStyle->SetPadColor(10);
  protoDUNEStyle->SetTitleFillColor(10);
  protoDUNEStyle->SetStatColor(10);

  //dont put a colored frame around the plots
  protoDUNEStyle->SetFrameBorderMode(0);
  protoDUNEStyle->SetCanvasBorderMode(0);
  protoDUNEStyle->SetPadBorderMode(0);

  //use the primary color palette
  protoDUNEStyle->SetPalette(1,0);

  //set the default line color for a histogram to be black
  protoDUNEStyle->SetHistLineColor(kBlack);

  //set the default line color for a fit function to be red
  protoDUNEStyle->SetFuncColor(kRed);

  //make the axis labels black
  protoDUNEStyle->SetLabelColor(kBlack,"xyz");

  //set the default title color to be black
  protoDUNEStyle->SetTitleColor(kBlack);

  // Sizes

  //set the margins
  protoDUNEStyle->SetPadBottomMargin(0.125);
  protoDUNEStyle->SetPadTopMargin(0.075);
  protoDUNEStyle->SetPadLeftMargin(0.12);
  protoDUNEStyle->SetPadRightMargin(0.1);

  //set axis label and title text sizes
  protoDUNEStyle->SetLabelSize(0.045,"xy");
  protoDUNEStyle->SetLabelSize(0.035,"z");
  protoDUNEStyle->SetLabelOffset(0.005,"xy");
  //protoDUNEStyle->SetLabelOffset(0.005,"z");
  protoDUNEStyle->SetTitleSize(0.05,"xyz");
  protoDUNEStyle->SetTitleOffset(1.15,"x");
  protoDUNEStyle->SetTitleOffset(1.15,"yz");
  protoDUNEStyle->SetStatFontSize(0.05);
  protoDUNEStyle->SetTextSize(0.05);
  protoDUNEStyle->SetTitleBorderSize(0);
  protoDUNEStyle->SetStatBorderSize(0);

  //set line widths
  protoDUNEStyle->SetHistLineWidth(3);
  protoDUNEStyle->SetFrameLineWidth(2);
  protoDUNEStyle->SetFuncWidth(2);

  // Misc

  //align the titles to be centered
  //protoDUNEStyle->SetTitleAlign(22);

  //set the number of divisions to show
  protoDUNEStyle->SetNdivisions(506, "xy");

  //turn off xy grids
  protoDUNEStyle->SetPadGridX(0);
  protoDUNEStyle->SetPadGridY(0);

  //set the tick mark style
  protoDUNEStyle->SetPadTickX(1);
  protoDUNEStyle->SetPadTickY(1);

  //show the fit parameters in a box
  protoDUNEStyle->SetOptFit(1111);

  //turn off all other stats
  //protoDUNEStyle->SetOptStat(0000000);
  protoDUNEStyle->SetStatBorderSize(1);
  protoDUNEStyle->SetStatFont      (62);
  protoDUNEStyle->SetOptStat       (111111);
  protoDUNEStyle->SetStatColor     (0);
  protoDUNEStyle->SetStatX         (0.93);
  protoDUNEStyle->SetStatY         (0.90);
  protoDUNEStyle->SetLegendBorderSize(0);

  //marker settings
  protoDUNEStyle->SetMarkerStyle(20);
  protoDUNEStyle->SetMarkerSize(0.9);

  // Fonts
  const int kProtoDUNEFont = 62;

  protoDUNEStyle->SetStatFont(kProtoDUNEFont);
  protoDUNEStyle->SetLabelFont(kProtoDUNEFont,"xyz");
  protoDUNEStyle->SetTitleFont(kProtoDUNEFont,"xyz");
  protoDUNEStyle->SetTitleFont(kProtoDUNEFont,"");
  protoDUNEStyle->SetTextFont(kProtoDUNEFont);

  // ---------------------------------------------------------
  // Additions from George 26/06/2019
  protoDUNEStyle->SetCanvasBorderSize(0);
  protoDUNEStyle->SetFrameBorderSize(0);
  protoDUNEStyle->SetDrawBorder(0);
  protoDUNEStyle->SetTitleBorderSize(0);

  // Set the size (in pixels) of the small lines drawn at the end of the error bars
  protoDUNEStyle->SetEndErrorSize(4);

  // Set option to strip decimals when drawing axis labels.
  protoDUNEStyle->SetStripDecimals(kFALSE);

  protoDUNEStyle->SetLegendBorderSize(0);
  protoDUNEStyle->SetLegendFont(kProtoDUNEFont);

  protoDUNEStyle->SetLabelOffset(0.015, "x");
  protoDUNEStyle->SetLabelOffset(0.015, "y");
  protoDUNEStyle->SetLabelOffset(0.01, "z");

  protoDUNEStyle->SetTitleStyle(0);
  protoDUNEStyle->SetTitleFont(kProtoDUNEFont, "pad");
  protoDUNEStyle->SetTitleX(0.1f);
  protoDUNEStyle->SetTitleY(.98f);
  protoDUNEStyle->SetTitleW(0.8f);
  protoDUNEStyle->SetLineStyleString(2, "[12 12]"); // postscript dashes

  // get rid of X error bars and y error bar caps
  protoDUNEStyle->SetErrorX(0.001);

  protoDUNEStyle->SetNumberContours(255);
  protoDUNEStyle->SetPalette(kBird);


  // ---------------------------------------------------------
  //done
  protoDUNEStyle->cd();

  gROOT->ForceStyle();
  gStyle->ls();

/*
  if (gROOT->GetVersionInt()>51200) {
    TColor::InitializeColors();
    const Int_t NRGBs = 5;
    const Int_t NCont = 255;
    Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
    Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
    Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
    Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
    TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
    gStyle->SetNumberContours(NCont);
  }
*/


  // Avoid too many decimal places in the axis labels
  //TGaxis::SetMaxDigits(4);
}

TFile f_pid("/dune/data/users/yinrui/HadronXSfiles/dEdxrestemplates.root");
TProfile* dedx_range_pro = (TProfile*)f_pid.Get("dedx_range_pro");
double chi2pid(std::vector<double> &trkdedx, std::vector<double> &trkres) { // from Heng-Ye
        int npt = 0;
        double chi2pro = 0;
        for (size_t i = 0; i<trkdedx.size(); ++i){ //hits
                //ignore the first and the last point
                if (i==0 || i==trkdedx.size()-1) continue;
                if (trkdedx[i]>1000) continue; //protect against large pulse height

                int bin = dedx_range_pro->FindBin(trkres[i]);
                //    std::cout<<"bin proton "<<bin<<std::endl;
                if (bin>=1&&bin<=dedx_range_pro->GetNbinsX()){
                        double bincpro = dedx_range_pro->GetBinContent(bin);
                        if (bincpro<1e-6){//for 0 bin content, using neighboring bins
                                bincpro = (dedx_range_pro->GetBinContent(bin-1)+dedx_range_pro->GetBinContent(bin+1))/2;
                        }
                        double binepro = dedx_range_pro->GetBinError(bin);
                        if (binepro<1e-6){
                                binepro = (dedx_range_pro->GetBinError(bin-1)+dedx_range_pro->GetBinError(bin+1))/2;
                        }
                        double errdedx = 0.04231+0.0001783*trkdedx[i]*trkdedx[i]; //resolution on dE/dx
                        errdedx *= trkdedx[i];
                        chi2pro += pow((trkdedx[i]-bincpro)/std::sqrt(pow(binepro,2)+pow(errdedx,2)),2);
                        ++npt;
                }
        }
        if (npt>0) return (chi2pro/npt);
        else return 9999;
}
