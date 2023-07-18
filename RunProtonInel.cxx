#include "ProtonInel.h"
#include "Unfold.h"
#include "EventType.h"
#include "EventSelection.h"
#include "TChain.h"
#include "json/json.h"
#include <fstream>
#include <iostream>
#include <string>

using namespace std;

int main(int argc, char ** argv){

  bool found_config = false;

  string config_file;

  for (int iArg = 1; iArg < argc; iArg++) {
    if (!strcasecmp(argv[iArg],"-c")) {
     config_file = argv[++iArg];
     found_config = true;
    }
    if (!strcasecmp(argv[iArg],"-h")) {
      std::cout << "Usage: RunCrossSection " <<
                   "-c config.json " << std::endl;
      return 1;
    }
  }

  if (!found_config){
    cout<<"Error: no configuration file was provided! Please provide with '-i'" <<endl;
    return 1;
  }

  Json::Value root;
  ifstream file(config_file);
  file >> root;
  cout<<root<<endl;
  //cout<<root["mcfile"].asString()<<endl;

  TChain *mcchain = new TChain();

  //chain->Add("/data/tjyang/dune/pduneana_Prod4.1_5_11_21.root/pduneana/beamana");
  //chain->Add("/data/tjyang/dune/pduneana_Prod4_1GeV_5_8_21.root/pduneana/beamana");

  //mcchain->Add("pduneana_mc.root/pduneana/beamana"); // test
  mcchain->Add(Form("%s/beamana", root["mcfile"].asString().c_str()));

  TChain *datachain = new TChain();
  //datachain->Add("pduneana.root/pduneana/beamana"); // test
  datachain->Add(Form("%s/beamana", root["datafile"].asString().c_str()));

  TH3D* hist_reco = new TH3D("hist_reco","hist_reco", pi::reco_nbins, pi::reco_bins, pi::reco_nbins, pi::reco_bins, pi::reco_nbins, pi::reco_bins);
  TH3D* hist_true = new TH3D("hist_true","hist_true", pi::true_nbins, pi::true_bins, pi::true_nbins, pi::true_bins, pi::true_nbins, pi::true_bins);
  Unfold uf(hist_reco, hist_true);

  anavar mcevt(mcchain);

  ProtonInel mcths;
  mcths.SetOutputFileName(root["mcoutfile"].asString());
  mcths.Run(mcevt, uf, -1, false, false);

  anavar dataevt(datachain);

  ProtonInel dataths;
  dataths.SetOutputFileName(root["dataoutfile"].asString());
  dataths.Run(dataevt, uf, -1, false, false);

  return 0;

}
