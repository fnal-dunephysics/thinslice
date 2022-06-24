#include "ProtonInel.h"
#include "Unfold.h"
#include "EventType.h"
#include "EventSelection.h"
#include "TChain.h"
#include <iostream>

int main(){


  TChain *mcchain = new TChain();

  //chain->Add("/data/tjyang/dune/pduneana_Prod4.1_5_11_21.root/pduneana/beamana");
  //chain->Add("/data/tjyang/dune/pduneana_Prod4_1GeV_5_8_21.root/pduneana/beamana");

  mcchain->Add("/dune/data/users/calcuttj/pduneana_Prod4a_1GeV_5_14_21.root/pduneana/beamana");


  TChain *datachain = new TChain();
  datachain->Add("/dune/data/users/calcuttj/pduneana_Prod4_1GeV_5387_5_12_21.root/pduneana/beamana");

  //Unfold uf(p::nthinslices+2, -1, p::nthinslices+1, p::nthinslices+2, -1, p::nthinslices+1);
  TH2D* hist_reco = new TH2D("hist_reco","hist_reco", pi::reco_nbins, pi::reco_bins, pi::reco_nbins, pi::reco_bins);
  TH2D* hist_true = new TH2D("hist_true","hist_true", pi::true_nbins, pi::true_bins, pi::true_nbins, pi::true_bins);
  Unfold uf(hist_reco, hist_true);

  anavar mcevt(mcchain);

  ProtonInel mcths;
  mcths.SetOutputFileName("mcprod4a_proton.root");
  mcths.Run(mcevt, uf);

  anavar dataevt(datachain);

  ProtonInel dataths;
  dataths.SetOutputFileName("data_proton.root");
  dataths.Run(dataevt, uf);

  return 0;

}
