#include "ThinSlice.h"
#include "EventType.h"
#include "EventSelection.h"
#include "TChain.h"
#include <iostream>

int main(){

  int nevents[nParTypes] = {0};

  TChain *chain = new TChain();

  chain->Add("/data/tjyang/dune/pduneana_Prod4.1_1GeV_3_26_21.root/pduneana/beamana");

  HadAna ana(chain);
  ana.AddTruePDG(-13);
  ana.AddTruePDG(211);
//  ana.SetPandoraSlicePDG(13);
//  ana.SetOutputFileName("hadana.root");
//  ana.BookHistograms();

  ThinSlice ths;
  TFile f("hadana.root","recreate");
  ths.BookHistograms();

  Long64_t nentries = ana.fChain->GetEntriesFast();
  
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    if (jentry%10000==0) std::cout<<jentry<<"/"<<nentries<<std::endl;
    Long64_t ientry = ana.LoadTree(jentry);
    if (ientry < 0) break;
    nb = ana.fChain->GetEntry(jentry);   nbytes += nb;
    // if (Cut(ientry) < 0) continue;
    //std::cout<<ana.run<<" "<<ana.event<<" "<<ana.MC<<" "<<ana.reco_beam_true_byE_matched<<" "<<ana.true_beam_PDG<<" "<<(*ana.true_beam_endProcess)<<std::endl;
    //std::cout<<GetParType(ana)<<std::endl;
    if (!ana.isTrueSelectedPart()) continue;

    ths.ProcessEvent(ana);

//    ana.ProcessEvent();
//    int partype = ana.GetParType();
//    ++nevents[0];
//    if (partype<nParTypes+1){
//      ++nevents[partype];
//    }
//    ana.FillHistograms(kNocut);
//    if (!ana.PassPandoraSliceCut()) continue;
//    ana.FillHistograms(kPandoraSlice);
//    if (!ana.PassBeamQualityCut()) continue;
//    ana.FillHistograms(kBeamQuality);
//    if (!ana.PassAPA3Cut()) continue;
//    ana.FillHistograms(kAPA3);
//    if (!ana.PassCaloSizeCut()) continue;
//    ana.FillHistograms(kCaloSize);
//    if (!ana.PassMichelScoreCut()) continue;
//    ana.FillHistograms(kMichelScore);
//    if (!ana.PassMediandEdxCut()) continue;
//    ana.FillHistograms(kMediandEdx);
  }

//  for (int i = 0; i<nParTypes; ++i){
//    std::cout<<i<<" "<<nevents[i]<<std::endl;
//  }
//
//  ana.SaveHistograms();

  ths.CalcXS();

  f.Write();
  return 0;

}