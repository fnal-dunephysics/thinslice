#ifndef EVENTSELECTION_H
#define EVENTSELECTION_H

namespace pi{
  const unsigned int nCuts = 7;

  // can change order of cuts
  const char cutName[nCuts][100] = {
    "Nocut",
    "PandoraSlice",
    "CaloSize",
    "BeamQuality",
    "ProtonCut",
    "MichelScore",
    "APA3",
  };
  
  enum cut{
    kNocut = 0,
    kPandoraSlice,
    kCaloSize,
    kBeamQuality,
    kProtonCut,
    kMichelScore,
    kAPA3,
  };
}

namespace p{
  const unsigned int nCuts = 5;

  const char cutName[nCuts][100] = {"Nocut",
                                    "PandoraSlice",
                                    "CaloSize",
                                    "BeamQuality",
                                    "StoppingProtonCut"};

  enum cut{
    kNocut = 0,
    kPandoraSlice,
    kCaloSize,
    kBeamQuality,
    kStoppingProtonCut
  };
}

const double beam_startX_data = -28.722;
const double beam_startY_data = 424.216;
const double beam_startZ_data = 3.17107;
const double beam_startX_rms_data = 3.86992;
const double beam_startY_rms_data = 4.59281;
const double beam_startZ_rms_data = 1.22898;

const double beam_startX_mc = -30.734;
const double beam_startY_mc = 422.495;
const double beam_startZ_mc = 0.0494442;
const double beam_startX_rms_mc = 4.75219;
const double beam_startY_rms_mc = 4.23009;
const double beam_startZ_rms_mc = 0.206941;

const double beam_angleX_data = 100.834;
const double beam_angleY_data = 104.119;
const double beam_angleZ_data = 18.2308;

const double beam_angleX_mc = 101.667;
const double beam_angleY_mc = 101.14;
const double beam_angleZ_mc = 16.5398;

const double beam_startX_data_inst = -31.2908;
const double beam_startY_data_inst = 422.105;
const double beam_startX_rms_data_inst = 3.79071;
const double beam_startY_rms_data_inst = 3.43891;

const double beam_startX_mc_inst = -29.0375;
const double beam_startY_mc_inst = 421.82;
const double beam_startX_rms_mc_inst = 4.61462;
const double beam_startY_rms_mc_inst = 3.78259;

// fiducial volume
const double fidvol_low = 0; //cm
const double fidvol_upp = 220; //cm

#endif
