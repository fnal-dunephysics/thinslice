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

const double beam_startX_data = -28.3483;
const double beam_startY_data = 424.553;
const double beam_startZ_data = 3.19841;
const double beam_startX_rms_data = 4.63594;
const double beam_startY_rms_data = 5.21649;
const double beam_startZ_rms_data = 1.2887;

const double beam_startX_mc = -30.6692;
const double beam_startY_mc = 422.263;
const double beam_startZ_mc = 0.1106;
const double beam_startX_rms_mc = 5.172;
const double beam_startY_rms_mc = 4.61689;
const double beam_startZ_rms_mc = 0.212763;

const double beam_angleX_data = 100.464;
const double beam_angleY_data = 103.442;
const double beam_angleZ_data = 17.6633;

const double beam_angleX_mc = 101.547;
const double beam_angleY_mc = 101.247;
const double beam_angleZ_mc = 16.5864;

const double beam_startX_data_inst = -30.9033;
const double beam_startY_data_inst = 422.406;
const double beam_startX_rms_data_inst = 4.17987;
const double beam_startY_rms_data_inst = 3.73181;

const double beam_startX_mc_inst = -28.8615;
const double beam_startY_mc_inst = 421.662;
const double beam_startX_rms_mc_inst = 4.551;
const double beam_startY_rms_mc_inst = 3.90137;

// fiducial volume
const double fidvol_low = 0; //cm
const double fidvol_upp = 220; //cm

#endif
