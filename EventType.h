#ifndef EVENTTYPE_H
#define EVENTTYPE_H

#include "anavar.h"

namespace pi{

  const unsigned int nIntTypes = 9;

  const char intTypeName[nIntTypes+1][100] = {"Data",
                                              "PiInel",
                                              "PiDecay",
                                              "Muon",
                                              "misID:cosmic",
                                              "misID:p",
                                              "misID:pi",
                                              "misID:mu",
                                              "misID:e/#gamma",
                                              "misID:other"};
  
  enum intType{
    kData,
    kPiInel,
    kPiElas,
    kMuon,
    kMIDcosmic,
    kMIDp,
    kMIDpi,
    kMIDmu,
    kMIDeg,
    kMIDother
  };
}

namespace p{

  const unsigned int nIntTypes = 8;

  const char intTypeName[nIntTypes+1][100] = {"Data",
                                              "PInel",
                                              "PElas",
                                              "misID:cosmic",
                                              "misID:p",
                                              "misID:pi",
                                              "misID:mu",
                                              "misID:e/#gamma",
                                              "misID:other"};
  
  enum intType{
    kData,
    kPInel,
    kPElas,
    kMIDcosmic,
    kMIDp,
    kMIDpi,
    kMIDmu,
    kMIDeg,
    kMIDother
  };
}


#endif
