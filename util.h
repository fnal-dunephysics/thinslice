#ifndef UTIL_H
#define UTIL_H

#include "TH1D.h"
#include "TH2D.h"
#include "EventType.h"
#include "TProfile.h"

double CalWeight(const anavar & evt, const int &partype);

double CalG4RW(const anavar & evt);

void FillHistVec1D(TH1D *hist[pi::nIntTypes+1], const double &value, const int &partype, double weight=1, bool fill_underflow=true, bool fill_overflow=true);

void FillHistVec2D(TH2D *hist[pi::nIntTypes+1], const double &value1, const double &value2, const int &partype, double weight=1);

void FillProfVec(TProfile *profile[pi::nIntTypes+1], const double &value1, const double &value2, const int &partype);

void FillHist1D(TH1D *hist, const double &value, const double &wei);

double GetPionKE(double length);

double GetTheta(double x, double y, double z);

double GetPhi(double x, double y, double z);

double GetThetaxz(double x, double y, double z);

double GetThetayz(double x, double y, double z);

void SetProtoDUNEStyle();

#endif
