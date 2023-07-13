#ifndef UTIL_H
#define UTIL_H

#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "EventType.h"
#include "TProfile.h"

double CalWeight(const anavar & evt, const int &partype, double rdm_radius, double rdm_angle);
double CalBkgW(const anavar & evt, double mu_weight=1., double p_weight=1., double spi_weight=1.);
double CalG4RW_ori(const anavar & evt, const double w1=1., const double w2=1.);
double CalG4RW(const anavar & evt, double weight[]);

void FillHistVec1D(TH1D *hist[pi::nIntTypes+1], const double &value, const int &partype, double weight=1, bool fill_underflow=true, bool fill_overflow=true);

void FillHistVec2D(TH2D *hist[pi::nIntTypes+1], const double &value1, const double &value2, const int &partype, double weight=1);

void FillHistVec3D(TH3D *hist[pi::nIntTypes+1], const double &value1, const double &value2, const double &value3, const int &partype, double weight=1);

void FillProfVec(TProfile *profile[pi::nIntTypes+1], const double &value1, const double &value2, const int &partype);

void FillHist1D(TH1D *hist, const double &value, const double &wei);

double GetPionKE(double length);

double GetTheta(double x, double y, double z);

double GetPhi(double x, double y, double z);

double GetThetaxz(double x, double y, double z);

double GetThetayz(double x, double y, double z);

void SetProtoDUNEStyle();

// from Heng-Ye
double chi2pid(std::vector<double> &trkdedx, std::vector<double> &trkres);

#endif
