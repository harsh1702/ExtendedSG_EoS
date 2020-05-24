#ifndef CALIB_LIQVAPNASGREF_H
#define CALIB_LIQVAPNASGREF_H

#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <cmath>


double meanValue(std::vector<double> const& vec);
//void computemeanvalues(double *mvL, double *mT, double *mp, double *mhL);
double computeMeanTp(std::vector<double> const& psatExp, std::vector<double> const& Texp, double pinfL);
double computeMeanTp2(std::vector<double> const& psatExp, std::vector<double> const& Texp, double pinfL);
double computeHeatCapDiffL(std::vector<double> const& psatExp, std::vector<double> const& Texp, std::vector<double> const& vLexp, double mvL, double mTp, double pinfL);
double computeDheatCapDiffL(std::vector<double> const& psatExp, std::vector<double> const& Texp, std::vector<double> const& vLexp, double mvL, double mTp, double mTp2, double pinfL);
double computebL(double mvl, double mTp, double diffC);
double computeDbL(double mTp, double mTp2, double diffC, double dDiffC);
double computeCpL(std::vector<double> const& Texp, std::vector<double> const& hLexp, std::vector<double> const& psatExp, double mhL, double mp, double mT, double bL);
double computeDcpL(std::vector<double> const& Texp, std::vector<double> const& psatExp, double mp, double mT, double dbL);
double computePinfL(std::vector<double> const& psatExp, std::vector<double> const& Texp, std::vector<double> const& vLexp, std::vector<double> const& hLexp, double p0, double ro0, double c0);

// calib_liqNASGref *calib_liqNASGref;
//
// calib_liqNASGref *getcalib_liqNASGref()
// { return calib_liqNASGref; }

#endif // CALIB_LIQVAP_H
