#ifndef CALIB_LIQVAP_H
#define CALIB_LIQVAP_H

#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <cmath>

// Reading experimental data files
void readLiqVapInput(double &p0L, double &ro0L, double &c0L, double &p0G, double &ro0G, double &coG);

// --- Vapor phase ---
double computecpG(std::vector<double> const& hGexp, std::vector<double> const& ThGexp);
double computeQg(std::vector<double> const& hGexp, std::vector<double> const& ThGexp, double cpG);
double computecvG(std::vector<double> const& vGexp, std::vector<double> const& Texp, std::vector<double> const& psatExp, double cpG);
double computeQprimG(std::vector<double> p, std::vector<double> T, double cpL, double cpG, double cvL, double cvG, double qL, double qG, double pinfL, double pinfG, double bL);
//modifications to vapour phase
double computebG(double mvG, double mTp, double diffC);
double computeHeatCapDiffG(std::vector<double> const& psatExp, std::vector<double> const& Texp, std::vector<double> const& vLexp, double mvL, double mTp, double pinfG);
double computeDHeatCapDiffG(std::vector<double> const& psatExp, std::vector<double> const& Texp, std::vector<double> const& vGexp, double mvG, double mTp, double mTp2, double pinfG);
double computeCpG(std::vector<double> const& Texp, std::vector<double> const& hGexp, std::vector<double> const& psatExp, double mhG, double mp, double mT, double bG);
double computePinfG(std::vector<double> const& psatExp, std::vector<double> const& Texp, std::vector<double> const& vGexp, std::vector<double> const& hGexp, double p0, double ro0, double c0);
double computeDbG(double mTp, double mTp2, double diffC, double dDiffC);
double computeDcpG(std::vector<double> const& Texp, std::vector<double> const& psatExp, double mp, double mT, double dbG);


// --- Liquid phase ---
double computeMeanTp(std::vector<double> const& psatExp, std::vector<double> const& Texp, double pinfL);
double computeMeanTp2(std::vector<double> const& psatExp, std::vector<double> const& Texp, double pinfL);
double computeHeatCapDiffL(std::vector<double> const& psatExp, std::vector<double> const& Texp, std::vector<double> const& vLexp, double mvL, double mTp, double pinfL);
double computeDheatCapDiffL(std::vector<double> const& psatExp, std::vector<double> const& Texp, std::vector<double> const& vLexp, double mvL, double mTp, double mTp2, double pinfL);
double computebL(double mvl, double mTp, double diffC);
double computeDbL(double mTp, double mTp2, double diffC, double dDiffC);
double computeCpL(std::vector<double> const& Texp, std::vector<double> const& hLexp, std::vector<double> const& psatExp, double mhL, double mp, double mT, double bL);
double computeDcpL(std::vector<double> const& Texp, std::vector<double> const& psatExp, double mp, double mT, double dbL);
double computeQl(double mhL, double mT, double mp, double cpL, double bL);
double computePinfL(std::vector<double> const& psatExp, std::vector<double> const& Texp, std::vector<double> const& vLexp, std::vector<double> const& hLexp, double p0, double ro0, double c0);

#endif // CALIB_LIQVAP_H
