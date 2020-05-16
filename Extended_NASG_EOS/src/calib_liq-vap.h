#ifndef CALIB_LIQVAP_EDITS_H
#define CALIB_LIQVAP_EDITS_H

#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <cmath>

using namespace std;

// Reading experimental data files
void readLiqInput(double &p0L, double &T0L, double &ro0L, double &e0L, double &brefLG);
void readVapInput(double &p0G, double &T0G, double &ro0G, double &e0G, double &brefG);
void readCritInput(double &Pc, double &Tc, double &vc, double &bc, double &pInfPrimeCrit);
void readAtmInput(double &pAtm, double &cAtm, double &vAtm);


//------For both phases-------
double computeb1(double bCrit, double brefL, double vCrit, double vrefL);
double computeb0(double brefL, double b1L, double vrefL);
double computeB(double b0L, double b1L);
double computefunction(vector<double> const& Pexp, vector<double> const& Texp, vector<double> const& vLexp, vector<double> const& eLexp, double pAtm, double cAtm, double vAtm, double p0, double T0, double e0, double b1, double b0, double pInfPrimeCrit, double Tc, double C);
//we have seperate computeC functions for gases and liquids because the initial value of C is different
//try to change the guessing of the C value 
double computeCL(vector<double> const& Pexp, vector<double> const& Texp, vector<double> const& vLexp, vector<double> const& eLexp, double pAtm, double cAtm, double vAtm, double p0, double T0, double e0, double b1, double b0, double pInfPrimeCrit, double Tc);
double computeCG(vector<double> const& Pexp, vector<double> const& Texp, vector<double> const& vLexp, vector<double> const& eLexp, double pAtm, double cAtm, double vAtm, double p0, double T0, double e0, double b1, double b0, double pInfPrimeCrit, double Tc);
double computeA(double CL, double Tc, double pInfPrimeCrit);
double computeSv1(vector<double> const& Pexp, vector<double> const& Texp, vector<double> const& vLexp, double AL, double BL, double CL, double b1L);
double computeSv2(vector<double> const& Pexp, vector<double> const& Texp, double AL, double CL, double b1L);
double computeSe1(vector<double> const& Pexp, vector<double> const& Texp, vector<double> const& eLexp, double p0L, double T0L, double e0L, double b1L, double AL, double CL);
double computeSe2(vector<double> const& Pexp, vector<double> const& Texp, double p0L, double T0L, double b1L, double AL, double CL);
double computeSe3(vector<double> const& Pexp, vector<double> const& Texp, vector<double> const& eLexp, double p0L, double T0L, double e0L, double b1L, double AL, double CL);
double computeSe4(vector<double> const& Pexp, vector<double> const& Texp, double p0L, double T0L, double b1L, double AL, double CL);
double computeSe5(vector<double> const& Pexp, vector<double> const& Texp, double p0L, double T0L, double b1L, double AL, double CL);
double computeD(double pressure, double temperature, double b1L, double AL, double CL);
double computeE(double pressure, double temperature, double b1L, double AL, double CL);
double computeF(double pressure, double temperature, double AL, double CL);
//these are function of my own method
// double computecv(vector<double> const& Pexp, vector<double> const& Texp, vector<double> const& eLexp, double p0L, double T0L, double eL0, double Sv1, double Sv2, double b1L, double AL, double CL);
// double computegamma(double cvL, double Sv1, double Sv2);
double computegamma(double Sv1, double Sv2, double Se1, double Se2, double Se3, double Se4, double Se5);
double computecv(double gamma, double Sv1, double Sv2);
double computePinf1(double gammaL, double AL);
double computePinf0(double gammaL, double b1L, double CL);
double computeq(double cvL, double p0L, double T0L, double e0L, double AL, double CL, double b1L, double gammaL);
double computeqPrime(vector<double> const& Pexp, vector<double> const& Texp, vector<double> const& sLexp, double cvL, double Pinf1L, double Pinf0L, double b1L, double gammaL);
double computeThEnthalpy(double cvL, double b1L, double b0L, double qL, double Pinf1L, double Pinf0L, double gammaL, double T, double P);
double computeVkTh(double cvL, double b1L, double b0L, double Pinf1L, double Pinf0L, double gammaL, double T, double P);
double computeThInternalEnergy(double cvL, double b1L, double qL, double Pinf1L, double Pinf0L, double gammaL, double T, double P);


#endif // CALIB_LIQVAP_EDITS_H
