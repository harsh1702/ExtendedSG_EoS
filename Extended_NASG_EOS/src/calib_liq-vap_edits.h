// #ifndef CALIB_LIQVAP_EDITS_H
// #define CALIB_LIQVAP_EDITS_H
//
// #include <iostream>
// #include <string>
// #include <fstream>
// #include <vector>
// #include <cmath>
//
// // Reading experimental data files
// void readLiqInput(double &p0L, double &T0L, double &ro0L, double &e0L, double &brefLG);
// void readVapInput(double &p0G, double &T0G, double &ro0G, double &e0G, double &brefG);
// void readCritInput(double &Pc, double &Tc, double &vc, double &bc);
// void readAtmInput(double pAtm, double cAtm, double vAtm);
//
//
// //------For both phases-------
// double computeb1(double bCrit, double brefL, double vCrit, double vrefL);
// double computeb0(double brefL, double b1L, double vrefL);
// double computeB(double b0L, double b1L);
// double computeC(double p0, double c0, double v0, double Sv1, double Sv2, double b1, double b0, double Pinf0, double Pinf1, double pInfPrimeCrit);
// double computeA(double CL, double Tc);
// double computeSv1(vector<double> const& Pexp, vector<double> const& Texp, vector<double> const& vLexp, double AL, double BL, double CL);
// double computeSv2(vector<double> const& Pexp, vector<double> const& Texp, double AL, double BL, double CL);
// double computeD(double pressure, double temperature, double b1L, double AL, double CL);
// double computeE(double pressure, double temperature, double b1L, double AL, double CL);
// double computeF(double pressure, double temperature, double b1L, double AL, double CL);
// double computecv(vector<double> const& Pexp, vector<double> const& Texp, vector<double> const& vLexp, vector<double> const& eLexp, double p0L, double T0L, double ro0L, double c0L, double eL0, double Sv1, double Sv2, double b1L);
// double computegamma(double cvL, double Sv1, double Sv2);
// double computePinf1(double gammaL, double AL);
// double computePinf0(double gammaL, double b1L, double CL);
// double computeq(double cvL, double p0L, double e0L);
// double computeqPrime(vector<double> const& Pexp, vector<double> const& Texp, vector<double> const& sLexp, double p0L, double T0L, double ro0L, double c0L, double eL0, double Sv1, double Sv2);
// double computeThEnthalpy(double cvL, double b1L, double b0L, double qL, double Pinf1L, double Pinf0L , double T, double P);
// double computeVkTh(double cvL, double b1L, double b0L, double Pinf1L, double Pinf0L, double T, double P);
//
//
// #endif // CALIB_LIQVAP_EDITS_H
