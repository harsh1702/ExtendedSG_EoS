// #include "calib_liq-vap_edits.h"
// #include "tools.h"
// #include <math.h>
//
// using namespace std;
//
// //Note that the arguements passed for all the functions are slightly wierd
//
//     //3 files for inputing the refference state and critical values data
//     //one for liquid state
//     //one for vapurs state
//     //one for critical values
//     //need to input much more data
//
// void readLiqInput(double &p0L, double &T0L, double &ro0L, double &e0L, double &brefLG)
// {
//
//     ifstream strmRefStates("input/refStateLiq.txt");
//     string line("");
//     if (strmRefStates) {
//         for (int i=1; i<5; i++) {getline(strmRefStates,line);}
//         p0L = stod(line);
//         getline(strmRefStates,line); getline(strmRefStates,line);
//         T0L = stod(line);
//         getline(strmRefStates,line); getline(strmRefStates,line);
//         ro0L = stod(line);
//         getline(strmRefStates,line); getline(strmRefStates,line);
//         e0L = stod(line);
//         getline(strmRefStates,line); getline(strmRefStates,line);
//         brefL = stod(line);
//     }
//     else {
//         cout << "Error : reading refStateLiq.txt file\n"; exit(0);
//     }
// }
//
// // **************************************************
//
// void readVapInput(double &p0G, double &T0G, double &ro0G, double &e0G, double &brefG)
// {
//
//     ifstream strmRefStates("input/refStateVap.txt");
//     string line("");
//     if (strmRefStates) {
//         for (int i=1; i<5; i++) {getline(strmRefStates,line);}
//         p0G = stod(line);
//         getline(strmRefStates,line); getline(strmRefStates,line);
//         T0G = stod(line);
//         getline(strmRefStates,line); getline(strmRefStates,line);
//         ro0G = stod(line);
//         getline(strmRefStates,line); getline(strmRefStates,line);
//         e0G = stod(line);
//         getline(strmRefStates,line); getline(strmRefStates,line);
//         brefG = stod(line);
//     }
//     else {
//         cout << "Error : reading refStateVap.txt file\n"; exit(0);
//     }
// }
//
// // **************************************************
//
// void readCritInput(double &Pc, double &Tc, double &vc, double &bc)
// {
//
//     ifstream strmRefStates("input/refStateCrit.txt");
//     string line("");
//     if (strmRefStates) {
//         for (int i=1; i<5; i++) {getline(strmRefStates,line);}
//         Pc = stod(line);
//         getline(strmRefStates,line); getline(strmRefStates,line);
//         Tc = stod(line);
//         getline(strmRefStates,line); getline(strmRefStates,line);
//         vc = stod(line);
//         getline(strmRefStates,line); getline(strmRefStates,line);
//         bc = stod(line);
//         getline(strmRefStates,line); getline(strmRefStates,line);
//         pInfPrimeCrit = stod(line);
//     }
//     else {
//         cout << "Error : reading refStateCrit.txt file\n"; exit(0);
//     }
// }
//
// // **************************************************
//
// void readAtmInput(double &pAtm, double &cAtm, double &vAtm)
// {
//
//     ifstream strmRefStates("input/refStateAtm.txt");
//     string line("");
//     if (strmRefStates) {
//         for (int i=1; i<5; i++) {getline(strmRefStates,line);}
//         pAtm = stod(line);
//         getline(strmRefStates,line); getline(strmRefStates,line);
//         cAtm = stod(line);
//         getline(strmRefStates,line); getline(strmRefStates,line);
//         vAtm = stod(line);
//     }
//     else {
//         cout << "Error : reading refStateCrit.txt file\n"; exit(0);
//     }
// }
//
// // **************************************************
//
// double computeb1(double bCrit, double brefL, double vCrit, double vrefL)
// {
//   return ((bCrit-brefL)/(vCrit-vrefL));
// }
//
// // **************************************************
//
// double computeb0(double brefL, double b1L, double vrefL)
// {
//   return (brefL - b1L*vrefL);
// }
//
// // **************************************************
//
// double computeB(double b0L, double b1L)
// {
//   return (b0L/(1-b1L));
// }
//
// // **************************************************
//
// //here all 0 subscipted quantities are atmospheric conditions
// //following equation A30 in the paper
// double computeC(double p0, double c0, double v0, double Sv1, double Sv2, double b1, double b0, double Pinf0, double Pinf1, double pInfPrimeCrit)
// {
//   double count(0);
//   double gamma(0), f(0), A(0), C(0), term2product1(0), term2product2(0), term2product3(0), term3product1(0), term3product2(0);
//   double dgamma(0), df(0), dterm2product1(0), dterm2product2(0), dterm2product3(0), dterm3product1(0), dterm3product2(0);
//
//   do {
//     count++;
//
//     //have to input the value of pInfPrimeCrit in the readfile or make it zero
//     A = (pInfPrimeCrit - C)/Tc;
//     gamma = C*b1/(C - Pinf0*(1 - b1));
//
//     //term2product1 is A
//     term2product2 = p0 + C;
//     term2product3 = Sv2*gamma/Sv1;
//
//     term3product1 = ((p0 + C)/(Sv1/Sv2 - A*(v0 - b1*v0 - b0)));
//     term3product2 = (A*Sv1/Sv2 - (gamma - b1)*Sv1/(Sv2*(v0 - b1*v0 - b0)));
//
//     f = -(c0*c0) - A*v0*v0*(p0 + C)*Sv2*gamma/Sv1 - v0*v0*term3product1*term3product2;
//
//     //derivative of gamma wrt C
//     dgamma = (b1*(C - Pinf0*(1 - b1)) - C*b1)/((C - Pinf0*(1 - b1))*(C - Pinf0*(1 - b1)));
//     //derivative of term2product1 wrt C basically dA
//     dterm2product1 = -1/Tc;
//     //derivative of term2product2 wrt C
//     dterm2product2 = 1;
//     //derivative of term2product3 wrt C
//     dterm2product3 = (Sv2/Sv1)*dgamma;
//     //derivative of term3product1 wrt C
//     dterm3product1 = (Sv1/Sv2 - A*(v0 - b1*v0 - b0) - ((p0 + C)*(v0 - b1*v0 - b0)/Tc);
//     //derivative of term3product2 wrt C
//     dterm3product2 = -(gamma - 1)/Tc + A*dgamma - dgamma*Sv1/(Sv2*(v0 - b1*v0 - b0)
//     //derivative of f wrt C for newton rhapson method
//     df = -v0*v0*(dterm2product1*term2product2*term2product3 + term2product1*dterm2product2*term2product3 + term2product1*term2product2*dterm2product3) - v0*v0*(dterm3product1*term3product2 + term3product1*dterm3product2);
//
//     C -= f/df;
//   } while(abs(f) < 1e-5 && count < 100);
// }
//
// // **************************************************
//
// double computeA(double CL, double Tc)
// {
//   //we need the CL value to compute AL
//   return (CL/Tc);
// }
//
// // **************************************************
//
// double computeSv1(vector<double> const& Pexp, vector<double> const& Texp, vector<double> const& vLexp, double AL, double BL, double CL)
// {
//   double numerator(0), denom(0), tot(0);
//
//   for (unsigned int i = 0; i < Texp.size(); i++)
//   {
//     numerator = (vLexp[i] - BL)*Texp[i];
//     denom = (1-b1L)*(Pexp[i] + AL*Texp[i] + CL);
//     tot += numerator/denom;
//   }
//
//   return tot;
// }
//
// // **************************************************
//
// double computeSv2(vector<double> const& Pexp, vector<double> const& Texp, double AL, double BL, double CL)
// {
//   double numerator(0), denom(0), tot(0);
//
//   for (unsigned int i = 0; i < Texp.size(); i++)
//   {
//     numerator = Texp[i]*Texp[i];
//     denom = (1-b1L)*(1-b1L)*(Pexp[i] + AL*Texp[i] + CL)*(Pexp[i] + AL*Texp[i] + CL);
//     tot += numerator/denom;
//   }
//
//   return tot;
// }
//
// // **************************************************
//
// double computeD(double pressure, double temperature, double b1L, double AL, double CL)
// {
//   double numerator(0), denom(0);
//   numerator = CL*(pressure + CL);
//   denom = (1-b1L)*(pressure + AL*temperature + CL);
//   return (numerator/denom);
// }
//
// // **************************************************
//
// double computeE(double pressure, double temperature, double b1L, double AL, double CL)
// {
//   double term1(0), term2(0);
//   term1 = (-CL*b1L)/(1-b1L);
//   term2 = (AL*temperature*CL)/((1-b1L)*(pressure + AL*temperature + CL);
//   return (term1 + term2);
// }
//
// // **************************************************
//
// double computeF(double pressure, double temperature, double b1L, double AL, double CL)
// {
//   return ((pressure + AL*temperature + CL)/temperature);
// }
//
// // **************************************************
//
// double computecv(vector<double> const& Pexp, vector<double> const& Texp, vector<double> const& vLexp, vector<double> const& eLexp, double p0L, double T0L, double ro0L, double c0L, double eL0, double Sv1, double Sv2, double b1L)
// {
//   double DL(O), EL(0), FL(0), DL_ref(0), EL_ref(0), FL_ref(0);
//   double num1(0), num2(0), denom(0);
//   double cvL_Coeff1(0), double cvL_Coeff2(0);
//
//   DL_ref = computeD(p0L, T0L, b1L, AL, CL);
//   EL_ref = computeE(p0L, T0L, b1L, AL, CL);
//   FL_ref = computeF(p0L, T0L, b1L, AL, CL);
//
//   //use method of least squares to find cvL from equation A26
//   for (unsigned int i = 0; i < Texp.size(); i++)
//   {
//     //do we need a DL array or does a single variable with name DL also work?
//     //put b1L and the other input parameters here
//     DL = computeD(Pexp[i], Texp[i], b1L, AL, CL);
//     EL = computeE(Pexp[i], Texp[i], b1L, AL, CL);
//     FL = computeF(Pexp[i], Texp[i], b1L, AL, CL);
//     // DL[i] = computeDL();
//     // EL[i] = computeEL();
//     // FL[i] = computeFL();
//     //in notes this is K1
//     cvL_Coeff1 = ((Pexp[i] + DL[i] + EL[i])/FL[i] - (p0L + DL_ref + EL_ref)/FL_ref);
//     //in notes this is K2
//     cvL_Coeff2 = (Sv1/Sv2)*((DL/FL) - (DL_ref/FL_ref));
//     //first term in the numerator
//     num1 += (eLexp[i]*(cvL_Coeff1 + cvL_Coeff2));
//     //second term in the numerator eL0 is the reference internal energy
//     num2 += (eL0*(cvL_Coeff1 + cvL_Coeff2));
//     //denominator
//     denom += (cvL_Coeff1 + cvL_Coeff2)*(cvL_Coeff1 + cvL_Coeff2);
//   }
//
//   return ((num1 - num2)/denom);
//
// }
//
// // **************************************************
//
// double computegamma(double cvL, double Sv1, double Sv2)
// {
//   return (Sv1/(cvL*Sv2) + 1);
// }
//
// // **************************************************
//
// double computePinf1(double gammaL, double AL)
// {
//   //computing the coeffecient of the temperature term in Pinf
//   return (AL/gammaL);
// }
//
// // **************************************************
//
// double computePinf0(double gammaL, double b1L, double CL)
// {
//   double numerator(0), double denom(0);
//   numerator = CL*(gammaL - b1L);
//   denom = gammaL*(1 - b1L);
//
//   return (numerator/denom);
// }
//
// // **************************************************
//
// double computeq(double cvL, double p0L, double e0L)
// {
//   double DL_ref(0), EL_ref(0), FL_ref(0);
//
//   DL_ref = computeD(p0L, T0L, b1L, AL, CL);
//   EL_ref = computeE(p0L, T0L, b1L, AL, CL);
//   FL_ref = computeF(p0L, T0L, b1L, AL, CL);
//
//   return (e0L - cvL*(p0L + DL_ref + EL_ref)/FL_ref));
// }
//
// // **************************************************
//
// double computeqPrime(vector<double> const& Pexp, vector<double> const& Texp, vector<double> const& sLexp, double p0L, double T0L, double ro0L, double c0L, double eL0, double Sv1, double Sv2)
// {
//   double term(0), term1(0), term2(0), term3(0);
//   double pInfPrime(0);
//
//   term0 = (gammaL - 1)/(1 - b1L);
//
//   for (unsigned int i = 0; i < Texp.size(); i++)
//   {
//     pInfPrime = gammaL*Pinf1L*Texp[i] + (gammaL*Pinf0L*(1 - b1L)/(gammaL - b1L));
//     term1 += sLexp[i];
//     term2 += cvL*Log((pow(Texp[i],term1))/(pow((Pexp[i] + pInfPrime),term1)));
//     term3 += (gammaL*Pinf1L*term1*cvL*Texp[i])/(Pexp[i] + pInfPrime);
//   }
//
//   return ((term1 + term2 + term3)/Texp.size());
// }
//
// // **************************************************
//
// double computeThEnthalpy(double cvL, double b1L, double b0L, double qL, double Pinf1L, double Pinf0L , double T, double P)
// {
//     // Purpose : compute phasic theoric enthalpy
//     // See eq. (56)
//
//     double pInf(0);
//     double denom1(0), denom2(0), num1(0), num2(0);
//
//     PinfL = Pinf1L*T + Pinf0L;
//
//     num1 = cvL*T*(gammaL*(P + PinfL) - P*b1L - gammaL*b1L*PinfL);
//     denom1 = (P + gammaL*Pinf1L*Texp[i] + (gammaL*Pinf0L*(1 - b1L)/(gammaL - b1L)))*(1 - b1L);
//     num2 = P*b0L;
//     denom2 = 1 - b1L;
//
//     return (num1/denom1 + num2/denom2 + qL);
// }
//
// // **************************************************
//
// double computeVkTh(double cvL, double b1L, double b0L, double Pinf1L, double Pinf0L, double T, double P)
// {
//     // Purpose : compute theoric specific vol. of phase k
//     // See eq. (56)
//
//     double denom1(0), denom2(0), num1(0), num2(0);
//
//     num1 = (gammaL - 1)*cvL*T;
//     denom1 = (P + gammaL*Pinf1L*Texp[i] + (gammaL*Pinf0L*(1 - b1L)/(gammaL - b1L)))*(1 - b1L);
//     num2 = b0L;
//     denom2 = 1 - b1L;
//
//     return (num1/denom1 + num2/denom2);
// }
