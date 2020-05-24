#include "calib_liq-vap.h"
#include "tools.h"
#include <math.h>

using namespace std;

//Note that the arguements passed for all the functions are slightly wierd

    //4 files for inputing the refference state and critical values data
    //one for liquid state
    //one for vapurs state
    //one for critical values
    //one for atmospheric conditions

void readLiqInput(double &p0L, double &T0L, double &ro0L, double &e0L, double &c0L, double &brefL)
//void readLiqInput(double &p0L, double &T0L, double &ro0L, double &e0L, double &brefL)
{

    ifstream strmRefStates("input/refStateLiq.txt");
    string line("");
    if (strmRefStates) {
        for (int i=1; i<5; i++) {getline(strmRefStates,line);}
        p0L = stod(line);
        getline(strmRefStates,line); getline(strmRefStates,line);
        T0L = stod(line);
        getline(strmRefStates,line); getline(strmRefStates,line);
        ro0L = stod(line);
        getline(strmRefStates,line); getline(strmRefStates,line);
        e0L = stod(line);
        getline(strmRefStates,line); getline(strmRefStates,line);
        c0L = stod(line);
        getline(strmRefStates,line); getline(strmRefStates,line);
        brefL = stod(line);
    }
    else {
        cout << "Error : reading refStateLiq.txt file\n"; exit(0);
    }
}

// **************************************************

void readVapInput(double &p0G, double &T0G, double &ro0G, double &e0G, double &c0G, double &brefG)
//void readVapInput(double &p0G, double &T0G, double &ro0G, double &e0G, double &brefG)
{

    ifstream strmRefStates("input/refStateVap.txt");
    string line("");
    if (strmRefStates) {
        for (int i=1; i<5; i++) {getline(strmRefStates,line);}
        p0G = stod(line);
        getline(strmRefStates,line); getline(strmRefStates,line);
        T0G = stod(line);
        getline(strmRefStates,line); getline(strmRefStates,line);
        ro0G = stod(line);
        getline(strmRefStates,line); getline(strmRefStates,line);
        e0G = stod(line);
        getline(strmRefStates,line); getline(strmRefStates,line);
        c0G = stod(line);
        getline(strmRefStates,line); getline(strmRefStates,line);
        brefG = stod(line);
    }
    else {
        cout << "Error : reading refStateVap.txt file\n"; exit(0);
    }
}

// **************************************************

void readCritInput(double &Pc, double &Tc, double &vc, double &bc, double &pInfPrimeCrit)
{

    ifstream strmRefStates("input/refStateCrit.txt");
    string line("");
    if (strmRefStates) {
        for (int i=1; i<5; i++) {getline(strmRefStates,line);}
        Pc = stod(line);
        getline(strmRefStates,line); getline(strmRefStates,line);
        Tc = stod(line);
        getline(strmRefStates,line); getline(strmRefStates,line);
        vc = stod(line);
        getline(strmRefStates,line); getline(strmRefStates,line);
        bc = stod(line);
        getline(strmRefStates,line); getline(strmRefStates,line);
        pInfPrimeCrit = stod(line);
    }
    else {
        cout << "Error : reading refStateCrit.txt file\n"; exit(0);
    }
}

// **************************************************

void readAtmInput(double &pAtm, double &cAtm, double &vAtm)
{

    ifstream strmRefStates("input/refStateAtm.txt");
    string line("");
    if (strmRefStates) {
        for (int i=1; i<5; i++) {getline(strmRefStates,line);}
        pAtm = stod(line);
        getline(strmRefStates,line); getline(strmRefStates,line);
        cAtm = stod(line);
        getline(strmRefStates,line); getline(strmRefStates,line);
        vAtm = stod(line);
    }
    else {
        cout << "Error : reading refStateCrit.txt file\n"; exit(0);
    }
}

// **************************************************

double computeb1(double bCrit, double brefL, double vCrit, double vrefL)
{
  return ((bCrit-brefL)/(vCrit-vrefL));
}

// **************************************************

double computeb0(double brefL, double b1L, double vrefL)
{
  return (brefL - b1L*vrefL);
}

// **************************************************

double computeB(double b0L, double b1L)
{
  return (b0L/(1-b1L));
}

// **************************************************

double computefunction(vector<double> const& Pexp, vector<double> const& Texp, vector<double> const& vLexp, vector<double> const& eLexp, double pAtm, double cAtm, double vAtm, double p0, double T0, double e0, double b1, double b0, double pInfPrimeCrit, double Tc, double C)
{
  double gamma(0), cv(0), A(0), B(b0/(1-b1));
  double term3product1(0), term3product2(0);
  double Sv1(0), Sv2(0), Se1(0), Se2(0), Se3(0), Se4(0), Se5(0);

  //have to input the value of pInfPrimeCrit in the readfile or make it zero
  A = (pInfPrimeCrit - C)/Tc;
  //have to find gamma in terms of Se1, Se2....Se5 as Pinf0 is unknown
  //gamma = C*b1/(C - Pinf0*(1 - b1)); (this formula can't be used)
  Sv1 = computeSv1(Pexp, Texp, vLexp, A, B, C, b1);
  Sv2 = computeSv2(Pexp, Texp, A, C, b1);
  Se1 = computeSe1(Pexp, Texp, eLexp, p0, T0, e0, b1, A, C);
  Se2 = computeSe2(Pexp, Texp, p0, T0, b1, A, C);
  Se3 = computeSe3(Pexp, Texp, eLexp, p0, T0, e0, b1, A, B);
  Se4 = computeSe4(Pexp, Texp, p0, T0, b1, A, C);
  Se5 = computeSe5(Pexp, Texp, p0, T0, b1, A, B);
  gamma = computegamma(Sv1, Sv2, Se1, Se2, Se3, Se4, Se5);
  // std::cout << "gamma " <<gamma<< '\n';
  // if (isnan(gamma)==1)
  // {std::cout << "change input C value " << '\n';  return 0;}
  //gamma = 1.0147;

  //cv value in terms of cl will include Sv1 and Sv2
  //add a formula for cv involving the first equation in A27 with the other relations of A27
  cv = computecv(gamma, Sv1, Sv2);
  //std::cout << "cv " <<cv<< '\n';

//std::cout << "cv value from computeC " <<cv<< '\n';
//std::cout << "gamma from computeC " <<gamma<< '\n';
  //cv = 4014;

  //term2product1 is A

  term3product1 = ((pAtm + C)/((gamma - 1)*cv - A*(vAtm - b1*vAtm - b0)));
  term3product2 = (A*(gamma - 1) - ((gamma - b1)*((gamma - 1)*cv)/(vAtm - b1*vAtm - b0)));

  return (-(cAtm*cAtm) - A*vAtm*vAtm*(pAtm + C)/cv - vAtm*vAtm*term3product1*term3product2);
}

//here all 0 subscipted quantities are atmospheric conditions
//following equation A30 in the paper
// double computeCL(vector<double> const& Pexp, vector<double> const& Texp, vector<double> const& vLexp, vector<double> const& eLexp, double pAtm, double cAtm, double vAtm, double p0, double T0, double e0, double b1, double b0, double pInfPrimeCrit, double Tc)
double computeCL(vector<double> const& Pexp, vector<double> const& Texp, vector<double> const& vLexp, vector<double> const& eLexp, double pAtm, double cAtm, double vAtm, double p0, double T0, double e0, double b1, double b0, double pInfPrimeCrit, double Tc, double NASGpinfL)
{
  double count(0), f(0), f2(0), C(NASGpinfL), C2(0);
  double f3(0), C3(0);

  do {
    count++;

    f = computefunction(Pexp, Texp, vLexp, eLexp, pAtm, cAtm, vAtm, p0, T0, e0, b1, b0, pInfPrimeCrit, Tc, C);
    //using steffensens iterative method

    if (isnan(C)) {break;}
    if (isnan(f)) {break;}

    std::cout << "the value of C in CL " <<C<< '\n';
    std::cout << "the value of f " <<f<< '\n';


    if (abs(f) < 7e-3)  //the lowest observed value for f is 3.04682e-05 with C as 1.96225e+08
    {break;}

    //C2 is xn + f(xn)
    //f2 is f(xn + f(xn))
    C2 = C + f;
    f2 = computefunction(Pexp, Texp, vLexp, eLexp, pAtm, cAtm, vAtm, p0, T0, e0, b1, b0, pInfPrimeCrit, Tc, C2);
//for central difference scheme it gives better convergence
//without central difference scheme the value of C goes upto infity
    C3 = C - f;
    f3 = computefunction(Pexp, Texp, vLexp, eLexp, pAtm, cAtm, vAtm, p0, T0, e0, b1, b0, pInfPrimeCrit, Tc, C3);

    C -= 2*f*f/(f2 - f3);

  } while(1);//count < 100);

  std::cout << "the value of C in CL" <<C<< '\n';

  return C;
}

// **************************************************

double computeCG(vector<double> const& Pexp, vector<double> const& Texp, vector<double> const& vLexp, vector<double> const& eLexp, double pAtm, double cAtm, double vAtm, double p0, double T0, double e0, double b1, double b0, double pInfPrimeCrit, double Tc, double NASGpinfG)
//double computeCG(vector<double> const& Pexp, vector<double> const& Texp, vector<double> const& vLexp, vector<double> const& eLexp, double pAtm, double cAtm, double vAtm, double p0, double T0, double e0, double b1, double b0, double pInfPrimeCrit, double Tc)
{
  //will check the convergence of C for the given initial value here then apply it for the CL function
  double count(0), f(0), f2(0), C(NASGpinfG), C2(0);
  double f3(0), C3(0);

  do {
    count++;

    f = computefunction(Pexp, Texp, vLexp, eLexp, pAtm, cAtm, vAtm, p0, T0, e0, b1, b0, pInfPrimeCrit, Tc, C);
    //using steffensens iterative method

    if (isnan(C)) {break;}
    if (isnan(f)) {break;}

    if (abs(f) < 1e-5)  //the lowest observed value for f is 3.04682e-05 with C as 1.96225e+08
    {
      std::cout << "the value of C in CG " <<C<< '\n';
      std::cout << "the value of f in CG" <<f<< '\n';
      break;
    }

    // std::cout << "the value of C in CG " <<C<< '\n';
    // std::cout << "the value of f " <<f<< '\n';

    //C2 is xn + f(xn)
    //f2 is f(xn + f(xn))
    C2 = C + f;
    f2 = computefunction(Pexp, Texp, vLexp, eLexp, pAtm, cAtm, vAtm, p0, T0, e0, b1, b0, pInfPrimeCrit, Tc, C2);
//for central difference scheme it gives better convergence
//without central difference scheme the value of C goes upto infity
    C3 = C - f;
    f3 = computefunction(Pexp, Texp, vLexp, eLexp, pAtm, cAtm, vAtm, p0, T0, e0, b1, b0, pInfPrimeCrit, Tc, C3);
    --C;
    //C -= 2*f*f/(f2 - f3);

  } while(1);//count < 100);
  //convergence is not exactly being obtained in a hundred iterations, hence will keep it infinite for now

  std::cout << "the value of C in CG" <<C<< '\n';
  //C should be in aroound 308.76e6 in value

  return C;
}

// **************************************************

double computeA(double CL, double Tc, double pInfPrimeCrit)
{
  //we need the CL value to compute AL
  //AL seems fine
  return ((pInfPrimeCrit - CL)/Tc);
}

// **************************************************

double computeSv1(vector<double> const& Pexp, vector<double> const& Texp, vector<double> const& vLexp, double AL, double BL, double CL, double b1L)
{
  double numerator(0), denom(0), tot(0);

  for (unsigned int i = 0; i < Texp.size(); i++)
  {
    numerator = (vLexp[i] - BL)*Texp[i];
    denom = (1-b1L)*(Pexp[i] + AL*Texp[i] + CL);
    tot += numerator/denom;
  }

  return tot;
}

// **************************************************

double computeSv2(vector<double> const& Pexp, vector<double> const& Texp, double AL, double CL, double b1L)
{
  double numerator(0), denom(0), tot(0);

  for (unsigned int i = 0; i < Texp.size(); i++)
  {
    numerator = Texp[i]*Texp[i];
    denom = (1-b1L)*(1-b1L)*(Pexp[i] + AL*Texp[i] + CL)*(Pexp[i] + AL*Texp[i] + CL);
    tot += numerator/denom;
  }

  return tot;
}

// **************************************************

double computeSe1(vector<double> const& Pexp, vector<double> const& Texp, vector<double> const& eLexp, double p0L, double T0L, double e0L, double b1L, double AL, double CL)
{
  double numerator1(0), numerator2(0), denom1(0), denom2(0), tot(0);
  double E_ref(0), F_ref(0);
  double E(0), F(0);

  E_ref = computeE(p0L, T0L, b1L, AL, CL);
  F_ref = computeF(p0L, T0L, AL, CL);

  for (unsigned int i = 0; i < Texp.size(); i++)
  {
    E = computeE(Pexp[i], Texp[i], b1L, AL, CL);
    F = computeF(Pexp[i], Texp[i], AL, CL);

    numerator1 = (eLexp[i] - e0L)*(Pexp[i] + E);
    denom1 = F;
    numerator2 = (eLexp[i] - e0L)*(p0L + E_ref);
    denom2 = F_ref;
    tot += (numerator1/denom1 - numerator2/denom2);
  }

  return tot;
}

// **************************************************

double computeSe2(vector<double> const& Pexp, vector<double> const& Texp, double p0L, double T0L, double b1L, double AL, double CL)
{
  double numerator1(0), numerator2(0), denom1(0), denom2(0), tot(0);
  double E_ref(0), F_ref(0);
  double E(0), F(0);

  E_ref = computeE(p0L, T0L, b1L, AL, CL);
  F_ref = computeF(p0L, T0L, AL, CL);

  for (unsigned int i = 0; i < Texp.size(); i++)
  {
    E = computeE(Pexp[i], Texp[i], b1L, AL, CL);
    F = computeF(Pexp[i], Texp[i], AL, CL);

    numerator1 = (Pexp[i] + E);
    denom1 = F;
    numerator2 = (p0L + E_ref);
    denom2 = F_ref;
    tot += ((numerator1/denom1 - numerator2/denom2)*(numerator1/denom1 - numerator2/denom2));
  }

  return tot;
}

// **************************************************

double computeSe3(vector<double> const& Pexp, vector<double> const& Texp, vector<double> const& eLexp, double p0L, double T0L, double e0L, double b1L, double AL, double CL)
{
  double numerator1(0), numerator2(0), denom1(0), denom2(0), tot(0);
  double D_ref(0), F_ref(0);
  double D(0), F(0);

  D_ref = computeD(p0L, T0L, b1L, AL, CL);
  F_ref = computeF(p0L, T0L, AL, CL);

  for (unsigned int i = 0; i < Texp.size(); i++)
  {
    D = computeD(Pexp[i], Texp[i], b1L, AL, CL);
    F = computeF(Pexp[i], Texp[i], AL, CL);

    numerator1 = (eLexp[i] - e0L)*D;
    denom1 = F;
    numerator2 = (eLexp[i] - e0L)*D_ref;
    denom2 = F_ref;
    tot += (numerator1/denom1 - numerator2/denom2);
  }

  return tot;
}

// **************************************************

double computeSe4(vector<double> const& Pexp, vector<double> const& Texp, double p0L, double T0L, double b1L, double AL, double CL)
{
  double numerator1(0), numerator2(0), numerator3(0), denom1(0), denom2(0), denom3(0), tot(0);
  double D_ref(0), E_ref(0), F_ref(0);
  double D(0), E(0), F(0);

  D_ref = computeD(p0L, T0L, b1L, AL, CL);
  E_ref = computeE(p0L, T0L, b1L, AL, CL);
  F_ref = computeF(p0L, T0L, AL, CL);

  for (unsigned int i = 0; i < Texp.size(); i++)
  {
    D = computeD(Pexp[i], Texp[i], b1L, AL, CL);
    E = computeE(Pexp[i], Texp[i], b1L, AL, CL);
    F = computeF(Pexp[i], Texp[i], AL, CL);

    numerator1 = 2*D*(Pexp[i] + E);
    denom1 = F*F;
    numerator2 = 2*D*(Pexp[i] + E) + 2*D_ref*(p0L + E_ref);
    denom2 = F*F_ref;
    numerator3 = 2*D_ref*(p0L + E_ref);
    denom3 = F_ref*F_ref;
    tot += (numerator1/denom1 - numerator2/denom2 + numerator3/denom3);
  }

  return tot;
}

// **************************************************

double computeSe5(vector<double> const& Pexp, vector<double> const& Texp, double p0L, double T0L, double b1L, double AL, double CL)
{
  double tot(0);
  double D_ref(0), F_ref(0);
  double D(0), F(0);

  D_ref = computeD(p0L, T0L, b1L, AL, CL);
  F_ref = computeF(p0L, T0L, AL, CL);

  for (unsigned int i = 0; i < Texp.size(); i++)
  {
    D = computeD(Pexp[i], Texp[i], b1L, AL, CL);
    F = computeF(Pexp[i], Texp[i], AL, CL);

    tot += (D/F - D_ref/F_ref)*(D/F - D_ref/F_ref);
  }

  return tot;
}

// **************************************************

double computeD(double pressure, double temperature, double b1L, double AL, double CL)
{
  double numerator(0), denom(0);
  numerator = CL*(pressure + CL);
  denom = (1-b1L)*(pressure + AL*temperature + CL);
  return (numerator/denom);
}

// **************************************************

double computeE(double pressure, double temperature, double b1L, double AL, double CL)
{
  double term1(0), term2(0);
  term1 = (-CL*b1L)/(1-b1L);
  term2 = (AL*temperature*CL)/((1-b1L)*(pressure + AL*temperature + CL));
  return (term1 + term2);
}

// **************************************************

double computeF(double pressure, double temperature, double AL, double CL)
{
  //basically F is independant of AL
  return (((pressure + AL*temperature + CL)/temperature) - AL);
}

// **************************************************

// double computecv(vector<double> const& Pexp, vector<double> const& Texp, vector<double> const& eLexp, double p0L, double T0L, double eL0, double Sv1, double Sv2, double b1L, double AL, double CL)
// {
//   //done using equation A26
//   double DL(0), EL(0), FL(0), DL_ref(0), EL_ref(0), FL_ref(0);
//   double num1(0), num2(0), denom(0);
//   double cvL_Coeff1(0), cvL_Coeff2(0);
//
//   DL_ref = computeD(p0L, T0L, b1L, AL, CL);
//   EL_ref = computeE(p0L, T0L, b1L, AL, CL);
//   FL_ref = computeF(p0L, T0L, AL, CL);
//
//   //use method of least squares to find cvL from equation A26
//   for (unsigned int i = 0; i < Texp.size(); i++)
//   {
//     //do we need a DL array or does a single variable with name DL also work?
//     //put b1L and the other input parameters here
//     DL = computeD(Pexp[i], Texp[i], b1L, AL, CL);
//     EL = computeE(Pexp[i], Texp[i], b1L, AL, CL);
//     FL = computeF(Pexp[i], Texp[i], AL, CL);
//     // DL[i] = computeDL();
//     // EL[i] = computeEL();
//     // FL[i] = computeFL();
//     //in notes this is K1
//     cvL_Coeff1 = ((Pexp[i] + DL + EL)/FL - (p0L + DL_ref + EL_ref)/FL_ref);
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
// std::cout << "cv from computecv" <<((num1 - num2)/denom)<< '\n';
//
//   return ((num1 - num2)/denom);
//
// }
//
// // **************************************************
//
// double computegamma(double cvL, double Sv1, double Sv2)
// {
// std::cout << "gamma from computegamma" <<(Sv1/(cvL*Sv2) + 1)<< '\n';
//   return (Sv1/(cvL*Sv2) + 1);
// }

// **************************************************

double computegamma(double Sv1, double Sv2, double Se1, double Se2, double Se3, double Se4, double Se5)
{
  double gamma(0), sqrt_val(0), non_sqrt_val(0);
  // gamma = ((-Sv2*Se1 + Sv2*Se3 + Sv1*Se4) + sqrt((-Sv2*Se1 + Sv2*Se3 + Sv1*Se4)*(-Sv2*Se1 + Sv2*Se3 + Sv1*Se4) + 4*(Sv2*Se1 + Sv1*Se2)*(Sv2*Se3 - Sv1*Se5)))/(2*Sv2*Se3 - 2*Sv1*Se5);
  // std::cout << "gammaL value 1 " <<gamma<< '\n';
  gamma = ((-Sv2*Se1 + Sv2*Se3 + Sv1*Se4) - sqrt((-Sv2*Se1 + Sv2*Se3 + Sv1*Se4)*(-Sv2*Se1 + Sv2*Se3 + Sv1*Se4) + 4*(Sv2*Se1 + Sv1*Se2)*(Sv2*Se3 - Sv1*Se5)))/(2*Sv2*Se3 - 2*Sv1*Se5);
//  std::cout << "gammaL value 2 " <<gamma<< '\n';
  return gamma;
}

// **************************************************

double computecv(double gamma, double Sv1, double Sv2)
{
  return (Sv1/((gamma - 1)*Sv2));
}

// **************************************************

double computePinf1(double gammaL, double AL)
{
  //computing the coeffecient of the temperature term in Pinf
  return (AL/gammaL);
}

// **************************************************

double computePinf0(double gammaL, double b1L, double CL)
{
  double numerator(0),  denom(0);
  numerator = CL*(gammaL - b1L);
  denom = gammaL*(1 - b1L);

  return (numerator/denom);
}

// **************************************************

double computeq(double cvL, double p0L, double T0L, double e0L, double AL, double CL, double b1L, double gammaL)
{
  double DL_ref(0), EL_ref(0), FL_ref(0);

  DL_ref = computeD(p0L, T0L, b1L, AL, CL);
  EL_ref = computeE(p0L, T0L, b1L, AL, CL);
  FL_ref = computeF(p0L, T0L, AL, CL);

  return (e0L - cvL*(p0L + gammaL*DL_ref + EL_ref)/FL_ref);
}

// **************************************************

double computeqPrime(vector<double> const& Pexp, vector<double> const& Texp, vector<double> const& sLexp, double cvL, double Pinf1L, double Pinf0L, double b1L, double gammaL)
{
  double term1(0), term2(0), term3(0);
  double pInfPrime(0);

  for (unsigned int i = 0; i < Texp.size(); i++)
  {
    pInfPrime = gammaL*Pinf1L*Texp[i] + (gammaL*Pinf0L*(1 - b1L)/(gammaL - b1L));
    term1 += sLexp[i];
    term2 += cvL*log((pow(Texp[i],term1))/(pow((Pexp[i] + pInfPrime),term1)));
    term3 += (gammaL*Pinf1L*term1*cvL*Texp[i])/(Pexp[i] + pInfPrime);
  }

  return ((term1 + term2 + term3)/Texp.size());
}


// **************************************************
double computeThInternalEnergy(double cvL, double b1L, double qL, double Pinf1L, double Pinf0L, double gammaL, double T, double P)
{
  double pInf(0);
  double num1(0), denom1(0);

  pInf = Pinf1L*T + Pinf0L;

  num1 = (P + gammaL*pInf)*cvL*T;
  denom1 = P + gammaL*Pinf1L*T + (gammaL*Pinf0L*(1 - b1L)/(gammaL - b1L));

  return (num1/denom1 + qL);
}

// **************************************************

double computeThEnthalpy(double cvL, double b1L, double b0L, double qL, double Pinf1L, double Pinf0L, double gammaL, double T, double P)
{
    // Purpose : compute phasic theoric enthalpy
    // See eq. (56)

    double pInf(0);
    double denom1(0), denom2(0), num1(0), num2(0);

    pInf = Pinf1L*T + Pinf0L;

    num1 = cvL*T*(gammaL*(P + pInf) - P*b1L - gammaL*b1L*pInf);
    denom1 = (1 - b1L)*(P + gammaL*Pinf1L*T + (gammaL*Pinf0L*(1 - b1L)/(gammaL - b1L)));
    num2 = P*b0L;
    denom2 = 1 - b1L;

    return (num1/denom1 + num2/denom2 + qL);
}

// **************************************************

double computeVkTh(double cvL, double b1L, double b0L, double Pinf1L, double Pinf0L, double gammaL, double T, double P)
{
    // Purpose : compute theoric specific vol. of phase k
    // See eq. (56)

    double denom1(0), denom2(0), num1(0), num2(0);

    num1 = (gammaL - 1)*cvL*T;
    denom1 = (1 - b1L)*(P + gammaL*Pinf1L*T + (gammaL*Pinf0L*(1 - b1L)/(gammaL - b1L)));
    num2 = b0L;
    denom2 = 1 - b1L;

    return (num1/denom1 + num2/denom2);
}
