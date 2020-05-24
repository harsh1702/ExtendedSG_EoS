#include "calib_liqNASGref.h"
#include "tools.h"

using namespace std;

double meanValue(vector<double> const& vec)
{
    // Purpose : compute the mean value of a vector
    double buf(0.);
    for (unsigned int i = 0; i < vec.size(); i++) {
        buf += vec[i];
    }
    return (buf/vec.size());
}

// **************************************************

// void computemeanvalues(double *mvL, double *mT, double *mp, double *mhL)
// {
//     mvL = meanValue(vLexp);
//     mT = meanValue(Texp);
//     mp = meanValue(psatExp);
//     mhL = meanValue(hLexp);
// }

// **************************************************

double computecvG(vector<double> const& vGexp, vector<double> const& Texp, vector<double> const& psatExp, double cpG)
{
    // Purpose : commpute heat capacity at constant volum of vapor phase with LSM
    // More : vGexp, Texp and psatExp are concording experimental points
    // See eq. (55)
    double num(0.), den(0.);
    for (unsigned int i = 0; i < Texp.size(); i++) {
        num += vGexp[i]*(Texp[i]/psatExp[i]);
        den += (Texp[i]/psatExp[i])*(Texp[i]/psatExp[i]);
    }
    return (cpG - num/den);
}

// **************************************************

double computeGammak(double cpk, double cvk)
{
    // Purpose : compute adiabatic index
    // See equation (25)
    return cpk/cvk;
}

// **************************************************

double computeMeanTp(vector<double> const& psatExp, vector<double> const& Texp, double pinfL)
{
	// pinfl will be assigned pinf1 so this is not specifically for liquid phase
    // Purpose : compute the mean value seen in (64) and (65)
    double mtp(0.);
    for (unsigned int i = 0; i < Texp.size(); i++) {
        mtp += Texp[i]/(psatExp[i]+pinfL);
    }
    mtp /= Texp.size();
    return mtp;
}

double computeMeanTp2(vector<double> const& psatExp, vector<double> const& Texp, double pinfL)
{
	// pinfl will be assigned pinf1 so this is not specifically for liquid phase
    // Purpose : compute a variation of the mean value seen in (64) and (65)
    double mtp(0.);
    for (unsigned int i = 0; i < Texp.size(); i++) {
        mtp += Texp[i]/((psatExp[i]+pinfL)*(psatExp[i]+pinfL));
    }
    mtp /= Texp.size();
    return mtp;
}

double computeHeatCapDiffL(vector<double> const& psatExp, vector<double> const& Texp, vector<double> const& vLexp, double mvL, double mTp, double pinfL)
{
    // Purpose : compute heat capacity cpL - cvL of liquid phase for pinfL iterative process with LSM
    // More : this fn is also used to solve eq. (68) with Newton-Raphson
    // Details : input values m'x' are the mean of the corresponding vector 'x'
    // See eq. (64)
    double num(0.), den(0.), bf(0.);
    for (unsigned int i = 0; i < Texp.size(); i++) {
        bf = Texp[i]/(psatExp[i]+pinfL);
        num += bf*(vLexp[i]-mvL);
        den += bf*(bf-mTp);
    }
    return num/den;
}

double computeDHeatCapDiffL(vector<double> const& psatExp, vector<double> const& Texp, vector<double> const& vLexp, double mvL, double mTp, double mTp2, double pinfL)
{
    // Purpose : compute pinfL derivative of fn computeHeatCapDiffL() for computePinfL Newton-Raphson process
    double m1(0.), m2(0.), m3(0.), m4(0.), bf1(0.), bf2(0.);
    for (unsigned int i = 0; i < Texp.size(); i++) {
        bf1 = Texp[i]/(psatExp[i]+pinfL);
        bf2 = bf1/(psatExp[i]+pinfL);
        m1 += -bf2*(vLexp[i]*mvL);
        m2 += bf1*(bf1 - (1./Texp.size())*mTp);
        m3 += bf1*(vLexp[i]-mvL);
        m4 += -bf2*(bf1-mTp)+bf1*(bf2+mTp2);
    }
    return (m1*m2-m3*m4)/(m2*m2);
}

// **************************************************

double computebL(double mvl, double mTp, double diffC)
{
    // Purpose : compute bl parameter for liquid phase
    // More : input values mvl, mT and mp are (m)ean values and diffC = cpL - cvL
    // Details : input values m'x' are the mean of the corresponding vector 'x'
    // See eq. (65)
    return (mvl-diffC*mTp);
}

double computeDbL(double mTp, double mTp2, double diffC, double dDiffC)
{
    // Purpose : compute pinfL derivative of fn computebL() for computePinfL Newton-Raphson process
    return (-dDiffC*mTp + diffC*mTp2);
}

// **************************************************

double computeCpL(vector<double> const& Texp, vector<double> const& hLexp, vector<double> const& psatExp, double mhL, double mp, double mT, double bL)
{
    // Purpose : compute heat cap. liquid phase at p = cst with LSM
    // More : this fn is used to solve eq. (68) with Newton-Raphson
    // Details : input values m'x' are the mean of the corresponding vector 'x'
    // See eq. (60)
    double num1(0.), num2(0.), num(0.), den(0.);
    for (unsigned int i = 0; i < Texp.size(); i++) {
        num1 += Texp[i]*(hLexp[i]-mhL);
        num2 += Texp[i]*(psatExp[i]-mp);
        den += Texp[i]*(Texp[i]-mT);
    }
    num = num1 - bL*num2;
    return num/den;
}

double computeDcpL(vector<double> const& Texp, vector<double> const& psatExp, double mp, double mT, double dbL)
{
    // Purpose : compute pinfL derivative fn of computeCpL() for computePinfL Newton-Raphson process
    double s1(0.), s2(0.);
    for (unsigned int i = 0; i < Texp.size(); i++) {
        s1 += Texp[i]*(psatExp[i]-mp);
        s2 += Texp[i]*(Texp[i]-mT);
    }
    return -dbL*(s1/s2);
}

// **************************************************

double computePinfL(vector<double> const& psatExp, vector<double> const& Texp, vector<double> const& vLexp, vector<double> const& hLexp, double p0L, double ro0L, double c0L)
{
    // Purpose : compute pinfL parameter of liquid phase using the Newton-Raphson procedure
    // More : liquid reference state and experimental data are used
    // See eq. (68)
    double fp, dfp, dfp1, dfp2, pinf1(2.e5), pinf2(0.), err(1.);
    double mp, mT, mvL, mhL, mTp, mTp2;
    double diffC, bL, cpL, dDiffC, dbL, dcpL;
    int count(0);

    mp = meanValue(psatExp);
    mT = meanValue(Texp);
    mvL = meanValue(vLexp);
    mhL = meanValue(hLexp);

    while (err > 1.e-5 && count < 100) {
        mTp = computeMeanTp(psatExp,Texp,pinf1);
        mTp2 = computeMeanTp2(psatExp,Texp,pinf1);

        diffC = computeHeatCapDiffL(psatExp,Texp,vLexp,mvL,mTp,pinf1);
        bL = computebL(mvL,mTp,diffC);
        cpL = computeCpL(Texp,hLexp,psatExp,mhL,mp,mT,bL);

        dDiffC = computeDHeatCapDiffL(psatExp,Texp,vLexp,mvL,mTp,mTp2,pinf1);
        dbL = computeDbL(mTp,mTp2,diffC,dDiffC);
        dcpL = computeDcpL(Texp,psatExp,mp,mT,dbL);

        dfp1 = -ro0L*dbL;
        dfp2 = (dDiffC*cpL-diffC*dcpL)/(cpL*cpL);

        fp = p0L + pinf1 - (1.-diffC/cpL)*ro0L*c0L*c0L*(1.-bL*ro0L);
        dfp = 1. - ro0L*c0L*c0L*dfp1 + dfp2*ro0L*c0L*c0L*(1.-bL*ro0L) + (diffC/cpL)*ro0L*c0L*c0L*dfp1;

        pinf2 = pinf1 - fp/dfp;
        err = fabs(pinf2-pinf1)/(0.5*(pinf1+pinf2));
        pinf1 = pinf2;
        count++;
        if (count >= 50) {
            cout << "Warning : newton-raphson of Psat(T) function not converged in PinfL\n"; exit(0);
        }
    }
    // cout << "Number of iteration NR pinfL : " << count << endl;
    return pinf1;
}
