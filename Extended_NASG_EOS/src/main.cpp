#include <iostream>
#include <string>

#include "tools.h"
#include "calib_liq-vap.h"
#include "calib_liqNASGref.h"

using namespace std;

int main()
{
    string run("res/");
    //here it is assumed that pinfprime at the critical point is 0
    //0 represents the refernce state
    double p0L, ro0L, e0L, T0L, c0L, brefL;
    double p0G, ro0G, e0G, T0G, c0G, brefG;
    double Pc, Tc, vc, bc, pInfPrimeCrit;

    //cP is not required
    double cpL(0.), qL(0.), Pinf1L(0.), Pinf0L(0.), cvL(0.), gammaL(0.), qPrimeL(0.), b1L(0.), b0L(0.);
    double cpG(0.), qG(0.), Pinf1G(0.), Pinf0G(0.), cvG(0.), gammaG(0.), qPrimeG(0.), b1G(0.), b0G(0.);
    //double mvL(0.), mT(0.), mp(0.), mhL(0.), mTp(0.), diffCl(0.), mvG(0.), mhG(0.),diffCg(0.);
    vector<double> Texp, Pexp, vGexp, vLexp, hGexp, hLexp, LvExp, eGexp, eLexp, sGexp, sLexp;

    double AL(0), BL(0), CL(0);
    double Sv1L(0), Sv2L(0);
    double Se1L(0), Se2L(0), Se3L(0), Se4L(0), Se5L(0);

    double AG(0), BG(0), CG(0);
    double Sv1G(0), Sv2G(0);
    double Se1G(0), Se2G(0), Se3G(0), Se4G(0), Se5G(0);

    //atmospheric conditions used for finding CL value iteratively
    double pAtm(0), vAtm(0), cAtm(0);

    //NASG parameters
    double NASGpinfL(0), NASGpinfG(0);

	//currently reading the reference state for liquid and gas properties
  //add the input of critical properties also
  //add the input of atmospheric propeties also
    readLiqInput(p0L,T0L,ro0L,e0L,c0L,brefL);
    readVapInput(p0G,T0G,ro0G,e0G,c0G,brefG);
    // readLiqInput(p0L,T0L,ro0L,e0L,brefL);
    // readVapInput(p0G,T0G,ro0G,e0G,brefG);
    readCritInput(Pc,Tc,vc,bc,pInfPrimeCrit);
    readAtmInput(pAtm,cAtm,vAtm);

    double v0L(1/ro0L), v0G(1/ro0G);

    //reading the saturation properties of the substance
    //have to change this also
    readExpData("input/expData.txt",Texp,Pexp,vGexp,vLexp,hGexp,hLexp,LvExp,eGexp,eLexp,sGexp,sLexp);


    //finding the Pinf for the NASG case to imporve the guess for C
    //change the input values of P and T to make it saturation properties
    //edit the reference states speed of sound must be added

    // --- Liquid phase ---
    NASGpinfL = computePinfL(Pexp,Texp,vLexp,hLexp,p0L,ro0L,c0L);
    // --- Gaseous phase ---
    NASGpinfG = computePinfL(Pexp,Texp,vGexp,hGexp,p0G,ro0G,c0G);


    // --- Liquid phase ---
    b1L = computeb1(bc, brefL, vc, v0L);
    b0L = computeb0(brefL, b1L, v0L);
    BL = computeB(b0L, b1L);
    //CL function
    //Sv1 cant be an arguement
    CL = computeCL(Pexp,Texp,vLexp,eLexp,pAtm,cAtm,vAtm,p0L,T0L,e0L,b1L,b0L,pInfPrimeCrit,Tc,NASGpinfL);
    //CL = computeCL(Pexp,Texp,vLexp,eLexp,pAtm,cAtm,vAtm,p0L,T0L,e0L,b1L,b0L,pInfPrimeCrit,Tc);
    //CL = 208.76e6;
    AL = computeA(CL,Tc,pInfPrimeCrit);
    Sv1L = computeSv1(Pexp, Texp, vLexp, AL, BL, CL, b1L);
    Sv2L = computeSv2(Pexp, Texp, AL, CL, b1L);
    Se1L = computeSe1(Pexp, Texp, eLexp, p0L, T0L, e0L, b1L, AL, CL);
    Se2L = computeSe2(Pexp, Texp, p0L, T0L, b1L, AL, CL);
    Se3L = computeSe3(Pexp, Texp, eLexp, p0L, T0L, e0L, b1L, AL, BL);
    Se4L = computeSe4(Pexp, Texp, p0L, T0L, b1L, AL, CL);
    Se5L = computeSe5(Pexp, Texp, p0L, T0L, b1L, AL, BL);

    // cvL = computecv(Pexp, Texp, eLexp, p0L, T0L, e0L, Sv1L, Sv2L, b1L, AL, CL);
    // gammaL = computegamma(cvL, Sv1L, Sv2L);
    gammaL = computegamma(Sv1L, Sv2L, Se1L, Se2L, Se3L, Se4L, Se5L);
    //std::cout << "gamma L value " <<gammaL<< '\n';
    cvL = computecv(gammaL, Sv1L, Sv2L);
    cpL = gammaL*cvL;
    Pinf1L = computePinf1(gammaL, AL);
    Pinf0L = computePinf0(gammaL, b1L, CL);

    qL = computeq(cvL, p0L, T0L, e0L, AL, CL, b1L, gammaL);
    qPrimeL = computeqPrime(Pexp, Texp, sLexp, cvL, Pinf1L, Pinf0L, b1L, gammaL);

    //mTp = computeMeanTp(psatExp,Texp,pinfL);
    //diffCl = computeHeatCapDiffL(psatExp,Texp,vLexp,mvL,mTp,pinfL);

    // --- Vapor phase ---
    b1G = computeb1(bc, brefG, vc, v0G);
    b0G = computeb0(brefG, b1G, v0G);
    BG = computeB(b0G, b1G);
    //CL function
    //Sv1 cant be an arguement
    CG = computeCG(Pexp,Texp,vGexp,eGexp,pAtm,cAtm,vAtm,p0G,T0G,e0G,b1G,b0G,pInfPrimeCrit,Tc,NASGpinfG);
    //CG = computeCG(Pexp,Texp,vGexp,eGexp,pAtm,cAtm,vAtm,p0G,T0G,e0G,b1G,b0G,pInfPrimeCrit,Tc);

    AG = computeA(CG,Tc,pInfPrimeCrit);
    Sv1G = computeSv1(Pexp, Texp, vGexp, AG, BG, CG, b1G);
    Sv2G = computeSv2(Pexp, Texp, AG, CG, b1G);
    Se1G = computeSe1(Pexp, Texp, eGexp, p0G, T0G, e0G, b1G, AG, CG);
    Se2G = computeSe2(Pexp, Texp, p0G, T0G, b1G, AG, CG);
    Se3G = computeSe3(Pexp, Texp, eGexp, p0G, T0G, e0G, b1G, AG, BG);
    Se4G = computeSe4(Pexp, Texp, p0G, T0G, b1G, AG, CG);
    Se5G = computeSe5(Pexp, Texp, p0G, T0G, b1G, AG, BG);

    // cvL = computecv(Pexp, Texp, eLexp, p0L, T0L, e0L, Sv1L, Sv2L, b1L, AL, CL);
    // gammaL = computegamma(cvL, Sv1L, Sv2L);
    gammaG = computegamma(Sv1G, Sv2G, Se1G, Se2G, Se3G, Se4G, Se5G);
    //std::cout << "gamma G value " <<gammaG<< '\n';
    cvG = computecv(gammaG, Sv1G, Sv2G);
    cpG = gammaG*cvG;
    Pinf1G = computePinf1(gammaG, AG);
    Pinf0G = computePinf0(gammaG, b1G, CG);

    qG = computeq(cvG, p0G, T0G, e0G, AG, CG, b1G, gammaG);
    qPrimeG = computeqPrime(Pexp, Texp, sGexp, cvG, Pinf1G, Pinf0G, b1G, gammaG);


    // --- Results ---
    cout << "--- Liquid (L) ---\n";
    cout << "cpL    (J.kg-1.K-1)  : " << cpL << endl;
    cout << "qL     (J.kg-1)      : " << qL << endl;
    cout << "Pinf1L  (Pa)          : " << Pinf1L << endl;
    cout << "Pinf0L  (Pa)          : " << Pinf0L << endl;
    cout << "cvL    (J.kg-1.K-1)  : " << cvL << endl;
    cout << "gammaL (-)           : " << gammaL << endl;
    cout << "q'L    (J.kg-1)      : " << qPrimeL << endl;
    cout << "b1L     (m3/kg)       : " << b1L << endl;
    cout << "b0L     (m3/kg)       : " << b0L << endl;
    cout << "\n";

    cout << "--- Gas (G) ---\n";
    cout << "cpG    (J.kg-1.K-1)  : " << cpG << endl;
    cout << "qG     (J.kg-1)      : " << qG << endl;
    cout << "Pinf1G  (Pa)          : " << Pinf1G << endl;
    cout << "Pinf0G  (Pa)          : " << Pinf0G << endl;
    cout << "cvG    (J.kg-1.K-1)  : " << cvG << endl;
    cout << "gammaG (-)           : " << gammaG << endl;
    cout << "q'G    (J.kg-1)      : " << qPrimeG << endl;
    cout << "b1G     (m3/kg)       : " << b1G << endl;
    cout << "b0G     (m3/kg)       : " << b0G << endl;
    cout << "\n";

    //writeResults(cpL,qL,pinfL,cvL,gammaL,qPrimL,bL,cpG,qG,PinfG,cvG,gammaG,qPrimG,bG);

    // --- Write theoric curves ---
    vector<double> hLth, hGth, LvTh, PsatTh, vlTh, vgTh, elTh, egTh;
//    double A(0.),B(0.),C(0.),D(0.),E(0.);

//    coeffPsatTh(cpG,cpL,cvG,cvL,qG,qL,qPrimG,qPrimL,bG,bL,A,B,C,D,E);

    for (unsigned int i = 0; i < Texp.size(); i++) {
        //calculates the Psat values only before the critical temperature
        //yet to implement the Psat and Tsat relationship
//        if (Texp[i] < Tc) {PsatTh.push_back(computePsatTh(A,B,C,D,E,pinfL,Texp[i]));}
        hLth.push_back(computeThEnthalpy(cvL,b1L,b0L,qL,Pinf1L,Pinf0L,gammaL,Texp[i],Pexp[i]));
        hGth.push_back(computeThEnthalpy(cvG,b1G,b0G,qG,Pinf1L,Pinf0L,gammaG,Texp[i],Pexp[i]));
        LvTh.push_back(hGth[i]-hLth[i]);
        vlTh.push_back(computeVkTh(cvL,b1L,b0L,Pinf1L,Pinf0L,gammaL,Texp[i],Pexp[i]));
        vgTh.push_back(computeVkTh(cvG,b1G,b0L,Pinf1L,Pinf0L,gammaG,Texp[i],Pexp[i]));
        elTh.push_back(computeThInternalEnergy(cvL,b1L,qL,Pinf1L,Pinf0L,gammaL,Texp[i],Pexp[i]));
        egTh.push_back(computeThInternalEnergy(cvG,b1G,qG,Pinf1G,Pinf0G,gammaG,Texp[i],Pexp[i]));
    }
//    writePlotFile("res/Psat_th.txt",Texp,PsatTh);
    writePlotFile("res/hL_th.txt",Texp,hLth);
    writePlotFile("res/hG_th.txt",Texp,hGth);
    writePlotFile("res/Lv_th.txt",Texp,LvTh);
    writePlotFile("res/vL_th.txt",Texp,vlTh);
    writePlotFile("res/vG_th.txt",Texp,vgTh);
    writePlotFile("res/eL_th.txt",Texp,elTh);
    writePlotFile("res/eG_th.txt",Texp,egTh);

    return 0;
}
