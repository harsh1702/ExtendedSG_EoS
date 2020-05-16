#include "tools.h"

using namespace std;

// **************************************************

void readExpData(string file, vector<double> &Texp, vector<double> &Pexp, vector<double> &vGexp, vector<double> &vLexp, vector<double> &hGexp, vector<double> &hLexp, vector<double> &LvExp, vector<double> &eLexp, vector<double> &eGexp, std::vector<double> &sLexp, std::vector<double> &sGexp)
{
    // Purpose : read the experimental data with 6 columns and an ignored header
    //hGexp and hLexp are only required to verify the final results
    ifstream strmIn(file.c_str());
    string line;
    double dat1, dat2, dat3, dat4, dat5, dat6, dat7, dat8, dat9, dat10, dat11;
    if (strmIn) {
        while (getline(strmIn,line)) {
            strmIn >> dat1 >> dat2 >> dat3 >> dat4 >> dat5 >> dat6 >> dat7>>dat8>>dat9>>dat10>>dat11;
            Texp.push_back(dat1);
            Pexp.push_back(dat2);
            vGexp.push_back(dat3);
            vLexp.push_back(dat4);
            hGexp.push_back(dat5);
            hLexp.push_back(dat6);
            LvExp.push_back(dat7);
            eGexp.push_back(dat8);
            eLexp.push_back(dat9);
            sGexp.push_back(dat10);
            sLexp.push_back(dat11);
        }
    }
    else {
        cout << "Error : reading experimental data file " << file << "\n"; exit(0);
    }
}

// **************************************************

void readFile(string const &file, vector<double> &tab_x, vector<double> &tab_y)
{
    // Purpose : read experimental data file with two columns and a header (ignored)
    ifstream streamFile(file);
    string line;
    double data1(0.), data2(0.);

    if (streamFile) {
        while (getline(streamFile,line)) {
            streamFile >> data1 >> data2;
            tab_x.push_back(data1);
            tab_y.push_back(data2);
            cout << data1 << " " << data2 << endl;
        }
    }
    else {
        cout << "Error : reading experimental data file " << file << ".txt\n"; exit(0);
    }
}

// **************************************************

void writePlotFile(string file, vector<double> &x, vector<double> &y)
{
    // Purpose : write file with two columns
    ofstream strm(file.c_str());
    for (unsigned int i = 0; i < x.size(); i++) {
        strm << x[i] << " " << y[i] << endl;
    }
}

// **************************************************

// double residual(std::vector<double> &tabKnow, std::vector<double> &tabEstimated)
// {
//     // Purpose : compute the residual of an Ordinary Lest Squares (OLS) procedure
//     // Linear fn y=f(x) - vector are y values
//     double scr(0.), sct(0.);
//     if (tabKnow.size() == tabEstimated.size()) {
//         for (unsigned int i = 0; i < tabKnow.size(); i++) {
//         scr += (tabKnow[i]-tabEstimated[i])*(tabKnow[i]-tabEstimated[i]);
//         sct += (tabKnow[i]-meanValue(tabKnow))*(tabKnow[i]-meanValue(tabKnow));
//         }
//     }
//     else {
//         cout << "Error : residual can not be computed due to different vector size\n";
//     }
//     return (1.-(scr/sct));
// }

// **************************************************

void writeResults(double cpL, double qL, double pinfL, double cvL, double gammaL, double qPrimL, double bL, double cpG, double qG, double pinfG, double cvG, double gammaG, double qPrimG, double bG)
{
    ofstream strmRes("res/res.txt");
    strmRes << "* Results *\n";
    strmRes << "--- Liquid (L) ---\n";
    strmRes << "cpL    (J.kg-1.K-1)  : " << cpL << endl;
    strmRes << "qL     (J.kg-1)      : " << qL << endl;
    strmRes << "pinfL  (Pa)          : " << pinfL << endl;
    strmRes << "cvL    (J.kg-1.K-1)  : " << cvL << endl;
    strmRes << "gammaL (-)           : " << gammaL << endl;
    strmRes << "q'L    (J.kg-1)      : " << qPrimL << endl;
    strmRes << "bL     (m3/kg)       : " << bL << endl;
    strmRes << "\n";

    strmRes << "--- Gas (G) ---\n";
    strmRes << "cpG    (J.kg-1.K-1)  : " << cpG << endl;
    strmRes << "qG     (J.kg-1)      : " << qG << endl;
    strmRes << "pinfG  (Pa)          : " << pinfG << endl;
    strmRes << "cvG    (J.kg-1.K-1)  : " << cvG << endl;
    strmRes << "gammaG (-)           : " << gammaG << endl;
    strmRes << "q'G    (J.kg-1)      : " << qPrimG << endl;
    strmRes << "bG     (m3/kg)       : " << bG << endl;
}
