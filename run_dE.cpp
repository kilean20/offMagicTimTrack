#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>

#include "efodo.h"
#include "line.h"
#include "global.h"
//#include <armadillo>


using namespace std;


int main(int argc, char *argv[])
{
    cout.precision(15);
    //Define the ring
    LINE FODO;    efodo(FODO);

    unsigned nTurn;
    stringstream(argv[3]) >> nTurn;

    unsigned nSkip;
    stringstream(argv[4]) >> nSkip;

    vector<double> x(10); std::fill(x.begin(), x.end(), 0.0);
    stringstream(argv[2]) >>  x[dE_];

    cout << "deltaE = " << x[dE_] << endl;

    ofstream xData,zData,sData,spinData;
    stringstream ss;
    ss<<"xData_dE"<<stringstream(argv[1]).str()<<".dat";
    xData.open(ss.str().c_str());

    ss.clear();ss.str("");
    ss<<"zData_dE"<<stringstream(argv[1]).str()<<".dat";
    zData.open(ss.str().c_str());

    ss.clear();ss.str("");
    ss<<"sData_dE"<<stringstream(argv[1]).str()<<".dat";
    sData.open(ss.str().c_str());

    ss.clear();ss.str("");
    ss<<"spinData_dE"<<stringstream(argv[1]).str()<<".dat";
    spinData.open(ss.str().c_str());

//    mat xData(nTurn/nSkip+1,2);xData.zeros();
//    mat zData(nTurn/nSkip+1,2);zData.zeros();
//    mat sData(nTurn/nSkip+1,2);sData.zeros();
//    mat spinData(nTurn/nSkip+1,2);spinData.zeros();


    for(unsigned i=0;i<nTurn;i++){
        if (i%nSkip==0)
        {
//            xData(i/nSkip,0)=x[x_];
//            xData(i/nSkip,1)=x[px_];
//            zData(i/nSkip,0)=x[z_];
//            zData(i/nSkip,1)=x[pz_];
//            sData(i/nSkip,0)=x[vt_];
//            sData(i/nSkip,1)=x[dE_];
//            spinData(i/nSkip,0)=x[Sx_];
//            spinData(i/nSkip,1)=x[Ss_];
            xData << x[x_] << " " << x[px_] << endl;
            zData << x[z_] << " " << x[pz_] << endl;
            sData << x[vt_] << " " << x[dE_] << endl;
            spinData << x[Sx_] << " " << x[Ss_] << endl;
        }
        for(unsigned j=0;j!=FODO.Ncell;j++) FODO.Cell[j].Pass(x);
    }
    zData << x[z_] << " " << x[pz_] << endl;
    sData << x[vt_] << " " << x[dE_] << endl;
    spinData << x[Sx_] << " " << x[Ss_] << endl;

//    stringstream ss;
//    ss<<"xData_dE"<<stringstream(argv[1])<<".dat";
//    xData.save(ss.str().c_str(),raw_ascii);
//    ss.clear();ss.str("");
//    ss<<"zData_dE"<<stringstream(argv[1])<<".dat";
//    zData.save(ss.str().c_str(),raw_ascii);
//    ss.clear();ss.str("");
//    ss<<"sData_dE"<<stringstream(argv[1])<<".dat";
//    sData.save(ss.str().c_str(),raw_ascii);
//    ss.clear();ss.str("");
//    ss<<"spinData_dE"<<stringstream(argv[1])<<".dat";
//    spinData.save(ss.str().c_str(),raw_ascii);
    return 0;
}
