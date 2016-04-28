#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <nlopt.hpp>
#include "efodo.h"
#include "line.h"
#include "global.h"
//#include <armadillo>


using namespace std;


double objFunc(const vector<double> &x, std::vector<double> &grad, void* objFunc_data)
{
    vector<double> x0(10);  std::fill(x0.begin(), x0.end(), 0.0);
    x0[x_]=x[0];x0[px_]=x[1];x0[dE_]=x[2];x0[vt_]=x[3];

    LINE * FODO = (LINE *)objFunc_data;
    for(unsigned j=0;j!=FODO->Ncell;j++) FODO->Cell[j].Pass(x0);

    return (x0[x_]-x[0])*(x0[x_]-x[0]) +(x0[px_]-x[1])*(x0[px_]-x[1]) +(x0[dE_]-x[2])*(x0[dE_]-x[2])+(x0[vt_]-x[3])*(x0[vt_]-x[3]);
}


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

    cout << "LastElemType = " << FODO.Cell[FODO.Ncell-1].Type << endl;
    cout << "beforeWrf = " << FODO.Cell[FODO.Ncell-1].RFcav.Wrf << endl;
    FODO.Cell[FODO.Ncell-1].RFcav.Wrf*= 1.0-(1.0E-5)/3.0;
    cout << "afterWrf = " << FODO.Cell[FODO.Ncell-1].RFcav.Wrf << endl;
    cout << "DesignLegnth = " << FODO.Length << endl;
    const double cSpeed=299792458.0;
    cout << "RFLength = " << 6.0*M_PI*BETA*cSpeed/FODO.Cell[FODO.Ncell-1].RFcav.Wrf << endl;

    nlopt::algorithm algo;
    algo=nlopt::LN_COBYLA;
    algo=nlopt::LN_SBPLX;
    nlopt::opt opt(algo, 4);
    vector<double> lb(4); std::fill(lb.begin(), lb.end(), -0.01);
    vector<double> ub(4); std::fill(ub.begin(), ub.end(), 0.01);
    opt.set_lower_bounds(lb);
    opt.set_upper_bounds(ub);
    opt.set_min_objective(objFunc, &FODO);
    opt.set_xtol_rel(1e-8);
    opt.set_ftol_abs(1e-25);


    double objResult=0;
    vector<double> obj_x;
    obj_x.push_back(0.0);obj_x.push_back(0.0);obj_x.push_back(x[dE_]);obj_x.push_back(x[vt_]);
    opt.optimize(obj_x, objResult);


    cout << "closedOrbit deltaE = " << obj_x[2] << ",  vt = " << obj_x[3] << ",  x = " << obj_x[0] << ",   px = " << obj_x[1] << endl;
    x[x_]=obj_x[0];x[px_]=obj_x[1];x[dE_]=obj_x[2];x[vt_]=obj_x[3];

    for(unsigned j=0;j!=FODO.Ncell;j++) FODO.Cell[j].Pass(x);
    cout << "after one turn deltaE = " << x[dE_] << ",  vt = " << x[vt_] << ",  x = " << x[x_] << ",   px = " << x[px_] << endl;

    stringstream(argv[2]) >>  x[dE_];
    x[x_]=obj_x[0];x[px_]=obj_x[1];x[dE_]+=obj_x[2];x[vt_]=obj_x[3];
    //x[x_]=0;x[px_]=0;x[dE_]+=obj_x[2];x[vt_]=0;
    cout << "Design deltaE = " << x[dE_] << endl;


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


    for(unsigned i=0;i<nTurn;i++){
        if (i%nSkip==0)
        {
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

    return 0;
}
