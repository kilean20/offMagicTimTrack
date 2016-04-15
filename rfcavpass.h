#ifndef RFCAVPASS_H
#define RFCAVPASS_H

#include <cmath>
#include <vector>

using namespace std;

#include "global.h"


inline void
RFcavPass(vector<double> &x, double VRF, double WRF, double PHASE){
    const double cSpeed=299792458.0;
    const double t = x[vt_]/(BETA*cSpeed);
    double delta_dE = -VRF*1.0E-6/(BETA2*ENERGY)*sin( WRF*t + PHASE);
    x[vt_]-= 6.0*M_PI*BETA*cSpeed/WRF;
    // sign corresponds to the negative charge
    x[dE_] += delta_dE;
}

#endif // RFCAVPASS_H

