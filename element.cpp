#include "element.h"
#include "driftpass.h"
#include "ebendpass.h"
#include "equadpass.h"
#include "rfcavpass.h"
#include "global.h"

//#include <armadillo>

//=========================================================================
//
//                            Element classes
//
//=========================================================================
//-------------------------------constructor0------------------------------
ELEMENT::ELEMENT():Type(0),L(0){Norder=2; FlagSpinTrack=1;}
//-------------------------------constructor1------------------------------
ELEMENT::ELEMENT(unsigned type, double l):Type(type),L(l){
    Norder=2;   FlagSpinTrack=1;
    switch (type){
    case DRIFT_:
        DRIFT d;
        Drift=d;
        break;
    default:
        break;
    }
}
//------------------------------constructor2--------------------------------
ELEMENT::ELEMENT(unsigned type, double l,double param):Type(type),L(l){
    Norder=2;   FlagSpinTrack=1;
    switch (type){
    case eBEND_:
        eBEND b;  b.Angle=param;
        eBend=b;
        Nint = ceil(param/0.02); //default Nint
        break;
    case eQUAD_:
        eQUAD q;  q.K1=param;
        eQuad=q;
        Nint = (L*L*abs(param)/0.02); //default Nint
        break;
    default:
        break;
    }
}
//------------------------------constructor3--------------------------------
ELEMENT::ELEMENT(unsigned type, double v,double f,double phi):Type(type),L(0){
    Norder=2;   FlagSpinTrack=1;
    switch (type){
    case RFcav_:
        RFCAV rf;   rf.Vrf=v;   rf.Wrf=2*M_PI*f; rf.Phase=phi;
        RFcav=rf;
    default:
        break;
    }
}
//=================================SetElem=================================
//---------------------------------SetElem1--------------------------------
void
ELEMENT::SetElem(unsigned type, double l){
    Type=type; L=l;
    switch (type){
    case DRIFT_:
        DRIFT d;
        Drift=d;
        break;
    default:
        break;
    }
}
//---------------------------------SetElem2--------------------------------
void
ELEMENT::SetElem(unsigned type, double l, double param){
    Type=type; L=l;
    switch (type){
    case eBEND_:
        eBEND b;  b.Angle=param;
        eBend=b;
        Nint = ceil(param/0.02); //default Nint
        break;
    case eQUAD_:
        eQUAD q;  q.K1=param;
        eQuad=q;
        Nint = (L*L*abs(param)/0.02); //default Nint
        break;
    default:
        break;
    }
}
//---------------------------------SetElem3--------------------------------
void
ELEMENT::SetElem(unsigned type, double v,double f,double phi){
    Type=type;
    switch (type){
    case RFcav_:
        RFCAV rf;   rf.Vrf=v;   rf.Wrf=2*M_PI*f; rf.Phase=phi;
        RFcav=rf;
        L=0;
    default:
        break;
    }
}
//=========================================================================
//                              Pass method
//=========================================================================
void
ELEMENT::Pass (vector<double> &x)
{
    switch (Type) {
    case DRIFT_:
        DriftPass(x, L);
        break;
    case eBEND_:
        if (FlagSpinTrack)
            if (Norder==2)
                eBendSpinPass(x, L, eBend.Angle, Nint);
            else
                eBendSpinPass(x, L, eBend.Angle, Nint, Norder);
        else
            if (Norder==2)
                eBendOrbitPass(x, L, eBend.Angle, Nint);
            else
                eBendOrbitPass(x, L, eBend.Angle, Nint, Norder);
        break;
    case eQUAD_:
        if (FlagSpinTrack)
            if (Norder==2)
                eQuadSpinPass(x, L, eQuad.K1, Nint);
            else
                eQuadSpinPass(x, L, eQuad.K1, Nint, Norder);
        else
            if (Norder==2)
                eQuadOrbitPass(x, L, eQuad.K1, Nint);
            else
                eQuadOrbitPass(x, L, eQuad.K1, Nint, Norder);
        break;
    case RFcav_:
        RFcavPass(x, RFcav.Vrf, RFcav.Wrf, RFcav.Phase);
        break;
    default:
        break;
    }
}
