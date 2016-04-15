#ifndef ELEMENT_H
#define ELEMENT_H
#include <vector>
#include <cmath>
#include "global.h"

using namespace std;
//=========================================================================
//                            electric elements
//=========================================================================
//--------------------------------DRIFT------------------------------------
struct DRIFT
{
};
//--------------------------------eBEND------------------------------------
struct eBEND
{
    double Angle;
};
//--------------------------------eQUAD------------------------------------
struct eQUAD
{
    double K1;
};
//=========================================================================
//                             RF elements
//=========================================================================
//--------------------------------RFCAV------------------------------------
struct RFCAV
{
    double Vrf, Wrf, Phase;
};

//=========================================================================
//
//      Base class of elements - static polymorphism using union
//
//=========================================================================
class ELEMENT
{
public:
    union {
        DRIFT Drift;
        eBEND eBend;
        eQUAD eQuad;
        RFCAV RFcav;
    };
    unsigned Type;


    ELEMENT();
    ELEMENT(unsigned Type, double l);
    ELEMENT(unsigned Type, double l, double param);
    ELEMENT(unsigned Type, double param1, double param2, double param3);

    void SetElem(unsigned Type, double l);
    void SetElem(unsigned Type, double l, double param);
    void SetElem(unsigned Type, double param1, double param2, double param3);

    unsigned Nint, Norder, FlagSpinTrack;
    double L, S;
    void Pass (vector<double> &x);
};


#endif
