#include "element.h"
#include "line.h"
#include <cmath>

using namespace std;

void efodo(LINE & FODO)
{
    double ang = 2*M_PI/12/2;
    double quadK = 5.0;

    double lBend = 0.4;
    double lQuad = 0.1;
    double lDrift = 0.2;

  ELEMENT temp;
  for(int i=0;i<12;i++){
      temp.SetElem(DRIFT_,0.5*lDrift);
      FODO.Append(temp);
      temp.SetElem(eBEND_,lBend, ang);temp.Nint=12;temp.FlagSpinTrack=1;
      FODO.Append(temp);
      temp.SetElem(DRIFT_,lDrift);
      FODO.Append(temp);
      temp.SetElem(eQUAD_,lQuad, quadK);temp.Nint=3;temp.FlagSpinTrack=1;
      FODO.Append(temp);
      temp.SetElem(DRIFT_,lDrift);
      FODO.Append(temp);
      temp.SetElem(eBEND_,lBend, ang);temp.Nint=12;temp.FlagSpinTrack=1;
      FODO.Append(temp);
      temp.SetElem(DRIFT_,lDrift);
      FODO.Append(temp);
      temp.SetElem(eQUAD_,lQuad, -quadK);temp.Nint=3;temp.FlagSpinTrack=1;
      FODO.Append(temp);
      temp.SetElem(DRIFT_,0.5*lDrift);
      FODO.Append(temp);
  }
  const double cSpeed=299792458.0;
  const double revFreq = (BETA*cSpeed)/FODO.Length;

  temp.SetElem(RFcav_,1000.0, 3.0*revFreq, 0.0);
  FODO.Append(temp);
}
