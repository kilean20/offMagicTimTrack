#ifndef EBENDPASS_H
#define EBENDPASS_H

#include <cmath>
#include <vector>
#include <iostream>
using namespace std;

#include "global.h"

//=============================================================================
//                       eBendSpinPass 2nd order method
//=============================================================================
inline void
eBendSpinPass(vector<double> &x, double L, double ANGLE, unsigned Nint){
    x[s_]=0.0; // initialize s. s=0 is the element entrance
    const double iRHO=ANGLE/L;
    double hs = (1.0 + x[x_]*iRHO );
    double dummy = x[dE_]*BETA2 +1.0 -log(hs)*BETA2;

    //spin precession vector and dummy variables
    vector<double> W(3); double dummy2, Cos, Sin;

    // initialize p_s : hard edge fringe field effect
    x[ps_] = hs*sqrt(dummy*dummy/BETA2-1.0/BETAGAMMA2-x[px_]*x[px_]-x[pz_]*x[pz_]);

    //time step
    double lambda = 0.5*L/x[ps_]/(double)Nint;

    // main loop
    for(unsigned i=0; i<Nint ;i++)
    {
        if(i==Nint-1)
        {
            lambda = 0.5*( (hs*hs)/x[ps_]*(L-x[s_])
                            +x[px_]*hs*hs*hs*(L-x[s_])*(L-x[s_])*iRHO/(x[ps_]*x[ps_]) );
        }
        if(i==0 || i==Nint-1)
        {
            //kick
            x[px_] += lambda*iRHO*( x[ps_]*x[ps_]/(hs*hs) -dummy)/hs;
            x[s_] += x[ps_]/(hs*hs) *lambda;
            x[vt_] += dummy *lambda;
        }
        //transverse drift
        x[x_] += x[px_] *lambda;
        x[z_] += x[pz_] *lambda;
        //spin kick
        hs = (1.0 + x[x_]*iRHO );
        dummy = 1.0 +BETA2*x[dE_] -BETA2*log(hs); //electric potential
        W[0] = -0.5*BETAGAMMA*EDM*dummy*iRHO/hs;
        //W[2] = -x[ps_]*iRHO/(hs*hs)*(1.0-dummy);
        W[2] = x[ps_]*iRHO/(hs*hs)*( offMagicDelta/magicGAMMA*(1.0+1.0/(1.0+offMagicDelta*magicBETA2))  -BETA2*(GAMMA-1.0)*(x[dE_]-log(hs))/(1.0+GAMMA*dummy) );
        dummy = BETA2*GAMMA*(MDM+1.0/(1.0+dummy*GAMMA));
        W[1] = -x[pz_]*dummy*iRHO/hs;
        //dummy=norm(W);
        dummy=sqrt(W[0]*W[0]+W[1]*W[1]+W[2]*W[2]);
        dummy2=sqrt(1.0-x[Sx_]*x[Sx_]-x[Ss_]*x[Ss_]);
        if(dummy>1.0e-150){
            for(vector<double>::iterator it = W.begin(); it != W.end(); ++it)  *it/=dummy;
            //W[0]/=dummy;W[1]/=dummy;W[2]/=dummy;;
            dummy*=lambda;
            Cos=cos(dummy);
            Sin=sin(dummy);
            dummy = 2.0*Sin*( x[Sx_]*(W[0]*W[0]-1.0)*Sin +dummy2*(W[2]*W[0]*Sin+W[1]*Cos) +x[Ss_]*(W[0]*W[1]*Sin-W[2]*Cos) );
            x[Ss_]+= 2.0*Sin*( x[Ss_]*(W[1]*W[1]-1.0)*Sin +dummy2*(W[2]*W[1]*Sin-W[0]*Cos) +x[Sx_]*(W[0]*W[1]*Sin+W[2]*Cos) );
            x[Sx_]+=dummy;
        }
        else{
            dummy=2.0*lambda*(W[1]*dummy2-W[2]*x[Ss_]);
            x[Ss_]+=2.0*lambda*(W[2]*x[Sx_]-W[0]*dummy2);
            x[Sx_]+=dummy;
        }
        //transverse drift
        x[x_] += x[px_] *lambda;
        x[z_] += x[pz_] *lambda;
        //kick
        hs = (1.0 + x[x_]*iRHO );
        dummy = 1.0 +BETA2*x[dE_] -BETA2*log(hs); //electric potential
        if(i>Nint-3)
        {
            x[px_] += lambda*iRHO*( x[ps_]*x[ps_]/(hs*hs) -dummy)/hs;
            x[s_] += x[ps_]/(hs*hs) *lambda;
            x[vt_] += dummy *lambda;
        }
        else
        {
            x[px_] += 2.0*lambda*iRHO*( x[ps_]*x[ps_]/(hs*hs) -dummy)/hs;
            x[s_] += x[ps_]/(hs*hs) *2.0*lambda;
            x[vt_] += dummy *2.0*lambda;
        }
    }
    //------------correction on extra or remaining legnth to the edge-----------------
    lambda = 0.5*(  (hs*hs)/x[ps_]*(L-x[s_])
                    + x[px_]*hs*hs*hs*(L-x[s_])*(L-x[s_])*iRHO/( x[ps_]*x[ps_] )  );
    //kick
    dummy = 1.0 +BETA2*x[dE_] -BETA2*log(hs); //electric potential
    x[px_] += lambda*iRHO*( x[ps_]*x[ps_]/(hs*hs) -dummy)/hs;
    //x[s_] += x[ps_]/(hs*hs) *lambda; //uncomment if you want to check the correction worked fine
    x[vt_] += dummy *lambda;
    //transverse drift
    x[x_] += x[px_] *lambda;
    x[z_] += x[pz_] *lambda;
    //spin kick
    hs = (1.0 + x[x_]*iRHO );
    dummy = 1.0 +BETA2*x[dE_] -BETA2*log(hs); //electric potential
    W[0] = -0.5*BETAGAMMA*EDM*dummy*iRHO/hs;
    //W[2] = -x[ps_]*iRHO/(hs*hs)*(1.0-dummy);
    W[2] = x[ps_]*iRHO/(hs*hs)*( offMagicDelta/magicGAMMA*(1.0+1.0/(1.0+offMagicDelta*magicBETA2))  -BETA2*(GAMMA-1.0)*(x[dE_]-log(hs))/(1.0+GAMMA*dummy) );
    dummy = BETA2*GAMMA*(MDM+1.0/(1.0+dummy*GAMMA));
    W[1] = -x[pz_]*dummy*iRHO/hs;
    //dummy=norm(W);
    dummy=sqrt(W[0]*W[0]+W[1]*W[1]+W[2]*W[2]);
    dummy2=sqrt(1.0-x[Sx_]*x[Sx_]-x[Ss_]*x[Ss_]);
    if(dummy>1.0e-150){
        for(vector<double>::iterator it = W.begin(); it != W.end(); ++it)  *it/=dummy;
        //W[0]/=dummy;W[1]/=dummy;W[2]/=dummy;;
        dummy*=lambda;
        Cos=cos(dummy);
        Sin=sin(dummy);
        dummy = 2.0*Sin*( x[Sx_]*(W[0]*W[0]-1.0)*Sin +dummy2*(W[2]*W[0]*Sin+W[1]*Cos) +x[Ss_]*(W[0]*W[1]*Sin-W[2]*Cos) );
        x[Ss_]+= 2.0*Sin*( x[Ss_]*(W[1]*W[1]-1.0)*Sin +dummy2*(W[2]*W[1]*Sin-W[0]*Cos) +x[Sx_]*(W[0]*W[1]*Sin+W[2]*Cos) );
        x[Sx_]+=dummy;
    }
    else{
        dummy=2.0*lambda*(W[1]*dummy2-W[2]*x[Ss_]);
        x[Ss_]+=2.0*lambda*(W[2]*x[Sx_]-W[1]*dummy2);
        x[Sx_]+=dummy;
    }
    //transverse drift
    x[x_] += x[px_] *lambda;
    x[z_] += x[pz_] *lambda;
    //kick
    hs = (1.0 + x[x_]*iRHO );
    dummy = 1.0 +BETA2*x[dE_] -BETA2*log(hs); //electric potential
    x[px_] += lambda*iRHO*( x[ps_]*x[ps_]/(hs*hs) -dummy)/hs;
    //x[s_] += x[ps_]/(hs*hs) *lambda;//uncomment if you want to check the correction worked fine
    x[vt_] += dummy * lambda;  // -L : deviation from the reference particle
}
//=============================================================================
//                 eBendSpinPass High order composition method
//=============================================================================
inline void
eBendSpinPass(vector<double> &x, double L, double ANGLE, unsigned Nint, unsigned Norder){
    // initialize step size
    vector<double> R;
    unsigned nR;
    switch (Norder) {
    case 2:
        R.push_back(1.0);
        nR = 1;
        break;
    case 4:
        for(unsigned i=0;i<3;i++)
        R.push_back(R4[i]);
        nR = 3;
        break;
    case 6:
        for(unsigned i=0;i<7;i++)
        R.push_back(R6[i]);
        nR = 7;
        break;
    default:
        for(unsigned i=0;i<3;i++)
        R.push_back(R4[i]);
        nR = 3;
        break;
    }

    x[s_]=0.0; // initialize s. s=0 is the element entrance
    const double iRHO=ANGLE/L;
    double hs = (1.0 + x[x_]*iRHO );
    double dummy = x[dE_]*BETA2 +1.0 -log(hs)*BETA2;

    //spin precession vector
    vector<double> W(3); double dummy2, Cos, Sin;

    // initialize p_s : hard edge fringe field effect
    x[ps_] = hs*sqrt(dummy*dummy/BETA2-1.0/BETAGAMMA2-x[px_]*x[px_]-x[pz_]*x[pz_]);

    //time step
    double lambda = 0.5*L/x[ps_]/(double)Nint;

    // main loop
    for(unsigned i=0; i<Nint ;i++)
    {
        if(i==Nint-1)
        {
            lambda = 0.5*(  (hs*hs)/x[ps_]*(L-x[s_])
                            + x[px_]*hs*hs*hs*(L-x[s_])*(L-x[s_])*iRHO/( x[ps_]*x[ps_] )  );
        }
        if(i==0 || i==Nint-1)
        {
            //kick
            x[px_] += R[0]*lambda*iRHO*( x[ps_]*x[ps_]/(hs*hs) -dummy)/hs;
            x[s_] += x[ps_]/(hs*hs) *R[0]*lambda;
            x[vt_] += dummy *R[0]*lambda;
        }
        //transverse drift
        x[x_] += x[px_] *R[0]*lambda;
        x[z_] += x[pz_] *R[0]*lambda;
        //spin kick
        hs = (1.0 + x[x_]*iRHO );
        dummy = 1.0 +BETA2*x[dE_] -BETA2*log(hs); //electric potential
        W[0] = -0.5*BETAGAMMA*EDM*dummy*iRHO/hs;
        //W[2] = -x[ps_]*iRHO/(hs*hs)*(1.0-dummy);
        W[2] = x[ps_]*iRHO/(hs*hs)*( offMagicDelta/magicGAMMA*(1.0+1.0/(1.0+offMagicDelta*magicBETA2))  -BETA2*(GAMMA-1.0)*(x[dE_]-log(hs))/(1.0+GAMMA*dummy) );
        dummy = BETA2*GAMMA*(MDM+1.0/(1.0+dummy*GAMMA));
        W[1] = -x[pz_]*dummy*iRHO/hs;
        //dummy=norm(W);
        dummy=sqrt(W[0]*W[0]+W[1]*W[1]+W[2]*W[2]);
        dummy2=sqrt(1.0-x[Sx_]*x[Sx_]-x[Ss_]*x[Ss_]);
        if(dummy>1.0e-150){
            for(vector<double>::iterator it = W.begin(); it != W.end(); ++it)  *it/=dummy;
            //W[0]/=dummy;W[1]/=dummy;W[2]/=dummy;;
            dummy*=R[0]*lambda;
            Cos=cos(dummy);
            Sin=sin(dummy);
            dummy = 2.0*Sin*( x[Sx_]*(W[0]*W[0]-1.0)*Sin +dummy2*(W[2]*W[0]*Sin+W[1]*Cos) +x[Ss_]*(W[0]*W[1]*Sin-W[2]*Cos) );
            x[Ss_]+= 2.0*Sin*( x[Ss_]*(W[1]*W[1]-1.0)*Sin +dummy2*(W[2]*W[1]*Sin-W[0]*Cos) +x[Sx_]*(W[0]*W[1]*Sin+W[2]*Cos) );
            x[Sx_]+=dummy;
        }
        else{
            dummy=2.0*R[0]*lambda*(W[1]*dummy2-W[2]*x[Ss_]);
            x[Ss_]+=2.0*R[0]*lambda*(W[2]*x[Sx_]-W[0]*dummy2);
            x[Sx_]+=dummy;
        }
        //transverse drift
        x[x_] += x[px_] *R[0]*lambda;
        x[z_] += x[pz_] *R[0]*lambda;
        for(unsigned r=1; r<nR; r++)
        {
            //kick
            hs = (1.0 + x[x_]*iRHO );
            dummy = 1.0 +BETA2*x[dE_] -BETA2*log(hs); //electric potential
            x[px_] += (R[r]+R[r-1])*lambda*iRHO*( x[ps_]*x[ps_]/(hs*hs) -dummy)/hs;
            x[s_] += x[ps_]/(hs*hs) *(R[r]+R[r-1])*lambda;
            x[vt_] += dummy * (R[r]+R[r-1])*lambda;
            //transverse drift
            x[x_] += x[px_] *R[r]*lambda;
            x[z_] += x[pz_] *R[r]*lambda;
            //spin kick
            hs = (1.0 + x[x_]*iRHO );
            dummy = 1.0 +BETA2*x[dE_] -BETA2*log(hs); //electric potential
            W[0] = -0.5*BETAGAMMA*EDM*dummy*iRHO/hs;
            //W[2] = -x[ps_]*iRHO/(hs*hs)*(1.0-dummy);
            W[2] = x[ps_]*iRHO/(hs*hs)*( offMagicDelta/magicGAMMA*(1.0+1.0/(1.0+offMagicDelta*magicBETA2))  -BETA2*(GAMMA-1.0)*(x[dE_]-log(hs))/(1.0+GAMMA*dummy) );
            dummy = BETA2*GAMMA*(MDM+1.0/(1.0+dummy*GAMMA));
            W[1] = -x[pz_]*dummy*iRHO/hs;
            //dummy=norm(W);
            dummy=sqrt(W[0]*W[0]+W[1]*W[1]+W[2]*W[2]);
            dummy2=sqrt(1.0-x[Sx_]*x[Sx_]-x[Ss_]*x[Ss_]);
            if(dummy>1.0e-150){
                for(vector<double>::iterator it = W.begin(); it != W.end(); ++it)  *it/=dummy;
                //W[0]/=dummy;W[1]/=dummy;W[2]/=dummy;;
                dummy*=R[r]*lambda;
                Cos=cos(dummy);
                Sin=sin(dummy);
                dummy = 2.0*Sin*( x[Sx_]*(W[0]*W[0]-1.0)*Sin +dummy2*(W[2]*W[0]*Sin+W[1]*Cos) +x[Ss_]*(W[0]*W[1]*Sin-W[2]*Cos) );
                x[Ss_]+= 2.0*Sin*( x[Ss_]*(W[1]*W[1]-1.0)*Sin +dummy2*(W[2]*W[1]*Sin-W[0]*Cos) +x[Sx_]*(W[0]*W[1]*Sin+W[2]*Cos) );
                x[Sx_]+=dummy;
            }
            else{
                dummy=2.0*R[r]*lambda*(W[1]*dummy2-W[2]*x[Ss_]);
                x[Ss_]+=2.0*R[r]*lambda*(W[2]*x[Sx_]-W[0]*dummy2);
                x[Sx_]+=dummy;
            }
            //transverse drift
            x[x_] += x[px_] *R[r]*lambda;
            x[z_] += x[pz_] *R[r]*lambda;
        }
        //kick
        hs = (1.0 + x[x_]*iRHO );
        dummy = 1.0 +BETA2*x[dE_] -BETA2*log(hs); //electric potential
        if(i>Nint-3)
        {
            x[px_] += R[nR-1]*lambda*iRHO*( x[ps_]*x[ps_]/(hs*hs) -dummy)/hs;
            x[s_] += x[ps_]/(hs*hs) *R[nR-1]*lambda;
            x[vt_] += dummy *R[nR-1]*lambda;
        }
        else
        {
            x[px_] += R[nR-1]*2.0*lambda*iRHO*( x[ps_]*x[ps_]/(hs*hs) -dummy)/hs;
            x[s_] += x[ps_]/(hs*hs) *R[nR-1]*2.0*lambda;
            x[vt_] += dummy *R[nR-1]*2.0*lambda;
        }
    }

    //------------correction on extra or remaining legnth to the edge-----------------
    lambda = 0.5*(  (hs*hs)/x[ps_]*(L-x[s_])
                    + x[px_]*hs*hs*hs*(L-x[s_])*(L-x[s_])*iRHO/( x[ps_]*x[ps_] )  );
    //kick
    dummy = 1.0 +BETA2*x[dE_] -BETA2*log(hs); //electric potential
    x[px_] += lambda*iRHO*( x[ps_]*x[ps_]/(hs*hs) -dummy)/hs;
    //x[s_] += x[ps_]/(hs*hs) *lambda; // uncomment when you like to check if correction worked fine
    x[vt_] += dummy *lambda;
    //transverse drift
    x[x_] += x[px_] *lambda;
    x[z_] += x[pz_] *lambda;
    //spin kick
    hs = (1.0 + x[x_]*iRHO );
    dummy = 1.0 +BETA2*x[dE_] -BETA2*log(hs); //electric potential
    W[0] = -0.5*BETAGAMMA*EDM*dummy*iRHO/hs;
    //W[2] = -x[ps_]*iRHO/(hs*hs)*(1.0-dummy);
    W[2] = x[ps_]*iRHO/(hs*hs)*( offMagicDelta/magicGAMMA*(1.0+1.0/(1.0+offMagicDelta*magicBETA2))  -BETA2*(GAMMA-1.0)*(x[dE_]-log(hs))/(1.0+GAMMA*dummy) );
    dummy = BETA2*GAMMA*(MDM+1.0/(1.0+dummy*GAMMA));
    W[1] = -x[pz_]*dummy*iRHO/hs;
    //dummy=norm(W);
    dummy=sqrt(W[0]*W[0]+W[1]*W[1]+W[2]*W[2]);
    dummy2=sqrt(1.0-x[Sx_]*x[Sx_]-x[Ss_]*x[Ss_]);
    if(dummy>1.0e-150){
        for(vector<double>::iterator it = W.begin(); it != W.end(); ++it)  *it/=dummy;
        //W[0]/=dummy;W[1]/=dummy;W[2]/=dummy;;
        dummy*=lambda;
        Cos=cos(dummy);
        Sin=sin(dummy);
        dummy = 2.0*Sin*( x[Sx_]*(W[0]*W[0]-1.0)*Sin +dummy2*(W[2]*W[0]*Sin+W[1]*Cos) +x[Ss_]*(W[0]*W[1]*Sin-W[2]*Cos) );
        x[Ss_]+= 2.0*Sin*( x[Ss_]*(W[1]*W[1]-1.0)*Sin +dummy2*(W[2]*W[1]*Sin-W[0]*Cos) +x[Sx_]*(W[0]*W[1]*Sin+W[2]*Cos) );
        x[Sx_]+=dummy;
    }
    else{
        dummy=2.0*lambda*(W[1]*dummy2-W[2]*x[Ss_]);
        x[Ss_]+=2.0*lambda*(W[2]*x[Sx_]-W[1]*dummy2);
        x[Sx_]+=dummy;
    }
    //transverse drift
    x[x_] += x[px_] *lambda;
    x[z_] += x[pz_] *lambda;
    //kick
    hs = (1.0 + x[x_]*iRHO );
    dummy = 1.0 +BETA2*x[dE_] -BETA2*log(hs); //electric potential
    x[px_] += lambda*iRHO*( x[ps_]*x[ps_]/(hs*hs) -dummy)/hs;
    //x[s_] += x[ps_]/(hs*hs) *lambda; // uncomment when you like to check if correction worked fine
    x[vt_] += dummy * lambda;
}
//=============================================================================
//                       eBendOrbitPass 2nd order method
//=============================================================================
inline void
eBendOrbitPass(vector<double> &x, double L, double ANGLE, unsigned Nint){
    x[s_]=0.0; // initialize s. s=0 is the element entrance
    const double iRHO=ANGLE/L;
    double hs = (1.0 + x[x_]*iRHO );
    double dummy = x[dE_]*BETA2 +1.0 -log(hs)*BETA2;

    // initialize p_s : hard edge fringe field effect
    x[ps_] = hs*sqrt(dummy*dummy/BETA2-1.0/BETAGAMMA2-x[px_]*x[px_]-x[pz_]*x[pz_]);

    //time step
    double lambda = 0.5*L/x[ps_]/(double)Nint;

    // main loop
    for(unsigned i=0; i<Nint ;i++)
    {
        if(i==Nint-1)
        {
            lambda = 0.5*( (hs*hs)/x[ps_]*(L-x[s_])
                            +x[px_]*hs*hs*hs*(L-x[s_])*(L-x[s_])*iRHO/(x[ps_]*x[ps_]) );
        }
        if(i==0 || i==Nint-1)
        {
            //kick
            x[px_] += lambda*iRHO*( x[ps_]*x[ps_]/(hs*hs) -dummy)/hs;
            x[s_] += x[ps_]/(hs*hs) *lambda;
            x[vt_] += dummy *lambda;
        }
        //transverse drift
        x[x_] += x[px_] *2.0*lambda;
        x[z_] += x[pz_] *2.0*lambda;
        //kick
        hs = (1.0 + x[x_]*iRHO );
        dummy = 1.0 +BETA2*x[dE_] -BETA2*log(hs); //electric potential
        if(i>Nint-3)
        {
            x[px_] += lambda*iRHO*( x[ps_]*x[ps_]/(hs*hs) -dummy)/hs;
            x[s_] += x[ps_]/(hs*hs) *lambda;
            x[vt_] += dummy *lambda;
        }
        else
        {
            x[px_] += 2.0*lambda*iRHO*( x[ps_]*x[ps_]/(hs*hs) -dummy)/hs;
            x[s_] += x[ps_]/(hs*hs) *2.0*lambda;
            x[vt_] += dummy *2.0*lambda;
        }
    }

    //------------correction on extra or remaining legnth to the edge-----------------
    lambda = 0.5*(  (hs*hs)/x[ps_]*(L-x[s_])
                    + x[px_]*hs*hs*hs*(L-x[s_])*(L-x[s_])*iRHO/( x[ps_]*x[ps_] )  );
    //kick
    dummy = 1.0 +BETA2*x[dE_] -BETA2*log(hs); //electric potential
    x[px_] += lambda*iRHO*( x[ps_]*x[ps_]/(hs*hs) -dummy)/hs;
    //x[s_] += x[ps_]/(hs*hs) *lambda; // uncomment when you like to check if correction worked fine
    x[vt_] += dummy *lambda;
    //transverse drift
    x[x_] += x[px_] *2.0*lambda;
    x[z_] += x[pz_] *2.0*lambda;
    //kick
    hs = (1.0 + x[x_]*iRHO );
    dummy = 1.0 +BETA2*x[dE_] -BETA2*log(hs); //electric potential
    x[px_] += lambda*iRHO*( x[ps_]*x[ps_]/(hs*hs) -dummy)/hs;
    //x[s_] += x[ps_]/(hs*hs) *lambda; // uncomment when you like to check if correction worked fine
    x[vt_] += dummy * lambda;
    //x[vt_] -= L; //deviation from the reference particle
}
//=============================================================================
//                eBendOrbitPass High order composition method
//=============================================================================
inline void
eBendOrbitPass(vector<double> &x, double L, double ANGLE, unsigned Nint, unsigned Norder){

    // initialize step size
    vector<double> R;
    unsigned nR;
    switch (Norder) {
    case 2:
        R.push_back(1.0);
        nR = 1;
        break;
    case 4:
        for(unsigned i=0;i<3;i++)
        R.push_back(R4[i]);
        nR = 3;
        break;
    case 6:
        for(unsigned i=0;i<7;i++)
        R.push_back(R6[i]);
        nR = 7;
        break;
    default:
        for(unsigned i=0;i<3;i++)
        R.push_back(R4[i]);
        nR = 3;
        break;
    }

    x[s_]=0.0; // initialize s. s=0 is the element entrance
    const double iRHO=ANGLE/L;
    double hs = (1.0 + x[x_]*iRHO );
    double dummy = x[dE_]*BETA2 +1.0 -log(hs)*BETA2;

    // initialize p_s : hard edge fringe field effect
    x[ps_] = hs*sqrt(dummy*dummy/BETA2-1.0/BETAGAMMA2-x[px_]*x[px_]-x[pz_]*x[pz_]);

    //time step
    double lambda = 0.5*L/x[ps_]/(double)Nint;

    // main loop
    for(unsigned i=0; i<Nint ;i++)
    {
        if(i==Nint-1)
        {
            lambda = 0.5*(  (hs*hs)/x[ps_]*(L-x[s_])
                            + x[px_]*hs*hs*hs*(L-x[s_])*(L-x[s_])*iRHO/( x[ps_]*x[ps_] )  );
        }
        if(i==0 || i==Nint-1)
        {
            //kick
            x[px_] += R[0]*lambda*iRHO*( x[ps_]*x[ps_]/(hs*hs) -dummy)/hs;
            x[s_] += x[ps_]/(hs*hs) *R[0]*lambda;
            x[vt_] += dummy *R[0]*lambda;
        }
        //transverse drift
        x[x_] += x[px_] *R[0]*2.0*lambda;
        x[z_] += x[pz_] *R[0]*2.0*lambda;
        for(unsigned r=1; r<nR; r++)
        {
            //kick
            hs = (1.0 + x[x_]*iRHO );
            dummy = 1.0 +BETA2*x[dE_] -BETA2*log(hs); //electric potential
            x[px_] += (R[r]+R[r-1])*lambda*iRHO*( x[ps_]*x[ps_]/(hs*hs) -dummy)/hs;
            x[s_] += x[ps_]/(hs*hs) *(R[r]+R[r-1])*lambda;
            x[vt_] += dummy * (R[r]+R[r-1])*lambda;
            //transverse drift
            x[x_] += x[px_] *R[r]*2.0*lambda;
            x[z_] += x[pz_] *R[r]*2.0*lambda;
        }
        //kick
        hs = (1.0 + x[x_]*iRHO );
        dummy = 1.0 +BETA2*x[dE_] -BETA2*log(hs); //electric potential
        if(i>Nint-3)
        {
            x[px_] += R[nR-1]*lambda*iRHO*( x[ps_]*x[ps_]/(hs*hs) -dummy)/hs;
            x[s_] += x[ps_]/(hs*hs) *R[nR-1]*lambda;
            x[vt_] += dummy *R[nR-1]*lambda;
        }
        else
        {
            x[px_] += R[nR-1]*2.0*lambda*iRHO*( x[ps_]*x[ps_]/(hs*hs) -dummy)/hs;
            x[s_] += x[ps_]/(hs*hs) *R[nR-1]*2.0*lambda;
            x[vt_] += dummy *R[nR-1]*2.0*lambda;
        }
    }

    //------------correction on extra or remaining legnth to the edge-----------------
    lambda = 0.5*(  (hs*hs)/x[ps_]*(L-x[s_])
                    + x[px_]*hs*hs*hs*(L-x[s_])*(L-x[s_])*iRHO/( x[ps_]*x[ps_] )  );
    //kick
    dummy = 1.0 +BETA2*x[dE_] -BETA2*log(hs); //electric potential
    x[px_] += lambda*iRHO*( x[ps_]*x[ps_]/(hs*hs) -dummy)/hs;
    //x[s_] += x[ps_]/(hs*hs) *lambda; // uncomment when you like to check if correction worked fine
    x[vt_] += dummy *lambda;
    //transverse drift
    x[x_] += x[px_] *2.0*lambda;
    x[z_] += x[pz_] *2.0*lambda;
    //kick
    hs = (1.0 + x[x_]*iRHO );
    dummy = 1.0 +BETA2*x[dE_] -BETA2*log(hs); //electric potential
    x[px_] += lambda*iRHO*( x[ps_]*x[ps_]/(hs*hs) -dummy)/hs;
    //x[s_] += x[ps_]/(hs*hs) *lambda; // uncomment when you like to check if correction worked fine
    x[vt_] += dummy * lambda;
}

#endif // EBENDPASS_H

