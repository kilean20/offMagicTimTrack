#ifndef EQUADPASS_H
#define EQUADPASS_H

#include <cmath>
#include <vector>
#include <iostream>
using namespace std;

#include "global.h"

//=============================================================================
//                       eQuadSpinPass 2nd order method
//=============================================================================
inline void
eQuadSpinPass(vector<double> &x, double L, double K1, unsigned Nint){
    const double betaGammaK1 = BETAGAMMA*K1;
    double dummy = x[dE_]*BETA2 + 1.0 - 0.5*K1*(x[x_]*x[x_]-x[z_]*x[z_])*BETA2;

    //spin precession vector and dummy variables
    vector<double> W(3); double dummy2, Cos, Sin;

    // initialize p_s : hard edge fringe field effect
    x[ps_]=sqrt( dummy*dummy/BETA2-1.0/BETAGAMMA2 - x[px_]*x[px_] - x[pz_]*x[pz_] );

    //timeStep
    double lambda = 0.5*L/x[ps_]/(double)Nint;

    // main loop
    for(unsigned i=0; i<Nint ;i++)
    {
        if(i==0)
        {
            //kick
            x[px_] -= dummy*K1*x[x_] *lambda;
            x[pz_] += dummy*K1*x[z_] *lambda;
            x[vt_] += dummy *lambda;
        }
        //transverse drift
        x[x_] += x[px_] *lambda;
        x[z_] += x[pz_] *lambda;
        //spin kick
        dummy = x[dE_]*BETA2 +1.0 -0.5*K1*(x[x_]*x[x_]-x[z_]*x[z_])*BETA2;
        dummy2 = BETA*(MDM+1.0/(1.0+GAMMA*dummy));
        W[0] = betaGammaK1*( -0.5*x[x_]*EDM*dummy + x[ps_]*x[z_]*dummy2  );
        W[1] = -betaGammaK1*( x[pz_]*x[x_] + x[px_]*x[z_]  )*dummy2;
        W[2] = betaGammaK1*( 0.5*x[z_]*EDM*dummy + x[ps_]*x[x_]*dummy2  );
        //dummy=norm(W);
        dummy=sqrt(W[0]*W[0]+W[1]*W[1]+W[2]*W[2]);
        dummy2=sqrt(1.0-x[Sx_]*x[Sx_]-x[Ss_]*x[Ss_]);
        if(dummy>1.0e-150){
            W[0]/=dummy;W[1]/=dummy;W[2]/=dummy;
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
        dummy = x[dE_]*BETA2 +1.0 - 0.5*K1*BETA2*(x[x_]*x[x_]-x[z_]*x[z_]);
        if(i==Nint-1)
        {
            x[px_] -= dummy*K1*x[x_] *lambda;
            x[pz_] += dummy*K1*x[z_] *lambda;
            x[vt_] += dummy *lambda;
        }else
        {
            x[px_] -= dummy*K1*x[x_] *2.0*lambda;
            x[pz_] += dummy*K1*x[z_] *2.0*lambda;
            x[vt_] += dummy *2.0*lambda;
        }
    }
    //deviation from the reference particle
    //x[vt_] -= L;
}
//=============================================================================
//                 eQuadSpinPass High order composition method
//=============================================================================
inline void
eQuadSpinPass(vector<double> &x, double L, double K1, unsigned Nint, unsigned Norder){
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

    const double betaGammaK1 = BETAGAMMA*K1;
    double dummy = x[dE_]*BETA2 + 1.0 - 0.5*K1*(x[x_]*x[x_]-x[z_]*x[z_])*BETA2;

    //spin precession vector and dummy variables
    vector<double> W(3); double dummy2, Cos, Sin;

    // initialize p_s : hard edge fringe field effect
    x[ps_]=sqrt( dummy*dummy/BETA2-1.0/BETAGAMMA2 - x[px_]*x[px_] - x[pz_]*x[pz_] );

    //timeStep
    double lambda = 0.5*L/x[ps_]/(double)Nint;

    // main loop
    for(unsigned i=0; i<Nint ;i++)
    {
        if(i==0)
        {
            //kick
            x[px_] -= dummy*K1*x[x_] *R[0]*lambda;
            x[pz_] += dummy*K1*x[z_] *R[0]*lambda;
            x[vt_] += dummy *R[0]*lambda;
        }
        //drift
        x[x_] += x[px_] *R[0]*lambda;
        x[z_] += x[pz_] *R[0]*lambda;
        //spin kick
        dummy = x[dE_]*BETA2 +1.0 -0.5*K1*(x[x_]*x[x_]-x[z_]*x[z_])*BETA2;
        dummy2 = BETA*(MDM+1.0/(1.0+GAMMA*dummy));
        W[0] = betaGammaK1*( -0.5*x[x_]*EDM*dummy + x[ps_]*x[z_]*dummy2  );
        W[1] = -betaGammaK1*( x[pz_]*x[x_] + x[px_]*x[z_]  )*dummy2;
        W[2] = betaGammaK1*( 0.5*x[z_]*EDM*dummy + x[ps_]*x[x_]*dummy2  );
        //dummy=norm(W);
        dummy=sqrt(W[0]*W[0]+W[1]*W[1]+W[2]*W[2]);
        dummy2=sqrt(1.0-x[Sx_]*x[Sx_]-x[Ss_]*x[Ss_]);
        if(dummy>1.0e-150){
            W[0]/=dummy;W[1]/=dummy;W[2]/=dummy;;
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
        //drift
        x[x_] += x[px_] *R[0]*lambda;
        x[z_] += x[pz_] *R[0]*lambda;
        for(unsigned r=1; r<nR; r++)
        {
            //kick
            dummy = x[dE_]*BETA2 +1.0 - 0.5*K1*BETA2*(x[x_]*x[x_]-x[z_]*x[z_]);
            x[px_] -= dummy*K1*x[x_] *(R[r]+R[r-1])*lambda;
            x[pz_] += dummy*K1*x[z_] *(R[r]+R[r-1])*lambda;
            x[vt_] += dummy *(R[r]+R[r-1])*lambda;
            //drift
            x[x_] += x[px_] *R[r]*lambda;
            x[z_] += x[pz_] *R[r]*lambda;
            //spin kick
            dummy = x[dE_]*BETA2 +1.0 -0.5*K1*(x[x_]*x[x_]-x[z_]*x[z_])*BETA2;
            dummy2 = BETA*(MDM+1.0/(1.0+GAMMA*dummy));
            W[0] = betaGammaK1*( -0.5*x[x_]*EDM*dummy + x[ps_]*x[z_]*dummy2  );
            W[1] = -betaGammaK1*( x[pz_]*x[x_] + x[px_]*x[z_]  )*dummy2;
            W[2] = betaGammaK1*( 0.5*x[z_]*EDM*dummy + x[ps_]*x[x_]*dummy2  );
            //dummy=norm(W);
            dummy=sqrt(W[0]*W[0]+W[1]*W[1]+W[2]*W[2]);
            dummy2=sqrt(1.0-x[Sx_]*x[Sx_]-x[Ss_]*x[Ss_]);
            if(dummy>1.0e-150){
                W[0]/=dummy;W[1]/=dummy;W[2]/=dummy;;
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
            //drift
            x[x_] += x[px_] *R[r]*lambda;
            x[z_] += x[pz_] *R[r]*lambda;
        }
        //kick
        dummy = x[dE_]*BETA2 +1.0 - 0.5*K1*BETA2*(x[x_]*x[x_]-x[z_]*x[z_]);
        if(i==Nint-1)
        {
            x[px_] -= dummy*K1*x[x_] *R[nR-1]*lambda;
            x[pz_] += dummy*K1*x[z_] *R[nR-1]*lambda;
            x[vt_] += dummy *R[nR-1]*lambda;
        }else
        {
            x[px_] -= dummy*K1*x[x_] *R[nR-1]*2.0*lambda;
            x[pz_] += dummy*K1*x[z_] *R[nR-1]*2.0*lambda;
            x[vt_] += dummy *R[nR-1]*2.0*lambda;
        }
    }
    //deviation from the reference particle
    //x[vt_] -= L;
}
//=============================================================================
//                       eQuadOrbitPass 2nd order method
//=============================================================================
inline void
eQuadOrbitPass(vector<double> &x, double L, double K1, unsigned Nint){
    double dummy = x[dE_]*BETA2 + 1.0 - 0.5*K1*(x[x_]*x[x_]-x[z_]*x[z_])*BETA2;

    // initialize p_s : hard edge fringe field effect
    x[ps_]=sqrt( dummy*dummy/BETA2-1.0/BETAGAMMA2 - x[px_]*x[px_] - x[pz_]*x[pz_] );

    //timeStep
    double lambda = 0.5*L/x[ps_]/(double)Nint;

    // main loop
    for(unsigned i=0; i<Nint ;i++)
    {
        if(i==0)
        {
            //kick
            x[px_] -= dummy*K1*x[x_] *lambda;
            x[pz_] += dummy*K1*x[z_] *lambda;
            x[vt_] += dummy *lambda;
        }
        //transverse drift
        x[x_] += x[px_] *2.0*lambda;
        x[z_] += x[pz_] *2.0*lambda;
        //kick
        dummy = x[dE_]*BETA2 +1.0 - 0.5*K1*BETA2*(x[x_]*x[x_]-x[z_]*x[z_]);
        if(i==Nint-1)
        {
            x[px_] -= dummy*K1*x[x_] *lambda;
            x[pz_] += dummy*K1*x[z_] *lambda;
            x[vt_] += dummy *lambda;
        }else
        {
            x[px_] -= dummy*K1*x[x_] *2.0*lambda;
            x[pz_] += dummy*K1*x[z_] *2.0*lambda;
            x[vt_] += dummy *2.0*lambda;
        }
    }
    //deviation from the reference particle
    //x[vt_] -= L;
}
//=============================================================================
//                eQuadOrbitPass High order composition method
//=============================================================================
inline void
eQuadOrbitPass(vector<double> &x, double L, double K1, unsigned Nint, unsigned Norder){
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

    double dummy = x[dE_]*BETA2 + 1.0 - 0.5*K1*(x[x_]*x[x_]-x[z_]*x[z_])*BETA2;

    // initialize p_s : hard edge fringe field effect
    x[ps_]=sqrt( dummy*dummy/BETA2-1.0/BETAGAMMA2 - x[px_]*x[px_] - x[pz_]*x[pz_] );

    //timeStep
    double lambda = 0.5*L/x[ps_]/(double)Nint;

    // main loop
    for(unsigned i=0; i<Nint ;i++)
    {
        if(i==0)
        {
            //kick
            x[px_] -= dummy*K1*x[x_] *R[0]*lambda;
            x[pz_] += dummy*K1*x[z_] *R[0]*lambda;
            x[vt_] += dummy *R[0]*lambda;
        }
        //drift
        x[x_] += x[px_] *R[0]*2.0*lambda;
        x[z_] += x[pz_] *R[0]*2.0*lambda;
        for(unsigned r=1; r<nR; r++)
        {
            //kick
            dummy = x[dE_]*BETA2 +1.0 - 0.5*K1*BETA2*(x[x_]*x[x_]-x[z_]*x[z_]);
            x[px_] -= dummy*K1*x[x_] *(R[r]+R[r-1])*lambda;
            x[pz_] += dummy*K1*x[z_] *(R[r]+R[r-1])*lambda;
            x[vt_] += dummy *(R[r]+R[r-1])*lambda;
            //drift
            x[x_] += x[px_] *R[r]*2.0*lambda;
            x[z_] += x[pz_] *R[r]*2.0*lambda;
        }
        //kick
        dummy = x[dE_]*BETA2 +1.0 - 0.5*K1*BETA2*(x[x_]*x[x_]-x[z_]*x[z_]);
        if(i==Nint-1)
        {
            x[px_] -= dummy*K1*x[x_] *R[nR-1]*lambda;
            x[pz_] += dummy*K1*x[z_] *R[nR-1]*lambda;
            x[vt_] += dummy *R[nR-1]*lambda;
        }else
        {
            x[px_] -= dummy*K1*x[x_] *R[nR-1]*2.0*lambda;
            x[pz_] += dummy*K1*x[z_] *R[nR-1]*2.0*lambda;
            x[vt_] += dummy *R[nR-1]*2.0*lambda;
        }
    }
    //deviation from the reference particle
    //x[vt_] -= L;
}

#endif  //EQUADPASS

