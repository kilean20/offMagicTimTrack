#ifndef GLOBAL_H
#define GLOBAL_H

#include <cmath>
//=========================================================================
//                           Global  constants 
//=========================================================================
//const double PI      = 3.14159265;
//const double light_speed =  299792458;
//--------------------------spin-parameters--------------------------------
const double MDM = 0.00115965218;  // anormalous magnetic momentum
//const double EDM = 0.0;  // electric momentum g-factor
//const double EDM = 14.6E-17;  // electric momentum g-factor
const double EDM = 7.3E-17;  // electric momentum g-factor

//-----------------------composition method step size---------------------
const double R4[3]={1.3512071919596578, -1.7024143839193155,
                    1.3512071919596578};
const double R6[7]={0.78451361047755726382, 0.23557321335935813368,
                   -1.17767998417887100695, 1.31518632068391122204,
                   -1.17767998417887100695, 0.23557321335935813368,
                    0.78451361047755726382};
//-------------------------state-vector-index------------------------------
enum x_index {x_=0, px_=1, s_= 2, ps_=3, z_=4, pz_=5, vt_=6, dE_=7, Sx_=8, Ss_=9};

//-------------------------state-vector-index------------------------------
enum type_index {DRIFT_=1, eBEND_=2, eQUAD_= 3, RFcav_=4};

//-----------relativistic parameters for electron magic energy-------------
const double magicGAMMA=29.38243572993826;
const double magicBETA=0.9994206777202302;
const double magicBETA2=magicBETA*magicBETA;
const double eMASS = 0.51099891013; // MeV/c^2
const double magicENERGY = magicGAMMA*eMASS; //MeV

//const double offMagicDelta = 1.0E-6;
const double offMagicDelta = 0.0;

const double GAMMA=magicGAMMA*(1.0+offMagicDelta)-offMagicDelta/magicGAMMA;
const double BETA=std::sqrt(1.0-1.0/GAMMA/GAMMA);
const double BETA2=BETA*BETA;
const double BETAGAMMA=GAMMA*BETA; // beta_0 gamma_0
const double BETAGAMMA2=BETAGAMMA*BETAGAMMA;
const double ENERGY = GAMMA*eMASS; //MeV
//-----------
#endif
                
