#include <iostream>
#include <fstream>
#include <string>

using namespace std;

ROOT::Math::Integrator ig;

Double_t PI = TMath::Pi();

// KM09b parameters ---------------------------------------------

// valence H
Double_t nval = 1.35;
//Double_t Mval = sqrt(0.64);
Double_t pval = 1.;
//Double_t rval = 1.11;
//Double_t bval = 2.4;

// sea H
Double_t nsea = 1.5;
Double_t rsea = 1.;
//Double_t Msea = sqrt(0.5);
Double_t psea = 2.;
Double_t bsea = 4.6;

// subtraction constant (H, E)
//Double_t C0 = 6.;
//Double_t Msub = 1.5;

// valence Htilde
//Double_t ntval = 3.;
//Double_t rtval = 1.;
//Double_t btval = 3./2.;
//Double_t Mtval = Mval;

// KM15 parameters ----------------------------------------------

// valence H
Double_t Mval = 0.789;
Double_t rval = 0.918;
Double_t bval = 0.4;

// subtraction constant (H, E)
Double_t C0 = 2.768;
Double_t Msub = 1.204;

// valence Htilde
Double_t Mtval = 3.993;
Double_t rtval = 0.881;
Double_t btval = 0.4;
Double_t ntval = 0.6;

// sea H
Double_t Msea = sqrt(0.482);

// "pion pole" Etilde
Double_t rpi = 2.646;
Double_t Mpi = 4.;


void ModKM15_CFFs(Double_t *kine, Double_t &ReH, Double_t &ImH, Double_t &ReE, Double_t &ReHtilde, Double_t &ImHtilde, Double_t &ReEtilde) {

   Double_t QQ = kine[0];     //Q^2 value
   Double_t xB = kine[1];      //Bjorken x
   Double_t t = kine[2];      //momentum transfer squared
   Double_t k = kine[3];      //Electron Beam energy

   Double_t alpha_val = 0.43 + 0.85 * t;
   Double_t alpha_sea = 1.13 + 0.15 * t;

   Double_t Ct = C0 / pow(1. - t / Msub / Msub, 2.);

   Double_t xi = xB / ( 2. - xB ); // Generalized Bjorken variable KM09
   //Double_t xi = xB * ( 1. + t / 2. / QQ ) / ( 2. - xB + xB * t / QQ ); // Generalized Bjorken variable

   // Imaginary part of the CFFs (H, Htilde)

   auto fHval = [&](double x){ return (nval * rval) / ( 1. + x ) * pow( ( 2. * x ) / ( 1. + x ) , -alpha_val ) * pow( ( 1. - x ) / ( 1. + x ), bval ) * 1. / pow( 1.- ( ( 1. - x ) / ( 1. + x ) ) * ( t / Mval / Mval ), pval); };

   auto fHsea = [&](double x){ return (nsea * rsea) / ( 1. + x ) * pow( ( 2. * x ) / ( 1. + x ) , -alpha_sea ) * pow( ( 1. - x ) / ( 1. + x ), bsea ) * 1. / pow( 1.- ( ( 1. - x ) / ( 1. + x ) ) * ( t / Msea / Msea ), psea); };

   auto fHtval = [&](double x){ return (ntval * rtval) / ( 1. + x ) * pow( ( 2. * x ) / ( 1. + x ) , -alpha_val ) * pow( ( 1. - x ) / ( 1. + x ), btval ) * 1. / ( 1.- ( ( 1. - x ) / ( 1. + x ) ) * ( t / Mtval / Mtval ) ); };

   auto fImH = [&](double x){ return PI * ( ( 2. * ( 4. / 9. ) + 1. / 9. ) * fHval(x) + 2. / 9. * fHsea(x)); };

   auto fImHt = [&](double x){ return PI * ( 2. * ( 4. / 9. ) + 1. / 9. ) * fHtval(x); };

   // Real part of the CFFs (H, Htilde)

   auto fPV_ReH = [&](double x){ return -2.* x / ( x + xi ) * fImH(x); };
   auto fPV_ReHt = [&](double x){ return -2.* xi / ( x + xi ) * fImHt(x); };

   ig.SetFunction(fPV_ReH);
   double DR_ReH = ig.IntegralCauchy(1.e-279, 1., xi);

   ig.SetFunction(fPV_ReHt);
   double DR_ReHt = ig.IntegralCauchy(1.e-323, 1., xi);

   // Evaluate the CFFs

   ImH = fImH(xi);
   ReH = 1. / PI * DR_ReH - Ct;
   ReE = Ct;
   ImHtilde = fImHt(xi);
   ReHtilde = 1. / PI * DR_ReHt;
   ReEtilde = rpi / xi * 2.164 / ( ( 0.0196 - t ) * pow( 1. - t / Mpi / Mpi, 2.) );

}
