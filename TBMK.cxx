/////////////////////////////////////////////////////////////////////////////////
//                                                                             //
//  Calculation of DVCS cross section using the BMK formulation of 2002 paper: //
//                                                                             //
//  arXiv:hep-ph/0112108v2                                                     //
//                                                                             //
//  Calculation is done up to Twist 2 approximation                            //
//                                                                             //
//  Written by: Liliet Calero Diaz                                             //
//                                                                             //
//  Email: lc2fc@virginia.edu                                                  //
//                                                                             //
/////////////////////////////////////////////////////////////////////////////////


#include "TBMK.h"

using namespace std;  	// std namespace: so you can do things like 'cout'

ClassImp(TBMK)				// classimp: necessary for root


//____________________________________________________________________________________
TBMK::TBMK() {
	// Default Constructor
}
//____________________________________________________________________________________
TBMK::~TBMK() {
	// Default Destructor
}
//___________________________________________________________________________________________________
TComplex TBMK::cdstar( TComplex c, TComplex d ){ // ( C D* ) product

  TComplex dstar = TComplex::Conjugate(d);

  return ( c.Re() * dstar.Re() - c.Im() * dstar.Im() ) + ( c.Re() * dstar.Im() + c.Im() * dstar.Re() ) * TComplex::I();
}
//____________________________________________________________________________________
void TBMK::SetKinematics(Double_t _QQ, Double_t _x, Double_t _t, Double_t _k){
  QQ = _QQ; //Q^2 value
  x = _x;   //Bjorken x
  t = _t;   //momentum transfer squared
  k = _k;   //Electron Beam energy
  //cout<<"Settings: "<<k<<", "<<QQ<<", "<<x<<", "<<t<<endl;
  ee = 4. * M2 * x * x / QQ; // epsilon squared
  y = sqrt(QQ) / ( sqrt(ee) * k );  // lepton energy fraction
  xi = x * ( 1. + t / 2. / QQ ) / ( 2. - x + x * t / QQ ); // Generalized Bjorken variable
  //Gamma1 = x * y / ALP_INV / ALP_INV / ALP_INV / PI / PI / 16. / QQ / sqrt( 1. + ee ); // factor in front of the cross section, eq. (22)
  Gamma1 = x * y * y / ALP_INV / ALP_INV / ALP_INV / PI / 8. / QQ / QQ / sqrt( 1. + ee ); // factor in front of the cross section, eq. (22)
  s = 2. * M * k + M2;
  Gamma2 = 1. / ALP_INV / ALP_INV / ALP_INV / PI / PI / 16. / ( s - M2 ) / ( s - M2 ) / sqrt( 1. + ee ) / x; // this was for testing Defurne way to write gamma
  tmin = - QQ * ( 2. * ( 1. - x ) * ( 1. - sqrt(1. + ee) ) + ee ) / ( 4. * x * ( 1. - x ) + ee ); // eq. (31)
  K2 = - ( t / QQ ) * ( 1. - x ) * ( 1. - y - y * y * ee / 4.) * ( 1. - tmin / t ) * ( sqrt(1. + ee) + ( ( 4. * x * ( 1. - x ) + ee ) / ( 4. * ( 1. - x ) ) ) * ( ( t - tmin ) / QQ )  ); // eq. (30)

  // Defurne's Jacobian
  jcob = 1./ ( 2. * M * x * k ) * 2. * PI * 2.;// not necessary

}
//___________________________________________________________________________________
void TBMK::BHLeptonPropagators(Double_t phi) {

  // K*D 4-vector product (phi-dependent)
  KD = - QQ / ( 2. * y * ( 1. + ee ) ) * ( 1. + 2. * sqrt(K2) * cos( PI - ( phi * RAD ) ) - t / QQ * ( 1. - x * ( 2. - y ) + y * ee / 2. ) + y * ee / 2.  ); // eq. (29)

  // lepton BH propagators P1 and P2 (contaminating phi-dependence)
  P1 = 1. + 2. * KD / QQ;
  P2 = t / QQ - 2. * KD / QQ;

  //cout << "k * D at "<<phi<<" degres: "<<KD<<endl;
  //cout<<"P1 = "<< P1<<" P2 = "<<P2<<endl;
  }
//___________________________________________________________________________________
Double_t TBMK::BHUU(Double_t phi, Double_t F1, Double_t F2) { // BH Unpolarized Cross Section

  // Get BH propagators
  BHLeptonPropagators(phi);

  // BH unpolarized Fourier harmonics eqs. (35 - 37)
  c0_BH = 8. * K2 * ( ( 2. + 3. * ee ) * ( QQ / t ) * ( F1 * F1  - F2 * F2 * t / ( 4. * M2 ) ) + 2. * x * x * ( F1 + F2 ) * ( F1 + F2 ) ) +
          ( 2. - y ) * ( 2. - y ) * ( ( 2. + ee ) * ( ( 4. * x * x * M2 / t ) * ( 1. + t / QQ ) * ( 1. + t / QQ ) + 4. * ( 1. - x ) * ( 1. + x * t / QQ ) ) * ( F1 * F1 - F2 * F2 * t / ( 4. * M2 ) ) +
          4. * x * x * ( x + ( 1. - x + ee / 2. ) * ( 1. - t / QQ ) * ( 1. - t / QQ ) - x * ( 1. - 2. * x ) * t * t / ( QQ * QQ ) ) * ( F1 + F2 ) * ( F1 + F2 ) ) +
          8. * ( 1. + ee ) * ( 1. - y - ee * y * y / 4. ) * ( 2. * ee * ( 1. - t / ( 4. * M2 ) ) * ( F1 * F1 - F2 * F2 * t / ( 4. * M2 ) ) - x * x * ( 1. - t / QQ ) * ( 1. - t / QQ ) * ( F1 + F2 ) * ( F1 + F2 ) );

  c1_BH = 8. * sqrt(K2) * ( 2. - y ) * ( ( 4. * x * x * M2 / t - 2. * x - ee ) * ( F1 * F1 - F2 * F2 * t / ( 4. * M2 ) ) + 2. * x * x * ( 1. - ( 1. - 2. * x ) * t / QQ ) * ( F1 + F2 ) * ( F1 + F2 ) );

  c2_BH = 8. * x * x * K2 * ( ( 4. * M2 / t ) * ( F1 * F1 - F2 * F2 * t / ( 4. * M2 ) ) + 2. * ( F1 + F2 ) * ( F1 + F2 ) );

  // BH squared amplitude eq (25) divided by e^6
  Amp2_BH = 1. / ( x * x * y * y * ( 1. + ee ) * ( 1. + ee ) * t * P1 * P2 ) * ( c0_BH + c1_BH * cos( PI - (phi * RAD) ) + c2_BH * cos( 2. * ( PI - ( phi * RAD ) ) )  );

  Amp2_BH = GeV2nb * Amp2_BH; // convertion to nb

  return dsigma_BH = Gamma1 * Amp2_BH;
}
//___________________________________________________________________________________________________
Double_t TBMK::DVCSUU( Double_t phi, TComplex t2cffs[4], Double_t &a_dvcs, Double_t &b_dvcs ) { // Pure DVCS Unpolarized Cross Section

  /* t2cffs = { H, E , Htilde, Etilde } Twist-2 Compton Form Factors*/
  TComplex H = t2cffs[0];
  TComplex E = t2cffs[1];
  TComplex Htilde = t2cffs[2];
  TComplex Etilde = t2cffs[3];

  // c coefficients (eq. 66) for pure DVCS .
  Double_t c_dvcs = 1./(2. - x)/(2. - x) * ( 4. * ( 1 - x ) * ( H.Rho2() + Htilde.Rho2() ) - x * x * ( cdstar(H, E) + cdstar(E, H) + cdstar(Htilde, Etilde) + cdstar(Etilde, Htilde) ) -
                    ( x * x + (2. - x) * (2. - x) * t / 4. / M2 ) * E.Rho2() - ( x * x * t / 4. / M2 ) * Etilde.Rho2() );

  // Pure DVCS unpolarized Fourier harmonics eqs. (43 - 44)
  c0_dvcs = 2. * ( 2. - 2. * y + y * y ) * c_dvcs;
  c1_dvcs = - ( 2. * xi / ( 1. + xi ) ) * ( 8. * sqrt(K2) / ( 2. - x ) ) * ( 2. - y ) * c_dvcs;

  // DVCS squared amplitude eq (26) divided by e^6
  Amp2_DVCS = 1. / ( y * y * QQ ) * ( c0_dvcs + c1_dvcs * cos( PI - (phi * RAD) ) );

  Amp2_DVCS = GeV2nb * Amp2_DVCS; // convertion to nb

  // DVCS has the form a_dvcs + b_dvcs cos(PI-phi)
  a_dvcs = GeV2nb * Gamma1 / ( y * y * QQ ) * c0_dvcs;
  b_dvcs = GeV2nb * Gamma1 / ( y * y * QQ ) * c1_dvcs;

  return dsigma_DVCS = Gamma1 * Amp2_DVCS;
}
//___________________________________________________________________________________________________
Double_t TBMK::DVCSUU2( TComplex t2cffs[4] ) { // Pure DVCS Unpolarized Cross Section with just considering thre c0 term (as supposely done by liuti)

  /* t2cffs = { H, E , Htilde, Etilde } Twist-2 Compton Form Factors*/
  TComplex H = t2cffs[0];
  TComplex E = t2cffs[1];
  TComplex Htilde = t2cffs[2];
  TComplex Etilde = t2cffs[3];

  // c coefficients (eq. 66) for pure DVCS .
  Double_t c_dvcs = 1./(2. - x)/(2. - x) * ( 4. * ( 1 - x ) * ( H.Rho2() + Htilde.Rho2() ) - x * x * ( cdstar(H, E) + cdstar(E, H) + cdstar(Htilde, Etilde) + cdstar(Etilde, Htilde) ) -
                    ( x * x + (2. - x) * (2. - x) * t / 4. / M2 ) * E.Rho2() - ( x * x * t / 4. / M2 ) * Etilde.Rho2() );

  // Pure DVCS unpolarized Fourier harmonics eqs. (43 - 44)
  c0_dvcs = 2. * ( 2. - 2. * y + y * y ) * c_dvcs;
  //c1_dvcs = - ( 2. * xi / ( 1. + xi ) ) * ( 8. * sqrt(K2) / ( 2. - x ) ) * ( 2. - y ) * c_dvcs;

  // DVCS squared amplitude eq (26) divided by e^6
  Amp2_DVCS = 1. / ( y * y * QQ ) *  c0_dvcs ;

  Amp2_DVCS = GeV2nb * Amp2_DVCS; // convertion to nb

  return dsigma_DVCS = Gamma1 * Amp2_DVCS;
}
//___________________________________________________________________________________________________
Double_t TBMK::IUU(Double_t phi, Double_t F1, Double_t F2, TComplex t2cffs[4]) { // Interference Unpolarized Cross Section (not easy to extract CFFs if written in this way, but used for cross check)

  // Get BH propagators
  BHLeptonPropagators(phi);

  /* t2cffs_I = { H, E , Htilde } Twist-2 Compton Form Factors present in the interference term*/
  TComplex H = t2cffs[0];
  TComplex E = t2cffs[1];
  TComplex Htilde = t2cffs[2];
  TComplex Etilde = t2cffs[3]; // This CFF does not appear in the interference

  // c (eq. 69) and deltac (eq. 72) interference coefficients.
  TComplex c_intf, dc_intf;
  c_intf = F1 * H + ( x / (2. - x) ) * ( F1 + F2 ) * Htilde - ( t / 4. / M2 ) * F2 * E; // eq. 69
  dc_intf = - ( x / (2. - x) ) * ( F1 + F2 ) * ( ( x / (2. - x) ) * ( H + E ) + Htilde ); // eq. 72

  // Interference unpolarized Fourier harmonics eqs. (53 - 55)
  c0_I = - 8. * ( 2. - y ) * ( ( 2. - y ) * ( 2. - y ) / ( 1. - y ) * K2 * c_intf.Re() + t / QQ * ( 1. - y ) * ( 2. - x ) * ( c_intf.Re() + dc_intf.Re() ) );
  c1_I = - 8. * sqrt(K2) * ( 2. - 2. * y + y * y ) * c_intf.Re();
  s1_I = 8. * sqrt(K2) * y * ( 2. - y ) * c_intf.Im(); // here lambda(lepton helicity) was set to 1 (= +1 if the spin is aligned with the direction of the lepton 3-momentum)
  c2_I = - 16. * K2 * ( 2. - y ) / ( 2. - x ) * ( - ( 2. * xi / ( 1. + xi ) ) * c_intf.Re() );
  s2_I = 16. * K2 * y / ( 2. - x ) * ( - ( 2. * xi / ( 1. + xi ) ) * c_intf.Im() ); // here lambda(lepton helicity) was set to 1 (= +1 if the spin is aligned with the direction of the lepton 3-momentum)

  // BH-DVCS interference squared amplitude eq (27) divided by e^6
  Amp2_I = 1. / ( x * y * y * y * t * P1 * P2 ) * ( c0_I + c1_I * cos( PI - (phi * RAD) ) + s1_I * sin( PI - (phi * RAD) ) + c2_I * cos( 2. * ( PI - ( phi * RAD ) ) ) + s2_I * sin( 2. * ( PI - ( phi * RAD ) ) )  );

  Amp2_I = GeV2nb * Amp2_I; // convertion to nb

  //cout<<"c0 = "<<c0_I<<", s2 = "<<s2_I<<endl;

  return dsigma_I = Gamma1 * Amp2_I;
}
//___________________________________________________________________________________________________
Double_t TBMK::IUU2(Double_t phi, Double_t F1, Double_t F2, TComplex t2cffs[4]) { // Interference Unpolarized Cross Section (writting as Liuti's style, easier extraction)

  // Get BH propagators
  BHLeptonPropagators(phi);

  /* t2cffs_I = { H, E , Htilde } Twist-2 Compton Form Factors present in the interference term*/
  TComplex H = t2cffs[0];
  TComplex E = t2cffs[1];
  TComplex Htilde = t2cffs[2];
  TComplex Etilde = t2cffs[3]; // This CFF does not appear in the interference

  Double_t A, B, C, D; //Coefficients in from to the CFFs

  A = - 8. * K2 * ( 2. - y ) * ( 2. - y ) * ( 2. - y ) / ( 1. - y ) - 8. * ( 2. - y ) * ( 1. - y ) * ( 2. - x ) * t / QQ - 8. * sqrt(K2) * ( 2. - 2. * y + y * y ) * cos( PI - (phi * RAD) )
      + 32. * K2 * xi * ( 2. - y ) / ( 2. - x ) / ( 1. + xi ) * cos( 2. * ( PI - ( phi * RAD ) ) );
  // no c2 term
  //A = - 8. * K2 * ( 2. - y ) * ( 2. - y ) * ( 2. - y ) / ( 1. - y ) - 8. * ( 2. - y ) * ( 1. - y ) * ( 2. - x ) * t / QQ - 8. * sqrt(K2) * ( 2. - 2. * y + y * y ) * cos( PI - (phi * RAD) );
  B = 8. * x * x * ( 2. - y ) * (1 - y ) / ( 2. - x ) * t / QQ;
  C = x / ( 2. - x ) * ( A + ( 2. - x ) * ( 2. - x) / x / x * B );
  D = 8. * sqrt(K2) * y * ( 2. - y ) * sin( PI - (phi * RAD) ) - 32. * K2 * y * xi / ( 2. - x ) / ( 1 + xi ) * sin( 2. * ( PI - ( phi * RAD ) ) );
  //D = 0.;

  // BH-DVCS interference squared amplitude eq (27) divided by e^6
  Amp2_I2 = 1. / ( x * y * y * y * t * P1 * P2 ) * ( A * ( F1 * H.Re() - t / 4. / M2 * F2 * E.Re() ) + B * ( F1 + F2 ) * ( H.Re() + E.Re() ) + C * ( F1 + F2 ) * Htilde.Re() + D * ( F1 * H.Im() - t / 4. / M2 * F2 * E.Im() + x / ( 2 - x ) * ( F1 + F2 ) * Htilde.Im() ) );

  Amp2_I2 = GeV2nb * Amp2_I2; // convertion to nb

  return dsigma_I2 = Gamma1 * Amp2_I2;

}
//___________________________________________________________________________________________________
Double_t TBMK::IUU3(Double_t phi, Double_t F1, Double_t F2, TComplex t2cffs[4]) { // Interference Unpolarized Cross Section (writting it as in Liuti's comparison paper, i.e just c0 and c1 terms)

  // Get BH propagators
  BHLeptonPropagators(phi);

  /* t2cffs_I = { H, E , Htilde } Twist-2 Compton Form Factors present in the interference term*/
  TComplex H = t2cffs[0];
  TComplex E = t2cffs[1];
  TComplex Htilde = t2cffs[2];
  TComplex Etilde = t2cffs[3]; // This CFF does not appear in the interference

  Double_t A, B, C, D; //Coefficients in from to the CFFs

  A = - 8. * K2 * ( 2. - y ) * ( 2. - y ) * ( 2. - y ) / ( 1. - y ) - 8. * ( 2. - y ) * ( 1. - y ) * ( 2. - x ) * t / QQ - 8. * sqrt(K2) * ( 2. - 2. * y + y * y ) * cos( PI - (phi * RAD) );
  B = 8. * x * x * ( 2. - y ) * (1 - y ) / ( 2. - x ) * t / QQ;
  C =  x / ( 2. - x ) * ( - 8. * K2 * ( 2. - y ) * ( 2. - y ) * ( 2. - y ) / ( 1. - y ) - 8. * sqrt(K2) * ( 2. - 2. * y + y * y ) * cos( PI - (phi * RAD) ) );
  //D = 8. * sqrt(K2) * y * ( 2. - y ) * sin( PI - (phi * RAD) ) - 32. * K2 * y * xi / ( 2. - x ) / ( 1 + xi ) * sin( 2. * ( PI - ( phi * RAD ) ) );
  D = 0.;

  // BH-DVCS interference squared amplitude eq (27) divided by e^6
  Amp2_I2 = 1. / ( x * y * y * y * t * P1 * P2 ) * ( A * ( F1 * H.Re() - t / 4. / M2 * F2 * E.Re() ) + B * ( F1 + F2 ) * ( H.Re() + E.Re() ) + C * ( F1 + F2 ) * Htilde.Re() + D * ( F1 * H.Im() - t / 4. / M2 * F2 * E.Im() + x / ( 2 - x ) * ( F1 + F2 ) * Htilde.Im() ) );

  Amp2_I2 = GeV2nb * Amp2_I2; // convertion to nb

  return dsigma_I2 = Gamma1 * Amp2_I2;

}
