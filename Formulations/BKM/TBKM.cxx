//////////////////////////////////////////////////////////////////////////////////
//                                                                              //
//  Calculation of unpolarized DVCS cross section using the BKM formulation:    //
//                                                                              //
//  - BKM 2002 cross sections (arXiv:hep-ph/0112108v2)                          //
//      -- BH, DVCDS_UU_02, I_UU_02                                             //
//                                                                              //
//  - BKM 2010 cross sections (arXiv:hep-ph/1005.5209v1)                        //
//      -- DVCS_UU_10, I_UU_10                                                  //
//                                                                              //
//  Written by: Liliet Calero Diaz                                              //
//                                                                              //
//  Email: lc2fc@virginia.edu                                                   //
//                                                                              //
//////////////////////////////////////////////////////////////////////////////////


#include "TBKM.h"

using namespace std;  	// std namespace: so you can do things like 'cout'

ClassImp(TBKM)				// classimp: necessary for root


//____________________________________________________________________________________
TBKM::TBKM() {
	// Default Constructor
}
//____________________________________________________________________________________
TBKM::~TBKM() {
	// Default Destructor
}
//___________________________________________________________________________________________________
TComplex TBKM::cdstar( TComplex c, TComplex d ){ // ( C D* ) product

    TComplex dstar = TComplex::Conjugate(d);

    return ( c.Re() * dstar.Re() - c.Im() * dstar.Im() ) + ( c.Re() * dstar.Im() + c.Im() * dstar.Re() ) * TComplex::I();
}
//____________________________________________________________________________________
void TBKM::SetKinematics( Double_t *kine ) {

    QQ = kine[0];     //Q^2 value
    x = kine[1];      //Bjorken x
    t = kine[2];      //momentum transfer squared
    k = kine[3];      //Electron Beam energy

    ee = 4. * M2 * x * x / QQ; // epsilon squared
    y = sqrt(QQ) / ( sqrt(ee) * k );  // lepton energy fraction
    xi = x * ( 1. + t / 2. / QQ ) / ( 2. - x + x * t / QQ ); // Generalized Bjorken variable
    s = 2. * M * k + M2; // Mandelstan variable
    Gamma = x * y * y / ALP_INV / ALP_INV / ALP_INV / PI / 8. / QQ / QQ / sqrt( 1. + ee ); // factor in front of the cross section, eq. (22)
    tmin = - QQ * ( 2. * ( 1. - x ) * ( 1. - sqrt(1. + ee) ) + ee ) / ( 4. * x * ( 1. - x ) + ee ); // eq. (31)
    Ktilde_10 = sqrt( tmin - t ) * sqrt( ( 1. - x ) * sqrt( 1. + ee ) + ( ( t - tmin ) * ( ee + 4. * x * ( 1. - x ) ) / 4. / QQ ) ) * sqrt( 1. - y - y * y * ee / 4. ) / sqrt( 1. - y + y * y * ee / 4.); // K tilde from 2010 paper
    K = sqrt( 1. - y + y * y * ee / 4.) * Ktilde_10 / sqrt(QQ);

}
//___________________________________________________________________________________
void TBKM::BHLeptonPropagators(Double_t *kine, Double_t phi) {

    SetKinematics(kine);
    // K*D 4-vector product (phi-dependent)
    KD = - QQ / ( 2. * y * ( 1. + ee ) ) * ( 1. + 2. * K * cos( PI - ( phi * RAD ) ) - t / QQ * ( 1. - x * ( 2. - y ) + y * ee / 2. ) + y * ee / 2.  ); // eq. (29)

    // lepton BH propagators P1 and P2 (contaminating phi-dependence)
    P1 = 1. + 2. * KD / QQ;
    P2 = t / QQ - 2. * KD / QQ;

    //cout << "k * D at "<<phi<<" degres: "<<KD<<endl;
    cout<<"P1 = "<< P1<<" P2 = "<<P2<<endl;
}
//___________________________________________________________________________________
Double_t TBKM::BH_UU(Double_t *kine, Double_t phi, Double_t F1, Double_t F2) { // BH Unpolarized Cross Section

    // Sets the kinematics and gets BH propagators
    BHLeptonPropagators(kine, phi);

    // BH unpolarized Fourier harmonics eqs. (35 - 37)
    c0_BH = 8. * K * K * ( ( 2. + 3. * ee ) * ( QQ / t ) * ( F1 * F1  - F2 * F2 * t / ( 4. * M2 ) ) + 2. * x * x * ( F1 + F2 ) * ( F1 + F2 ) ) +
          ( 2. - y ) * ( 2. - y ) * ( ( 2. + ee ) * ( ( 4. * x * x * M2 / t ) * ( 1. + t / QQ ) * ( 1. + t / QQ ) + 4. * ( 1. - x ) * ( 1. + x * t / QQ ) ) * ( F1 * F1 - F2 * F2 * t / ( 4. * M2 ) ) + 4. * x * x * ( x + ( 1. - x + ee / 2. ) * ( 1. - t / QQ ) * ( 1. - t / QQ ) - x * ( 1. - 2. * x ) * t * t / ( QQ * QQ ) ) * ( F1 + F2 ) * ( F1 + F2 ) ) +
          8. * ( 1. + ee ) * ( 1. - y - ee * y * y / 4. ) * ( 2. * ee * ( 1. - t / ( 4. * M2 ) ) * ( F1 * F1 - F2 * F2 * t / ( 4. * M2 ) ) - x * x * ( 1. - t / QQ ) * ( 1. - t / QQ ) * ( F1 + F2 ) * ( F1 + F2 ) );

    c1_BH = 8. * K * ( 2. - y ) * ( ( 4. * x * x * M2 / t - 2. * x - ee ) * ( F1 * F1 - F2 * F2 * t / ( 4. * M2 ) ) + 2. * x * x * ( 1. - ( 1. - 2. * x ) * t / QQ ) * ( F1 + F2 ) * ( F1 + F2 ) );

    c2_BH = 8. * x * x * K * K * ( ( 4. * M2 / t ) * ( F1 * F1 - F2 * F2 * t / ( 4. * M2 ) ) + 2. * ( F1 + F2 ) * ( F1 + F2 ) );

    // BH squared amplitude eq (25) divided by e^6
    Amp2_BH = 1. / ( x * x * y * y * ( 1. + ee ) * ( 1. + ee ) * t * P1 * P2 ) * ( c0_BH + c1_BH * cos( PI - (phi * RAD) ) + c2_BH * cos( 2. * ( PI - ( phi * RAD ) ) )  );

    Amp2_BH = GeV2nb * Amp2_BH; // convertion to nb

    return dsigma_BH = Gamma * Amp2_BH;
}
//===================================================================================================
//                   B   K   M   -   2   0   0   2
//===================================================================================================
//___________________________________________________________________________________________________
Double_t TBKM::DVCS_UU_02(Double_t *kine, Double_t phi, TComplex t2cffs[4], TString twist = "T2") { // Pure DVCS Unpolarized Cross Section

    SetKinematics(kine);

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
    c1_dvcs = - ( 2. * xi / ( 1. + xi ) ) * ( 8. * K / ( 2. - x ) ) * ( 2. - y ) * c_dvcs;

    // DVCS squared amplitude eq (26) divided by e^6
    if( twist == "T2") // F_eff = 0
        Amp2_DVCS = 1. / ( y * y * QQ ) *  c0_dvcs ;

    if( twist == "T3")  // F_eff = -2xi/(1+xi) F
        Amp2_DVCS = 1. / ( y * y * QQ ) * ( c0_dvcs + c1_dvcs * cos( PI - (phi * RAD) ) );

    Amp2_DVCS = GeV2nb * Amp2_DVCS; // convertion to nb

    return dsigma_DVCS = Gamma * Amp2_DVCS;
}
//___________________________________________________________________________________________________
Double_t TBKM::I_UU_02(Double_t *kine, Double_t phi, Double_t F1, Double_t F2, TComplex t2cffs[4], TString twist = "T2") { // Interference Unpolarized Cross Section (Liuti's style)

    // Get BH propagators and set the kinematics
    BHLeptonPropagators(kine, phi);

    /* t2cffs_I = { H, E , Htilde, Etilde } Twist-2 Compton Form Factors ( Etilde does not appear in the interference term ) */
    TComplex H = t2cffs[0];
    TComplex E = t2cffs[1];
    TComplex Htilde = t2cffs[2];
    TComplex Etilde = t2cffs[3]; // This CFF does not appear in the interference

    Double_t A, B, C; // Coefficients in from to the CFFs

    if( twist == "T2") // F_eff = 0, no c2 term (no cos(2phi))
        A = - 8. * K * K * ( 2. - y ) * ( 2. - y ) * ( 2. - y ) / ( 1. - y ) - 8. * ( 2. - y ) * ( 1. - y ) * ( 2. - x ) * t / QQ - 8. * K * ( 2. - 2. * y + y * y ) * cos( PI - (phi * RAD) );

    if( twist == "T3") // F_eff = -2xi/(1+xi) F
        A = - 8. * K * K * ( 2. - y ) * ( 2. - y ) * ( 2. - y ) / ( 1. - y ) - 8. * ( 2. - y ) * ( 1. - y ) * ( 2. - x ) * t / QQ - 8. * K * ( 2. - 2. * y + y * y ) * cos( PI - (phi * RAD) )
            + 32. * K * K * xi * ( 2. - y ) / ( 2. - x ) / ( 1. + xi ) * cos( 2. * ( PI - ( phi * RAD ) ) );

    B = 8. * x * x * ( 2. - y ) * (1 - y ) / ( 2. - x ) * t / QQ;
    // C =  x / ( 2. - x ) * ( - 8. * K * K * ( 2. - y ) * ( 2. - y ) * ( 2. - y ) / ( 1. - y ) - 8. * K * ( 2. - 2. * y + y * y ) * cos( PI - (phi * RAD) ) );
    C = x / ( 2. - x ) * ( A + ( 2. - x ) * ( 2. - x) / x / x * B );

    // BH-DVCS interference squared amplitude eq (27) divided by e^6
    I = 1. / ( x * y * y * y * t * P1 * P2 ) * ( A * ( F1 * H.Re() - t / 4. / M2 * F2 * E.Re() ) + B * ( F1 + F2 ) * ( H.Re() + E.Re() ) + C * ( F1 + F2 ) * Htilde.Re() );

    I = GeV2nb * I; // convertion to nb

    return dsigma_I = Gamma * I;
}
//===================================================================================================
//                   B   K   M   -   2   0   1   0
//===================================================================================================
