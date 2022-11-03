//////////////////////////////////////////////////////////////////////////////////
//                                                                              //
//  Based on unpolarized DVCS cross section using the BKM formulation by        //
//  Liliet Calero Diaz  -  lc2fc@virginia.edu                                   //                             
//                                                                              //
//////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////
//                                                                              //
//  Calculation of polarized DVCS cross section using the BKM formulation:      //
//                                                                              //
//  - BKM 2002 cross sections (arXiv:hep-ph/0112108v2)                          //
//      -- BH, DVCDS_LP_02, I_LP_02                                             //
//                                                                              //
//  - BKM 2010 cross sections (arXiv:hep-ph/1005.5209v1)                        //
//      -- DVCS_LP_10, I_LP_10                                                  //
//                                                                              //
//  Written by: Duncan Beauch and Wyndham White                                 //
//                                                                              //
//  Email: drb5wqd@virginia.edu & wrw2ztk@virginia.edu                          //
//                                                                              //
//////////////////////////////////////////////////////////////////////////////////

#include "TBKM.h"


using namespace std;  	// std namespace: so you can do things like 'cout'

ClassImp(TBKM)			// classimp: necessary for root


//_______________________________________________________________________________________________________________________________
TBKM::TBKM() {
	// Default Constructor
}
//_______________________________________________________________________________________________________________________________
TBKM::~TBKM() {
	// Default Destructor
}
//_______________________________________________________________________________________________________________________________
TComplex TBKM::cdstar( TComplex c, TComplex d ){ // ( C D* ) product

    TComplex dstar = TComplex::Conjugate(d);

    return ( c.Re() * dstar.Re() - c.Im() * dstar.Im() ) + ( c.Re() * dstar.Im() + c.Im() * dstar.Re() ) * TComplex::I();
}
//_______________________________________________________________________________________________________________________________
void TBKM::SetCFFs( TComplex *t2cffs ) { // Twist-2 Compton Form Factors

     H = t2cffs[0];
     E = t2cffs[1];
     Htilde = t2cffs[2];
     Etilde = t2cffs[3];
}
//_______________________________________________________________________________________________________________________________
void TBKM::SetKinematics( Double_t *kine ) {

    QQ = kine[0];     //Q^2 value
    x = kine[1];      //Bjorken x
    t = kine[2];      //momentum transfer squared
    k = kine[3];      //Electron Beam energy

    ee = 4. * M2 * x * x / QQ; // epsilon squared
    sqrtOnePlusEE = sqrt(1 + ee); // square root of one plus epsilon squared
    y = sqrt(QQ) / ( sqrt(ee) * k );  // lepton energy fraction
    xi = x * ( 1. + t / 2. / QQ ) / ( 2. - x + x * t / QQ ); // Generalized Bjorken variable
    s = 2. * M * k + M2; // Mandelstan variable
    Gamma = x * y * y / ALP_INV / ALP_INV / ALP_INV / PI / 8. / QQ / QQ / sqrt( 1. + ee ); // factor in front of the cross section, eq. (22)
    tmin = - QQ * ( 2. * ( 1. - x ) * ( 1. - sqrt(1. + ee) ) + ee ) / ( 4. * x * ( 1. - x ) + ee ); // eq. (31)
    Ktilde_10 = sqrt( tmin - t ) * sqrt( ( 1. - x ) * sqrt( 1. + ee ) + ( ( t - tmin ) * ( ee + 4. * x * ( 1. - x ) ) / 4. / QQ ) ) * sqrt( 1. - y - y * y * ee / 4. )
                / sqrt( 1. - y + y * y * ee / 4.); // K tilde from 2010 paper
    K = sqrt( 1. - y + y * y * ee / 4.) * Ktilde_10 / sqrt(QQ);
}
//_______________________________________________________________________________________________________________________________
void TBKM::BHLeptonPropagators(Double_t *kine, Double_t phi) {

    SetKinematics(kine);
    // K*D 4-vector product (phi-dependent)
    KD = - QQ / ( 2. * y * ( 1. + ee ) ) * ( 1. + 2. * K * cos( PI - ( phi * RAD ) ) - t / QQ * ( 1. - x * ( 2. - y ) + y * ee / 2. ) + y * ee / 2.  ); // eq. (29)

    // lepton BH propagators P1 and P2 (contaminating phi-dependence)
    P1 = 1. + 2. * KD / QQ;
    P2 = t / QQ - 2. * KD / QQ;
}
//_______________________________________________________________________________________________________________________________
Double_t TBKM::BH_LP(Double_t *kine, Double_t phi, Double_t F1, Double_t F2) { // BH Unpolarized Cross Section

    // Sets the kinematics and gets BH propagators
    BHLeptonPropagators(kine, phi);

    // BH unpolarized Fourier harmonics eqs. (35 - 37)
    c0_BH = 8. * K * K * ( ( 2. + 3. * ee ) * ( QQ / t ) * ( F1 * F1  - F2 * F2 * t / ( 4. * M2 ) ) + 2. * x * x * ( F1 + F2 ) * ( F1 + F2 ) )
            + ( 2. - y ) * ( 2. - y ) * ( ( 2. + ee ) * ( ( 4. * x * x * M2 / t ) * ( 1. + t / QQ ) * ( 1. + t / QQ ) + 4. * ( 1. - x ) * ( 1. + x * t / QQ ) ) * ( F1 * F1 - F2 * F2 * t / ( 4. * M2 ) )
            + 4. * x * x * ( x + ( 1. - x + ee / 2. ) * ( 1. - t / QQ ) * ( 1. - t / QQ ) - x * ( 1. - 2. * x ) * t * t / ( QQ * QQ ) ) * ( F1 + F2 ) * ( F1 + F2 ) )
            + 8. * ( 1. + ee ) * ( 1. - y - ee * y * y / 4. ) * ( 2. * ee * ( 1. - t / ( 4. * M2 ) ) * ( F1 * F1 - F2 * F2 * t / ( 4. * M2 ) ) - x * x * ( 1. - t / QQ ) * ( 1. - t / QQ ) * ( F1 + F2 ) * ( F1 + F2 ) );

    c1_BH = 8. * K * ( 2. - y ) * ( ( 4. * x * x * M2 / t - 2. * x - ee ) * ( F1 * F1 - F2 * F2 * t / ( 4. * M2 ) ) + 2. * x * x * ( 1. - ( 1. - 2. * x ) * t / QQ ) * ( F1 + F2 ) * ( F1 + F2 ) );

    c2_BH = 8. * x * x * K * K * ( ( 4. * M2 / t ) * ( F1 * F1 - F2 * F2 * t / ( 4. * M2 ) ) + 2. * ( F1 + F2 ) * ( F1 + F2 ) );

    // BH squared amplitude eq (25) divided by e^6
    Amp2_BH = 1. / ( x * x * y * y * ( 1. + ee ) * ( 1. + ee ) * t * P1 * P2 ) * ( c0_BH + c1_BH * cos( PI - (phi * RAD) ) + c2_BH * cos( 2. * ( PI - ( phi * RAD ) ) )  );

    Amp2_BH = GeV2nb * Amp2_BH; // convertion to nb

    return dsigma_BH = Gamma * Amp2_BH;
}
//===============================================================================================================================
//                   B   K   M   -   2   0   0   2
//===============================================================================================================================
//_______________________________________________________________________________________________________________________________
Double_t TBKM::DVCS_LP_02(Double_t *kine, Double_t phi, TComplex *t2cffs, TString twist = "t2") { // Pure DVCS Unpolarized Cross Section

    SetKinematics(kine);

    SetCFFs(t2cffs);

    // c coefficients (BKM02 eqs. [66]) for pure DVCS
    c_dvcs = 1./(2. - x)/(2. - x) * ( 4. * ( 1 - x ) * ( H.Rho2() + Htilde.Rho2() ) - x * x * ( cdstar(H, E) + cdstar(E, H) + cdstar(Htilde, Etilde) + cdstar(Etilde, Htilde) )
             - ( x * x + (2. - x) * (2. - x) * t / 4. / M2 ) * E.Rho2() - ( x * x * t / 4. / M2 ) * Etilde.Rho2() );

    // Pure DVCS unpolarized Fourier harmonics (BKM02 eqs. [43, 44])
    c0_dvcs = 2. * ( 2. - 2. * y + y * y ) * c_dvcs;
    c1_dvcs = - ( 2. * xi / ( 1. + xi ) ) * ( 8. * K / ( 2. - x ) ) * ( 2. - y ) * c_dvcs;

    // DVCS squared amplitude eq (26) divided by e^6
    if( twist == "t2") // F_eff = 0
        Amp2_DVCS = 1. / ( y * y * QQ ) *  c0_dvcs ;

    if( twist == "t3")  // F_eff = -2xi/(1+xi) F
        Amp2_DVCS = 1. / ( y * y * QQ ) * ( c0_dvcs + c1_dvcs * cos( PI - (phi * RAD) ) );

    Amp2_DVCS = GeV2nb * Amp2_DVCS; // convertion to nb

    return dsigma_DVCS = Gamma * Amp2_DVCS;
}
//_______________________________________________________________________________________________________________________________
Double_t TBKM::I_LP_02(Double_t *kine, Double_t phi, Double_t F1, Double_t F2, TComplex *t2cffs, TString twist = "t2") { // Interference Unpolarized Cross Section (Liuti's style)

    // Get BH propagators and set the kinematics
    BHLeptonPropagators(kine, phi);

    SetCFFs(t2cffs); // Etilde CFF does not appear in the interference

    if( twist == "t2") // F_eff = 0, no c2 term (no cos(2phi))
        A_02 = - 8. * K * K * ( 2. - y ) * ( 2. - y ) * ( 2. - y ) / ( 1. - y ) - 8. * ( 2. - y ) * ( 1. - y ) * ( 2. - x ) * t / QQ - 8. * K * ( 2. - 2. * y + y * y ) * cos( PI - (phi * RAD) );

    if( twist == "t3") // F_eff = -2xi/(1+xi) F
        A_02 = - 8. * K * K * ( 2. - y ) * ( 2. - y ) * ( 2. - y ) / ( 1. - y ) - 8. * ( 2. - y ) * ( 1. - y ) * ( 2. - x ) * t / QQ - 8. * K * ( 2. - 2. * y + y * y ) * cos( PI - (phi * RAD) )
               + 32. * K * K * xi * ( 2. - y ) / ( 2. - x ) / ( 1. + xi ) * cos( 2. * ( PI - ( phi * RAD ) ) );

    B_02 = 8. * x * x * ( 2. - y ) * (1 - y ) / ( 2. - x ) * t / QQ;
    // C =  x / ( 2. - x ) * ( - 8. * K * K * ( 2. - y ) * ( 2. - y ) * ( 2. - y ) / ( 1. - y ) - 8. * K * ( 2. - 2. * y + y * y ) * cos( PI - (phi * RAD) ) );
    C_02 = x / ( 2. - x ) * ( A_02 + ( 2. - x ) * ( 2. - x) / x / x * B_02 );

    // BH-DVCS interference squared amplitude eq (27) divided by e^6
    I = 1. / ( x * y * y * y * t * P1 * P2 ) * ( A_02 * ( F1 * H.Re() - t / 4. / M2 * F2 * E.Re() ) + B_02 * ( F1 + F2 ) * ( H.Re() + E.Re() ) + C_02 * ( F1 + F2 ) * Htilde.Re() );

    I = GeV2nb * I; // convertion to nb

    return dsigma_I = Gamma * I;
}
//===============================================================================================================================
//                   B   K   M   -   2   0   1   0
//===============================================================================================================================
//_______________________________________________________________________________________________________________________________

//Wyndham done with this
Double_t TBKM::DVCS_LP_10(Double_t *kine, Double_t phi, TComplex *t2cffs, TString twist = "t2") { // Pure DVCS Polarized Cross Section

    SetKinematics(kine);

    SetCFFs(t2cffs);

    // F_eff = f * F
    if(twist == "t2") f = 0; // F_eff = 0 --> DVCS is constant
    if(twist == "t3") f = - 2. * xi / ( 1. + xi );
    if(twist == "t3ww") f = 2. / ( 1. + xi );

    //MODIFIED FROM HERE
    // c_dvcs_polarized(F,F*) coefficients (BKM10 eqs. [2.23]) for pure DVCS
    c_dvcs_ffs = QQ * ( QQ + x * t ) / (sqrtOnePlusEE * pow( ( ( 2. - x ) * QQ + x * t ), 2)) * ( 4. * ( 1. - x  + ((3. - 2*x) * QQ + t)/(QQ + x * t) * ee / 4) * (cdstar(H, Htilde) + cdstar(Htilde, H))
                - (QQ - x * (1 - 2*x) * t) / (QQ + x * t) * x * x * (cdstar(H, Etilde) + cdstar(Etilde, H) + cdstar(Htilde, E) + cdstart(E, Htilde))
                - (4. * (1 - x) * (QQ + x * t) * t + pow( QQ + t, 2 ) * ee)/() * x * (cdstar(Htilde, E) + cdstar(E, Htilde))
                - ((2. - x) * QQ + x * t) / (QQ + x * t) * ((x*x * pow( QQ + t, 2 )  / (2.*QQ * ((2 - x) * QQ + x * t)) + t / (4.*M2)) * x * (cdstar(E, Etilde) + cdstar(Etilde, E)));


    //ASK ABOUT THIS
    // c_dvcs_unp(Feff,F*)
    c_dvcs_efffs = f * c_dvcs_ffs;
    

    // dvcs c_n coefficients (BKM10 eqs. [2.20], [2.21])
    c0_dvcs_10 = 2. * lambda * bigLambda * y * (2 - y) / sqrtOnePlusEE  * c_dvcs_ffs; 
    c1_dvcs_10 = -1. * 8. * bigLambda * K / ( 2. - x ) / ( 1. + ee ) * ( -1. * lambda * y * sqrtOnePlusEE ) * c_dvcs_efffs;
    

    //Is this [2.17] we have sine coeff.
    Amp2_DVCS_10 = 1. / ( y * y * QQ ) * ( c0_dvcs_10 + c1_dvcs_10 * cos( PI - (phi * RAD) ) );

    //TO HERE

    Amp2_DVCS_10 = GeV2nb * Amp2_DVCS_10; // convertion to nb

    return dsigma_DVCS_10 = Gamma * Amp2_DVCS_10;
}
//_______________________________________________________________________________________________________________________________
Double_t TBKM::I_LP_10(Double_t *kine, Double_t phi, Double_t F1, Double_t F2, TComplex *t2cffs, TString twist = "t2") { // Interference Unpolarized Cross Section (Liuti's style)

    // Get BH propagators and set the kinematics
    BHLeptonPropagators(kine, phi);

    // Set the CFFs
    SetCFFs(t2cffs); // Etilde CFF does not appear in the interference

    // Get A_LP_I, B_LP_I and C_LP_I interference coefficients
    ABC_LP_I_10(kine, phi, A_U_I, B_U_I, C_U_I, twist);

    // BH-DVCS interference squared amplitude [at 2.34 for unpolarized?] (where is the polarized analog of this)? Does it not change?
    I_10 = 1. / ( x * y * y * y * t * P1 * P2 ) * ( A_U_I * ( F1 * H.Re() - t / 4. / M2 * F2 * E.Re() ) + B_U_I * ( F1 + F2 ) * ( H.Re() + E.Re() ) + C_U_I * ( F1 + F2 ) * Htilde.Re() );

    I_10 = GeV2nb * I_10; // convertion to nb

    return dsigma_I_10 = Gamma * I_10;
}
//_______________________________________________________________________________________________________________________________
void TBKM::ABC_LP_I_10(Double_t *kine, Double_t phi, Double_t &A_U_I, Double_t &B_U_I, Double_t &C_U_I, TString twist = "t2") { // Get A_UU_I, B_UU_I and C_UU_I interference coefficients (BKM10)

    SetKinematics(kine);

    // F_eff = f * F
    if(twist == "t2") f = 0; // F_eff = 0 ( pure twist 2)
    if(twist == "t3") f = - 2. * xi / ( 1. + xi );
    if(twist == "t3ww") f = 2. / ( 1. + xi );

    // Interference coefficients  (BKM10 Appendix A.2)
    // n = 0 -----------------------------------------
    // helicity - conserving (F) DONE
    C_110 = - 4. * lambda * bigLambda * y * (1 + sqrt(1 + ee)) / pow(1 + ee, 5./2.) * (pow(2. - y, 2) * Ktilde_10*Ktilde_10 / QQ + (1. - y + ee/4. * y*y) * (x*t/QQ - ee/2. * (1. - t / QQ)) * 1. + (sqrt(1. + ee) - 1. + 2.*x)/(1 + sqrt(1. + ee)) * t/QQ);
    C_110_V = 4. * lambda * bigLambda * y * (1 + sqrt(1 + ee)) / pow(1 + ee, 5./2.) * t/QQ * (pow(2. - y, 2) * (1. + sqrt(1. + ee) - 2.*x) / (1. + sqrt(1. + ee)) * Ktilde_10*Ktilde_10/QQ + (1. - y - ee/4.*y*y) * (2. - x + 3.*ee/2.) * (1. + (4.*(1.-x)*x + ee) / (4. - 2.*x + 3.*ee) * t/QQ) * (1. + (sqrt(1. + ee) - 1. + 2.*x) / (1. + sqrt(1. + ee)) * t/QQ)); 
    C_110_A = 4. * lambda * bigLambda * y / pow(1 + ee, 5./2.) * x*t/QQ * (2. * pow(2. - y, 2) * Ktilde_10*Ktilde_10/QQ + (1. - y - ee/4.*y*y) * (1. + sqrt(1. + ee)) * (1. - (1. - 2.*x) * t/QQ) * (1. + (sqrt(1. + ee) - 1 + 2.*x) / (1. + sqrt(1. + ee)) * t/QQ));
    // helicity - changing (F_eff) by one unit
    C_010 = 8. * sqrt(2.) * lambda * bigLambda * K * ( 1. - x ) * y * sqrt( 1. - y - ee / 4. * y * y ) / pow( 1. + ee, 2) * t / QQ;
    C_010_V = 8. * sqrt(2.) * lambda * bigLambda * K * y * sqrt( 1. - y - ee / 4. * y * y ) / pow(1. + ee, 2) * t / QQ * (x - t/QQ * (1. - 2.*x));
    C_010_A = 8. * sqrt(2.) * lambda * bigLambda * K * y * sqrt( 1. - y - ee / 4. * y * y ) / pow(1. + ee, 2) * t*x / QQ * (1. + t/QQ);
    // helicity - changing (F_eff) by two units
    C_MP0 = 4. * lambda * bigLambda * y / pow(1 + ee, 5./2.) * (Ktilde_10*Ktilde_10/QQ * pow(2. - y, 2) * (1. - sqrt(1. + ee)) + 0.5 * (1. - y - y*y*ee/4.) * (2.*x*t/QQ - (1. - t/QQ)) * (1. - sqrt(1. + ee) - t/QQ * (1. - 2.*x + sqrt(1. + ee))));
    C_MP0_V = 2. * lambda * bigLambda * y / pow(1. + ee, 5./2.) * t/QQ * ((4. - 2.*x + 3.*ee) * (1. - y - y*y*ee/4.) * (1 + t/QQ * (4. * x * (1. - x) + ee) / (4. - 2.*x + 3*ee)) * (sqrt(1. + ee) - 1. + t/QQ * (1. - 2. * x + sqrt(1. + ee))) + 2. * pow(2. - y, 2) * (sqrt(1. + ee) - 1. + 2*x) * Ktilde_10*Ktilde_10/QQ);

    // n = 1
    //tPrime is t - tmin
    C_111 = -4. * lambda * bigLambda * K * y * (2 - y) / pow(1 + ee, 5./2.) * (1 + sqrtOnePlusEE - ee) * (1 - (1 - 2*x * (2 + sqrtOnePlusEE) / (1 + sqrtOnePlusEE - ee)) * t / QQ);
    C_111_V = 8. * lambda * bigLambda * K * (2-y) * y / pow(1 + ee, 2) * (sqrtOnePlusEE + 2*(1 - x)) * t / QQ * (1 - (1 + (1 - ee) / (sqrtOnePlusEE) - 2*x * (1 + (4*(1 - x)) / (sqrtOnePlusEE))) / (2*(sqrtOnePlusEE + 2*(1 - x))) * (t - tmin) / QQ);
    C_111_A = 16. * lambda * bigLambda * K * (2 - y) * y / pow(1 + ee, 5./2.) * x * t / QQ * (1 - (1 - 2*x) * t / QQ);

    C_011 = -8. * sqrt(2.) * lambda * bigLambda * K * y * (1 - y) * sqrt(1 - y - y*y * ee / 4.) / pow(1 + ee, 2) * t / QQ;
    C_011_V = 8. * sqrt(2.) * lambda * bigLambda * y * (2 - y) * sqrt(1 - y - y*y * ee / 4.) / pow(1 + ee, 2) * t * Ktilde_10*K_tilde_10 / (QQ * QQ);
    //No C_011_A

    C_MP1 = 4. * lambda * bigLambda * K * y * (2 - y) / pow(1 + ee, 5./2.) * (1 - ee - sqrtOnePlusEE - t / QQ * (1 - ee - sqrtOnePlusEE - 2*x * (2 - sqrtOnePlusEE)));
    C_MP1_V = -4. * lambda * bigLambda * y * (2 - y) / pow(1 + ee, 5./2.) * t / QQ * (5 - 4*x + 3*ee - sqrtOnePlusEE 
            - t / QQ * (1 - ee - sqrtOnePlusEE - 2*x * (4 - 4*x - sqrtOnePlusEE)));
    C_MP1_A = -16. * lambda * bigLambda * x * y * (2 - y) / pow(1 + ee, 5./2.) * t / QQ * (1 - (1 - 2*x) * t / QQ);

    // n = 2 -----------------------------------------
    // helicity - conserving (F)
    //Wyndham done here
    C_112 = -1. * 4. * lambda * bigLambda * y * (1 - y - ee / 4 * y*y) / pow(1 + ee, 5./2.) * (x * t / QQ - (1 - t / QQ) * ee / 2)
                                                                                                            * (1 - sqrtOnePlusEE - (1 + sqrtOnePlusEE - 2*x) * t / QQ);
    C_112_V = -1. * 2. * lambda * bigLambda * y * (1 - y - ee / 4 * y*y) / pow(1 + ee, 5./2.) * (4 - 2*x + 3*ee) * t / QQ * (1 + (4 * (1 - x) * x + ee) / (4 - 2*x + 3*ee) * t / QQ)
                                                                                                            * (sqrtOnePlusEE - 1 + (1 + sqrtOnePlusEE - 2*x) * t / QQ);
    C_112_A = 4. * lambda * bigLambda * (1 - y - ee / 4 * y*y) / pow(1 + ee, 5./2.) * x * t / QQ * (1 - (1 - 2*x) * t / QQ)
                                                                                                            * (1 - sqrtOnePlusEE - (1 + sqrtOnePlusEE - 2*x) * t / QQ);
                                                                         
    // helicity - changing ONLY 1 UNIT (F_eff) [page 18]
    C_012 = -8. * sqrt(2.) * lambda * bigLambda * K * y * sqrt(1 - y - y*y * ee / 4.) / pow(1 + ee, 2) * (1 + x * t / QQ);
    
    C_012_V = 8. * sqrt(2.) * lambda * bigLambda * K * y * (1 - x) * sqrt(1 - y - y*y * ee / 4.) / pow(1 + ee, 5./2.) * t / QQ;
    
    C_012_A = 8. * sqrt(2.) * lambda * bigLambda * k * sqrt(1 - y - y*y * ee / 4) / pow(1 + ee, 2) * e * t / QQ * (1 + t / QQ);

    // helicity - changing 2 units [page 19]
    //MP meaning minus plus, might change version
    C_MP2 = -2. * lambda * bigLambda * y * (1 - y - y*y * ee / 4.) / pow(1 + ee, 5./2.) * (ee * (1 + sqrtOnePlusEE)
            - 2 * t / QQ * ((1 - x) * ee + x * (1 + sqrtOnePlusEE)) + t*t / (QQ*QQ) * (2*x + ee) * (1 - 2*x - sqrtOnePlusEE));
   
    C_MP2_V = -2. * lambda * bigLambda * y * (1 - y - y*y * ee / 4.) / pow(1 + ee, 5./2.) * t / QQ * (4 - 2*x + 3*ee + t / QQ * (4*x * (1 - x) + ee))
            * (1 + sqrtOnePlusEE - t / QQ * (1 - sqrtOnePlusEE - 2*x));
    
    C_MP2_A = -4. * lambda * bigLambda * x * y * (1 - y - y*y * ee / 4.) / pow(1 + ee, 5./2.)
            * t / QQ * (1 - (1 - 2*x) * t / QQ) * (1 + sqrtOnePlusEE - t / QQ * (1 - sqrtOnePlusEE - 2*x));
    
    //No n = 3 for polarized case as far as I can tell    

    // A_U_I, B_U_I and C_U_I
    A_U_I = C_110 + sqrt(2) / ( 2. - x ) * Ktilde_10 / sqrt(QQ) * f * C_010 + ( C_111 + sqrt(2) / ( 2. - x ) * Ktilde_10 / sqrt(QQ) * f * C_011 ) * cos( PI - (phi * RAD) )
            + ( C_112 + sqrt(2) / ( 2. - x ) * Ktilde_10 / sqrt(QQ) * f * C_012 ) * cos( 2. * ( PI - (phi * RAD) ) ) + C_113 * cos( 3. * ( PI - (phi * RAD) ) );
    B_U_I = xi / ( 1. + t / 2. / QQ ) * ( C_110_V + sqrt(2) / ( 2. - x ) * Ktilde_10 / sqrt(QQ) * f * C_010_V + ( C_111_V + sqrt(2) / ( 2. - x ) * Ktilde_10 / sqrt(QQ) * f * C_011_V ) * cos( PI - (phi * RAD) )
            + ( C_112_V + sqrt(2) / ( 2. - x ) * Ktilde_10 / sqrt(QQ) * f * C_012_V ) * cos( 2. * ( PI - (phi * RAD) ) ) + C_113_V * cos( 3. * ( PI - (phi * RAD) ) ) );
    C_U_I = xi / ( 1. + t / 2. / QQ ) * ( C_110 + C_110_A + sqrt(2) / ( 2. - x ) * Ktilde_10 / sqrt(QQ) * f * ( C_010 + C_010_A ) + ( C_111 + C_111_A + sqrt(2) / ( 2. - x ) * Ktilde_10 / sqrt(QQ) * f
            * ( C_011 + C_011_A ) ) * cos( PI - (phi * RAD) ) + ( C_112 + C_112_A + sqrt(2) / ( 2. - x ) * Ktilde_10 / sqrt(QQ) * f * ( C_012 + C_012_A ) ) * cos( 2. * ( PI - (phi * RAD) ) )
            + ( C_113 + C_113_A ) * cos( 3. * ( PI - (phi * RAD) ) ) );
}
