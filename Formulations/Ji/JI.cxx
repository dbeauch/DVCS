/* 
    Ji Formulation Implementation
    Made By Wyndham White wrw2ztk@virginia.edu

*/

#include "JI.h"

using namespace std;  	// std namespace: so you can do things like 'cout'

ClassImp(JI)			// classimp: necessary for root


//_______________________________________________________________________________________________________________________________
JI::JI() {
	// Default Constructor
}
//_______________________________________________________________________________________________________________________________
JI::~JI() {
	// Default Destructor
}
//_______________________________________________________________________________________________________________________________
TComplex TBKM::cdstar( TComplex c, TComplex d ){ // ( C D* ) product

    TComplex dstar = TComplex::Conjugate(d);

    return ( c.Re() * dstar.Re() - c.Im() * dstar.Im() ) + ( c.Re() * dstar.Im() + c.Im() * dstar.Re() ) * TComplex::I();
}
//_______________________________________________________________________________________________________________________________
void JI::SetCFFs( TComplex *t2cffs ) { // Twist-2 Compton Form Factors

     H = t2cffs[0];
     E = t2cffs[1];
     Htilde = t2cffs[2];
     Etilde = t2cffs[3];
}
//_______________________________________________________________________________________________________________________________
void JI::SetKinematics( Double_t *kine ) {

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
    Ktilde_10 = sqrt( tmin - t ) * sqrt( ( 1. - x ) * sqrt( 1. + ee ) + ( ( t - tmin ) * ( ee + 4. * x * ( 1. - x ) ) / 4. / QQ ) ) * sqrt( 1. - y - y * y * ee / 4. )
                / sqrt( 1. - y + y * y * ee / 4.); // K tilde from 2010 paper
    K = sqrt( 1. - y + y * y * ee / 4.) * Ktilde_10 / sqrt(QQ);
    lambda = 0.;
    bigLambda = 0.;
}

// EQUATION 3.17 - not finished
Double_t JI::DVCS_squaredAmplitude(){
    Double_t L_DVCS = 1.0;
    Double_t H_DVCS = H_DVCS_Calculation();
    return L_DVCS * H_DVCS / QQ / QQ;
}

// EQUATION 3.29 - not finished 
Double_t JI:H_DVCS_Calculation() {
    Double_t H_U = H_DVCS_U_Calculation();

    //Not relevant for Lambda_L = 0
    Double_t H_L;

    //Neither of these relevant at the moment because considering Lambda_T = 0
    Double_t H_T_in;
    Double_t H_T_out;


    //Other terms need to be added into this
    return H_U;
}

// EQUATION 3.34 - not finished
Double_t JI:: H_DVCS_U_Calculation() {

    //NOT sure what script H is 
    Double_t scriptH = 1.0;
    Double_t scriptHTilde = 1.0;
    return  4.0 * ( (1 - xi*xi) * (scriptH*cdstar(H, H) + scriptHTilde * cdstar(Htilde, Htilde)) 
                    - t / 4*M2 * (scriptH * cdstar(E, E) + xi*xi * scriptHTilde * cdstar(Etilde, Etilde))
                    - xi*xi * (scriptH * cdstar(E, E) + scriptH * (cdstar(H, E) + cdstair(E, H)) + scriptHTilde(cdstar(Htilde, Etilde) + cdstar(Etilde, Htilde))) );

}



