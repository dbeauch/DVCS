/* 
    Ji Formulation Implementation
    Made By Wyndham White wrw2ztk@virginia.edu

*/

#include "JI.h"

using namespace std; 	// std namespace: so you can do things like 'cout'

//ClassImp(JI);			// classimp: necessary for root


//_______________________________________________________________________________________________________________________________
JI::JI() {
	// Default Constructor
}
//_______________________________________________________________________________________________________________________________
JI::~JI() {
	// Default Destructor
}
//_______________________________________________________________________________________________________________________________

//_______________________________________________________________________________________________________________________________
void JI::SetCFFs( complexDouble *t2cffs ) { // Twist-2 Compton Form Factors

     H = t2cffs[0];
     E = t2cffs[1];
     Htilde = t2cffs[2];
     Etilde = t2cffs[3];
}
//_______________________________________________________________________________________________________________________________
void JI::SetKinematics( double *kine ) {

    QQ = kine[0];     //Q^2 value
    q = sqrt(-1.0*QQ);
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

void JI::setAlphaBeta(double newAlpha, double newBeta){
    alpha = newAlpha;
    beta = newBeta;
}

void JI::setValuesOf_h(double phi, int twist){
    
    //convert phi to radians
    phi *= RAD;
    
    //h_+U always 0 regardless of twist
    h_plusU = 0.0;

    for(int i = 2; i <= twist; i++){
        switch(i) {
            case 2:
                h_U = (2 - 2*y + y*y) 
                        / (2*y*y);
                h_minusL = (y - 2) 
                            / (2*y);
                continue;
            case 3:
                h_U += 2*x * (y - 2) * sqrt((t * (x - 1) - M2 * x*x) * (1 - y)) * (beta - 1) * cos(phi) 
                                        / ((1 + (beta - 1) * x) * (y*y));
                h_minusL += 2*x * sqrt((t * (x - 1) - M2 * x*x) * (1 - y)) * (beta - 1) * cos(phi) 
                                            / ((1 + (beta - 1) * x) * y);
                continue;
            case 4:
                h_U +=                   1 
                    / (2 * pow(1 + (beta - 1) * x, 2) * y*y)
                    * (8*pow(beta - 1, 2) * x*x * (y - 1) * pow(cos(phi), 2) * (M2 * x*x - t * x + t)
                    - x*x * pow(y - 2, 2) * (M2 * (pow(2*(beta - 1) * x + 1, 2) + 1) - 2*pow(beta - 1, 2) * t * (x - 1)));
                h_minusL += x*x * (y - 2) * (M2 * (2*pow(beta - 1, 2) * x*x + 2*(beta - 1) * x + 1) - pow(beta - 1, 2) * t * (x - 1))
                                                / (y * pow((beta - 1) * x + 1, 2));
                continue;
        }
    }
    //h_U~ always h_U regardless of twist
    double h_U_tilde = h_U;

    return;

    /* THIS IS HOW IT USED BE STRUCTURED. LEAVING HERE
    IN CASE WE NEED TO COMBINE FOR SOME REASON
    double h_U_twist2 = (2 - 2*y + y*y) 
                           / (2*y*y);
    double h_U_tilde_twist2 = h_U_twist2;

    double h_U_twist3 = 2*x * (y - 2) * sqrt((t * (x - 1) - M2 * x*x) * (1 - y)) * (beta - 1) * cos(phi) 
                                            / ((1 + (beta - 1) * x) * (y*y));
    double h_U_tilde_twist3 = h_U_twist3;

    double h_U_twist4 =           1 
                / (2 * pow(1 + (beta - 1) * x, 2) * y*y)
                * (8*pow(beta - 1, 2) * x*x * (y - 1) * pow(cos(phi), 2) * (M2 * x*x - t * x + t)
                - x*x * pow(y - 2, 2) * (M2 * (pow(2*(beta - 1) * x + 1, 2) + 1) - 2*pow(beta - 1, 2) * t * (x - 1)));
    double h_U_tilde_twist4 = h_U_twist4;

    //h -L
    double h_minusL_twist2 = (y - 2) 
                              / (2*y);
    double h_minusL_twist3 = 2*x * sqrt((t * (x - 1) - M2 * x*x) * (1 - y)) * (beta - 1) * cos(phi) 
                                                / ((1 + (beta - 1) * x) * y);
    double h_minusL_twist4 = x*x * (y - 2) * (M2 * (2*pow(beta - 1, 2) * x*x + 2*(beta - 1) * x + 1) - pow(beta - 1, 2) * t * (x - 1))
                                                        / (y * pow((beta - 1) * x + 1, 2));
    //h +U
    double h_plusU_twist2 = 0.0;
    double h_plusU_twist3 = 0.0;
    double h_plusU_twist4 = 0.0;
    */
}

double JI::get_h_U(){
    return h_U;
}


// ALL OLD: PRIOR TO SWITCHING TO DIFFERENT MODEL OF SAME FORMULATION

// //Using 5.8
// void JI::set_h_values() {

//     h_U = 2 * (k * q) * (k * qPrime) / (q * qPrime) - (q*q) * pow((k * qPrime), 2) / pow(q * qPrime, 2) - k * q;
//     h_tilde_U = h_U;
//     h_minus_L = -1.0 * q*q * (k * qPrime) / (q * qPrime) + (k * q);
//     h_plus_U = 0;

// }

// // EQUATION 3.17 - not finished - stores squared amplitude in Amp2_DVCS
// void JI::DVCS_squaredAmplitude(){
//     Double_t L_DVCS = 1.0;
//     Double_t H_DVCS = H_DVCS_Calculation();
//     Amp2_DVCS = L_DVCS * H_DVCS / QQ / QQ;
// }

// // EQUATION 3.29 - not finished 
// void JI:H_DVCS_Calculation() {
//     Double_t H_U = H_DVCS_U_Calculation();

//     //Not relevant for Lambda_L = 0
//     H_L_DVCS = 0.0;

//     //Neither of these relevant at the moment because considering Lambda_T = 0
//     H_T_in_DVCS = 0.0;
//     H_T_out_DVCS = 0.0;


//     //Other terms need to be added into this
//     H_DVCS = H_U;
// }

// // EQUATION 3.34 - not finished
// void JI:: H_DVCS_U_Calculation() {

//     //NOT sure what script H is 
//     Double_t scriptH = 1.0;
//     Double_t scriptHTilde = 1.0;
//     return  4.0 * ( (1 - xi*xi) * (scriptH*cdstar(H, H) + scriptHTilde * cdstar(Htilde, Htilde)) 
//                     - t / 4*M2 * (scriptH * cdstar(E, E) + xi*xi * scriptHTilde * cdstar(Etilde, Etilde))
//                     - xi*xi * (scriptH * cdstar(E, E) + scriptH * (cdstar(H, E) + cdstair(E, H)) + scriptHTilde(cdstar(Htilde, Etilde) + cdstar(Etilde, Htilde))) );

// }

