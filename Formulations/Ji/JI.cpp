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
//NOTE: TAKEN FROM BKM SO ANY DIFFERENT COEFFICIENTS COULD BE DIFFERENT
void JI::SetKinematics( double *kine ) {

    QQ = kine[0];     //Q^2 value
    q = sqrt(-1.0*QQ);
    x = kine[1];      //Bjorken x
    t = kine[2];      //momentum transfer squared
    k = kine[3];      //Electron Beam energy

    ee = 4. * M2 * x * x / QQ; // epsilon squared
    y = sqrt(QQ) / ( sqrt(ee) * k );  // lepton energy fraction
    //DIFFERENT XI FOR JI
    // Generalized Bjorken variable
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

void JI::updateXi(){
    xi = (x / (2 - x)) - ( (2*x * (2 * alpha * M2 * x * x + (x - 1) * (1 - 2*alpha + (beta - 1) * x) * t)) 
                                    / (pow(x - 2, 2) * (1 + (beta - 1) * x) * QQ) );
    N = sqrt(-4.0 * M2 * xi*xi - t * (1 - xi*xi)) / M;
}

void JI::setValuesOf_h(double phi, int twist){
    
    //convert phi to radians
    phi *= RAD;
    
    //h_+U always 0 regardless of twist
    h_plusU = 0.0;

    for(int i = 2; i <= twist; i++){
        switch(i) {
            case 2:
                h_U = QQ * (2 - 2*y + y*y) 
                        / (2*y*y);
                h_minusL = (y - 2) 
                            / (2*y);
                continue;
            case 3:
                h_U += sqrt(QQ) * (2*x * (y - 2) * sqrt((t * (x - 1) - M2 * x*x) * (1 - y)) * (beta - 1) * cos(phi) 
                                        / ((1 + (beta - 1) * x) * (y*y)));
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
    
}

//B.1-B.6
void setValuesOf_F(complexDouble H, complexDouble H_tilde, complexDouble E, complexDouble E_tilde){
    //F_UU
    complexDouble firstTerm = (cdstar(H, H) * h_U) + (cdstar(H_tilde, H_tilde) * h_U_tilde);
        firstTerm = firstTerm * (1 - xi*xi);
    complexDouble secondTerm = (cdstar(E, E) * h_U) + (cdstar(E_tilde, E_tilde) * (xi*xi * h_U_tilde));
        secondTerm = secondTerm * (-1.0 * t / (4*M2));
    complexDouble thirdTerm = (cdstar(E, E) * h_U) + ((cdstar(H, E) + cdstar(E, H)) * h_U) + (cdstar(H_tilde, E_tilde) + cdstar(E_tilde, H_tilde));
        thirdTerm = thirdTerm * (-1.0 * xi*xi);
    complexDouble sumTerm;

    F_UU = firstTerm + secondTerm + thirdTerm;
    F_UU = F_UU * 4;

    //F_UT_out
    firstTerm = (cdstar(E, H) * h_U);
    secondTerm = cdstar(E, H_tilde) * (-1.0 * h_U_tilde * xi);
    sumTerm = firstTerm + secondTerm;

    F_UT_out.imaginaryPart = N * 4 * sumTerm.imaginaryPart;
    F_UT_out.realPart = 0.0;

    //F_UL & F_LL
    firstTerm = cdstar(H, H_tilde) * (1 - xi*xi);
    secondTerm = (cdstar(E, H_tilde) + cdstar(H, E_tilde)) * (-1.0 * xi*xi);
    thirdTerm =  cdstar(E, E_tilde) * (-1.0 * xi * ((xi * xi / (1 + xi)) + (t / (4*M2))));
    sumTerm = firstTerm + secondTerm + thirdTerm;
    F_UL.imaginaryPart = sumTerm.imaginaryPart * 8 * h_plusU;
    F_UL.realPart = 0.0;
    F_LL.realPart = sumTerm.realPart * 8 * h_minusL;
    F_LL.imaginaryPart = 0.0;

    //F_UT & F_LL
    firstTerm = cdstar(E, H_tilde);
    secondTerm = cdstar(H, E_tilde) * (-1.0 * xi);
    thirdTerm = cdstar(E, E_tilde) * (-1.0 * xi * xi / (1 + xi));
    sumTerm = firstTerm + secondTerm + thirdTerm;
    F_UT_in.imaginaryPart = sumTerm.imaginaryPart * (4 * N * h_plusU);
    F_UT_in.realPart = 0;
    F_LT_in.imaginaryPart = 0;
    F_LT_in.realPart = sumTerm.realPart * (4 * N * h_minusL);
}

double JI::get_h_U(){
    return h_U;
}

double JI::get_h_U_tilde(){
    return h_U_tilde;
}

double JI::get_h_minusL(){
    return h_minusL;
}

double JI::get_h_plusU(){
    return h_plusU;
}


