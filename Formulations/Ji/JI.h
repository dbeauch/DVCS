#ifndef JI_H
#define JI_H

#include <cmath>
#include "complexDouble.hpp"

class JI {


private:

	double ALP_INV = 137.0359998; // 1 / Electromagnetic Fine Structure Constant
	double PI = M_PI;
	double RAD = M_PI / 180.;
	double M = 0.938272; //Mass of the proton in GeV
	double M2 = M*M; //Mass of the proton  squared in GeV
	double GeV2nb = .389379*1000000; // Conversion from 1/GeV2 to NanoBarn

	double QQ, x, t, k, q; // kinematics
	double ee, y, xi, tmin, s, Gamma, sqrtOnePlusEE;
	double qPrime, alpha, beta;
	double K, Ktilde_10, KD;
	double lambda, bigLambda;
	double P1, P2; // lepton propagators
	complexDouble H, E, Htilde, Etilde; // Twist-2 CFFs

	double Amp2_DVCS; //Eqn. 3.17
	//Double_t H_DVCS, H_U_DVCS, H_L_DVCS, H_T_in_DVCS, H_T_out_DVCS; //3.29 OLD MODEL
	double h_U, h_U_tilde, h_plusU, h_minusL; // B.13-B.15

public:

	JI();  // Constructor
	~JI(); // Destructor

	void SetCFFs( complexDouble *t2cffs ); // t2cffs = { H, E , Htilde, Etilde } Twist-2 Compton Form Factors
	void SetKinematics(double *kine);
	void setAlphaBeta(double newAlpha, double newBeta);
	//void BHLeptonPropagators(Double_t *kine, Double_t phi);
	
	void setValuesOf_h(double phi, int twist); //PHI IN DEGREES FOR THIS
	double get_h_U();

	// WORK FROM BEFORE SWITCHING MODEL
	//void set_h_values();
	//void DVCS_squaredAmplitude(); 	// EQUATION 3.17 - not finished
	//void H_DVCS_Calculation();    	// EQUATION 3.29 - not finished 
	//void H_DVCS_U_Calculation()		// EQUATION 3.34 - not finished

	//ClassDef(JI,1);


};

#endif
