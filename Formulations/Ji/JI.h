#ifndef JI_H
#define JI_H


class JI {


private:

	Double_t ALP_INV = 137.0359998; // 1 / Electromagnetic Fine Structure Constant
	Double_t PI = TMath::Pi();
	Double_t RAD = PI / 180.;
	Double_t M = 0.938272; //Mass of the proton in GeV
	Double_t M2 = M*M; //Mass of the proton  squared in GeV
	Double_t GeV2nb = .389379*1000000; // Conversion from 1/GeV2 to NanoBarn

	Double_t QQ, x, t, k; // kinematics
	Double_t ee, y, xi, tmin, s, Gamma, sqrtOnePlusEE;
	Double_t  K, Ktilde_10, KD;
	Double_t lambda, bigLambda;
	Double_t P1, P2; // lepton propagators
	TComplex H, E, Htilde, Etilde; // Twist-2 CFFs

	Double_t A_02, B_02, C_02; // A_UU_I, B_UU_I, C_UU_I interference coefficients in BKM02
	Double_t Amp2_DVCS,; // squared amplitudes BKM02

public:

	JI();  // Constructor
	~JI(); // Destructor

	TComplex cdstar( TComplex c, TComplex d ); // complex and complex conjugate numbers product
	void SetCFFs( TComplex *t2cffs ); // t2cffs = { H, E , Htilde, Etilde } Twist-2 Compton Form Factors
	void SetKinematics(Double_t *kine);
	void BHLeptonPropagators(Double_t *kine, Double_t phi);

	Double_t DVCS_squaredAmplitude(); 	// EQUATION 3.17 - not finished
	Double_t H_DVCS_Calculation();    	// EQUATION 3.29 - not finished 
	Double_t H_DVCS_U_Calculation()		// EQUATION 3.34 - not finished

	ClassDef(JI,1);


};

#endif
