#ifndef TBKM_H
#define TBKM_H


class TBKM {


private:

	Double_t ALP_INV = 137.0359998; // 1 / Electromagnetic Fine Structure Constant
	Double_t PI = TMath::Pi();
	Double_t RAD = PI / 180.;
	Double_t M = 0.938272; //Mass of the proton in GeV
	Double_t M2 = M*M; //Mass of the proton  squared in GeV
	Double_t GeV2nb = .389379*1000000; // Conversion from 1/GeV2 to NanoBarn

	Double_t QQ, x, t, k; // kinematics
	Double_t ee, y, xi, tmin, s, Gamma;
	Double_t  K, Ktilde_10, KD;
	Double_t P1, P2; // lepton propagators
	TComplex H, E, Htilde, Etilde; // Twist-2 CFFs

	Double_t c0_BH, c1_BH, c2_BH; // BH unpolarized c coefficients (BKM02 eqs. [35, 37])
	Double_t c0_dvcs, c1_dvcs; // DVCS unpolarized Fourier harmonics (BKM02 eqs. [43, 44])
	Double_t Amp2_BH, Amp2_DVCS, I; // squared amplitudes
	Double_t dsigma_BH, dsigma_DVCS, dsigma_I; // 4-fold differential cross sections

public:

	TBKM();  // Constructor
	~TBKM(); // Destructor

	TComplex cdstar( TComplex c, TComplex d ); // complex and complex conjugate numbers product
	void SetCFFs( TComplex *t2cffs ); // t2cffs = { H, E , Htilde, Etilde } Twist-2 Compton Form Factors
	void SetKinematics(Double_t *kine);
	void BHLeptonPropagators(Double_t *kine, Double_t phi);
	//BKM02
	Double_t BH_UU(Double_t *kine, Double_t phi, Double_t F1, Double_t F2);
	Double_t DVCS_UU_02(Double_t *kine, Double_t phi, TComplex *t2cffs, TString twist);
	Double_t I_UU_02(Double_t *kine, Double_t phi, Double_t F1, Double_t F2, TComplex *t2cffs, TString twist);
	//BKM10


	ClassDef(TBKM,1);


};

#endif
