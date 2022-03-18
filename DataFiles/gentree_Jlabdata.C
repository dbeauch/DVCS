// *************************************************************************************************************************
// Converts the .csv datafile of all Jlab data into a .root file   																			*
// There are 195 kinematic set for the unpolarized data; 85 sets from Hall-A and 110 sets from Hall-B  							*
// Every TTree entry corresponds to one of the kimenatic sets																					*
// Data references:																																			*
// E00-110 Hall-A sets 1 - 20 (20 sets)	--> https://arxiv.org/abs/1504.05453															*
// E12-06-114 Hall-A sets 21 - 65 (45 sets) -->  https://arxiv.org/pdf/2201.03714.pdf 													*
// E07–007 HAll-A sets 176 - 195 (20 sets) --> https://arxiv.org/pdf/1703.09442.pdf														*
// Hall-B sets 66 - 175 (110 sets) -->	https://arxiv.org/pdf/1504.02009.pdf									  							*
// *************************************************************************************************************************

#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include "TList.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"

using namespace std;

//_________________________________________________________________________________________________________
void gentree_Jlabdata() {

	const Int_t kMaxNumOfPts = 24;
	Int_t npoints;

	// File to save the TTree
	TFile fout("dvcs_Jlabdata.root","RECREATE");

	TCanvas *c1;

	struct kin_t {
		Double_t k;
		Double_t QQ;
		Double_t xB;
		Double_t t;
   };
	kin_t kin;
	Double_t phi[kMaxNumOfPts];
	Double_t F[kMaxNumOfPts], var_F[kMaxNumOfPts];
	Double_t sys_high_F[kMaxNumOfPts], sys_low_F[kMaxNumOfPts], sys_corr_F[kMaxNumOfPts];

	TTree *t3 = new TTree("Jlab_dvcs_unpol","dvcs");
	t3->Branch("kinematics",&kin.k,"k/D:QQ:xB:t");
	t3->Branch("npoints",&npoints,"npoints/I");
	t3->Branch("phi",phi,"phi[npoints]/D");
	t3->Branch("F",F,"F[npoints]/D");
	t3->Branch("var_F",var_F,"var_F[npoints]/D");

   TTree *T1 = new TTree("T1","tree data from ascii file");
   Long64_t nlines_T1 = T1->ReadFile("AllJlabData.csv","f_set/I:f_index:f_k/D:f_QQ:f_xB:f_t:f_phi:f_F:f_sigmaF:f_varF:f_F1:f_F2");
   printf(" found %lld points\n",nlines_T1);

	TTree *T2 = new TTree("T2","tree data from ascii file");
   Long64_t nlines_T2 = T2->ReadFile("AllJlabData.csv","f_set/I:f_index:f_k/D:f_QQ:f_xB:f_t:f_phi:f_F:f_sigmaF:f_varF:f_F1:f_F2");
   //T2->Write();

	Int_t f_set_T1, f_set_T2, f_index;
	Double_t  f_k, f_QQ, f_xB, f_t, f_phi, f_F, f_sigmaF, f_varF, f_F1, f_F2;

   T1->SetBranchAddress("f_set",&f_set_T1);
	T2->SetBranchAddress("f_set",&f_set_T2);
	T1->SetBranchAddress("f_index",&f_index);
	T1->SetBranchAddress("f_k",&f_k);
	T1->SetBranchAddress("f_QQ",&f_QQ);
	T1->SetBranchAddress("f_xB",&f_xB);
	T1->SetBranchAddress("f_t",&f_t);
	T1->SetBranchAddress("f_phi",&f_phi);
	T1->SetBranchAddress("f_F",&f_F);
	T1->SetBranchAddress("f_sigmaF",&f_sigmaF);

	Long_t nentries = T1->GetEntries();

	//Long_t lOneTenthOfNData = ((double)(nentries) / 10. );

	Int_t iset = 0, ipt = 0;

	for(Long_t i = 0; i<nentries; i++) { // tree entry loop

		T1->GetEntry(i);
		//if( icand % lOneTenthOfNData == 0 )
		//cout<<" Currently at data point........: "<<icand<<" / "<<lNDataPoints<<" ( "<<(long)(((double)(icand)/(double)(lNDataPoints))*(100.+1e-3))<<"% )"<<endl;
		if (i != nentries-1)  T2->GetEntry(i+1);
		else f_set_T2 = -1; // dummy value to be able to select the last set of the file

		cout<<"entry "<<i<<": "<<"set T1 = "<<f_set_T1<<", set T2 = "<<f_set_T2<<endl;
		cout<<"index "<<f_index<<", QQ = "<<f_QQ<<", xB = "<<f_xB<<", t = "<<f_t<<endl;
		cout<<"phi = "<<f_phi<<", F = "<<f_F<<", var_F = "<<f_sigmaF<<endl;

		if ( (f_F == 0.) &&  (f_set_T1 == f_set_T2) ) continue;		// if it is not the last point of the set, then continue
		if ( (f_F == 0.) &&  (f_set_T1 != f_set_T2) ) {					// if it is the last point of the set, then fill the tree with last values where F was not zero and continue

			cout<<"=============== Set "<<iset<<" ======================="<<endl;
			cout<<"k = "<< f_k<<endl;
			cout<<"QQ = "<< f_QQ<<endl;
			cout<<"xB = "<< f_xB<<endl;
			cout<<"t = "<< f_t<<endl;

			kin.k =  f_k;
			kin.QQ = f_QQ;
			kin.xB = f_xB;
			kin.t = f_t;

			iset++;

			npoints = ipt;

			ipt = 0;

			t3 ->Fill();

			continue;
		}

		phi[ipt] = f_phi;
		F[ipt] = f_F;
		var_F[ipt] = f_sigmaF;

		ipt++;

		if( f_set_T1 != f_set_T2 ){

			cout<<"=============== Set "<<iset<<" ======================="<<endl;
			cout<<"k = "<< f_k<<endl;
			cout<<"QQ = "<< f_QQ<<endl;
			cout<<"xB = "<< f_xB<<endl;
			cout<<"t = "<< f_t<<endl;

			kin.k =  f_k;
			kin.QQ = f_QQ;
			kin.xB = f_xB;
			kin.t = f_t;

			iset++;

			npoints = ipt;
			cout<<"npoints = "<<npoints<<endl;
			ipt = 0;

			t3 ->Fill();
		}
	}

	t3->Print();
	t3->Show(0);
	t3->Show(1);
	t3->Write();

}
