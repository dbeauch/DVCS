#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include "Fit/Fitter.h"
#include "Fit/BinData.h"
#include "Fit/Chi2FCN.h"
#include "TList.h"
#include "Math/WrappedMultiTF1.h"
#include "HFitInterface.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "/media/lily/Data/GPDs/DVCS/Formulations/BKM/TBKM.h"
#include "/media/lily/Data/GPDs/DVCS/Formulations/BKM/TBKM.cxx"
#include "/media/lily/Data/GPDs/DVCS/FFs/TFormFactors.cxx"
#include "/media/lily/Data/GPDs/DVCS/FFs/TFormFactors.h"
#include "/media/lily/Data/GPDs/DVCS/GPD_Models/TGPDModels.h"

using namespace std;

TBKM *bkm = new TBKM;

TFormFactors *ff = new TFormFactors;

// user defined functions
Double_t BH(Double_t *angle, Double_t *par);
Double_t TotalUUXS_BKM02(Double_t *angle, Double_t *par);
Double_t TotalUUXS_BKM10(Double_t *angle, Double_t *par);
Double_t TotalUUXS_BKM02_Fit(Double_t *angle, Double_t *par);
Double_t TotalUUXS_BKM10_Fit(Double_t *angle, Double_t *par);
void DoROOTFit(TGraphErrors* gGenData, TF1* ffit );

void genpseudoBKM( Bool_t IsModel = 1, TString kFormulation = "BKM10", TString kKinematics = "Jlab_all", TString kGPDModel = "KM15" )
{
   Int_t seed = 400;
   Double_t var = 0.05;
   const Int_t kMaxNumOfDataPts = 24;

   // Get Jlab kinematics and phi values from the data file
   TFile* dataFile = TFile::Open("/media/lily/Data/GPDs/DataFiles/dvcs_Jlabdata.root");

   // Create a TTreeReader for the tree, by passing the TTree's name and the  TFile it is in.
   TTreeReader dataReader("Jlab_dvcs_unpol", dataFile);
   // read data tree
   TTreeReaderValue<Double_t> k_jlab(dataReader, "kinematics.k");
   TTreeReaderValue<Double_t> QQ_jlab(dataReader, "kinematics.QQ");
   TTreeReaderValue<Double_t> xB_jlab(dataReader, "kinematics.xB");
   TTreeReaderValue<Double_t> t_jlab(dataReader, "kinematics.t");
   TTreeReaderArray<Double_t> phi_jlab(dataReader, "phi");
   TTreeReaderValue<Int_t> npts_jlab(dataReader, "npoints");

   // Generated CFFs 2D graphs vs t and xB
   TGraph2D *g2D_ReH = new TGraph2D();
      g2D_ReH ->SetName("ReH");
      g2D_ReH ->SetTitle("; t [GeV^{2}]; x_{B};ReH");
      g2D_ReH ->SetMarkerSize(2);
      g2D_ReH ->SetMarkerStyle(22);
      g2D_ReH ->SetLineColor(1);
      g2D_ReH ->SetLineWidth(1);
   TGraph2D *g2D_ReE = new TGraph2D();
      g2D_ReE ->SetName("ReE");
      g2D_ReE ->SetTitle("; t [GeV^{2}]; x_{B};ReE");
      g2D_ReE ->SetMarkerSize(2);
      g2D_ReE ->SetMarkerStyle(22);
      g2D_ReE ->SetLineColor(1);
      g2D_ReE ->SetLineWidth(1);
   TGraph2D *g2D_ReHtilde = new TGraph2D();
      g2D_ReHtilde ->SetName("ReHtilde");
      g2D_ReHtilde ->SetTitle("; t [GeV^{2}]; x_{B};ReHtilde");
      g2D_ReHtilde ->SetMarkerSize(2);
      g2D_ReHtilde ->SetMarkerStyle(22);
      g2D_ReHtilde ->SetLineColor(1);
      g2D_ReHtilde ->SetLineWidth(1);
   TGraph2D *g2D_dvcs = new TGraph2D();
      g2D_dvcs ->SetName("dvcs");
      g2D_dvcs ->SetTitle("; t [GeV^{2}]; x_{B};dvcs");
      g2D_dvcs ->SetMarkerSize(2);
      g2D_dvcs ->SetMarkerStyle(22);
      g2D_dvcs ->SetLineColor(1);
      g2D_dvcs ->SetLineWidth(1);

   // ----- CFFs results -----
   TGraphErrors* gReH_true = new TGraphErrors(); // Generated ReH
   TGraphErrors* gReE_true = new TGraphErrors(); // Generated ReE
   TGraphErrors* gReHt_true = new TGraphErrors(); // Generated ReHt
   TGraphErrors* gdvcs_true = new TGraphErrors(); // Generated ReH
   TGraphErrors *gReH_fit = new TGraphErrors();
   TGraphErrors *gReE_fit = new TGraphErrors();
   TGraphErrors *gReHt_fit = new TGraphErrors();
   TGraphErrors *gdvcs_fit = new TGraphErrors();

   gReH_fit ->SetName("phi_fit");
   gReH_fit ->SetMarkerColor(kBlue);
   gReH_fit ->SetLineColor(kBlue);
   gReH_fit ->SetMarkerStyle(22);
   gReH_fit ->SetMarkerSize(1.3);

   gReE_fit ->SetName("phi_fit");
   gReE_fit ->SetMarkerColor(kBlue);
   gReE_fit ->SetLineColor(kBlue);
   gReE_fit ->SetMarkerStyle(22);
   gReE_fit ->SetMarkerSize(1.3);

   gReHt_fit ->SetName("phi_fit");
   gReHt_fit ->SetMarkerColor(kBlue);
   gReHt_fit ->SetLineColor(kBlue);
   gReHt_fit ->SetMarkerStyle(22);
   gReHt_fit ->SetMarkerSize(1.3);

   gdvcs_fit ->SetName("phi_fit");
   gdvcs_fit ->SetMarkerColor(kBlue);
   gdvcs_fit ->SetLineColor(kBlue);
   gdvcs_fit ->SetMarkerStyle(22);
   gdvcs_fit ->SetMarkerSize(1.3);

   // Read file containing CFFs from ANN Global Fit
   TTree *tree_globalCFFs = new TTree("tree_globalCFFs","tree data from ascii file");
   Double_t  f_global_ReH, f_global_ReE, f_global_ReHtilde, f_global_dvcs;

   if ( !IsModel ) {
      if(kFormulation == "BKM02")
         tree_globalCFFs ->ReadFile("./ANN_GlobalFit_CFFs/BKM02_ModelFromData.txt","f_set/I:f_k/D:f_QQ:f_xB:f_t:f_global_ReH:f_global_ReE:f_global_ReHtilde:f_global_dvcs");

      if(kFormulation == "BKM10")
         tree_globalCFFs ->ReadFile("./ANN_GlobalFit_CFFs/BKM10_ModelFromData.txt","f_set/I:f_k/D:f_QQ:f_xB:f_t:f_global_ReH:f_global_ReE:f_global_ReHtilde:f_global_dvcs");

      tree_globalCFFs->SetBranchAddress("f_global_ReH",&f_global_ReH);
      tree_globalCFFs->SetBranchAddress("f_global_ReE",&f_global_ReE);
      tree_globalCFFs->SetBranchAddress("f_global_ReHtilde",&f_global_ReHtilde);
      tree_globalCFFs->SetBranchAddress("f_global_dvcs",&f_global_dvcs);
   }

   // Generating function initialization
   TF1 *fgen;
   if ( IsModel ) {
      if(kFormulation == "BKM02")   fgen = new TF1("fgen", TotalUUXS_BKM02, 0, 360, 12);
      if(kFormulation == "BKM10")   fgen = new TF1("fgen", TotalUUXS_BKM10, 0, 360, 12);
   }
   if ( !IsModel ) {
      if(kFormulation == "BKM02")   fgen = new TF1("fgen", TotalUUXS_BKM02_Fit, 0, 360, 8);
      if(kFormulation == "BKM10")   fgen = new TF1("fgen", TotalUUXS_BKM10_Fit, 0, 360, 8);
   }
   // Fit function initialization
   TF1* ffit;
      if(kFormulation == "BKM02")   ffit = new TF1("ffit", TotalUUXS_BKM02_Fit, 0, 360, 8);
      if(kFormulation == "BKM10")   ffit = new TF1("ffit", TotalUUXS_BKM10_Fit, 0, 360, 8);
   ffit ->SetParNames("k", "QQ", "xB", "t", "ReH", "ReE", "ReHtilde", "dvcs");
   // BH contribution
   TF1 *fBH = new TF1("fBH", BH, 0, 360, 4);
   fBH ->SetParNames("k", "QQ", "xB", "t");

   // Pseudo-data cross section
   TGraphErrors* gGenDVCS; // Total
   TGraphErrors* gGenIntfDVCS;   // F - BH
   TGraphErrors *gSelPseudoF[200];
   TGraphErrors *gSelPseudoIntfDVCS[200];

   TCanvas *c_F;
   TMultiGraph *mgr_pseudo[200];

   ofstream myfile; //output file
   if(IsModel) myfile.open (Form("./pseudo_%s_%s_%s_t2_5pct.csv", kGPDModel.Data(), kFormulation.Data(), kKinematics.Data()));
   if(!IsModel) myfile.open (Form("./pseudo_ANN-global_%s_%s_t2_5pct.csv", kFormulation.Data(), kKinematics.Data()));
   myfile<<"#Set,index,k,QQ,x_b,t,phi_x,F,sigmaF,varF,ReH,ReE,ReHTilde,dvcs"<<endl;

   // File to save the TTree
   TString outputFileName;
   if(IsModel) outputFileName = Form("./pseudo_%s_%s_%s_t2_5pct.root", kGPDModel.Data(), kFormulation.Data(), kKinematics.Data());
   if(!IsModel) outputFileName = Form("./pseudo_ANN-global_%s_%s_t2_5pct.root", kFormulation.Data(), kKinematics.Data());

   TFile fout(outputFileName.Data(),"RECREATE");

   struct kin_t {
      Double_t k;
      Double_t QQ;
      Double_t xB;
      Double_t t;
   };
   kin_t kin;
   Int_t npoints;
   Double_t phi[kMaxNumOfDataPts];
   Double_t F[kMaxNumOfDataPts], errF[kMaxNumOfDataPts];
   Double_t varF;
   Double_t gdvcs, dvcs, e_dvcs;
   Double_t gReH, gReE, gReHtilde, gReEtilde, ReH, ReE, ReHtilde;
   Double_t e_ReH, e_ReE, e_ReHtilde;
   Double_t gImH, gImE, gImHtilde, gImEtilde;

   TTree *t3 = new TTree("dvcs","generated dvcs");
   t3->Branch("kinematics",&kin.k,"k/D:QQ:xB:t");
   t3->Branch("npoints",&npoints,"npoints/I");
   t3->Branch("phi",phi,"phi[npoints]/D");
   t3->Branch("F",F,"F[npoints]/D");
   t3->Branch("errF",errF,"errF[npoints]/D");
   t3->Branch("varF",&varF,"varF/D");
   t3->Branch("gReH",&gReH,"gReH/D");
   t3->Branch("gReE",&gReE,"gReE/D");
   t3->Branch("gReHtilde",&gReHtilde,"gReHtilde/D");
   t3->Branch("gdvcs",&gdvcs,"gdvcs/D");
   t3->Branch("ReH",&ReH,"ReH/D");
   t3->Branch("ReE",&ReE,"ReE/D");
   t3->Branch("ReHtilde",&ReHtilde,"ReHtilde/D");
   t3->Branch("e_ReH",&e_ReH,"e_ReH/D");
   t3->Branch("e_ReE",&e_ReE,"e_ReE/D");
   t3->Branch("e_ReHtilde",&e_ReHtilde,"e_ReHtilde/D");
   t3->Branch("dvcs",&dvcs,"dvcs/D");
   t3->Branch("e_dvcs",&e_dvcs,"e_dvcs/D");

   TRandom *r = new TRandom3();
   Int_t iset = 0, ig = 0;
   TComplex cffs[4] = {0};
   Double_t k, QQ, xB, t;

   if (kKinematics == "Jlab_all")   dataReader.SetEntriesRange(0, -1); // Will load all 195 entries.
   if (kKinematics == "HallA_E00-110") dataReader.SetEntriesRange(0, 20); // Will load entries from 0 to 19.
   if (kKinematics == "HallA_E12-06-114") dataReader.SetEntriesRange(20, 65); // Will load entries from 20 to 64.
   if (kKinematics == "HallB")   dataReader.SetEntriesRange(65, 175); // Will load entries from 65 to 174.
   if (kKinematics == "HallA_E07–007")  dataReader.SetEntriesRange(175, -1); // Will load entries from 175 to 194.

   while (dataReader.Next()) {

      r->SetSeed(seed+iset);

      k = *k_jlab;
      QQ = *QQ_jlab;
      xB = *xB_jlab;
      t = *t_jlab;
      npoints = *npts_jlab;

      cout<<"No of points: "<<npoints<<endl;

      kin.k = k;
      kin.QQ = QQ;
      kin.xB = xB;
      kin.t = t;

      Double_t kine[4] = {QQ, xB, t, k};

      // Set CFFs values  ------------------------------------------------------
      if ( IsModel ) {
         if(kGPDModel == "model1"){ // out of the hat model 1
            gReH = -45.* pow( t , 4 ) + 10.* t + 1.45 / xB / xB - 7.;
            gReE = -1./ t / t + 40.* xB;
            gReHtilde =  40.* t + 5./ xB;
            gReEtilde = -15./ t + 5. / xB;
            gImH = 0;
            gImE = 0;
            gImHtilde = 0;
            gImEtilde = 0;
         }
         if(kGPDModel == "KM15"){   // CFFs from KM15 Modified model
            ModKM15_CFFs(kine, gReH, gImH, gReE, gReHtilde, gImHtilde, gReEtilde);
            gImE = 0;
            gImEtilde = 0;
         }
      }
      if ( !IsModel ) {
         tree_globalCFFs ->GetEntry(iset);
         gReH = f_global_ReH;
         gReE = f_global_ReE;
         gReHtilde = f_global_ReHtilde;
         gdvcs = f_global_dvcs;
      }

      // Fill graphs with true CFFs and dvcs values ----------------------------

      gReH_true ->SetPoint( iset, iset+1, gReH );
      gReH_true ->SetPointError( iset, 0.45, 0 );

      gReE_true ->SetPoint( iset, iset+1, gReE );
      gReE_true ->SetPointError( iset, 0.45, 0 );

      gReHt_true ->SetPoint( iset, iset+1, gReHtilde );
      gReHt_true ->SetPointError( iset, 0.45, 0 );

      // true dvcs
      if ( IsModel ) {
         cffs[0] = TComplex(gReH,gImH);   // H
         cffs[1] = TComplex(gReE,gImE);   // E
         cffs[2] = TComplex(gReHtilde,gImHtilde);  // Htilde
         cffs[3] = TComplex(gReEtilde,gImEtilde);  // Etilde

         if(kFormulation == "BKM02") gdvcs = bkm ->DVCS_UU_02(kine, 0, cffs, "t2"); // dvcs
         if(kFormulation == "BKM10") gdvcs = bkm ->DVCS_UU_10(kine, 0, cffs, "t2"); // dvcs
      }
      gdvcs_true ->SetPoint( iset, iset+1, gdvcs );
      gdvcs_true ->SetPointError( iset, 0.45, 0 );

      // Fix gen function parameters -------------------------------------------
      if ( IsModel ) {
         fgen->FixParameter(0, k); //k
         fgen->FixParameter(1, QQ); //QQ
         fgen->FixParameter(2, xB); //xB
         fgen->FixParameter(3, t); //t
         fgen->FixParameter(4, gReH);
         fgen->FixParameter(5, gReE);
         fgen->FixParameter(6, gReHtilde);
         fgen->FixParameter(7, gReEtilde);
         fgen->FixParameter(8, gImH);
         fgen->FixParameter(9, gImE);
         fgen->FixParameter(10, gImHtilde);
         fgen->FixParameter(11, gImEtilde);
      }
      if ( !IsModel ) {
         fgen->FixParameter(0, k); //k
         fgen->FixParameter(1, QQ); //QQ
         fgen->FixParameter(2, xB); //xB
         fgen->FixParameter(3, t); //t
         fgen->FixParameter(4, gReH);
         fgen->FixParameter(5, gReE);
         fgen->FixParameter(6, gReHtilde);
         fgen->FixParameter(7, gdvcs);
      }
      fBH->FixParameter(0, k); //k
      fBH->FixParameter(1, QQ); //QQ
      fBH->FixParameter(2, xB); //xB
      fBH->FixParameter(3, t); //t

      // Fix fit function parameters -------------------------------------------
      ffit->FixParameter(0, k); //k
      ffit->FixParameter(1, QQ); //QQ
      ffit->FixParameter(2, xB); //xB
      ffit->FixParameter(3, t); //t

      // Fill 2D graphs for ReH, ReE, ReHtilde, dvcs
      g2D_ReH->SetPoint(iset,t, xB, gReH );
      g2D_ReE->SetPoint(iset,t, xB, gReE );
      g2D_ReHtilde->SetPoint(iset,t, xB, gReHtilde );
      g2D_dvcs->SetPoint(iset,t, xB, gdvcs );

      // Cross section graphs
      gGenDVCS = new TGraphErrors(npoints);
      gGenDVCS ->SetName("total xs (F)");
      gGenDVCS ->SetMarkerStyle(20);
      gGenDVCS ->SetMarkerSize(1.8);
      gGenIntfDVCS = new TGraphErrors(npoints);
      gGenIntfDVCS ->SetName("Intf+dvcs");
      gGenIntfDVCS ->SetMarkerStyle(22);
      gGenIntfDVCS ->SetMarkerSize(1.8);

      //Generate cross section from kinematics
      for (Int_t i=0; i<npoints; i++) {

         phi[i] = phi_jlab[i];

         F[i] = r->Gaus(fgen->Eval(phi[i]),var*fgen->Eval(phi[i]));// find cross section with 5% variation in output
         //Double_t f1 = r->Uniform(F[i],0.4);// Generate Uniform variation
         errF[i] = var*F[i];
         //Double_t f2 = 0.1*F + r->Gaus(F*.01,0.001); // A simulated physical error (with absolute and relative contributions)

         // Fill pseudo data graph ---------------
         // Total xs
         gGenDVCS ->SetPoint( i, phi[i], F[i] );
         gGenDVCS ->SetPointError( i, 0, errF[i] );
         // F-BH = I+dvcs xs
         gGenIntfDVCS ->SetPoint( i, phi[i], F[i]-(fBH->Eval(phi[i])) );
         gGenIntfDVCS ->SetPointError( i, 0, errF[i] );

      }// end phi loop

      if ( F[0] < 0 ) { cout<<" \"Not Valid Kinematic Set\" --> Skipped"<<endl; continue; }

      // Fit raw generated pseudo-data
      //gGenDVCS ->Fit(ffit,"MQRB");  // This fit option gives larger errors
      //gGenDVCS ->Fit("ffit","FQR");
      DoROOTFit(gGenDVCS, ffit);

      // with pre-fit
      // gGenDVCS ->Fit(ffit, "QWRN"); // "initial pre-fit"
      // gGenDVCS ->Fit(ffit, "QMR+"); // "final fit"

      // Get fit parameters from total xs
      ReH = ffit ->GetParameter(4); //ReH
      ReE = ffit ->GetParameter(5); //ReE
      ReHtilde = ffit ->GetParameter(6); //ReHtilde
      e_ReH = ffit ->GetParError(4); //ReH fit error
      e_ReE = ffit ->GetParError(5); //ReE fit error
      e_ReHtilde = ffit ->GetParError(6); //ReHtilde fit error
      dvcs = ffit ->GetParameter(7); //ReHtilde
      e_dvcs = ffit ->GetParError(7); //ReH fit error

      gReH_fit ->SetPoint( iset, iset+1, ReH );
      gReH_fit ->SetPointError( iset, 0, e_ReH );

      gReE_fit ->SetPoint( iset, iset+1, ReE );
      gReE_fit ->SetPointError( iset, 0, e_ReE );

      gReHt_fit ->SetPoint( iset, iset+1, ReHtilde );
      gReHt_fit ->SetPointError( iset, 0, e_ReHtilde );

      gdvcs_fit ->SetPoint( iset, iset+1, dvcs );
      gdvcs_fit ->SetPointError( iset, 0, e_dvcs );

      // Get percent change
      Double_t pct_change_ReH = 100. * (ReH - gReH) / gReH ;
      Double_t pct_change_ReE = 100. * (ReE - gReE) / gReE ;
      Double_t pct_change_ReHtilde = 100. * (ReHtilde - gReHtilde) / gReHtilde;

      cout<<"set: "<<iset+1<<" QQ = "<<QQ<<", x_B = "<<xB<<", t = "<<t<<endl;
      cout<<"********* CFF extraction *********"<<endl;
      cout<<"   ReH  	%       |	ReE  	   %        | ReHtilde     %        |"<<endl;
      cout<<ReH<<" ( "<<gReH<<" )  "<<TMath::Abs(pct_change_ReH)<<" | "	<<ReE<<" ( "<<gReE<<" ) "<<TMath::Abs(pct_change_ReE)<<" | "<<ReHtilde<<" ( "<<gReHtilde<<" ) "<<TMath::Abs(pct_change_ReHtilde)<<" | "<<endl;
      cout<<endl;
      cout<<endl;

      // Print results to a .csv file
      for (Int_t i=0; i<npoints; i++) {
         myfile<<ig+1<<","<<i<<","<<k<<","<<QQ<<","<<xB<<","<<t<<","<<phi[i]<<","<<F[i]<<","<<errF[i]<<","<<var<<","<<gReH<<","
             <<gReE<<","<<gReHtilde<<","<<gdvcs<<endl;
      }

      // If selections are to be done the selected sets can be save in these graphs
      // Save the xs graphs for the "selected" kinematics settings
      gSelPseudoF[ig] = (TGraphErrors*)gGenDVCS->Clone(Form("gSelPseudoF_%d ", ig));
      gSelPseudoF[ig] ->SetTitle(Form("set %d: k = %.2f, Q^{2} = %.2f, xB = %.2f, t = %.2f; #phi [deg];d^{4}#sigma [nb/GeV^{4}]", iset+1, k, QQ, xB, t));
      gSelPseudoF[ig] ->SetName("total xs (F)");
      gSelPseudoIntfDVCS[ig] = (TGraphErrors*)gGenIntfDVCS->Clone(Form("gSelPseudoIntfDVCS_%d ", ig));
      gSelPseudoIntfDVCS[ig] ->SetTitle(Form("set %d: k = %.2f, Q^{2} = %.2f, xB = %.2f, t = %.2f; #phi [deg];d^{4}#sigma [nb/GeV^{4}]", iset+1, k, QQ, xB, t));
      gSelPseudoIntfDVCS[ig] ->SetName("Intf+dvcs");

      mgr_pseudo[ig] = new TMultiGraph();
      mgr_pseudo[ig] ->SetTitle(Form("set %d: k = %.2f, Q^{2} = %.2f, xB = %.2f, t = %.2f; #phi [deg];d^{4}#sigma [nb/GeV^{4}]", iset+1, k, QQ, xB, t));
      mgr_pseudo[ig] ->Add(gSelPseudoF[ig]);
      mgr_pseudo[ig] ->Add(gSelPseudoIntfDVCS[ig]);

      ig++;
      t3 ->Fill();

      iset++;

   }// end kinematics loop

   gROOT->SetBatch(kTRUE);
   for(Int_t i = 0; i<iset; i=i+5){

      c_F = new TCanvas(Form("c_F_%d",i),"test", 552, 274, 2198, 1710);
      c_F ->Divide(3,2);
      c_F ->cd(1);
      mgr_pseudo[i] ->Draw("ap");
      c_F ->cd(2);
      mgr_pseudo[i+1] ->Draw("ap");
      c_F ->cd(3);
      mgr_pseudo[i+2] ->Draw("ap");
      c_F ->cd(4);
      mgr_pseudo[i+3] ->Draw("ap");
      c_F ->cd(5);
      mgr_pseudo[i+4] ->Draw("ap");

      if(IsModel){
         if (i == 0) c_F->Print(Form("./Graphs/xs_pseudo_%s_%s_%s_t2_5pct.pdf(", kGPDModel.Data(), kFormulation.Data(), kKinematics.Data()),"pdf");
         else c_F->Print(Form("./Graphs/xs_pseudo_%s_%s_%s_t2_5pct.pdf", kGPDModel.Data(), kFormulation.Data(), kKinematics.Data()),"pdf");
      }
      if(!IsModel){
         if (i == 0) c_F->Print(Form("./Graphs/xs_pseudo_ANN-global_%s_%s_t2_5pct.pdf(", kFormulation.Data(), kKinematics.Data()),"pdf");
         else c_F->Print(Form("./Graphs/xs_pseudo_ANN-global_%s_%s_t2_5pct.pdf", kFormulation.Data(), kKinematics.Data()),"pdf");
      }
   }

   if(IsModel) c_F->Print(Form("./Graphs/xs_pseudo_%s_%s_%s_t2_5pct.pdf]", kGPDModel.Data(), kFormulation.Data(), kKinematics.Data()),"pdf");
   if(!IsModel) c_F->Print(Form("./Graphs/xs_pseudo_ANN-global_%s_%s_t2_5pct.pdf]", kFormulation.Data(), kKinematics.Data()),"pdf");

   gROOT->SetBatch(kFALSE);

   TCanvas *c_cffs = new TCanvas("c_cffs", "c_cffs",97,274,2473,690);
   c_cffs->Divide(4,1);
   c_cffs->cd(1);
   g2D_ReH->Draw("TRI1 P0");
   c_cffs->cd(2);
   g2D_ReE->Draw("TRI1 P0");
   c_cffs->cd(3);
   g2D_ReHtilde->Draw("TRI1 P0");
   c_cffs->cd(4);
   g2D_dvcs->Draw("TRI1 P0");

   if(IsModel) c_cffs ->Print(Form("./Graphs/trueCFFs_%s_%s_t2_5pct.pdf)", kGPDModel.Data(), kKinematics.Data()),"pdf");
   if(!IsModel) c_cffs ->Print(Form("./Graphs/trueCFFs_ANN-global_%s_%s_t2_5pct.pdf)", kFormulation.Data(), kKinematics.Data()),"pdf");

   fout.cd();
   t3->Write();

   // Draw extracted CFFs and true values
   gReH_true->SetName("true");
   gReE_true->SetName("true");
   gReHt_true->SetName("true");
   gdvcs_true->SetName("true");
   gReH_true ->SetLineColor(2);
   gReH_true ->SetMarkerSize(0);
   gReE_true ->SetLineColor(2);
   gReE_true ->SetMarkerSize(0);
   gReHt_true ->SetLineColor(2);
   gReHt_true ->SetMarkerSize(0);
   gdvcs_true ->SetLineColor(2);
   gdvcs_true ->SetMarkerSize(0);

   TMultiGraph *mgr2 = new TMultiGraph();
   mgr2 ->SetTitle("ReH - 5\% F - t2; set;ReH");
   mgr2 ->Add(gReH_true);
   mgr2 ->Add(gReH_fit);

   TMultiGraph *mgr3 = new TMultiGraph();
   mgr3 ->SetTitle("ReE- 5\% F - t2; set;ReE");
   mgr3 ->Add(gReE_true);
   mgr3 ->Add(gReE_fit);

   TMultiGraph *mgr4 = new TMultiGraph();
   mgr4 ->SetTitle("ReHtilde - 5\% F - t2; set;ReHtilde");
   mgr4 ->Add(gReHt_true);
   mgr4 ->Add(gReHt_fit);

   TMultiGraph *mgr5 = new TMultiGraph();
   mgr5 ->SetTitle("dvcs - 5\% F - t2; set;dvcs");
   mgr5 ->Add(gdvcs_true);
   mgr5 ->Add(gdvcs_fit);

   TCanvas * c2 = new TCanvas("c2","ReH", 4446,332,1685,899);
   mgr2 ->Draw("ap");
   gPad->BuildLegend(0.19, 0.8, 0.3, 0.9);
   if(IsModel) c2 ->Print(Form("./Graphs/ReH_%s_%s_%s_t2_5pct.pdf)", kGPDModel.Data(), kFormulation.Data(), kKinematics.Data()),"pdf");
   if(!IsModel) c2 ->Print(Form("./Graphs/ReH_ANN-global_%s_%s_t2_5pct.pdf)", kFormulation.Data(), kKinematics.Data()),"pdf");

   TCanvas * c3 = new TCanvas("c3","ReE", 4446,332,1685,899);
   mgr3 ->Draw("ap");
   gPad->BuildLegend(0.19, 0.8, 0.3, 0.9);
   if(IsModel) c3 ->Print(Form("./Graphs/ReE_%s_%s_%s_t2_5pct.pdf)", kGPDModel.Data(), kFormulation.Data(), kKinematics.Data()),"pdf");
   if(!IsModel) c3 ->Print(Form("./Graphs/ReE_ANN-global_%s_%s_t2_5pct.pdf)", kFormulation.Data(), kKinematics.Data()),"pdf");

   TCanvas * c4 = new TCanvas("c4","ReHtilde", 4446,332,1685,899);
   mgr4 ->Draw("ap");
   gPad->BuildLegend(0.19, 0.8, 0.3, 0.9);
   if(IsModel) c4 ->Print(Form("./Graphs/ReHtilde_%s_%s_%s_t2_5pct.pdf)", kGPDModel.Data(), kFormulation.Data(), kKinematics.Data()),"pdf");
   if(!IsModel) c4 ->Print(Form("./Graphs/ReHtilde_ANN-global_%s_%s_t2_5pct.pdf)", kFormulation.Data(), kKinematics.Data()),"pdf");

   TCanvas * c5 = new TCanvas("c5","dvcs", 4446,332,1685,899);
   mgr5 ->Draw("ap");
   gPad->BuildLegend(0.19, 0.8, 0.3, 0.9);
   if(IsModel) c5 ->Print(Form("./Graphs/dvcs_%s_%s_%s_t2_5pct.pdf)", kGPDModel.Data(), kFormulation.Data(), kKinematics.Data()),"pdf");
   if(!IsModel) c5 ->Print(Form("./Graphs/dvcs_ANN-global_%s_%s_t2_5pct.pdf)", kFormulation.Data(), kKinematics.Data()),"pdf");

}
//____________________________________________________________________________________________________
Double_t BH(Double_t *angle, Double_t *par)	// Total Cross Section - user defined function
{
    Double_t _phi = angle[0];
    Double_t _k = par[0];
    Double_t _QQ = par[1];
    Double_t _xB = par[2];
    Double_t _t = par[3];

    Double_t kine[4] = {_QQ, _xB, _t, _k};

    Double_t _F1 = ff->ffF1_K(_t);
    Double_t _F2 = ff->ffF2_K(_t);

    // Set QQ, xB, t and k
    Double_t xsbhuu = bkm ->BH_UU(kine, _phi, _F1, _F2); // BH cross section

    return xsbhuu;
}
//____________________________________________________________________________________________________
Double_t TotalUUXS_BKM02(Double_t *angle, Double_t *par)	// Total Cross Section - user defined function
{
    Double_t _phi = angle[0];
    Double_t _k = par[0];
    Double_t _QQ = par[1];
    Double_t _xB = par[2];
    Double_t _t = par[3];

    Double_t kine[4] = {_QQ, _xB, _t, _k};

    /* F = { H, E , Htilde, Etilde} Twist-2 Compton Form Factors*/
    TComplex _F[4] = {0};
    _F[0] = TComplex(par[4],par[8]);   // H
    _F[1] = TComplex(par[5],par[9]);   // E
    _F[2] = TComplex(par[6],par[10]);  // Htilde
    _F[3] = TComplex(par[7],par[11]);  // Etilde

    Double_t _F1 = ff->ffF1_K(_t);
    Double_t _F2 = ff->ffF2_K(_t);

    // Set QQ, xB, t and k
    Double_t xsbhuu = bkm ->BH_UU(kine, _phi, _F1, _F2); // BH cross section
    Double_t xsiuu = bkm ->I_UU_02(kine, _phi, _F1, _F2, _F, "t2"); // Interference
    Double_t xsdvcs = bkm ->DVCS_UU_02(kine, _phi, _F, "t2"); // dvcs

    Double_t tot_sigma_uu = xsbhuu + xsiuu + xsdvcs; // Constant added to account for DVCS contribution

    return tot_sigma_uu;
}
//____________________________________________________________________________________________________
Double_t TotalUUXS_BKM10(Double_t *angle, Double_t *par)	// Total Cross Section - user defined function
{
    Double_t _phi = angle[0];
    Double_t _k = par[0];
    Double_t _QQ = par[1];
    Double_t _xB = par[2];
    Double_t _t = par[3];

    Double_t kine[4] = {_QQ, _xB, _t, _k};

    /* F = { H, E , Htilde, Etilde} Twist-2 Compton Form Factors*/
    TComplex _F[4] = {0};
    _F[0] = TComplex(par[4],par[8]);   // H
    _F[1] = TComplex(par[5],par[9]);   // E
    _F[2] = TComplex(par[6],par[10]);  // Htilde
    _F[3] = TComplex(par[7],par[11]);  // Etilde

    Double_t _F1 = ff->ffF1_K(_t);
    Double_t _F2 = ff->ffF2_K(_t);

    // Set QQ, xB, t and k
    Double_t xsbhuu = bkm ->BH_UU(kine, _phi, _F1, _F2); // BH cross section
    Double_t xsiuu = bkm ->I_UU_10(kine, _phi, _F1, _F2, _F, "t2"); // Interference
    Double_t xsdvcs = bkm ->DVCS_UU_10(kine, _phi, _F, "t2"); // dvcs

    Double_t tot_sigma_uu = xsbhuu + xsiuu + xsdvcs; // Constant added to account for DVCS contribution



    return tot_sigma_uu;
}
//____________________________________________________________________________________________________
Double_t TotalUUXS_BKM02_Fit(Double_t *angle, Double_t *par)	// Total Cross Section - user defined function
{
    Double_t _phi = angle[0];
    Double_t _k = par[0];
    Double_t _QQ = par[1];
    Double_t _xB = par[2];
    Double_t _t = par[3];

    Double_t kine[4] = {_QQ, _xB, _t, _k};

    /* F = { H, E , Htilde, Etilde} Twist-2 Compton Form Factors*/
    TComplex _F[4] = {0};
    _F[0] = TComplex(par[4],0.);
    _F[1] = TComplex(par[5],0.);
    _F[2] = TComplex(par[6],0.);
    //_F[3] = TComplex(par[7],par[11]);

    Double_t xsdvcs = par[7];

    Double_t _F1 = ff->ffF1_K(_t);
    Double_t _F2 = ff->ffF2_K(_t);

    // Set QQ, xB, t and k
    Double_t xsbhuu	 = bkm ->BH_UU(kine, _phi, _F1, _F2); // BH cross section
    Double_t xsiuu = bkm ->I_UU_02(kine, _phi, _F1, _F2, _F, "t2"); // Interference

    Double_t tot_sigma_uu = xsbhuu + xsiuu + xsdvcs;

    return tot_sigma_uu;
}
//____________________________________________________________________________________________________
Double_t TotalUUXS_BKM10_Fit(Double_t *angle, Double_t *par)	// Total Cross Section - user defined function
{
    Double_t _phi = angle[0];
    Double_t _k = par[0];
    Double_t _QQ = par[1];
    Double_t _xB = par[2];
    Double_t _t = par[3];

    Double_t kine[4] = {_QQ, _xB, _t, _k};

    /* F = { H, E , Htilde, Etilde} Twist-2 Compton Form Factors*/
    TComplex _F[4] = {0};
    _F[0] = TComplex(par[4],0.);
    _F[1] = TComplex(par[5],0.);
    _F[2] = TComplex(par[6],0.);
    //_F[3] = TComplex(par[7],par[11]);

    Double_t xsdvcs = par[7];

    Double_t _F1 = ff->ffF1_K(_t);
    Double_t _F2 = ff->ffF2_K(_t);

    // Set QQ, xB, t and k
    Double_t xsbhuu	 = bkm ->BH_UU(kine, _phi, _F1, _F2); // BH cross section
    Double_t xsiuu = bkm ->I_UU_10(kine, _phi, _F1, _F2, _F, "t2"); // Interference

    Double_t tot_sigma_uu = xsbhuu + xsiuu + xsdvcs;

    return tot_sigma_uu;
}
//___________________________________________________________________________________________________
void DoROOTFit(TGraphErrors* gGenData, TF1* ffit ) { //Fit with ROOT::Fit method

   ROOT::Fit::DataOptions opt;
   ROOT::Fit::DataRange range;
   // set the data range
   range.SetRange(0,360);

   ROOT::Fit::BinData data(opt,range);
   ROOT::Fit::FillData(data, gGenData);

   ROOT::Math::WrappedMultiTF1 fitfunc(*ffit,1);
   ROOT::Fit::Fitter fitter;
   fitter.SetFunction(fitfunc, false);

   // fix parameters
   fitter.Config().ParSettings(0).Fix();
   fitter.Config().ParSettings(1).Fix();
   fitter.Config().ParSettings(2).Fix();
   fitter.Config().ParSettings(3).Fix();
   fitter.Config().ParSettings(7).SetLowerLimit(0);
   //fitter.Config().ParSettings(7).SetLimits(0,0.5);
   //fitter.Config().ParSettings(7).SetStepSize(1E-10);

   ROOT::Math::MinimizerOptions(mopt);
   mopt.SetMinimizerType("Minuit2");
   mopt.SetStrategy(0);
   mopt.SetTolerance(1);
   mopt.SetMaxFunctionCalls(1E5);
   mopt.SetPrecision(1E-10);
   mopt.SetPrintLevel(0);

   // print the default minimizer option values
   //mopt.Print();
   fitter.Config().SetMinimizerOptions(mopt);
   fitter.Fit(data); // chi2 fit'
   //fitter.LikelihoodFit(data);
   ROOT::Fit::FitResult result = fitter.Result();
   result.Print(std::cout);
   ffit ->SetFitResult(result);
   gGenData ->GetListOfFunctions()->Add(ffit);
}
