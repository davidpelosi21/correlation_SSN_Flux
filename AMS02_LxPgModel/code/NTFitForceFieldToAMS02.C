#include "TROOT.h"
#include "TMath.h"
#include "TGraphErrors.h"
#include "TGraph.h"
#include "TH2D.h"
#include "TAxis.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TApplication.h"
#include "TMinuit.h"
#include <iostream>

#include "NTCosmicRayData.cpp"
#include "NTAMS02Data.cpp"

using namespace std;


// --- graph style ---
void SetStyleHistoVSEkn(TH2D* hh);
void SetStyleHistoVSTime(TH2D* hh);
void SetTimeGraphStyle(TGraphErrors* gr);

// --- graph manipulation ---
void FlattenGraph(TGraph* gr, double OldFlatIndex, double NewFlatIndex); // for LIS
void ScaleGraph(TGraph* gr, double ScaleFactor);
void MakeDataToModel(TGraphErrors* grData, TGraph* grModel, TGraphErrors* grDataToModel);

// --- modulation ---
void ModulateGraph(TGraph* gr, double Phi, int Z, int A); // DEPRECATED!
void ModulateGraph(TGraph* grLIS, TGraph* grMMOD, double Phi, int Z, int A);
void SolarModulation(double *Ekin, double *cr_density, int NEKN, int z, int a, double phi, int kOpt);

// --- minimization ---
double GetChiSquare(double Phi, int iTime, int iExperiment);
void FCN(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);


// ---- constants ----
static const int nTimePAMELA_Proton  = 47;
static const int nTimeSOHO_Proton    = 20;
static const int nTimeBESS_Proton    = 2; //Polar I/II [no 00/TeV]
static const int nTimeAMSPRL2015_Proton   = 1;
static const int nTimeAMS02_Proton    = 81;
static const int nTimeAMS02_Helium    = 81;

const int NITERATIONS = 10000;

double EMIN= 0.07;
double EMAX= 100.;



double bestChi2_PAMELA[nTimePAMELA_Proton];
double bestChi2_SOHO[nTimeSOHO_Proton];
double bestChi2_BESS[nTimeBESS_Proton];
double bestChi2_AMSPRL2015[nTimeAMSPRL2015_Proton];
double bestChi2_AMS02[nTimeAMS02_Proton];
double bestChi2_AMS02Helium[nTimeAMS02_Helium];

int NDF_PAMELA[nTimePAMELA_Proton];
int NDF_SOHO[nTimeSOHO_Proton];
int NDF_BESS[nTimeBESS_Proton];
int NDF_AMSPRL2015[nTimeAMSPRL2015_Proton];
int NDF_AMS02[nTimeAMS02_Proton];
int NDF_AMS02Helium[nTimeAMS02_Helium];


// ---- external objects ----
TGraph* grLISFluxVSEkn_Proton;
TGraph* grLISFluxVSEkn_Helium;

NTCosmicRayData* crdata;
NTAMS02Data* amsdata;


// ---- results ----
TGraphErrors* grPHIvsTime_PAMELA;
TGraphErrors* grPHIvsTime_SOHO;
TGraphErrors* grPHIvsTime_BESS;
TGraphErrors* grPHIvsTime_AMSPRL2015;
TGraphErrors* grPHIvsTime_AMS02;
TGraphErrors* grPHIvsTime_AMS02Helium;

bool kSTORE= false; // store PHIvsTime graphs
bool kForceFieldSpectraSTORE = false;


// --- LIS options ---
// 1: based on TimeLag ApJL-2017 paper
// 2: based on BCUnc PRD-2017 paper
// 3: based on pHe ASR-2017 paper

int kLIS = 1; // [ApJL | PRD | ASR]

double NormProton= 1.;
double NormHelium= 1.;

// LIS1 w.r.t. paper [but it is already scaled]
// double NormProton= 0.97;
// double NormHelium= 1.03; 

// LIS2


int main(int argc, char **argv){ //TOP
  TApplication theApp("App",&argc,argv);

  gROOT->ProcessLine(".x $NTBASEDIR/THEORY/00crdata/style/MyRootStyle.C");  
  gROOT->ProcessLine(".L ../Matisse/lib/libNTMatisseLib.so");


  // ---- set cosmic-ray data ----
  crdata= new NTCosmicRayData();
  crdata->SetProtonData();

  // ---- set AMS02 data ----
  amsdata= new NTAMS02Data();
  amsdata->SetProtonData();
  amsdata->SetHeliumData();



  // ---- graph with results ----
  grPHIvsTime_PAMELA= new TGraphErrors();
  grPHIvsTime_PAMELA->SetName("grPhiVSTime_PAMELA");
  grPHIvsTime_PAMELA->SetMarkerStyle( crdata->grPAMELA_ProtonFluxVSTime[0]->GetMarkerStyle() );
  grPHIvsTime_PAMELA->SetMarkerColor( crdata->grPAMELA_ProtonFluxVSTime[0]->GetMarkerColor() );
  grPHIvsTime_PAMELA->SetMarkerSize( crdata->grPAMELA_ProtonFluxVSTime[0]->GetMarkerSize() );
  grPHIvsTime_PAMELA->SetLineWidth( crdata->grPAMELA_ProtonFluxVSTime[0]->GetLineWidth() );
  grPHIvsTime_PAMELA->SetLineColor( crdata->grPAMELA_ProtonFluxVSTime[0]->GetLineColor() );
  SetTimeGraphStyle(grPHIvsTime_PAMELA);

  grPHIvsTime_SOHO= new TGraphErrors();
  grPHIvsTime_SOHO->SetName("grPhiVSTime_SOHO");
  grPHIvsTime_SOHO->SetMarkerStyle( crdata->grSOHO_ProtonFluxVSTime[0]->GetMarkerStyle() );
  grPHIvsTime_SOHO->SetMarkerColor( crdata->grSOHO_ProtonFluxVSTime[0]->GetMarkerColor() );
  grPHIvsTime_SOHO->SetMarkerSize( crdata->grSOHO_ProtonFluxVSTime[0]->GetMarkerSize() );
  grPHIvsTime_SOHO->SetLineWidth( crdata->grSOHO_ProtonFluxVSTime[0]->GetLineWidth() );
  grPHIvsTime_SOHO->SetLineColor( crdata->grSOHO_ProtonFluxVSTime[0]->GetLineColor() );
  SetTimeGraphStyle(grPHIvsTime_SOHO);

  grPHIvsTime_BESS= new TGraphErrors();
  grPHIvsTime_BESS->SetName("grPhiVSTime_BESS");
  grPHIvsTime_BESS->SetMarkerStyle( crdata->grBESS_ProtonFluxVSTime[0]->GetMarkerStyle() );
  grPHIvsTime_BESS->SetMarkerColor( crdata->grBESS_ProtonFluxVSTime[0]->GetMarkerColor() );
  grPHIvsTime_BESS->SetMarkerSize( crdata->grBESS_ProtonFluxVSTime[0]->GetMarkerSize() );
  grPHIvsTime_BESS->SetLineWidth( crdata->grBESS_ProtonFluxVSTime[0]->GetLineWidth() );
  grPHIvsTime_BESS->SetLineColor( crdata->grBESS_ProtonFluxVSTime[0]->GetLineColor() );
  SetTimeGraphStyle(grPHIvsTime_BESS);

  grPHIvsTime_AMSPRL2015= new TGraphErrors();
  grPHIvsTime_AMSPRL2015->SetName("grPhiVSTime_AMSPRL2015");
  // grPHIvsTime_AMSPRL2015->SetMarkerStyle( crdata->grAMS02_ProtonFluxVSTime[0]->GetMarkerStyle() );
  // grPHIvsTime_AMSPRL2015->SetMarkerColor( crdata->grAMS02_ProtonFluxVSTime[0]->GetMarkerColor() );
  // grPHIvsTime_AMSPRL2015->SetMarkerSize( crdata->grAMS02_ProtonFluxVSTime[0]->GetMarkerSize() );
  // grPHIvsTime_AMSPRL2015->SetLineWidth( crdata->grAMS02_ProtonFluxVSTime[0]->GetLineWidth() );
  // grPHIvsTime_AMSPRL2015->SetLineColor( crdata->grAMS02_ProtonFluxVSTime[0]->GetLineColor() );

  // NT SET LARGE LINE
  grPHIvsTime_AMSPRL2015->SetLineWidth( 3 );
  grPHIvsTime_AMSPRL2015->SetLineColor( kRed+2 );
  grPHIvsTime_AMSPRL2015->SetMarkerSize( 0 );
  grPHIvsTime_AMSPRL2015->SetMarkerColor( kRed+2 );
  SetTimeGraphStyle(grPHIvsTime_AMSPRL2015);

  
  grPHIvsTime_AMS02= new TGraphErrors();
  grPHIvsTime_AMS02->SetName("grPhiVSTime_AMS02");
  grPHIvsTime_AMS02->SetMarkerStyle( amsdata->grAMS02_ProtonFluxVSTime[0]->GetMarkerStyle() );
  grPHIvsTime_AMS02->SetMarkerColor( amsdata->grAMS02_ProtonFluxVSTime[0]->GetMarkerColor() );
  grPHIvsTime_AMS02->SetMarkerSize( amsdata->grAMS02_ProtonFluxVSTime[0]->GetMarkerSize() );
  grPHIvsTime_AMS02->SetLineWidth( amsdata->grAMS02_ProtonFluxVSTime[0]->GetLineWidth() );
  grPHIvsTime_AMS02->SetLineColor( amsdata->grAMS02_ProtonFluxVSTime[0]->GetLineColor() );
  SetTimeGraphStyle(grPHIvsTime_AMS02);

  grPHIvsTime_AMS02Helium= new TGraphErrors();
  grPHIvsTime_AMS02Helium->SetName("grPhiVSTime_AMS02Helium");
  grPHIvsTime_AMS02Helium->SetMarkerStyle( amsdata->grAMS02_HeliumFluxVSTime[0]->GetMarkerStyle() );
  grPHIvsTime_AMS02Helium->SetMarkerColor( amsdata->grAMS02_HeliumFluxVSTime[0]->GetMarkerColor() );
  grPHIvsTime_AMS02Helium->SetMarkerSize( amsdata->grAMS02_HeliumFluxVSTime[0]->GetMarkerSize() );
  grPHIvsTime_AMS02Helium->SetLineWidth( amsdata->grAMS02_HeliumFluxVSTime[0]->GetLineWidth() );
  grPHIvsTime_AMS02Helium->SetLineColor( amsdata->grAMS02_HeliumFluxVSTime[0]->GetLineColor() );
  SetTimeGraphStyle(grPHIvsTime_AMS02Helium);

  

  // ---- pick up models of LIS fluxes ---- [NT due modelli da considerare!]

  // --- LIS from ApJL-2017 TimeLag ---
  if(kLIS==1){
    TFile *inProtonLIS= new TFile("../SimulatedData/grLIS_Proton_NT_ApJL_2017_TimeLag.root","READ");
    grLISFluxVSEkn_Proton= (TGraph*)(inProtonLIS->Get("grLISFluxVSEkn_Proton"));
    inProtonLIS->Close();

    TFile *inHeliumLIS= new TFile("../SimulatedData/grLIS_Helium_NT_ApJL_2017_TimeLag.root","READ");
    grLISFluxVSEkn_Helium= (TGraph*)(inHeliumLIS->Get("grLISFluxVSEkn_Helium"));
    inHeliumLIS->Close();

    NormProton= 1.;
    NormHelium= 1.;
  }

  // grLISFluxVSEkn_Proton=new TGraph("../SimulatedData/LIS_Proton.dat");
  // grLISFluxVSEkn_Helium=new TGraph("../SimulatedData/LIS_Helium.dat");

  
  // --- LIS from PRD-2017 BCUnc ---
  if(kLIS==2){
    TFile *inProtonLIS= new TFile("../SimulatedData/grLIS_Proton_NT_PRD_2017_BCUnc.root","READ");
    grLISFluxVSEkn_Proton= (TGraph*)(inProtonLIS->Get("grLISFluxVSEkn_Proton"));
    inProtonLIS->Close();

    TFile *inHeliumLIS= new TFile("../SimulatedData/grLIS_Helium_NT_PRD_2017_BCUnc.root","READ");
    grLISFluxVSEkn_Helium= (TGraph*)(inHeliumLIS->Get("grLISFluxVSEkn_Helium"));
    inHeliumLIS->Close();

    NormProton= 1.;
    NormHelium= 1.;
  }

  // --- LIS from ASR-2017 pHe ---
  if(kLIS==3){
    TFile *inProtonLIS= new TFile("../SimulatedData/grLIS_Proton_NT_ASR_2017_SolarPHe.root","READ");
    grLISFluxVSEkn_Proton= (TGraph*)(inProtonLIS->Get("grLISFluxVSEkn_Proton"));
    inProtonLIS->Close();

    TFile *inHeliumLIS= new TFile("../SimulatedData/grLIS_Helium_NT_ASR_2017_SolarPHe.root","READ");
    grLISFluxVSEkn_Helium= (TGraph*)(inHeliumLIS->Get("grLISFluxVSEkn_Helium"));
    inHeliumLIS->Close();

    NormProton= 1.;
    NormHelium= 1.;
  }
    
    // TFile *inPHeLIS= new TFile("../SimulatedData/grLIS_ProtonHelium_ASR2017.root","READ");
    // grLISFluxVSEkn_Proton= (TGraph*)(inPHeLIS->Get("grLISFluxH"));
    // FlattenGraph(grLISFluxVSEkn_Proton, 2.5, 0.0);
    // grLISFluxVSEkn_Helium= (TGraph*)(inPHeLIS->Get("grLISFluxHe"));
    // FlattenGraph(grLISFluxVSEkn_Helium, 2.5, 0.0);
    // inPHeLIS->Close();

  
  
  // --- LIS style graphs ---
  grLISFluxVSEkn_Proton->SetName("grLISFluxVSEkn_Proton");
  grLISFluxVSEkn_Proton->SetLineColor(kBlue+1);
  grLISFluxVSEkn_Proton->SetLineWidth(3);
  ScaleGraph(grLISFluxVSEkn_Proton, NormProton);

  grLISFluxVSEkn_Helium->SetName("grLISFluxVSEkn_Helium");
  grLISFluxVSEkn_Helium->SetLineColor(kBlue+1);
  grLISFluxVSEkn_Helium->SetLineWidth(3);
  ScaleGraph(grLISFluxVSEkn_Helium, NormHelium);

  
  // --- NTDEC2017 rivedere NORM! ---
  // NO: NORM IS APPLIED TO "MODULATED" FLUX INSIDE THE MODULATION ROUTINE [->NOW REMOVED!]
  // ScaleGraph(grLISFluxVSEkn_Proton, 0.97); // _NormProton!!

  
  // ---- best-fit modulated graphs ----
  TGraph* grFluxEknH_PAMELA[nTimePAMELA_Proton];
  TGraph* grFluxEknH_SOHO[nTimeSOHO_Proton];
  TGraph* grFluxEknH_BESS[nTimeBESS_Proton];
  TGraph* grFluxEknH_AMSPRL2015[nTimeAMSPRL2015_Proton];
  TGraph* grFluxEknH_AMS02[nTimeAMS02_Proton];
  TGraph* grFluxEknH_AMS02Helium[nTimeAMS02_Helium];

  for(int tt=0;tt<nTimePAMELA_Proton;tt++){
    grFluxEknH_PAMELA[tt]= (TGraph*)grLISFluxVSEkn_Proton->Clone(Form("grFluxEknH_PAMELA_T%d",tt));
    grFluxEknH_PAMELA[tt]->SetName(Form("grFluxEknH_PAMELA_T%d",tt));
    grFluxEknH_PAMELA[tt]->SetLineColor(kRed+1);
    grFluxEknH_PAMELA[tt]->SetLineWidth(2);
  }

  for(int tt=0;tt<nTimeSOHO_Proton;tt++){
    grFluxEknH_SOHO[tt]= (TGraph*)grLISFluxVSEkn_Proton->Clone(Form("grFluxEknH_SOHO_T%d",tt));
    grFluxEknH_SOHO[tt]->SetName(Form("grFluxEknH_SOHO_T%d",tt));
    grFluxEknH_SOHO[tt]->SetLineColor(kRed+1);
    grFluxEknH_SOHO[tt]->SetLineWidth(2);
  }

  for(int tt=0;tt<nTimeBESS_Proton;tt++){
    grFluxEknH_BESS[tt]= (TGraph*)grLISFluxVSEkn_Proton->Clone(Form("grFluxEknH_BESS_T%d",tt));
    grFluxEknH_BESS[tt]->SetName(Form("grFluxEknH_BESS_T%d",tt));
    grFluxEknH_BESS[tt]->SetLineColor(kRed+1);
    grFluxEknH_BESS[tt]->SetLineWidth(2);
  }

  for(int tt=0;tt<nTimeAMSPRL2015_Proton;tt++){
    grFluxEknH_AMSPRL2015[tt]= (TGraph*)grLISFluxVSEkn_Proton->Clone(Form("grFluxEknH_AMSPRL2015_T%d",tt));
    grFluxEknH_AMSPRL2015[tt]->SetName(Form("grFluxEknH_AMSPRL2015_T%d",tt));
    grFluxEknH_AMSPRL2015[tt]->SetLineColor(kRed+1);
    grFluxEknH_AMSPRL2015[tt]->SetLineWidth(2);
  }

  for(int tt=0;tt<nTimeAMS02_Proton;tt++){
    grFluxEknH_AMS02[tt]= (TGraph*)grLISFluxVSEkn_Proton->Clone(Form("grFluxEknH_AMS02_T%d",tt));
    grFluxEknH_AMS02[tt]->SetName(Form("grFluxEknH_AMS02_T%d",tt));
    grFluxEknH_AMS02[tt]->SetLineColor(kRed+1);
    grFluxEknH_AMS02[tt]->SetLineWidth(2);
  }

  for(int tt=0;tt<nTimeAMS02_Proton;tt++){
    grFluxEknH_AMS02Helium[tt]= (TGraph*)grLISFluxVSEkn_Helium->Clone(Form("grFluxEknH_AMS02Helium_T%d",tt));
    grFluxEknH_AMS02Helium[tt]->SetName(Form("grFluxEknHe_AMS02Helium_T%d",tt));
    grFluxEknH_AMS02Helium[tt]->SetLineColor(kRed+1);
    grFluxEknH_AMS02Helium[tt]->SetLineWidth(2);
  }


  for(int tt=0;tt<nTimePAMELA_Proton;tt++) bestChi2_PAMELA[tt]= 1.e+9;
  for(int tt=0;tt<nTimeSOHO_Proton;tt++) bestChi2_SOHO[tt]= 1.e+9;
  for(int tt=0;tt<nTimeBESS_Proton;tt++) bestChi2_BESS[tt]= 1.e+9;
  for(int tt=0;tt<nTimeAMSPRL2015_Proton;tt++) bestChi2_AMSPRL2015[tt]= 1.e+9;
  for(int tt=0;tt<nTimeAMS02_Proton;tt++) bestChi2_AMS02[tt]= 1.e+9;
  for(int tt=0;tt<nTimeAMS02_Helium;tt++) bestChi2_AMS02Helium[tt]= 1.e+9;


  // ----- minimization ----
  const int NPAR = 3; // PHI, iExp, iTime, [ Emin-Emax from ABOVE]
  TMinuit* gMinuit = new TMinuit(NPAR);  
  gMinuit->SetPrintLevel(-1); // -1 to 3
  gMinuit->SetFCN(FCN);

  Double_t arglist[10];
  Int_t    ierflg = 0;
  arglist[0] = 1;
  gMinuit->mnexcm("SET ERR", arglist ,1,ierflg);

  double minPhi= 0.10;
  double maxPhi= 1.90;
  double stepPhi=1.e-5;
  double startPhi=0.5;
  
  // [ iExp | 1:PAMELA | 2:SOHO | 3:BESS | 4:AMS02 vs T  | 5:AMS02-He vs T | 9:AMSPRL2015-Full ]
  
  // ---- fitting loop: PAMELA ---
  cout<<"**** Fit to PAMELA ****"<<endl;

  for(int tt=0;tt<nTimePAMELA_Proton;tt++){

    gMinuit->DefineParameter(0, "Phi", startPhi, stepPhi, minPhi, maxPhi);
    gMinuit->DefineParameter(1, "Experiment", 1, 0, 0, 0);
    gMinuit->DefineParameter(2, "Time Index", tt, 0, 0, 0);
    gMinuit->FixParameter(1); 
    gMinuit->FixParameter(2); 
    
    arglist[0] = NITERATIONS;
    arglist[1] = 1.;
    gMinuit->mnexcm("MIGRAD", arglist ,2,ierflg);
    
    double bestPhi= 0;
    double errPhi = 0.;
    gMinuit->GetParameter(0,bestPhi,errPhi);

    // update starting value for next fit
    startPhi= bestPhi; 
    
    // convert GV to MV
    bestPhi *= 1000.;
    errPhi  *= 1000.;
    errPhi  *= 2.;
    
    // ---- put results into graph ----
    grPHIvsTime_PAMELA->SetPoint(tt, crdata->xTimePAMELA_ProtonFlux[tt], bestPhi);
    grPHIvsTime_PAMELA->SetPointError(tt, crdata->eTimePAMELA_ProtonFlux, errPhi);
    cout <<"iTime: "<<tt<<"   PHI= " << bestPhi << " +- "<<errPhi<<"      GV"<<"    Chi2: "<<bestChi2_PAMELA[tt]<< " / "<<NDF_PAMELA[tt]<<endl;
    
    // ---- set best-fit flux vs ekn calculations ----
    ModulateGraph(grFluxEknH_PAMELA[tt], bestPhi/1000., 1, 1);
  }
  
  
  
  // ---- fitting loop: SOHO ---
  cout<<"**** Fit to SOHO ****"<<endl;
  for(int tt=0;tt<nTimeSOHO_Proton;tt++){
    
    gMinuit->DefineParameter(0, "Phi", startPhi, stepPhi, minPhi, maxPhi);
    gMinuit->DefineParameter(1, "Experiment", 2, 0, 0, 0);
    gMinuit->DefineParameter(2, "Time Index", tt, 0, 0, 0);
    gMinuit->FixParameter(1); 
    gMinuit->FixParameter(2); 
    
    arglist[0] = NITERATIONS;
    arglist[1] = 1.;
    gMinuit->mnexcm("MIGRAD", arglist ,2,ierflg);

    double bestPhi= 0.;
    double errPhi = 0.;
    gMinuit->GetParameter(0,bestPhi,errPhi);

    // update starting value for next fit
    startPhi= bestPhi; 
    
    // convert GV to MV
    bestPhi *= 1000.;
    errPhi  *= 1000.;
    errPhi  *= 2.;
    
    // ---- put results into graph ----
    grPHIvsTime_SOHO->SetPoint(tt, crdata->xTimeSOHO_ProtonFlux[tt], bestPhi);
    grPHIvsTime_SOHO->SetPointError(tt, crdata->eTimeSOHO_ProtonFlux, errPhi);
    cout <<"iTime: "<<tt<<"   PHI= " << bestPhi << " +- "<<errPhi<<"      GV"<<"    Chi2: "<<bestChi2_SOHO[tt]<< " / "<<NDF_SOHO[tt]<<endl;
    
    // ---- set best-fit flux vs ekn calculations ----
    ModulateGraph(grFluxEknH_SOHO[tt], bestPhi/1000., 1, 1);
  }



  // ---- fitting loop: BESS ---
  cout<<"**** Fit to BESS ****"<<endl;
  for(int tt=0;tt<nTimeBESS_Proton;tt++){
    
      gMinuit->DefineParameter(0, "Phi", startPhi, stepPhi, minPhi, maxPhi);
      gMinuit->DefineParameter(1, "Experiment", 3, 0, 0, 0);
      gMinuit->DefineParameter(2, "Time Index", tt, 0, 0, 0);
      gMinuit->FixParameter(1); 
      gMinuit->FixParameter(2); 

      arglist[0] = NITERATIONS;
      arglist[1] = 1.;
      gMinuit->mnexcm("MIGRAD", arglist ,2,ierflg);

      double bestPhi= 0.;
      double errPhi = 0.;
      gMinuit->GetParameter(0,bestPhi,errPhi);

      
      // update starting value for next fit
      startPhi= bestPhi; 

      // convert GV to MV
      bestPhi *= 1000.;
      errPhi  *= 1000.;
      errPhi  *= 2.;

      // ---- put results into graph ----
      grPHIvsTime_BESS->SetPoint(tt, crdata->xTimeBESS_ProtonFlux[tt], bestPhi);
      grPHIvsTime_BESS->SetPointError(tt, crdata->eTimeBESS_ProtonFlux[tt], errPhi);
      cout <<"iTime: "<<tt<<"   PHI= " << bestPhi << " +- "<<errPhi<<"      GV"<<"    Chi2: "<<bestChi2_BESS[tt]<< " / "<<NDF_BESS[tt]<<endl;
      
      // ---- set best-fit flux vs ekn calculations ----
      ModulateGraph(grFluxEknH_BESS[tt], bestPhi/1000., 1, 1);
  }


    
  // ---- fitting loop: AMS02 ---
  cout<<"**** Fit to AMS02 Proton ****"<<endl;
  int iPoint=0;
  for(int tt=0;tt<nTimeAMS02_Proton;tt++){
    int indTime = tt;
    if(tt==amsdata->SkipThis || tt==amsdata->SkipThis+1) continue;
    // cout<<"MAKE FIT TO AMS.... T"<<indTime<<endl;
    
      gMinuit->DefineParameter(0, "Phi", startPhi, stepPhi, minPhi, maxPhi);
      gMinuit->DefineParameter(1, "Experiment", 4, 0, 0, 0);
      gMinuit->DefineParameter(2, "Time Index", indTime, 0, 0, 0);
      gMinuit->FixParameter(1); 
      gMinuit->FixParameter(2); 

      arglist[0] = NITERATIONS;
      arglist[1] = 1.;
      gMinuit->mnexcm("MIGRAD", arglist ,2,ierflg);

      double bestPhi= 0.;
      double errPhi = 0.;
      gMinuit->GetParameter(0,bestPhi,errPhi);

      
      // update starting value for next fit
      startPhi= bestPhi; 

      // convert GV to MV
      bestPhi *= 1000.;
      errPhi  *= 1000.;
      errPhi  *= 2.;

      // ---- put results into graph ----
      grPHIvsTime_AMS02->SetPoint(iPoint, amsdata->xTimeAMS02_ProtonFlux[indTime], bestPhi);
      grPHIvsTime_AMS02->SetPointError(iPoint, amsdata->eTimeAMS02_ProtonFlux, errPhi);
      cout <<"iTime: "<<indTime<<"   PHI= " << bestPhi << " +- "<<errPhi<<"      GV"<<"    Chi2: "<<bestChi2_AMS02[indTime]<< " / "<<NDF_AMS02[indTime]<<endl;
      iPoint++;

      // ---- set best-fit flux vs ekn calculations ----
      ModulateGraph(grFluxEknH_AMS02[indTime], bestPhi/1000., 1, 1);
  }


  

  // ---- fitting loop: AMS02 Helium --- PER ORA NO
  cout<<"**** Fit to AMS02 Helium ****"<<endl;
  iPoint=0;
  for(int tt=0;tt<nTimeAMS02_Helium;tt++){
    int indTime = tt;
    if(tt==amsdata->SkipThis || tt==amsdata->SkipThis+1) continue;
    // cout<<"MAKE FIT TO AMS.... T"<<indTime<<endl;
    
      gMinuit->DefineParameter(0, "Phi", startPhi, stepPhi, minPhi, maxPhi);
      gMinuit->DefineParameter(1, "Experiment", 5, 0, 0, 0);
      gMinuit->DefineParameter(2, "Time Index", indTime, 0, 0, 0);
      gMinuit->FixParameter(1); 
      gMinuit->FixParameter(2); 

      arglist[0] = NITERATIONS;
      arglist[1] = 1.;
      gMinuit->mnexcm("MIGRAD", arglist ,2,ierflg);

      double bestPhi= 0.;
      double errPhi = 0.;
      gMinuit->GetParameter(0,bestPhi,errPhi);

      
      // update starting value for next fit
      startPhi= bestPhi; 

      // convert GV to MV
      bestPhi *= 1000.;
      errPhi  *= 1000.;
      errPhi  *= 2.;

      // ---- put results into graph ----
      grPHIvsTime_AMS02Helium->SetPoint(iPoint, amsdata->xTimeAMS02_HeliumFlux[indTime], bestPhi);
      grPHIvsTime_AMS02Helium->SetPointError(iPoint, amsdata->eTimeAMS02_HeliumFlux, errPhi);
      cout <<"iTime: "<<indTime<<"   PHI= " << bestPhi << " +- "<<errPhi<<"      GV"<<"    Chi2: "<<bestChi2_AMS02Helium[indTime]<< " / "<<NDF_AMS02Helium[indTime]<<endl;
      iPoint++;

      // ---- set best-fit flux vs ekn calculations ----
      ModulateGraph(grFluxEknH_AMS02Helium[indTime], bestPhi/1000., 2, 4); // Z-A for Helium4 !!
  }

  

  // ---- fitting loop: AMS02 PRL-2015: one point ---
  cout<<"**** Fit to AMSPRL2015 ****"<<endl;
  for(int tt=0;tt<nTimeAMSPRL2015_Proton;tt++){
    
      gMinuit->DefineParameter(0, "Phi", startPhi, stepPhi, minPhi, maxPhi);
      gMinuit->DefineParameter(1, "Experiment", 9, 0, 0, 0);
      gMinuit->DefineParameter(2, "Time Index", tt, 0, 0, 0);
      gMinuit->FixParameter(1); 
      gMinuit->FixParameter(2); 

      arglist[0] = NITERATIONS;
      arglist[1] = 1.;
      gMinuit->mnexcm("MIGRAD", arglist ,2,ierflg);

      double bestPhi= 0.;
      double errPhi = 0.;
      gMinuit->GetParameter(0,bestPhi,errPhi);

      
      // update starting value for next fit
      startPhi= bestPhi; 

      // convert GV to MV
      bestPhi *= 1000.;
      errPhi  *= 1000.;
      errPhi  *= 2.;

      // ---- put results into graph ----
      grPHIvsTime_AMSPRL2015->SetPoint(tt, crdata->xTimeAMS02_ProtonFlux[tt], bestPhi);
      grPHIvsTime_AMSPRL2015->SetPointError(tt, crdata->eTimeAMS02_ProtonFlux, errPhi);
      cout <<"iTime: "<<tt<<"   PHI= " << bestPhi << " +- "<<errPhi<<"      GV"<<"    Chi2: "<<bestChi2_AMSPRL2015[tt]<< " / "<<NDF_AMSPRL2015[tt]<<endl;

      // ---- set best-fit flux vs ekn calculations ----
      // bestPhi=575.;
      ModulateGraph(grFluxEknH_AMSPRL2015[tt], bestPhi/1000., 1, 1);
  }



  // ---- draw results ---
  TDatime Date1( 2000, 01, 01, 0, 0, 0);
  TDatime Date2( 2018, 01, 01, 0, 0, 0);
  double UTMIN= (double)Date1.Convert();
  double UTMAX= (double)Date2.Convert();
  double FMIN = -10.;
  double FMAX = 2200.;

  TH2D* hFramePHIvsTime= new TH2D("hFramePHIvsTime","PHI vs Time",200, UTMIN, UTMAX, 500, 0.20, 1200.0);
  SetStyleHistoVSTime(hFramePHIvsTime);
  hFramePHIvsTime->GetYaxis()->SetNdivisions(506);
  hFramePHIvsTime->GetYaxis()->SetTitleSize(0.07);
  hFramePHIvsTime->GetYaxis()->SetTitle("#phi (MV)");
  hFramePHIvsTime->GetYaxis()->SetTitleOffset(0.50);

  // ---- plot fluxes ----
  TCanvas* ccBestFitPhiVSTime= new TCanvas(Form("ccBestFitPhiVSTime_N%d",0), Form("BEST PHI FROM FIT TO CR DATA | N%d",0), 1200,580 );

  ccBestFitPhiVSTime->cd();  
  gPad->SetBottomMargin(0.15);
  gPad->SetTopMargin(0.05);
  gPad->SetLeftMargin(0.150);
  gPad->SetRightMargin(0.06);
  
  gPad->SetTicky(1);
  gPad->SetTickx(1); 
  hFramePHIvsTime->Draw();

  SetTimeGraphStyle(grPHIvsTime_PAMELA);
  SetTimeGraphStyle(grPHIvsTime_SOHO);
  SetTimeGraphStyle(grPHIvsTime_BESS);
  SetTimeGraphStyle(grPHIvsTime_AMSPRL2015);
  
  grPHIvsTime_SOHO->Draw("pZ");
  grPHIvsTime_PAMELA->Draw("pZ");
  grPHIvsTime_BESS->Draw("pZ");
  grPHIvsTime_AMSPRL2015->Draw("pZ");
  grPHIvsTime_AMS02->Draw("pZl");
  grPHIvsTime_AMS02Helium->Draw("pZl");

  gPad->SetGridy();
  ccBestFitPhiVSTime->Update();


  // ---- plot LIS and energy spectra ----

  
  // ---- frames ----
  TH2D* hFrameProtonFluxVSEkn= new TH2D("hFrameProtonFluxVSEkn","hFrameProtonFluxVSEkn",200, 0.08, 50.0, 200, 0.3, 3.e+4);
  SetStyleHistoVSEkn(hFrameProtonFluxVSEkn);
  hFrameProtonFluxVSEkn->GetXaxis()->SetTitle("kinetic energy (GeV)");
  hFrameProtonFluxVSEkn->GetYaxis()->SetNdivisions(506);
  hFrameProtonFluxVSEkn->GetYaxis()->SetTitle("J ( GeV^{ -1} m^{ -2} s^{ -1} sr^{ -1} )");
  
  TH2D* hFrameHeliumFluxVSEkn= new TH2D("hFrameHeliumFluxVSEkn","hFrameHeliumFluxVSEkn",200, 0.08, 50.0, 200, 0.1, 1.e+4);
  SetStyleHistoVSEkn(hFrameHeliumFluxVSEkn);
  hFrameHeliumFluxVSEkn->GetXaxis()->SetTitle("kinetic energy (GeV/n)");
  hFrameHeliumFluxVSEkn->GetYaxis()->SetNdivisions(506);
  hFrameHeliumFluxVSEkn->GetYaxis()->SetTitle("J ( (GeV/n)^{ -1} m^{ -2} s^{ -1} sr^{ -1} )");

  // ---- plot fluxes ----
  TCanvas* ccProtonFluxEnergySpectra= new TCanvas(Form("ccProtonFluxEnergySpectr_N%d",2), Form("Proton Flux Energy Spectrum | N%d",0), 1200,580 );
  ccProtonFluxEnergySpectra->Divide(2,1,0.01,0.01);

  // proton
  ccProtonFluxEnergySpectra->cd(1);  
  gPad->SetBottomMargin(0.15);
  gPad->SetTopMargin(0.05);
  gPad->SetLeftMargin(0.150);
  gPad->SetRightMargin(0.06);
  gPad->SetTicky(1);
  gPad->SetTickx(1); 
  hFrameProtonFluxVSEkn->Draw();
  gPad->SetLogx();
  gPad->SetLogy();    

  crdata->grPAMELA_ProtonFluxVSEkn[2]->Draw("pZ");
  crdata->grPAMELA_ProtonFluxVSEkn[25]->Draw("pZ"); 
  crdata->grSOHO_ProtonFluxVSEkn[18]->Draw("pZ");
  amsdata->grAMS02_ProtonFluxVSEkn[4]->Draw("pZ");
  amsdata->grAMS02_ProtonFluxVSEkn[40]->Draw("pZ");

  // grFluxEknH_PAMELA[2]->Draw("c");
  // grFluxEknH_PAMELA[25]->Draw("c");
  // grFluxEknH_SOHO[18]->Draw("c");
  grFluxEknH_AMS02[4]->Draw("c");
  grFluxEknH_AMS02[40]->Draw("c");

  amsdata->grVOYAGER1_ProtonFluxVSEkn->Draw("pZ");
  grLISFluxVSEkn_Proton->Draw("c");

  // helium
  ccProtonFluxEnergySpectra->cd(2);  
  gPad->SetBottomMargin(0.15);
  gPad->SetTopMargin(0.05);
  gPad->SetLeftMargin(0.150);
  gPad->SetRightMargin(0.06);
  gPad->SetTicky(1);
  gPad->SetTickx(1); 
  hFrameProtonFluxVSEkn->Draw();
  gPad->SetLogx();
  gPad->SetLogy();    

  amsdata->grAMS02_HeliumFluxVSEkn[4]->Draw("pZ");
  amsdata->grAMS02_HeliumFluxVSEkn[40]->Draw("pZ");

  grFluxEknH_AMS02Helium[4]->Draw("c");
  grFluxEknH_AMS02Helium[40]->Draw("c");

  amsdata->grVOYAGER1_HeliumFluxVSEkn->Draw("pZ");
  grLISFluxVSEkn_Helium->Draw("c");

  ccProtonFluxEnergySpectra->Update();

  
  // ---- plot DATA/LIS ratios ----

  // AMS 78th BARTED
  TGraphErrors* grAMS02ToLIS_Proton= new TGraphErrors();
  grAMS02ToLIS_Proton->SetName("grAMS02ToLIS_Proton_T78");  
  MakeDataToModel(amsdata->grAMS02_ProtonFluxVSEkn[78], grLISFluxVSEkn_Proton, grAMS02ToLIS_Proton);

  TGraphErrors* grAMS02ToLIS_Helium= new TGraphErrors();
  grAMS02ToLIS_Helium->SetName("grAMS02ToLIS_Helium_T78");  
  MakeDataToModel(amsdata->grAMS02_HeliumFluxVSEkn[78], grLISFluxVSEkn_Helium, grAMS02ToLIS_Helium);

  // AMS PRL2015
  TGraphErrors* grAMSPRL2015ToLIS_Proton= new TGraphErrors();
  grAMSPRL2015ToLIS_Proton->SetName("grAMSPRL2015ToLIS_Proton_T78");  
  MakeDataToModel(crdata->grAMS02_ProtonFluxVSEkn, grLISFluxVSEkn_Proton, grAMSPRL2015ToLIS_Proton);
  grAMSPRL2015ToLIS_Proton->SetMarkerColor(kRed+2);
  grAMSPRL2015ToLIS_Proton->SetLineColor(kRed+2);
  
  // PAMELA  
  TGraphErrors* grPAMELAToLIS_Proton= new TGraphErrors();
  grPAMELAToLIS_Proton->SetName("grPAMELAToLIS_Proton_T42");  
  MakeDataToModel(crdata->grPAMELA_ProtonFluxVSEkn[42], grLISFluxVSEkn_Proton, grPAMELAToLIS_Proton);

  
  // ---- frames ----
  TH2D* hFrameDataToLIS= new TH2D("hFrameDataToLIS","hFrameDataToLIS",200, 0.08, 1000.0, 200, 0.15, 1.25);
  SetStyleHistoVSEkn(hFrameDataToLIS);
  hFrameDataToLIS->GetYaxis()->SetNdivisions(506);
  hFrameDataToLIS->GetXaxis()->SetTitle("kinetic energy (GeV/n)");
  hFrameDataToLIS->GetYaxis()->SetTitle("DATA / LIS ");

  // ---- plot fluxes ----
  TCanvas* ccDataToLIS= new TCanvas("ccDataToLIS","DATA / LIS ratios", 1200, 450 );
  ccDataToLIS->Divide(2,1,0.01,0.01);

  ccDataToLIS->cd(1);  
  gPad->SetBottomMargin(0.15);
  gPad->SetTopMargin(0.05);
  gPad->SetLeftMargin(0.150);
  gPad->SetRightMargin(0.06);
  gPad->SetTicky(1);
  gPad->SetTickx(1); 
  hFrameDataToLIS->Draw();
  gPad->SetLogx();
  grPAMELAToLIS_Proton->Draw("pZ");
  grAMS02ToLIS_Proton->Draw("pZ");
  grAMSPRL2015ToLIS_Proton->Draw("pZ");
  gPad->SetGridy();
  
  ccDataToLIS->cd(2);  
  gPad->SetBottomMargin(0.15);
  gPad->SetTopMargin(0.05);
  gPad->SetLeftMargin(0.150);
  gPad->SetRightMargin(0.06);
  gPad->SetTicky(1);
  gPad->SetTickx(1); 
  hFrameDataToLIS->Draw();
  gPad->SetLogx();
  grAMS02ToLIS_Helium->Draw("pZ");
  gPad->SetGridy();

  
  
  // --- the end ---
  theApp.Run();
  return 1;


  
  // ---- draw energy spectra at three epochs ----
  // ---- index ----
  int indT1= 0;
  int indT2= 35; // 21;
  int indT3= 43; // NOT USED
  int indTA= 0; // AMS02 NOT USED
  int indTB= 0; // BESS
  int indTS= 18; // SOHO LAST!

  // ---- time ----
  double T1 = crdata->xTimePAMELA_ProtonFlux[indT1]; // for BESS-Polar-I
  double T2 = crdata->xTimePAMELA_ProtonFlux[indT2]; // for BESS-Polar-II ind=1 | PAMELA ind=18
  double T3 = crdata->xTimePAMELA_ProtonFlux[indT3]; // for PAMELA T-index=43
  double TA = crdata->xTimeAMS02_ProtonFlux[0]; // one point
  double TB = crdata->xTimeBESS_ProtonFlux[0]; 
  double TS = crdata->xTimeSOHO_ProtonFlux[indTS]; 


  // ----- define legend ----
  crdata->grPAMELA_ProtonFluxVSEkn[indT2]->SetMarkerStyle(26);
  crdata->grPAMELA_ProtonFluxVSEkn[indT2]->SetMarkerColor(kMagenta-5);
  TLegend* legProtonData= new TLegend(0.21,0.22,0.56,0.47);
  legProtonData->SetBorderSize(0);
  legProtonData->SetTextSize(0.038);
  legProtonData->SetTextFont(42);
  legProtonData->SetFillColor(0);
  legProtonData->SetFillStyle(0);
  legProtonData->SetNColumns(1);
  legProtonData->AddEntry(crdata->grVOYAGER1_ProtonFluxVSEkn,"Voyager-1 ","p");
  legProtonData->AddEntry(crdata->grPAMELA_ProtonFluxVSEkn[indT1],"PAMELA July 2006","p");
  legProtonData->AddEntry(crdata->grPAMELA_ProtonFluxVSEkn[indT2],"PAMELA July 2008","p");
  legProtonData->AddEntry(crdata->grSOHO_ProtonFluxVSEkn[1],"EPHIN/SOHO 2014","p");
  //legProtonData->AddEntry(crdata->grBESS_ProtonFluxVSEkn[1],"BESS-Polar ","p");


  // ---- plot fluxes ----
  TCanvas* ccProtonFluxEnergySpectraBis= new TCanvas(Form("ccProtonFluxEnergySpectr_N%d",2), Form("Proton Flux Energy Spectrum | N%d",0), 600,580 );

  ccProtonFluxEnergySpectraBis->cd();  
  gPad->SetBottomMargin(0.15);
  gPad->SetTopMargin(0.05);
  gPad->SetLeftMargin(0.150);
  gPad->SetRightMargin(0.06);
  
  gPad->SetTicky(1);
  gPad->SetTickx(1); 
  hFrameProtonFluxVSEkn->Draw();
  gPad->SetLogx();
  gPad->SetLogy();    
  legProtonData->Draw();

  crdata->grVOYAGER1_ProtonFluxVSEkn->Draw("pZ");
  // crdata->grBESS_ProtonFluxVSEkn[1]->Draw("pZ"); // come si accorda com pamela T2?
  crdata->grPAMELA_ProtonFluxVSEkn[indT1]->Draw("pZ");
  crdata->grPAMELA_ProtonFluxVSEkn[indT2]->Draw("pZ"); // 23 o 25???
  //crdata->grPAMELA_ProtonFluxVSEkn[indT3]->Draw("pZ");
  // crdata->grAMS02_ProtonFluxVSEkn->Draw("pZ");
  // crdata->grBESS_ProtonFluxVSEkn[indTB]->Draw("pZ"); 
  crdata->grSOHO_ProtonFluxVSEkn[indTS]->Draw("pZ");

  grLISFluxVSEkn_Proton->Draw("c");
  grFluxEknH_PAMELA[indT1]->Draw("c");
  grFluxEknH_PAMELA[indT2]->Draw("c");
  grFluxEknH_SOHO[indTS]->Draw("c");

  // grFluxEknH_BESS[indTB]->Draw("c");

  // grFluxEknH_PAMELA[indT3]->Draw("c");
  // grFluxEknH_AMSPRL2015[indTA]->Draw("c");
  ccProtonFluxEnergySpectraBis->Update();



  

  // ---- LAST: plot fluxes after FF-modulation with NM-Phi! ----
  // NT: NOT USED. WE USE DrawForceFieldFitsToPamela which makes this! 
  if(1 == 0){
    // 1. get phi from NEWK NMs
    
    TFile* inFileNMPhi= new TFile("grPhiVSTime_Matisse_Correct.root");
    inFileNMPhi->cd();
    
    // NT: MA USOSKIN è UN PO' MEGLIO
    TGraphErrors* grNMPhiVSTime= (TGraphErrors*)inFileNMPhi->Get("grPhiVSTime_NEWK"); // MV
    
    // phi in GV
    double PHI_T1= (1.e-3)*grNMPhiVSTime->Eval(T1);
    double PHI_T2= (1.e-3)*grNMPhiVSTime->Eval(T2);
    double PHI_TS= (1.e-3)*grNMPhiVSTime->Eval(TS);
    inFileNMPhi->Close();
    
    // 2. make NM-driven modulated graphs [ NTDEC2017: CAMBIARE IL LIS come fatto SOPRA!!!!!]
    TGraph* grNMDriven_T1=new TGraph("../SimulatedData/LIS_Proton.dat");
    grNMDriven_T1->SetName("grNMDriven_T1");
    grNMDriven_T1->SetLineColor(kGreen+1);
    grNMDriven_T1->SetLineWidth(2);
    ModulateGraph(grNMDriven_T1, PHI_T1, 1, 1);
    
    TGraph* grNMDriven_T2=new TGraph("../SimulatedData/LIS_Proton.dat");
    grNMDriven_T2->SetName("grNMDriven_T2");
    grNMDriven_T2->SetLineColor(kGreen+1);
    grNMDriven_T2->SetLineWidth(2);
    ModulateGraph(grNMDriven_T2, PHI_T2, 1, 1);
    
    TGraph* grNMDriven_TS=new TGraph("../SimulatedData/LIS_Proton.dat");
    grNMDriven_TS->SetName("grNMDriven_TS");
    grNMDriven_TS->SetLineColor(kGreen+1);
    grNMDriven_TS->SetLineWidth(2);
    ModulateGraph(grNMDriven_TS, PHI_TS, 1, 1);
    
    
    // 3. plot 
    ccProtonFluxEnergySpectraBis->cd();
    grNMDriven_T1->Draw("c");
    grNMDriven_T2->Draw("c");
    grNMDriven_TS->Draw("c");
    
  }
  // ------- end plot NM-driven fluxes ----------
  
  
  ccProtonFluxEnergySpectraBis->Update();




  // ---- store a few FF-modulated energy spectra ---
  if(kForceFieldSpectraSTORE){
    grFluxEknH_PAMELA[indT1]->SetName("grProtonFluxVSEkn_PAMELA2006");
    grFluxEknH_PAMELA[indT2]->SetName("grProtonFluxVSEkn_PAMELA2008");
    grFluxEknH_SOHO[indTS]->SetName("grProtonFluxVSEkn_SOHO2014");

    TFile* outFileFF= new TFile("grTestFFNC_N2.root","RECREATE");
    outFileFF->cd();
    grFluxEknH_PAMELA[indT1]->Write();
    grFluxEknH_PAMELA[indT2]->Write();
    grFluxEknH_SOHO[indTS]->Write();
    outFileFF->Close();
  }



  


  // ---- store all PHI reconstructions ----
  if(kSTORE){
    TFile* outFile= new TFile("grPhiVSTime_BestFitToCRData.root","RECREATE");
    outFile->cd();
    grPHIvsTime_PAMELA->Write();
    grPHIvsTime_SOHO->Write();
    grPHIvsTime_BESS->Write();
    grPHIvsTime_AMSPRL2015->Write();
    outFile->Write();
    outFile->Close();
  }


  theApp.Run();
  return 1;

}

// BOTTOM


// modulation applied to the same input graph
void ModulateGraph(TGraph* gr, double Phi, int Z, int A){
  int NP= gr->GetN();
  double* xEkn= gr->GetX();
  double* xFlux= gr->GetY();
  SolarModulation(xEkn, xFlux, NP, Z, A, Phi, 1); // log-log
  for(int ee=0;ee<NP;ee++){
    gr->SetPoint(ee, xEkn[ee], xFlux[ee]);
  }
}


// from LIS graph, modulation applied to MOD graph
void ModulateGraph(TGraph* grLIS, TGraph* grMOD, double Phi, int Z, int A){
  int NP= grLIS->GetN();

  // NO: if we modify arrays xEkn/xFlux, the graph gets modified
  // double* xEkn= grLIS->GetX();
  // double* xFlux= grLIS->GetY();

  // YES: make a replica of arrays
  double* xEkn= new double[NP];
  double* xFlux= new double[NP];
  for(int ee=0;ee<NP;ee++){
    xEkn[ee]= grLIS->GetX()[ee];
    xFlux[ee]=grLIS->GetY()[ee];
  }

  SolarModulation(xEkn, xFlux, NP, Z, A, Phi, 1); // log-log

  // reset MOD graph 
  for(int ee=grMOD->GetN()-1;ee>=0;ee--) grMOD->RemovePoint(ee);

  for(int ee=0;ee<NP;ee++){
    grMOD->SetPoint(ee, xEkn[ee], xFlux[ee]);
  }

}



void SolarModulation(double *Ekin, double *cr_density, int NEKN, int z, int a, double phi, int kOpt){
  if(phi<0.05) return; // no modulation

   int j;
   double *density,B,C,T,y;

   // NT CLARIFY!!
   //double Mp=939.; // Mp = (mp+mn)/2 ~939 MeV [o forse  0.939 GeV?????????????]
   double Mp=0.939; // Mp = (mp+mn)/2 ~939 MeV [o forse  0.939 GeV?????????????]

   density=new double[NEKN];

   int key= kOpt;
   if(key<0 || key>2)key=0;
 
   for(int i=0; i<NEKN; i++)
   {

     T = Ekin[i]+fabs(z)*phi/fabs(a); // NT: quindi l'input è Ekin, NON Ekn!!! Giusto?
     // T = Ekin[i]+abs(z)*phi/a; // NT: quindi l'input è Ekin, NON Ekn!!! Giusto?
     for(j=0; j<NEKN; j++) if(T<Ekin[j]) break;
     if(j==NEKN) { density[i] = cr_density[i]; break; }
      

     if(key==0){ // linear interpolation in energy
       y = cr_density[j-1] +(T-Ekin[j-1]) *(cr_density[j]-cr_density[j-1]) /(Ekin[j]-Ekin[j-1]);
     }
     
     if(key==1){ // log-log interpolation
       double logy= log10(cr_density[j-1]) + (log10(T)- log10(Ekin[j-1])) *( log10(cr_density[j])-log10(cr_density[j-1])) /(log10(Ekin[j])-log10(Ekin[j-1]));
       y=TMath::Power(10, logy);
     }
     
     if(key==2){ // power-law interpolation
       if(cr_density[j-1]*cr_density[j]>0.)
	 {
	   B = log(cr_density[j-1]/cr_density[j])/log(Ekin[j-1]/Ekin[j]);
	   C = cr_density[j]*pow(Ekin[j],-B);
	   y = C*pow(T,B);
	 }
     }
     
     // ---- IMOS: remove last factor if working with flux ----
     // density[i]= y * Ekin[i] *(Ekin[i]+2*Mp) /T /(T+2*Mp) * pow(Ekin[i]/T,2);

     // ---- NTSEP2017: SO WE REMOVED LAST FACTOR ----
     density[i]= y * Ekin[i] *(Ekin[i]+2*Mp) /T /(T+2*Mp); // *pow(Ekin[i]/T,2);
   }
   
   // APPLY NORMALIZATION FACTOR [NTDEC2017: NormProton NOW REMOVED from HERE!!]
   for(j=0;j<NEKN;j++) cr_density[j] = density[j];
   delete[] density;
}



void FCN(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag){
  double Phi= par[0]; // GV
  int iExperiment   = (int)par[1]; // PAMELA / SOHO / BESS / AMS02
  int iTime         = (int)par[2];
  f= GetChiSquare(Phi, iTime, iExperiment);
}



// chisquare at fixed epoch
double GetChiSquare(double Phi, int iTime, int iExperiment){

  double ChiSquare = 0.;
  double Emin = EMIN;
  double Emax = EMAX;

  int Z=1;
  int A=1;
  if(iExperiment==5){ Z=2; A=4; } // AMS-He: PER IL MOMENTO NO

  // ---- modulate ----
  TGraph* grMOD= new TGraph();
  if(Z==1 && A==1) ModulateGraph(grLISFluxVSEkn_Proton,grMOD, Phi, Z, A);
  if(Z==2 && A==4) ModulateGraph(grLISFluxVSEkn_Helium,grMOD, Phi, Z, A);

    
  // --- EXP1: PAMELA ----
  if( iExperiment==1){
    NDF_PAMELA[iTime]=0;
    for(int ee=0;ee<crdata->nEknPAMELA_Proton;ee++){ // ekn loop

      // check RANGE
      double Energy= crdata->grPAMELA_ProtonFluxVSEkn[iTime]->GetX()[ee];
      if(Energy<Emin || Energy>Emax) continue; 

      
      // get DATA 
      double Flux  = crdata->grPAMELA_ProtonFluxVSEkn[iTime]->GetY()[ee];
      double eFlux = crdata->grPAMELA_ProtonFluxVSEkn[iTime]->GetEY()[ee];
      
      // get CALCULATIONS
      double FluxModel = grMOD->Eval(Energy);

      // calc CHISQUARE      
      double Delta = (Flux - FluxModel)/eFlux;
      ChiSquare += (Delta*Delta);
      //cout<<"Time: "<<iTime<<"     Chi: "<<ChiSquare<<"     PHI: "<<Phi<<endl;      
    
      NDF_PAMELA[iTime]++;
    }
    if( bestChi2_PAMELA[iTime] > ChiSquare ) bestChi2_PAMELA[iTime] = ChiSquare;
  }


  // --- EXP2: SOHO ----
  if( iExperiment==2){
    NDF_SOHO[iTime]=0;
    for(int ee=0;ee<crdata->nEknSOHO_Proton;ee++){ // ekn loop
      
      // check RANGE
      double Energy= crdata->grSOHO_ProtonFluxVSEkn[iTime]->GetX()[ee];
      if(Energy<Emin || Energy>Emax) continue; 
      
      // get DATA 
      double Flux  = crdata->grSOHO_ProtonFluxVSEkn[iTime]->GetY()[ee];
      double eFlux = crdata->grSOHO_ProtonFluxVSEkn[iTime]->GetEY()[ee];
      
      // get CALCULATIONS
      double FluxModel = grMOD->Eval(Energy);
      
      // calc CHISQUARE      
      double Delta = (Flux - FluxModel)/eFlux;
      ChiSquare += Delta*Delta;

      NDF_SOHO[iTime]++;
    }
    if( bestChi2_SOHO[iTime] > ChiSquare ) bestChi2_SOHO[iTime] = ChiSquare;
  }


  // --- EXP3: BESS ----
  if( iExperiment==3){
    NDF_BESS[iTime]=0;
    for(int ee=0;ee<crdata->nEknBESS_Proton;ee++){ // ekn loop
      
      // check RANGE
      double Energy= crdata->grBESS_ProtonFluxVSEkn[iTime]->GetX()[ee];
      if(Energy<Emin || Energy>Emax) continue; 
      
      // get DATA 
      double Flux  = crdata->grBESS_ProtonFluxVSEkn[iTime]->GetY()[ee];
      double eFlux = crdata->grBESS_ProtonFluxVSEkn[iTime]->GetEY()[ee];
      
      // get CALCULATIONS
      double FluxModel = grMOD->Eval(Energy);
      
      // calc CHISQUARE      
      double Delta = (Flux - FluxModel)/eFlux;
      ChiSquare += Delta*Delta;

      NDF_BESS[iTime]++;
    }
    if( bestChi2_BESS[iTime] > ChiSquare ) bestChi2_BESS[iTime] = ChiSquare;
  }



  // --- EXP4: AMS02 PROTON ----
  if( iExperiment==4){
    NDF_AMS02[iTime]=0;

    for(int ee=0;ee<amsdata->nEknAMS02_Proton;ee++){ // ekn loop
      
      double Energy= amsdata->grAMS02_ProtonFluxVSEkn[iTime]->GetX()[ee];
      if(Energy<Emin || Energy>Emax) continue; 
      
      // --- re-check range for AMS02 ---
      // if(Energy>60.0) continue; // Range should be R=1-60 GV [actual range 1-100 GV]
      
      // get DATA 
      double Flux  = amsdata->grAMS02_ProtonFluxVSEkn[iTime]->GetY()[ee];
      double eFlux = amsdata->grAMS02_ProtonFluxVSEkn[iTime]->GetEY()[ee];
      
      // get CALCULATIONS
      double FluxModel = grMOD->Eval(Energy);
      
      // calc CHISQUARE      
      double Delta = (Flux - FluxModel)/eFlux;
      ChiSquare += Delta*Delta;

      NDF_AMS02[iTime]++;
    }
    if( bestChi2_AMS02[iTime] > ChiSquare ) bestChi2_AMS02[iTime] = ChiSquare;
  }


  // --- EXP5: AMS02 HELIUM ---- 
  if( iExperiment==5){
    NDF_AMS02Helium[iTime]=0;

    for(int ee=0;ee<amsdata->nEknAMS02_Helium;ee++){ // ekn loop
      
      double Energy= amsdata->grAMS02_HeliumFluxVSEkn[iTime]->GetX()[ee];
      if(Energy<Emin || Energy>Emax) continue; 
      
      // --- re-check range for AMS02 ---
      // if(Energy>60.0) continue; // Range should be R=1-60 GV [actual range 1-100 GV]
      
      // get DATA 
      double Flux  = amsdata->grAMS02_HeliumFluxVSEkn[iTime]->GetY()[ee];
      double eFlux = amsdata->grAMS02_HeliumFluxVSEkn[iTime]->GetEY()[ee];
      // if( !(eFlux>0) ) continue; // NT SOME PROBLEMS WITH Helium ERRORS!
    
      // get CALCULATIONS
      double FluxModel = grMOD->Eval(Energy);

	    
      // calc CHISQUARE      
      double Delta = (Flux - FluxModel)/eFlux;
      ChiSquare += Delta*Delta;

      /// DUMP
      // cout<<"T: "<<iTime<< "      E: "<<ee<<  "       DATA: "<< Flux<<"    Err: "<<eFlux<< "          MODEL: "<<FluxModel<<"    Delta: "<<Delta<<"    Chi2: "<<ChiSquare<<endl;
      
      
      NDF_AMS02Helium[iTime]++;
    }
    if( bestChi2_AMS02Helium[iTime] > ChiSquare ) bestChi2_AMS02Helium[iTime] = ChiSquare;
  }

  


  // --- EXP9: AMSPRL2015 PRL-2015 ONE POINT ---- [no iTime index!]
  if( iExperiment==9){
    NDF_AMSPRL2015[iTime]=0;
    for(int ee=0;ee<crdata->nEknAMS02_Proton;ee++){ // ekn loop
      
      // check RANGE
      double Energy= crdata->grAMS02_ProtonFluxVSEkn->GetX()[ee];
      if(Energy<Emin || Energy>Emax) continue; 
      
      // get DATA 
      double Flux  = crdata->grAMS02_ProtonFluxVSEkn->GetY()[ee];
      double eFlux = crdata->grAMS02_ProtonFluxVSEkn->GetEY()[ee];
      
      // get CALCULATIONS
      double FluxModel = grMOD->Eval(Energy);
      
      // calc CHISQUARE      
      double Delta = (Flux - FluxModel)/eFlux;
      ChiSquare += Delta*Delta;

      NDF_AMSPRL2015[iTime]++;
    }
    if( bestChi2_AMSPRL2015[iTime] > ChiSquare ) bestChi2_AMSPRL2015[iTime] = ChiSquare;
  }


  delete grMOD;
  return ChiSquare;
}








void SetStyleHistoVSEkn(TH2D* hh){
  hh->SetTitle(0);
  hh->SetLabelFont(42,"X");
  hh->SetLabelFont(42,"Y");
  hh->SetTitleFont(42,"X");
  hh->SetTitleFont(42,"Y");
  hh->GetXaxis()->SetTitle("kinetic energy (GeV)");


  //hh->GetYaxis()->CenterTitle();
  hh->GetXaxis()->SetNdivisions(513);
  hh->GetYaxis()->SetNdivisions(508);

  hh->GetXaxis()->SetLabelSize(0.05);
  hh->GetYaxis()->SetLabelSize(0.05);
  hh->GetYaxis()->SetTitleSize(0.05);
  hh->GetXaxis()->SetTitleSize(0.05);

  hh->GetXaxis()->SetTitleOffset(1.40);
  hh->GetYaxis()->SetTitleOffset(1.40);

  hh->GetXaxis()->SetLabelOffset(0.00);
  hh->GetYaxis()->SetLabelOffset(0.01);


}


void SetStyleHistoVSTime(TH2D* hh){
  hh->SetTitle(0);
  hh->GetXaxis()->SetTitle(0);
  //hh->GetYaxis()->SetTitle("J(R) [ GV^{ -1} m^{ -2} s^{ -1} sr^{ -1} ]");
  hh->GetXaxis()->SetTimeDisplay(1);
  hh->GetXaxis()->SetTimeFormat("#splitline{%Y}{%b,%d}%F1970-01-01 00:00:00s0");

  hh->GetXaxis()->SetNdivisions(513);
  hh->GetXaxis()->SetLabelSize(0.05);
  hh->GetYaxis()->SetNdivisions(508);
  hh->GetYaxis()->SetLabelSize(0.05);
  hh->GetYaxis()->SetTitleSize(0.05);
  hh->GetYaxis()->CenterTitle();
  hh->GetZaxis()->SetLabelSize(0.05);
  //hh->GetXaxis()->SetLabelOffset(0.06);
  //hh->GetYaxis()->SetTitleOffset(1.30);
}

void SetTimeGraphStyle(TGraphErrors* gr){
    gr->SetTitle(0);
    gr->GetXaxis()->SetTitle(0);
    gr->GetYaxis()->SetTitle("#phi (GV) ");
    gr->GetXaxis()->SetTimeDisplay(1);
    gr->GetXaxis()->SetTimeFormat("#splitline{%Y}{%b,%d}%F1970-01-01 00:00:00s0");    
    gr->GetXaxis()->SetNdivisions(513);
    gr->GetXaxis()->SetLabelSize(0.06);
    gr->GetXaxis()->SetLabelOffset(0.06);
    gr->GetYaxis()->SetNdivisions(508);
    gr->GetYaxis()->SetLabelSize(0.06);
    gr->GetYaxis()->SetTitleOffset(0.65);
    gr->GetYaxis()->SetTitleSize(0.06);
    gr->GetYaxis()->CenterTitle();
}


// // ---- scale graph by input factor ----
void ScaleGraph(TGraph* gr, double ScaleFactor){
   for(int ii=0;ii<gr->GetN();ii++){
     double x, y;
     gr->GetPoint(ii, x, y);
     y *= ScaleFactor;
     gr->SetPoint(ii, x, y);
   }
}




void FlattenGraph(TGraph* gr, double OldFlatIndex, double NewFlatIndex){
  int NP= gr->GetN();
  for(int ii=0;ii<NP;ii++){
    double Xxx, Val; // ErrX, ErrY;
    gr->GetPoint(ii, Xxx, Val);
    // de-flattening AND re-flattening
    Val  *= pow( Xxx, NewFlatIndex - OldFlatIndex );
    gr->SetPoint(ii, Xxx, Val);  
  }
}


void MakeDataToModel(TGraphErrors* grData, TGraph* grModel, TGraphErrors* grDataToModel){
  for(int ii=0;ii<grData->GetN();ii++){
    double dX= grData->GetX()[ii];
    double dY= grData->GetY()[ii];
    double eY= grData->GetEY()[ii];
    double mY= grModel->Eval(dX);
    double D2M = dY/mY;
    double eD2M= eY/mY;

    grDataToModel->SetPoint(ii, dX, D2M);
    grDataToModel->SetPointError(ii, 0.0, eD2M);
  }

  grDataToModel->SetMarkerStyle( grData->GetMarkerStyle() );
  grDataToModel->SetMarkerSize( grData->GetMarkerSize() );
  grDataToModel->SetMarkerColor( grData->GetMarkerColor() );
  grDataToModel->SetLineColor( grData->GetLineColor() );
  grDataToModel->SetLineWidth( grData->GetLineWidth() );
}




  // ---- TEST PLOT ----
  /*

  TGraph* grMOD= new TGraph();
  grMOD->SetLineColor(kPink+2);
  grMOD->SetLineStyle(7);
  grMOD->SetLineWidth(3);
  ModulateGraph(grLISFluxVSEkn_Proton,grMOD, 0.4, 1, 1);


  TH2D* hFrameProtonFluxVSEkn= new TH2D("hFrameProtonFluxVSEkn","hFrameProtonFluxVSEkn",200, 0.1, 30.0, 200, 5.e-1, 3.e+4);
  SetStyleHistoVSEkn(hFrameProtonFluxVSEkn);
  hFrameProtonFluxVSEkn->GetXaxis()->SetTitle("kinetic energy (GeV)");
  hFrameProtonFluxVSEkn->GetYaxis()->SetNdivisions(506);
  hFrameProtonFluxVSEkn->GetYaxis()->SetTitle("J ( GeV^{ -1} m^{ -2} s^{ -1} sr^{ -1} )");


  TCanvas* ccProtonFluxEnergySpectra= new TCanvas(Form("ccProtonFluxEnergySpectr_N%d",0), Form("Proton Flux Energy Spectrum | N%d",0), 600,580 );

  ccProtonFluxEnergySpectra->cd();  
  gPad->SetBottomMargin(0.15);
  gPad->SetTopMargin(0.05);
  gPad->SetLeftMargin(0.150);
  gPad->SetRightMargin(0.06);
  
  gPad->SetTicky(1);
  gPad->SetTickx(1); 
  hFrameProtonFluxVSEkn->Draw();
  gPad->SetLogx();
  gPad->SetLogy();    

  crdata->grVOYAGER1_ProtonFluxVSEkn->Draw("pZ");
  // crdata->grBESS_ProtonFluxVSEkn[1]->Draw("pZ"); // come si accorda com pamela T2?
  crdata->grPAMELA_ProtonFluxVSEkn[20]->Draw("pZ");
  grLISFluxVSEkn_Proton->Draw("c");
  grMOD->Draw("c");
  ccProtonFluxEnergySpectra->Update();    
  
  theApp.Run();
  return 1;
*/
