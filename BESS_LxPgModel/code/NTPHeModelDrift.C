#include "TROOT.h"
#include "TMath.h"
#include "TGraphErrors.h"
#include "TGraph.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TAxis.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TApplication.h"
#include "TMinuit.h"
#include <iostream>
#include "NTCosmicRayData.cpp"
#include "NTAMS02Data.cpp"


// --------- NTPHeModel -----------
// -- Make K0-Fits to CR proton data: PAMELA-SOHO-BESS-AMS02
// -- Make Predictions to CR helium data at same epochs
// -- Mix 3He + 4He templates. Vary the composition
// -- Compute the p/He ratio
// ---------------------------------

// TOPTOP

int kLISModel = 0;
bool kReNorm  = false;
bool kDriftON = false;

// ---- constants ----
const int NITERATIONS = 10000;

double EMIN= 0.02; //0.07
double EMAX= 20.; 
int thisRIG1 = 1; // rig bin for p/He plots Low-R
int thisRIG2 = 4; // rig bin for p/He plots High-R
int thisRIG3 = 8; // rig bin for p/He plots High-R

// --- histo limits ---
// double hLnEmin;
// double hLnEmax;
// double hLnRmin;
// double hLnRmax;
// double hXiDmin;
// double hXiDmax;
// double hLnK0min;
// double hLnK0max;


static const int nTimePAMELA_Proton     = 47;
static const int nTimeSOHO_Proton       = 20;
static const int nTimeBESS_Proton       = 2; // Polar I/II
static const int nTimeBESSTeV_Proton    = 1; 
static const int nTimeBESS00_Proton     = 1; 
static const int nTimeAMS02_Proton      = 81;
static const int nTimeAMSPRL2015_Proton = 1;


// ---- experimental index ----
const int iExp_PAMELA      = 1;
const int iExp_SOHO        = 2;
const int iExp_BESS        = 3;
const int iExp_BESSTeV     = 4;
const int iExp_BESS00      = 5;
const int iExp_AMS02       = 6;
const int iExp_AMSPRL2015  = 7;

// ---- renormalization factor ----
double  startN0_PAMELA   = 0.955; // PAMELA
double  startN0_SOHO     = 0.955; // SOHO
double  startN0_BESSPolar= 0.96;  // BESSPolar
double  startN0_BESSTeV  = 1.0;   // BESSTeV
double  startN0_BESS00   = 1.0;   // BESS00
double  startN0_AMS02    = 1.0;   // AMS02


static const int nTimeAMS02_Helium = 81;

double bestChi2_PAMELA[nTimePAMELA_Proton];
double bestChi2_SOHO[nTimeSOHO_Proton];
double bestChi2_BESS[nTimeBESS_Proton];
double bestChi2_BESSTeV[nTimeBESSTeV_Proton];
double bestChi2_BESS00[nTimeBESS00_Proton];
double bestChi2_AMSPRL2015[nTimeAMSPRL2015_Proton];
double bestChi2_AMS02[nTimeAMS02_Proton];
// double bestChi2_AMS02Helium[nTimeAMS02_Helium];

int NDF_PAMELA[nTimePAMELA_Proton];
int NDF_SOHO[nTimeSOHO_Proton];
int NDF_BESS[nTimeBESS_Proton];
int NDF_BESSTeV[nTimeBESSTeV_Proton];
int NDF_BESS00[nTimeBESS00_Proton];
int NDF_AMSPRL2015[nTimeAMSPRL2015_Proton];
int NDF_AMS02[nTimeAMS02_Proton];
// int NDF_AMS02Helium[nTimeAMS02_Helium];




// --- graph style ---
void SetStyleHistoVSEkn(TH2F* hh);
void SetStyleHistoVSTime(TH2F* hh);
void SetTimeGraphStyle(TGraphErrors* gr);

void SetGraphK0Style(); 
void SetGraphN0Style(); 
void SetGraphXiDStyle();


// --- graph manipulation ---
void FlattenGraph(TGraph* gr, double OldFlatIndex, double NewFlatIndex); // for LIS
void ScaleGraph(TGraph* gr, double ScaleFactor);
void MakeDataToModel(TGraphErrors* grData, TGraph* grModel, TGraphErrors* grDataToModel);


// --- minimization ---
double GetChiSquare(double K0, int iTime, int iExperiment, double XiD, double Norm);
void FCN(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);
double Interpolation(TH3F* hh, double x, double y, double z);
double Interpolation(TH2F* hh, double x, double y);

// ---- external objects ----

// --- cosmic-ray data ---
NTCosmicRayData* crdata;
NTAMS02Data* amsdata;

// --- model output ---

// 3D HISTOS: DRIFT ON
TH3F* hLn_Flux_vs_XiD_vs_KScale_vs_KEnergy_Z1_A1; // PROTON vs Ekn: for FIT
TH3F* hLn_Flux_vs_XiD_vs_KScale_vs_KEnergy_Z2_A3; // He3 vs Ekn: for PLOT
TH3F* hLn_Flux_vs_XiD_vs_KScale_vs_KEnergy_Z2_A4; // He4 vs Ekn: for PLOT
TH3F* hLn_Flux_vs_XiD_vs_KScale_vs_Rigidity_Z1_A1; // vs Rig: for PREDICTIONS
TH3F* hLn_Flux_vs_XiD_vs_KScale_vs_Rigidity_Z2_A3;
TH3F* hLn_Flux_vs_XiD_vs_KScale_vs_Rigidity_Z2_A4;

// 2D HISTOS: DRIFT OFF [NTJAN2018]
TH2F* hLn_Flux_vs_KScale_vs_KEnergy_Z1_A1; // PROTON vs Ekn: for FIT
TH2F* hLn_Flux_vs_KScale_vs_KEnergy_Z2_A3; // He3 vs Ekn: for PLOT
TH2F* hLn_Flux_vs_KScale_vs_KEnergy_Z2_A4; // He4 vs Ekn: for PLOT
TH2F* hLn_Flux_vs_KScale_vs_Rigidity_Z1_A1; // vs Rig: for PREDICTIONS
TH2F* hLn_Flux_vs_KScale_vs_Rigidity_Z2_A3;
TH2F* hLn_Flux_vs_KScale_vs_Rigidity_Z2_A4;



// aggiungere TH2 per il caso drift-off?

// --- fit results ---

// K0 diffusion scaling
TGraphErrors* grK0vsTime_PAMELA;
TGraphErrors* grK0vsTime_SOHO;
TGraphErrors* grK0vsTime_BESS;
TGraphErrors* grK0vsTime_BESSTeV;
TGraphErrors* grK0vsTime_BESS00;
TGraphErrors* grK0vsTime_AMSPRL2015;
TGraphErrors* grK0vsTime_AMS02;

// XiD drift velocity component
TGraphErrors* grXiDvsTime_PAMELA;
TGraphErrors* grXiDvsTime_SOHO;
TGraphErrors* grXiDvsTime_BESS;
TGraphErrors* grXiDvsTime_BESSTeV;
TGraphErrors* grXiDvsTime_BESS00;
TGraphErrors* grXiDvsTime_AMSPRL2015;
TGraphErrors* grXiDvsTime_AMS02;


// N0 flux normalization
TGraphErrors* grN0vsTime_PAMELA;
TGraphErrors* grN0vsTime_SOHO;
TGraphErrors* grN0vsTime_BESS;
TGraphErrors* grN0vsTime_BESSTeV;
TGraphErrors* grN0vsTime_BESS00;
TGraphErrors* grN0vsTime_AMSPRL2015;
TGraphErrors* grN0vsTime_AMS02;


int main(int argc, char **argv){ //TOP

  TApplication theApp("App",&argc,argv);
  gROOT->ProcessLine(".x $NTBASEDIR/THEORY/00crdata/style/MyRootStyle.C");  
  //gROOT->ProcessLine(".L ../Matisse/lib/libNTMatisseLib.so");
  
  // ---- set cosmic-ray data ----
  crdata= new NTCosmicRayData();
  crdata->SetProtonData();
  
  // ---- set AMS02 data ----
  amsdata= new NTAMS02Data();
  amsdata->SetProtonData();
  amsdata->SetHeliumData();
  
  // --- set BESS data --- NO!
  // BESSTeV & BESS00 now added under CRData
  // bessdata = new NTBESSData();
  // bessdata->SetProtonData();

  
  // ---- set model grids ----
  TFile* gridFile;
  if(kDriftON){
    gridFile= new TFile(Form("../outroot/hPHeModel_DriftON_LIS%d.root",kLISModel),"READ");

  // OLD
  // TFile* gridFile= new TFile(Form("../hHeIsoFluxTablesDrift_LIS%d.root",kLISModel),"READ");
  // TFile* gridFile= new TFile(Form("../outroot/hHeIsoFluxTables_Norm.root"),"READ");
  
  hLn_Flux_vs_XiD_vs_KScale_vs_KEnergy_Z1_A1= (TH3F*)gridFile->Get("hLn_Flux_vs_XiD_vs_KScale_vs_KEnergy_Z1_A1")->Clone("hhLn_Flux_vs_XiD_vs_KScale_vs_KEnergy_Z1_A1");
  hLn_Flux_vs_XiD_vs_KScale_vs_KEnergy_Z1_A1->SetDirectory(0);
  
  hLn_Flux_vs_XiD_vs_KScale_vs_KEnergy_Z2_A3= (TH3F*)gridFile->Get("hLn_Flux_vs_XiD_vs_KScale_vs_KEnergy_Z2_A3")->Clone("hhLn_Flux_vs_XiD_vs_KScale_vs_KEnergy_Z2_A3"); // A3 -> A4!!!
  hLn_Flux_vs_XiD_vs_KScale_vs_KEnergy_Z2_A3->SetDirectory(0);

  hLn_Flux_vs_XiD_vs_KScale_vs_KEnergy_Z2_A4= (TH3F*)gridFile->Get("hLn_Flux_vs_XiD_vs_KScale_vs_KEnergy_Z2_A4")->Clone("hhLn_Flux_vs_XiD_vs_KScale_vs_KEnergy_Z2_A4");
  hLn_Flux_vs_XiD_vs_KScale_vs_KEnergy_Z2_A4->SetDirectory(0);
  
  hLn_Flux_vs_XiD_vs_KScale_vs_Rigidity_Z1_A1= (TH3F*)gridFile->Get("hLn_Flux_vs_XiD_vs_KScale_vs_Rigidity_Z1_A1")->Clone("hhLn_Flux_vs_XiD_vs_KScale_vs_Rigidity_Z1_A1");
  hLn_Flux_vs_XiD_vs_KScale_vs_Rigidity_Z1_A1->SetDirectory(0);

  hLn_Flux_vs_XiD_vs_KScale_vs_Rigidity_Z2_A3= (TH3F*)gridFile->Get("hLn_Flux_vs_XiD_vs_KScale_vs_Rigidity_Z2_A3")->Clone("hhLn_Flux_vs_XiD_vs_KScale_vs_Rigidity_Z2_A3"); // A3 -> A4!!!
  hLn_Flux_vs_XiD_vs_KScale_vs_Rigidity_Z2_A3->SetDirectory(0);
  
  hLn_Flux_vs_XiD_vs_KScale_vs_Rigidity_Z2_A4= (TH3F*)gridFile->Get("hLn_Flux_vs_XiD_vs_KScale_vs_Rigidity_Z2_A4")->Clone("hhLn_Flux_vs_XiD_vs_KScale_vs_Rigidity_Z2_A4");
  hLn_Flux_vs_XiD_vs_KScale_vs_Rigidity_Z2_A4->SetDirectory(0);  
  gridFile->Close();
  }

  if(!kDriftON){ // 2D HISTOS | only LIS0 currently?
    gridFile= new TFile(Form("../outroot/hPHeModel_DriftOFF_LIS%dc_AllTH2_HeTOT.root",kLISModel),"READ");

    hLn_Flux_vs_KScale_vs_KEnergy_Z1_A1= (TH2F*)gridFile->Get("hLn_Flux_vs_KScale_vs_KEnergy_Z1_A1");
    hLn_Flux_vs_KScale_vs_KEnergy_Z1_A1->SetDirectory(0);

    hLn_Flux_vs_KScale_vs_KEnergy_Z2_A3= (TH2F*)gridFile->Get("hLn_Flux_vs_KScale_vs_KEnergy_Z2_A3");
    hLn_Flux_vs_KScale_vs_KEnergy_Z2_A3->SetDirectory(0);
    
    hLn_Flux_vs_KScale_vs_KEnergy_Z2_A4= (TH2F*)gridFile->Get("hLn_Flux_vs_KScale_vs_KEnergy_Z2_A4");
    hLn_Flux_vs_KScale_vs_KEnergy_Z2_A4->SetDirectory(0);
    
    hLn_Flux_vs_KScale_vs_Rigidity_Z1_A1= (TH2F*)gridFile->Get("hLn_Flux_vs_KScale_vs_Rigidity_Z1_A1");
    hLn_Flux_vs_KScale_vs_Rigidity_Z1_A1->SetDirectory(0);
    
    hLn_Flux_vs_KScale_vs_Rigidity_Z2_A3= (TH2F*)gridFile->Get("hLn_Flux_vs_KScale_vs_Rigidity_Z2_A3");
    hLn_Flux_vs_KScale_vs_Rigidity_Z2_A3->SetDirectory(0);
    
    hLn_Flux_vs_KScale_vs_Rigidity_Z2_A4= (TH2F*)gridFile->Get("hLn_Flux_vs_KScale_vs_Rigidity_Z2_A4");
    hLn_Flux_vs_KScale_vs_Rigidity_Z2_A4->SetDirectory(0);  
    gridFile->Close();
  }


  // ---- graph with K0-results ---- 
  grK0vsTime_PAMELA= new TGraphErrors();
  grK0vsTime_SOHO= new TGraphErrors();
  grK0vsTime_BESS= new TGraphErrors();
  grK0vsTime_BESSTeV= new TGraphErrors();
  grK0vsTime_BESS00= new TGraphErrors();
  grK0vsTime_AMS02= new TGraphErrors();
  grK0vsTime_AMSPRL2015= new TGraphErrors();
  SetGraphK0Style();

  grN0vsTime_PAMELA= new TGraphErrors();
  grN0vsTime_SOHO= new TGraphErrors();
  grN0vsTime_BESS= new TGraphErrors();
  grN0vsTime_BESSTeV= new TGraphErrors();
  grN0vsTime_BESS00= new TGraphErrors();
  grN0vsTime_AMS02= new TGraphErrors();
  grN0vsTime_AMSPRL2015= new TGraphErrors();
  SetGraphN0Style();


  grXiDvsTime_PAMELA= new TGraphErrors();
  grXiDvsTime_SOHO= new TGraphErrors();
  grXiDvsTime_BESS= new TGraphErrors();
  grXiDvsTime_BESSTeV= new TGraphErrors();
  grXiDvsTime_BESS00= new TGraphErrors();
  grXiDvsTime_AMS02= new TGraphErrors();
  grXiDvsTime_AMSPRL2015= new TGraphErrors();
  SetGraphXiDStyle();

  
  // ---- initializations ----
  for(int tt=0;tt<nTimePAMELA_Proton;tt++) bestChi2_PAMELA[tt]= 1.e+9;
  for(int tt=0;tt<nTimeSOHO_Proton;tt++) bestChi2_SOHO[tt]= 1.e+9;
  for(int tt=0;tt<nTimeBESS_Proton;tt++) bestChi2_BESS[tt]= 1.e+9;
  for(int tt=0;tt<nTimeBESSTeV_Proton;tt++) bestChi2_BESSTeV[tt]= 1.e+9;
  for(int tt=0;tt<nTimeBESS00_Proton;tt++) bestChi2_BESS00[tt]= 1.e+9;
  for(int tt=0;tt<nTimeAMSPRL2015_Proton;tt++) bestChi2_AMSPRL2015[tt]= 1.e+9;
  for(int tt=0;tt<nTimeAMS02_Proton;tt++) bestChi2_AMS02[tt]= 1.e+9;


  // ----- minimization ----
  const int NPAR = 5; // K0, iExp, iTime, XiD, NORM [ Emin-Emax from ABOVE]
  TMinuit* gMinuit = new TMinuit(NPAR);  
  gMinuit->SetPrintLevel(-1);
  gMinuit->SetFCN(FCN);

  Double_t arglist[10];
  Int_t    ierflg = 0;
  arglist[0] = 1;
  gMinuit->mnexcm("SET ERR", arglist ,1,ierflg);
  
  double minK0= 1.5;
  double maxK0= 9.; //25.;
  double stepK0=1.e-2;
  double startK0=4.5;

  double minXiD= -1.4;
  double maxXiD=  1.4;
  double stepXiD=1.e-2;
  double startXiD=0.0;

  double minN0   = 0.95;
  double maxN0   = 1.05;
  double stepN0  = 1.e-3;
  double startN0 =  1.0;


  // ---- fitting loop: PAMELA ---
  cout<<"**** Fit to PAMELA ****"<<endl;

  // Dedicated NORM for PAMELA:  startN0= 0.955;
  startN0=  startN0_PAMELA;


  for(int tt=0;tt<nTimePAMELA_Proton;tt++){
    gMinuit->DefineParameter(0, "K_{0}", startK0, stepK0, minK0, maxK0);
    gMinuit->DefineParameter(1, "Experiment", iExp_PAMELA, 0, 0, 0);
    gMinuit->DefineParameter(2, "Time Index", tt, 0, 0, 0);
    gMinuit->DefineParameter(3, "XiD", startXiD, stepXiD, minXiD, maxXiD);
    gMinuit->DefineParameter(4, "NORM", startN0, stepN0, minN0, maxN0);
    gMinuit->FixParameter(1); // FIX TO WHAT? 
    gMinuit->FixParameter(2); 
    if(!kDriftON) gMinuit->FixParameter(3); // set Drift OFF
    if(!kReNorm) gMinuit->FixParameter(4); // No Renormalization parameter

    arglist[0] = NITERATIONS;
    arglist[1] = 1.;
    gMinuit->mnexcm("MIGRAD", arglist ,2,ierflg);
    
    double bestK0= 0;
    double errK0 = 0.;
    double bestN0= 0;
    double errN0 = 0.;
    double bestXiD= 0;
    double errXiD = 0.;
    gMinuit->GetParameter(0,bestK0,errK0);
    gMinuit->GetParameter(3,bestXiD,errXiD);
    gMinuit->GetParameter(4,bestN0,errN0);
    startK0= bestK0; // update for next fit
    startXiD= bestXiD; // update for next fit

    cout <<"iTime: "<<tt<<"   K0= " << bestK0 << " +- "<<errK0<<" | "<<"    N0: "<<bestN0<<"     Chi2: "<<bestChi2_PAMELA[tt]<< " / "<<NDF_PAMELA[tt]<<endl;
    
    // ---- put results into graph ----
    grK0vsTime_PAMELA->SetPoint(tt, crdata->xTimePAMELA_ProtonFlux[tt], bestK0);
    grK0vsTime_PAMELA->SetPointError(tt, crdata->eTimePAMELA_ProtonFlux, errK0);
    grN0vsTime_PAMELA->SetPoint(tt, crdata->xTimePAMELA_ProtonFlux[tt], bestN0);
    grN0vsTime_PAMELA->SetPointError(tt, crdata->eTimePAMELA_ProtonFlux, errN0);
    grXiDvsTime_PAMELA->SetPoint(tt, crdata->xTimePAMELA_ProtonFlux[tt], bestXiD);
    grXiDvsTime_PAMELA->SetPointError(tt, crdata->eTimePAMELA_ProtonFlux, errXiD);
}


  

  // ---- fitting loop: SOHO ---
  cout<<"**** Fit to SOHO ****"<<endl;

  // Dedicated NORM for SOHO startN0= 0.955;
  startN0=  startN0_SOHO;

  
  for(int tt=0;tt<nTimeSOHO_Proton;tt++){
    gMinuit->DefineParameter(0, "K_{0}", startK0, stepK0, minK0, maxK0);
    gMinuit->DefineParameter(1, "Experiment", iExp_SOHO, 0, 0, 0);
    gMinuit->DefineParameter(2, "Time Index", tt, 0, 0, 0);
    gMinuit->DefineParameter(3, "XiD", startXiD, stepXiD, minXiD, maxXiD);
    gMinuit->DefineParameter(4, "NORM", startN0, stepN0, minN0, maxN0);
    gMinuit->FixParameter(1); // FIX TO WHAT? 
    gMinuit->FixParameter(2); 
    if(!kDriftON) gMinuit->FixParameter(3); // set Drift OFF
    if(!kReNorm) gMinuit->FixParameter(4); // No Renormalization parameter
    
    arglist[0] = NITERATIONS;
    arglist[1] = 1.;
    gMinuit->mnexcm("MIGRAD", arglist ,2,ierflg);
    
    double bestK0= 0;
    double errK0 = 0.;
    double bestN0= 0;
    double errN0 = 0.;
    double bestXiD= 0;
    double errXiD = 0.;
    gMinuit->GetParameter(0,bestK0,errK0);
    gMinuit->GetParameter(3,bestXiD,errXiD);
    gMinuit->GetParameter(4,bestN0,errN0);
    startK0= bestK0; // update for next fit
    startXiD= bestXiD; // update for next fit

    cout <<"iTime: "<<tt<<"   K0= " << bestK0 << " +- "<<errK0<<" | "<<"    N0: "<<bestN0<<"     Chi2: "<<bestChi2_SOHO[tt]<< " / "<<NDF_SOHO[tt]<<endl;
    
    // ---- put results into graph ----
    grK0vsTime_SOHO->SetPoint(tt, crdata->xTimeSOHO_ProtonFlux[tt], bestK0);
    grK0vsTime_SOHO->SetPointError(tt, crdata->eTimeSOHO_ProtonFlux, errK0);
    grN0vsTime_SOHO->SetPoint(tt, crdata->xTimeSOHO_ProtonFlux[tt], bestN0);
    grN0vsTime_SOHO->SetPointError(tt, crdata->eTimeSOHO_ProtonFlux, errN0);
    grXiDvsTime_SOHO->SetPoint(tt, crdata->xTimeSOHO_ProtonFlux[tt], bestXiD);
    grXiDvsTime_SOHO->SetPointError(tt, crdata->eTimeSOHO_ProtonFlux, errXiD);
}



  // ---- fitting loop: BESS ---
  cout<<"**** Fit to BESS ****"<<endl;

  // Dedicated NORM for BESS   startN0= 0.96;
  startN0=   startN0_BESSPolar;

  
  for(int tt=0;tt<nTimeBESS_Proton;tt++){
    gMinuit->DefineParameter(0, "K_{0}", startK0, stepK0, minK0, maxK0);
    gMinuit->DefineParameter(1, "Experiment", iExp_BESS, 0, 0, 0);
    gMinuit->DefineParameter(2, "Time Index", tt, 0, 0, 0);
    gMinuit->DefineParameter(3, "XiD", startXiD, stepXiD, minXiD, maxXiD);
    gMinuit->DefineParameter(4, "NORM", startN0, stepN0, minN0, maxN0);
    gMinuit->FixParameter(1); // FIX TO WHAT? 
    gMinuit->FixParameter(2); 
    if(!kDriftON) gMinuit->FixParameter(3); // set Drift OFF
    if(!kReNorm) gMinuit->FixParameter(4); // No Renormalization parameter
    
    arglist[0] = NITERATIONS;
    arglist[1] = 1.;
    gMinuit->mnexcm("MIGRAD", arglist ,2,ierflg);
    
    double bestK0= 0;
    double errK0 = 0.;
    double bestN0= 0;
    double errN0 = 0.;
    double bestXiD= 0;
    double errXiD = 0.;
    gMinuit->GetParameter(0,bestK0,errK0);
    gMinuit->GetParameter(3,bestXiD,errXiD);
    gMinuit->GetParameter(4,bestN0,errN0);
    startK0= bestK0; // update for next fit
    startXiD= bestXiD; // update for next fit

    cout <<"iTime: "<<tt<<"   K0= " << bestK0 << " +- "<<errK0<<" | "<<"    N0: "<<bestN0<<"     Chi2: "<<bestChi2_BESS[tt]<< " / "<<NDF_BESS[tt]<<endl;
    
    // ---- put results into graph ----
    grK0vsTime_BESS->SetPoint(tt, crdata->xTimeBESS_ProtonFlux[tt], bestK0);
    grK0vsTime_BESS->SetPointError(tt, crdata->eTimeBESS_ProtonFlux[tt], errK0);
    grN0vsTime_BESS->SetPoint(tt, crdata->xTimeBESS_ProtonFlux[tt], bestN0);
    grN0vsTime_BESS->SetPointError(tt, crdata->eTimeBESS_ProtonFlux[tt], errN0);
    grXiDvsTime_BESS->SetPoint(tt, crdata->xTimeBESS_ProtonFlux[tt], bestXiD);
    grXiDvsTime_BESS->SetPointError(tt, crdata->eTimeBESS_ProtonFlux[tt], errXiD);
  }




  // ---- fitting loop: BESSTeV ---
  cout<<"**** Fit to BESSTeV ****"<<endl;

  // Dedicated NORM for BESS startN0= 1.0;
  startN0=   startN0_BESSTeV;

  for(int tt=0;tt<nTimeBESSTeV_Proton;tt++){
    gMinuit->DefineParameter(0, "K_{0}", startK0, stepK0, minK0, maxK0);
    gMinuit->DefineParameter(1, "Experiment", iExp_BESSTeV, 0, 0, 0);
    gMinuit->DefineParameter(2, "Time Index", tt, 0, 0, 0);
    gMinuit->DefineParameter(3, "XiD", startXiD, stepXiD, minXiD, maxXiD);
    gMinuit->DefineParameter(4, "NORM", startN0, stepN0, minN0, maxN0);
    gMinuit->FixParameter(1); // FIX TO WHAT? 
    gMinuit->FixParameter(2); 
    if(!kDriftON) gMinuit->FixParameter(3); // set Drift OFF
    if(!kReNorm) gMinuit->FixParameter(4); // No Renormalization parameter
    
    arglist[0] = NITERATIONS;
    arglist[1] = 1.;
    gMinuit->mnexcm("MIGRAD", arglist ,2,ierflg);
    
    double bestK0= 0;
    double errK0 = 0.;
    double bestN0= 0;
    double errN0 = 0.;
    double bestXiD= 0;
    double errXiD = 0.;
    gMinuit->GetParameter(0,bestK0,errK0);
    gMinuit->GetParameter(3,bestXiD,errXiD);
    gMinuit->GetParameter(4,bestN0,errN0);
    startK0= bestK0; // update for next fit
    startXiD= bestXiD; // update for next fit

    cout <<"iTime: "<<tt<<"   K0= " << bestK0 << " +- "<<errK0<<" | "<<"    N0: "<<bestN0<<"     Chi2: "<<bestChi2_BESSTeV[tt]<< " / "<<NDF_BESSTeV[tt]<<endl;
    
    // ---- put results into graph ----
    grK0vsTime_BESSTeV->SetPoint(tt, crdata->xTimeBESSTeV_ProtonFlux[tt], bestK0);
    grK0vsTime_BESSTeV->SetPointError(tt, crdata->eTimeBESSTeV_ProtonFlux[tt], errK0);
    grN0vsTime_BESSTeV->SetPoint(tt, crdata->xTimeBESSTeV_ProtonFlux[tt], bestN0);
    grN0vsTime_BESSTeV->SetPointError(tt, crdata->eTimeBESSTeV_ProtonFlux[tt], errN0);
    grXiDvsTime_BESSTeV->SetPoint(tt, crdata->xTimeBESSTeV_ProtonFlux[tt], bestXiD);
    grXiDvsTime_BESSTeV->SetPointError(tt, crdata->eTimeBESSTeV_ProtonFlux[tt], errXiD);
  }


  
  // ---- fitting loop: BESS00 ---
  cout<<"**** Fit to BESS00 ****"<<endl;

  // Dedicated NORM for BESS00 startN0= 1.0;
  startN0= startN0_BESS00;
  
  for(int tt=0;tt<nTimeBESS00_Proton;tt++){
    gMinuit->DefineParameter(0, "K_{0}", startK0, stepK0, minK0, maxK0);
    gMinuit->DefineParameter(1, "Experiment", iExp_BESS00, 0, 0, 0);
    gMinuit->DefineParameter(2, "Time Index", tt, 0, 0, 0);
    gMinuit->DefineParameter(3, "XiD", startXiD, stepXiD, minXiD, maxXiD);
    gMinuit->DefineParameter(4, "NORM", startN0, stepN0, minN0, maxN0);
    gMinuit->FixParameter(1); // FIX TO WHAT? 
    gMinuit->FixParameter(2); 
    if(!kDriftON) gMinuit->FixParameter(3); // set Drift OFF
    if(!kReNorm) gMinuit->FixParameter(4); // No Renormalization parameter
    
    arglist[0] = NITERATIONS;
    arglist[1] = 1.;
    gMinuit->mnexcm("MIGRAD", arglist ,2,ierflg);
    
    double bestK0= 0;
    double errK0 = 0.;
    double bestN0= 0;
    double errN0 = 0.;
    double bestXiD= 0;
    double errXiD = 0.;
    gMinuit->GetParameter(0,bestK0,errK0);
    gMinuit->GetParameter(3,bestXiD,errXiD);
    gMinuit->GetParameter(4,bestN0,errN0);
    startK0= bestK0; // update for next fit
    startXiD= bestXiD; // update for next fit

    cout <<"iTime: "<<tt<<"   K0= " << bestK0 << " +- "<<errK0<<" | "<<"    N0: "<<bestN0<<"     Chi2: "<<bestChi2_BESS00[tt]<< " / "<<NDF_BESS00[tt]<<endl;
    
    // ---- put results into graph ----
    grK0vsTime_BESS00->SetPoint(tt, crdata->xTimeBESS00_ProtonFlux[tt], bestK0);
    grK0vsTime_BESS00->SetPointError(tt, crdata->eTimeBESS00_ProtonFlux[tt], errK0);
    grN0vsTime_BESS00->SetPoint(tt, crdata->xTimeBESS00_ProtonFlux[tt], bestN0);
    grN0vsTime_BESS00->SetPointError(tt, crdata->eTimeBESS00_ProtonFlux[tt], errN0);
    grXiDvsTime_BESS00->SetPoint(tt, crdata->xTimeBESS00_ProtonFlux[tt], bestXiD);
    grXiDvsTime_BESS00->SetPointError(tt, crdata->eTimeBESS00_ProtonFlux[tt], errXiD);
}


  // ---- fitting loop: AMS02 ---
  cout<<"**** Fit to AMS02 ****"<<endl;

  // Dedicated NORM for AMS02 startN0= 1.0;
  startN0 =  startN0_AMS02;

  int indTime=0;
  for(int tt=0;tt<nTimeAMS02_Proton;tt++){
    // int indTime = tt; // Setta TUTTI i punti (anche vuoti; commentare per eliminarli)
    if(tt==amsdata->SkipThis || tt==amsdata->SkipThis+1) continue;
    
    gMinuit->DefineParameter(0, "K_{0}", startK0, stepK0, minK0, maxK0);
    gMinuit->DefineParameter(1, "Experiment", iExp_AMS02, 0, 0, 0);
    gMinuit->DefineParameter(2, "Time Index", tt, 0, 0, 0);
    gMinuit->DefineParameter(3, "XiD", startXiD, stepXiD, minXiD, maxXiD);
    gMinuit->DefineParameter(4, "NORM", startN0, stepN0, minN0, maxN0);
    gMinuit->FixParameter(1); // FIX TO WHAT? 
    gMinuit->FixParameter(2); 
    if(!kDriftON) gMinuit->FixParameter(3); // set Drift OFF
    if(!kReNorm) gMinuit->FixParameter(4); // No Renormalization parameter
    
    arglist[0] = NITERATIONS;
    arglist[1] = 1.;
    gMinuit->mnexcm("MIGRAD", arglist ,2,ierflg);
    
    double bestK0= 0;
    double errK0 = 0.;
    double bestN0= 0;
    double errN0 = 0.;
    double bestXiD= 0;
    double errXiD = 0.;
    gMinuit->GetParameter(0,bestK0,errK0);
    gMinuit->GetParameter(3,bestXiD,errXiD);
    gMinuit->GetParameter(4,bestN0,errN0);
    startK0= bestK0; // update for next fit
    startXiD= bestXiD; // update for next fit

    cout <<"iTime: "<<tt<<"   K0= " << bestK0 << " +- "<<errK0<<" | "<<"    N0: "<<bestN0<<"     Chi2: "<<bestChi2_AMS02[tt]<< " / "<<NDF_AMS02[tt]<<endl;
    
    // ---- put results into graph ----
    grK0vsTime_AMS02->SetPoint(indTime, amsdata->xTimeAMS02_ProtonFlux[tt], bestK0);
    grK0vsTime_AMS02->SetPointError(indTime, amsdata->eTimeAMS02_ProtonFlux, errK0);
    grN0vsTime_AMS02->SetPoint(indTime, amsdata->xTimeAMS02_ProtonFlux[tt], bestN0);
    grN0vsTime_AMS02->SetPointError(indTime, amsdata->eTimeAMS02_ProtonFlux, errN0);
    grXiDvsTime_AMS02->SetPoint(indTime, amsdata->xTimeAMS02_ProtonFlux[tt], bestXiD);
    grXiDvsTime_AMS02->SetPointError(indTime, amsdata->eTimeAMS02_ProtonFlux, errXiD);

    indTime++; // increment time index
  }

  // cout<<"TIME-POINTS OF AMS02: "<<grK0vsTime_AMS02->GetN()<<endl;


  
  // ----- plot results ----

  // ---- draw results ---
  TDatime Date1( 2000, 01, 01, 0, 0, 0);
  TDatime Date2( 2018, 01, 01, 0, 0, 0);
  double UTMIN= (double)Date1.Convert();
  double UTMAX= (double)Date2.Convert();

  double FMIN = 0.;
  double FMAX = 10.;
  TH2F* hFrameK0vsTime= new TH2F("hFrameK0vsTime","K0 vs Time",200, UTMIN, UTMAX, 500, FMIN, FMAX);
  SetStyleHistoVSTime(hFrameK0vsTime);
  hFrameK0vsTime->GetYaxis()->SetNdivisions(506);
  hFrameK0vsTime->GetYaxis()->SetTitleSize(0.07);
  hFrameK0vsTime->GetYaxis()->SetTitle("k_{0} ");
  hFrameK0vsTime->GetYaxis()->SetTitleOffset(0.50);

  
  FMIN= -5;
  FMAX=  5;
  TH2F* hFrameXiDvsTime= new TH2F("hFrameXiDvsTime","XiD vs Time",200, UTMIN, UTMAX, 500, FMIN, FMAX);
  SetStyleHistoVSTime(hFrameXiDvsTime);
  hFrameXiDvsTime->GetYaxis()->SetNdivisions(506);
  hFrameXiDvsTime->GetYaxis()->SetTitleSize(0.07);
  hFrameXiDvsTime->GetYaxis()->SetTitle("#xi_{d} ");
  hFrameXiDvsTime->GetYaxis()->SetTitleOffset(0.50);

  FMIN=  0.80;
  FMAX=  1.20;
  TH2F* hFrameN0vsTime= new TH2F("hFrameN0vsTime","N0 vs Time",200, UTMIN, UTMAX, 500, FMIN, FMAX);
  SetStyleHistoVSTime(hFrameN0vsTime);
  hFrameN0vsTime->GetYaxis()->SetNdivisions(506);
  hFrameN0vsTime->GetYaxis()->SetTitleSize(0.07);
  hFrameN0vsTime->GetYaxis()->SetTitle("N_{0} ");
  hFrameN0vsTime->GetYaxis()->SetTitleOffset(0.50);


  // ---- make results & predictions  ----  

  // --- time profiles on p-He and p/He ratios at given Ekn/Rig ---
  

  // ---- pick up AMS-02 data ---- [from Claudio]
  TGraphErrors* grDataTimeProfileProton_AMS02_Rig[45];
  TGraphErrors* grDataTimeProfileHelium_AMS02_Rig[45];
  TGraphErrors* grDataTimeProfilePHeRatio_AMS02_Rig[45];

  TFile* inAMSDataRig= new TFile("$NTBASEDIR/THEORY/00matisse/ANALYSIS/ExperimentalData/grAMS02_pHe_monthly_BR2426_BR2506_OFFICIAL_Rigidity_HI.root","READ");
  inAMSDataRig->cd();

  for(int rr=0;rr<45;rr++){
    grDataTimeProfileProton_AMS02_Rig[rr]= (TGraphErrors*)inAMSDataRig->Get(Form("g_pflux_bin%02d",rr));
    grDataTimeProfileProton_AMS02_Rig[rr]->SetMarkerSize(0.5);
  }

for(int rr=0;rr<40;rr++){
      grDataTimeProfileHelium_AMS02_Rig[rr]= (TGraphErrors*)inAMSDataRig->Get(Form("g_heflux_bin%02d",rr));
      grDataTimeProfileHelium_AMS02_Rig[rr]->SetMarkerSize(0.5);
 }  

 for(int rr=0;rr<40;rr++){
   grDataTimeProfilePHeRatio_AMS02_Rig[rr]= (TGraphErrors*)inAMSDataRig->Get(Form("g_pohe_bin%02d",rr));
   grDataTimeProfilePHeRatio_AMS02_Rig[rr]->SetMarkerSize(0.5);
 } 


 
 // grDataTimeProfileProton_AMS02_Rig[0]->Draw("apZ");
 // grDataTimeProfilePHeRatio_AMS02_Rig[0]->Draw("apZ");

 // --- pick up AMS-02 RIGIDITY grids for P and He ---
 TH1F* hAMS02_RigidityGrid_Proton= (TH1F*)inAMSDataRig->Get("h_pflux_BR2426"); 
 hAMS02_RigidityGrid_Proton->SetDirectory(0);
 
 TH1F* hAMS02_RigidityGrid_Helium= (TH1F*)inAMSDataRig->Get("h_heflux_BR2426"); 
 hAMS02_RigidityGrid_Helium->SetDirectory(0);
 
 inAMSDataRig->Close();

 double pRigAMS02Proton[45];
 double pRigAMS02Helium[40]; 

 for(int pp=0;pp<45;pp++) pRigAMS02Proton[pp]=hAMS02_RigidityGrid_Proton->GetXaxis()->GetBinCenterLog(pp+1);
for(int pp=0;pp<40;pp++) pRigAMS02Helium[pp]=hAMS02_RigidityGrid_Helium->GetXaxis()->GetBinCenterLog(pp+1);
 
 
  // ---- COMPUTE CALCULATIONS FOR AMS-02 ----
  TGraph* grResultsTimeProfileProton_AMS02_Rig[45];
  TGraph* grResultsTimeProfileHelium3_AMS02_Rig[40];  
  TGraph* grResultsTimeProfileHelium4_AMS02_Rig[40];  
  TGraph* grResultsTimeProfilePHeRatio3_AMS02_Rig[40];  
  TGraph* grResultsTimeProfilePHeRatio4_AMS02_Rig[40];  



  // ---- compute proton flux over AMS-02 points ----
  int indRig  = 0;
  for(int rr=0;rr<45;rr++){
    indRig = rr;
    double RigGeV= pRigAMS02Proton[rr]; 
    double LnRigMeV= log( (RigGeV*1.e+3) ); // get ln-rigidity in MV 

    grResultsTimeProfileProton_AMS02_Rig[indRig]= new TGraph();
    grResultsTimeProfileProton_AMS02_Rig[indRig]->SetName(Form("grTimeProfile_Proton_R%2.2f",RigGeV));
    grResultsTimeProfileProton_AMS02_Rig[indRig]->SetLineColor(kRed+1);
    grResultsTimeProfileProton_AMS02_Rig[indRig]->SetLineWidth(3);
    
    int indTime = 0;
    for(int tt=0;tt<nTimeAMS02_Proton;tt++){
      if(tt==amsdata->SkipThis || tt==amsdata->SkipThis+1) continue;
    
      // get K0 from proton fit
      double Time  = grK0vsTime_AMS02->GetX()[indTime];
      double LnK0  = log( grK0vsTime_AMS02->GetY()[indTime] );
      double XiD   = ( grXiDvsTime_AMS02->GetY()[indTime] );
      // double LnFluxMeV= hLn_Flux_vs_XiD_vs_KScale_vs_Rigidity_Z1_A1->Interpolate( LnRigMeV, LnK0, XiD );
      double LnFluxMeV;
      if(kDriftON) LnFluxMeV= Interpolation( hLn_Flux_vs_XiD_vs_KScale_vs_Rigidity_Z1_A1, LnRigMeV, LnK0, XiD );
      if(!kDriftON)LnFluxMeV = Interpolation( hLn_Flux_vs_KScale_vs_Rigidity_Z1_A1, LnRigMeV, LnK0); //NTJAN2018 arrivati qui
      double FluxGeV = exp(LnFluxMeV)*1.e+3;
      if(kReNorm) FluxGeV*= grN0vsTime_AMS02->GetY()[indTime]; // re-normalize by N0
      
      //cout<<"Time: "<<tt<<" / "<<Time<<"   Rig: "<<RigGeV <<"    lnK0: "<<LnK0<<"   Flux: "<<FluxGeV<<endl;
      grResultsTimeProfileProton_AMS02_Rig[rr]->SetPoint(indTime, Time, FluxGeV);
      indTime++;
    } 
  }


  // ---- compute helium 3/4 fluxes over AMS-02 points ----
  indRig  = 0;
  for(int rr=0;rr<40;rr++){
    indRig = rr;
    double RigGeV= pRigAMS02Helium[rr]; 
    double LnRigMeV= log( (RigGeV*1.e+3) ); // get ln-rigidity in MV 
    
    grResultsTimeProfileHelium3_AMS02_Rig[indRig]= new TGraph();
    grResultsTimeProfileHelium3_AMS02_Rig[indRig]->SetName(Form("grTimeProfile_Helium4_R%2.2f",RigGeV));
    grResultsTimeProfileHelium3_AMS02_Rig[indRig]->SetLineColor(kOrange+1);
    grResultsTimeProfileHelium3_AMS02_Rig[indRig]->SetLineWidth(3);

    grResultsTimeProfileHelium4_AMS02_Rig[indRig]= new TGraph();
    grResultsTimeProfileHelium4_AMS02_Rig[indRig]->SetName(Form("grTimeProfile_Helium4_R%2.2f",RigGeV));
    grResultsTimeProfileHelium4_AMS02_Rig[indRig]->SetLineColor(kGreen+1);
    grResultsTimeProfileHelium4_AMS02_Rig[indRig]->SetLineWidth(3);

    grResultsTimeProfilePHeRatio3_AMS02_Rig[indRig]= new TGraph();
    grResultsTimeProfilePHeRatio3_AMS02_Rig[indRig]->SetName(Form("grTimeProfile_PHeRatio3_R%2.2f",RigGeV));
    grResultsTimeProfilePHeRatio3_AMS02_Rig[indRig]->SetLineColor(kOrange+1);
    grResultsTimeProfilePHeRatio3_AMS02_Rig[indRig]->SetLineWidth(3);
    
    grResultsTimeProfilePHeRatio4_AMS02_Rig[indRig]= new TGraph();
    grResultsTimeProfilePHeRatio4_AMS02_Rig[indRig]->SetName(Form("grTimeProfile_PHeRatio4_R%2.2f",RigGeV));
    grResultsTimeProfilePHeRatio4_AMS02_Rig[indRig]->SetLineColor(kGreen+1);
    grResultsTimeProfilePHeRatio4_AMS02_Rig[indRig]->SetLineWidth(3);



    int indTime = 0;
    for(int tt=0;tt<nTimeAMS02_Helium;tt++){
      if(tt==amsdata->SkipThis || tt==amsdata->SkipThis+1) continue;
    
      // get K0 from *proton* fit
      double Time  = grK0vsTime_AMS02->GetX()[indTime];
      double LnK0  = log( grK0vsTime_AMS02->GetY()[indTime] );
      double XiD   = grXiDvsTime_AMS02->GetY()[indTime];

      // double LnFluxMeV_He3= hLn_Flux_vs_XiD_vs_KScale_vs_Rigidity_Z2_A3->Interpolate( LnRigMeV, LnK0, XiD );
      // double LnFluxMeV_He4= hLn_Flux_vs_XiD_vs_KScale_vs_Rigidity_Z2_A4->Interpolate( LnRigMeV, LnK0, XiD );

      double LnFluxMeV_He3=0.;
      double LnFluxMeV_He4=0.;

      if(kDriftON){
	LnFluxMeV_He3= Interpolation(hLn_Flux_vs_XiD_vs_KScale_vs_Rigidity_Z2_A3, LnRigMeV, LnK0, XiD );
	LnFluxMeV_He4= Interpolation(hLn_Flux_vs_XiD_vs_KScale_vs_Rigidity_Z2_A4, LnRigMeV, LnK0, XiD );
      }

      if(!kDriftON){
	LnFluxMeV_He3= Interpolation(hLn_Flux_vs_KScale_vs_Rigidity_Z2_A3, LnRigMeV, LnK0 );
	LnFluxMeV_He4= Interpolation(hLn_Flux_vs_KScale_vs_Rigidity_Z2_A4, LnRigMeV, LnK0 );
      }
      
      
      double FluxGeV_He3 = exp(LnFluxMeV_He3)*1.e+3;
      double FluxGeV_He4 = exp(LnFluxMeV_He4)*1.e+3;

      if(kReNorm){
	FluxGeV_He3*= grN0vsTime_AMS02->GetY()[indTime]; // re-normalize by N0
	FluxGeV_He4*= grN0vsTime_AMS02->GetY()[indTime]; // re-normalize by N0
      }
      
      // cout<<FluxGeV_He3<<"   --  "<<FluxGeV_He4<<endl;

      // compute He fluxes 
      grResultsTimeProfileHelium3_AMS02_Rig[rr]->SetPoint(indTime, Time, FluxGeV_He3 + FluxGeV_He4);
      grResultsTimeProfileHelium4_AMS02_Rig[rr]->SetPoint(indTime, Time, FluxGeV_He4);
      
      // compute p/He ratios
      double FluxGeV_Proton= grResultsTimeProfileProton_AMS02_Rig[rr+5]->GetY()[indTime];
      double RatioPHe3= FluxGeV_Proton/(FluxGeV_He3 + FluxGeV_He4);
      double RatioPHe4= FluxGeV_Proton/FluxGeV_He4;
      grResultsTimeProfilePHeRatio3_AMS02_Rig[rr]->SetPoint(indTime, Time, RatioPHe3);
      grResultsTimeProfilePHeRatio4_AMS02_Rig[rr]->SetPoint(indTime, Time, RatioPHe4);
      
      indTime++;
    }
  }



  
  // **** results on energy/rigidity specta : servono a qualcosa? ****
  
  // ---- results: flux vs ekn ---
  TGraph* grResultsFluxVSEkn_PAMELA[nTimePAMELA_Proton];
  TGraph* grResultsFluxVSEkn_SOHO[nTimeSOHO_Proton];
  TGraph* grResultsFluxVSEkn_BESS[nTimeBESS_Proton];
  TGraph* grResultsFluxVSEkn_BESSTeV[nTimeBESSTeV_Proton];
  TGraph* grResultsFluxVSEkn_BESS00[nTimeBESS00_Proton];
  TGraph* grResultsFluxVSEkn_AMS02[nTimeAMS02_Proton];
  TGraph* grResultsFluxVSEkn_He4_AMS02[nTimeAMS02_Proton];
  TGraph* grResultsFluxVSEkn_He3_AMS02[nTimeAMS02_Proton];
  
  // ---- results: flux vs rig ---
  TGraph* grResultsFluxVSRig_PAMELA[nTimePAMELA_Proton];
  TGraph* grResultsFluxVSRig_SOHO[nTimeSOHO_Proton];
  TGraph* grResultsFluxVSRig_BESS[nTimeBESS_Proton];
  TGraph* grResultsFluxVSRig_BESSTeV[nTimeBESSTeV_Proton];
  TGraph* grResultsFluxVSRig_BESS00[nTimeBESS00_Proton];
  TGraph* grResultsFluxVSRig_AMS02[nTimeAMS02_Proton];

  
  // --- N energy points for model ---
  int NEKN=0;
  if(kDriftON) NEKN= hLn_Flux_vs_XiD_vs_KScale_vs_KEnergy_Z1_A1->GetNbinsX();
  if(!kDriftON) NEKN= hLn_Flux_vs_KScale_vs_KEnergy_Z1_A1->GetNbinsX();

  
  // --- PAMELA ---
  for(int tt=0;tt<nTimePAMELA_Proton;tt++){
    grResultsFluxVSEkn_PAMELA[tt]= new TGraphErrors();
    grResultsFluxVSEkn_PAMELA[tt]->SetName(Form("grResultsFluxVSEkn_PAMELA_T%d",tt));
    grResultsFluxVSEkn_PAMELA[tt]->SetLineColor(kRed+1);
    grResultsFluxVSEkn_PAMELA[tt]->SetLineWidth(3);
    double LnK0  = log( grK0vsTime_PAMELA->GetY()[tt] );    
    double XiD   =  grXiDvsTime_PAMELA->GetY()[tt];
    
    for(int ee=0;ee<NEKN-1;ee++){
      double LnE=0;
      double LnFlux=0.;
      if(kDriftON)  LnE=hLn_Flux_vs_XiD_vs_KScale_vs_KEnergy_Z1_A1->GetXaxis()->GetBinCenter(ee+1);
      if(!kDriftON) LnE=hLn_Flux_vs_KScale_vs_KEnergy_Z1_A1->GetXaxis()->GetBinCenter(ee+1);      
      if(kDriftON)  LnFlux= Interpolation(hLn_Flux_vs_XiD_vs_KScale_vs_KEnergy_Z1_A1, LnE, LnK0, XiD );
      if(!kDriftON) LnFlux= Interpolation(hLn_Flux_vs_KScale_vs_KEnergy_Z1_A1, LnE, LnK0 );
      double EMeV= exp( LnE );
      double EGeV= EMeV*1.e-3;
      double FluxMeV= exp( LnFlux );
      double FluxGeV= FluxMeV*1.e+3;
      if(kReNorm) FluxGeV*= grN0vsTime_PAMELA->GetY()[tt]; // renormalize by N0
      grResultsFluxVSEkn_PAMELA[tt]->SetPoint(ee, EGeV, FluxGeV);
    }
  }

  
  // --- SOHO ---
  for(int tt=0;tt<nTimeSOHO_Proton;tt++){
    grResultsFluxVSEkn_SOHO[tt]= new TGraphErrors();
    grResultsFluxVSEkn_SOHO[tt]->SetName(Form("grResultsFluxVSEkn_SOHO_T%d",tt));
    grResultsFluxVSEkn_SOHO[tt]->SetLineColor(kRed+1);
    grResultsFluxVSEkn_SOHO[tt]->SetLineWidth(3);
    double LnK0  = log( grK0vsTime_SOHO->GetY()[tt] );    
    double XiD   = grXiDvsTime_SOHO->GetY()[tt];

    for(int ee=0;ee<NEKN-1;ee++){
      double LnE=0.;
      double LnFlux=0.;
      if(kDriftON)  LnE=hLn_Flux_vs_XiD_vs_KScale_vs_KEnergy_Z1_A1->GetXaxis()->GetBinCenter(ee+1);
      if(!kDriftON) LnE=hLn_Flux_vs_KScale_vs_KEnergy_Z1_A1->GetXaxis()->GetBinCenter(ee+1);      
      if(kDriftON)  LnFlux= Interpolation( hLn_Flux_vs_XiD_vs_KScale_vs_KEnergy_Z1_A1, LnE, LnK0, XiD );
      if(!kDriftON) LnFlux= Interpolation( hLn_Flux_vs_KScale_vs_KEnergy_Z1_A1, LnE, LnK0 );
      double EMeV= exp( LnE );
      double EGeV= EMeV*1.e-3;
      double FluxMeV= exp( LnFlux );
      double FluxGeV= FluxMeV*1.e+3;
      if(kReNorm) FluxGeV*= grN0vsTime_SOHO->GetY()[tt]; // renormalize by N0
      grResultsFluxVSEkn_SOHO[tt]->SetPoint(ee, EGeV, FluxGeV);
    }
  }


  // --- BESS ---
  for(int tt=0;tt<nTimeBESS_Proton;tt++){
    grResultsFluxVSEkn_BESS[tt]= new TGraphErrors();
    grResultsFluxVSEkn_BESS[tt]->SetName(Form("grResultsFluxVSEkn_BESS_T%d",tt));
    grResultsFluxVSEkn_BESS[tt]->SetLineColor(kRed+1);
    grResultsFluxVSEkn_BESS[tt]->SetLineWidth(3);
    double LnK0  = log( grK0vsTime_BESS->GetY()[tt] );    
    double XiD   = grXiDvsTime_BESS->GetY()[tt];
    
    for(int ee=0;ee<NEKN-1;ee++){
      double LnE=0.;
      double LnFlux=0.;
      if(kDriftON)  LnE=hLn_Flux_vs_XiD_vs_KScale_vs_KEnergy_Z1_A1->GetXaxis()->GetBinCenter(ee+1);
      if(!kDriftON) LnE=hLn_Flux_vs_KScale_vs_KEnergy_Z1_A1->GetXaxis()->GetBinCenter(ee+1);
      if(kDriftON)  LnFlux= Interpolation(hLn_Flux_vs_XiD_vs_KScale_vs_KEnergy_Z1_A1, LnE, LnK0, XiD);
      if(!kDriftON) LnFlux= Interpolation(hLn_Flux_vs_KScale_vs_KEnergy_Z1_A1, LnE, LnK0);
       
      double EMeV= exp( LnE );
      double EGeV= EMeV*1.e-3;
      double FluxMeV= exp( LnFlux );
      double FluxGeV= FluxMeV*1.e+3;
      if(kReNorm) FluxGeV*= grN0vsTime_BESS->GetY()[tt]; // renormalize by N0
      grResultsFluxVSEkn_BESS[tt]->SetPoint(ee, EGeV, FluxGeV);
    }
  }


  // --- BESSTeV ---
  for(int tt=0;tt<nTimeBESSTeV_Proton;tt++){
    grResultsFluxVSEkn_BESSTeV[tt]= new TGraphErrors();
    grResultsFluxVSEkn_BESSTeV[tt]->SetName(Form("grResultsFluxVSEkn_BESSTeV_T%d",tt));
    grResultsFluxVSEkn_BESSTeV[tt]->SetLineColor(kRed+1);
    grResultsFluxVSEkn_BESSTeV[tt]->SetLineWidth(3);
    double LnK0  = log( grK0vsTime_BESSTeV->GetY()[tt] );    
    double XiD   = grXiDvsTime_BESSTeV->GetY()[tt];

    for(int ee=0;ee<NEKN-1;ee++){
      double LnE=0.;
      double LnFlux=0.;
      if(kDriftON)  LnE=hLn_Flux_vs_XiD_vs_KScale_vs_KEnergy_Z1_A1->GetXaxis()->GetBinCenter(ee+1);
      if(!kDriftON) LnE=hLn_Flux_vs_KScale_vs_KEnergy_Z1_A1->GetXaxis()->GetBinCenter(ee+1);
      if(kDriftON)  LnFlux= Interpolation( hLn_Flux_vs_XiD_vs_KScale_vs_KEnergy_Z1_A1, LnE, LnK0, XiD );
      if(!kDriftON) LnFlux= Interpolation( hLn_Flux_vs_KScale_vs_KEnergy_Z1_A1, LnE, LnK0 );
      double EMeV= exp( LnE );
      double EGeV= EMeV*1.e-3;
      double FluxMeV= exp( LnFlux );
      double FluxGeV= FluxMeV*1.e+3;
      if(kReNorm) FluxGeV*= grN0vsTime_BESSTeV->GetY()[tt]; // renormalize by N0
      grResultsFluxVSEkn_BESSTeV[tt]->SetPoint(ee, EGeV, FluxGeV);
    }
  }

  
  // --- BESS00 ---
  for(int tt=0;tt<nTimeBESS00_Proton;tt++){
    grResultsFluxVSEkn_BESS00[tt]= new TGraphErrors();
    grResultsFluxVSEkn_BESS00[tt]->SetName(Form("grResultsFluxVSEkn_BESS00_T%d",tt));
    grResultsFluxVSEkn_BESS00[tt]->SetLineColor(kRed+1);
    grResultsFluxVSEkn_BESS00[tt]->SetLineWidth(3);
    double LnK0  = log( grK0vsTime_BESS00->GetY()[tt] );    
    double XiD   = grXiDvsTime_BESS00->GetY()[tt];

    for(int ee=0;ee<NEKN-1;ee++){
      double LnE=0.;
      double LnFlux=0.;
      if(kDriftON)  LnE=hLn_Flux_vs_XiD_vs_KScale_vs_KEnergy_Z1_A1->GetXaxis()->GetBinCenter(ee+1);
      if(!kDriftON) LnE=hLn_Flux_vs_KScale_vs_KEnergy_Z1_A1->GetXaxis()->GetBinCenter(ee+1);
      if(kDriftON)   LnFlux= Interpolation(hLn_Flux_vs_XiD_vs_KScale_vs_KEnergy_Z1_A1, LnE, LnK0, XiD );
      if(!kDriftON)  LnFlux= Interpolation(hLn_Flux_vs_KScale_vs_KEnergy_Z1_A1, LnE, LnK0 );
      double EMeV= exp( LnE );
      double EGeV= EMeV*1.e-3;
      double FluxMeV= exp( LnFlux );
      double FluxGeV= FluxMeV*1.e+3;
      if(kReNorm) FluxGeV*= grN0vsTime_BESS00->GetY()[tt]; // renormalize by N0
      grResultsFluxVSEkn_BESS00[tt]->SetPoint(ee, EGeV, FluxGeV);
    }
  }
 

  // --- AMS02 ---
  indTime=0; // per gestire i due buchi
  for(int tt=0;tt<nTimeAMS02_Proton;tt++){
    grResultsFluxVSEkn_AMS02[tt]= new TGraphErrors();
    grResultsFluxVSEkn_AMS02[tt]->SetName(Form("grResultsFluxVSEkn_AMS02_T%d",tt));
    grResultsFluxVSEkn_AMS02[tt]->SetLineColor(kRed+1);
    grResultsFluxVSEkn_AMS02[tt]->SetLineWidth(3);

    grResultsFluxVSEkn_He4_AMS02[tt]= new TGraphErrors();
    grResultsFluxVSEkn_He4_AMS02[tt]->SetName(Form("grResultsFluxVSEkn_He4_AMS02_T%d",tt));
    grResultsFluxVSEkn_He4_AMS02[tt]->SetLineColor(kGreen+1);
    grResultsFluxVSEkn_He4_AMS02[tt]->SetLineWidth(3);

    // He3 [actually HeTOT?]
    grResultsFluxVSEkn_He3_AMS02[tt]= new TGraphErrors();
    grResultsFluxVSEkn_He3_AMS02[tt]->SetName(Form("grResultsFluxVSEkn_He3_AMS02_T%d",tt));
    grResultsFluxVSEkn_He3_AMS02[tt]->SetLineColor(kOrange+1);
    grResultsFluxVSEkn_He3_AMS02[tt]->SetLineWidth(3);

    
    if(tt==amsdata->SkipThis || tt==amsdata->SkipThis+1) continue;
    double LnK0  = log( grK0vsTime_AMS02->GetY()[indTime] );    
    double XiD   = grXiDvsTime_AMS02->GetY()[indTime];

    // NB: energy grid from proton | but for helium it is different
    for(int ee=0;ee<NEKN-1;ee++){

      // energy
      double LnE=0.;
      if(kDriftON)  LnE=hLn_Flux_vs_XiD_vs_KScale_vs_KEnergy_Z1_A1->GetXaxis()->GetBinCenter(ee+1);
      if(!kDriftON) LnE=hLn_Flux_vs_KScale_vs_KEnergy_Z1_A1->GetXaxis()->GetBinCenter(ee+1);
      double EMeV= exp( LnE );
      double EGeV= EMeV*1.e-3;
      if(EGeV > 50) continue;
      
      // proton
      double LnFlux= 0.;
      if(kDriftON) LnFlux = Interpolation(hLn_Flux_vs_XiD_vs_KScale_vs_KEnergy_Z1_A1, LnE, LnK0, XiD );
      if(!kDriftON) LnFlux = Interpolation(hLn_Flux_vs_KScale_vs_KEnergy_Z1_A1, LnE, LnK0 );
      
      double FluxMeV= exp( LnFlux );
      double FluxGeV= FluxMeV*1.e+3;
      if(kReNorm) FluxGeV*= grN0vsTime_AMS02->GetY()[indTime]; // renormalize by N0
      grResultsFluxVSEkn_AMS02[tt]->SetPoint(ee, EGeV, FluxGeV);

      // helium
      double LnFluxHe4= 0.;
      if(kDriftON) LnFluxHe4= Interpolation(hLn_Flux_vs_XiD_vs_KScale_vs_KEnergy_Z2_A4, LnE, LnK0, XiD );
      if(!kDriftON) LnFluxHe4= Interpolation(hLn_Flux_vs_KScale_vs_KEnergy_Z2_A4, LnE, LnK0 );
      double FluxHe4GeV= exp(LnFluxHe4)*1.e+3;
      if(kReNorm) FluxHe4GeV*= grN0vsTime_AMS02->GetY()[indTime]; // renormalize also He4 by N0!!!
      grResultsFluxVSEkn_He4_AMS02[tt]->SetPoint(ee, EGeV, FluxHe4GeV);

      // helium3 [actually HeTOT: accounting for He3 NTJAN2018]
      double LnFluxHe3= 0.;
      if(kDriftON) LnFluxHe3= Interpolation(hLn_Flux_vs_XiD_vs_KScale_vs_KEnergy_Z2_A3, LnE, LnK0, XiD );
      if(!kDriftON) LnFluxHe3= Interpolation(hLn_Flux_vs_KScale_vs_KEnergy_Z2_A3, LnE, LnK0 );
      double FluxHe3GeV= exp(LnFluxHe3)*1.e+3;
      if(kReNorm) FluxHe3GeV*= grN0vsTime_AMS02->GetY()[indTime]; // renormalize also He3 by N0
      grResultsFluxVSEkn_He3_AMS02[tt]->SetPoint(ee, EGeV, FluxHe4GeV+FluxHe3GeV); // HeTOT FLUX!!! NTJAN2018

      
    }
    indTime++;
  }




  // ---- draw some results ----


  // ---- 1: Plot best-fit K0 and XiD vs Time ----
  TCanvas* ccBestFitK0VSTime= new TCanvas(Form("ccBestFitK0VSTime_N%d",0), Form("BEST-FIT PARAMETERS | N%d",0), 900, 730 );
  ccBestFitK0VSTime->Divide(1,3,0.01,0.01);  

  // --- K0 vs Time fit results ---
  ccBestFitK0VSTime->cd(1);  
  gPad->SetBottomMargin(0.15);
  gPad->SetTopMargin(0.05);
  gPad->SetLeftMargin(0.150);
  gPad->SetRightMargin(0.06);
  gPad->SetTicky(1);
  gPad->SetTickx(1); 
  hFrameK0vsTime->Draw();
  SetTimeGraphStyle(grK0vsTime_PAMELA);
  SetTimeGraphStyle(grK0vsTime_SOHO);
  SetTimeGraphStyle(grK0vsTime_BESS);
  SetTimeGraphStyle(grK0vsTime_BESSTeV);
  SetTimeGraphStyle(grK0vsTime_BESS00);
  SetTimeGraphStyle(grK0vsTime_AMS02);
  // SetTimeGraphStyle(grK0vsTime_AMSPRL2015);
  
  grK0vsTime_SOHO->Draw("pZ");
  grK0vsTime_PAMELA->Draw("pZ");
  grK0vsTime_BESS->Draw("pZ");
  grK0vsTime_BESSTeV->Draw("pZ");
  grK0vsTime_BESS00->Draw("pZ");
  // grK0vsTime_AMSPRL2015->Draw("pZ");
  grK0vsTime_AMS02->Draw("pZ");
  // grK0vsTime_AMS02Helium->Draw("pZl");
  gPad->SetGridy();


  // --- XiD vs Time fit results ---
  ccBestFitK0VSTime->cd(2);  
  gPad->SetBottomMargin(0.15);
  gPad->SetTopMargin(0.05);
  gPad->SetLeftMargin(0.150);
  gPad->SetRightMargin(0.06);
  gPad->SetTicky(1);
  gPad->SetTickx(1); 
  hFrameXiDvsTime->Draw();
  SetTimeGraphStyle(grXiDvsTime_PAMELA);
  SetTimeGraphStyle(grXiDvsTime_SOHO);
  SetTimeGraphStyle(grXiDvsTime_BESS);
  SetTimeGraphStyle(grXiDvsTime_BESSTeV);
  SetTimeGraphStyle(grXiDvsTime_BESS00);
  SetTimeGraphStyle(grXiDvsTime_AMS02);
  // SetTimeGraphStyle(grXiDvsTime_AMSPRL2015);
  
  grXiDvsTime_SOHO->Draw("pZ");
  grXiDvsTime_PAMELA->Draw("pZ");
  grXiDvsTime_BESS->Draw("pZ");
  grXiDvsTime_BESSTeV->Draw("pZ");
  grXiDvsTime_BESS00->Draw("pZ");
  // grXiDvsTime_AMSPRL2015->Draw("pZ");
  grXiDvsTime_AMS02->Draw("pZ");
  // grXiDvsTime_AMS02Helium->Draw("pZl");
  gPad->SetGridy();


  // --- N0 vs Time fit results ---
  ccBestFitK0VSTime->cd(3);  
  gPad->SetBottomMargin(0.15);
  gPad->SetTopMargin(0.05);
  gPad->SetLeftMargin(0.150);
  gPad->SetRightMargin(0.06);
  gPad->SetTicky(1);
  gPad->SetTickx(1); 
  hFrameN0vsTime->Draw();
  SetTimeGraphStyle(grN0vsTime_PAMELA);
  SetTimeGraphStyle(grN0vsTime_SOHO);
  SetTimeGraphStyle(grN0vsTime_BESS);
  SetTimeGraphStyle(grN0vsTime_BESSTeV);
  SetTimeGraphStyle(grN0vsTime_BESS00);
  SetTimeGraphStyle(grN0vsTime_AMS02);
  // SetTimeGraphStyle(grN0vsTime_AMSPRL2015);
  
  grN0vsTime_SOHO->Draw("pZ");
  grN0vsTime_PAMELA->Draw("pZ");
  grN0vsTime_BESS->Draw("pZ");
  grN0vsTime_BESSTeV->Draw("pZ");
  grN0vsTime_BESS00->Draw("pZ");
  // grN0vsTime_AMSPRL2015->Draw("pZ");
  grN0vsTime_AMS02->Draw("pZ");
  // grN0vsTime_AMS02Helium->Draw("pZl");
  gPad->SetGridy();


  ccBestFitK0VSTime->cd();
  ccBestFitK0VSTime->Update();



  // ---- fluxes vs time and p/He ratio ----

  // ---- draw results ---
  TDatime DateFlux1( 2011, 01, 01, 0, 0, 0);
  TDatime DateFlux2( 2018, 01, 01, 0, 0, 0);
  UTMIN= (double)DateFlux1.Convert();
  UTMAX= (double)DateFlux2.Convert();
  FMIN = 5.8;
  FMAX = 8.4;

  TH2F* hFramePHeRatiovsTime= new TH2F("hFramePHeRatiovsTime","K0 vs Time",200, UTMIN, UTMAX, 500, FMIN, FMAX);
  SetStyleHistoVSTime(hFramePHeRatiovsTime);
  hFramePHeRatiovsTime->GetYaxis()->SetNdivisions(506);
  hFramePHeRatiovsTime->GetYaxis()->SetTitleSize(0.07);
  hFramePHeRatiovsTime->GetYaxis()->SetTitle("p/He ratio");
  hFramePHeRatiovsTime->GetYaxis()->SetTitleOffset(0.50);


  // ---- plot p-He and p/He at selected rigidities ----

  int thisRig1 = thisRIG1;
  int thisRig2 = thisRIG2;
  int thisRig3 = thisRIG3;
  cout<<"TIME-PROFILES at Rlow: "<<  pRigAMS02Proton[thisRig1+5]<< "    Rmid: "<<pRigAMS02Helium[thisRig2]<<"    Rhigh: "<<pRigAMS02Helium[thisRig3]<<endl;


  // --- plot p-He fluxes and p/He ratio vs time ---
  TCanvas* ccResultsTimeProfiles= new TCanvas(Form("ccResultsTimeProfiles_N%d",0), Form("Flux Time Profiles | N%d",0), 1300, 730 );
  ccResultsTimeProfiles->Divide(3,3,0.01,0.01);  


  // up-left: proton LowRig
  ccResultsTimeProfiles->cd(1);
  gPad->SetBottomMargin(0.15);
  gPad->SetTopMargin(0.05);
  gPad->SetLeftMargin(0.150);
  gPad->SetRightMargin(0.06);

  grDataTimeProfileProton_AMS02_Rig[thisRig1+5]->Draw("apZ");
  grResultsTimeProfileProton_AMS02_Rig[thisRig1+5]->Draw("l");

  // up-mid: proton MidRig
  ccResultsTimeProfiles->cd(2);
  gPad->SetBottomMargin(0.15);
  gPad->SetTopMargin(0.05);
  gPad->SetLeftMargin(0.150);
  gPad->SetRightMargin(0.06);

  grDataTimeProfileProton_AMS02_Rig[thisRig2+5]->Draw("apZ");
  grResultsTimeProfileProton_AMS02_Rig[thisRig2+5]->Draw("l");

  // up-righ: proton HigRig
  ccResultsTimeProfiles->cd(3);
  gPad->SetBottomMargin(0.15);
  gPad->SetTopMargin(0.05);
  gPad->SetLeftMargin(0.150);
  gPad->SetRightMargin(0.06);

  grDataTimeProfileProton_AMS02_Rig[thisRig3+5]->Draw("apZ");
  grResultsTimeProfileProton_AMS02_Rig[thisRig3+5]->Draw("l");

  // mid-left: helium LowRig
  ccResultsTimeProfiles->cd(4);
  gPad->SetBottomMargin(0.15);
  gPad->SetTopMargin(0.05);
  gPad->SetLeftMargin(0.150);
  gPad->SetRightMargin(0.06);
  
  grDataTimeProfileHelium_AMS02_Rig[thisRig1]->Draw("apZ");
  grResultsTimeProfileHelium3_AMS02_Rig[thisRig1]->Draw("l");
  grResultsTimeProfileHelium4_AMS02_Rig[thisRig1]->Draw("l");

  // mid-mid: helium MidRig
  ccResultsTimeProfiles->cd(5);
  gPad->SetBottomMargin(0.15);
  gPad->SetTopMargin(0.05);
  gPad->SetLeftMargin(0.150);
  gPad->SetRightMargin(0.06);
  
  grDataTimeProfileHelium_AMS02_Rig[thisRig2]->Draw("apZ");
  grResultsTimeProfileHelium3_AMS02_Rig[thisRig2]->Draw("l");
  grResultsTimeProfileHelium4_AMS02_Rig[thisRig2]->Draw("l");

  // mid-right: helium HigRig
  ccResultsTimeProfiles->cd(6);
  gPad->SetBottomMargin(0.15);
  gPad->SetTopMargin(0.05);
  gPad->SetLeftMargin(0.150);
  gPad->SetRightMargin(0.06);
  
  grDataTimeProfileHelium_AMS02_Rig[thisRig3]->Draw("apZ");
  grResultsTimeProfileHelium3_AMS02_Rig[thisRig3]->Draw("l");
  grResultsTimeProfileHelium4_AMS02_Rig[thisRig3]->Draw("l");


  // bot-left phe ratio RigLow
  ccResultsTimeProfiles->cd(7);
  gPad->SetBottomMargin(0.15);
  gPad->SetTopMargin(0.05);
  gPad->SetLeftMargin(0.150);
  gPad->SetRightMargin(0.06);

  hFramePHeRatiovsTime->Draw();
  grDataTimeProfilePHeRatio_AMS02_Rig[thisRig1]->Draw("pZ");
  grResultsTimeProfilePHeRatio3_AMS02_Rig[thisRig1]->Draw("l");
  grResultsTimeProfilePHeRatio4_AMS02_Rig[thisRig1]->Draw("l");

  // bot-mid phe ratio RigMid
  ccResultsTimeProfiles->cd(8);
  gPad->SetBottomMargin(0.15);
  gPad->SetTopMargin(0.05);
  gPad->SetLeftMargin(0.150);
  gPad->SetRightMargin(0.06);

  hFramePHeRatiovsTime->Draw();
  grDataTimeProfilePHeRatio_AMS02_Rig[thisRig2]->Draw("pZ");
  grResultsTimeProfilePHeRatio3_AMS02_Rig[thisRig2]->Draw("l");
  grResultsTimeProfilePHeRatio4_AMS02_Rig[thisRig2]->Draw("l");

  // bot-right phe ratio RigHig
  ccResultsTimeProfiles->cd(9);
  gPad->SetBottomMargin(0.15);
  gPad->SetTopMargin(0.05);
  gPad->SetLeftMargin(0.150);
  gPad->SetRightMargin(0.06);

  hFramePHeRatiovsTime->Draw();
  grDataTimeProfilePHeRatio_AMS02_Rig[thisRig3]->Draw("pZ");
  grResultsTimeProfilePHeRatio3_AMS02_Rig[thisRig3]->Draw("l");
  grResultsTimeProfilePHeRatio4_AMS02_Rig[thisRig3]->Draw("l");

  /*
  TF1* f1= new TF1("f1","[2]+[0]*sin(([3]+[1]*x)/1.e+10)",grResultsTimeProfilePHeRatio4_AMS02_Rig[0]->GetX()[0],grResultsTimeProfilePHeRatio4_AMS02_Rig[0]->GetX()[79]);
  f1->SetLineColor(kRed+1);
  grDataTimeProfilePHeRatio_AMS02_Rig[thisRig1]->Fit(f1);

  TF1* f2= new TF1("f2","[2]+[0]*sin(([3]+[1]*x)/1.e+10)",grResultsTimeProfilePHeRatio4_AMS02_Rig[0]->GetX()[0],grResultsTimeProfilePHeRatio4_AMS02_Rig[0]->GetX()[79]);
  f2->SetLineColor(kRed+1);
  grDataTimeProfilePHeRatio_AMS02_Rig[thisRig2]->Fit(f2);
  */

  ccResultsTimeProfiles->cd();
  ccResultsTimeProfiles->Update();




  // --- plot for CC: just p/He ---
  TCanvas* ccPHeRatioResultsVSTime= new TCanvas(Form("ccPHeRatioResultsVSTime_N%d",0), Form("p/He Time Profiles | N%d",0), 800, 500);
  ccPHeRatioResultsVSTime->cd();
  gPad->SetBottomMargin(0.15);
  gPad->SetTopMargin(0.05);
  gPad->SetLeftMargin(0.150);
  gPad->SetRightMargin(0.06);

  hFramePHeRatiovsTime->Draw();
  grDataTimeProfilePHeRatio_AMS02_Rig[thisRig1]->Draw("pZ");
  grResultsTimeProfilePHeRatio4_AMS02_Rig[thisRig1]->Draw("l");
  grResultsTimeProfilePHeRatio3_AMS02_Rig[thisRig1]->Draw("l");

  grDataTimeProfilePHeRatio_AMS02_Rig[thisRig2]->Draw("pZ");
  grResultsTimeProfilePHeRatio4_AMS02_Rig[thisRig2]->Draw("l");
  grResultsTimeProfilePHeRatio3_AMS02_Rig[thisRig2]->Draw("l");

  grDataTimeProfilePHeRatio_AMS02_Rig[thisRig3]->Draw("pZ");
  grResultsTimeProfilePHeRatio4_AMS02_Rig[thisRig3]->Draw("l");
  grResultsTimeProfilePHeRatio3_AMS02_Rig[thisRig3]->Draw("l");


  // --- the end ---
  // theApp.Run();
  // return 1;
  


  // ---- frames ----
  TH2F* hFrameProtonFluxVSEkn= new TH2F("hFrameProtonFluxVSEkn","hFrameProtonFluxVSEkn",200, 0.01, 50.0, 200, 0.1, 1.e+4);
  SetStyleHistoVSEkn(hFrameProtonFluxVSEkn);
  hFrameProtonFluxVSEkn->GetXaxis()->SetTitle("kinetic energy (GeV)");
  hFrameProtonFluxVSEkn->GetYaxis()->SetNdivisions(506);
  hFrameProtonFluxVSEkn->GetYaxis()->SetTitle("J ( GeV^{ -1} m^{ -2} s^{ -1} sr^{ -1} )");

  
  // ---- plot all AMS-02 spectra ----
  TCanvas* ccResultsEnergySpectraAMS02= new TCanvas(Form("ccResultsEnergySpectraAMS02_N%d",0), "AMS-02 ENERGY SPECTRA", 1850, 1000 );
  ccResultsEnergySpectraAMS02->Divide(7,5,0.01,0.01);  

  indTime=0; // per gestire i due buchi
  for(int tt=0;tt<35;tt++){ // time loop
    if(tt==amsdata->SkipThis || tt==amsdata->SkipThis+1) continue;
    
    ccResultsEnergySpectraAMS02->cd(indTime+1);
    gPad->SetBottomMargin(0.10);
    gPad->SetTopMargin(0.02);
    gPad->SetLeftMargin(0.08);
    gPad->SetRightMargin(0.02);

    // proton
    hFrameProtonFluxVSEkn->Draw();
    gPad->SetLogx();
    gPad->SetLogy();
    amsdata->grAMS02_ProtonFluxVSEkn[tt]->SetMarkerSize(0.6);
    amsdata->grAMS02_ProtonFluxVSEkn[tt]->Draw("pZ"); // DATA
    grResultsFluxVSEkn_AMS02[tt]->Draw("l"); // MODEL

    // helium4
    amsdata->grAMS02_HeliumFluxVSEkn[tt]->SetMarkerSize(0.6);
    amsdata->grAMS02_HeliumFluxVSEkn[tt]->SetMarkerStyle(24);
    amsdata->grAMS02_HeliumFluxVSEkn[tt]->Draw("pZ"); // DATA converted to He(A=4)
    grResultsFluxVSEkn_He4_AMS02[tt]->Draw("l"); // MODEL He4
    grResultsFluxVSEkn_He3_AMS02[tt]->Draw("l"); // MODEL HeTOT?
    indTime++;
  }

  ccResultsEnergySpectraAMS02->cd();
  ccResultsEnergySpectraAMS02->Update();


  // ---- plot all [many] PAMELA spectra ----
  TCanvas* ccResultsEnergySpectraPAMELA= new TCanvas(Form("ccResultsEnergySpectraPAMELA_N%d",0), "PAMELA ENERGY SPECTRA", 1850, 1000 );
  ccResultsEnergySpectraPAMELA->Divide(7,5,0.01,0.01);  

  for(int tt=0;tt<35;tt++){ // time loop
    ccResultsEnergySpectraPAMELA->cd(tt+1);
    gPad->SetBottomMargin(0.10);
    gPad->SetTopMargin(0.02);
    gPad->SetLeftMargin(0.08);
    gPad->SetRightMargin(0.02);

    hFrameProtonFluxVSEkn->Draw();
    gPad->SetLogx();
    gPad->SetLogy();
    crdata->grPAMELA_ProtonFluxVSEkn[tt]->SetMarkerSize(0.6);
    crdata->grPAMELA_ProtonFluxVSEkn[tt]->Draw("pZ"); // DATA
    grResultsFluxVSEkn_PAMELA[tt]->Draw("l"); // MODEL
  }

  ccResultsEnergySpectraPAMELA->cd();
  ccResultsEnergySpectraPAMELA->Update();
  

  // --- the end ---
  theApp.Run();
  return 1;

  // YYYYYYYYY BOTTOM YYYYYYYYYYY
	  

  
}




void FCN(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag){
  double K0    = par[0]; 
  double XiD   = par[3]; 
  double N0    = par[4];

  // iExperiment: PAMELA / SOHO / BESS / BESSTeV / BESS00 / AMS02 / AMS02PRL2015
  int iExperiment   = (int)par[1]; 
  int iTime         = (int)par[2];

  f= GetChiSquare(K0, iTime, iExperiment, XiD, N0);

}



// ---- chisquare of spectrum | fixed epoch | given experiment | input K0 ----
double GetChiSquare(double K0, int iTime, int iExperiment, double XiD, double N0){

  double ChiSquare = 0.;
  double Emin = EMIN;
  double Emax = EMAX;
  double LnK0= log(K0);
  // double XiD = XiD
  
  // get data vs ekn of given exp at given time
  // get spectrum from TH2F-interpolation f(K0,E)
  // compute chisquare between the two

  // --- EXP1: PAMELA ----
  if( iExperiment==iExp_PAMELA){
    NDF_PAMELA[iTime]=0;
    for(int ee=0;ee<crdata->nEknPAMELA_Proton;ee++){ // ekn loop

      
      // get DATA
      double Energy= crdata->grPAMELA_ProtonFluxVSEkn[iTime]->GetX()[ee];
      if(Energy<Emin || Energy>Emax) continue;  // Check Range
	    
      double Ekn   = (1.e+3)*Energy;  // convert GeV->MeV
      double Flux  = (1.e-3)*crdata->grPAMELA_ProtonFluxVSEkn[iTime]->GetY()[ee];
      double eFlux = (1.e-3)*crdata->grPAMELA_ProtonFluxVSEkn[iTime]->GetEY()[ee];
      
      // get model
      double LnEkn= log(Ekn);
      double LnFluxModel= 0.;
      if(kDriftON)  LnFluxModel= hLn_Flux_vs_XiD_vs_KScale_vs_KEnergy_Z1_A1->Interpolate(LnEkn, LnK0, XiD);
      if(!kDriftON) LnFluxModel= hLn_Flux_vs_KScale_vs_KEnergy_Z1_A1->Interpolate(LnEkn, LnK0);
      double FluxModel= N0*exp(LnFluxModel);
      
      // get CHI2
      double Delta = (Flux - FluxModel)/eFlux;
      ChiSquare += (Delta*Delta);
      
      //cout<<"Time: "<<iTime<<"     Chi: "<<ChiSquare<<"     K0: "<<K0<<endl;      
    
      NDF_PAMELA[iTime]++;
    }
    if( bestChi2_PAMELA[iTime] > ChiSquare ) bestChi2_PAMELA[iTime] = ChiSquare;
  }

  
  // --- EXP2: SOHO ----
  if( iExperiment==iExp_SOHO){
    NDF_SOHO[iTime]=0;
    for(int ee=0;ee<crdata->nEknSOHO_Proton;ee++){ // ekn loop

      // get DATA
      double Energy= crdata->grSOHO_ProtonFluxVSEkn[iTime]->GetX()[ee];
      if(Energy<Emin || Energy>Emax) continue;  // Check Range

      double Ekn   = (1.e+3)*Energy;  // convert GeV->MeV
      double Flux  = (1.e-3)*crdata->grSOHO_ProtonFluxVSEkn[iTime]->GetY()[ee];
      double eFlux = (1.e-3)*crdata->grSOHO_ProtonFluxVSEkn[iTime]->GetEY()[ee];
      
      // get model
      double LnEkn= log(Ekn);
      double LnFluxModel= 0;
      if(kDriftON)  LnFluxModel= hLn_Flux_vs_XiD_vs_KScale_vs_KEnergy_Z1_A1->Interpolate(LnEkn, LnK0, XiD);
      if(!kDriftON) LnFluxModel= hLn_Flux_vs_KScale_vs_KEnergy_Z1_A1->Interpolate(LnEkn, LnK0);
      double FluxModel= N0*exp(LnFluxModel);
      
      // get CHI2
      double Delta = (Flux - FluxModel)/eFlux;
      ChiSquare += (Delta*Delta);
      
      //cout<<"Time: "<<iTime<<"     Chi: "<<ChiSquare<<"     K0: "<<K0<<endl;      
    
      NDF_SOHO[iTime]++;
    }
    if( bestChi2_SOHO[iTime] > ChiSquare ) bestChi2_SOHO[iTime] = ChiSquare;
  }


  // --- EXP3: BESS [Polar I/II] ----
  if( iExperiment==iExp_BESS){
    NDF_BESS[iTime]=0;
    for(int ee=0;ee<crdata->nEknBESS_Proton;ee++){ // ekn loop

      // get DATA
      double Energy= crdata->grBESS_ProtonFluxVSEkn[iTime]->GetX()[ee];
      if(Energy<Emin || Energy>Emax) continue;  // Check Range

      double Ekn   = (1.e+3)*Energy;  // convert GeV->MeV
      double Flux  = (1.e-3)*crdata->grBESS_ProtonFluxVSEkn[iTime]->GetY()[ee];
      double eFlux = (1.e-3)*crdata->grBESS_ProtonFluxVSEkn[iTime]->GetEY()[ee];
      
      // get model
      double LnEkn= log(Ekn);
      double LnFluxModel= 0.;
      if(kDriftON) LnFluxModel= hLn_Flux_vs_XiD_vs_KScale_vs_KEnergy_Z1_A1->Interpolate(LnEkn, LnK0, XiD);
      if(!kDriftON) LnFluxModel= hLn_Flux_vs_KScale_vs_KEnergy_Z1_A1->Interpolate(LnEkn, LnK0);
      double FluxModel= N0*exp(LnFluxModel);
      
      // get CHI2
      double Delta = (Flux - FluxModel)/eFlux;
      ChiSquare += (Delta*Delta);
      
      //cout<<"Time: "<<iTime<<"     Chi: "<<ChiSquare<<"     K0: "<<K0<<endl;      
    
      NDF_BESS[iTime]++;
    }
    if( bestChi2_BESS[iTime] > ChiSquare ) bestChi2_BESS[iTime] = ChiSquare;
  }



  
  // --- EXP4: BESSTeV ----
  if( iExperiment==iExp_BESSTeV){
    NDF_BESSTeV[iTime]=0;
    for(int ee=0;ee<crdata->nEknBESSTeV_Proton;ee++){ // ekn loop

      // get DATA
      double Energy= crdata->grBESSTeV_ProtonFluxVSEkn[iTime]->GetX()[ee];
      if(Energy<Emin || Energy>Emax) continue;  // Check Range

      double Ekn   = (1.e+3)*Energy;  // convert GeV->MeV
      double Flux  = (1.e-3)*crdata->grBESSTeV_ProtonFluxVSEkn[iTime]->GetY()[ee];
      double eFlux = (1.e-3)*crdata->grBESSTeV_ProtonFluxVSEkn[iTime]->GetEY()[ee];
      
      // get model
      double LnEkn= log(Ekn);
      double LnFluxModel= 0.;
      if(kDriftON) LnFluxModel= hLn_Flux_vs_XiD_vs_KScale_vs_KEnergy_Z1_A1->Interpolate(LnEkn, LnK0, XiD);
      if(!kDriftON)LnFluxModel= hLn_Flux_vs_KScale_vs_KEnergy_Z1_A1->Interpolate(LnEkn, LnK0);
      double FluxModel= N0*exp(LnFluxModel);
      
      // get CHI2
      double Delta = (Flux - FluxModel)/eFlux;
      ChiSquare += (Delta*Delta);
      
      //cout<<"Time: "<<iTime<<"     Chi: "<<ChiSquare<<"     K0: "<<K0<<endl;      
    
      NDF_BESSTeV[iTime]++;
    }
    if( bestChi2_BESSTeV[iTime] > ChiSquare ) bestChi2_BESSTeV[iTime] = ChiSquare;
  }


  // --- EXP5: BESS00 ----
  if( iExperiment==iExp_BESS00){
    NDF_BESS00[iTime]=0;
    for(int ee=0;ee<crdata->nEknBESS00_Proton;ee++){ // ekn loop

      // get DATA
      double Energy= crdata->grBESS00_ProtonFluxVSEkn[iTime]->GetX()[ee];
      if(Energy<Emin || Energy>Emax) continue;  // Check Range

      double Ekn   = (1.e+3)*Energy;  // convert GeV->MeV
      double Flux  = (1.e-3)*crdata->grBESS00_ProtonFluxVSEkn[iTime]->GetY()[ee];
      double eFlux = (1.e-3)*crdata->grBESS00_ProtonFluxVSEkn[iTime]->GetEY()[ee];
      
      // get model
      double LnEkn= log(Ekn);
      double LnFluxModel= 0.;
      if(kDriftON) LnFluxModel= hLn_Flux_vs_XiD_vs_KScale_vs_KEnergy_Z1_A1->Interpolate(LnEkn, LnK0, XiD);
      if(!kDriftON)LnFluxModel= hLn_Flux_vs_KScale_vs_KEnergy_Z1_A1->Interpolate(LnEkn, LnK0);
      double FluxModel= N0*exp(LnFluxModel);
      
      // get CHI2
      double Delta = (Flux - FluxModel)/eFlux;
      ChiSquare += (Delta*Delta);
      
      //cout<<"Time: "<<iTime<<"     Chi: "<<ChiSquare<<"     K0: "<<K0<<endl;      
    
      NDF_BESS00[iTime]++;
    }
    if( bestChi2_BESS00[iTime] > ChiSquare ) bestChi2_BESS00[iTime] = ChiSquare;
  }

  

  // --- EXP6: AMS02 ----
  if( iExperiment==iExp_AMS02){
    NDF_AMS02[iTime]=0;
    for(int ee=0;ee<amsdata->nEknAMS02_Proton;ee++){ // ekn loop

      // get DATA
      double Energy= amsdata->grAMS02_ProtonFluxVSEkn[iTime]->GetX()[ee];
      if(Energy<Emin || Energy>Emax) continue;  // Check Range

      double Ekn   = (1.e+3)*Energy;  // convert GeV->MeV
      double Flux  = (1.e-3)*amsdata->grAMS02_ProtonFluxVSEkn[iTime]->GetY()[ee];
      double eFlux = (1.e-3)*amsdata->grAMS02_ProtonFluxVSEkn[iTime]->GetEY()[ee];
      
      // get model
      double LnEkn= log(Ekn);
      double LnFluxModel= 0.;
      if(kDriftON)  LnFluxModel= hLn_Flux_vs_XiD_vs_KScale_vs_KEnergy_Z1_A1->Interpolate(LnEkn, LnK0, XiD);
      if(!kDriftON) LnFluxModel= hLn_Flux_vs_KScale_vs_KEnergy_Z1_A1->Interpolate(LnEkn, LnK0);
      double FluxModel= N0*exp(LnFluxModel);

      
      // get CHI2
      double Delta = (Flux - FluxModel)/eFlux;
      ChiSquare += (Delta*Delta);
      
      //cout<<"Time: "<<iTime<<"     Chi: "<<ChiSquare<<"     K0: "<<K0<<endl;      
    
      NDF_AMS02[iTime]++;
    }
    if( bestChi2_AMS02[iTime] > ChiSquare ) bestChi2_AMS02[iTime] = ChiSquare;
  }

  
  // --- EXP7: AMS02 PRL-2015 ---- ONE TIME POINT
  if( iExperiment==iExp_AMSPRL2015){
    NDF_AMSPRL2015[iTime]=0;
    for(int ee=0;ee<crdata->nEknAMS02_Proton;ee++){ // ekn loop

      // get DATA
      double Energy= crdata->grAMS02_ProtonFluxVSEkn->GetX()[ee];
      if(Energy<Emin || Energy>Emax) continue;  // Check Range

      double Ekn   = (1.e+3)*Energy;  // convert GeV->MeV
      double Flux  = (1.e-3)*crdata->grAMS02_ProtonFluxVSEkn->GetY()[ee];
      double eFlux = (1.e-3)*crdata->grAMS02_ProtonFluxVSEkn->GetEY()[ee];
      
      // get model
      double LnEkn= log(Ekn);
      double LnFluxModel= 0.;
      if(kDriftON) LnFluxModel= hLn_Flux_vs_XiD_vs_KScale_vs_KEnergy_Z1_A1->Interpolate(LnEkn, LnK0, XiD);
      if(!kDriftON)LnFluxModel= hLn_Flux_vs_KScale_vs_KEnergy_Z1_A1->Interpolate(LnEkn, LnK0);
      double FluxModel= N0*exp(LnFluxModel);
      
      // get CHI2
      double Delta = (Flux - FluxModel)/eFlux;
      ChiSquare += (Delta*Delta);
      
      //cout<<"Time: "<<iTime<<"     Chi: "<<ChiSquare<<"     K0: "<<K0<<endl;      
    
      NDF_AMSPRL2015[iTime]++;
    }
    if( bestChi2_AMSPRL2015[iTime] > ChiSquare ) bestChi2_AMSPRL2015[iTime] = ChiSquare;
  }

  
  return ChiSquare;
}


double Interpolation(TH3F* hh, double x, double y, double z){

  int NX = hh->GetXaxis()->GetNbins();
  int NY = hh->GetYaxis()->GetNbins();
  int NZ = hh->GetZaxis()->GetNbins();

  double Xmin = hh->GetXaxis()->GetBinCenter(1);
  double Xmax = hh->GetXaxis()->GetBinCenter(NX);
  double Ymin = hh->GetYaxis()->GetBinCenter(1);
  double Ymax = hh->GetYaxis()->GetBinCenter(NX);
  double Zmin = hh->GetZaxis()->GetBinCenter(1);
  double Zmax = hh->GetZaxis()->GetBinCenter(NX);

  double X = x;
  double Y = y;
  double Z = z;

  if(X < Xmin ) X = Xmin;
  if(Y < Ymin ) Y = Ymin;
  if(Z < Zmin ) Z = Zmin;

  if(X > Xmax ) X = Xmax;
  if(Y > Ymax ) Y = Ymax;
  if(Z > Zmax ) Z = Zmax;

  return hh->Interpolate(X, Y, Z);
}


double Interpolation(TH2F* hh, double x, double y){

  int NX = hh->GetXaxis()->GetNbins();
  int NY = hh->GetYaxis()->GetNbins();

  double Xmin = hh->GetXaxis()->GetBinCenter(1);
  double Xmax = hh->GetXaxis()->GetBinCenter(NX);
  double Ymin = hh->GetYaxis()->GetBinCenter(1);
  double Ymax = hh->GetYaxis()->GetBinCenter(NX);

  double X = x;
  double Y = y;

  if(X < Xmin ) X = Xmin;
  if(Y < Ymin ) Y = Ymin;

  if(X > Xmax ) X = Xmax;
  if(Y > Ymax ) Y = Ymax;

  return hh->Interpolate(X, Y);
}


void SetStyleHistoVSEkn(TH2F* hh){
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


void SetStyleHistoVSTime(TH2F* hh){
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
    gr->GetYaxis()->SetTitle("k_{0} ");
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










void SetGraphK0Style(){
  grK0vsTime_PAMELA->SetName("grK0VSTime_PAMELA");
  grK0vsTime_PAMELA->SetMarkerStyle( crdata->grPAMELA_ProtonFluxVSTime[0]->GetMarkerStyle() );
  grK0vsTime_PAMELA->SetMarkerColor( crdata->grPAMELA_ProtonFluxVSTime[0]->GetMarkerColor() );
  grK0vsTime_PAMELA->SetMarkerSize( crdata->grPAMELA_ProtonFluxVSTime[0]->GetMarkerSize() );
  grK0vsTime_PAMELA->SetLineWidth( crdata->grPAMELA_ProtonFluxVSTime[0]->GetLineWidth() );
  grK0vsTime_PAMELA->SetLineColor( crdata->grPAMELA_ProtonFluxVSTime[0]->GetLineColor() );
  SetTimeGraphStyle(grK0vsTime_PAMELA);

  grK0vsTime_SOHO->SetName("grK0VSTime_SOHO");
  grK0vsTime_SOHO->SetMarkerStyle( crdata->grSOHO_ProtonFluxVSTime[0]->GetMarkerStyle() );
  grK0vsTime_SOHO->SetMarkerColor( crdata->grSOHO_ProtonFluxVSTime[0]->GetMarkerColor() );
  grK0vsTime_SOHO->SetMarkerSize( crdata->grSOHO_ProtonFluxVSTime[0]->GetMarkerSize() );
  grK0vsTime_SOHO->SetLineWidth( crdata->grSOHO_ProtonFluxVSTime[0]->GetLineWidth() );
  grK0vsTime_SOHO->SetLineColor( crdata->grSOHO_ProtonFluxVSTime[0]->GetLineColor() );
  SetTimeGraphStyle(grK0vsTime_SOHO);

  grK0vsTime_BESS->SetName("grK0VSTime_BESS");
  grK0vsTime_BESS->SetMarkerStyle( crdata->grBESS_ProtonFluxVSTime[0]->GetMarkerStyle() );
  grK0vsTime_BESS->SetMarkerColor( crdata->grBESS_ProtonFluxVSTime[0]->GetMarkerColor() );
  grK0vsTime_BESS->SetMarkerSize( crdata->grBESS_ProtonFluxVSTime[0]->GetMarkerSize() );
  grK0vsTime_BESS->SetLineWidth( crdata->grBESS_ProtonFluxVSTime[0]->GetLineWidth() );
  grK0vsTime_BESS->SetLineColor( crdata->grBESS_ProtonFluxVSTime[0]->GetLineColor() );
  SetTimeGraphStyle(grK0vsTime_BESS);
  
  grK0vsTime_BESSTeV->SetName("grK0VSTime_BESSTeV");
  grK0vsTime_BESSTeV->SetMarkerStyle( crdata->grBESSTeV_ProtonFluxVSTime[0]->GetMarkerStyle() );
  grK0vsTime_BESSTeV->SetMarkerColor( crdata->grBESSTeV_ProtonFluxVSTime[0]->GetMarkerColor() );
  grK0vsTime_BESSTeV->SetMarkerSize( crdata->grBESSTeV_ProtonFluxVSTime[0]->GetMarkerSize() );
  grK0vsTime_BESSTeV->SetLineWidth( crdata->grBESSTeV_ProtonFluxVSTime[0]->GetLineWidth() );
  grK0vsTime_BESSTeV->SetLineColor( crdata->grBESSTeV_ProtonFluxVSTime[0]->GetLineColor() );
  SetTimeGraphStyle(grK0vsTime_BESSTeV);

  grK0vsTime_BESS00->SetName("grK0VSTime_BESS00");
  grK0vsTime_BESS00->SetMarkerStyle( crdata->grBESS00_ProtonFluxVSTime[0]->GetMarkerStyle() );
  grK0vsTime_BESS00->SetMarkerColor( crdata->grBESS00_ProtonFluxVSTime[0]->GetMarkerColor() );
  grK0vsTime_BESS00->SetMarkerSize( crdata->grBESS00_ProtonFluxVSTime[0]->GetMarkerSize() );
  grK0vsTime_BESS00->SetLineWidth( crdata->grBESS00_ProtonFluxVSTime[0]->GetLineWidth() );
  grK0vsTime_BESS00->SetLineColor( crdata->grBESS00_ProtonFluxVSTime[0]->GetLineColor() );
  SetTimeGraphStyle(grK0vsTime_BESS00);

  grK0vsTime_AMS02->SetName("grK0VSTime_AMS02");
  grK0vsTime_AMS02->SetMarkerStyle( amsdata->grAMS02_ProtonFluxVSTime[0]->GetMarkerStyle() );
  grK0vsTime_AMS02->SetMarkerColor( amsdata->grAMS02_ProtonFluxVSTime[0]->GetMarkerColor() );
  grK0vsTime_AMS02->SetMarkerSize( amsdata->grAMS02_ProtonFluxVSTime[0]->GetMarkerSize() );
  grK0vsTime_AMS02->SetLineWidth( amsdata->grAMS02_ProtonFluxVSTime[0]->GetLineWidth() );
  grK0vsTime_AMS02->SetLineColor( amsdata->grAMS02_ProtonFluxVSTime[0]->GetLineColor() );
  SetTimeGraphStyle(grK0vsTime_AMS02);

  grK0vsTime_AMSPRL2015->SetName("grK0VSTime_AMSPRL2015");
  grK0vsTime_AMSPRL2015->SetLineWidth( 3 );
  grK0vsTime_AMSPRL2015->SetLineColor( kRed+2 );
  grK0vsTime_AMSPRL2015->SetMarkerSize( 0 );
  grK0vsTime_AMSPRL2015->SetMarkerColor( kRed+2 );
  SetTimeGraphStyle(grK0vsTime_AMSPRL2015);
}



void SetGraphN0Style(){
  grN0vsTime_PAMELA->SetName("grN0VSTime_PAMELA");
  grN0vsTime_PAMELA->SetMarkerStyle( crdata->grPAMELA_ProtonFluxVSTime[0]->GetMarkerStyle() );
  grN0vsTime_PAMELA->SetMarkerColor( crdata->grPAMELA_ProtonFluxVSTime[0]->GetMarkerColor() );
  grN0vsTime_PAMELA->SetMarkerSize( crdata->grPAMELA_ProtonFluxVSTime[0]->GetMarkerSize() );
  grN0vsTime_PAMELA->SetLineWidth( crdata->grPAMELA_ProtonFluxVSTime[0]->GetLineWidth() );
  grN0vsTime_PAMELA->SetLineColor( crdata->grPAMELA_ProtonFluxVSTime[0]->GetLineColor() );
  SetTimeGraphStyle(grN0vsTime_PAMELA);

  grN0vsTime_SOHO->SetName("grN0VSTime_SOHO");
  grN0vsTime_SOHO->SetMarkerStyle( crdata->grSOHO_ProtonFluxVSTime[0]->GetMarkerStyle() );
  grN0vsTime_SOHO->SetMarkerColor( crdata->grSOHO_ProtonFluxVSTime[0]->GetMarkerColor() );
  grN0vsTime_SOHO->SetMarkerSize( crdata->grSOHO_ProtonFluxVSTime[0]->GetMarkerSize() );
  grN0vsTime_SOHO->SetLineWidth( crdata->grSOHO_ProtonFluxVSTime[0]->GetLineWidth() );
  grN0vsTime_SOHO->SetLineColor( crdata->grSOHO_ProtonFluxVSTime[0]->GetLineColor() );
  SetTimeGraphStyle(grN0vsTime_SOHO);

  grN0vsTime_BESS->SetName("grN0VSTime_BESS");
  grN0vsTime_BESS->SetMarkerStyle( crdata->grBESS_ProtonFluxVSTime[0]->GetMarkerStyle() );
  grN0vsTime_BESS->SetMarkerColor( crdata->grBESS_ProtonFluxVSTime[0]->GetMarkerColor() );
  grN0vsTime_BESS->SetMarkerSize( crdata->grBESS_ProtonFluxVSTime[0]->GetMarkerSize() );
  grN0vsTime_BESS->SetLineWidth( crdata->grBESS_ProtonFluxVSTime[0]->GetLineWidth() );
  grN0vsTime_BESS->SetLineColor( crdata->grBESS_ProtonFluxVSTime[0]->GetLineColor() );
  SetTimeGraphStyle(grN0vsTime_BESS);
  
  grN0vsTime_BESSTeV->SetName("grN0VSTime_BESSTeV");
  grN0vsTime_BESSTeV->SetMarkerStyle( crdata->grBESSTeV_ProtonFluxVSTime[0]->GetMarkerStyle() );
  grN0vsTime_BESSTeV->SetMarkerColor( crdata->grBESSTeV_ProtonFluxVSTime[0]->GetMarkerColor() );
  grN0vsTime_BESSTeV->SetMarkerSize( crdata->grBESSTeV_ProtonFluxVSTime[0]->GetMarkerSize() );
  grN0vsTime_BESSTeV->SetLineWidth( crdata->grBESSTeV_ProtonFluxVSTime[0]->GetLineWidth() );
  grN0vsTime_BESSTeV->SetLineColor( crdata->grBESSTeV_ProtonFluxVSTime[0]->GetLineColor() );
  SetTimeGraphStyle(grN0vsTime_BESSTeV);

  grN0vsTime_BESS00->SetName("grN0VSTime_BESS00");
  grN0vsTime_BESS00->SetMarkerStyle( crdata->grBESS00_ProtonFluxVSTime[0]->GetMarkerStyle() );
  grN0vsTime_BESS00->SetMarkerColor( crdata->grBESS00_ProtonFluxVSTime[0]->GetMarkerColor() );
  grN0vsTime_BESS00->SetMarkerSize( crdata->grBESS00_ProtonFluxVSTime[0]->GetMarkerSize() );
  grN0vsTime_BESS00->SetLineWidth( crdata->grBESS00_ProtonFluxVSTime[0]->GetLineWidth() );
  grN0vsTime_BESS00->SetLineColor( crdata->grBESS00_ProtonFluxVSTime[0]->GetLineColor() );
  SetTimeGraphStyle(grN0vsTime_BESS00);

  grN0vsTime_AMS02->SetName("grN0VSTime_AMS02");
  grN0vsTime_AMS02->SetMarkerStyle( amsdata->grAMS02_ProtonFluxVSTime[0]->GetMarkerStyle() );
  grN0vsTime_AMS02->SetMarkerColor( amsdata->grAMS02_ProtonFluxVSTime[0]->GetMarkerColor() );
  grN0vsTime_AMS02->SetMarkerSize( amsdata->grAMS02_ProtonFluxVSTime[0]->GetMarkerSize() );
  grN0vsTime_AMS02->SetLineWidth( amsdata->grAMS02_ProtonFluxVSTime[0]->GetLineWidth() );
  grN0vsTime_AMS02->SetLineColor( amsdata->grAMS02_ProtonFluxVSTime[0]->GetLineColor() );
  SetTimeGraphStyle(grN0vsTime_AMS02);

  grN0vsTime_AMSPRL2015->SetName("grN0VSTime_AMSPRL2015");
  grN0vsTime_AMSPRL2015->SetLineWidth( 3 );
  grN0vsTime_AMSPRL2015->SetLineColor( kRed+2 );
  grN0vsTime_AMSPRL2015->SetMarkerSize( 0 );
  grN0vsTime_AMSPRL2015->SetMarkerColor( kRed+2 );
  SetTimeGraphStyle(grN0vsTime_AMSPRL2015);
}



void SetGraphXiDStyle(){
  grXiDvsTime_PAMELA->SetName("grXiDVSTime_PAMELA");
  grXiDvsTime_PAMELA->SetMarkerStyle( crdata->grPAMELA_ProtonFluxVSTime[0]->GetMarkerStyle() );
  grXiDvsTime_PAMELA->SetMarkerColor( crdata->grPAMELA_ProtonFluxVSTime[0]->GetMarkerColor() );
  grXiDvsTime_PAMELA->SetMarkerSize( crdata->grPAMELA_ProtonFluxVSTime[0]->GetMarkerSize() );
  grXiDvsTime_PAMELA->SetLineWidth( crdata->grPAMELA_ProtonFluxVSTime[0]->GetLineWidth() );
  grXiDvsTime_PAMELA->SetLineColor( crdata->grPAMELA_ProtonFluxVSTime[0]->GetLineColor() );
  SetTimeGraphStyle(grXiDvsTime_PAMELA);

  grXiDvsTime_SOHO->SetName("grXiDVSTime_SOHO");
  grXiDvsTime_SOHO->SetMarkerStyle( crdata->grSOHO_ProtonFluxVSTime[0]->GetMarkerStyle() );
  grXiDvsTime_SOHO->SetMarkerColor( crdata->grSOHO_ProtonFluxVSTime[0]->GetMarkerColor() );
  grXiDvsTime_SOHO->SetMarkerSize( crdata->grSOHO_ProtonFluxVSTime[0]->GetMarkerSize() );
  grXiDvsTime_SOHO->SetLineWidth( crdata->grSOHO_ProtonFluxVSTime[0]->GetLineWidth() );
  grXiDvsTime_SOHO->SetLineColor( crdata->grSOHO_ProtonFluxVSTime[0]->GetLineColor() );
  SetTimeGraphStyle(grXiDvsTime_SOHO);

  grXiDvsTime_BESS->SetName("grXiDVSTime_BESS");
  grXiDvsTime_BESS->SetMarkerStyle( crdata->grBESS_ProtonFluxVSTime[0]->GetMarkerStyle() );
  grXiDvsTime_BESS->SetMarkerColor( crdata->grBESS_ProtonFluxVSTime[0]->GetMarkerColor() );
  grXiDvsTime_BESS->SetMarkerSize( crdata->grBESS_ProtonFluxVSTime[0]->GetMarkerSize() );
  grXiDvsTime_BESS->SetLineWidth( crdata->grBESS_ProtonFluxVSTime[0]->GetLineWidth() );
  grXiDvsTime_BESS->SetLineColor( crdata->grBESS_ProtonFluxVSTime[0]->GetLineColor() );
  SetTimeGraphStyle(grXiDvsTime_BESS);
  
  grXiDvsTime_BESSTeV->SetName("grXiDVSTime_BESSTeV");
  grXiDvsTime_BESSTeV->SetMarkerStyle( crdata->grBESSTeV_ProtonFluxVSTime[0]->GetMarkerStyle() );
  grXiDvsTime_BESSTeV->SetMarkerColor( crdata->grBESSTeV_ProtonFluxVSTime[0]->GetMarkerColor() );
  grXiDvsTime_BESSTeV->SetMarkerSize( crdata->grBESSTeV_ProtonFluxVSTime[0]->GetMarkerSize() );
  grXiDvsTime_BESSTeV->SetLineWidth( crdata->grBESSTeV_ProtonFluxVSTime[0]->GetLineWidth() );
  grXiDvsTime_BESSTeV->SetLineColor( crdata->grBESSTeV_ProtonFluxVSTime[0]->GetLineColor() );
  SetTimeGraphStyle(grXiDvsTime_BESSTeV);

  grXiDvsTime_BESS00->SetName("grXiDVSTime_BESS00");
  grXiDvsTime_BESS00->SetMarkerStyle( crdata->grBESS00_ProtonFluxVSTime[0]->GetMarkerStyle() );
  grXiDvsTime_BESS00->SetMarkerColor( crdata->grBESS00_ProtonFluxVSTime[0]->GetMarkerColor() );
  grXiDvsTime_BESS00->SetMarkerSize( crdata->grBESS00_ProtonFluxVSTime[0]->GetMarkerSize() );
  grXiDvsTime_BESS00->SetLineWidth( crdata->grBESS00_ProtonFluxVSTime[0]->GetLineWidth() );
  grXiDvsTime_BESS00->SetLineColor( crdata->grBESS00_ProtonFluxVSTime[0]->GetLineColor() );
  SetTimeGraphStyle(grXiDvsTime_BESS00);

  grXiDvsTime_AMS02->SetName("grXiDVSTime_AMS02");
  grXiDvsTime_AMS02->SetMarkerStyle( amsdata->grAMS02_ProtonFluxVSTime[0]->GetMarkerStyle() );
  grXiDvsTime_AMS02->SetMarkerColor( amsdata->grAMS02_ProtonFluxVSTime[0]->GetMarkerColor() );
  grXiDvsTime_AMS02->SetMarkerSize( amsdata->grAMS02_ProtonFluxVSTime[0]->GetMarkerSize() );
  grXiDvsTime_AMS02->SetLineWidth( amsdata->grAMS02_ProtonFluxVSTime[0]->GetLineWidth() );
  grXiDvsTime_AMS02->SetLineColor( amsdata->grAMS02_ProtonFluxVSTime[0]->GetLineColor() );
  SetTimeGraphStyle(grXiDvsTime_AMS02);

  grXiDvsTime_AMSPRL2015->SetName("grXiDVSTime_AMSPRL2015");
  grXiDvsTime_AMSPRL2015->SetLineWidth( 3 );
  grXiDvsTime_AMSPRL2015->SetLineColor( kRed+2 );
  grXiDvsTime_AMSPRL2015->SetMarkerSize( 0 );
  grXiDvsTime_AMSPRL2015->SetMarkerColor( kRed+2 );
  SetTimeGraphStyle(grXiDvsTime_AMSPRL2015);
}



  
