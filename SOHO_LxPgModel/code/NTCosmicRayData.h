#ifndef NTCosmicRayData_h
#define NTCosmicRayData_h

#include "TObject.h"
#include "TROOT.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TStyle.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"

#include "TMultiGraph.h"
#include "TMath.h"
#include <iostream>


using namespace std;

class NTCosmicRayData{

 public:


  // ---- number of TIME points ----
  static const int nTimePAMELA_Proton = 47;
  static const int nTimeSOHO_Proton   = 20;
  static const int nTimeBESS_Proton   = 2; //Polar I/II [no 00/TeV]
  static const int nTimeAMS02_Proton  = 1; 
  static const int nTimePAMELA_Electron     = 7;
  static const int nTimePAMELA_LeptonRatio  = 35;
  

  // ---- number of EKN points ----
  static const int nEknPAMELA_Proton   = 78;
  static const int nEknSOHO_Proton     = 13;
  static const int nEknBESS_Proton     = 89;
  static const int nEknAMS02_Proton    = 72;
  static const int nEknVOYAGER1_Proton = 15;

  
  static const int nEknVOYAGER1_AllElectron = 4; 
  static const int nEknPAMELA_Electron      = 26;
  static const int nEknAMS02_Electron       = 73;
  static const int nEknAMS02_Positron       = 73;

  static const int nEknPAMELA_LeptonRatio      =  3;
  static const int nEknPAMELA_PositronFraction = 16;
  static const int nEknAMS02_PositronFraction  = 66;
  static const int nEknAMS01_PositronFraction  =  9;

  

  // ---- time arrays ----
  double xTimePAMELA_ProtonFlux[nTimePAMELA_Proton];
  double xTimeSOHO_ProtonFlux[nTimeSOHO_Proton];
  double xTimeBESS_ProtonFlux[nTimeBESS_Proton];
  double xTimeAMS02_ProtonFlux[nTimeAMS02_Proton]; // one point

  double xTimePAMELA_ElectronFlux[nTimePAMELA_Electron];
  double xTimePAMELA_LeptonRatio[nTimePAMELA_LeptonRatio];
  double xTimeAMS02_ElectronFlux[1];
  double xTimeAMS02_PositronFlux[1];

  double xTimePAMELA_PositronFraction[1];
  double xTimeAMS02_PositronFraction[1];
  double xTimeAMS01_PositronFraction[1];

  double xTimeAMS02_PbarPRatio[1];
  double xTimePAMELA_PbarPRatio[1];
  double xTimeBESSPolarI_PbarPRatio[1];
  double xTimeBESSPolarII_PbarPRatio[1];
  double xTimeBESSTeV_PbarPRatio[1];

  // ---- time errors ----

  double eTimePAMELA_ProtonFlux;
  double eTimeBESS_ProtonFlux[2];
  double eTimeSOHO_ProtonFlux;
  double eTimeAMS02_ProtonFlux;

  double eTimePAMELA_ElectronFlux;
  double eTimeAMS02_ElectronFlux;
  double eTimeAMS02_PositronFlux;

  double eTimePAMELA_PositronFraction;
  double eTimeAMS02_PositronFraction;
  double eTimeAMS01_PositronFraction;
  double eTimePAMELA_LeptonRatio;

  double eTimeAMS02_PbarPRatio;
  double eTimePAMELA_PbarPRatio;
  double eTimeBESSPolarI_PbarPRatio;
  double eTimeBESSPolarII_PbarPRatio;
  double eTimeBESSTeV_PbarPRatio;
  
  
  // ---- energy arrays ----
  double xEknPAMELA_ProtonFlux[nEknPAMELA_Proton];
  double xEknSOHO_ProtonFlux[nEknSOHO_Proton];
  double xEknBESS_ProtonFlux[nEknBESS_Proton];
  double xEknAMS02_ProtonFlux[nEknAMS02_Proton];

  double xEknPAMELA_ElectronFlux[nEknPAMELA_Electron];
  double xEknAMS02_ElectronFlux[nEknAMS02_Electron];
  double xEknVOYAGER1_AllElectronFlux[nEknVOYAGER1_AllElectron];

  double xEknPAMELA_LeptonRatio[nEknPAMELA_Electron];
  double xEknPAMELA_PositronFraction[nEknPAMELA_PositronFraction];
  double xEknAMS02_PositronFraction[nEknAMS02_PositronFraction];
  double xEknAMS01_PositronFraction[nEknAMS01_PositronFraction];
  double xEknAMS02_PositronFlux[nEknAMS02_Positron];

  
  // NT2017 OTHER BESS PROTON SPECTRA
  static const int nTimeBESSTeV_Proton = 1;
  static const int nEknBESSTeV_Proton = 47;
  double xTimeBESSTeV_ProtonFlux[nTimeBESSTeV_Proton];
  double eTimeBESSTeV_ProtonFlux[nTimeBESSTeV_Proton];
  double xEknBESSTeV_ProtonFlux[nEknBESSTeV_Proton];
  
  static const int nTimeBESS00_Proton = 1;
  static const int nEknBESS00_Proton  = 30;
  double xTimeBESS00_ProtonFlux[nTimeBESS00_Proton];
  double eTimeBESS00_ProtonFlux[nTimeBESS00_Proton];
  double xEknBESS00_ProtonFlux[nEknBESS00_Proton];

  
 public:

  // ---- proton flux VS ekn ----
  TGraphErrors* grPAMELA_ProtonFluxVSEkn[nTimePAMELA_Proton];
  TGraphErrors* grSOHO_ProtonFluxVSEkn[nTimeSOHO_Proton];
  TGraphErrors* grBESS_ProtonFluxVSEkn[nTimeBESS_Proton];
  TGraphErrors* grAMS02_ProtonFluxVSEkn;
  TGraphErrors* grVOYAGER1_ProtonFluxVSEkn;

  TGraphErrors* grBESSTeV_ProtonFluxVSEkn[nTimeBESSTeV_Proton];
  TGraphErrors* grBESS00_ProtonFluxVSEkn[nTimeBESS00_Proton];

  TMultiGraph* grALL_ProtonFluxVSEkn;


  // ---- proton flux VS time ----
  TGraphErrors* grPAMELA_ProtonFluxVSTime[nEknPAMELA_Proton];
  TGraphErrors* grSOHO_ProtonFluxVSTime[nEknSOHO_Proton];
  TGraphErrors* grBESS_ProtonFluxVSTime[nEknBESS_Proton];
  TGraphErrors* grAMS02_ProtonFluxVSTime[nEknAMS02_Proton]; //1-point graph

  TGraphErrors* grBESSTeV_ProtonFluxVSTime[nEknBESSTeV_Proton]; //1-point graph
  TGraphErrors* grBESS00_ProtonFluxVSTime[nEknBESS00_Proton]; //1-point graph

  TMultiGraph* grALL_ProtonFluxVSTime;


  // ---- electron flux VS ekn ----
  TGraphErrors* grPAMELA_ElectronFluxVSEkn[nTimePAMELA_Electron];
  TGraphErrors* grAMS02_ElectronFluxVSEkn; 
  TGraphErrors* grVOYAGER1_AllElectronFluxVSEkn;
  TMultiGraph* grALL_ElectronFluxVSEkn;

  // ---- antiproton/proton ratio VS ekn ----
  TGraphErrors* grAMS02_PbarPRatioVSEkn;     // Aguilar et al. 2016 PRL
  TGraphErrors* grPAMELA_PbarPRatioVSEkn;    // Adriani et al. 2010 PRL 121101 
  TGraphErrors* grBESSPolarI_PbarPRatioVSEkn; // Abe et al., PLB 2008 670, 103
  TGraphErrors* grBESSPolarII_PbarPRatioVSEkn; // Abe et al., PRL 2012
  TGraphErrors* grBESSTeV_PbarPRatioVSEkn;   // Haino et al., 2005 ICRC 3, 13
  TMultiGraph* grALL_PbarPRatioVSEkn;
  
  
  //TGraphAsymmErrors* grHEATPbar_PbarPRatioVSEkn_2001;    // Beach et al., 2001 PRL 87, 271101  

  
  // ---- positron fraction VS ekn ----
  TGraphErrors* grPAMELA_PositronFractionVSEkn;
  TGraphErrors* grAMS02_PositronFractionVSEkn;
  TGraphErrors* grAMS01_PositronFractionVSEkn;
  TMultiGraph*  grALL_PositronFractionVSEkn;

  // ---- positron flux ----
  TGraphErrors* grAMS02_PositronFluxVSEkn; 


  // ---- leptons VS time ----
  TGraphErrors* grPAMELA_ElectronFluxVSTime[nEknPAMELA_Electron];
  TGraphErrors* grPAMELA_LeptonRatioVSTime[nEknPAMELA_LeptonRatio];

  NTCosmicRayData();

  void SetAllData();
  void SetProtonData();
  void SetElectronData();
  void SetPositronData();
  void SetAntiprotonData();

  void SymmetrizeGraph(TGraphAsymmErrors* grA, TGraphErrors* grS);
  void ScaleGraphErrors(TGraphErrors* gr, double ScaleValue);
  void CleanXErrors(TGraphErrors* gr);

  //SetStyleHistoVSEkn(TH2F* hh);
  //SetStyleHistoVSTime(TH2D* hh);

};


#endif
