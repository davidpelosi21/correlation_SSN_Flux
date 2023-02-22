#ifndef NTAMS02Data_h
#define NTAMS02Data_h

#include "TObject.h"
#include "TROOT.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TStyle.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h" // USEd?

#include "TMultiGraph.h"
#include "TMath.h"
#include <iostream>


using namespace std;

class NTAMS02Data{

 public:

  // ---- number of TIME points ----
  static const int nTimeAMS02_Proton  = 81; // 79+2 with UTTPS off: 81
  static const int nTimeAMS02_Helium  = 81; // 79+2 with UTTPS off: 81

  // ---- number of EKN points ----
  static const int nEknAMS02_Proton   = 45; // R>1 GV [up 60GV: 45 | up 100GV: 52 | up 2TV: 74];


  static const int nEknAMS02_Helium   = 40; // R>2 GV [up 60GV: 40 | up 100GV: 47 | up 3TV: 71]; 

  static const int nEknVOYAGER1_Proton = 15;
  static const int nEknVOYAGER1_Helium = 16;

  // ---- time arrays ----
  double xTimeAMS02_ProtonFlux[nTimeAMS02_Proton]; 
  double xTimeAMS02_HeliumFlux[nTimeAMS02_Helium]; 

  // ---- time errors ----
  double eTimeAMS02_ProtonFlux;
  double eTimeAMS02_HeliumFlux;
  
  // ---- energy arrays ----
  double xEknAMS02_ProtonFlux[nEknAMS02_Proton];
  double xEknAMS02_HeliumFlux[nEknAMS02_Helium];

  double xEknVOYAGER1_ProtonFlux[nEknVOYAGER1_Proton];
  double xEknVOYAGER1_HeliumFlux[nEknVOYAGER1_Helium];

  // ----- skip bartel ----
  int SkipThis; // same for p and He
  
 public:
  // ---- proton flux VS ekn ----
  TGraphErrors* grAMS02_ProtonFluxVSEkn[nTimeAMS02_Proton];
  TGraphErrors* grAMS02_HeliumFluxVSEkn[nTimeAMS02_Helium];

  TGraphErrors* grVOYAGER1_ProtonFluxVSEkn;
  TGraphErrors* grVOYAGER1_HeliumFluxVSEkn;

  
  // ---- proton flux VS time -------
  TGraphErrors* grAMS02_ProtonFluxVSTime[nEknAMS02_Proton]; //1-point graph
  TGraphErrors* grAMS02_HeliumFluxVSTime[nEknAMS02_Helium]; //1-point graph


  NTAMS02Data();
  void SetAllData();
  void SetProtonData();
  void SetHeliumData();


  //// void SymmetrizeGraph(TGraphAsymmErrors* grA, TGraphErrors* grS); // USED?
  void ScaleGraphErrors(TGraphErrors* gr, double ScaleValue);
  void CleanXErrors(TGraphErrors* gr);

};


#endif
