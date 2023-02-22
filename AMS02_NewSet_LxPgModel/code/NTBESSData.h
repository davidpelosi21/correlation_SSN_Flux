#ifndef NTBESSData_h
#define NTBESSData_h

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

class NTBESSData{

 public:

  // ---- number of TIME points ----
  static const int nTimeBESSTeV_Proton   = 1; 
  static const int nTimeBESS00_Proton    = 1; 
  static const int nTimeBESSPolar_Proton = 2; 
  
  // ---- number of EKN points ----
  static const int nEknBESSTeV_Proton   = 47; 
  static const int nEknBESS00_Proton    = 30;
  static const int nEknBESSPolar_Proton = 89;


  // ---- time arrays ----
  double xTimeBESSTeV_ProtonFlux[nTimeBESSTeV_Proton]; 
  double xTimeBESS00_ProtonFlux[nTimeBESS00_Proton]; 
  double xTimeBESSPolar_ProtonFlux[nTimeBESSPolar_Proton]; 

  // ---- time errors ----
  double eTimeBESSTeV_ProtonFlux[nTimeBESSTeV_Proton]; // NT array!
  double eTimeBESS00_ProtonFlux[nTimeBESS00_Proton];  // NT array!
  double eTimeBESSPolar_ProtonFlux[nTimeBESSPolar_Proton];
  
  // ---- energy arrays ----
  double xEknBESSTeV_ProtonFlux[nEknBESSTeV_Proton];
  double xEknBESS00_ProtonFlux[nEknBESS00_Proton];
  double xEknBESSPolar_ProtonFlux[nEknBESSPolar_Proton];

  
 public:

  // ---- proton flux VS ekn ----
  TGraphErrors* grBESSTeV_ProtonFluxVSEkn[nTimeBESSTeV_Proton];
  TGraphErrors* grBESS00_ProtonFluxVSEkn[nTimeBESS00_Proton];
  TGraphErrors* grBESSPolar_ProtonFluxVSEkn[nTimeBESSPolar_Proton];

  
  // ---- proton flux VS time ----
  TGraphErrors* grBESSTeV_ProtonFluxVSTime[nEknBESSTeV_Proton]; //1-point graph
  TGraphErrors* grBESS00_ProtonFluxVSTime[nEknBESS00_Proton]; //1-point graph
  TGraphErrors* grBESSPolar_ProtonFluxVSTime[nEknBESSPolar_Proton]; //2 points


  NTBESSData();
  void SetAllData();
  void SetProtonData();
  void SetHeliumData();


  void SymmetrizeGraph(TGraphAsymmErrors* grA, TGraphErrors* grS); 
  void ScaleGraphErrors(TGraphErrors* gr, double ScaleValue);
  void CleanXErrors(TGraphErrors* gr);

};


#endif
