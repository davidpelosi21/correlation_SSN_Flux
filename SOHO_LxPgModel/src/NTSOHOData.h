#ifndef NTSOHOData_h
#define NTSOHOData_h

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

class NTSOHOData{

 public:

  int kTimeUnits;
  
  // ---- number of TIME points ----
  static const int nTimeSOHO_Proton  = 20;

  // ---- number of EKN points ----
  static const int nEknSOHO_Proton   = 13;


  // ---- time arrays ----
  double xTimeSOHO_ProtonFlux[nTimeSOHO_Proton]; 

  // ---- time errors ----
  double eTimeSOHO_ProtonFlux[nTimeSOHO_Proton]; 
  
  // ---- energy arrays ----
  double xEknSOHO_ProtonFlux[nEknSOHO_Proton];

  // ---- proton flux VS ekn ----
  TGraphErrors* grSOHO_ProtonFluxVSEkn[nTimeSOHO_Proton];
  
  // ---- proton flux VS time ----
  TGraphErrors* grSOHO_ProtonFluxVSTime[nEknSOHO_Proton]; 


  // ---- basic class ----
  NTSOHOData();
  NTSOHOData(int TimeUnits);
  void SetAllData();
  void SetProtonData();

  void ScaleGraphErrors(TGraphErrors* gr, double ScaleValue);
  void CleanXErrors(TGraphErrors* gr);
  void RemoveZero(TGraphErrors* gr);
  void ConvertUnixTime2FractYear(TGraphErrors* gr);

};


#endif
