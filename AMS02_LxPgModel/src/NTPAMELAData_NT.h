#ifndef NTPAMELAData_h
#define NTPAMELAData_h

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

class NTPAMELAData{

 public:

  // ---- number of TIME points ----
  static const int nTimePAMELA_Proton  = 46; // 79+2 with UTTPS off: 81

  // ---- number of EKN points ----
  static const int nEknPAMELA_Proton   = 78; // R>1 GV [up 60GV: 45 | up 100GV: 52 | up 2TV: 74];

  // ---- time arrays ----
  double xTimePAMELA_ProtonFlux[nTimePAMELA_Proton]; 

  // ---- time errors ----
  double eTimePAMELA_ProtonFlux;

  
  // ---- energy arrays ----
  double xEknPAMELA_ProtonFlux[nEknPAMELA_Proton];


 public:
  // ---- proton flux VS ekn ----
  TGraphErrors* grPAMELA_ProtonFluxVSEkn[nTimePAMELA_Proton];

  // ---- proton flux VS time -------
  TGraphErrors* grPAMELA_ProtonFluxVSTime[nEknPAMELA_Proton]; //1-point graph


  NTPAMELAData();
  void SetAllData();
  void SetProtonData();


  //// void SymmetrizeGraph(TGraphAsymmErrors* grA, TGraphErrors* grS); // USED?
  void ScaleGraphErrors(TGraphErrors* gr, double ScaleValue);
  void CleanXErrors(TGraphErrors* gr);

};


#endif
