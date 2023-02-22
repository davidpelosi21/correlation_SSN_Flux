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
  static const int nTimeBESS_Proton  = 8;  //bess + polar 

  // ---- number of EKN points ----
  //int nEknBESS_Proton[10]  = {14,14,14,30,41,30,30,47,89,89}; // R>1 GV [up 60GV: 45 | up 100GV: 52 | up 2TV: 74];
int nEknBESS_Proton[8]  = {14,30,41,30,30,47,89,89}; 


  // ---- time arrays ----
  
   double xTimeBESS_ProtonFlux[nTimeBESS_Proton]; 


  // ---- time errors ----
  double eTimeBESS_ProtonFlux;

  
  // ---- energy arrays ----
  //double xEknBESS_ProtonFlux[4][nEknBESS_Proton];



 public:
  // ---- proton flux VS ekn ----
  TGraphErrors* grBESS_ProtonFluxVSEkn[nTimeBESS_Proton];


  TGraphErrors* grVOYAGER1_ProtonFluxVSEkn;


  
  // ---- proton flux VS time -------
  TGraphErrors *grBESS_ProtonFluxVSTime = new TGraphErrors(); //1-point graph



  NTBESSData();
  void SetAllData();
  void SetProtonData();



  //// void SymmetrizeGraph(TGraphAsymmErrors* grA, TGraphErrors* grS); // USED?
  void ScaleGraphErrors(TGraphErrors* gr, double ScaleValue);
  void CleanXErrors(TGraphErrors* gr);

};


#endif
