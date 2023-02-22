#include "NTAMS02Data.h"

//TString PATH_02 = "/Users/davidpelosi/Desktop/TESI-Master-current/Fitting_LxPg_2";
TString PATH_02 = "./";

NTAMS02Data::NTAMS02Data(){
}


void NTAMS02Data::SetAllData(){
  SetProtonData();
  SetHeliumData();
}


void NTAMS02Data::SetProtonData(){
  SkipThis = 46; // skip this and this+1

  // ---- pick up data files ----
TFile* inFileAMS02= new TFile(PATH_02+"/data/grAMS02_ProtonHelium_vs_Ekn_vs_Time_DEC2017.root","READ");

inFileAMS02->cd();
inFileAMS02->ls();
 
  // --- get fluxes VS ekn ---
  for(int tt=0;tt<nTimeAMS02_Proton;tt++){
    if(tt==SkipThis || tt==SkipThis+1){
      grAMS02_ProtonFluxVSEkn[tt]= new TGraphErrors(); // empty
    }
    else{    
      grAMS02_ProtonFluxVSEkn[tt]=(TGraphErrors*)inFileAMS02->Get(Form("grProton_vs_Ekn_T%d",tt));
      grAMS02_ProtonFluxVSEkn[tt]->SetName(Form("grAMS02_ProtonFluxVSEkn_T%d",tt));
      grAMS02_ProtonFluxVSEkn[tt]->SetMarkerColor(kBlack);
      grAMS02_ProtonFluxVSEkn[tt]->SetLineColor(kBlack);
      grAMS02_ProtonFluxVSEkn[tt]->SetMarkerStyle(20);
      grAMS02_ProtonFluxVSEkn[tt]->SetMarkerSize(0.9);
      grAMS02_ProtonFluxVSEkn[tt]->SetLineWidth(1);
      //ScaleGraphErrors(grAMS02_ProtonFluxVSEkn[tt], 1.0); 
    }
    
    // cout<<"Proton Spectrum N. "<<tt<<" ... done!"<<endl;
  }

  cout<<"PROTON EKN.. SPECTRA DONE!"<<endl;
  
   TFile* inFileAMS02_1= new TFile(PATH_02+"/data/grAMS02_ProtonHelium_vs_Ekn_vs_Time_DEC2017.root","READ");
   
  // ---- get fluxes vs time ---- | NT BUT they are RIGIDITY fluxes!!
  for(int ee=0;ee<nEknAMS02_Proton;ee++){
    grAMS02_ProtonFluxVSTime[ee]=(TGraphErrors*)inFileAMS02_1->Get(Form("grProton_vs_Time_R%d",ee+1));
    grAMS02_ProtonFluxVSTime[ee]->SetName(Form("grAMS02_ProtonFluxVSTime_E%d",ee));
    grAMS02_ProtonFluxVSTime[ee]->SetMarkerColor(kBlack);
    grAMS02_ProtonFluxVSTime[ee]->SetLineColor(kBlack);
    grAMS02_ProtonFluxVSTime[ee]->SetMarkerStyle(20);
    grAMS02_ProtonFluxVSTime[ee]->SetMarkerSize(0.9);
    grAMS02_ProtonFluxVSTime[ee]->SetLineWidth(2);
    //ScaleGraphErrors(grAMS02_ProtonFluxVSTime[ee], 1.0);

    // cout<<"Proton profile N. "<<ee<<" ... done!"<<endl;
  }
  cout<<" PROTON TIME-PROFILES DONE!"<<endl;

  
  inFileAMS02->Close();
  delete inFileAMS02;

  // ---- get time array ---- 
  double* xTimeAMS02= (double*)grAMS02_ProtonFluxVSTime[0]->GetX();
  for(int tt=0;tt<nTimeAMS02_Proton;tt++)xTimeAMS02_ProtonFlux[tt]=xTimeAMS02[tt];
  eTimeAMS02_ProtonFlux= 0.5*(xTimeAMS02_ProtonFlux[1]-xTimeAMS02_ProtonFlux[0]);
  
  // ---- get ekn array ----

  double* xEknAMS02= (double*)grAMS02_ProtonFluxVSEkn[0]->GetX();
  for(int ee=0;ee<nEknAMS02_Proton;ee++)xEknAMS02_ProtonFlux[ee]=xEknAMS02[ee];


  
  // ----- SET VOYAGER-1 proton data ---- [ONLY VS EKN]

  double Voyager1_Hydrogen_Ekn[15]={3.8,5.4,6.9,10.2,15.4,23.9,39.0,52.0,79.1,143.9,164.9,181.3,204.1,245.3,308.0};
  double Voyager1_Hydrogen_Flux[15]={1.908e+01,2.116e+01,2.369e+01,2.665e+01,2.968e+01,2.909e+01,2.829e+01,2.914e+01,2.443e+01,1.927e+01,1.820e+01,1.652e+01,1.492e+01,1.289e+01,1.042e+01};
  double stat_Voyager1_Hydrogen_Flux[15]={9.572e-01,1.061e+00,1.188e+00,1.345e+00,1.496e+00,1.460e+00,1.419e+00,1.468e+00,1.263e+00,9.723e-01,9.197e-01,8.410e-01,7.519e-01,6.484e-01,5.236e-01};
  double syst_Voyager1_Hydrogen_Flux[15]={9.540e-01,1.058e+00,1.184e+00,1.333e+00,1.484e+00,1.454e+00,1.415e+00,1.457e+00,1.221e+00,9.635e-01,9.100e-01,8.260e-01,7.460e-01,6.445e-01,5.210e-01};
  
  for(int ee=0;ee<15;ee++) { // unit conversion
    Voyager1_Hydrogen_Ekn[ee] *= 1.e-3;
    Voyager1_Hydrogen_Flux[ee] *= 1.e+3;
    stat_Voyager1_Hydrogen_Flux[ee] *= 1.e+3;
    syst_Voyager1_Hydrogen_Flux[ee] *= 1.e+3;
  }

    double eVoyager1_Hydrogen_Ekn[15];
    double eVoyager1_Hydrogen_Flux[15];
    for(int ee=0;ee<15;ee++){
      eVoyager1_Hydrogen_Ekn[ee] = 0.;
      eVoyager1_Hydrogen_Flux[ee] = TMath::Sqrt(stat_Voyager1_Hydrogen_Flux[ee]*stat_Voyager1_Hydrogen_Flux[ee] + syst_Voyager1_Hydrogen_Flux[ee]*syst_Voyager1_Hydrogen_Flux[ee]);
    }
    
    grVOYAGER1_ProtonFluxVSEkn = new TGraphErrors(15, Voyager1_Hydrogen_Ekn, Voyager1_Hydrogen_Flux, eVoyager1_Hydrogen_Ekn, eVoyager1_Hydrogen_Flux);
    grVOYAGER1_ProtonFluxVSEkn->SetName("grVOYAGER1_ProtonFluxVSEkn");
    grVOYAGER1_ProtonFluxVSEkn->SetMarkerStyle(21);
    grVOYAGER1_ProtonFluxVSEkn->SetMarkerSize(1.1);
    grVOYAGER1_ProtonFluxVSEkn->SetMarkerColor(kBlack);
    grVOYAGER1_ProtonFluxVSEkn->SetLineColor(kBlack);


}



void NTAMS02Data::SetHeliumData(){
  SkipThis = 46; // skip this and this+1

  // ---- pick up data files ---- SAME FILE AS PROTON
  // TFile* inFileAMS02= new TFile("$NTBASEDIR/THEORY/00matisse/ANALYSIS/ExperimentalData/grAMS02_ProtonHelium_Fluxes_VS_Ekn_2017_RangePaper.root","READ");
          
  TFile* inFileAMS02= new TFile(PATH_02+"/data/grAMS02_ProtonHelium_vs_Ekn_vs_Time_DEC2017.root","READ");

  inFileAMS02->cd();
  
  
  // --- get fluxes VS ekn ---
  for(int tt=0;tt<nTimeAMS02_Helium;tt++){
    if(tt==SkipThis || tt==SkipThis+1){
      grAMS02_HeliumFluxVSEkn[tt]= new TGraphErrors(); // empty
    }
    else{
      grAMS02_HeliumFluxVSEkn[tt]=(TGraphErrors*)inFileAMS02->Get(Form("grHelium_vs_Ekn_T%d",tt));
      grAMS02_HeliumFluxVSEkn[tt]->SetName(Form("grAMS02_HeliumFluxVSEkn_T%d",tt));
      grAMS02_HeliumFluxVSEkn[tt]->SetMarkerColor(kBlack);
      grAMS02_HeliumFluxVSEkn[tt]->SetLineColor(kBlack);
      grAMS02_HeliumFluxVSEkn[tt]->SetMarkerStyle(20);
      grAMS02_HeliumFluxVSEkn[tt]->SetMarkerSize(0.9);
      grAMS02_HeliumFluxVSEkn[tt]->SetLineWidth(1);
      //ScaleGraphErrors(grAMS02_HeliumFluxVSEkn[tt], 1.0); 
    }
  }

  cout<<"HELIUM EKN SPECTRA DONE!"<<endl;

  // ---- get fluxes vs time ---- | NT BUT they are RIGIDITY fluxes!!
  for(int ee=0;ee<nEknAMS02_Helium;ee++){
    grAMS02_HeliumFluxVSTime[ee]=(TGraphErrors*)inFileAMS02->Get(Form("grHelium_vs_Time_R%d",ee+6));
    grAMS02_HeliumFluxVSTime[ee]->SetName(Form("grAMS02_HeliumFluxVSTime_E%d",ee));
    grAMS02_HeliumFluxVSTime[ee]->SetMarkerColor(kBlack);
    grAMS02_HeliumFluxVSTime[ee]->SetLineColor(kBlack);
    grAMS02_HeliumFluxVSTime[ee]->SetMarkerStyle(24);
    grAMS02_HeliumFluxVSTime[ee]->SetMarkerSize(0.9);
    grAMS02_HeliumFluxVSTime[ee]->SetLineWidth(2);
    //ScaleGraphErrors(grAMS02_HeliumFluxVSTime[ee], 1.0);
  }
  cout<<"PROTON TIME-PROFILES DONE!"<<endl;

  inFileAMS02->Close();
  delete inFileAMS02;
  
  // ---- get time array ---- 
  double* xTimeAMS02= (double*)grAMS02_HeliumFluxVSTime[0]->GetX();
  for(int tt=0;tt<nTimeAMS02_Helium;tt++)xTimeAMS02_HeliumFlux[tt]=xTimeAMS02[tt];
  eTimeAMS02_HeliumFlux= 0.5*(xTimeAMS02_HeliumFlux[1]-xTimeAMS02_HeliumFlux[0]);
  
  // ---- get ekn array ----
  double* xEknAMS02= (double*)grAMS02_HeliumFluxVSEkn[0]->GetX();
  for(int ee=0;ee<nEknAMS02_Helium;ee++)xEknAMS02_HeliumFlux[ee]=xEknAMS02[ee];

  

  // ----- SET VOYAGER-1 helium data ---- [ONLY VS EKN]
  double Voyager1_Helium_Ekn[16]={3.8,5.4,6.9,10.2,15.4,23.9,39.0,52.0,122.0,142.8,172.4,216.0,278.3,349.5,430.9,571.0};
  double Voyager1_Helium_Flux[16]={1.530e+00,1.767e+00,1.936e+00,2.228e+00,2.475e+00,2.563e+00,2.370e+00,2.241e+00,1.746e+00,1.585e+00,1.409e+00,1.184e+00,9.229e-01,7.844e-01,5.954e-01,4.235e-01};
  double stat_Voyager1_Helium_Flux[16]={7.964e-02,9.149e-02,1.002e-01,1.233e-01,1.358e-01,1.336e-01,1.226e-01,1.226e-01,8.831e-02,8.059e-02,7.103e-02,5.973e-02,4.645e-02,3.960e-02,3.001e-02,2.131e-02};
  double syst_Voyager1_Helium_Flux[16]={7.650e-02,8.835e-02,9.680e-02,1.114e-01,1.238e-01,1.282e-01,1.185e-01,1.120e-01,8.730e-02,7.925e-02,7.045e-02,5.920e-02,4.614e-02,3.922e-02,2.977e-02,2.118e-02};

  for(int ee=0;ee<16;ee++){ // unit conversion
    Voyager1_Helium_Ekn[ee] *= 1.e-3;
    Voyager1_Helium_Flux[ee] *= 1.e+3;
    stat_Voyager1_Helium_Flux[ee] *= 1.e+3;
    syst_Voyager1_Helium_Flux[ee] *= 1.e+3;
  }
  
  double eVoyager1_Helium_Ekn[16];
  double eVoyager1_Helium_Flux[16];
  for(int ee=0;ee<16;ee++){
    eVoyager1_Helium_Ekn[ee] = 0.;
    eVoyager1_Helium_Flux[ee] = TMath::Sqrt(stat_Voyager1_Helium_Flux[ee]*stat_Voyager1_Helium_Flux[ee] + syst_Voyager1_Helium_Flux[ee]*syst_Voyager1_Helium_Flux[ee]);
  }
  
  grVOYAGER1_HeliumFluxVSEkn = new TGraphErrors(16, Voyager1_Helium_Ekn, Voyager1_Helium_Flux, eVoyager1_Helium_Ekn, eVoyager1_Helium_Flux);
  grVOYAGER1_HeliumFluxVSEkn->SetName("grVOYAGER1_HeliumFluxVSEkn");
  grVOYAGER1_HeliumFluxVSEkn->SetMarkerStyle(21);
  grVOYAGER1_HeliumFluxVSEkn->SetMarkerSize(1.1);
  grVOYAGER1_HeliumFluxVSEkn->SetMarkerColor(kBlack);
  grVOYAGER1_HeliumFluxVSEkn->SetLineColor(kBlack);
  
}




void NTAMS02Data::ScaleGraphErrors(TGraphErrors* gr, double ScaleValue){
  int NP= gr->GetN();
  for(int ii=0;ii<NP;ii++){
    double Ekn, Val, ErrX, ErrY;
    gr->GetPoint(ii, Ekn, Val);
    ErrX=gr->GetErrorX(ii);
    ErrY=gr->GetErrorY(ii);
    Val  *= ScaleValue;
    ErrY *= ScaleValue;
    gr->SetPoint(ii, Ekn, Val);  
    gr->SetPointError(ii, ErrX, ErrY);  
  }
}



void NTAMS02Data::CleanXErrors(TGraphErrors* gr){
  int NP=gr->GetN();
  for(int ii=0;ii<NP;ii++){
    double e_ekn=0.;
    double e_val=gr->GetErrorY(ii);
    gr->SetPointError(ii, e_ekn,e_val); // SOMETHING WRONG HERE...?
  }
}


/*
void NTAMS02Data::RemoveZero(TGraphErrors* gr){  
  for(int ii=gr->GetN()-1;ii>=0;ii--){
    if( !(gr->GetY()[ii] > 0) ){
      gr->RemovePoint(ii);
    }
  }
}
*/
