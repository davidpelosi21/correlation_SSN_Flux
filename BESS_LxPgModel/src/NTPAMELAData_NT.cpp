#include "NTPAMELAData.h"

//TString PATH_02 = "/Users/davidpelosi/Desktop/TESI-Master-current/Fitting_LxPg_2";
TString PATH_PAMELA = "./";

NTPAMELAData::NTPAMELAData(){
}


void NTPAMELAData::SetAllData(){
  SetProtonData();

}


void NTPAMELAData::SetProtonData(){
  // ---- pick up data files ----
  // TFile* inFileAMS02= new TFile("$NTBASEDIR/THEORY/00matisse/ANALYSIS/ExperimentalData/grAMS02_ProtonHelium_Fluxes_VS_Ekn_2017_RangePaper.root","READ");
  TFile* inFilePAMELA= new TFile(PATH_PAMELA+"/data/grPAMELA_Proton_VS_Ekn_VS_Time.root","READ");
  
  inFilePAMELA->cd();
inFilePAMELA->ls();
 
  // --- get fluxes VS ekn ---
  for(int tt=0;tt<nTimePAMELA_Proton;tt++){
  
      grPAMELA_ProtonFluxVSEkn[tt]=(TGraphErrors*)inFilePAMELA->Get(Form("grPAMELA_ProtonFluxVSEkn_T%d",tt));
      grPAMELA_ProtonFluxVSEkn[tt]->SetName(Form("grPAMELA_ProtonFluxVSEkn_T%d",tt));
      grPAMELA_ProtonFluxVSEkn[tt]->SetMarkerColor(kBlack);
      grPAMELA_ProtonFluxVSEkn[tt]->SetLineColor(kBlack);
      grPAMELA_ProtonFluxVSEkn[tt]->SetMarkerStyle(20);
      grPAMELA_ProtonFluxVSEkn[tt]->SetMarkerSize(0.9);
      grPAMELA_ProtonFluxVSEkn[tt]->SetLineWidth(1);
      //ScaleGraphErrors(grAMS02_ProtonFluxVSEkn[tt], 1.0); 

    // cout<<"Proton Spectrum N. "<<tt<<" ... done!"<<endl;
  }

  cout<<"PROTON EKN.. SPECTRA DONE!"<<endl;
  
  // ---- get fluxes vs time ---- | NT BUT they are RIGIDITY fluxes!!
  for(int ee=0;ee<nEknPAMELA_Proton;ee++){
    grPAMELA_ProtonFluxVSTime[ee]=(TGraphErrors*)inFilePAMELA->Get(Form("grPamelaProtonVSTime_E%d",ee+1));
    
    /*grPAMELA_ProtonFluxVSTime[ee]->SetName(Form("grPAMELA_ProtonFluxVSTime_E%d",ee));
    grPAMELA_ProtonFluxVSTime[ee]->SetMarkerColor(kBlack);
    grPAMELA_ProtonFluxVSTime[ee]->SetLineColor(kBlack);
    grPAMELA_ProtonFluxVSTime[ee]->SetMarkerStyle(20);
    grPAMELA_ProtonFluxVSTime[ee]->SetMarkerSize(0.9);
    grPAMELA_ProtonFluxVSTime[ee]->SetLineWidth(2);
    */
    //ScaleGraphErrors(grAMS02_ProtonFluxVSTime[ee], 1.0);

    // cout<<"Proton profile N. "<<ee<<" ... done!"<<endl;
  }
  cout<<" PROTON TIME-PROFILES DONE!"<<endl;

  
  inFilePAMELA->Close();
  delete inFilePAMELA;

  // ---- get time array ---- 
  double* xTimePAMELA= (double*)grPAMELA_ProtonFluxVSTime[0]->GetX();
  for(int tt=0;tt<nTimePAMELA_Proton;tt++)xTimePAMELA_ProtonFlux[tt]=xTimePAMELA[tt];
  eTimePAMELA_ProtonFlux= 0.5*(xTimePAMELA_ProtonFlux[1]-xTimePAMELA_ProtonFlux[0]);
  
  // ---- get ekn array ----

  double* xEknPAMELA = (double*)grPAMELA_ProtonFluxVSEkn[0]->GetX();
  for(int ee=0;ee<nEknPAMELA_Proton;ee++) xEknPAMELA_ProtonFlux[ee]=xEknPAMELA[ee];

cout<<" import pamela completed"<<endl;
}



void NTPAMELAData::ScaleGraphErrors(TGraphErrors* gr, double ScaleValue){
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



void NTPAMELAData::CleanXErrors(TGraphErrors* gr){
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
