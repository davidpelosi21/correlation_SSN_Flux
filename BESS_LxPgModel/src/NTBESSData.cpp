#include "NTBESSData.h"

//TString PATH_02 = "/Users/davidpelosi/Desktop/TESI-Master-current/Fitting_LxPg_2";
TString PATH_bess = "./";

NTBESSData::NTBESSData(){
}


void NTBESSData::SetAllData(){
  SetProtonData();

}




void NTBESSData::SetProtonData(){

  // ---- pick up data files ----
TFile* inFileBESS= new TFile(PATH_bess+"/data/grBESS_ProtonFluxes.root","READ");

TFile* inFileBESS_extend= new TFile(PATH_bess+"/data/crdb_extracted.root","READ");

inFileBESS->cd();
inFileBESS->ls();

inFileBESS_extend->cd();
inFileBESS_extend->ls();

 
// --- get fluxes VS ekn ---
  for(int tt=0;tt<nTimeBESS_Proton;tt++){
      grBESS_ProtonFluxVSEkn[tt]= new TGraphErrors(); // empty

  /*  
    if(tt==0){
      grBESS_ProtonFluxVSEkn[tt]=(TGraphErrors*)inFileBESS_extend->Get(Form("gr_exp3_errtot",tt));
      grBESS_ProtonFluxVSEkn[tt]->SetName(Form("grBESS_ProtonFluxVSEkn_T%d",tt));
      grBESS_ProtonFluxVSEkn[tt]->SetMarkerColor(kBlack);
      grBESS_ProtonFluxVSEkn[tt]->SetLineColor(kBlack);
      grBESS_ProtonFluxVSEkn[tt]->SetMarkerStyle(20);
      grBESS_ProtonFluxVSEkn[tt]->SetMarkerSize(0.9);
      grBESS_ProtonFluxVSEkn[tt]->SetLineWidth(1);
}


    if(tt==1){
      grBESS_ProtonFluxVSEkn[tt]=(TGraphErrors*)inFileBESS_extend->Get(Form("gr_exp4_errtot",tt));
      grBESS_ProtonFluxVSEkn[tt]->SetName(Form("grBESS_ProtonFluxVSEkn_T%d",tt));
      grBESS_ProtonFluxVSEkn[tt]->SetMarkerColor(kBlack);
      grBESS_ProtonFluxVSEkn[tt]->SetLineColor(kBlack);
      grBESS_ProtonFluxVSEkn[tt]->SetMarkerStyle(20);
      grBESS_ProtonFluxVSEkn[tt]->SetMarkerSize(0.9);
      grBESS_ProtonFluxVSEkn[tt]->SetLineWidth(1);
}

*/

    if(tt==0){
      grBESS_ProtonFluxVSEkn[tt]=(TGraphErrors*)inFileBESS_extend->Get(Form("gr_exp5_errtot",tt));
      grBESS_ProtonFluxVSEkn[tt]->SetName(Form("grBESS_ProtonFluxVSEkn_T%d",tt));
      grBESS_ProtonFluxVSEkn[tt]->SetMarkerColor(kBlack);
      grBESS_ProtonFluxVSEkn[tt]->SetLineColor(kBlack);
      grBESS_ProtonFluxVSEkn[tt]->SetMarkerStyle(20);
      grBESS_ProtonFluxVSEkn[tt]->SetMarkerSize(0.9);
      grBESS_ProtonFluxVSEkn[tt]->SetLineWidth(1);
}


    if(tt==1){
      grBESS_ProtonFluxVSEkn[tt]=(TGraphErrors*)inFileBESS_extend->Get(Form("gr_exp6_errtot",tt));
      grBESS_ProtonFluxVSEkn[tt]->SetName(Form("grBESS_ProtonFluxVSEkn_T%d",tt));
      grBESS_ProtonFluxVSEkn[tt]->SetMarkerColor(kBlack);
      grBESS_ProtonFluxVSEkn[tt]->SetLineColor(kBlack);
      grBESS_ProtonFluxVSEkn[tt]->SetMarkerStyle(20);
      grBESS_ProtonFluxVSEkn[tt]->SetMarkerSize(0.9);
      grBESS_ProtonFluxVSEkn[tt]->SetLineWidth(1);
}



    if(tt==2){
      grBESS_ProtonFluxVSEkn[tt]=(TGraphErrors*)inFileBESS_extend->Get(Form("gr_exp7_errtot",tt));
      grBESS_ProtonFluxVSEkn[tt]->SetName(Form("grBESS_ProtonFluxVSEkn_T%d",tt));
      grBESS_ProtonFluxVSEkn[tt]->SetMarkerColor(kBlack);
      grBESS_ProtonFluxVSEkn[tt]->SetLineColor(kBlack);
      grBESS_ProtonFluxVSEkn[tt]->SetMarkerStyle(20);
      grBESS_ProtonFluxVSEkn[tt]->SetMarkerSize(0.9);
      grBESS_ProtonFluxVSEkn[tt]->SetLineWidth(1);
}

    if(tt==3){
      grBESS_ProtonFluxVSEkn[tt]=(TGraphErrors*)inFileBESS_extend->Get(Form("gr_exp8_errtot",tt));
      grBESS_ProtonFluxVSEkn[tt]->SetName(Form("grBESS_ProtonFluxVSEkn_T%d",tt));
      grBESS_ProtonFluxVSEkn[tt]->SetMarkerColor(kBlack);
      grBESS_ProtonFluxVSEkn[tt]->SetLineColor(kBlack);
      grBESS_ProtonFluxVSEkn[tt]->SetMarkerStyle(20);
      grBESS_ProtonFluxVSEkn[tt]->SetMarkerSize(0.9);
      grBESS_ProtonFluxVSEkn[tt]->SetLineWidth(1);
}


    if(tt==4){
      grBESS_ProtonFluxVSEkn[tt]=(TGraphErrors*)inFileBESS_extend->Get(Form("gr_exp2_errtot",tt));
      grBESS_ProtonFluxVSEkn[tt]->SetName(Form("grBESS_ProtonFluxVSEkn_T%d",tt));
      grBESS_ProtonFluxVSEkn[tt]->SetMarkerColor(kBlack);
      grBESS_ProtonFluxVSEkn[tt]->SetLineColor(kBlack);
      grBESS_ProtonFluxVSEkn[tt]->SetMarkerStyle(20);
      grBESS_ProtonFluxVSEkn[tt]->SetMarkerSize(0.9);
      grBESS_ProtonFluxVSEkn[tt]->SetLineWidth(1);
}



    if(tt==5){
      grBESS_ProtonFluxVSEkn[tt]=(TGraphErrors*)inFileBESS_extend->Get(Form("gr_exp1_errtot",tt));
      grBESS_ProtonFluxVSEkn[tt]->SetName(Form("grBESS_ProtonFluxVSEkn_T%d",tt));
      grBESS_ProtonFluxVSEkn[tt]->SetMarkerColor(kBlack);
      grBESS_ProtonFluxVSEkn[tt]->SetLineColor(kBlack);
      grBESS_ProtonFluxVSEkn[tt]->SetMarkerStyle(20);
      grBESS_ProtonFluxVSEkn[tt]->SetMarkerSize(0.9);
      grBESS_ProtonFluxVSEkn[tt]->SetLineWidth(1);
}



     if( tt == 6){
      grBESS_ProtonFluxVSEkn[tt]=(TGraphErrors*)inFileBESS->Get("grBESSPolarI_ProtonFlux_2016");
      grBESS_ProtonFluxVSEkn[tt]->SetName(Form("grBESS_ProtonFluxVSEkn_T%d",tt));
      grBESS_ProtonFluxVSEkn[tt]->SetMarkerColor(kBlack);
      grBESS_ProtonFluxVSEkn[tt]->SetLineColor(kBlack);
      grBESS_ProtonFluxVSEkn[tt]->SetMarkerStyle(20);
      grBESS_ProtonFluxVSEkn[tt]->SetMarkerSize(0.9);
      grBESS_ProtonFluxVSEkn[tt]->SetLineWidth(1);
      //ScaleGraphErrors(grBESS_ProtonFluxVSEkn[tt], 1.0); 
  }



     if( tt == 7){
      grBESS_ProtonFluxVSEkn[tt]=(TGraphErrors*)inFileBESS->Get("grBESSPolarII_ProtonFlux_2016");
      grBESS_ProtonFluxVSEkn[tt]->SetName(Form("grBESS_ProtonFluxVSEkn_T%d",tt));
      grBESS_ProtonFluxVSEkn[tt]->SetMarkerColor(kBlack);
      grBESS_ProtonFluxVSEkn[tt]->SetLineColor(kBlack);
      grBESS_ProtonFluxVSEkn[tt]->SetMarkerStyle(20);
      grBESS_ProtonFluxVSEkn[tt]->SetMarkerSize(0.9);
      grBESS_ProtonFluxVSEkn[tt]->SetLineWidth(1);
      //ScaleGraphErrors(grBESS_ProtonFluxVSEkn[tt], 1.0); 
  }


  }


  cout<<"PROTON EKN.. SPECTRA DONE!"<<endl;
  // ---- get fluxes vs time 


    
   grBESS_ProtonFluxVSTime->SetName("grBESS_ProtonFluxVSTime");
    grBESS_ProtonFluxVSTime->SetMarkerColor(kBlack);
    grBESS_ProtonFluxVSTime->SetLineColor(kBlack);
    grBESS_ProtonFluxVSTime->SetMarkerStyle(20);
    grBESS_ProtonFluxVSTime->SetMarkerSize(0.9);
    grBESS_ProtonFluxVSTime->SetLineWidth(2);
 
  
  
  //TDatime Date0( 1995, 7, 01, 0, 0, 0);
  //TDatime Date1( 1995, 7, 01, 0, 0, 0);

  //ok 
  TDatime Date2( 1995, 7, 01, 0, 0, 0);
  TDatime Date3( 1997, 7, 01, 0, 0, 0);
  TDatime Date4( 1998, 7, 01, 0, 0, 0);
  TDatime Date5( 1999, 8, 01, 0, 0, 0);
  TDatime Date6( 2000, 8, 01, 0, 0, 0);
  TDatime Date7( 2002, 8, 01, 0, 0, 0);
  TDatime Date8( 2004, 12, 01, 0, 0, 0);
  TDatime Date9( 2008, 1, 01, 0, 0, 0);


  //double t0= (double)Date2.Convert() - 2*31622400; //unsigned int 10 cifre
  //double t1= (double)Date2.Convert() - 31622400;  //unsigned int 10 cifre


  double t2= (double)Date2.Convert();  //unsigned int 10 cifre
  double t3= (double)Date3.Convert();  //unsigned int 10 cifre
  double t4= (double)Date4.Convert();  //unsigned int 10 cifre
  double t5= (double)Date5.Convert();  //unsigned int 10 cifre
  double t6= (double)Date6.Convert();  //unsigned int 10 cifre
  double t7= (double)Date7.Convert();  //unsigned int 10 cifre
  double t8= (double)Date8.Convert();  //unsigned int 10 cifre
  double t9= (double)Date9.Convert();  //unsigned int 10 cifre


  //xTimeBESS_ProtonFlux[0]=t0;
  //xTimeBESS_ProtonFlux[1]=t1;
  xTimeBESS_ProtonFlux[0]=t2;
  xTimeBESS_ProtonFlux[1]=t3;
   xTimeBESS_ProtonFlux[2]=t4;
    xTimeBESS_ProtonFlux[3]=t5;
     xTimeBESS_ProtonFlux[4]=t6;
      xTimeBESS_ProtonFlux[5]=t7;
       xTimeBESS_ProtonFlux[6]=t8;
        xTimeBESS_ProtonFlux[7]=t9;



 eTimeBESS_ProtonFlux = ((3.154e+7)/2.0001); // time error: 6 months

  cout<<" PROTON TIME-PROFILES DONE!"<<endl;

  
  inFileBESS->Close();
  delete inFileBESS;


   
  // ---- get ekn array ----

  /*double* xEknBESS= (double*)grBESS_ProtonFluxVSEkn[0]->GetX();

  for(int ee=0;ee<nEknBESS_Proton;ee++)  {

        xEknBESS_ProtonFlux[ee]=xEknBESS[ee];
  }
*/

}









void NTBESSData::ScaleGraphErrors(TGraphErrors* gr, double ScaleValue){
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



void NTBESSData::CleanXErrors(TGraphErrors* gr){
  int NP=gr->GetN();
  for(int ii=0;ii<NP;ii++){
    double e_ekn=0.;
    double e_val=gr->GetErrorY(ii);
    gr->SetPointError(ii, e_ekn,e_val); // SOMETHING WRONG HERE...?
  }
}


/*
void NTBESSData::RemoveZero(TGraphErrors* gr){  
  for(int ii=gr->GetN()-1;ii>=0;ii--){
    if( !(gr->GetY()[ii] > 0) ){
      gr->RemovePoint(ii);
    }
  }
}
*/
