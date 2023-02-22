//  open a root file and extract th2f, for each x bin create the projection on y and save it in a tgraph

//import root libraries
#include "../src/Conversion.h"
#include "TFile.h"
#include "TH2F.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TStyle.h"
//include  TDatime
#include "TDatime.h"
#include "TLegend.h"
#include "TMath.h"
#include "TString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TColor.h"
//stdio libraries
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

//namespace std
//using namespace std;

void ams() {

TFile *f = new TFile("ProtonMonthlyRebin.root");

TH2F *h= (TH2F*)f->Get("MonthlyFlux_Rebin");
TH2F *h2= (TH2F*)f->Get("Monthly_Stat_errors_Rebin");
TH2F *h3= (TH2F*)f->Get("Monthly_TDep_Syst_errors_Rebin");
TH2F *h4= (TH2F*)f->Get("Monthly_TIndep_Syst_errors_Rebin");

//for each time bin

TH1F *flux[h->GetNbinsX()];
TH1F *rel_stat_err[h2->GetNbinsX()];
TH1F *rel_sys_Dep_err[h3->GetNbinsX()];
TH1F *rel_sys_InDep_err[h4->GetNbinsX()];

for (int i=0; i<h->GetNbinsX(); i++) {
  flux[i] = (TH1F*)h->ProjectionY(Form("g_N%d",i),i+1,i+1);
}

for (int i=0; i<h2->GetNbinsX(); i++) {
  rel_stat_err[i] = (TH1F*)h2->ProjectionY(Form("rel_stat_err_N%d",i),i+1,i+1);
}

for (int i=0; i<h3->GetNbinsX(); i++) {
  rel_sys_Dep_err[i] = (TH1F*)h3->ProjectionY(Form("rel_sys_Dep_err_N%d",i),i+1,i+1);
}

for(int i=0; i<h4->GetNbinsX(); i++) {
   rel_sys_InDep_err[i] = (TH1F*)h4->ProjectionY(Form("rel_sys_InDep_err_N%d",i),i+1,i+1);
}

//for each rigidty bin 

TH1F *flux_R[h->GetNbinsY()];
TH1F *rel_stat_err_R[h2->GetNbinsY()];
TH1F *rel_sys_Dep_err_R[h3->GetNbinsY()];
TH1F *rel_sys_InDep_err_R[h4->GetNbinsY()];

for (int i=0; i<h->GetNbinsY(); i++) {
  flux_R[i] = (TH1F*)h->ProjectionX(Form("R_g_N%d",i),i+1,i+1);
}

for (int i=0; i<h2->GetNbinsY(); i++) {
  rel_stat_err_R[i] = (TH1F*)h2->ProjectionX(Form("R_rel_sys_err_N%d",i),i+1,i+1);
}

for (int i=0; i<h3->GetNbinsY(); i++) {
  rel_sys_Dep_err_R[i] = (TH1F*)h3->ProjectionX(Form("R_rel_sys_Dep_err_N%d",i),i+1,i+1);
}

for (int i=0; i<h4->GetNbinsY(); i++) {
  rel_sys_InDep_err_R[i] = (TH1F*)h4->ProjectionX(Form("R_rel_sys_InDep__err_N%d",i),i+1,i+1);
}


TGraphErrors *flux_err_stat = new TGraphErrors[h->GetNbinsX()];
TGraphErrors *flux_err_stat_rig = new TGraphErrors[h->GetNbinsX()];

TH1F *time = new TH1F();
double date_root;
TDatime * t2 = new TDatime();
time = (TH1F*)h->ProjectionX();

for (int j=0; j<h->GetNbinsX(); j++) {
   date_root = time->GetBinCenter(j+1);
   t2->Set(date_root);

   flux_err_stat[j].SetName(Form("grProton_vs_Ekn_T%d",j));
   flux_err_stat_rig[j].SetName(Form("grProton_vs_Rig_T%d",j));
   


   flux_err_stat[j].SetTitle(Form("%d/%d/%d",t2->GetDay(),t2->GetMonth(),t2->GetYear()));
   //flux_err_stat[j].SetTitle(Form("%f",date_root));


   for (int i=0; i<h->GetNbinsY(); i++)   {
      //ekn
      flux_err_stat[j].SetPoint(i,RigToEkn(sqrt(flux[j]->GetBinLowEdge(i+2)  * flux[j]->GetBinLowEdge(i+2+1))), flux[j]->GetBinContent(i+2)*dRdE(RigToEkn(flux[j]->GetBinCenter(i+2))) );
      flux_err_stat[j].SetPointError(i,0,RigToEkn( flux[j]->GetBinContent(i+2)*   sqrt(  pow(rel_stat_err[j]->GetBinContent(i+2),2)  +  pow(rel_sys_Dep_err[j]->GetBinContent(i+2),2) + pow(rel_sys_InDep_err[j]->GetBinContent(i+2),2)   )        ) );

   //rig
   flux_err_stat_rig[j].SetPoint(i,sqrt(flux[j]->GetBinLowEdge(i+2)  * flux[j]->GetBinLowEdge(i+2+1)),flux[j]->GetBinContent(i+2));
   flux_err_stat_rig[j].SetPointError(i,0, flux[j]->GetBinContent(i+2)*   sqrt(  pow(rel_stat_err[j]->GetBinContent(i+2),2)  +  pow(rel_sys_Dep_err[j]->GetBinContent(i+2),2) + pow(rel_sys_InDep_err[j]->GetBinContent(i+2),2)          ) );

   }
}


//time graph
TGraphErrors *flux_err_stat_R = new TGraphErrors[h->GetNbinsY()];

for (int j=0; j<h->GetNbinsY(); j++) {
   flux_err_stat_R[j].SetName(Form("grProton_vs_Time_R%d",j+1));
   

   for (int i=0; i<114; i++) {

     flux_err_stat_R[j].SetPoint(i,sqrt(flux_R[j]->GetBinLowEdge(i+2)*flux_R[j]->GetBinLowEdge(i+2+1)), flux_R[j]->GetBinContent(i+2) );
     
      flux_err_stat_R[j].SetPointError( i, 0, flux_R[j]->GetBinContent(i+2)* sqrt(pow(rel_stat_err_R[j]->GetBinContent(i+2),2)  +  pow(rel_sys_Dep_err_R[j]->GetBinContent(i+2),2) + pow(rel_sys_InDep_err_R[j]->GetBinContent(i+2),2)  ) );
   }
}


TFile *f1 = new TFile("AMS02.root","RECREATE");

f1->cd();
for (int i=0; i<h->GetNbinsX(); i++) {
  //flux[i]->Write();
  //stat_err[i]->Write();
  //if (flux[i]->GetEntries() != 0 &&  stat_err[i]->GetEntries() != 0) {
   if (flux[i]->GetEntries() != 0) {
      flux_err_stat[i].SetMarkerColor(4);
      flux_err_stat[i].SetMarkerStyle(21);

      flux_err_stat_rig[i].SetMarkerColor(4);
      flux_err_stat_rig[i].SetMarkerStyle(21);


      flux_err_stat[i].Write();
      flux_err_stat_rig[i].Write();


  }  
}

for (int i=0; i<h->GetNbinsY(); i++) {
  //flux[i]->Write();
  //stat_err[i]->Write();
  //if (flux[i]->GetEntries() != 0 &&  stat_err[i]->GetEntries() != 0) {
   if (flux_R[i]->GetEntries() != 0) {
      flux_err_stat_R[i].SetMarkerColor(4);
      flux_err_stat_R[i].SetMarkerStyle(21);
      flux_err_stat_R[i].Write();


  }  
}


//log scale on y axis of a tgraph
//gPad->SetLogy();
flux_err_stat->Write();

f1->Close();
f->Close();


cout << "End" << endl;

}







