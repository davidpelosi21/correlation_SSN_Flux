#include "TROOT.h"
#include "TMath.h"
#include "TGraphErrors.h"
#include "TGraph.h"
#include "TH2F.h"
#include "TAxis.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TApplication.h"
#include "TMinuit.h"
#include <iostream>
#include <iomanip>
#include <string.h>
#include <tgmath.h>

// *******   MAIN  *********

void cross(){ //TOP

    cout<<"********************"<<endl;
    cout<<"Cross correlation Diffusion Coefficient vs Solar Proxy (A polarity label)"<<endl;
    cout<<"********************"<<endl;
//read K0_vs time from fitting ams data -> time_scale BR
    

  //TFile* gridFile= new TFile(Form("K0_vs_time_Potgeiter3GV.root"),"READ");
    
 TFile* gridFile= new TFile(Form("K0_vs_time_Potgeiter.root"),"READ");

   TGraphErrors* grK0vsTime = (TGraphErrors*)gridFile->Get("grK0VSTime_AMS02");


    gridFile->Close();


    // intervallo temporale start date
    TDatime *t= new TDatime();
    t->Set(grK0vsTime->GetX()[0]);
    int year_start = t->GetYear();
    int month_start = t->GetMonth();
    int day_start = t->GetDay();
    
      

    cout<<"K_{0} data from: "<<day_start<<"/"<<month_start<<"/"<<year_start;
    
    //end date
    t->Set(grK0vsTime->GetX()[grK0vsTime->GetN()-1]);
    int year_end = t->GetYear();
    int month_end = t->GetMonth();
    int day_end = t->GetDay();
    
    cout<<"   to: "<<day_end<<"/"<<month_end<<"/"<<year_end<<endl;
// read solar proxy
    

TFile *f = new TFile("SSN_convert.root");

TGraphErrors *SSN_daily = (TGraphErrors*)f->Get("SSN_convert");
    //SSN_daily->SetDirectory(0);

    
//--------------------  CROSS CORRELATION  ---------------
//colori diversi per polarità
    TGraphErrors* Cross_K0_Proxy =  new TGraphErrors(); //general A    

    int ii=0;

//lag  11.4 mesi
    int nSeconds = (int)(11.4 * 30.44 * 24 * 60 * 60);

    for (int i = 0; i<grK0vsTime->GetN(); i++) {

                    TTimeStamp timestamp(grK0vsTime->GetX()[i]);
                    timestamp.Add(-nSeconds);
                    TDatime datetime2(timestamp.GetSec());
        
                    
                    //Cross_K0_Proxy->SetPoint(ii, SSN_daily->Eval(grK0vsTime->GetX()[i]),grK0vsTime->GetY()[i]);
                    Cross_K0_Proxy->SetPoint(ii, SSN_daily->Eval(datetime2.Convert()),grK0vsTime->GetY()[i]);
                    Cross_K0_Proxy->SetPointError(ii,0,grK0vsTime->GetEY()[i]);

                    ii++;
         }

    TCanvas* c1= new TCanvas("c1", "Cross correlation K_{0} vs solar proxy");
    
    c1->cd();
        
    Cross_K0_Proxy->SetMarkerStyle(21);

    //setmarkersize
    Cross_K0_Proxy->SetMarkerSize(0.5);
//setmarker color
    Cross_K0_Proxy->SetMarkerColor(kBlack);
//setline color
    Cross_K0_Proxy->SetLineColor(kBlack);

Cross_K0_Proxy->GetXaxis()->SetTitle("SSN");
Cross_K0_Proxy->GetYaxis()->SetTitle("K_{0} [10^{22} cm^{2} s^{-1}]");
//Cross_K0_Proxy->GetYaxis()->SetTitle("a(t)");

   Cross_K0_Proxy->Draw("ap");
  


//save TGraphErrors in a root file
    TFile* f1 = new TFile("Cross_K0_ProxySSN.root","RECREATE");
    Cross_K0_Proxy->SetName("Cross_K0_ProxySSN");
    Cross_K0_Proxy->Write();
    f1->Close();

    cout<<"********************"<<endl;
    cout<<"Cross correlation K_{0} vs solar proxy"<<endl;
    cout<<"********************"<<endl;
    cout<<"Cross correlation saved in Cross_K0_Proxy.root"<<endl;
    cout<<"********************"<<endl;
    
    return;

  return 1;
    
  // YYYYYYYYY end  YYYYYYYYYYY
}



