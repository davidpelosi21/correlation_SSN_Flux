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

void cross_Polarity(){ //TOP

    cout<<"********************"<<endl;
    cout<<"Cross correlation Diffusion Coefficient vs Solar Proxy (A polarity label)"<<endl;
    cout<<"********************"<<endl;
//read K0_vs time from fitting ams data -> time_scale BR
    

  //TFile* gridFile= new TFile(Form("K0_vs_time_Potgeiter3GV.root"),"READ");
    
 TFile* gridFile= new TFile(Form("K0_vs_time_Potgeiter.root"),"READ");
    
   //TGraphErrors* grK0vsTime_PAMELA = (TGraphErrors*)gridFile->Get("grK0VSTime_PAMELA");
   TGraphErrors* grK0vsTime_PAMELA = (TGraphErrors*)gridFile->Get("grK0VSTime_PAMELA");
    //grK0vsTime_PAMELA->SetDirectory(0);

    gridFile->Close();


    // intervallo temporale start date
    TDatime *t= new TDatime();
    t->Set(grK0vsTime_PAMELA->GetX()[0]);
    int year_start = t->GetYear();
    int month_start = t->GetMonth();
    int day_start = t->GetDay();
    
      

    cout<<"K_{0} data from: "<<day_start<<"/"<<month_start<<"/"<<year_start;
    
    //end date
    t->Set(grK0vsTime_PAMELA->GetX()[grK0vsTime_PAMELA->GetN()-1]);
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
   
    TGraphErrors* Cross_K0_Proxy_A_pos =  new TGraphErrors();
    TGraphErrors* Cross_K0_Proxy_A_neg =  new TGraphErrors();
    TGraphErrors* Cross_K0_Proxy_A_un =  new TGraphErrors();

    TGraphErrors* Cross_K0_Proxy =  new TGraphErrors(); //general A    
   /* for (int i = 1; i<10; i++) {
        
        cout<<"Proxy x value: "<<std::setprecision(10)<<SSN_daily->GetX()[i]<<endl;
    }*/

    
    TDatime * t1 = new TDatime();
    TDatime * t2 = new TDatime();


    t1->Set(2012, 9, 1,0,0,0);
    t2->Set(2014, 3, 1,0,0,0);

    
    int cpos=0;
    int cneg=0;
    int cundef=0;
    
    int ii,ij,ik;
    ii=0;
    ij=0;
    ik=0;
//lag  11.4 mesi
    int nSeconds = (int)(11.4 * 30.44 * 24 * 60 * 60);

cout<<t1->Convert()<<"  "<<t2->Convert()<<"  |"<<grK0vsTime_PAMELA->GetX()[4]<<endl;

    for (int i = 0; i<grK0vsTime_PAMELA->GetN(); i++) {
       // cout<<"K0 x value: "<<std::setprecision(10)<<grK0vsTime_PAMELA->GetX()[i]<<"  --- "<<decimalyear<<"    "<<day<<"/"<<month<<"/"<<year<<endl;
               
                    TTimeStamp timestamp(grK0vsTime_PAMELA->GetX()[i]);
                    timestamp.Add(-nSeconds);
                    TDatime datetime2(timestamp.GetSec());
            
      
                 if (grK0vsTime_PAMELA->GetX()[i] <= t1->Convert()) {
                    //cout<<"   A<0"<<endl;
                    
                    Cross_K0_Proxy_A_neg->SetPoint(ii, SSN_daily->Eval(grK0vsTime_PAMELA->GetX()[i]),grK0vsTime_PAMELA->GetY()[i]);
                    //Cross_K0_Proxy_A_neg->SetPoint(ii, SSN_daily->Eval(datetime2.Convert()),grK0vsTime_PAMELA->GetY()[i]);
                    Cross_K0_Proxy_A_neg->SetPointError(ii,0,grK0vsTime_PAMELA->GetEY()[i]);
                    cneg++; ii++;
                }
                else if (grK0vsTime_PAMELA->GetX()[i] >= t2->Convert()) {
                    //cout<<"   A>0"<<endl;

                    Cross_K0_Proxy_A_pos->SetPoint(ij, SSN_daily->Eval(grK0vsTime_PAMELA->GetX()[i]),grK0vsTime_PAMELA->GetY()[i]);
                    //Cross_K0_Proxy_A_pos->SetPoint(ij, SSN_daily->Eval(datetime2.Convert()),grK0vsTime_PAMELA->GetY()[i]);
                    Cross_K0_Proxy_A_pos->SetPointError(ij,0,grK0vsTime_PAMELA->GetEY()[i]);
                    cpos++; ij++;
                }
                
                else {
                    //cout<<"  A undefined"<<endl;
                
                    Cross_K0_Proxy_A_un->SetPoint(ik,SSN_daily->Eval(grK0vsTime_PAMELA->GetX()[i]),grK0vsTime_PAMELA->GetY()[i]);
                    //Cross_K0_Proxy_A_un->SetPoint(ik,SSN_daily->Eval(datetime2.Convert()),grK0vsTime_PAMELA->GetY()[i]);
                    Cross_K0_Proxy_A_un->SetPointError(ik,0,grK0vsTime_PAMELA->GetEY()[i]);
                    cundef++; ik++;
                }  
                

            //Cross_K0_Proxy->SetPoint(i,SSN_daily->Eval(grK0vsTime_PAMELA->GetX()[i]),grK0vsTime_PAMELA->GetY()[i]);
            //Cross_K0_Proxy->SetPointError(i,0,grK0vsTime_PAMELA->GetEY()[i]);
            

        }

    
    
cout<<"Cross correlation data "<<cpos+cneg+cundef<<"  == K0 data  "<<grK0vsTime_PAMELA->GetN()<<endl;

    
//cout<<"Date reverse  "<<revert(2017.001)<<endl;
    TCanvas* c1= new TCanvas("c1", "Cross correlation K_{0} vs solar proxy");
    
    c1->cd();
        
    Cross_K0_Proxy_A_un->SetMarkerStyle(21);
    Cross_K0_Proxy_A_pos->SetMarkerStyle(21);
    Cross_K0_Proxy_A_neg->SetMarkerStyle(21);
    Cross_K0_Proxy->SetMarkerStyle(21);

    //setmarkersize
    Cross_K0_Proxy_A_un->SetMarkerSize(0.5);
    Cross_K0_Proxy_A_pos->SetMarkerSize(0.5);
    Cross_K0_Proxy_A_neg->SetMarkerSize(0.5);
    Cross_K0_Proxy->SetMarkerSize(0.5);

//setmarker color
    Cross_K0_Proxy_A_un->SetMarkerColor(kGreen);
    Cross_K0_Proxy_A_pos->SetMarkerColor(kRed);
    Cross_K0_Proxy_A_neg->SetMarkerColor(kBlue);
    Cross_K0_Proxy->SetMarkerColor(kBlack);
//setline color
    Cross_K0_Proxy_A_un->SetLineColor(kGreen);
    Cross_K0_Proxy_A_pos->SetLineColor(kRed);
    Cross_K0_Proxy_A_neg->SetLineColor(kBlue);
    Cross_K0_Proxy->SetLineColor(kBlack);


    Cross_K0_Proxy_A_un->SetTitle("A undefined");
    Cross_K0_Proxy_A_pos->SetTitle("A >0 ");
    Cross_K0_Proxy_A_neg->SetTitle("A < 0 ");

    
    
    TMultiGraph *mg = new TMultiGraph();
    mg->Add(Cross_K0_Proxy_A_un,"p");
    mg->Add(Cross_K0_Proxy_A_pos,"p");
    mg->Add(Cross_K0_Proxy_A_neg,"p");
    
    mg->SetTitle("Cross correlation K_{0} vs solar proxy; SSN(t) ;K_{0} [10^{22} cm^{2} s^{-1}]");
    
    //mg->GetXaxis()->SetTitle(" SSN(t)");
    //mg->GetYaxis()->SetTitle("K_{0} [10^{22} cm^{2} s^{-1}]");
    //Cross_K0_Proxy_A_neg->Draw("ap");
    mg->Draw("ap");
   //Cross_K0_Proxy->Draw("ap");

    //grK0vsTime_PAMELA->Draw("ap");
    c1->BuildLegend();
                                

//save TGraphErrors in a root file
    TFile* f1 = new TFile("Cross_K0_ProxySSN_Polarity.root","RECREATE");
    Cross_K0_Proxy_A_un->SetName("Cross_K0_ProxySSN_A_un");
    Cross_K0_Proxy_A_un->Write();
    Cross_K0_Proxy_A_pos->SetName("Cross_K0_ProxySSN_A_pos");
    Cross_K0_Proxy_A_pos->Write();
    Cross_K0_Proxy_A_neg->SetName("Cross_K0_ProxySSN_A_neg");
    Cross_K0_Proxy_A_neg->Write();

    f1->Close();

    cout<<"********************"<<endl;
    cout<<"Cross correlation K_{0} vs solar proxy"<<endl;
    cout<<"********************"<<endl;
    cout<<"Cross correlation saved in Cross_K0_Proxy_Polarity.root"<<endl;
    cout<<"********************"<<endl;


  return 1;
    
  // YYYYYYYYY end  YYYYYYYYYYY
}



