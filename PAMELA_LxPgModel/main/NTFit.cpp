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
#include <chrono>
//#include "NTCosmicRayData.cpp"
#include "NTAMS02Data.h"
#include "NTPAMELAData.h"
#include "Constants.h"
#include "SolarModulation.h"


TString PATHm = "./";

// --------- NTFit -----------
// -- Make K0 e N0 -Fits to CR proton data: -AMS02
// ---------------------------------

//int Kmodel =  1; //standard model
//int Kmodel =  2; //power law model
int Kmodel =  3; // Potgeiter  model
//int Kmodel =  4; // Potgeiter Model with Rk free

double Z = 1;
double A = 1;
int kLISModel = 4;  //serve solo per scegliere i th2f per i vari LIS
int iXiD = 49; // zero drift [49]
//scelgo indice  int iK0  = 35; // avg diffusion level [35]

double K0=1.0;
double a = 1.0;
double b = 1.0;
double rk = 1793.96; // MV; //3 GV
int iK0 = 35; // 1AU should be = mult = 5

int  run=0;

bool kReNorm = true;
// ---- constants ----
const int NITERATIONS = 10000;

double EMIN= 0.02; //0.07
double EMAX= 15.;
int thisRIG1 = 0; // rig bin for p/He plots Low-R
int thisRIG2 = 3; // rig bin for p/He plots High-R

static const int nTimePAMELA_Proton      = 83; //81
const int iExp_PAMELA       = 6;

double bestChi2_PAMELA[nTimePAMELA_Proton];

int NDF_PAMELA[nTimePAMELA_Proton];

// --- graph style ---
void SetStyleHistoVSEkn(TH2F* hh);
void SetStyleHistoVSTime(TH2F* hh);

void SetTimeGraphStyle(TGraphErrors* gr);


void SetGraphK0Style(); // init K0vsTime graphs
void SetGraphAStyle(); // init K0vsTime graphs
void SetGraphBStyle(); // init K0vsTime graphs
void SetGraphRkStyle(); // init K0vsTime graphs
void SetGraphChis2Style(); // init K0vsTime graphs
void SetGraphN0Style(); // 
// --- graph manipulation ---
void FlattenGraph(TGraph* gr, double OldFlatIndex, double NewFlatIndex); // for LIS
void ScaleGraph(TGraph* gr, double ScaleFactor);
void MakeDataToModel(TGraphErrors* grData, TGraph* grModel, TGraphErrors* grDataToModel);


// --- minimization ---
double GetChiSquare(SolarModulation* NTSolution, double K0, int iTime, int iExperiment, double Norm,double a,double b, double rk);
void FCN(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);


// ---- external objects ----
NTAMS02Data* amsdata;  //def puntantore alla classe amsdata
NTPAMELAData* pameladata;  //def puntantore alla classe amsdata


// --- model output ---
TH2F* hLn_Flux_vs_KScale_vs_KEnergy_Z1_A1; // PROTON vs Ekn: for FIT

TH2F* hLn_Flux_vs_KScale_vs_Rigidity_Z1_A1; // vs Rig: for PREDICTIONS

TGraphErrors* grK0vsTime_PAMELA;
TGraphErrors* grAvsTime_PAMELA;
TGraphErrors* grBvsTime_PAMELA;
TGraphErrors* grRkvsTime_PAMELA;
TGraphErrors* grChis2vsTime_PAMELA;
TGraphErrors* grN0vsTime_PAMELA;

TGraph* Chi2_K0 = new TGraph();
int c=0;


//def globale della classe SolarModulation
SolarModulation* NTSolution= new SolarModulation(Z, A, kLISModel, iK0,K0,a,b,rk,iXiD, Kmodel);

using namespace std::chrono;
auto start= chrono::high_resolution_clock::now();
auto stop = chrono::high_resolution_clock::now();
auto startM= chrono::high_resolution_clock::now();
auto stopM = chrono::high_resolution_clock::now();


//*****+  MAIN ************
int main(int argc, char **argv){
    //TApplication theApp("App",&argc,argv);
    startM = chrono::high_resolution_clock::now();
    
    NTSolution->InitJLis();
    NTSolution->SetJLis();
    NTSolution->KscaleGrid();
    NTSolution->InitModulation();
    NTSolution->Solve_BC();

    // ---- set AMS02 data ----
    amsdata= new NTAMS02Data();
    pameladata= new NTPAMELAData();

    amsdata->SetProtonData();
    pameladata->SetProtonData();

    // importo histo per algoritmo interpolazione
    TFile* gridFile= new TFile(Form("../LxPgModel1D_set_TH2F/out/Model_DriftOFF/histo_FluxResults_LIS%d_Z1_A1.root",kLISModel),"READ");
    
    hLn_Flux_vs_KScale_vs_KEnergy_Z1_A1= (TH2F*)gridFile->Get("hLn_Flux_vs_KScale_vs_KEnergy_Z1_A1");
    hLn_Flux_vs_KScale_vs_KEnergy_Z1_A1->SetDirectory(0);

    hLn_Flux_vs_KScale_vs_Rigidity_Z1_A1= (TH2F*)gridFile->Get("hLn_Flux_vs_KScale_vs_Rigidity_Z1_A1");
    hLn_Flux_vs_KScale_vs_Rigidity_Z1_A1->SetDirectory(0);

    gridFile->Close();

  
    grK0vsTime_PAMELA= new TGraphErrors();

    SetGraphK0Style();

    grAvsTime_PAMELA= new TGraphErrors();
    
    SetGraphAStyle();

  grBvsTime_PAMELA= new TGraphErrors();
    
    SetGraphBStyle();

   grRkvsTime_PAMELA= new TGraphErrors();
    SetGraphRkStyle();

    grChis2vsTime_PAMELA= new TGraphErrors();
    
    SetGraphChis2Style();
    
    grN0vsTime_PAMELA= new TGraphErrors();
    
    SetGraphN0Style();
    
// ---- initializations ----
  for(int tt=0;tt<nTimePAMELA_Proton;tt++) bestChi2_PAMELA[tt]= 1.e+9;
  // ----- minimization ----
    //const int NPAR = 4;  //K0, iExp, iTime, NORM [ Emin-Emax from ABOVE]

    const int NPAR = 7; // K0, iExp, iTime, NORM, a,b,rk [ Emin-Emax from ABOVE]
    TMinuit* gMinuit = new TMinuit(NPAR);
    gMinuit->SetPrintLevel(-1);
    gMinuit->SetFCN(FCN);

    Double_t arglist[10];
    Int_t    ierflg = 0;
    arglist[0] = 1;
    gMinuit->mnexcm("SET ERR", arglist ,1,ierflg);

    double minK0= 1.0;
    double maxK0= 25.;
    double stepK0=1.e-2;
    double startK0=4.5; //4.5

  double minN0   = 0.9;
  double maxN0   = 2.;
  double stepN0  = 1.e-3;
  double startN0 =  0.98079;

    double mina= 0;
    double maxa= 2.;
    double stepa=1.e-2;
    double starta=1. ; // 1

    double minb= 0;
    double maxb= 2.;
    double stepb=1.e-2;
    double startb= 1; // 1

       double minrk= 100;
    double maxrk= 10.e3;
    double steprk=1;
    double startrk=3000;


  // ---- fitting loop: AMS02 ---
    cout<<"**** Fit to PAMELA ****"<<endl;

    int indTime=0;
    
    
   for(int tt=0;tt<nTimePAMELA_Proton;tt++){
    
//select BR o range di BR

//if( tt == 73) {
    //if( tt > -1 && tt < 33) {
      if( tt > -1) {
            cout<<endl<<endl;
            cout<<"**** Current BR ************************************** "<<tt<<endl;

            
            gMinuit->DefineParameter(0, "K_{0}", startK0, stepK0, minK0, maxK0);
            gMinuit->DefineParameter(1, "Experiment", iExp_PAMELA, 0, 0, 0);
            gMinuit->DefineParameter(2, "Time Index", tt, 0, 0, 0);
            gMinuit->DefineParameter(3, "NORM", startN0, stepN0,  minN0, maxN0);
            gMinuit->DefineParameter(4, "index", starta, stepa, mina, maxa);
            gMinuit->DefineParameter(5, "indexb", startb, stepb, minb, maxb);
              gMinuit->DefineParameter(6, "rk", startrk, steprk, minrk, maxrk);
            //gMinuit->DefineParameter(4, "a", starta, stepK0, minK0, maxK0);
            gMinuit->FixParameter(1); // FIX TO WHAT?
            gMinuit->FixParameter(2);

            gMinuit->FixParameter(3);
            //gMinuit->FixParameter(0);


          if(Kmodel == 1){
                gMinuit->FixParameter(4);
                gMinuit->FixParameter(5);
                gMinuit->FixParameter(6);
          }

             if(Kmodel == 2){
                gMinuit->FixParameter(5);
                gMinuit->FixParameter(6);
          }
            // se scelgo standard model il Parametro 4 (i.e. a) lo fisso esplicitamente per risparmiare costo computazionale

           if(Kmodel == 3){
                gMinuit->FixParameter(6);

           }
           
           // se scelgo standard model il Parametro 4 (i.e. a) lo fisso esplicitamente per risparmiare costo computazionale


            if(!kReNorm) gMinuit->FixParameter(3); // No Renormalization parameter
            
            arglist[0] = NITERATIONS;
            arglist[1] = 1.;
            gMinuit->mnexcm("MIGRAD", arglist ,2,ierflg);
            
            double bestK0= 0;
            double errK0 = 0.;

            double bestN0= 0;
            double errN0 = 0.;

      
            double besta = 0;
            double erra = 0.;
            double bestb = 0;
            double errb = 0.;

            double bestrk = 0;
            double errrk = 0.; 

            gMinuit->GetParameter(0,bestK0,errK0);
            gMinuit->GetParameter(3,bestN0,errN0);
            gMinuit->GetParameter(4,besta,erra);
            gMinuit->GetParameter(5,bestb,errb);
             gMinuit->GetParameter(6,bestrk,errrk);

             gMinuit->mnmatu(1);

            startK0= bestK0; // update for next fit
            startN0= bestN0; // update for next fit
            starta= besta; // update for next fit
            startb= bestb; // update for next fit
            
            cout <<" iTime: "<<tt<<"   K0= " << bestK0 << " +- "<<errK0<<" | "<<"    a: "<<besta<< " +- "<<erra <<" | b: "<<bestb<<" +- "<<errb<<" | N0 = "<<bestN0<<" +- "<<errN0<<" | rk: "<<bestrk<<" +- "<<errrk<<" |    Chi2: "<<bestChi2_PAMELA[tt]<<" / "<<NDF_PAMELA[tt]<<endl;
  
            // ---- put results into graph ----
            grK0vsTime_PAMELA->SetPoint(indTime, pameladata->xTimePAMELA_ProtonFlux[tt], bestK0);
            grK0vsTime_PAMELA->SetPointError(indTime, pameladata->eTimePAMELA_ProtonFlux, errK0);
            
   
            grAvsTime_PAMELA->SetPoint(indTime, pameladata->xTimePAMELA_ProtonFlux[tt], besta);
            grAvsTime_PAMELA->SetPointError(indTime, pameladata->eTimePAMELA_ProtonFlux, erra);
           
            grBvsTime_PAMELA->SetPoint(indTime, pameladata->xTimePAMELA_ProtonFlux[tt], bestb);
            grBvsTime_PAMELA->SetPointError(indTime, pameladata->eTimePAMELA_ProtonFlux, errb);
            
            grRkvsTime_PAMELA->SetPoint(indTime, pameladata->xTimePAMELA_ProtonFlux[tt], bestrk);
            grRkvsTime_PAMELA->SetPointError(indTime, pameladata->eTimePAMELA_ProtonFlux, errrk);
            
            grN0vsTime_PAMELA->SetPoint(indTime, pameladata->xTimePAMELA_ProtonFlux[tt], bestN0);
            grN0vsTime_PAMELA->SetPointError(indTime, pameladata->eTimePAMELA_ProtonFlux, errN0);
            
          
            grChis2vsTime_PAMELA->SetPoint(indTime, pameladata->xTimePAMELA_ProtonFlux[tt], bestChi2_PAMELA[tt]/ NDF_PAMELA[tt]);
            grChis2vsTime_PAMELA->SetPointError(indTime, pameladata->eTimePAMELA_ProtonFlux, 0);

          


            indTime++; // increment time index

            //NTSolution->PlotSolution_BR(tt);
        
       }
           
    }
    
    
  // ----- plot results ----
  TDatime Date1( 2000, 01, 01, 0, 0, 0);
  TDatime Date2( 2018, 01, 01, 0, 0, 0);
  double UTMIN= (double)Date1.Convert(); //unsigned int 10 cifre
  double UTMAX= (double)Date2.Convert();  //unsigned int 10 cifre
  double FMIN = 0.;
  double FMAX = 10.;


    
  TH2F* hFrameK0vsTime= new TH2F("hFrameK0vsTime","K0 vs Time",200, UTMIN, UTMAX, 500, FMIN, FMAX);
  SetStyleHistoVSTime(hFrameK0vsTime);
  hFrameK0vsTime->GetYaxis()->SetNdivisions(506);
  hFrameK0vsTime->GetYaxis()->SetTitleSize(0.07);
  hFrameK0vsTime->GetYaxis()->SetTitle("k_{0} ");
  hFrameK0vsTime->GetYaxis()->SetTitleOffset(0.50);
    
  



    
    TFile* outFile;
    if (Kmodel==1)
    {
       outFile= new TFile(Form(PATHm+"/out/K0_vs_time_StandardModel.root"),"recreate");
    }
    
       if (Kmodel==2){  outFile= new TFile(Form(PATHm+"/out/K0_vs_time_PowerLaw.root"),"recreate");}
         if (Kmodel==3){  outFile= new TFile(Form(PATHm+"/out/K0_vs_time_Potgeiter.root"),"recreate");}
         if (Kmodel==4){  outFile= new TFile(Form(PATHm+"/out/K0_vs_time_Potgeiter_Rk_free.root"),"recreate");}
    outFile->cd();
    
    grK0vsTime_PAMELA->Write();
    outFile->Write();
    outFile->Close();

    delete outFile;
    
  

    TFile* outFile2;
    if (Kmodel==1)
    {
       outFile2= new TFile(Form(PATHm+"/out/A_vs_time_StandardModel.root"),"recreate");
    }
    
    if (Kmodel==2){  outFile2= new TFile(Form(PATHm+"/out/A_vs_time_PowerLaw.root"),"recreate");}
    if (Kmodel==3){  outFile2= new TFile(Form(PATHm+"/out/A_vs_time_Potgeiter.root"),"recreate");}
 
    if (Kmodel==4){  outFile2= new TFile(Form(PATHm+"/out/A_vs_time_Potgeiter_Rk_free.root"),"recreate");}
    outFile2->cd();
    
    grAvsTime_PAMELA->Write();
    outFile2->Write();
    outFile2->Close();

    delete outFile2;
    

       TFile* outFile4;
    if (Kmodel==1)
    {
       outFile4= new TFile(Form(PATHm+"/out/B_vs_time_StandardModel.root"),"recreate");
    }
    
 if (Kmodel==2){  outFile4 = new TFile(Form(PATHm+"/out/B_vs_time_PowerLaw.root"),"recreate");}
if (Kmodel==3){  outFile4 = new TFile(Form(PATHm+"/out/B_vs_time_Potgeiter.root"),"recreate");}
  if (Kmodel==4){  outFile4 = new TFile(Form(PATHm+"/out/B_vs_time_Potgeiter_Rk_free.root"),"recreate");}
    outFile4->cd();
    
    grBvsTime_PAMELA->Write();
    outFile4->Write();
    outFile4->Close();

    delete outFile4;



    TFile* outFile3;
    if (Kmodel==1) {
       outFile3 = new TFile(Form(PATHm+"/out/Chis2_vs_time_StandardModel.root"),"recreate");
       }
    if (Kmodel==2){  outFile3= new TFile(Form(PATHm+"/out/Chis2_vs_time_PowerLaw.root"),"recreate");}
    if (Kmodel==3){  outFile3= new TFile(Form(PATHm+"/out/Chis2_vs_time_Potgeiter.root"),"recreate");}
   if (Kmodel==4){  outFile3= new TFile(Form(PATHm+"/out/Chis2_vs_time_Potgeiter_Rk_free.root"),"recreate");}
    outFile3->cd();
    
    grChis2vsTime_PAMELA->Write();
    outFile3->Write();
    outFile3->Close();

    delete outFile3;

  cout<<"first step"<<endl;


    TFile* outFile5;
    if (Kmodel==1) {
       outFile5 = new TFile(Form(PATHm+"/out/Rk_vs_time_StandardModel.root"),"recreate");
       }

    if (Kmodel==2){  outFile5= new TFile(Form(PATHm+"/out/Rk_vs_time_PowerLaw.root"),"recreate");}
    if (Kmodel==3){  outFile5= new TFile(Form(PATHm+"/out/Rk_vs_time_Potgeiter.root"),"recreate");}
 if (Kmodel==4){  outFile5= new TFile(Form(PATHm+"/out/Rk_vs_time_Potgeiter_Rk_free.root"),"recreate");}

    outFile5->cd();
    
    grRkvsTime_PAMELA->Write();
    outFile5->Write();
    outFile5->Close();

    delete outFile5;



  TFile* outFileN0;
    if (Kmodel==1) {
       outFileN0 = new TFile(Form(PATHm+"/out/N0_vs_time_StandardModel.root"),"recreate");
       }
    if (Kmodel==2){  outFileN0= new TFile(Form(PATHm+"/out/N0_vs_time_PowerLaw.root"),"recreate");}
    if (Kmodel==3){  outFileN0= new TFile(Form(PATHm+"/out/N0_vs_time_Potgeiter.root"),"recreate");}
   if (Kmodel==4){  outFileN0= new TFile(Form(PATHm+"/out/N0_vs_time_Potgeiter_Rk_free.root"),"recreate");}
    outFileN0->cd();
    
    grN0vsTime_PAMELA->Write();
    outFileN0->Write();
    outFileN0->Close();

    delete outFileN0;



  cout<<"first step"<<endl;


    
    // run App per single Draw
    Chi2_K0->SetMarkerStyle(8);
    Chi2_K0->SetMarkerSize(1);
   // Chi2_K0->Draw();
    
    
    // --- the end ---
    //stopM = high_resolution_clock::now();
    //auto durationM = duration_cast<seconds>(stopM - startM);
    stopM = chrono::high_resolution_clock::now();
    auto durationM = chrono::duration_cast<chrono::seconds>(stopM - startM);
    
    
    cout << "Gloabl TIme execution: " << durationM.count() <<"  seconds"<< endl;
    
    // theApp.Run();
    return 1;

}



void FCN(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag){
    double K0    = par[0];
    double N0  = par[3];
    double a = par[4];
    double b = par[5];
    double rk = par[6];
    
    // iExperiment: PAMELA / SOHO / BESS / BESSTeV / BESS00 / AMS02 / AMS02PRL2015
    int iExperiment   = (int)par[1];
    int iTime         = (int)par[2];
    
    cout<<endl<<endl;
    
   cout<<run<<" Current run :  K0=  "<<K0<<"| a="<<a<<" | b="<<b<<" | N0= "<<N0<<" | rk= "<<rk<<endl;
  
    run++;
    
    // --- EXP6: AMS02 ----
    if( iExperiment==iExp_PAMELA){
        NDF_PAMELA[iTime]=0;
     
      f= GetChiSquare(NTSolution,K0, iTime, iExperiment, N0,a,b,rk);
       // delete NTSolution;

        Chi2_K0->SetPoint(c,K0,f);
        c++;
    }
}


// ---- chisquare of spectrum | fixed epoch | given experiment | input K0 ----


double GetChiSquare(SolarModulation* NTSolution ,double K0, int iTime, int iExperiment, double N0,double a,double b,double rk) {
    
    NTSolution->Iterate_ThisKScale(K0,a,b,rk);
    

  cout<<"****************"<<endl;

      NTSolution->Iter_Modulation(); //set K0 in SolarModulation class
        
    
      start = chrono::high_resolution_clock::now();
      NTSolution->Solve();  //metodo lento
      stop = chrono::high_resolution_clock::now();
      auto duration = chrono::duration_cast<chrono::milliseconds>(stop - start);
      cout << "Solve TIme execution: " << duration.count() << " milliseconds"<<endl;
 
  
  double ChiSquare = 0.;
  double Emin = EMIN;
  double Emax = EMAX;
  double LnK0= log(K0);
    
  // --- EXP6: AMS02 ----
    if( iExperiment==iExp_PAMELA){
        NDF_PAMELA[iTime]=0;
        for(int ee=0;ee<pameladata->nEknPAMELA_Proton;ee++){ // ekn loop
            // get DATA
            
            double Energy= pameladata->grPAMELA_ProtonFluxVSEkn[iTime]->GetX()[ee];

 

            if(Energy<Emin || Energy>Emax) continue;  // Check Range
            
           
            double Ekn   = (1.e+3)*Energy;  // convert GeV->MeV
 
            //double Flux  = (1.e-3)*pameladata->grPAMELA_ProtonFluxVSEkn[iTime]->GetY()[ee];
            //double eFlux = (1.e-3)*pameladata->grPAMELA_ProtonFluxVSEkn[iTime]->GetEY()[ee];
     
            double Flux  = (1.e-3)*pameladata->grPAMELA_ProtonFluxVSEkn[iTime]->GetY()[ee];

            double eFlux = (1.e-3)*0.5*(pameladata->grPAMELA_ProtonFluxVSEkn[iTime]->GetEYhigh()[ee] +pameladata->grPAMELA_ProtonFluxVSEkn[iTime]->GetEYlow()[ee]) ;
     

    
            // get model
            double LnEkn= log(Ekn);
            //cout<<"Intrpolation****"<<endl;
            
            // Per interpolazione
            // double LnFluxModel= hLn_Flux_vs_KScale_vs_KEnergy_Z1_A1->Interpolate(LnEkn, LnK0);
         
            //per iterazione
            double LnFluxModel =  log(NTSolution->CurrentFlux(Ekn));
            

            double FluxModel= N0*exp(LnFluxModel);
            // get CHI2
            double Delta = (Flux - FluxModel)/eFlux;
            ChiSquare += (Delta*Delta);
            NDF_PAMELA[iTime]++;
    

    
        }

        cout<<" K0 = "<<K0<<"|  a = "<<a<<" | b= "<<b<<" | rk= "<<rk<<" | Chis2 = "<<ChiSquare<<endl;
      
  
        if( bestChi2_PAMELA[iTime] > ChiSquare ) bestChi2_PAMELA[iTime] = ChiSquare;
    }
    
    return ChiSquare;
}




// funzione per aspetto grafico dei plot finali

void SetStyleHistoVSEkn(TH2F* hh){
  hh->SetTitle(0);
  hh->SetLabelFont(42,"X");
  hh->SetLabelFont(42,"Y");
  hh->SetTitleFont(42,"X");
  hh->SetTitleFont(42,"Y");
  hh->GetXaxis()->SetTitle("kinetic energy (GeV)");

  //hh->GetYaxis()->CenterTitle();
  hh->GetXaxis()->SetNdivisions(513);
  hh->GetYaxis()->SetNdivisions(508);
  hh->GetXaxis()->SetLabelSize(0.05);
  hh->GetYaxis()->SetLabelSize(0.05);
  hh->GetYaxis()->SetTitleSize(0.05);
  hh->GetXaxis()->SetTitleSize(0.05);
  hh->GetXaxis()->SetTitleOffset(1.40);
  hh->GetYaxis()->SetTitleOffset(1.40);
  hh->GetXaxis()->SetLabelOffset(0.00);
  hh->GetYaxis()->SetLabelOffset(0.01);


}


void SetStyleHistoVSTime(TH2F* hh){
  hh->SetTitle(0);
  hh->GetXaxis()->SetTitle(0);
  //hh->GetYaxis()->SetTitle("J(R) [ GV^{ -1} m^{ -2} s^{ -1} sr^{ -1} ]");
  hh->GetXaxis()->SetTimeDisplay(1);
  hh->GetXaxis()->SetTimeFormat("#splitline{%Y}{%b,%d}%F1970-01-01 00:00:00s0");

  hh->GetXaxis()->SetNdivisions(513);
  hh->GetXaxis()->SetLabelSize(0.05);
  hh->GetYaxis()->SetNdivisions(508);
  hh->GetYaxis()->SetLabelSize(0.05);
  hh->GetYaxis()->SetTitleSize(0.05);
  hh->GetYaxis()->CenterTitle();
  hh->GetZaxis()->SetLabelSize(0.05);
  //hh->GetXaxis()->SetLabelOffset(0.06);
  //hh->GetYaxis()->SetTitleOffset(1.30);
}

void SetTimeGraphStyle(TGraphErrors* gr){
    gr->SetTitle(Form("Iteration From Lis model %d",kLISModel));
    gr->GetXaxis()->SetTitle(0);
    gr->GetYaxis()->SetTitle("K_{0} [10^{22} cm^{2} s^{-1}] ");
    
    gr->GetXaxis()->SetTimeDisplay(1);
    
    gr->GetXaxis()->SetTimeFormat("#splitline{%Y}{%b,%d}%F1970-01-01 00:00:00s0");
    gr->GetXaxis()->SetNdivisions(513);
    gr->GetXaxis()->SetLabelSize(0.06);
    gr->GetXaxis()->SetLabelOffset(0.06);
    gr->GetYaxis()->SetNdivisions(508);
    gr->GetYaxis()->SetLabelSize(0.06);
    gr->GetYaxis()->SetTitleOffset(0.65);
    gr->GetYaxis()->SetTitleSize(0.06);
    gr->GetYaxis()->CenterTitle();
}



void SetGraphK0Style(){
  /*grK0vsTime_PAMELA->SetName("grK0VSTime_PAMELA");
  grK0vsTime_PAMELA->SetMarkerStyle( crdata->grPAMELA_ProtonFluxVSTime[0]->GetMarkerStyle() );
  grK0vsTime_PAMELA->SetMarkerColor( crdata->grPAMELA_ProtonFluxVSTime[0]->GetMarkerColor() );
  grK0vsTime_PAMELA->SetMarkerSize( crdata->grPAMELA_ProtonFluxVSTime[0]->GetMarkerSize() );
  grK0vsTime_PAMELA->SetLineWidth( crdata->grPAMELA_ProtonFluxVSTime[0]->GetLineWidth() );
  grK0vsTime_PAMELA->SetLineColor( crdata->grPAMELA_ProtonFluxVSTime[0]->GetLineColor() );
  SetTimeGraphStyle(grK0vsTime_PAMELA);

  grK0vsTime_SOHO->SetName("grK0VSTime_SOHO");
  grK0vsTime_SOHO->SetMarkerStyle( crdata->grSOHO_ProtonFluxVSTime[0]->GetMarkerStyle() );
  grK0vsTime_SOHO->SetMarkerColor( crdata->grSOHO_ProtonFluxVSTime[0]->GetMarkerColor() );
  grK0vsTime_SOHO->SetMarkerSize( crdata->grSOHO_ProtonFluxVSTime[0]->GetMarkerSize() );
  grK0vsTime_SOHO->SetLineWidth( crdata->grSOHO_ProtonFluxVSTime[0]->GetLineWidth() );
  grK0vsTime_SOHO->SetLineColor( crdata->grSOHO_ProtonFluxVSTime[0]->GetLineColor() );
  SetTimeGraphStyle(grK0vsTime_SOHO);

  grK0vsTime_BESS->SetName("grK0VSTime_BESS");
  grK0vsTime_BESS->SetMarkerStyle( crdata->grBESS_ProtonFluxVSTime[0]->GetMarkerStyle() );
  grK0vsTime_BESS->SetMarkerColor( crdata->grBESS_ProtonFluxVSTime[0]->GetMarkerColor() );
  grK0vsTime_BESS->SetMarkerSize( crdata->grBESS_ProtonFluxVSTime[0]->GetMarkerSize() );
  grK0vsTime_BESS->SetLineWidth( crdata->grBESS_ProtonFluxVSTime[0]->GetLineWidth() );
  grK0vsTime_BESS->SetLineColor( crdata->grBESS_ProtonFluxVSTime[0]->GetLineColor() );
  SetTimeGraphStyle(grK0vsTime_BESS);
  
  grK0vsTime_BESSTeV->SetName("grK0VSTime_BESSTeV");
  grK0vsTime_BESSTeV->SetMarkerStyle( crdata->grBESSTeV_ProtonFluxVSTime[0]->GetMarkerStyle() );
  grK0vsTime_BESSTeV->SetMarkerColor( crdata->grBESSTeV_ProtonFluxVSTime[0]->GetMarkerColor() );
  grK0vsTime_BESSTeV->SetMarkerSize( crdata->grBESSTeV_ProtonFluxVSTime[0]->GetMarkerSize() );
  grK0vsTime_BESSTeV->SetLineWidth( crdata->grBESSTeV_ProtonFluxVSTime[0]->GetLineWidth() );
  grK0vsTime_BESSTeV->SetLineColor( crdata->grBESSTeV_ProtonFluxVSTime[0]->GetLineColor() );
  SetTimeGraphStyle(grK0vsTime_BESSTeV);

  grK0vsTime_BESS00->SetName("grK0VSTime_BESS00");
  grK0vsTime_BESS00->SetMarkerStyle( crdata->grBESS00_ProtonFluxVSTime[0]->GetMarkerStyle() );
  grK0vsTime_BESS00->SetMarkerColor( crdata->grBESS00_ProtonFluxVSTime[0]->GetMarkerColor() );
  grK0vsTime_BESS00->SetMarkerSize( crdata->grBESS00_ProtonFluxVSTime[0]->GetMarkerSize() );
  grK0vsTime_BESS00->SetLineWidth( crdata->grBESS00_ProtonFluxVSTime[0]->GetLineWidth() );
  grK0vsTime_BESS00->SetLineColor( crdata->grBESS00_ProtonFluxVSTime[0]->GetLineColor() );
  SetTimeGraphStyle(grK0vsTime_BESS00);
*/
    
  grK0vsTime_PAMELA->SetName("grK0VSTime_PAMELA");
  grK0vsTime_PAMELA->SetMarkerStyle( pameladata->grPAMELA_ProtonFluxVSTime[0]->GetMarkerStyle() );
  grK0vsTime_PAMELA->SetMarkerColor( pameladata->grPAMELA_ProtonFluxVSTime[0]->GetMarkerColor() );
  grK0vsTime_PAMELA->SetMarkerSize( pameladata->grPAMELA_ProtonFluxVSTime[0]->GetMarkerSize() );
  grK0vsTime_PAMELA->SetLineWidth( pameladata->grPAMELA_ProtonFluxVSTime[0]->GetLineWidth() );
  grK0vsTime_PAMELA->SetLineColor( pameladata->grPAMELA_ProtonFluxVSTime[0]->GetLineColor() );
  SetTimeGraphStyle(grK0vsTime_PAMELA);

    /*
  grK0vsTime_AMSPRL2015->SetName("grK0VSTime_AMSPRL2015");
  grK0vsTime_AMSPRL2015->SetLineWidth( 3 );
  grK0vsTime_AMSPRL2015->SetLineColor( kRed+2 );
  grK0vsTime_AMSPRL2015->SetMarkerSize( 0 );
  grK0vsTime_AMSPRL2015->SetMarkerColor( kRed+2 );
  SetTimeGraphStyle(grK0vsTime_AMSPRL2015);
     */
}



void SetGraphRkStyle(){

    
  grRkvsTime_PAMELA->SetName("grRkvsTime_PAMELA");

  grRkvsTime_PAMELA->SetMarkerStyle( pameladata->grPAMELA_ProtonFluxVSTime[0]->GetMarkerStyle() );
  grRkvsTime_PAMELA->SetMarkerColor( pameladata->grPAMELA_ProtonFluxVSTime[0]->GetMarkerColor() );
  grRkvsTime_PAMELA->SetMarkerSize( pameladata->grPAMELA_ProtonFluxVSTime[0]->GetMarkerSize() );
  grRkvsTime_PAMELA->SetLineWidth( pameladata->grPAMELA_ProtonFluxVSTime[0]->GetLineWidth() );
  grRkvsTime_PAMELA->SetLineColor( pameladata->grPAMELA_ProtonFluxVSTime[0]->GetLineColor() );
  SetTimeGraphStyle(grRkvsTime_PAMELA);
  grRkvsTime_PAMELA->GetYaxis()->SetTitle("Rk");

}


void SetGraphN0Style(){

    
  grN0vsTime_PAMELA->SetName("grN0vsTime_PAMELA");

  grN0vsTime_PAMELA->SetMarkerStyle( pameladata->grPAMELA_ProtonFluxVSTime[0]->GetMarkerStyle() );
  grN0vsTime_PAMELA->SetMarkerColor( pameladata->grPAMELA_ProtonFluxVSTime[0]->GetMarkerColor() );
  grN0vsTime_PAMELA->SetMarkerSize( pameladata->grPAMELA_ProtonFluxVSTime[0]->GetMarkerSize() );
  grN0vsTime_PAMELA->SetLineWidth( pameladata->grPAMELA_ProtonFluxVSTime[0]->GetLineWidth() );
  grN0vsTime_PAMELA->SetLineColor( pameladata->grPAMELA_ProtonFluxVSTime[0]->GetLineColor() );
  SetTimeGraphStyle(grN0vsTime_PAMELA);
  grN0vsTime_PAMELA->GetYaxis()->SetTitle("N0");

}



void SetGraphAStyle(){
/*
  grAvsTime_PAMELA->SetName("grAvsTime_PAMELA");
  grAvsTime_PAMELA->SetMarkerStyle( crdata->grPAMELA_ProtonFluxVSTime[0]->GetMarkerStyle() );
  grAvsTime_PAMELA->SetMarkerColor( crdata->grPAMELA_ProtonFluxVSTime[0]->GetMarkerColor() );
  grAvsTime_PAMELA->SetMarkerSize( crdata->grPAMELA_ProtonFluxVSTime[0]->GetMarkerSize() );
  grAvsTime_PAMELA->SetLineWidth( crdata->grPAMELA_ProtonFluxVSTime[0]->GetLineWidth() );
  grAvsTime_PAMELA->SetLineColor( crdata->grPAMELA_ProtonFluxVSTime[0]->GetLineColor() );
  SetTimeGraphStyle(grAvsTime_PAMELA);

  grN0vsTime_SOHO->SetName("grN0VSTime_SOHO");
  grN0vsTime_SOHO->SetMarkerStyle( crdata->grSOHO_ProtonFluxVSTime[0]->GetMarkerStyle() );
  grN0vsTime_SOHO->SetMarkerColor( crdata->grSOHO_ProtonFluxVSTime[0]->GetMarkerColor() );
  grN0vsTime_SOHO->SetMarkerSize( crdata->grSOHO_ProtonFluxVSTime[0]->GetMarkerSize() );
  grN0vsTime_SOHO->SetLineWidth( crdata->grSOHO_ProtonFluxVSTime[0]->GetLineWidth() );
  grN0vsTime_SOHO->SetLineColor( crdata->grSOHO_ProtonFluxVSTime[0]->GetLineColor() );
  SetTimeGraphStyle(grN0vsTime_SOHO);

  grN0vsTime_BESS->SetName("grN0VSTime_BESS");
  grN0vsTime_BESS->SetMarkerStyle( crdata->grBESS_ProtonFluxVSTime[0]->GetMarkerStyle() );
  grN0vsTime_BESS->SetMarkerColor( crdata->grBESS_ProtonFluxVSTime[0]->GetMarkerColor() );
  grN0vsTime_BESS->SetMarkerSize( crdata->grBESS_ProtonFluxVSTime[0]->GetMarkerSize() );
  grN0vsTime_BESS->SetLineWidth( crdata->grBESS_ProtonFluxVSTime[0]->GetLineWidth() );
  grN0vsTime_BESS->SetLineColor( crdata->grBESS_ProtonFluxVSTime[0]->GetLineColor() );
  SetTimeGraphStyle(grN0vsTime_BESS);
  
  grN0vsTime_BESSTeV->SetName("grN0VSTime_BESSTeV");
  grN0vsTime_BESSTeV->SetMarkerStyle( crdata->grBESSTeV_ProtonFluxVSTime[0]->GetMarkerStyle() );
  grN0vsTime_BESSTeV->SetMarkerColor( crdata->grBESSTeV_ProtonFluxVSTime[0]->GetMarkerColor() );
  grN0vsTime_BESSTeV->SetMarkerSize( crdata->grBESSTeV_ProtonFluxVSTime[0]->GetMarkerSize() );
  grN0vsTime_BESSTeV->SetLineWidth( crdata->grBESSTeV_ProtonFluxVSTime[0]->GetLineWidth() );
  grN0vsTime_BESSTeV->SetLineColor( crdata->grBESSTeV_ProtonFluxVSTime[0]->GetLineColor() );
  SetTimeGraphStyle(grN0vsTime_BESSTeV);

  grN0vsTime_BESS00->SetName("grN0VSTime_BESS00");
  grN0vsTime_BESS00->SetMarkerStyle( crdata->grBESS00_ProtonFluxVSTime[0]->GetMarkerStyle() );
  grN0vsTime_BESS00->SetMarkerColor( crdata->grBESS00_ProtonFluxVSTime[0]->GetMarkerColor() );
  grN0vsTime_BESS00->SetMarkerSize( crdata->grBESS00_ProtonFluxVSTime[0]->GetMarkerSize() );
  grN0vsTime_BESS00->SetLineWidth( crdata->grBESS00_ProtonFluxVSTime[0]->GetLineWidth() );
  grN0vsTime_BESS00->SetLineColor( crdata->grBESS00_ProtonFluxVSTime[0]->GetLineColor() );
  SetTimeGraphStyle(grN0vsTime_BESS00);
*/
    
  grAvsTime_PAMELA->SetName("grAvsTime_PAMELA");
  grAvsTime_PAMELA->SetMarkerStyle( pameladata->grPAMELA_ProtonFluxVSTime[0]->GetMarkerStyle() );
  grAvsTime_PAMELA->SetMarkerColor( pameladata->grPAMELA_ProtonFluxVSTime[0]->GetMarkerColor() );
  grAvsTime_PAMELA->SetMarkerSize( pameladata->grPAMELA_ProtonFluxVSTime[0]->GetMarkerSize() );
  grAvsTime_PAMELA->SetLineWidth( pameladata->grPAMELA_ProtonFluxVSTime[0]->GetLineWidth() );
  grAvsTime_PAMELA->SetLineColor( pameladata->grPAMELA_ProtonFluxVSTime[0]->GetLineColor() );
  SetTimeGraphStyle(grAvsTime_PAMELA);
  grAvsTime_PAMELA->GetYaxis()->SetTitle("a");

    /*
  grN0vsTime_AMSPRL2015->SetName("grN0VSTime_AMSPRL2015");
  grN0vsTime_AMSPRL2015->SetLineWidth( 3 );
  grN0vsTime_AMSPRL2015->SetLineColor( kRed+2 );
  grN0vsTime_AMSPRL2015->SetMarkerSize( 0 );
  grN0vsTime_AMSPRL2015->SetMarkerColor( kRed+2 );
  SetTimeGraphStyle(grN0vsTime_AMSPRL2015);*/
}



  


void SetGraphBStyle(){

    
  grBvsTime_PAMELA->SetName("grBvsTime_PAMELA");

  grBvsTime_PAMELA->SetMarkerStyle( pameladata->grPAMELA_ProtonFluxVSTime[0]->GetMarkerStyle() );
  grBvsTime_PAMELA->SetMarkerColor( pameladata->grPAMELA_ProtonFluxVSTime[0]->GetMarkerColor() );
  grBvsTime_PAMELA->SetMarkerSize( pameladata->grPAMELA_ProtonFluxVSTime[0]->GetMarkerSize() );
  grBvsTime_PAMELA->SetLineWidth( pameladata->grPAMELA_ProtonFluxVSTime[0]->GetLineWidth() );
  grBvsTime_PAMELA->SetLineColor( pameladata->grPAMELA_ProtonFluxVSTime[0]->GetLineColor() );
  SetTimeGraphStyle(grBvsTime_PAMELA);
  grBvsTime_PAMELA->GetYaxis()->SetTitle("b");

}





void SetGraphChis2Style(){

  grChis2vsTime_PAMELA->SetName("grChis2vsTime_PAMELA");

  grChis2vsTime_PAMELA->SetMarkerStyle( pameladata->grPAMELA_ProtonFluxVSTime[0]->GetMarkerStyle() );
  grChis2vsTime_PAMELA->SetMarkerColor( pameladata->grPAMELA_ProtonFluxVSTime[0]->GetMarkerColor() );
  grChis2vsTime_PAMELA->SetMarkerSize( pameladata->grPAMELA_ProtonFluxVSTime[0]->GetMarkerSize() );
  grChis2vsTime_PAMELA->SetLineWidth( pameladata->grPAMELA_ProtonFluxVSTime[0]->GetLineWidth() );
  grChis2vsTime_PAMELA->SetLineColor( pameladata->grPAMELA_ProtonFluxVSTime[0]->GetLineColor() );
  SetTimeGraphStyle(grChis2vsTime_PAMELA);
  grChis2vsTime_PAMELA->GetYaxis()->SetTitle("#chi^{2}/d.o.f.");

}