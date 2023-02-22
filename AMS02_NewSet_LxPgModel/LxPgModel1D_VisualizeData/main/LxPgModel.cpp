#include "SolarModulation.h"
/*
 ****************************************** ****************************************** ******************************************
 ****************************************** ****************************************** ******************************************
 ****************************************** ****************************************** ******************************************
 
 Run solo per graficare dati di voyager e ams (p o He) e flussi LIS (0-5), Flusso modulato (grid numerica) e ForceField Solution
 Flusso in Ekn o Rig
 
 ****************************************** **************************************** ******************************************
 
 
 */

#include "SolarModulation.h"

int main(int argc, char **argv){

    //nel caso multigraph
    int mult = 5;
    int ar = 90*mult;
    int ap = 200; //90*mult;
    /*
   // // distance [AU]
    double rmin = 0;
    double rmax = 90*AU;
   // // rigidity [MV]
    double pmin = 1.0;
    double pmax = 500.e+3; // 500 GV 10E6;
    */
    vector<int> radius;
    radius.push_back(1*mult);
    radius.push_back(10*mult);
    radius.push_back(40*mult);
    radius.push_back(60*mult);
    radius.push_back(80*mult);


  // --- specify nucleus and LIS model option ----
  double Z = 1;
  double A = 1;
  int  LIS = 2; // [0 to 4] [CC2016ApJ | NT2017ApJ | NT2017PRD | NT2017ASR | REINA2020]
  // Run mode [0: runID as input | 1: int indices of diffusion as input]
  int RunMode=0;

  // --- input ---
  int RID1; // RunID input
  int RID2; // RunID2 is used for making a scan between RunID-1 and RunID-2
  int iXiD; // index of Xi parameter [for drift, if enabled]
  int iK0;  // index of K0 parameter [for diffusion]
  int Kmodel;  // broken line or standard model
  //  double Alfa=1.6;
  //  double Beta=1;
  // the run dependes on the input arguments, as follows
  // ---- zero arguments: use default indices ----
    //se al momento del run non passo stringhe, argc=1, e.g. ./LxPgModel
  if(argc==1){
    RunMode=1;
    iXiD = 49; // zero drift [49]
    iK0  = 35; // avg diffusion level [35]
      // funziona se nel main è cosi:
        // // grid sizebl
  }
// ---- one argument: single run with RunID as input ----
//se passo una stringa al momento del run " ./LxPgModel 2 "  -> argc = 2  e argv[1] = 2, atoi serve a convertire la stringa di input in numero
  if(argc==2){
    RunMode= 0;
    RID1= atoi(argv[1]);
    RID2= RID1+1;
  }
  // ---- two arguments: scan between Run1 and Run2 ----
  if(argc==3){
    RunMode= 0;
    RID1= atoi(argv[1]);
    RID2= atoi(argv[2]);
  }
  // ---- three arguments: single run with indices XiD and K0 ----
  if(argc==5){
    RunMode= 1;
    iXiD= atoi(argv[1]);
    iK0 = atoi(argv[2]);
    Kmodel = atoi(argv[3]);
    // argv[3]) select DC model (standard or borken line)
    // the third argument-> plot vs E or plot vs Rig
  }
  // --- TO PLOT OR STORE SOLUTION, MODIFY BELOW ---
  // PLOT MAKES SENSE IS FOR SINGLE RUNS. FOR SCAN, YOU STORE RESULTS
  if(RunMode==0){ // input is run ID: RID
    for(int RID=RID1; RID<RID2; RID++){
      SolarModulation* NTSolution= new SolarModulation(Z, A, LIS, RID);
      NTSolution->InitJLis();
      NTSolution->SetJLis();
      NTSolution->Solve();
      NTSolution->ForceFieldSolution();
      NTSolution->SetSolutions(); // NT this is also called inside PlotSolution
      NTSolution->PlotSolution();  // PLOT SINGLE RUN
      //NTSolution->StoreResults();  // STORE FOR SCANS
      delete NTSolution;
    }
  }

  // --- by Transport indices ----
  if(RunMode==1){ // input is indices for CR transport: iK0 and iXiD
      
    SolarModulation* NTSolution= new SolarModulation(Z, A, LIS, iK0, iXiD,Kmodel);
    NTSolution->InitJLis();
    NTSolution->SetJLis();
    NTSolution->Solve();
    NTSolution->ForceFieldSolution();
    NTSolution->SetSolutions();
      
    if (string(argv[4]) == 'R') {
          NTSolution->PlotSolution_Rig();
    }
    
    if (string(argv[4]) == 'E') {
            NTSolution->PlotSolution();
    }
    
      if (string(argv[4]) == 'S') {
          NTSolution->StoreResults();  // STORE FOR SCANS
      }
    //NTSolution->PlotSolution(); // PLOT SINGLE RUN
    //NTSolution->StoreResults();  // STORE FOR SCANS
    //NTSolution->GetDifferentialIntensity(radius,mult);
    delete NTSolution;
  }

    
    // completo gli istogrammi
    if(RunMode==2){ // input is indices for CR transport: iK0 and iXiD
        
      SolarModulation* NTSolution= new SolarModulation(Z, A, LIS, iK0, iXiD,Kmodel);
      NTSolution->InitJLis();
      NTSolution->SetJLis();
      NTSolution->Solve();
      NTSolution->ForceFieldSolution();
      NTSolution->SetSolutions();
        
      if (string(argv[4]) == 'R') {
            NTSolution->PlotSolution_Rig();
      }
      
      if (string(argv[4]) == 'E') {
              NTSolution->PlotSolution();
      }
      
        if (string(argv[4]) == 'S') {
            NTSolution->StoreResults();  // STORE FOR SCANS
        }
      //NTSolution->PlotSolution(); // PLOT SINGLE RUN
      //NTSolution->StoreResults();  // STORE FOR SCANS
      //NTSolution->GetDifferentialIntensity(radius,mult);
      delete NTSolution;
    }
  return 1;

}





