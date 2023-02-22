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
//#include "TMinuit.h"

int main(int argc, char **argv){

  // --- specify nucleus and LIS model option ----
  double Z = 1;
  double A = 1;
  int  LIS = 5; // [0 to 5] [CC2016ApJ | NT2017ApJ | NT2017PRD | NT2017ASR | REINA2020] [BOSCHINI]
  // Run mode [0: runID as input | 1: int indices of diffusion as input]
  int RunMode=0;

    // --- input ---
  int RID1; // RunID input
  int RID2; // RunID2 is used for making a scan between RunID-1 and RunID-2
  int iXiD; // index of Xi parameter [for drift, if enabled]
  int iK0;  // index of K0 parameter [for diffusion]
  int Kmodel;  // broken line or standard model for K function
  char plotmode;


//plot singolo con K0 e Drift selezionati dai 2 indici iXiD e iK0;
// Kmodel
// Plot vs E o R
  if(argc==1){
    RunMode=1;
    iXiD = 49; // zero drift [49]
    iK0  = 35; // avg diffusion level [35]
    Kmodel =  1; //standard model
    //Kmodel = 2; // // broken line model
    //plotmode = 'E'; // vs Ekn
     plotmode = 'E'; // vs Rig
      
    // argv[2] -> plot vs 'E' or plot vs 'R'
  }
    
// ---- three arguments: single run with indices XiD and K0 --> run tipico: ./bin/LxPgModel.exe 49 35 1
  if(argc==4){
    RunMode= 2;
    iXiD= atoi(argv[1]);
    iK0 = atoi(argv[2]);
    Kmodel = atoi(argv[3]);     // argv[3]) select DC model (standard or broken line)
    // argv[4] -> plot vs 'E' or plot vs 'R'
  }
      

// RunMode = 0 con definizione degli indici di Run
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

//***********    RUNMODE 0    ************************
//****************************************************
    
  if(RunMode==0){ // input is run ID: RID
    for(int RID=RID1; RID<RID2; RID++){
      SolarModulation* NTSolution= new SolarModulation(Z, A, LIS, RID);
      NTSolution->InitJLis();
      NTSolution->SetJLis();
      NTSolution->Solve();
      NTSolution->ForceFieldSolution();
      NTSolution->SetSolutions();
      NTSolution->PlotSolution();  // PLOT SINGLE RUN
      //NTSolution->StoreResults();  // STORE FOR SCANS
      delete NTSolution;
    }
  }

    
//***********    RUNMODE 1    ************************
//****************************************************


if(RunMode==1){ // input is indices for CR transport: iK0 and iXiD
    cout<<"§§§§§§§§§§§§§§   RunMode 1  §§§§§§§§§§§§§§§§"<<endl;
    SolarModulation* NTSolution= new SolarModulation(Z, A, LIS, iK0, iXiD, Kmodel);
    NTSolution->KscaleGrid();
    NTSolution->MakeGridsForHistograms();
    NTSolution->InitModulation();
    NTSolution->InitJLis();
    NTSolution->SetJLis();
    NTSolution->Solve();  
    NTSolution->ForceFieldSolution();
    NTSolution->SetSolutions();

    if (plotmode == 'R') {
          NTSolution->PlotSolution_Rig();
    }
    
    if (plotmode == 'E') {
            NTSolution->PlotSolution();
    }
    delete NTSolution;
  
}

    
//***********    RUNMODE 2    ************************
//****************************************************
// creation TH2F for fitting procedure
    
if(RunMode==2){
    
         SolarModulation* Solution_Grid= new SolarModulation(Z, A, LIS,iK0, iXiD,Kmodel);
         Solution_Grid->KscaleGrid();
         Solution_Grid->MakeGridsForHistograms();

        
         for (int ik=0; ik<Solution_Grid->nk; ik++) {
        //for (int ik=0; ik<6; ik++) {
             SolarModulation* NTSolution= new SolarModulation(Z, A, LIS,iK0, iXiD,Kmodel);
             cout<<"Start....Iteration    "<<ik<<endl;
             NTSolution->SetJLis(); //initjlis e setjlis
             NTSolution->KscaleGrid();
            
             //metodo che modifica ThisKscale ad ogni ciclo
             NTSolution->Iterate_ThisKScale(ik);
             NTSolution->InitModulation();
             NTSolution->Solve();
             NTSolution->ForceFieldSolution();
             NTSolution->SetSolutions();
             delete NTSolution;
            
        }
        
            Solution_Grid->StoreResults();  // STORE FOR SCANS

            //delete Solution_Grid
}
        

  return 1;

}

//end





