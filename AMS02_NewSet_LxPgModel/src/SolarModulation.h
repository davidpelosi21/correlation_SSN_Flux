#ifndef __SolarModulation__
#define __SolarModulation__

#include <stdio.h>
#include <cmath>
#include <cstdio>
#include "Constants.h"
#include <vector>

#include "Vec.h"

#include "FCMatrixBanded.h"
#include "EqSolver.h"
#include "TFormula.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TMath.h" // NT
#include "TApplication.h"
#include "TCanvas.h"
#include "TSystem.h"
#include "TGraph.h"

#include "TGraphErrors.h" // NT
#include "TGraphAsymmErrors.h"
#include "TH2F.h" // NT
#include "TH3F.h" // NT
#include "TFile.h" // NT

#include "TAxis.h"
#include "TF1.h"
#include "TLegend.h"
#include "TH1F.h"
#include "DataPoints.h"
#include "Spline3Interpolator.h"
#include "TRandom.h"
#include <ctime>

using namespace std;

class SolarModulation {

 public:
  SolarModulation(double Z, double A, int LIS,  int iKScale, double ThisKScale2,double a, double b,double rk,int iXiD, int Kmodel);
    
    //SolarModulation(double Z, double A, int LIS,  int iKScale, double ThisKScale2, int iXiD, int Kmodel);
    
  ~SolarModulation();
  void InitModulation();
    void Iter_Modulation();

    
    // definition of local interstellar spectrum
    void SetJLis(TF1*); // set LIS from input function
    void SetJLis();     // set LIS inside the class
    void InitJLis();    // initialize LIS model
    void Iterate_ThisKScale(double iKScale, double a,double b,double rk);    // initialize LIS model
    void KscaleGrid();    // define kscale grid

    // distances in A.U, rigidity in GeV
    void SetGridMatrices(); // defines, for each grid, two matrices that lead to the solution
    void Solve();
    void Solve_BC();

    double* GetRadiusGrid() {return radius;}
    double* GetlnRigidityGrid() {return lnRigidity;}
    double** GetSolution(){return f;}
    

    // NT handle solution
    void SetSolutions();
    void PlotSolution();
    void PlotSolution_Rig();

      void PlotSolution_BR(int br);

    void StoreResults();
    
    // print 4 different radius in AU
    void GetDifferentialIntensity(vector<int> R,string filename="", int mult = 6);

    void ForceFieldSolution(double phi); // phi imposed
    void ForceFieldSolution(); // phi calculated internally
    bool IsForceField = false;

    void ImportFile(string,string,string,string);
    void ImportFile(vector<string>);
    void BestParameter(double KCoeffInitial, int i);
    double GetBestParameter(){return BestPar;}

    void MakeGridsForHistograms();
    
    int nk;   // NT2017  kscale!
    
    double CurrentFlux(double E_c);
    double ThisKScale2 = -1.;
    double index = -1.;
    double indexb = -1.;
    double indexrk = -1.;
    double ThisKScale = -1.;
    
    
 private:
    
    EqSolver *S;
    Vec sol;
    Vec B;
    
    int RunID= -1;

    double ThisPhi    = -1.;
    double ThisXiD    = -1.;

    // servono?
    int indexKScale = -1;
    int indexPhi    = -1; // NOT same as KScale
    int indexXiD    = -1;
    int indexKmodel = 1; // select standard or broken-line model
    
    //np, nr number of grid points
    double** f;
    double** flux;
    int nr;   // radius
    int np;   // rigidity
//    int nk;   // NT2017  kscale!
    int nphi; // NT2017
    int nxid; // NT2018 drift

    double rmin, rmax;
    double pmin, pmax;
    double emin, emax; // NT ekn grid!
    double kmin, kmax; // NT2017 kscale!
    double phimin, phimax; // NT2017 phi
    double xidmin, xidmax; // NT2018 drift
    
    double deltar;
    double deltalogp;
    double deltalogk; // NT2017 kscale!
    double deltaphi;  // NT2016 phi
    double deltaxid;  // NT2018 drift
    
    double* radius=NULL;
    double* Rigidity=NULL; // NT2017
    double* KEnergy=NULL;
    double* KScale=NULL; // NT2017 kscale!

    double* lnRigidity=NULL;
    double* lnKEnergy=NULL;
    double* lnKScale=NULL; // NT2017 kscale!
    double* Phi=NULL; // NT2017
    double* XiD=NULL; // NT2018 drift
    
    FCmatrixFull** PA;

  //  EqSolver* S;
    
    TF1* SolarWind=NULL;
    TF1* JLisT=NULL;
    TF1* JLisRig=NULL;
    TF1* JLis_Scaled=NULL;
    TF1* KCoeff=NULL;
    TF1* ForceField=NULL;

    vector<TGraph*> DataSaver;
    double BestPar;
    
    // --- adimensional variables ---
    double Armin, Armax;
    double Apmin, Apmax;
    double deltaAr;

    // --- particle identity ---
    double pZ; // nuclear charge: Z
    double pA; // n of nucleons: A
    double pM; // particle mass: M = A*Mp

    
    // --- histogram grids ---
/*
    TH1F* hLIS_vs_Rigidity;
    TH1F* hLIS_vs_KEnergy;
    TH1F* hFlux_vs_Rigidity;
    TH1F* hFlux_vs_KEnergy;
    TH1F* hFlux_vs_Radius;
    TH1F* hFlux_vs_KScale;

    TH1F* hFluxFF_vs_Rigidity;
    TH1F* hFluxFF_vs_KEnergy;
    TH2F* hFluxFF_vs_Phi_vs_Rigidity;
    TH2F* hFluxFF_vs_Phi_vs_KEnergy;

    // lin-lin-lin histograms
    TH2F* hFlux_vs_KScale_vs_Rigidity;
    TH2F* hFlux_vs_KScale_vs_KEnergy;
    TH2F* hFlux_vs_Radius_vs_Rigidity;
    TH2F* hFlux_vs_Radius_vs_KEnergy;
    
    // log-log-log histograms
    TH2F* hLn_Flux_vs_KScale_vs_Rigidity;
    TH2F* hLn_Flux_vs_KScale_vs_KEnergy;
    TH2F* hLn_FluxFF_vs_Phi_vs_Rigidity; // log-R and log-Flux only
    TH2F* hLn_FluxFF_vs_Phi_vs_KEnergy; // log-E and log-Flux only

    // histograms with drift NT2018
    TH3F* hFlux_vs_XiD_vs_KScale_vs_Rigidity;
    TH3F* hFlux_vs_XiD_vs_KScale_vs_KEnergy;
    TH3F* hLn_Flux_vs_XiD_vs_KScale_vs_Rigidity;
    TH3F* hLn_Flux_vs_XiD_vs_KScale_vs_KEnergy;
*/
    
    // redundant!
    int nRadius;
    int nRigidity;
    int nKEnergy;
    int nKScale;
    int nPhi;
    int nXiD; // NT2018
    
    double* xRadius   = NULL;
    double* xRigidity = NULL;
    double* xKEnergy  = NULL;
    double* xKScale   = NULL;
    double* xPhi      = NULL;
    double* xXiD      = NULL; // NT2018 drift

    double* xLnRigidity = NULL;
    double* xLnKEnergy  = NULL;
    double* xLnKScale   = NULL;

    bool IsGrid= false;
    bool IsSolution= false;
  
    // ----------------------
    

    
    // auxiliar function: lis vs rigidity
    double fJLisRigvalue(double* x,double* par){
      double R= x[0];
      double E = RigToEkn(R);
      return (JLisT->Eval(E)*dEdR(R))*1E2;  //scaled by a factor of 100 [NT: why?]
    }

    
    // auxiliary function: lis scaled
    double fJLis_ScaledValue(double* x,double* par){ //scaled by a factor of 100 [NT: why?]
      return JLisT->Eval(x[0])*1E2;
    }

    
    // transform rigidity R=p/Z into kinetic energy per nucleon E=T/A
    double RigToEkn(double R){return (sqrt( (R*pZ)*(R*pZ) + pM*pM )-pM)/pA;}
    
    // transform kinetic energy per nucleon E=T/A into rigidity R=p/Z
    double EknToRig(double E){return (sqrt( (E*pA+pM)*(E*pA+pM)- (pM*pM)))/pZ;}
    
    // jacobian J(R) = dEkn/dRig vs Rig = beta
    double dEdR(double Rig){
      return (pZ/pA)*(pZ/pA)*Rig/sqrt((pZ/pA)*(pZ/pA)*(Rig*Rig)+ProtonMass*ProtonMass);
    }
    
    // jacobian J(E) = dRig/dEkn vs Ekn = 1/beta
    double dRdE(double Ekn){
      return (pA/pZ)*(ProtonMass + Ekn)/sqrt(Ekn*Ekn + 2.*ProtonMass*Ekn);
    }

    
    // flux array conversion rig->ekn [not used but OK]
    void FluxRig2Ekn(double* specRig, double* specEkn){ // Grid gia definite
      for(int ee=0;ee<np;ee++) specEkn[ee]= specRig[ee]*dRdE(KEnergy[ee]);
    }

    // flux array conversion ekn->rig [not used but OK]
    void FluxEkn2Rig(double* specEkn, double* specRig){ // Grid gia definite
      for(int ee=0;ee<np;ee++) specRig[ee]= specEkn[ee]*dEdR(Rigidity[ee]);
    }




    // ---- force field for all Z/A nuclei ----
    double functorFF(double* x, double* par){
      double phi= par[0];
      double PHI = (pZ/pA)*phi;

      // --- Standard FF: NT PRD 96 103005 2017 ---
      double E_TOA = x[0];
      double E_LIS = E_TOA + PHI;
      double factor = E_TOA*(E_TOA + 2*ProtonMass)/(E_LIS*(E_LIS + 2*ProtonMass));
      return JLisT->Eval(E_LIS)*factor;
    }


    double fSpline_1AU(double* x, double* par){
      Spline3Interpolator G(np,KEnergy,f[1]);
      return G.Interpolate(x[0]);
    }

    
    // --- stuff for LIS from external graphs ----
    bool IsLIS = false;
    int kLISModel; // default: Mod0!
    double* xLnEkn;
    double* xLnLIS;
    int     nLIS;

    // ---- initialize LIS models ----
    // void InitJLis(); // moved as public
    

    
    // --- functors for LIS ----
    double functorLISvsEKN(double* x, double* par){
      double funct= 0;
      if(kLISModel==0) funct= functorLISvsEKN_Mod0(x, par);
      if(kLISModel>0 ) funct= functorLISvsEKN_Mod123(x, par);
      return funct;
    }

    double functorLISvsRIG(double* x, double* par){
      double funct=0;
      if(kLISModel==0) funct= functorLISvsRIG_Mod0(x, par);
      if(kLISModel>0 ) funct= functorLISvsRIG_Mod123(x, par);
      return funct;
    }


    // ---- LIS vs EKN models from arrays/literature ----
    double functorLISvsEKN_Mod123(double* x, double* par){
      double LnE= log(x[0]);
      if(LnE< xLnEkn[0]) LnE=xLnEkn[0];
      if(LnE> xLnEkn[nLIS-1]) LnE=0.999*xLnEkn[nLIS-1];
      
      // find index
      int indE=-1;
      for(int ee=0;ee<nLIS-1;ee++){
    if(LnE>=xLnEkn[ee] && LnE<=xLnEkn[ee+1]){
      indE= ee;
      break;
    }
      }

      // interpolate linear
      double X1= xLnEkn[indE];
      double X2= xLnEkn[indE+1];
      double Y1= xLnLIS[indE];
      double Y2= xLnLIS[indE+1];
      double LnLIS= Y1+ ((Y2-Y1)/(X2-X1))*(LnE - X1);
      double LIS= exp( LnLIS );
      return LIS;
    }


    // ---- LIS vs RIG models from arrays/literature ----
    double functorLISvsRIG_Mod123(double* x, double* par) {
      double R = x[0];
      double E =  RigToEkn(R);
      double JacR = dEdR(R);
      double* EE = new double[1];
      EE[0]= E;
      double LIS_EKN= functorLISvsEKN(EE, par);
      double LIS_RIG = JacR*LIS_EKN;
      return LIS_RIG;
    }


    // --- Mod0: functor for LIS parametric function a la Corti ----
    double functorLISvsEKN_Mod0(double* x, double* par){
      double E= (1.e-3)*x[0]; // MeV->GeV
      int Z= (int)pZ;
      int A= (int)pA;
      
      // --- mass for LIS calculation purposes ----
      double rA= pA; // mass number we use for LIS
      if(Z==2) rA=4; // NT! for He3 hypothesis, we just apply a correction below.
      
      double Mp=0.938;
      double pM=Mp*rA;
      
      // rigidity [use rA not pA, for trick]
      double R= (sqrt( (E*rA+pM)*(E*rA+pM)- (pM*pM)))/pZ;
      double LnR= TMath::Log(R);
      
      // jacobians dR/dE and dE/dR
      // double dEdR = (pZ/rA)*(pZ/rA)*R/sqrt((pZ/rA)*(pZ/rA)*(R*R)+Mp*Mp);
      double _dRdE= (rA/pZ)*(Mp + E)/sqrt(E*E + 2.*Mp*E);
      
      // low-R parameterization
      double Norm= 11600; // 11740.;
      double mu= -0.559;
      double sigma= 0.563;
      double G1 = -2.4482;
      double nu= 0.431;
      
      if( Z==2 ){ // He4 or He3

    // STD: He un po troppo basso a LE [GV]
    //Norm= 1740.;
    //mu= -0.075;
    //sigma= 0.495;
    //G1 = -2.364;

    // ALT: He un po' piu alto a LE [GV] [ma basso a 10 MeV]
    Norm= 1740.;
    mu= -0.090;
    sigma= 0.44;
    G1 = -2.364;

      }

      
      double LowR= (TMath::Power( 1.+ TMath::Exp(-(LnR-mu)/sigma), -1./nu ))*TMath::Power(R, G1);
      
      // high-R parameterization
      double Rb1 = 6.2;
      double DG1 = -0.4227;
      double S1  = -0.108;
      double Rb2 = 545.;
      double DG2 = -0.6;
      double S2  = -0.4;
      
      double HighR= TMath::Power(1.+ TMath::Power(  (R/Rb1)*TMath::Power(1.+ TMath::Power(R/Rb2, DG2/S2), S2) , DG1/S1 ), S1);
            
      ///double LIS_RIG =  Norm*LowR*HighR;
      double Total = Norm*LowR*HighR;
      
      // He3/He4 ratio
      double Rb=1.2;
      double N0=0.18;

      // NOT BAD
      //double GammaHE= 0.33; // HE
      //double GammaLE= -0.1; // LE

      // TEST: STEEPER
      double GammaHE= 0.45; // HE
      double GammaLE= -0.07; // LE

      double S=10;
      double Ratio= N0*pow(R, -GammaLE)*(pow( (1. + pow(R/Rb, S)), -(GammaHE-GammaLE)/S));

      // ---- scale by 3He/4He if A==3 ----
      double LIS_RIG= Total; //
      if(A==3) LIS_RIG= Ratio*Total/(1.+Ratio); // tot 3He flux
      if(A==4) LIS_RIG= Total/(1.+Ratio); // tot 4He flux
      // ----------------------------------
      
      double LIS_EKN = _dRdE*LIS_RIG; //GeV^-1
      double LIS_MEV = LIS_EKN*(1.e-3); // MeV^-1
      
      return LIS_MEV;
    }


    // --- Mod0: LIS vs RIG parametric functor ---
    double functorLISvsRIG_Mod0(double* x, double* par){
      
      int Z= (int)pZ;
      double R= (1.e-3)*x[0]; // MeV->GeV
      double LnR= TMath::Log(R);
      
      // --- jacobians dR/dE and dE/dR ---
      // double dEdR = (pZ/pA)*(pZ/pA)*R/sqrt((pZ/pA)*(pZ/pA)*(R*R)+Mp*Mp);
      // double dRdE= (pA/pZ)*(Mp + E)/sqrt(E*E + 2.*Mp*E);
      

      // low-R parameterization
      double Norm= 11600; // 11740.;
      double mu= -0.559;
      double sigma= 0.563;
      double G1 = -2.4482;
      double nu= 0.431;
      
      if( Z==2 ){ // He4 or He3
    Norm= 1740.; // 1800.;
    mu= -0.075; //-0.08;
    sigma= 0.495; //0.52;
    G1 = -2.364; // -2.37
      }

      
      double LowR= (TMath::Power( 1.+ TMath::Exp(-(LnR-mu)/sigma), -1./nu ))*TMath::Power(R, G1);
      
      // high-R parameterization
      double Rb1 = 6.2;
      double DG1 = -0.4227;
      double S1  = -0.108;
      double Rb2 = 545.;
      double DG2 = -0.6;
      double S2  = -0.4;
  
      double HighR= TMath::Power(1.+ TMath::Power(  (R/Rb1)*TMath::Power(1.+ TMath::Power(R/Rb2, DG2/S2), S2) , DG1/S1 ), S1);
      
      
      double LIS_RIG =  Norm*LowR*HighR; // GV^-1
      double LIS_MV  = LIS_RIG*(1.e-3); // MV^-1
      return LIS_MV;
    }
   
    
    // --- other functors: not used but they could be included here ----
    void SetKCoeff(); // init TF1
    void SetKCoeffP(); // init TF1 broken line( potgeiter 2013)
     void SetKCoeffP2(); // init TF1 broken line( potgeiter 2013) rk free

      void SetKCoeffPower(); // init TF1 broken line( David 2023)
    void SetSolarWind(); // init TF1
    

    /* void SetKScale(); // set KCoeff scale from RunID */
    /* void SetPhi();    // compute Phi from KScale */
    /* void SetXiD();    // set drift level from RunID */

    
    
    // ---- wind profile. default: radially constant ----
    double functorSW(double* x, double* par){
      double V= par[0];
      if(V<0) V = (1.e+3)*Vw; // SCALE HERE // default from Constants.h
      return V; // [m/s]
    }


    // ---- diffusion coefficient: standard form ----
    double functorDC(double* x, double* par){
      double Norm     = (1.e+18)*par[0]; // norm factor
      if(Norm<0) Norm = 4.38e+18; // usare ThisKScale???
      if(indexKScale>-1) Norm = (1.e+18)*ThisKScale2; // SCALE HERE !!!
      double Z = pZ;   // charge, class member
      double A = pA;   // mass number, class member;
      double Rigidity  = x[0];           // rigidity MV
      double Momentum  = Rigidity*pZ;    // momentum MeV/c
      double Mass      = pA*ProtonMass;  // or pM [MeV/c^2]
      double Energy    = sqrt(Momentum*Momentum + Mass*Mass);
      double Beta      = Momentum/Energy;
      return Norm*(1.e-3)*Rigidity*Beta; // diffusion coeff [m^2/s]
      // NT: factor 1.e-3 is to have (R/1GV)?
    }

    // ---- diffusion coefficient: DAVID 2023 ----
    double functorDCPower(double* x, double* par){

      //double Norm  = (1.e+18)*par[0]; // norm factor
      double Norm  = (1.e+18)*ThisKScale2; // norm factor
      //cout<<par[0]<<"  "<<ThisKScale2<<endl;
      double a = index;

      //if(Norm<0) Norm = 4.38e+18; // usare ThisKScale???
      
      //if(indexKScale>-1) Norm = (1.e+18)*ThisKScale2; // SCALE HERE !!!

      double Z = pZ;   // charge, class member
      double A = pA;   // mass number, class member;
      double Rigidity  = x[0];           // rigidity MV
      double Momentum  = Rigidity*pZ;    // momentum MeV/c
      double Mass      = pA*ProtonMass;  // or pM [MeV/c^2]
      double Energy    = sqrt(Momentum*Momentum + Mass*Mass);
      double Beta      = Momentum/Energy;
       
      return Norm*pow((1.e-3)*Rigidity,a)*Beta; // diffusion coeff [m^2/s]
     
      
      // NT: factor 1.e-3 is to have (R/1GV)?


    }



    // --- diffusion coefficient Potgieter 2013 ---
    double functorDCP(double* x, double* par){

      double Norm     = (1.e+18)*par[0]; // norm factor [m2/s]
      if(Norm<0) Norm = 4.38e+18; // usare ThisKScale???
      if(indexKScale>-1) Norm = (1.e+18)*ThisKScale2; // SCALE HERE !!!

        double a = index;
        double b = indexb;

      double Rigidity  = x[0];           // rigidity MV */
      double Momentum  = Rigidity*pZ;    // momentum MeV/c
      double Mass      = pA*ProtonMass;  // or pM [MeV/c^2]
      double Energy    = sqrt(Momentum*Momentum + Mass*Mass);
      double Beta      = Momentum/Energy;

      double c=3; // 3.0

      //double Pk=4.e+3; // MV
      double Pk=2802.2; // MV

      double P0=1.e+3; // MV
      double K0=1.;
      double P = Rigidity; // MV
      double B=5.05; // nT
      double B0=1; // nT

      double K= Norm*Beta*pow(P/P0,a)*pow(((pow(P/P0,c) + pow(Pk/P0,c))/( 1.+pow(Pk/P0,c) ) ), (b-a)/c  );
      return K;
    }

    



    // --- diffusion coefficient Potgieter 2013 rk free ---
    double functorDCP2(double* x, double* par){

      double Norm     = (1.e+18)*par[0]; // norm factor [m2/s]
      if(Norm<0) Norm = 4.38e+18; // usare ThisKScale???
      if(indexKScale>-1) Norm = (1.e+18)*ThisKScale2; // SCALE HERE !!!

        double a = index;
        double b = indexb;
        double Pk = indexrk;

      double Rigidity  = x[0];           // rigidity MV */
      double Momentum  = Rigidity*pZ;    // momentum MeV/c
      double Mass      = pA*ProtonMass;  // or pM [MeV/c^2]
      double Energy    = sqrt(Momentum*Momentum + Mass*Mass);
      double Beta      = Momentum/Energy;

      double c=3; // 3.0
      //double Pk=4.e+3; // MV
      double P0=1.e+3; // MV
      double K0=1.;
      double P = Rigidity; // MV
      double B=5.05; // nT
      double B0=1; // nT

      double K= Norm*Beta*pow(P/P0,a)*pow(((pow(P/P0,c) + pow(Pk/P0,c))/( 1.+pow(Pk/P0,c) ) ), (b-a)/c  );
      return K;
    }









    // ---- phenomenological drift ----
    double GetXiD(){
      double XiD= ThisXiD;
      if(indexXiD<0 || DRIFT==0) XiD=0.;
      if(fabs(XiD)<1.e-6)XiD=0;
      return XiD;
    }



    // alternative constructor where Wind and KCoeff are defined internally
    // SolarModulation(double frmin, double frmax, int N_R, double fpmin, double fpmax, int N_P, double Z, double A);


    // --- beta calculations / not used ---
    /*
    double RigToBeta(double R){
      double E= RigToEkn(R);
      double Gamma = E/ProtonMass + 1.;
      double Beta= sqrt(1. - 1./(Gamma*Gamma));
      return Beta;
    }

    double BetaToRig(double B){
      double gamma = 1./sqrt(1. - (B*B));
      double E= (gamma-1.)*ProtonMass;
      double R= EknToRig(E);
      return R;
    }

    double EknToBeta(double E){
      double Gamma = E/ProtonMass + 1.;
      double Beta= sqrt(1. - 1./(Gamma*Gamma));
      return Beta;
    }

    double BetaToEkn(double B){
      double gamma = 1./sqrt(1. - (B*B));
      double E= (gamma-1.)*ProtonMass;
      return E;
    }
    */

    
};

#endif



