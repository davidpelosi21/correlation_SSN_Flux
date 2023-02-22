#include "SolarModulation.h"
#include "EqSolver.h"
#include <chrono>

TString PATH = "../LxPgModel1D_set_TH2F";

using namespace std::chrono;
//globali
// --- histogram grids ---
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


// ---- NT Constructor ---
SolarModulation::SolarModulation(double Z, double A, int LIS, int iKScale, double ThisKScale2, double a, double b,double rk,int iXiD, int Kmodel){

  S = new EqSolver();
  
  Vec* Par = new Vec(608);
  B = *Par;
    
  pZ = Z;
  pA = A;
  pM = ProtonMass*A;
  kLISModel = LIS;
  indexKScale= iKScale;
  indexXiD   = iXiD;
  indexKmodel = Kmodel;
  index = a;
  indexb = b;
  indexrk = rk;
   // indexAlfa = Alfa;
   // indexBeta = Beta;

  JLisT = NULL;
  KCoeff= NULL;
  SolarWind= NULL;

  rmax = RMAX;
  rmin = RMIN;
  nr   = NR;

  pmax = PMAX;
  pmin = PMIN;
  np   = NP;

  emin = RigToEkn(pmin);
  emax = RigToEkn(pmax);

  Armax = rmax/r0;
  Armin = rmin/r0;

  kmax= KMAX;
  kmin= KMIN;
  nk  = NK;

  phimax= PHIxKSCALE/kmin;
  phimin= PHIxKSCALE/kmax;
  nphi  = NK;

  xidmax= XMAX;
  xidmin= XMIN;
  nxid  = NX;

  // --- set RunID index ---
  if(indexKScale>-1 && indexXiD>-1){
    int indRID=-1;
    for(int dd=0;dd<NX;dd++){
      for(int kk=0;kk<NK;kk++){
          indRID++;
          if( kk==indexKScale && dd==indexXiD ) RunID=indRID;
      }
    }
      //printf("enter in first if \n");
  }
    if(indexXiD<0 || DRIFT==0) {RunID = indexKScale;}
}



void SolarModulation::InitModulation(){
    // --- set diffusion parameters ----
    ThisKScale = 6.;
    ThisPhi    = 407.;
    ThisXiD    = 0.;
    
    if(indexKScale>-1) ThisKScale= KScale[indexKScale];
    if(indexKScale>-1) ThisPhi= PHIxKSCALE/ThisKScale;
    if(indexXiD>-1 && DRIFT>0) ThisXiD = XiD[indexXiD];
    if( fabs(ThisXiD)< 1.e-6 ) ThisXiD=0.;  //fabs returns abs
    
    // --- init functions ----
    if ( indexKmodel == 1) {
        SetKCoeff();
    }
    if (indexKmodel == 2) {  SetKCoeffPower(); }
    if (indexKmodel == 3) {  SetKCoeffP(); }
    if (indexKmodel == 4) {  SetKCoeffP2(); }

    SetSolarWind();
    // SetXiD(); // no TF1 for XiD
    
    // --- creating empty matrix for saving solution f(r,lnp) ---
    f = new double*[nr];
    flux = new double*[nr];
    for(int i =0; i<nr; ++i){
        f[i] = new double[np];
        flux[i] = new double[np];
    }
    for(int i=0; i<nr; ++i){
        for(int j=0; j<np;++j){
            f[i][j] = 0.;
            flux[i][j] = 0.;
        }
    }

    // --- creating auxiliar matrix for setting the np-1 matrices ---
    double** aux = new double*[nr-2];   //array di puntatori
    for(int i=0; i<nr-2; ++i) aux[i] = new double[nr-2];
    for(int i=0; i<nr-2; ++i){
        for(int j=0; j<nr-2; ++j) aux[i][j] = 0.;
    }

    // --- equation to solve: A U(n+1) = B U(n) n: energy index ---
    // PA : array of matrices for every energy iteration equation

    PA = new FCmatrixFull*[np-1];
    for(int i=0; i<np-1; ++i) PA[i] = new FCmatrixFull(aux, nr-2, nr-2);

    for(int j=0; j<np-1;++j){
        double p = exp(lnRigidity[j]); //current rigidity

    // NT2018 use DRIFT via GetXiD
        for(int i=1; i<nr-1; ++i){ //loop on radius
            double r = radius[i]; //current radius
            double A = KCoeff->Eval(p)/(r0*r0); // Z/A set?
            double B = (-(1.+GetXiD())*SolarWind->Eval(r) + 2.*KCoeff->Eval(p)/r)/r0;
            double C = (2./3.)*SolarWind->Eval(r)/(r);

            // normal equations for i(radius) = [1,nr-2]
            if(i!=1 && i!=nr-2){
                (*(PA[j]))[i-1][i-2] = A/(2.*deltaAr*deltaAr) - B/(4.*deltaAr);
                (*(PA[j]))[i-1][i-1] = -A/(deltaAr*deltaAr) - C/deltalogp;
                (*(PA[j]))[i-1][i] =  A/(2.*deltaAr*deltaAr) + B/(4.*deltaAr);                 }

            // different for the first equation
            else if(i==1){
                (*(PA[j]))[i-1][i-1] = -A/(deltaAr*deltaAr) - C/deltalogp;
                (*(PA[j]))[i-1][i] =  A/(2.*deltaAr*deltaAr) + B/(4.*deltaAr);
            }

            // and for the last one
            else{
                (*(PA[j]))[i-1][i-2] = A/(2.*deltaAr*deltaAr) - B/(4.*deltaAr);
                (*(PA[j]))[i-1][i-1] = -A/(deltaAr*deltaAr) - C/deltalogp;
            }
        }
    }
    
    for(int i=0; i<nr-2; ++i) delete[] aux[i];
    delete[] aux;
    
    if ( indexKmodel == 1) {
        SetKCoeff();
    }
    if (indexKmodel == 2) {  SetKCoeffPower(); }
    if (indexKmodel == 3) {  SetKCoeffP(); }
    if (indexKmodel == 4) {  SetKCoeffP2(); }

    SetSolarWind();
    // SetXiD(); // no TF1 for XiD

}


void SolarModulation::Iter_Modulation(){
   
    for(int j=0; j<np-1;++j){
        double p = exp(lnRigidity[j]); //current rigidity

    // NT2018 use DRIFT via GetXiD
        for(int i=1; i<nr-1; ++i){ //loop on radius
            double r = radius[i]; //current radius
            double A = KCoeff->Eval(p)/(r0*r0); // Z/A set?
            double B = (-(1.+GetXiD())*SolarWind->Eval(r) + 2.*KCoeff->Eval(p)/r)/r0;
            double C = (2./3.)*SolarWind->Eval(r)/(r);

            // normal equations for i(radius) = [1,nr-2]
            if(i!=1 && i!=nr-2){
                (*(PA[j]))[i-1][i-2] = A/(2.*deltaAr*deltaAr) - B/(4.*deltaAr);
                (*(PA[j]))[i-1][i-1] = -A/(deltaAr*deltaAr) - C/deltalogp;
                (*(PA[j]))[i-1][i] =  A/(2.*deltaAr*deltaAr) + B/(4.*deltaAr);                 }

            // different for the first equation
            else if(i==1){
                (*(PA[j]))[i-1][i-1] = -A/(deltaAr*deltaAr) - C/deltalogp;
                (*(PA[j]))[i-1][i] =  A/(2.*deltaAr*deltaAr) + B/(4.*deltaAr);
            }

            // and for the last one
            else{
                (*(PA[j]))[i-1][i-2] = A/(2.*deltaAr*deltaAr) - B/(4.*deltaAr);
                (*(PA[j]))[i-1][i-1] = -A/(deltaAr*deltaAr) - C/deltalogp;
            }
        }
    }

}




SolarModulation::~SolarModulation(){   //distruttore
    //cout<<"Deleting SolarMOdulation Object"<<endl;
    for(int i=0; i<nr; ++i){
        delete[] f[i];
        delete[] flux[i];
    }
    
    delete[] flux;
    delete[] f;
    
    delete[] radius;
    delete[] lnRigidity;
    delete[] Rigidity;
    delete[] KEnergy;
    delete[] KScale;
    delete[] Phi;
    delete[] XiD;
    
    
    for(int i=0; i<np-1; ++i) delete PA[i];
    delete[] PA;
    
    delete SolarWind;
    delete JLisT;
    delete JLisRig;
    delete JLis_Scaled;
    delete KCoeff;
    delete ForceField;
    
    
    delete[] xRigidity;
    delete[] xKEnergy;
    delete[] xRadius;
    delete[] xKScale;
    delete[] xPhi;
    delete[] xXiD;
    
   delete S;
    
    DataSaver.clear();
}


void SolarModulation::MakeGridsForHistograms(){
  int zz= (int)pZ;
  int aa= (int)pA;

  // ---- 1D histograms ----
  hLIS_vs_Rigidity= new TH1F(Form("hLIS_vs_Rigidity_Z%d_A%d", zz, aa), Form("hLIS_vs_Rigidity_Z%d_A%d", zz, aa), nRigidity, xRigidity);
  hLIS_vs_KEnergy= new TH1F(Form("hLIS_vs_KEnergy_Z%d_A%d", zz, aa), Form("hLIS_vs_KEnergy_Z%d_A%d", zz, aa), nKEnergy, xKEnergy);

  hFlux_vs_Rigidity= new TH1F(Form("hFlux_vs_Rigidity_Z%d_A%d", zz, aa), Form("hFlux_vs_Rigidity_Z%d_A%d", zz, aa), nRigidity, xRigidity);
  hFlux_vs_KEnergy= new TH1F(Form("hFlux_vs_KEnergy_Z%d_A%d", zz, aa), Form("hFlux_vs_KEnergy_Z%d_A%d", zz, aa), nKEnergy, xKEnergy);

  hFlux_vs_Radius= new TH1F(Form("hFlux_vs_Radius_Z%d_A%d", zz, aa), Form("hFlux_vs_Radius_Z%d_A%d", zz, aa), nRadius, xRadius);

  hFlux_vs_KScale= new TH1F(Form("hFlux_vs_KScale_Z%d_A%d", zz, aa), Form("hFlux_vs_KScale_Z%d_A%d", zz, aa), nKScale, xKScale);

  hFluxFF_vs_Rigidity= new TH1F(Form("hFluxFF_vs_Rigidity_Z%d_A%d", zz, aa), Form("hFluxFF_vs_Rigidity_Z%d_A%d", zz, aa), nRigidity, xRigidity);
  hFluxFF_vs_KEnergy= new TH1F(Form("hFluxFF_vs_KEnergy_Z%d_A%d", zz, aa), Form("hFluxFF_vs_KEnergy_Z%d_A%d", zz, aa), nKEnergy, xKEnergy);
  

  // ---- 2D histograms: flux vs kscale vs rig/ekn ----
  hFlux_vs_KScale_vs_Rigidity= new TH2F(Form("hFlux_vs_KScale_vs_Rigidity_Z%d_A%d",zz,aa), Form("hFlux_vs_KScale_vs_Rigidity_Z%d_A%d",zz,aa), nRigidity, xRigidity, nKScale, xKScale);
  hFlux_vs_KScale_vs_KEnergy= new TH2F(Form("hFlux_vs_KScale_vs_KEnergy_Z%d_A%d",zz,aa), Form("hFlux_vs_KScale_vs_KEnergy_Z%d_A%d",zz,aa), nKEnergy, xKEnergy, nKScale, xKScale);

  // LOG HISTOS
  hLn_Flux_vs_KScale_vs_Rigidity= new TH2F(Form("hLn_Flux_vs_KScale_vs_Rigidity_Z%d_A%d",zz,aa), Form("hLn_Flux_vs_KScale_vs_Rigidity_Z%d_A%d",zz,aa), nRigidity, xLnRigidity, nKScale, xLnKScale);
  hLn_Flux_vs_KScale_vs_KEnergy= new TH2F(Form("hLn_Flux_vs_KScale_vs_KEnergy_Z%d_A%d",zz,aa), Form("hLn_Flux_vs_KScale_vs_KEnergy_Z%d_A%d",zz,aa), nKEnergy, xLnKEnergy, nKScale, xLnKScale);

  
  // ---- 2D histograms: flux vs radius vs Rig/Ekn ----
  hFlux_vs_Radius_vs_Rigidity= new TH2F(Form("hFlux_vs_Radius_vs_Rigidity_Z%d_A%d",zz,aa), Form("hFlux_vs_Radius_vs_Rigidity_Z%d_A%d",zz,aa), nRigidity, xRigidity, nRadius, xRadius);
  hFlux_vs_Radius_vs_KEnergy= new TH2F(Form("hFlux_vs_Radius_vs_KEnergy_Z%d_A%d",zz,aa), Form("hFlux_vs_Radius_vs_KEnergy_Z%d_A%d",zz,aa), nKEnergy, xKEnergy, nRadius, xRadius);


  // ---- 2D histograms: fluxFF vs Phi vs Rig/Ekn ----
  hFluxFF_vs_Phi_vs_Rigidity= new TH2F(Form("hFluxFF_vs_Phi_vs_Rigidity_Z%d_A%d",zz,aa), Form("hFluxFF_vs_Phi_vs_Rigidity_Z%d_A%d",zz,aa), nRigidity, xRigidity, nPhi, xPhi);
  hFluxFF_vs_Phi_vs_KEnergy= new TH2F(Form("hFluxFF_vs_Phi_vs_KEnergy_Z%d_A%d",zz,aa), Form("hFluxFF_vs_Phi_vs_KEnergy_Z%d_A%d",zz,aa), nKEnergy, xKEnergy, nPhi, xPhi);

  hLn_FluxFF_vs_Phi_vs_Rigidity= new TH2F(Form("hLn_FluxFF_vs_Phi_vs_Rigidity_Z%d_A%d",zz,aa), Form("hLn_FluxFF_vs_Phi_vs_Rigidity_Z%d_A%d",zz,aa), nRigidity, xLnRigidity, nPhi, xPhi);
  hLn_FluxFF_vs_Phi_vs_KEnergy= new TH2F(Form("hLn_FluxFF_vs_Phi_vs_KEnergy_Z%d_A%d",zz,aa), Form("hLn_FluxFF_vs_Phi_vs_KEnergy_Z%d_A%d",zz,aa), nKEnergy, xLnKEnergy, nPhi, xPhi);


  // ---- 3D histograms: flux vs XiD vs KScale vs Rig/Ekn ----

  hFlux_vs_XiD_vs_KScale_vs_Rigidity= new TH3F(Form("hFlux_vs_XiD_vs_KScale_vs_Rigidity_Z%d_A%d",zz,aa), Form("hFlux_vs_XiD_vs_KScale_vs_Rigidity_Z%d_A%d",zz,aa), nRigidity, xRigidity, nKScale, xKScale, nXiD, xXiD);
  
  hFlux_vs_XiD_vs_KScale_vs_KEnergy= new TH3F(Form("hFlux_vs_XiD_vs_KScale_vs_KEnergy_Z%d_A%d",zz,aa), Form("hFlux_vs_XiD_vs_KScale_vs_KEnergy_Z%d_A%d",zz,aa), nKEnergy, xKEnergy, nKScale, xKScale, nXiD, xXiD);

  // LOG HISTOS
  hLn_Flux_vs_XiD_vs_KScale_vs_Rigidity= new TH3F(Form("hLn_Flux_vs_XiD_vs_KScale_vs_Rigidity_Z%d_A%d",zz,aa), Form("hLn_Flux_vs_XiD_vs_KScale_vs_Rigidity_Z%d_A%d",zz,aa), nRigidity, xLnRigidity, nKScale, xLnKScale, nXiD, xXiD);
  
  hLn_Flux_vs_XiD_vs_KScale_vs_KEnergy= new TH3F(Form("hLn_Flux_vs_XiD_vs_KScale_vs_KEnergy_Z%d_A%d",zz,aa), Form("hLn_Flux_vs_XiD_vs_KScale_vs_KEnergy_Z%d_A%d",zz,aa), nKEnergy, xLnKEnergy, nKScale, xLnKScale, nXiD, xXiD);
  
  IsGrid=true;
}


// --- set LIS from an input function | DONT USE IT ---
void SolarModulation::SetJLis(TF1* aJLisT){
  int zz= (int)pZ;
  int aa= (int)pA;

  //printf("[Run # %d] Set  LIS spectrum from input TF1\n", RunID);
  JLisT=new TF1(*aJLisT);
  JLisRig = new TF1(Form("JLisRig_Z%d_A%d",zz,aa),this,&SolarModulation::fJLisRigvalue, 1,1000000, 0,Form("SolarModulation_Z%d_A%d",zz,aa),Form("JLisRig_Z%d_A%d",zz,aa));
}



void SolarModulation::InitJLis(){
  int Z= (int)pZ;
  // int A= (int)pA; // GESTIRE 3He! TBD

  // --- Mod0 Corti et al 2016 | do nothing  ---
  if( kLISModel==0) IsLIS = true;  // parametric TF1!

  // --- Mod1-Mod3 from NT: get arrays from graph ---
  // -- Mod4 from REINA
  if( kLISModel>0){
    TGraph* grLISFluxVSEkn;

    TFile* inLIS;
    // Mod1 NT et al ApJ 2017
    if( kLISModel==1) inLIS= new TFile(PATH+"/input/grLIS_pHe_Tomassetti_ApJ_2017.root","READ");

    // Mod2 NT PRD 2017
    if( kLISModel==2) inLIS= new TFile(PATH+"/input/grLIS_pHe_Tomassetti_PRD_2017.root","READ");

    // Mod3 NT ASR 2017
    if( kLISModel==3) inLIS= new TFile(PATH+"/input/grLIS_pHe_Tomassetti_ASR_2017.root","READ");

    // Mod4 REINA 2020 [p-He-BCO | RIGIDITY LIS!]
    if( kLISModel==4) inLIS= new TFile(PATH+"/input/grLIS_pHeBCO_Reina_2020.root","READ");
      
      // Mod5 BOSCHINI 2017 [p-He | RIGIDITY LIS!]
    if( kLISModel==5) inLIS= new TFile(PATH+"/input/grLIS_pHe_Boschini_ApJ_2017.root","READ");
      


    // --- get arrays ---  [good for Mod0 to Mod3. For Mod4 REINA2020, LIS are VS rigidity!]
    inLIS->cd();
   // cout<<"CONversion 1 GV ---> "<<pow(RigToEkn(1.e3),-1)<<"  MeV"<<endl;
    if (kLISModel < 4) {
      if(Z==1)grLISFluxVSEkn= (TGraph*)(inLIS->Get("grLISFluxVSEkn_Proton"));
      if(Z==2)grLISFluxVSEkn= (TGraph*)(inLIS->Get("grLISFluxVSEkn_Helium"));
      inLIS->Close();
      nLIS= grLISFluxVSEkn->GetN();
      xLnEkn= new double[nLIS];
      xLnLIS= new double[nLIS];
      for(int ee=0;ee<nLIS;ee++){
        xLnEkn[ee] = log( (1.e+3)*grLISFluxVSEkn->GetX()[ee] ); // MeV units
        xLnLIS[ee] = log( (1.e-3)*grLISFluxVSEkn->GetY()[ee] ); // MeV units
     
      }

    }
      
      
      if (kLISModel == 4)  {
      //in realta dip da R
      if(Z==1)grLISFluxVSEkn= (TGraph*)(inLIS->Get("LIS_R_Z01"));
      if(Z==2)grLISFluxVSEkn= (TGraph*)(inLIS->Get("LIS_R_Z02"));
      if(Z==5)grLISFluxVSEkn= (TGraph*)(inLIS->Get("LIS_R_Z05"));
      if(Z==6)grLISFluxVSEkn= (TGraph*)(inLIS->Get("LIS_R_Z06"));
      if(Z==8)grLISFluxVSEkn= (TGraph*)(inLIS->Get("LIS_R_Z08"));
      inLIS->Close();
      nLIS= grLISFluxVSEkn->GetN();
      xLnEkn= new double[nLIS];
      xLnLIS= new double[nLIS];
      for(int ee=0;ee<nLIS;ee++){
            xLnLIS[ee]= log( (1.e-3)*grLISFluxVSEkn->GetY()[ee]*dRdE((RigToEkn((1.e+3)*grLISFluxVSEkn->GetX()[ee]))));
            xLnEkn[ee] = log( RigToEkn((1.e+3)*grLISFluxVSEkn->GetX()[ee]) );
        }
      }
      
      if (kLISModel == 5)  {
      //in realta dip da R. Trasforma il graph di input da Rig -> E
      if(Z==1)grLISFluxVSEkn= (TGraph*)(inLIS->Get("grLISFluxVSRig_Proton"));
      if(Z==2)grLISFluxVSEkn= (TGraph*)(inLIS->Get("grLISFluxVSRig_Helium"));

      inLIS->Close();
      nLIS= grLISFluxVSEkn->GetN();
      xLnEkn= new double[nLIS];
      xLnLIS= new double[nLIS];
    //flusso Rig -> Flusso Ekn (ne prendo il log)
    // * Jac (dRdE) e conversione E -> (MeV)
    // flusso in input in GV
    
          for(int ee=0;ee<nLIS;ee++){
              xLnLIS[ee]= log( (1.e-3)*grLISFluxVSEkn->GetY()[ee]*dRdE((RigToEkn((1.e+3)*grLISFluxVSEkn->GetX()[ee]))));
              xLnEkn[ee] = log( RigToEkn((1.e+3)*grLISFluxVSEkn->GetX()[ee]) );
          }
    }



    IsLIS = true;
    delete grLISFluxVSEkn;
    delete inLIS;

  }
}



// --- set default LIS based on Z and A ---
void SolarModulation::SetJLis(){
  int zz= (int)pZ;
  int aa= (int)pA;
  //printf("Set LIS spectrum for nucleus [Z=%d | A=%d].  Use LIS Model N.%d\n",zz,aa,kLISModel);

  // ---- LIS MODELS ----
  // 0: Corti et al. 2016 ApJ [DEFAULT]
  // 1: Tomassetti et al. 2017 ApJ [OK]
  // 2: Tomassetti et al. 2017 PRD []
  // 3: Tomassetti et al. 2017 ASR
  
  
  // initialize LIS to given model
  if(!IsLIS) InitJLis();
  
  // ---- Set parametric LIS from Corti et al. ApJ 2016 ----
  
  // --- set LIS vs ekn ---
  TF1* JLISvsEKN = new TF1(Form("fLIS_vs_EKN_Z%d_A%d",zz,aa),this,&SolarModulation::functorLISvsEKN,
             emin, emax, 0, Form("fLIS_vs_EKN_Z%d_A%d",zz,aa), Form("fLIS_vs_EKN_Z%d_A%d",zz,aa));


  JLisT=new TF1(*JLISvsEKN); // prende anche il NOME?
  JLisT->SetNpx(1000);
  delete JLISvsEKN;

  
  // --- setLIS vs rig ---
  TF1* JLISvsRIG = new TF1(Form("fLIS_vs_RIG_Z%d_A%d",zz,aa),this,&SolarModulation::functorLISvsRIG,
             pmin, pmax, 0, Form("fLIS_vs_RIG_Z%d_A%d",zz,aa), Form("fLIS_vs_RIG_Z%d_A%d",zz,aa));

  JLisRig= new TF1(*JLISvsRIG);
  JLisRig->SetNpx(1000);
  delete JLISvsRIG;

}



// --- force field solver  ---
void SolarModulation::ForceFieldSolution(double phi){
  ForceField = new TF1(Form("fFF_Z%d_A%d",(int)pZ,(int)pA),this,&SolarModulation::functorFF, emin, emax, 1,"SM","functorFF");
//1 indica il numero di parametri da passare
  ForceField->SetParameter(0,phi);  //questo è par[0]
}

void SolarModulation::ForceFieldSolution(){
  if(!IsForceField){
    ForceField = new TF1(Form("fFF_Z%d_A%d",(int)pZ,(int)pA),this,&SolarModulation::functorFF, emin, emax, 1,"SM","functorFF");
  }

  ForceField->SetParameter(0,ThisPhi); // Phi from KScale ( i.e. par[0] )
  ForceField->SetLineColor(kGreen+1);
  ForceField->SetLineWidth(3);
  IsForceField= true;
}


// --- NT set diffusion coefficient TF1 ---
void SolarModulation::SetKCoeff(){
    //1: è il numero di parametri
  KCoeff = new TF1(Form("fDC_Z%d_A%d",(int)pZ,(int)pA),this,&SolarModulation::functorDC, pmin, pmax, 1,"DC","functorDC");
  //if(indexKScale<0) KCoeff->SetParameter(0, 6.0);
  //if(indexKScale>-1) KCoeff->SetParameter(0, ThisKScale2);  // useless
    KCoeff->SetParameter(0, ThisKScale2);
}

// --- DAVID set diffusion coefficient TF1 ---
void SolarModulation::SetKCoeffPower(){
    //1: è il numero di parametri
  KCoeff = new TF1(Form("fDC_Z%d_A%d",(int)pZ,(int)pA),this,&SolarModulation::functorDCPower, pmin, pmax, 2,"DC","functorDC");
  //if(indexKScale<0) KCoeff->SetParameter(0, 6.0);
  //if(indexKScale>-1) KCoeff->SetParameter(0, ThisKScale2);  // useless
    KCoeff->SetParameter(0, ThisKScale2);
    cout<<"thiskscale   "<<ThisKScale2<<endl;
    KCoeff->SetParameter(1, index);
}



//broken line model for DC
void SolarModulation::SetKCoeffP(){
    //1: è il numero di parametri
  KCoeff = new TF1(Form("fDC_Z%d_A%d",(int)pZ,(int)pA),this,&SolarModulation::functorDCP, pmin, pmax, 1,"DC","functorDC");
  //if(indexKScale<0) KCoeff->SetParameter(0, 6.0);
  //if(indexKScale>-1) KCoeff->SetParameter(0, ThisKScale2);  // useless
    KCoeff->SetParameter(0, ThisKScale2);
  //KCoeff->SetParameter(1, indexAlfa);
  //KCoeff->SetParameter(2, indexBeta);
}


//broken line model for DC rk free
void SolarModulation::SetKCoeffP2(){
    //1: è il numero di parametri
  KCoeff = new TF1(Form("fDC_Z%d_A%d",(int)pZ,(int)pA),this,&SolarModulation::functorDCP2, pmin, pmax, 1,"DC","functorDC");
  //if(indexKScale<0) KCoeff->SetParameter(0, 6.0);
  //if(indexKScale>-1) KCoeff->SetParameter(0, ThisKScale2);  // useless
    KCoeff->SetParameter(0, ThisKScale2);
  //KCoeff->SetParameter(1, indexAlfa);
  //KCoeff->SetParameter(2, indexBeta);
}


void SolarModulation::SetSolarWind(){
  SolarWind = new TF1(Form("fSW_Z%d_A%d",(int)pZ,(int)pA),this,&SolarModulation::functorSW, rmin, rmax, 1,"SW","functorSW");
  SolarWind->SetParameter(0, -1); // NEG -> Set Vw x 1.e-3
}

//solve LIS part (comune per tutti i K0 ), boundary conditions
void SolarModulation::Solve_BC(){
  // should we put a pA/pZ-factor?
  double NTF= 1.; // pA/pZ; NO

  // set boundary conditions
  double Tmax = RigToEkn(pmax);

  for(int i=0; i<np; ++i){ //loop on rigidity p
    double p = exp(lnRigidity[i]);
    double T = RigToEkn(p);
    f[0][i] = 0.;
    f[nr-1][i] = (JLisT->Eval(T)/(p*p*NTF))/f0;   //valori dei nodi della griglia psi al contorno
    // division by f0 is for adimensionalisation
    // NT Should an pA/pZ factor be here?
  }
    
  for(int i=0; i<nr; ++i) f[i][np-1] = (JLisT->Eval(Tmax)/(pmax*pmax*NTF))/f0;
    //valori dei nodi della griglia psi al contorno
    //f è la psi,
  // NT Should an pA/pZ factor be here?
}




// --- numerical solver ---
void SolarModulation::Solve(){
    double NTF= 1.;
    
  //--- Generate the solution ---
    for(int j=np-2; j>=0; --j){ //loop on rigidity j
        
        
        double p = exp(lnRigidity[j]); //current rigidity p [MV]
        double T = RigToEkn(p); // current kinetic energy per nucleon
        
        // creating aux
        double aux[nr-2];
        
        // setting coeficients | NT2018 use DRIFT via GetXiD
        for(int i=1; i<nr-1; ++i){ //loop on radius i
            double r = radius[i]; //current radius [m]
            double A = KCoeff->Eval(p)/(r0*r0);   //phi_1
            double B = (-(1.+GetXiD())*SolarWind->Eval(r) + 2.*KCoeff->Eval(p)/r)/r0; //phi_2
            double C = (2./3.)*SolarWind->Eval(r)/(r); //phi_3
            
            if(i!=nr-2) aux[i-1] = f[i-1][j+1]*(-A/(2.*deltaAr*deltaAr)+B/(4.*deltaAr)) + f[i][j+1]*(A/(deltaAr*deltaAr)-C/deltalogp) + f[i+1][j+1]*(-A/(2.*deltaAr*deltaAr)-B/(4.*deltaAr));
            else aux[i-1] = f[i-1][j+1]*(-A/(2.*deltaAr*deltaAr)+B/(4.*deltaAr)) + f[i][j+1]*(A/(deltaAr*deltaAr)-C/deltalogp) + f[i+1][j+1]*(-A/(2.*deltaAr*deltaAr)-B/(4.*deltaAr)) - f[nr-1][j]*(A/(2.*deltaAr*deltaAr)+B/(4.*deltaAr));
            
            
#ifdef DEBUG
            printf(" j = %d i = %d |  A = %lf , B = %lf, C = %lf  AUX = %lf\n",j,i,A,B,C,aux[i]);
#endif
            
        }
        // creating objects for solving the linear sistem
       // Vec B(nr-2, aux);

        B.Set(nr-2, aux);
        
        S->SetConstants( *(PA[j])  ,B);
        
         sol = S->RecursiveTriDiagonal();
        // filling the solution matrix
        for(int k=1; k<nr-1; ++k) f[k][j] = sol[k-1];
        
   }

    
  // --- get flux from distribution function's solution f ---
  for(int i=0; i<nr; ++i){
    for(int j=0; j<np; ++j){
      double p = exp(lnRigidity[j]);
      flux[i][j] = f[i][j]*(p*p*NTF)*f0; //NT: pA/pZ here?
    }
  }
}




void SolarModulation::SetSolutions(){
  if(!IsGrid) MakeGridsForHistograms();
  int zz= (int)pZ;
  int aa= (int)pA;

  // --- set LIS histograms ----
  for(int pp=0;pp<nKEnergy;pp++){
    double pFluxLISvsEkn=  JLisT->Eval( KEnergy[pp] );
    hLIS_vs_KEnergy->SetBinContent(pp+1, pFluxLISvsEkn);
  }

  for(int pp=0;pp<nRigidity;pp++){
    double pFluxLISvsRig=  JLisRig->Eval( Rigidity[pp] );
    hLIS_vs_Rigidity->SetBinContent(pp+1, pFluxLISvsRig);
  }


  // --- best radius index for local flux at R0=1AU ---
  int indR0 = -1; // 1AU should be = mult = 5
  double dist=1.e+9;
  for(int rr=0;rr<nRadius;rr++){
    double ddist = fabs( radius[rr]/AU - 1. );
    if(ddist<dist){
      dist=ddist;
      indR0= rr;
    }
  }
 
cout<<"Reference flux at radius R="<<radius[indR0]/AU<<" AU"<<" is at index "<<indR0<<endl;

  // --- best energy index for flux in heliosphere --- [~1 GeV/n]
  int indE0 = -1; // 1AU should be = mult = 5
  dist=1.e+20;
  for(int ee=0;ee<nKEnergy;ee++){
    double ddist = fabs( KEnergy[ee] - 1000. ); //MeV/n!
    if(ddist<dist){
      dist=ddist;
      indE0= ee;
    }
  }
  cout<<"Reference flux at energy E="<<KEnergy[indE0]<<" MeV/n"<<endl;

  // --- best rigidity index for flux in heliosphere --- [~1 GV]
  int indP0 = -1; // 1AU should be = mult = 5
  dist=1.e+20;
  for(int pp=0;pp<nRigidity;pp++){
    double ddist = fabs( Rigidity[pp] - 1000. ); // MV!
    if(ddist<dist){
      dist=ddist;
      indP0= pp;
    }
  }
  cout<<"Reference flux at rigidity P="<<Rigidity[indP0]<<" MV"<<endl;


  // --- kscale index for calculations ----
  int indK0= indexKScale; // from RunID;
  cout<<"Current KScale="<< ThisKScale <<" is at index: "<<indexKScale<<endl;

  // --- drift index for calculations ----
  int indXiD= indexXiD; // from RunID;
  cout<<"Current Xi-Drift="<< ThisXiD <<" is at index: "<<indexXiD<<endl;
  
    
    cout<<endl;
    cout<<"Model for Diffusion Coeff K(R)"<<endl;
    if (indexKmodel != 1) {
        cout<<"Model Potgeiter 2013"<<endl;
    }
    else{cout<<"Standard Model"<<endl;}
    
    cout<<endl;
    
    
  // --- phi index for binning ---
  int indPHI = -1; // 1AU should be = mult = 5
  dist=1.e+9;
  for(int rr=0;rr<nPhi;rr++){
    double ddist = fabs( Phi[rr] - ThisPhi );
    if(ddist<dist){
      dist=ddist;
      indPHI= rr;
    }
  }
  indexPhi= indPHI;
  cout<<"Current FF modulation Phi="<<Phi[indPHI]<<" is at index: "<<indPHI<<endl;


  
  // --- near-Earth flux vs Ekn | at R=1 AU ---
  for(int ee=0;ee<nKEnergy;ee++){ // set flux vs ekn
    double pFluxEkn= flux[indR0][ee];
    hFlux_vs_KEnergy->SetBinContent(ee+1, pFluxEkn);
      //cout<<pFluxEkn<<endl;
  }

  
  // --- near-Earth flux vs Rig | at R=1 AU ---
  for(int pp=0;pp<nRigidity;pp++){ // set flux vs rig
    double pRigidit = Rigidity[pp];
    double pFluxEkn = flux[indR0][pp];
    double pFluxRig = pFluxEkn*dEdR(pRigidit);
    hFlux_vs_Rigidity->SetBinContent(pp+1, pFluxRig);
  }


  // --- force-field flux vs Ekn ----
  if(!IsForceField) ForceFieldSolution();

  for(int ee=0;ee<nKEnergy;ee++){ // set flux vs ekn
    double pFluxEkn= ForceField->Eval(KEnergy[ee]);
    hFluxFF_vs_KEnergy->SetBinContent(ee+1, pFluxEkn);
    hFluxFF_vs_Phi_vs_KEnergy->SetBinContent(ee+1, indPHI+1, pFluxEkn);
    hLn_FluxFF_vs_Phi_vs_KEnergy->SetBinContent(ee+1, indPHI+1, log(pFluxEkn));
  }

  // --- force-field flux vs Rig ---
  for(int pp=0;pp<nRigidity;pp++){ // set flux vs rig
    double pRigidit = Rigidity[pp];
    double pFluxEkn = ForceField->Eval(KEnergy[pp]);
    double pFluxRig = pFluxEkn*dEdR(pRigidit);
    hFluxFF_vs_Rigidity->SetBinContent(pp+1, pFluxRig);
    hFluxFF_vs_Phi_vs_Rigidity->SetBinContent(pp+1, indPHI+1, pFluxRig);
    hLn_FluxFF_vs_Phi_vs_Rigidity->SetBinContent(pp+1, indPHI+1, log(pFluxRig));
  }


  // --- flux vs radius | at fixed energy 1 GeV ---
  for(int rr=0;rr<nRadius;rr++){
    double pFluxEkn = flux[rr][indE0]; // ekn flux
    hFlux_vs_Radius->SetBinContent(rr+1, pFluxEkn);
  }

  // --- flux vs kscale | at fixed E= 1 GeV and R=1 AU ---
  double pFluxEkn = flux[indR0][indE0]; // ekn flux
  hFlux_vs_KScale->SetBinContent(indK0+1, pFluxEkn); // ONE POINT PER RUN!
  
  

  // --- local flux vs kcoeff vs ekn | at 1 AU ---
  if(indK0>=0 && indK0<nKScale){

    // local flux
    for(int ee=0;ee<nKEnergy;ee++){
      double pFluxEkn= flux[indR0][ee];
      hFlux_vs_KScale_vs_KEnergy->SetBinContent(ee+1, indK0+1, pFluxEkn);
      hLn_Flux_vs_KScale_vs_KEnergy->SetBinContent(ee+1, indK0+1, 1.*log(pFluxEkn));
    }

    for(int pp=0;pp<nRigidity;pp++){
      double pRigidit = Rigidity[pp];
      double pFluxEkn = flux[indR0][pp];
      double pFluxRig = pFluxEkn*dEdR(pRigidit);
      hFlux_vs_KScale_vs_Rigidity->SetBinContent(pp+1, indK0+1, pFluxRig);
      hLn_Flux_vs_KScale_vs_Rigidity->SetBinContent(pp+1, indK0+1, 1.*log(pFluxRig));
    }

  }
    
    
  // ---- flux vs radius vs energy ----
  for(int rr=0;rr<nRadius;rr++){
    for(int ee=0;ee<nKEnergy;ee++){
      double pFluxEkn = flux[rr][ee];
      hFlux_vs_Radius_vs_KEnergy->SetBinContent(ee+1, rr+1, pFluxEkn);
    }
  }
  
  // ---- flux vs radius vs rigidity ----
  for(int rr=0;rr<nRadius;rr++){
    for(int pp=0;pp<nRigidity;pp++){
      double pRigidit = Rigidity[pp];
      double pFluxEkn = flux[rr][pp];
      double pFluxRig = pFluxEkn*dEdR(pRigidit);
      hFlux_vs_Radius_vs_Rigidity->SetBinContent(pp+1, rr+1, pFluxRig);
    }
  }


  
  // ---- DRIFT! local flux vs xi-drift vs kscale vs energy | at r=1AU ----

  if(indK0>=0 && indK0<nKScale && indXiD>=0 && indXiD<nXiD){

    for(int ee=0;ee<nKEnergy;ee++){
      double pFluxEkn= flux[indR0][ee];
      hFlux_vs_XiD_vs_KScale_vs_KEnergy->SetBinContent(ee+1, indK0+1, indXiD+1, pFluxEkn);
      hLn_Flux_vs_XiD_vs_KScale_vs_KEnergy->SetBinContent(ee+1, indK0+1, indXiD+1, log(pFluxEkn));
    }

    for(int pp=0;pp<nRigidity;pp++){
      double pRigidit = Rigidity[pp];
      double pFluxEkn = flux[indR0][pp];
      double pFluxRig = pFluxEkn*dEdR(pRigidit);
      hFlux_vs_XiD_vs_KScale_vs_Rigidity->SetBinContent(pp+1, indK0+1, indXiD+1, pFluxRig);
      hLn_Flux_vs_XiD_vs_KScale_vs_Rigidity->SetBinContent(pp+1, indK0+1, indXiD+1, log(pFluxRig));
    }

  }
    
  
  // ---- full solutions are set ----
  IsSolution= true;
}



void SolarModulation::StoreResults(){

  int zz= (int)pZ;
  int aa= (int)pA;
  int kk= indexKScale;
  int dd= indexXiD;
  
  Printf("Store results for particle [Z%d, A%d], Run_%d with indK0_%d and indXiD_%d \n",zz,aa,RunID,kk,dd);

  // store locally
  TFile* outFile;
  if(DRIFT==0){
    // outFile= new TFile(Form("$NTBASEDIR/THEORY/00matisse/PHeModel/outroot/hFluxResults_LIS%d_Z%d_A%d_RunID%d.root",kLISModel,zz,aa,RunID),"RECREATE");
    // outFile= new TFile(Form("./hrun_FluxResults_LIS%d_Z%d_A%d_RunID%d.root",kLISModel,zz,aa,RunID),"RECREATE");
    outFile= new TFile(Form(PATH+"/OUT/Model_DriftOFF/histo_FluxResults_LIS%d_Z%d_A%d.root",kLISModel,zz,aa),"RECREATE");
  }
  
  if(DRIFT==1){ // kk and dd indices
    // outFile= new TFile(Form("$NTBASEDIR/THEORY/00matisse/PHeModel/outroot/hFluxResults_LIS%d_Z%d_A%d_XiD%d_K%d_RunID%d.root",kLISModel,zz,aa,dd,kk,RunID),"RECREATE");
    // outFile= new TFile(Form("./hrun_FluxResults_LIS%d_Z%d_A%d_XiD%d_K%d_RunID%d.root",kLISModel,zz,aa,kk,dd,RunID),"RECREATE");
    outFile= new TFile(Form(PATH+"/OUT/Model_DriftON/hrun_FluxResults_LIS%d_Z%d_A%d_XiD%d_K%d_RunID%d.root",zz,aa,kLISModel,zz,aa,dd,kk,RunID),"RECREATE");
  }
  
  outFile->cd();

  // --- save 1D-histos: never ----

  // --- save 2D-histos: for zero-drift only ----
  if(DRIFT==0 || indexXiD==49 || ThisXiD==0){
    hFlux_vs_KScale_vs_Rigidity->Write();
    hFlux_vs_KScale_vs_KEnergy->Write();
    //hFluxFF_vs_Phi_vs_Rigidity->Write();
    //hFluxFF_vs_Phi_vs_KEnergy->Write();
    //hFlux_vs_Radius->Write();
    
    hLn_Flux_vs_KScale_vs_Rigidity->Write();
    hLn_Flux_vs_KScale_vs_KEnergy->Write();
    //hLn_FluxFF_vs_Phi_vs_Rigidity->Write();
    //hLn_FluxFF_vs_Phi_vs_KEnergy->Write();
  }

  // --- save 3D-histos: ONLY for DRIFT=1 ---
  if(DRIFT==1){
    hFlux_vs_XiD_vs_KScale_vs_Rigidity->Write();
    hFlux_vs_XiD_vs_KScale_vs_KEnergy->Write();
    hLn_Flux_vs_XiD_vs_KScale_vs_Rigidity->Write();
    hLn_Flux_vs_XiD_vs_KScale_vs_KEnergy->Write();
  }
  
  outFile->Write();
  outFile->Close();
  
  delete outFile;

  // Printf("Program executed ¯\\_(ツ)_/¯\n\n",zz,aa,RunID,kk,dd);
  Printf("Program executed ¯\\_(ツ)_/¯\n\n");

}



void SolarModulation::PlotSolution(){
  int zz= (int)pZ;
  int aa= (int)pA;
  TApplication* theApp = new TApplication(Form("App_Z%d_A%d",zz,aa), 0, 0);

  // ---- get Voyager-1 data ----


  // ---- Voyager-1 Hydrogen flux ---- A. C. Cumming et al. 2016, ApJ 831, 18
  double Voyager1_Hydrogen_Ekn[15]={3.8,5.4,6.9,10.2,15.4,23.9,39.0,52.0,79.1,143.9,164.9,181.3,204.1,245.3,308.0};
  double Voyager1_Hydrogen_Flux[15]={1.908e+01,2.116e+01,2.369e+01,2.665e+01,2.968e+01,2.909e+01,2.829e+01,2.914e+01,2.443e+01,1.927e+01,1.820e+01,1.652e+01,1.492e+01,1.289e+01,1.042e+01};
  double stat_Voyager1_Hydrogen_Flux[15]={9.572e-01,1.061e+00,1.188e+00,1.345e+00,1.496e+00,1.460e+00,1.419e+00,1.468e+00,1.263e+00,9.723e-01,9.197e-01,8.410e-01,7.519e-01,6.484e-01,5.236e-01};
  double syst_Voyager1_Hydrogen_Flux[15]={9.540e-01,1.058e+00,1.184e+00,1.333e+00,1.484e+00,1.454e+00,1.415e+00,1.457e+00,1.221e+00,9.635e-01,9.100e-01,8.260e-01,7.460e-01,6.445e-01,5.210e-01};

  double eVoyager1_Hydrogen_Ekn[15];
  double eVoyager1_Hydrogen_Flux[15];
  for(int ee=0;ee<15;ee++){
    eVoyager1_Hydrogen_Ekn[ee] = 0.;
    eVoyager1_Hydrogen_Flux[ee] = sqrt(stat_Voyager1_Hydrogen_Flux[ee]*stat_Voyager1_Hydrogen_Flux[ee] + syst_Voyager1_Hydrogen_Flux[ee]*syst_Voyager1_Hydrogen_Flux[ee]);
  }

  TGraphErrors* grVOYAGER1_Hydrogen_Flux_2016 = new TGraphErrors(15, Voyager1_Hydrogen_Ekn, Voyager1_Hydrogen_Flux,
                                 eVoyager1_Hydrogen_Ekn, eVoyager1_Hydrogen_Flux);
  grVOYAGER1_Hydrogen_Flux_2016->SetMarkerStyle(25);
  grVOYAGER1_Hydrogen_Flux_2016->SetMarkerSize(1);
  grVOYAGER1_Hydrogen_Flux_2016->SetMarkerColor(kBlack);
  grVOYAGER1_Hydrogen_Flux_2016->SetLineColor(kBlack);


  // ---- Voyager-1 Helium flux ---- A. C. Cumming et al. 2016, ApJ 831, 18
// sono i dati dentro in file Voyager in data
  double Voyager1_Helium_Ekn[16]={3.8,5.4,6.9,10.2,15.4,23.9,39.0,52.0,122.0,142.8,172.4,216.0,278.3,349.5,430.9,571.0};
  double Voyager1_Helium_Flux[16]={1.530e+00,1.767e+00,1.936e+00,2.228e+00,2.475e+00,2.563e+00,2.370e+00,2.241e+00,1.746e+00,1.585e+00,1.409e+00,1.184e+00,9.229e-01,7.844e-01,5.954e-01,4.235e-01};
  double stat_Voyager1_Helium_Flux[16]={7.964e-02,9.149e-02,1.002e-01,1.233e-01,1.358e-01,1.336e-01,1.226e-01,1.226e-01,8.831e-02,8.059e-02,7.103e-02,5.973e-02,4.645e-02,3.960e-02,3.001e-02,2.131e-02};
  double syst_Voyager1_Helium_Flux[16]={7.650e-02,8.835e-02,9.680e-02,1.114e-01,1.238e-01,1.282e-01,1.185e-01,1.120e-01,8.730e-02,7.925e-02,7.045e-02,5.920e-02,4.614e-02,3.922e-02,2.977e-02,2.118e-02};


  double eVoyager1_Helium_Ekn[16];
  double eVoyager1_Helium_Flux[16];
  for(int ee=0;ee<16;ee++){
    eVoyager1_Helium_Ekn[ee] = 0.;
    eVoyager1_Helium_Flux[ee] = sqrt(stat_Voyager1_Helium_Flux[ee]*stat_Voyager1_Helium_Flux[ee] + syst_Voyager1_Helium_Flux[ee]*syst_Voyager1_Helium_Flux[ee]);
  }


    TGraphErrors* grVOYAGER1_Helium_Flux_2016 = new TGraphErrors(16, Voyager1_Helium_Ekn, Voyager1_Helium_Flux,
                                 eVoyager1_Helium_Ekn, eVoyager1_Helium_Flux);
    grVOYAGER1_Helium_Flux_2016->SetMarkerStyle(25);
    grVOYAGER1_Helium_Flux_2016->SetMarkerSize(1);
    grVOYAGER1_Helium_Flux_2016->SetMarkerColor(kBlack);
    grVOYAGER1_Helium_Flux_2016->SetLineColor(kBlack);


    //Voyager1 Helium3 data (sono identici ai dati importati a mano sopra)

     //Helium AMS
    // ---- get AMS02 data vs Ekn ----
    TFile* inFileHelium_AMS= new TFile(PATH+"/data/grAMS02_Helium_Ekn.root","READ");
    inFileHelium_AMS->cd();
    TGraphAsymmErrors* grAMS02_Helium_Flux= (TGraphAsymmErrors*)inFileHelium_AMS->Get("grAMS02_Helium_FluxVSEkn");
    grAMS02_Helium_Flux->SetMarkerStyle(20);
    grAMS02_Helium_Flux->SetMarkerColor(kBlack);
    grAMS02_Helium_Flux->SetMarkerSize(0.8);
    inFileHelium_AMS->Close();

    // --- GeV->MeV conversion ---
    for(int ee=0;ee<grAMS02_Helium_Flux->GetN();ee++){
        grAMS02_Helium_Flux->GetX()[ee] *= 1.e+3;
        grAMS02_Helium_Flux->GetY()[ee] *= 1.e-3;
        //grAMS02_Helium_Flux->GetEXlow()[ee] *= 1.e+3;
        //grAMS02_Helium_Flux->GetEXhigh()[ee] *= 1.e+3;
        grAMS02_Helium_Flux->GetEY()[ee] *= 1.e-3;
        //grAMS02_Helium_Flux->GetEYhigh()[ee] *= 1.e-3;
    }
    
    
    // ---- get AMS02 data vs Ekn ----
    TFile* inFileProton= new TFile(PATH+"/data/grProtonData.root","READ");
    inFileProton->cd();
    TGraphAsymmErrors* grAMS02_Proton_Flux_2015= (TGraphAsymmErrors*)inFileProton->Get("gr_exp1");
    grAMS02_Proton_Flux_2015->SetMarkerStyle(20);
    grAMS02_Proton_Flux_2015->SetMarkerColor(kBlack);
    grAMS02_Proton_Flux_2015->SetMarkerSize(0.8);
    inFileProton->Close();

    // --- GeV->MeV conversion ---
    for(int ee=0;ee<grAMS02_Proton_Flux_2015->GetN();ee++){
      grAMS02_Proton_Flux_2015->GetX()[ee] *= 1.e+3;
      grAMS02_Proton_Flux_2015->GetY()[ee] *= 1.e-3;
      grAMS02_Proton_Flux_2015->GetEXlow()[ee] *= 1.e+3;
      grAMS02_Proton_Flux_2015->GetEXhigh()[ee] *= 1.e+3;
      grAMS02_Proton_Flux_2015->GetEYlow()[ee] *= 1.e-3;
      grAMS02_Proton_Flux_2015->GetEYhigh()[ee] *= 1.e-3;
    }

    // ---- graph with solution ----
    // 5 è indR0, non posso mettere una variabile
    //TGraph* grFluxSolution= new TGraph(np, KEnergy, flux[5]);
    //ricalcolo indr0;

    int indR0 = -1; // 1AU should be = mult = 5
    double dist=1.e+9;
    for(int rr=0;rr<nRadius;rr++){
      double ddist = fabs( radius[rr]/AU - 1. );
      if(ddist<dist){
        dist=ddist;
        indR0= rr;
      }
    }

    //cout<<"distance in m "<<radius[indR0]/AU <<endl;
    //cout<<"distance in m "<<dist*AU<<endl;
    //cout<<"index "<<indR0<<endl;
    //cout<<"Ekin-min "<<emin<<"  Ekin_max  "<<emax<<endl;
    //cout<<"Ekin-min "<<KEnergy[0]<<"  Ekin_max  "<<KEnergy[499]<<endl;

    TGraph* grFluxSolution= new TGraph(np, KEnergy, flux[indR0]);

    //TGraph* grFluxSolution= new TGraph(np, KEnergy, flux[5]);
    //se voglio calcore sempre near-earth
    //TGraph* grFluxSolution= new TGraph(np, KEnergy, flux_at_r);
    grFluxSolution->SetLineColor(kAzure+1);
    grFluxSolution->SetLineWidth(3);
      // mult= 5 -> flux[mult][pp]
    // for (int pp=0; pp<np;++pp){
    //   double flux=flux[r[5]][pp];

    // ---- some printout ----
    // for(int rr=0;rr<nr;rr++){
    //   cout<<rr<<"    Radius: "<<radius[rr]/AU<<endl;
    // }


    // ---- frames ----
    // double EMIN=1.; // 1 MeV
    // double EMAX=300000.; // 300 GeV
    //TH2F* hFrameProtonEkn= new TH2F("hFrameProtonEkn","hFrameProtonEkn",200, EMIN, EMAX, 200, 1.e-7, 1.e+2); // MeV

    double EMIN=100; // 100 MeV
    double EMAX=100000.; // 100 GeV

    TH2F* hFrameProtonEkn= new TH2F("hFrameProtonEkn","hFrameProtonEkn",200, EMIN, EMAX, 200, 1.e-6, 1.e+2); // MeV
    //TH2F* hFrameProtonEkn= new TH2F("hFrameProtonEkn","hFrameProtonEkn",300, emin, emax, 300, 1.e-7, 1.e4); // MeV


    hFrameProtonEkn->GetXaxis()->SetTitle("kinetic energy [MeV/n]");
    hFrameProtonEkn->GetYaxis()->SetTitle("J(E) [MeV^{ -1} m^{ -2} s^{ -1} sr^{ -1} ]");
    hFrameProtonEkn->SetTitle(0);
    hFrameProtonEkn->GetYaxis()->SetNdivisions(504);
    hFrameProtonEkn->SetLabelFont(42,"X");
    hFrameProtonEkn->SetLabelFont(42,"Y");
    hFrameProtonEkn->SetTitleFont(42,"X");
    hFrameProtonEkn->SetTitleFont(42,"Y");
    hFrameProtonEkn->GetXaxis()->SetLabelSize(0.05);
    hFrameProtonEkn->GetYaxis()->SetLabelSize(0.05);
    hFrameProtonEkn->GetXaxis()->SetTitleSize(0.05);
    hFrameProtonEkn->GetYaxis()->SetTitleSize(0.05);
    hFrameProtonEkn->GetXaxis()->SetTitleOffset(1.375);
    hFrameProtonEkn->GetYaxis()->SetTitleOffset(1.55);
    hFrameProtonEkn->SetStats(0);

    TCanvas* ccFluxSolution = new TCanvas(Form("ccFluxSolution_Z%d_A%d",zz,aa), Form("SOLUTION for Z=%d A=%d Model LIS=%d",zz,aa,kLISModel), 600, 600); //1400, 1000);
    ccFluxSolution->cd();
    gPad->SetBottomMargin(0.18);
    gPad->SetLeftMargin(0.18);
    gPad->SetRightMargin(0.05);
    gPad->SetTopMargin(0.05);
    
    hFrameProtonEkn->Draw(); //background th2f
    
    gPad->SetLogx();
    gPad->SetLogy();

   grVOYAGER1_Hydrogen_Flux_2016->Draw("pZ");
//grVOYAGER1_Helium_Flux_2016->Draw("pZ");
    
    grAMS02_Proton_Flux_2015->Draw("pZ");
  //   grAMS02_Helium_Flux->Draw("pZ"); // METTERE...?

    JLisT->SetLineWidth(3);
    grFluxSolution->SetLineWidth(3);

    JLisT->Draw("l same");
    ForceField->Draw("l same");
    grFluxSolution->Draw("l same");

    TCanvas* ccKCoeff = new TCanvas(Form("ccKCoeff"), Form("K vs Rigidity"), 600, 600); //1400, 1000);
    ccKCoeff->cd();
    TH2F* hFrameK= new TH2F("hFrameK","hFrameK",200, pmin, pmax, 200, KCoeff->Eval(pmin), KCoeff->Eval(pmax)); // MeV
    //TH2F* hFrameProtonEkn= new TH2F("hFrameProtonEkn","hFrameProtonEkn",300, emin, emax, 300, 1.e-7, 1.e4); // MeV

    hFrameK->GetXaxis()->SetTitle("Rigidity[MV]");
    hFrameK->GetYaxis()->SetTitle("K [m^{2} s^{-1}]");
    hFrameK->SetTitle(0);
    hFrameK->GetYaxis()->SetNdivisions(504);
    hFrameK->SetLabelFont(42,"X");
    hFrameK->SetLabelFont(42,"Y");
    hFrameK->SetTitleFont(42,"X");
    hFrameK->SetTitleFont(42,"Y");
    hFrameK->GetXaxis()->SetLabelSize(0.05);
    hFrameK->GetYaxis()->SetLabelSize(0.05);
    hFrameK->GetXaxis()->SetTitleSize(0.05);
    hFrameK->GetYaxis()->SetTitleSize(0.05);
    hFrameK->GetXaxis()->SetTitleOffset(1.375);
    hFrameK->GetYaxis()->SetTitleOffset(1.55);
    hFrameK->SetStats(0);
    gPad->SetBottomMargin(0.18);
    gPad->SetLeftMargin(0.18);
    gPad->SetRightMargin(0.05);
    gPad->SetTopMargin(0.05);
    hFrameK->Draw();
    gPad->SetLogx();
    gPad->SetLogy();
    KCoeff->SetLineWidth(3);
    KCoeff->Draw("l same");
    // store locally
    TFile* outFile;
    outFile= new TFile(Form(PATH+"/out/Model_DriftOFF/TGraph_Flux_data_LIS%d_Z%d_A%d.root",kLISModel,zz,aa),"recreate");
    outFile->cd();
    
    //salva grafici su file root
 
 // WARNING for TF1 saving to root file
   // ccFluxSolution->Write();
    //ForceField->Write();
    /*
 Instead TF1 based on real C/C++ functions (e.g. lambda or functors), cases 4,5,6 cannot be cloned and saved in a file, because they depend on some functional code.
 It is instead possible to copy them using the copy constructor or TF1::Copy, but this will work as far the function code on which they are based is still available and valid.
 */
    grFluxSolution->Write();
    outFile->Write();
    outFile->Close();

    delete outFile;

    //// ForceField->Draw("l same"); // FORCE FIELD

    // plot solution
    if(!IsSolution) SetSolutions();

    /*
    hLIS_vs_KEnergy->Draw("hist same");
    hFlux_vs_KEnergy->Draw("hist same");
    hFluxFF_vs_KEnergy->Draw("hist same");
    */

    // hLIS_vs_Rigidity->Draw("hist same");
    // hFlux_vs_Rigidity->Draw("hist same");

    
    ccFluxSolution->Update();

    Printf("FINAL plotting");
    theApp->Run();

}



void SolarModulation::GetDifferentialIntensity(vector<int> r,string filename, int mult){
  
  int zz= (int)pZ;
  int aa= (int)pA;

  //Creating drawing tools
  TApplication* theApp = new TApplication(Form("App_Z%d_A%d",zz,aa), 0, 0);
  
  TCanvas* ccSolution = new TCanvas(Form("ccSolution_Z%d_A%d",zz,aa), Form("SOLUTION for Z=%d A=%d",zz,aa), 1400, 1000);

  TPad* d1 = new TPad(Form("Draw_Z%d_A%d",zz,aa),"Pad1",0,0,0.5,1,0,1);
  TPad* d2 = new TPad(Form("Draw_Z%d_A%d",zz,aa),"Pad2",0.5,0.5,1,1,0,1);
  TPad* d3  = new TPad(Form("Draw_Z%d_A%d",zz,aa),"Pad3",0.5,0,1,0.5,0,1);
  
  d1->Draw();
  d2->Draw();
  d3->Draw();
  
  d1->cd();
  gPad->SetTopMargin(0.04);
  gPad->SetBottomMargin(0.1);
  float* xbins1 = new float[61];
  for(int i=0; i<61; ++i) xbins1[i] = 0.01*pow(10,i*1/10.);
  
  TH1F* H1 = new TH1F("" , "", 60, xbins1);
  gPad->SetLogx();
  gPad->SetLogy();
  gStyle->SetOptStat(0);
  
  H1->SetMaximum(1E4);
  H1->SetMinimum(.1);
  H1->Draw();
  H1->GetXaxis()->SetTitle("kinetic energy [MeV/n]");
  H1->GetYaxis()->SetTitle("J [ s^{-1} m^{2} sr (MeV/n)^{-1} ]");
  H1->GetXaxis()->SetTitleOffset(1.4);
  H1->GetYaxis()->SetTitleOffset(1.4);
  
    
  //creating graph containing the problem's solution
  vector<double> vf; //scale solution for better visibility by factores of 10^i (i=0,0.5,1,1,5..)
  
    for (int i = 0; i < r.size(); ++i){vf.push_back(pow(sqrt(10),i));}
    double* rigidita = new double[np];
    vector<double*> fF;
    double** fF1 = new double*[nr];
    vector<TGraph*> G;
    
    for(int i =0; i<nr; ++i){fF1[i] = new double[np];}
    for(int k = 0; k < r.size();++k){
      for (int j = 0; j < np; ++j){fF1[r[k]][j] = vf[k]*flux[r[k]][j];}
      fF.push_back(fF1[r[k]]);
    }
    
    for (int i = 0; i < r.size(); ++i){
      G.push_back(new TGraph(np,KEnergy,fF[i]));
      G[i]->SetLineColor(i+2);
      G[i]->SetLineWidth(3);
      G[i]->SetLineStyle(1);
      G[i]->Draw("C");
    }
    
    for(int i=0; i<np; ++i) rigidita[i] =  exp(lnRigidity[i]);
    
    // ---- drawing JLisT to compare data ----
    TF1* JLis_Scaled = new TF1(Form("JLisScaled_Z%d_A%d",zz,aa),this,&SolarModulation::fJLis_ScaledValue, 1,1000000,0,"SolarModulation","JLis");
    JLis_Scaled->SetLineColor(1);
    JLis_Scaled->SetLineWidth(3);
    JLis_Scaled->Draw("same");
    
    //draw the force field solution
    ForceField->SetLineColor(2);
    ForceField->SetLineWidth(3);
    ForceField->SetLineStyle(2);
    ForceField->Draw("same");
    
    //draw data points
    for (int i=0; i<(int)DataSaver.size(); ++i) {
        DataSaver[i]->SetMarkerStyle(23+3*i);
        DataSaver[i]->SetMarkerColor(30+4*i);
        DataSaver[i]->SetMarkerSize(1);
        DataSaver[i]->Draw("P");
    }
    
    //Legend of the TPad
    TLegend* legend = new TLegend(0.12,0.73,0.5,0.93);
    legend->SetHeader("");
    legend->AddEntry(JLis_Scaled,"JLisT(90 AU)","l");
    legend->AddEntry(ForceField,"Force Field with phi(r = 1AU) = 407 MV","l");
    legend->AddEntry(DataSaver[0],"Pamela 2006","p");
    legend->AddEntry(DataSaver[1],"AMS_01","p");
    legend->AddEntry(DataSaver[2],"IMP75","p");
    legend->AddEntry(DataSaver[3],"IMP87","p");
    for (int i = 0; i < r.size(); ++i){legend->AddEntry(G[i],Form("Solution: Flux, r = %d AU",r[i]/mult),"L");}
    legend->Draw("same");
    
    d2->cd();
    
    gPad->SetTopMargin(0.08);
    float* xbins2 = new float[31];
    for(int i=0; i<31; ++i) xbins2[i] = 10*pow(10,i*1/10.);
    
    TH1F* H2 = new TH1F("" , "", 30, xbins2);
    gPad->SetLogx();
    gPad->SetLogy();
    gStyle->SetOptStat(0);
    
    H2->SetMaximum(10);
    H2->SetMinimum(.01);
    H2->Draw();
    H2->GetXaxis()->SetTitle("kinetic energy [MeV/n]");
    H2->GetYaxis()->SetTitle("J [ s^{-1} m^{2} sr (MeV/n)^{-1} ]");
    H2->GetXaxis()->SetTitleOffset(1.4);
    H2->GetYaxis()->SetTitleOffset(1.4);
    
    G[0]->Draw("C");
    ForceField->Draw("same");

    JLisT->SetLineColor(1);
    JLisT->SetLineWidth(3);
    JLisT->Draw("same");

    
    for (int i=0; i<(int)DataSaver.size(); ++i) {
        DataSaver[i]->SetMarkerStyle(23+3*i);
        DataSaver[i]->SetMarkerColor(30+4*i);
        DataSaver[i]->SetMarkerSize(1);
        DataSaver[i]->Draw("P");
    }
    
    d3->cd();
    gPad->SetTopMargin(0.04);
    gPad->SetBottomMargin(0.2);
    float* xbins3 = new float[41];
    for(int i=0; i<41; ++i) xbins3[i] = 1*pow(10,i*1/10.);
    
    TH1F* H3 = new TH1F(Form("H3_Z%d_A%d",zz,aa), "", 40, xbins3);
    gPad->SetLogx();
    gPad->SetLogy();
    gStyle->SetOptStat(0);
    
    H3->SetMaximum(4E2);
    H3->SetMinimum(.01);
    H3->Draw();
    H3->GetXaxis()->SetTitle("kinetic energy [MeV/n]");
    H3->GetYaxis()->SetTitle("J [ s^{-1} m^{2} sr (MeV/n)^{-1} ]");
    H3->GetXaxis()->SetTitleOffset(1.4);
    H3->GetYaxis()->SetTitleOffset(1.4);
    
    vector<TGraph*> G1;
    
    for (int i = 0; i < r.size(); ++i){
        G1.push_back(new TGraph(np,KEnergy,flux[r[i]]));
        G1[i]->SetLineColor(i+2);
        G1[i]->SetLineWidth(3);
        G1[i]->SetLineStyle(1);
        G1[i]->Draw("C");
    }
    
    JLisT->SetLineColor(1);
    JLisT->SetLineWidth(3);
    JLisT->Draw("same");
    
    ccSolution->Modified();
    ccSolution->Update();

    // while(ccSolution->WaitPrimitive()) gSystem->ProcessEvents();
    // ccSolution->Print(filename.c_str());

    theApp->Run();

    
    //freeing memory
    delete[] xbins1;
    delete[] xbins2;
    delete[] xbins3;
    delete H1;
    delete H2;
    delete H3;
    for(int i=0; i<nr; ++i){
        delete[] fF1[i];
    }
    delete[] fF1;
    
    delete d1;
    delete d2;
    delete d3;
    delete ccSolution;
}


void SolarModulation::ImportFile(string FILE1, string FILE2,string FILE3,string FILE4){
    DataSaver.push_back(new TGraph(FILE1.c_str(),"%lg %lg",""));
    DataSaver.push_back(new TGraph(FILE2.c_str(),"%lg %lg",""));
    DataSaver.push_back(new TGraph(FILE3.c_str(),"%lg %lg",""));
    DataSaver.push_back(new TGraph(FILE4.c_str(),"%lg %lg",""));
    printf("[SolarModulation::ImportFile()] %d files imported!\n", (int)DataSaver.size());
}


void SolarModulation::ImportFile(vector<string> F) {
    for (int i=0; i<(int)F.size(); ++i) {
        DataSaver.push_back(new TGraph(F[i].c_str(),"%lg %lg",""));
    }
    printf("[SolarModulation::ImportFile()] %d files imported!\n", (int)DataSaver.size());
}


void SolarModulation::BestParameter(double KCoeffInitial, int i)
{

  int zz= (int)pZ;
  int aa= (int)pA;
  double eps = 1;

    Solve();

    TF1* f1 = new TF1(Form("1AU_Z%d_A%d",zz,aa), this,&SolarModulation::fSpline_1AU,0,10000000, 0, " ", " ");
    double Xi1 = DataSaver[i]->Chisquare(f1);
    double K1 = KCoeffInitial;
    double K2= KCoeffInitial + KCoeffInitial/10;
    KCoeff->SetParameter(0, K2);

    Solve();

    TF1* f2 = new TF1(Form("1AU_Z%d_A%d",zz,aa), this,&SolarModulation::fSpline_1AU,0,10000000, 0, " ", " ");
    double Xi2 = DataSaver[i]->Chisquare(f2);

    delete f1;
    delete f2;

    while(abs(Xi2-Xi1)>eps){
      K1 = K2;
      TRandom* aleatorio = new TRandom();
      aleatorio->SetSeed(time(NULL));
      if(aleatorio->Uniform(-1,1)>=0) K2+= K2/10;
      else K2-= K2/10;

      Solve();

        TF1* f1 = new TF1(Form("1AU_Z%d_A%d",zz,aa), this,&SolarModulation::fSpline_1AU,0,10000000, 0, " ", " ");
        Xi1 = DataSaver[i]->Chisquare(f1);
        KCoeff->SetParameter(0, K2);
        Solve();

        TF1* f2 = new TF1(Form("1AU_Z%d_A%d",zz,aa), this,&SolarModulation::fSpline_1AU,0,10000000, 0, " ", " ");
        Xi2 = DataSaver[i]->Chisquare(f2);
        delete f1;
        delete f2;
        BestPar = K2;
    }

}




void SolarModulation::PlotSolution_Rig(){
  int zz= (int)pZ;
  int aa= (int)pA;
  TApplication* theApp = new TApplication(Form("App_Z%d_A%d",zz,aa), 0, 0);

  //----- get Voyager-1 data -----

    TFile* inFileProton= new TFile(PATH+"/data/grVOYAGER1_pHe3_Spectra_Ekn_Rig.root","READ");
    inFileProton->cd();
    // Proton Flux vs GV
    TGraphAsymmErrors* grVoyager1_Proton_Flux= (TGraphAsymmErrors*)inFileProton->Get("grVoyager1_Proton_Flux_Rig");
    grVoyager1_Proton_Flux->SetMarkerStyle(20);
    grVoyager1_Proton_Flux->SetMarkerColor(kBlack);
    grVoyager1_Proton_Flux->SetMarkerSize(0.8);


    TFile* inFileRig4 = new TFile(PATH+"/data/grVOYAGER1_pHe4_Spectra_Ekn_Rig.root","READ");
    inFileRig4->cd();
    // He4 Flux vs GV
    TGraphAsymmErrors* grVoyager1_He4_Flux= (TGraphAsymmErrors*)inFileRig4->Get("grVoyager1_Helium4_Flux_Rig");
    grVoyager1_He4_Flux->SetMarkerStyle(20);
    grVoyager1_He4_Flux->SetMarkerColor(kBlack);
    grVoyager1_He4_Flux->SetMarkerSize(0.8);

    inFileRig4->Close();

    // --- GV->MV conversion ---
    for(int ee=0;ee<grVoyager1_Proton_Flux->GetN();ee++){
        grVoyager1_Proton_Flux->GetX()[ee] *= 1.e+3;
        grVoyager1_Proton_Flux->GetY()[ee] *= 1.e-3;
        grVoyager1_Proton_Flux->GetEY()[ee] *= 1.e-3;
        //grVoyager1_He3_Flux->GetEXlow()[ee] *= 1.e+3;
        //grVoyager1_He3_Flux->GetEXhigh()[ee] *= 1.e+3;
        //grVoyager1_He3_Flux->GetEYhigh()[ee] *= 1.e-3;
    }
        
    
    for(int ee=0;ee<grVoyager1_He4_Flux->GetN();ee++){
        grVoyager1_He4_Flux->GetX()[ee] *= 1.e+3;
        grVoyager1_He4_Flux->GetY()[ee] *= 1.e-3;
        grVoyager1_He4_Flux->GetEY()[ee] *= 1.e-3;
        //grVoyager1_He3_Flux->GetEXlow()[ee] *= 1.e+3;
        //grVoyager1_He3_Flux->GetEXhigh()[ee] *= 1.e+3;
        //grVoyager1_He3_Flux->GetEYhigh()[ee] *= 1.e-3;
    }
    
        
    
    //ricalcolo indr0;
    int indR0 = -1; // 1AU should be = mult = 5
    double dist=1.e+9;
    for(int rr=0;rr<nRadius;rr++){
      double ddist = fabs( radius[rr]/AU - 1. );
      if(ddist<dist){
        dist=ddist;
        indR0= rr;
      }
    }

    // --- near-Earth flux vs Rig | at R=1 AU ---
    for(int pp=0;pp<nRigidity;pp++){ // set flux vs rig
      double pRigidit = Rigidity[pp];
      double pFluxEkn = flux[indR0][pp];
      double pFluxRig = pFluxEkn*dEdR(pRigidit);
      flux[indR0][pp] = pFluxRig;
    }
    
 
    TGraph* grFluxSolution= new TGraph(np, Rigidity, flux[indR0]); // flux Rig
    grFluxSolution->SetLineColor(kAzure+1);
    grFluxSolution->SetLineWidth(3);


    double PMIN=100.; // 100 MV
    double PMAX=200000.; // 200 GV

    TH2F* hFrameRig= new TH2F("hFrameProtonEkn","hFrameProtonEkn",200, PMIN, PMAX, 200, 1.e-6, 1.e+3); // MeV
    //TH2F* hFrameProtonEkn= new TH2F("hFrameProtonEkn","hFrameProtonEkn",300, emin, emax, 300, 1.e-7, 1.e4); // MeV


    hFrameRig->GetXaxis()->SetTitle("Rigidity [MV]");
    hFrameRig->GetYaxis()->SetTitle("J(E) [MV^{ -1} m^{ -2} s^{ -1} sr^{ -1} ]");
    hFrameRig->SetTitle(0);
    hFrameRig->GetYaxis()->SetNdivisions(504);
    hFrameRig->SetLabelFont(42,"X");
    hFrameRig->SetLabelFont(42,"Y");
    hFrameRig->SetTitleFont(42,"X");
    hFrameRig->SetTitleFont(42,"Y");
    hFrameRig->GetXaxis()->SetLabelSize(0.05);
    hFrameRig->GetYaxis()->SetLabelSize(0.05);
    hFrameRig->GetXaxis()->SetTitleSize(0.05);
    hFrameRig->GetYaxis()->SetTitleSize(0.05);
    hFrameRig->GetXaxis()->SetTitleOffset(1.375);
    hFrameRig->GetYaxis()->SetTitleOffset(1.55);
    hFrameRig->SetStats(0);


    TCanvas* ccFluxSolution = new TCanvas(Form("ccFluxSolution_Z%d_A%d",zz,aa), Form("SOLUTION vs Rig for Z=%d A=%d Model LIS=%d",zz,aa,kLISModel), 600, 600); //1400, 1000);
    ccFluxSolution->cd();
    gPad->SetBottomMargin(0.18);
    gPad->SetLeftMargin(0.18);
    gPad->SetRightMargin(0.05);
    gPad->SetTopMargin(0.05);
    
    hFrameRig->Draw();
    
    gPad->SetLogx();
    gPad->SetLogy();
    
    grVoyager1_Proton_Flux->Draw("pZ");
    //grVoyager1_He4_Flux->Draw("pZ");
    
    JLisRig->SetLineWidth(3);
    grFluxSolution->SetLineWidth(3);

    JLisRig->Draw("l same");
    //ForceField->Draw("l same");
    grFluxSolution->Draw("l same");

    TCanvas* ccKCoeff = new TCanvas(Form("ccKCoeff"), Form("K vs Rigidity"), 600, 600); //1400, 1000);
    ccKCoeff->cd();
    TH2F* hFrameK= new TH2F("hFrameK","hFrameK",200, pmin, pmax, 200, KCoeff->Eval(pmin), KCoeff->Eval(pmax)); // MeV
    //TH2F* hFrameProtonEkn= new TH2F("hFrameProtonEkn","hFrameProtonEkn",300, emin, emax, 300, 1.e-7, 1.e4); // MeV

    hFrameK->GetXaxis()->SetTitle("Rigidity[MV]");
    hFrameK->GetYaxis()->SetTitle("K [m^{2} s^{-1}]");
    hFrameK->SetTitle(0);
    hFrameK->GetYaxis()->SetNdivisions(504);
    hFrameK->SetLabelFont(42,"X");
    hFrameK->SetLabelFont(42,"Y");
    hFrameK->SetTitleFont(42,"X");
    hFrameK->SetTitleFont(42,"Y");
    hFrameK->GetXaxis()->SetLabelSize(0.05);
    hFrameK->GetYaxis()->SetLabelSize(0.05);
    hFrameK->GetXaxis()->SetTitleSize(0.05);
    hFrameK->GetYaxis()->SetTitleSize(0.05);
    hFrameK->GetXaxis()->SetTitleOffset(1.375);
    hFrameK->GetYaxis()->SetTitleOffset(1.55);
    hFrameK->SetStats(0);
    gPad->SetBottomMargin(0.18);
    gPad->SetLeftMargin(0.18);
    gPad->SetRightMargin(0.05);
    gPad->SetTopMargin(0.05);
    hFrameK->Draw();
    gPad->SetLogx();
    gPad->SetLogy();
    KCoeff->SetLineWidth(3);
    KCoeff->Draw("l same");
    // store locally
    
    TFile* outFile;
    outFile= new TFile(Form(PATH+"/OUT/TGraph_FluxVSRig_data_LIS%d_Z%d_A%d.root",kLISModel,zz,aa),"recreate");
    outFile->cd();
    //ccFluxSolution->Write();
    //ForceField->Write();
    grFluxSolution->Write();
    outFile->Write();
    outFile->Close();

    delete outFile;

    //// ForceField->Draw("l same"); // FORCE FIELD

    // plot solution
    if(!IsSolution) SetSolutions();

    ccFluxSolution->Update();

    Printf("FINAL plotting");
    theApp->Run();

}



void SolarModulation::Iterate_ThisKScale(double iKScale, double a, double b,double rk){ // Grid gia definite
    //indexKScale= iKScale;
    ThisKScale2 = iKScale;
    index = a;
    indexb = b;
    indexrk = rk;

}


void SolarModulation::KscaleGrid(){
    
    //estratti da initmodulation
    // --- NT2017 create kscale grid | fixed ranges ---
    KScale  = new double[nk];
    lnKScale  = new double[nk];
    deltalogk = (log(kmax)-log(kmin))/(nk-1);
    
    
    for(int kk=0; kk<nk; ++kk) {
        lnKScale[kk] = log(kmin) + deltalogk*kk;
        //cout<<lnKScale[kk]<<endl;
    }
    
    
    
    for(int kk=0; kk<nk; ++kk) KScale[kk] = exp( lnKScale[kk] );
    // --- NT2017 create phi-scale grid ---
    Phi     = new double[nphi];
    for(int mm=0;mm<nphi;mm++) Phi[nphi-mm-1]=PHIxKSCALE/KScale[mm];

    // --- NT2018 create drift grid ---
    deltaxid= (xidmax-xidmin)/(nxid-1);
    XiD    = new double[nxid];
    for(int dd=0;dd<nxid;dd++) XiD[dd]= xidmin+ deltaxid*dd;

    // --- other grid vectors: radius, rigidity and energy ---
    deltar = (rmax-rmin)/(nr-1);
    deltalogp = (log(pmax)-log(pmin))/(np-1);
    deltaAr = (Armax - Armin)/(nr-1);  // adimensional grid

    radius = new double[nr];
    Rigidity = new double[np];
    KEnergy = new double [np];
    lnRigidity = new double[np];
    lnKEnergy = new double[np];

    for(int i=0; i<nr; ++i) radius[i] = rmin + deltar*i;
    for(int i=0; i<np; ++i) lnRigidity[i] = log(pmin)+deltalogp*i;
    for(int i=0; i<np; ++i) Rigidity[i] = exp(lnRigidity[i]);
    for(int i=0; i<np; ++i) KEnergy[i] = RigToEkn(Rigidity[i]); // (exp(lnRigidity[i]));
    for(int i=0; i<np; ++i) lnKEnergy[i] = log(KEnergy[i]);
// printf("Vectors created. Number of elements: %d and %d\n", nr, np);  //numeri della grid: r vs lnp
    
    
    
    
    
//estratti da makegrid istogram
    int zz= (int)pZ;
  int aa= (int)pA;

  //cout<<"Build grid for final histograms [ Z="<<zz<<", A="<<aa<<"] "<<endl;



  // ---- radius grid ----
  nRadius = nr;
  xRadius = new double[nRadius+1];
  for(int rr=0;rr<nRadius-1;rr++){
    double pRad1= radius[rr];
    double pRad2= radius[rr+1];
    double xRad = 0.5*(pRad1 + pRad2);
    xRadius[rr+1]= xRad; // AU UNITS!
  }
    

  // first and last values | LINEAR progression
  xRadius[0]= radius[0] - (xRadius[1]-radius[0]); // = 2*radius[0] - xRadius[1]
  xRadius[nRadius]= radius[nRadius-1] + (radius[nRadius-1]-xRadius[nRadius-1]);

    
  // NT2017: convert xRadius from Meter to AU
  for(int rr=0;rr<nRadius+1;rr++) xRadius[rr] /= AU;

  // ---- rigidity grid ----
  nRigidity = np;
  xRigidity = new double[nRigidity+1];
  xLnRigidity = new double[nRigidity+1];
  for(int pp=0;pp<nRigidity-1;pp++){
    double pR1 = exp( lnRigidity[pp] );
    double pR2 = exp( lnRigidity[pp+1] );
    double xR  = sqrt(pR1*pR2);
    xRigidity[pp+1]= xR;
    xLnRigidity[pp+1]= log(xR);
  }

  // first and last value | LOG PROGRESSION
  xRigidity[0] = Rigidity[0]*Rigidity[0]/xRigidity[1];
  xRigidity[nRigidity] = Rigidity[nRigidity-1]*Rigidity[nRigidity-1]/xRigidity[nRigidity-1];
  xLnRigidity[0]=log( xRigidity[0]);
  xLnRigidity[nRigidity]=log( xRigidity[nRigidity]);


  // ---- kinetic energy grid ----
  nKEnergy = nRigidity;
  xKEnergy = new double[nKEnergy+1];
  xLnKEnergy = new double[nKEnergy+1];
  for(int pp=0;pp<nRigidity+1;pp++){
    xKEnergy[pp] = RigToEkn( xRigidity[pp] );
    xLnKEnergy[pp]= log(xKEnergy[pp]);
  }

  // ---- kscale grid | LOG ----
  nKScale = nk;
  xKScale = new double[nKScale+1];
  xLnKScale = new double[nKScale+1];
  for(int pp=0;pp<nKScale-1;pp++){
    double pR1 = KScale[pp];
    double pR2 = KScale[pp+1];
    double xR  = sqrt(pR1*pR2);
    xKScale[pp+1]= xR;
    xLnKScale[pp+1]= log(xR);
  }

  // first and last value | LOG PROGRESSION
  xKScale[0] = (KScale[0]*KScale[0])/xKScale[1];
  xKScale[nKScale] = (KScale[nKScale-1]*KScale[nKScale-1])/xKScale[nKScale-1];
  xLnKScale[0]=log( xKScale[0]);
  xLnKScale[nKScale]=log( xKScale[nKScale]);

//    cout<<"Kscale grid min "<<xKScale[0]<<endl;
//    cout<<"Kscale grid max"<<xKScale[nKScale]<<endl;
    
  // ---- build PHI grid for TH1 ----
  nPhi = nphi;
  xPhi = new double[nPhi+1];
  for(int mm=0;mm<nPhi+1;mm++) xPhi[nPhi-mm]=PHIxKSCALE/xKScale[mm];

  // ---- xi-drift grid ----
  nXiD = nxid;
  xXiD = new double[nXiD+1];
  for(int dd=0;dd<nXiD-1;dd++){
    double pRad1= XiD[dd];
    double pRad2= XiD[dd+1];
    double xRad = 0.5*(pRad1 + pRad2);
    xXiD[dd+1]= xRad; // AU UNITS!
  }

  // first and last values | LINEAR progression
  xXiD[0]= XiD[0] - (xXiD[1]-XiD[0]);
  xXiD[nXiD]= XiD[nXiD-1] + (XiD[nXiD-1]-xXiD[nXiD-1]);
    
    IsGrid=true;
}



double SolarModulation::CurrentFlux(double E_c){

// --- best radius index for local flux at R0=1AU ---
  int indR0 = -1; // 1AU should be = mult = 5
  double dist=1.e+9;
    
//    cout<<"0 Inside Current FLux ..."<<endl;
    
  for(int rr=0;rr<nRadius;rr++){
    double ddist = fabs( radius[rr]/AU - 1. );
    if(ddist<dist){
      dist=ddist;
      indR0= rr;
    }
  }
    
    
// --- best energy index for flux in heliosphere --- [~1 GeV/n]
  int indE0 = -1; // 1AU should be = mult = 5
  dist=1.e+20;
  for(int ee=0;ee<nKEnergy;ee++){
    double ddist = fabs( KEnergy[ee] - E_c ); //MeV/n!
    if(ddist<dist){
        dist=ddist;
        indE0= ee;
     }
   }
   
    
   //trovare index di Ec:
    
    double flux_c = flux[indR0][indE0];
    
    return flux_c;
    
}






















void SolarModulation::PlotSolution_BR(int br){
  int zz= (int)pZ;
  int aa= (int)pA;
 

       
    TFile* inFileProton= new TFile("./data/SOHO.root","READ");
    
  //  TFile* inFileProton= new TFile("/Users/davidpelosi/Desktop/AMS-PAMELA_N0/PAMELA_LxPgModel/data/ssdc_canvas_Cut.root","READ");
   //TFile* inFileProton= new TFile("./data/ssdc_canvas.root","READ");

    inFileProton->cd();
    TGraphErrors* grAMS02_Proton_Flux_2015 = (TGraphErrors*)inFileProton->Get(Form("grFlux_SOHO_%d",br));
    
    //TGraphErrors* grPAMELA_Proton_Flux_2015 = (TGraphErrors*)inFileProton->Get(Form("graph_CUT%d",br));
    // TGraphErrors* grPAMELA_Proton_Flux_2015 = (TGraphErrors*)inFileProton->Get(Form("graph%d",br));
    

    grAMS02_Proton_Flux_2015->SetMarkerStyle(20);
    grAMS02_Proton_Flux_2015->SetMarkerColor(kBlack);
    grAMS02_Proton_Flux_2015->SetMarkerSize(0.8);
    inFileProton->Close();

    // --- GeV->MeV conversion ---
    for(int ee=0;ee<grAMS02_Proton_Flux_2015->GetN();ee++){
      grAMS02_Proton_Flux_2015->GetX()[ee] *= 1.e+3;
      grAMS02_Proton_Flux_2015->GetY()[ee] *= 1.e-3;
      
      grAMS02_Proton_Flux_2015->GetEX()[ee] *= 1.e+3;
      grAMS02_Proton_Flux_2015->GetEY()[ee] *= 1.e-3;


    }

    int indR0 = -1; // 1AU should be = mult = 5
    double dist=1.e+9;
    for(int rr=0;rr<nRadius;rr++){
      double ddist = fabs( radius[rr]/AU - 1. );
      if(ddist<dist){
        dist=ddist;
        indR0= rr;
      }
    }

    TGraph* grFluxSolution= new TGraph(np, KEnergy, flux[indR0]);

    //TGraph* grFluxSolution= new TGraph(np, KEnergy, flux[5]);
    //se voglio calcore sempre near-earth
    //TGraph* grFluxSolution= new TGraph(np, KEnergy, flux_at_r);
    grFluxSolution->SetLineColor(kAzure+1);
    grFluxSolution->SetLineWidth(3);
  
    double EMIN=100; // 100 MeV
    double EMAX=100000.; // 100 GeV

    TH2F* hFrameProtonEkn= new TH2F("hFrameProtonEkn","hFrameProtonEkn",200, EMIN, EMAX, 200, 1.e-6, 1.e+2); // MeV
    //TH2F* hFrameProtonEkn= new TH2F("hFrameProtonEkn","hFrameProtonEkn",300, emin, emax, 300, 1.e-7, 1.e4); // MeV


    hFrameProtonEkn->GetXaxis()->SetTitle("kinetic energy [MeV/n]");
    hFrameProtonEkn->GetYaxis()->SetTitle("J(E) [MeV^{ -1} m^{ -2} s^{ -1} sr^{ -1} ]");
    hFrameProtonEkn->SetTitle(0);
    hFrameProtonEkn->GetYaxis()->SetNdivisions(504);
    hFrameProtonEkn->SetLabelFont(42,"X");
    hFrameProtonEkn->SetLabelFont(42,"Y");
    hFrameProtonEkn->SetTitleFont(42,"X");
    hFrameProtonEkn->SetTitleFont(42,"Y");
    hFrameProtonEkn->GetXaxis()->SetLabelSize(0.05);
    hFrameProtonEkn->GetYaxis()->SetLabelSize(0.05);
    hFrameProtonEkn->GetXaxis()->SetTitleSize(0.05);
    hFrameProtonEkn->GetYaxis()->SetTitleSize(0.05);
    hFrameProtonEkn->GetXaxis()->SetTitleOffset(1.375);
    hFrameProtonEkn->GetYaxis()->SetTitleOffset(1.55);
    hFrameProtonEkn->SetStats(0);

    TCanvas* ccFluxSolution = new TCanvas(Form("ccFluxSolution_Z%d_A%d",zz,aa), Form("SOLUTION for Z=%d A=%d Model LIS=%d",zz,aa,kLISModel), 600, 600); //1400, 1000);
    ccFluxSolution->cd();
    gPad->SetBottomMargin(0.18);
    gPad->SetLeftMargin(0.18);
    gPad->SetRightMargin(0.05);
    gPad->SetTopMargin(0.05);
    
    hFrameProtonEkn->Draw(); //background th2f
    
    gPad->SetLogx();
    gPad->SetLogy();

   //grVOYAGER1_Hydrogen_Flux_2016->Draw("pZ");
//grVOYAGER1_Helium_Flux_2016->Draw("pZ");
    
    grAMS02_Proton_Flux_2015->Draw("pZ");
  //   grAMS02_Helium_Flux->Draw("pZ"); // METTERE...?

    JLisT->SetLineWidth(3);
    grFluxSolution->SetLineWidth(3);

    JLisT->Draw("l same");
    //ForceField->Draw("l same");
    grFluxSolution->Draw("l same");

    // store locally
    TFile* outFile;
    if (indexKmodel == 1)
    {
    outFile= new TFile(Form("./out/TGraph_Flux_data_LIS%d_Z%d_A%d_BR%d_DCStandard.root",kLISModel,zz,aa,br),"recreate");
    }
    
    if (indexKmodel == 2) { 
       outFile= new TFile(Form("./out/TGraph_Flux_data_LIS%d_Z%d_A%d_BR%d_DCPowerLaw.root",kLISModel,zz,aa,br),"recreate");
    }
 
   if (indexKmodel == 3) { 
       outFile= new TFile(Form("./out/TGraph_Flux_data_LIS%d_Z%d_A%d_BR%d_DCPotgeiter.root",kLISModel,zz,aa,br),"recreate");
    }
    if (indexKmodel == 4) { 
       outFile= new TFile(Form("./out/TGraph_Flux_data_LIS%d_Z%d_A%d_BR%d_DCPotgeiter_Rk_free.root",kLISModel,zz,aa,br),"recreate");
    }
   outFile->cd();
    

    grFluxSolution->Write();
    //outFile->Write();

    ccFluxSolution->SetName(Form("BR %d",br));
    ccFluxSolution->Write();
    JLisT->Write();
    grAMS02_Proton_Flux_2015->Write();

    outFile->Close();

    delete outFile;

    //// ForceField->Draw("l same"); // FORCE FIELD

    // plot solution
    //if(!IsSolution) SetSolutions();

    /*
    hLIS_vs_KEnergy->Draw("hist same");
    hFlux_vs_KEnergy->Draw("hist same");
    hFluxFF_vs_KEnergy->Draw("hist same");
    */
    // hLIS_vs_Rigidity->Draw("hist same");
    // hFlux_vs_Rigidity->Draw("hist same");    
    ccFluxSolution->Update();

    Printf("FINAL plotting");

   // theApp->Run();

}