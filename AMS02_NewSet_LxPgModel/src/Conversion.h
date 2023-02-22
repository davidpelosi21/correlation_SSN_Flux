#include <stdio.h>
#include <cmath>
#include <cstdio>
#include <vector>

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
#include "TRandom.h"
#include <ctime>

using namespace std;

int pZ = 1; // proton charge
double ProtonMass = 0.9382720813; // proton mass
int pA = 1; // proton mass number
double pM = pA*ProtonMass; // proton mass


    
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

 