//  open a root file and extract th2f, for each x bin create the projection on y and save it in a tgraph

//import root libraries
#include "../src/Conversion.h"
#include "TFile.h"
#include "TH2F.h"
#include "TGraphErrors.h"
//include tcollection
//include tmultigraph
#include "TMultiGraph.h"
#include "TCollection.h"
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
using namespace std;


/*int main() {

//extract 2 tgrapherrors from root file and merge them in a single tgrapherrors
//open root file
TFile *f = new TFile("pamela_Flux_time.root");
f->cd();

//merge two tgrapherrors in a single tgrapherrors
TGraphErrors *graph1 = (TGraphErrors*)f->Get("graph1");
TGraphErrors *graph2 = (TGraphErrors*)f->Get("graph2");
TGraphErrors *graph3 = new TGraphErrors();

//create a tcollection
//TCollection *list = new TCollection();
TCollection *list;
list->Add(graph1);
list->Add(graph2);

graph3->Merge(list);

f->Close();

//create a new root file
TFile *f1 = new TFile("pamela_Flux_time_merged.root","RECREATE");
graph1->Write();
f1->Close();
}
*/

int main() {
TFile *f = new TFile("pamela_Flux_time.root");
f->cd();

//merge two tgrapherrors in a single tgrapherrors
TGraphErrors *graph1 = (TGraphErrors*)f->Get("graph1");
TGraphErrors *graph2 = (TGraphErrors*)f->Get("graph2");


double* xEknPAMELA_1 = (double*)graph1->GetX();
double* xEknPAMELA_2 = (double*)graph2->GetX();

double* xEknPAMELA = new double[graph1->GetN()+graph2->GetN()];

for (int i=0; i<graph1->GetN()+graph2->GetN(); i++) {
  if (i < graph1->GetN())
  {
      xEknPAMELA[i] = xEknPAMELA_1[i];
  }
else
  {
      xEknPAMELA[i] = xEknPAMELA_2[i-graph1->GetN()];
  }
}

for (int i=0; i<graph1->GetN()+graph2->GetN(); i++) {
  cout << i <<"\t"<<xEknPAMELA[i] << endl;
}

f->Close();
}






