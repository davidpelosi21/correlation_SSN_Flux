//extract tgrapherrors from root file and merge them in a single tgrapherrors
//open root file

void merge() {

TFile *f = new TFile("pamela_Flux_time.root");
f->cd();

//merge two tgrapherrors in a single tgrapherrors
TGraphErrors *graph1 = (TGraphErrors*)f->Get("graph1");
TGraphErrors *graph2 = (TGraphErrors*)f->Get("graph2");

TGraphErrors *mergedGraph = new TGraphErrors();

int n1 = graph1->GetN();
int n2 = graph2->GetN();

for (int i = 0; i < n1; i++) {
  double x, y;
  graph1->GetPoint(i, x, y);
  double ex = graph1->GetErrorX(i);
  double ey = graph1->GetErrorY(i);
  mergedGraph->SetPoint(i, x, y);
  mergedGraph->SetPointError(i, ex, ey);
}

for (int i = 0; i < n2; i++) {
  double x, y;
  graph2->GetPoint(i, x, y);
  double ex = graph2->GetErrorX(i);
  double ey = graph2->GetErrorY(i);
  mergedGraph->SetPoint(i + n1, x, y);
  mergedGraph->SetPointError(i + n1, ex, ey);
}

//mergedGraph->Draw("AP");
//save tgrapherrors in a root file
TFile *f1 = new TFile("pamela_Flux_time_merged.root","RECREATE");

mergedGraph->SetName("MergedGraph");

mergedGraph->Write();
f1->Close();

}