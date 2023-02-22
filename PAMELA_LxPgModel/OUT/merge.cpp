
void merge() {

TFile *f = new TFile("K0_vs_timep29.root");
f->cd();

//merge two tgrapherrors in a single tgrapherrors
TGraphErrors *graph1 = (TGraphErrors*)f->Get("grK0VSTime_PAMELA");

TFile *f2 = new TFile("K0_vs_time.root");
f2->cd();
TGraphErrors *graph2 = (TGraphErrors*)f2->Get("grK0VSTime_PAMELA");


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
TFile *f1 = new TFile("grK0VSTime_PAMELA.root","RECREATE");

mergedGraph->SetName(f1->GetName());

//draw a TgraphErrors only points
mergedGraph->SetMarkerStyle(20);
mergedGraph->SetMarkerColor(kRed);
mergedGraph->SetMarkerSize(0.5);
//time conversion of the x axis
mergedGraph->GetXaxis()->SetTimeDisplay(1); 
mergedGraph->GetXaxis()->SetTimeFormat("%d/%m/%y %H:%M");
mergedGraph->GetXaxis()->SetTimeOffset(0,"gmt");
mergedGraph->GetXaxis()->SetNdivisions(-503);
mergedGraph->GetXaxis()->SetLabelSize(0.03);
mergedGraph->GetXaxis()->SetTitleSize(0.03);
mergedGraph->GetXaxis()->SetTitleOffset(1.2);
mergedGraph->GetXaxis()->SetTitle("Date");
mergedGraph->GetYaxis()->SetTitle("K0 flux");
mergedGraph->GetYaxis()->SetTitleOffset(1.2);
mergedGraph->GetYaxis()->SetTitleSize(0.03);
mergedGraph->GetYaxis()->SetLabelSize(0.03);
mergedGraph->Write();

f1->Close();

}




//create a tmultigraph

/*
void merge() {

TFile *f = new TFile("K0_vs_timep29.root");
f->cd();
TGraphErrors *graph1 = (TGraphErrors*)f->Get("grK0VSTime_PAMELA");


TFile *f2 = new TFile("K0_vs_time.root");
f2->cd();
TGraphErrors *graph2 = (TGraphErrors*)f2->Get("grK0VSTime_PAMELA");

f->Close();
f2->Close();

TMultiGraph *mg = new TMultiGraph();
//add the two tgrapherrors to the multigraph
mg->Add(graph1);
mg->Add(graph2);
//draw the multigraph
mg->Draw("AP");

//save the multigraph in a root file
TFile *f1 = new TFile("grK0VSTime_PAMELA.root","RECREATE");
//mg->SetName();
mg->Write();

f1->Close();

}
*/
