//extract two tgrapherrors from a root file and merge them in a single tgrapherrors

void N0(){

TString path = "./";
//graph2 k vs t for PAMELA
TFile *f1 = new TFile("/Users/davidpelosi/Desktop/Correlation-testing/PAMELA_LxPgModel/OUT/N0_vs_time_Potgeiter20.root");
f1->cd();
TGraphErrors *graph1 = (TGraphErrors*)f1->Get("grN0vsTime_PAMELA");

graph1->SetMarkerStyle(1);
graph1->SetMarkerColor(kRed);
graph1->SetLineColor(kRed);
graph1->SetLineWidth(2);

TFile *f2 = new TFile("/Users/davidpelosi/Desktop/Correlation-testing/PAMELA_LxPgModel/OUT/N0_vs_time_Potgeiter20-82.root");
f2->cd();
TGraphErrors *graph2 = (TGraphErrors*)f2->Get("grN0vsTime_PAMELA");

graph2->SetMarkerStyle(1);
graph2->SetMarkerColor(kRed);
graph2->SetLineColor(kRed);
graph2->SetLineWidth(2);

//merge graph1 and graph2 in a single graph
TGraphErrors *graph = new TGraphErrors();
int n1 = graph1->GetN();
int n2 = graph2->GetN();
int n = n1 + n2;
graph->Set(n);
for(int i=0; i<n1; i++){
double x,y;
graph1->GetPoint(i,x,y);
graph->SetPoint(i,x,y);

double ex,ey;
ex = graph1->GetErrorX(i);
ey = graph1->GetErrorY(i);
graph->SetPointError(i,ex,ey);
}

for(int i=0; i<n2; i++){
double x,y;
graph2->GetPoint(i,x,y);
graph->SetPoint(i+n1,x,y);

double ex,ey;
ex = graph2->GetErrorX(i);
ey = graph2->GetErrorY(i);
graph->SetPointError(i+n1,ex,ey);
}

graph->SetMarkerStyle(1);
graph->SetMarkerColor(kRed);
graph->SetLineColor(kRed);
graph->SetLineWidth(2);
graph->SetMarkerSize(0.2);


//x axis in time format
graph->GetXaxis()->SetTimeDisplay(1);
graph->GetXaxis()->SetTimeFormat("%d/%m/%y");
graph->GetXaxis()->SetTimeOffset(0,"gmt");


graph->Draw("AP");

f1->Close();
f2->Close();

//save in a root file
TFile *f = new TFile("N0_vs_time_PotgeiterALL.root","RECREATE");
f->cd();
graph->SetName("grN0vsTime_PAMELA");
graph->Write();
f->Close();



}