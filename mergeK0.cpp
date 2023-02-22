//create a tmultigraph
void merge() {
TString path = "./";
//graph2 k vs t for PAMELA
TFile *f1 = new TFile(path+"/PAMELA_LxPgModel/OUT/K0_vs_time_StandardModel.root");
f1->cd();
TGraphErrors *graph1 = (TGraphErrors*)f1->Get("grK0VSTime_PAMELA");

graph1->SetMarkerStyle(1);
graph1->SetMarkerColor(4);
graph1->SetLineColor(4);
graph1->SetLineWidth(2);
graph1->SetMarkerSize(0.2);

//cout<<"OK"<endl;

TFile *f2 = new TFile(path+"/AMS02_NewSet_LxPgModel/OUT/K0_vs_time_StandardModel.root");
f2->cd();
TGraphErrors *graph2 = (TGraphErrors*)f2->Get("grK0VSTime_AMS02");


graph2->SetMarkerStyle(1);
graph2->SetMarkerColor(kRed);
graph2->SetLineColor(kRed);
graph2->SetLineWidth(2);
graph2->SetMarkerSize(0.2);


TFile *f3 = new TFile(path+"BESS_LxPgModel/OUT/K0_vs_time_StandardModel.root");
f3->cd();
TGraphErrors *graph3 = (TGraphErrors*)f3->Get("grK0VSTime_BESS");
graph3->SetMarkerStyle(1);
graph3->SetMarkerColor(kGreen+1);
graph3->SetLineColor(kGreen+1);
graph3->SetLineWidth(2);
graph3->SetMarkerSize(0.2);


TFile *f4 = new TFile(path+"SOHO_LxPgModel/OUT/K0_vs_time_StandardModel.root");
f4->cd();
TGraphErrors *graph4 = (TGraphErrors*)f4->Get("grK0VSTime_SOHO");
graph4->SetMarkerStyle(1);
graph4->SetMarkerColor(kOrange+2);
graph4->SetLineColor(kOrange+2);
graph4->SetLineWidth(2);
graph4->SetMarkerSize(0.2);



f1->Close();
f2->Close();
f3->Close();


//add graph2 to tmultigraph
TMultiGraph *mg = new TMultiGraph();
mg->Add(graph1);
mg->Add(graph2);
mg->Add(graph3);
mg->Add(graph4);
//set multigraph options
//mg->SetTitle("A vs time");
//mg->GetXaxis()->SetTitle("Date");
//mg->GetYaxis()->SetTitle("K_{0} [10^{22} cm^{2} s^{-1}]");
mg->GetYaxis()->CenterTitle();
mg->GetYaxis()->SetTitle("K_{0} [10^{22} cm^{2} s^{-1}]");
//mg->GetYaxis()->SetTitle("b(t)");
mg->GetXaxis()->CenterTitle();
//mg->GetXaxis()->SetTitle("S(t)");

//se lag
mg->GetXaxis()->SetTitle("");



//create canvas
TCanvas *c1 = new TCanvas("c1","c1",800,600);
c1->cd();
//mg->GetYaxis()->SetRangeUser(0, 2);

//x axis time format
mg->GetXaxis()->SetTimeDisplay(1);
mg->GetXaxis()->SetTimeFormat("%d/%m/%y");
mg->GetXaxis()->SetTimeOffset(0,"gmt");


mg->Draw("AP");
//add legend
TLegend *leg = new TLegend(0.1,0.7,0.48,0.9);
leg->AddEntry(graph1,"PAMELA","l");
leg->AddEntry(graph2,"AMS02 ","l");
leg->AddEntry(graph3,"BESS ","l");
leg->AddEntry(graph4,"SOHO ","l");
leg->Draw();

//save the multigraph in a root file
TFile *f5 = new TFile("K0_Standard.root","RECREATE");
c1->Write();
graph1->Write();
graph2->Write();
graph3->Write();
graph4->Write();



f5->Close();


}


