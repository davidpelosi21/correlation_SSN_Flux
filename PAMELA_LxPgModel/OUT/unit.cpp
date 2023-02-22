void unit() {

  // Open the second root file and extract the TGraphErrors object
  TFile* file2 = TFile::Open("K0_vs_time_Potgeiter.root");
  TGraphErrors* graph2 = (TGraphErrors*)file2->Get("grK0VSTime_PAMELA");

  // Save the merged TGraphErrors to a new root file
TFile* file3 = TFile::Open("K0_vs_time_PotgeiterALL.root", "RECREATE");
graph2->SetName("grK0VSTime_PAMELA");
graph2->SetMarkerStyle(20);
graph2->SetMarkerSize(0.5);
graph2->SetMarkerColor(kRed);
graph2->SetLineColor(kRed);
graph2->SetLineWidth(2);
graph2->GetXaxis()->SetTitle("Date");
//mg->GetYaxis()->SetTitle("");


//x axis in time format
graph2->GetXaxis()->SetTimeDisplay(1);
graph2->GetXaxis()->SetTimeFormat("%d/%m/%y");

graph2->GetXaxis()->SetTimeOffset(0,"gmt");


  graph2->Write();

  file3->Close();

  delete graph2;
  delete file3;
  delete file2;
}