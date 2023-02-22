void unit2() {
  // Open the first root file and extract the TGraphErrors object
  TFile* file1 = TFile::Open("B_vs_time_Potgeiter23.root");
  TGraphErrors* graph1 = (TGraphErrors*)file1->Get("grBvsTime_PAMELA");

  // Open the second root file and extract the TGraphErrors object
  TFile* file2 = TFile::Open("B_vs_time_Potgeiter24.root");
  TGraphErrors* graph2 = (TGraphErrors*)file2->Get("grBvsTime_PAMELA");

  // Merge the two TGraphErrors into a new TGraphErrors object
  int nPoints = graph1->GetN() + graph2->GetN();
  double* x = new double[nPoints];
  double* y = new double[nPoints];
  double* ex = new double[nPoints];
  double* ey = new double[nPoints];
  for (int i = 0; i < graph1->GetN(); i++) {
    x[i] = graph1->GetX()[i];
    y[i] = graph1->GetY()[i];
    ex[i] = graph1->GetEX()[i];
    ey[i] = graph1->GetEY()[i];
  }
  for (int i = 0; i < graph2->GetN(); i++) {
    x[i + graph1->GetN()] = graph2->GetX()[i];
    y[i + graph1->GetN()] = graph2->GetY()[i];
    ex[i + graph1->GetN()] = graph2->GetEX()[i];
    ey[i + graph1->GetN()] = graph2->GetEY()[i];
  }
  TGraphErrors* graph3 = new TGraphErrors(nPoints, x, y, ex, ey);

  // Save the merged TGraphErrors to a new root file
  TFile* file3 = TFile::Open("B_vs_time_Potgeiter.root", "RECREATE");
  graph3->SetMarkerStyle(20);
graph3->SetMarkerSize(0.5);
graph3->SetMarkerColor(kRed);
graph3->SetLineColor(kRed);
graph3->SetLineWidth(2);
graph3->GetXaxis()->SetTitle("Date");
//mg->GetYaxis()->SetTitle("");


//x axis in time format
graph3->GetXaxis()->SetTimeDisplay(1);
graph3->GetXaxis()->SetTimeFormat("%d/%m/%y");

graph3->GetXaxis()->SetTimeOffset(0,"gmt");

  graph3->Write();
  file3->Close();

  // Clean up memory
  delete[] x;
  delete[] y;
  delete[] ex;
  delete[] ey;
  delete graph1;
  delete graph2;
  delete graph3;
  delete file1;
  delete file2;
}