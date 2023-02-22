void MergeGraphs() {
  // Open the first root file and extract the TGraphErrors object
  TFile* file1 = TFile::Open("/PAMELA_LxPgModel/OUT/Cross_K0_ProxySSN.root");
  TGraphErrors* graph1 = (TGraphErrors*)file1->Get("Cross_K0_ProxySSN");

  // Open the second root file and extract the TGraphErrors object
  TFile* file2 = TFile::Open("/AMS02_NewSet_LxPgModel/OUT/Cross_K0_ProxySSN.root");
  TGraphErrors* graph2 = (TGraphErrors*)file2->Get("Cross_K0_ProxySSN");

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
  TFile* file3 = TFile::Open("cross_AMS02_PAMELA.root", "RECREATE");
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