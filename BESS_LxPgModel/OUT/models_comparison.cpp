//create a macro to extract 3 graphs from a root file and add to a TmultiGraph


void compare(int br){

    TMultiGraph *mg = new TMultiGraph();
    TGraphErrors *model_standard;
    TGraphErrors *model_powerlaw;
    TGraphErrors *model_pot;

  TGraphErrors *ams;
   TF1 *lis = new TF1();



    TFile* model_standard_File= new TFile(Form("TGraph_Flux_data_LIS4_Z1_A1_BR%d_DCStandard.root",br),"READ");
 
  TFile* model_power_File= new TFile(Form("TGraph_Flux_data_LIS4_Z1_A1_BR%d_DCPowerLaw.root",br),"READ");

   TFile* model_pot_File= new TFile(Form("TGraph_Flux_data_LIS4_Z1_A1_BR%d_DCPotgeiter.root",br),"READ");
 

 //extract tf1 from a canvas saved a root file


lis = (TF1*)model_standard_File->Get("fLIS_vs_EKN_Z1_A1");
ams = (TGraphErrors*)model_standard_File->Get(Form("grProton_vs_Ekn_T%d",br));


    model_standard= (TGraphErrors*)model_standard_File->Get("Graph");
   model_powerlaw= (TGraphErrors*)model_power_File->Get("Graph");
   model_pot= (TGraphErrors*)model_pot_File->Get("Graph");

    model_standard->SetMarkerColor(2);
    model_standard->SetMarkerStyle(20);
    model_standard->SetMarkerSize(0.1);
    model_standard->SetLineColor(2);
    model_standard->SetLineWidth(1);

    model_powerlaw->SetMarkerColor(4);
    model_powerlaw->SetMarkerStyle(20);
    model_powerlaw->SetMarkerSize(0.1);
    model_powerlaw->SetLineColor(4);
    //set line width
    model_powerlaw->SetLineWidth(1);

    model_pot->SetMarkerColor(6);
    model_pot->SetMarkerStyle(20);
    model_pot->SetMarkerSize(0.1);
    model_pot->SetLineColor(6);
    model_pot->SetLineWidth(1);

//lis set green line color
lis->SetMarkerColor(kGreen+2);
lis->SetLineWidth(2);

ams->SetMarkerColor(1);
//ams->SetMarkerStyle(8);
ams->SetMarkerSize(0.1);
ams->SetLineColor(1);


    mg->Add(model_standard,"APL");
    mg->Add(model_powerlaw,"APL");
    mg->Add(model_pot,"APL");
    mg->Add(ams,"AP");

//log scale
    //gPad->SetLogy();
    //gPad->SetLogx();
    
//canvas
    TCanvas *c1 = new TCanvas("c1","c1",800,600);
    c1->cd();
    c1->SetLogy();
    c1->SetLogx();
    mg->GetXaxis()->SetTitle("Energy [MeV]");
    mg->GetYaxis()->SetTitle("Flux [cm^{-2} s^{-1} sr^{-1} MeV^{-1}]");
    mg->Draw("AP");
    lis->Draw("same");
 

//add tgraph to a multigraph with AP option


    //add mg legend
    TLegend *leg = new TLegend(0.1,0.7,0.48,0.9);
    leg->AddEntry(model_standard,"Standard","l");
    leg->AddEntry(model_powerlaw,"PowerLaw","l");
    leg->AddEntry(model_pot,"Potgeiter","l");
    leg->AddEntry(lis,"LIS","l");
    leg->AddEntry(ams,"AMS","l");
    leg->Draw();
    c1->Update();

//add canvas to a root file
    TFile* model_comparison_File= new TFile(Form("Flux_data_LIS4_Z1_A1_BR%d_DCComparison.root",br),"RECREATE");
    model_comparison_File->cd();
    c1->Write();
    //set log sclae axis
    model_comparison_File->Close();

}
