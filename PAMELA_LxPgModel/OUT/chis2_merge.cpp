//create a macro to extract two tgrapherrors from a root file and add to a TmultiGraph
//the macro is called by the script chis2_merge.cpp
//the macro is called by the script chis2_merge.cpp


void chi2_merge(){
    TMultiGraph *mg = new TMultiGraph();
    TGraphErrors *chis2_standard;
    TGraphErrors *chis2_powerlaw;

    TFile* chis2_standard_File= new TFile("Chis2_vs_time_StandardModel.root","READ");
    TFile* chis2_powerlaw_File= new TFile("Chis2_vs_time_PowerLaw.root","READ");
    

    chis2_standard= (TGraphErrors*)chis2_standard_File->Get("grChis2vsTime_PAMELA");
 
 chis2_standard->SetMarkerColor(2);
    chis2_standard->SetMarkerStyle(20);
    chis2_standard->SetMarkerSize(0.5);
    chis2_standard->SetLineColor(2);
    chis2_standard->SetLineWidth(2);

    chis2_powerlaw= (TGraphErrors*)chis2_powerlaw_File->Get("grChis2vsTime_PAMELA");


mg->Add(chis2_standard);
mg->Add(chis2_powerlaw);
//convert x axis in date fomat
mg->GetXaxis()->SetTimeDisplay(1);
mg->GetXaxis()->SetTimeFormat("%d/%m/%y");
mg->GetXaxis()->SetTimeOffset(0,"gmt");


//add axis
//create a canvas
TCanvas *c1 = new TCanvas("c1","c1",800,600);
c1->cd();
mg->GetXaxis()->SetTitle("Date");
mg->GetYaxis()->SetTitle("#chi^{2}/d.o.f");
mg->Draw("AP");

//add legend    
TLegend *leg = new TLegend(0.1,0.7,0.48,0.9);
leg->AddEntry(chis2_standard,"Standard","p");
leg->AddEntry(chis2_powerlaw,"PowerLaw","p");
leg->Draw();

//add mg to a root file
TFile* chis2_merge_File= new TFile("Chis2_Standard_vs_PowerLaw.root","RECREATE");
chis2_merge_File->cd();
c1->Write();



chis2_merge_File->Close();

}
