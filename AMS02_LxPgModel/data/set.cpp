void res() {
//extract an array of tgrapherror from a root file

//create an array of tgrapherrors
TFile *f = new TFile("grAMS02_ProtonHelium_vs_Ekn_vs_Time_DEC2017.root");
TGraphErrors *g1[81];
TGraphErrors *g[81];

for (int j=0; j<81; j++) {
   g[j] = new TGraphErrors();
   g[j] = (TGraphErrors*)f->Get(Form("grProton_vs_Ekn_T%d",j));

   g1[j] = new TGraphErrors();
   g1[j]->SetName(Form("grProton_vs_Ekn_T%d",j));

for (int i=0; i<30; i++) { 
   double x,y;
   g[j]->GetPoint(i,x,y);
   g1[j]->SetPoint(i,x,y);
   g1[j]->SetPointError(i,0,0.3*g[j]->GetErrorY(i));
   g1[j]->SetMarkerStyle(20);
    g1[j]->SetMarkerSize(0.5);

}

}

//save the array of tgrapherrors in a root file
TFile *f1 = new TFile("grAMS02_ProtonHelium_vs_Ekn_vs_Time_DEC2017_RES.root","RECREATE");
f1->cd();
for (int j=0; j<81; j++) {
   g1[j]->Write();
}
f1->Close();



}