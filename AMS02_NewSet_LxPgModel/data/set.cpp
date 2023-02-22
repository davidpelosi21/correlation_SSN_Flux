void res() {
//extract an array of tgrapherror from a root file
//create an array of tgrapherrors
TFile *f = new TFile("AMS02.root");

TGraphErrors *g1[115];
TGraphErrors *g[115];

for (int j=0; j<115; j++) {
   g[j] = new TGraphErrors();

    if( j != 46 && j != 47 ){

   g[j] = (TGraphErrors*)f->Get(Form("grProton_vs_Ekn_T%d",j));

   g1[j] = new TGraphErrors();
   g1[j]->SetName(Form("grProton_vs_Ekn_T%d",j));

    cout<<j<<" ok"<<endl;


   for (int i=0; i<23; i++)  { 
   double x,y;
   g[j]->GetPoint(i,x,y);
   g1[j]->SetPoint(i,x,y);
   g1[j]->SetPointError(i,0,g[j]->GetErrorY(i));
   g1[j]->SetMarkerStyle(20);
   g1[j]->SetMarkerSize(0.5);
   }
}

}

//save the array of tgrapherrors in a root file

TFile *f1 = new TFile("AMS02_RES23.root","RECREATE");
f1->cd();
for (int j=0; j<115; j++) {
    if (  j != 46 && j != 47 )
    {
       g1[j]->Write();
    }

}
f1->Close();



}