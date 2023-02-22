//import two th2f from a root file, take the projection on y for each th2f and combine the two projections in a tgrapherrors
//

void ams() {

TFile *f = new TFile("ProtonMonthlyRebinProj.root");

TGraphErrors *flux_err_stat[12];
TGraphErrors *flux_err_stat2[12];


for (int i=0; i<12; i++) {
   flux_err_stat[i] = (TGraphErrors*)f->Get(Form("flux_err_stat_%d",i));
   flux_err_stat2[i] = (TGraphErrors*)f->Get(Form("flux_err_stat_%d",i));
}

TGraphErrors *flux_err_stat_ams[12];

for (int i=0; i<12; i++) {
   flux_err_stat_ams[i] = new TGraphErrors();
   flux_err_stat_ams[i]->SetName(Form("flux_err_stat_ams_%d",i));
   for (int j=0; j<flux_err_stat[i]->GetN(); j++) {
      double x,y;
      flux_err_stat[i]->GetPoint(j,x,y);
      double x2,y2;
      flux_err_stat2[i]->GetPoint(j,x2,y2);
      flux_err_stat_ams[i]->SetPoint(j,x,y+y2);
      flux_err_stat_ams[i]->SetPointError(j,0,flux_err_stat[i]->GetErrorY(j)+flux_err_stat2[i]->GetErrorY(j));
   }
}

//modify point of a tgraph

                     
//save g in a root file
TFile *f1 = new TFile("ProtonMonthlyRebinProj_ams.root","RECREATE");
f1->cd();
for (int i=0; i<12; i++) {
   flux_err_stat_ams[i]->Write();
}

f1->Close();
f->Close();

cout << "End" << endl;

}
