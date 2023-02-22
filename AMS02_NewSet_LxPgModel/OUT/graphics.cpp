//extract tgrapherorrs from a root file

void graphics()
{
  TFile *f = new TFile("K0_vs_time_Potgeiter.root");
  TGraphErrors *gr = (TGraphErrors*)f->Get("grK0vsTime_AMS02");

//setpoint errors
    int np = gr->GetN();
   for (int i=0; i<np; i++) {
      gr->SetPointError(i,0,gr->GetErrorY(i));
   }



    gr->SetMarkerStyle(20);
    gr->SetMarkerColor(4);
    gr->SetLineColor(4);
    gr->SetLineWidth(2);
    gr->SetMarkerSize(0.2);
    gr->Draw("APL");
    
   


    //save the graph in a root file
    TFile *f1 = new TFile("graphics.root","RECREATE");
    gr->Write();
    f1->Close();

}