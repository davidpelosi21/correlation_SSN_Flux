void residuals2() {
    TFile *f = new TFile("/Users/davidpelosi/Desktop/AMS-PAMELA_SOHO_BESS/PAMELA_LxPgModel/OUT/TGraph_Flux_data_LIS4_Z1_A1_BR30_DCStandard.root");
    TGraphErrors *gr_data2 = (TGraphErrors*)f->Get("graph30");

    double x,y;
    gr_data2->GetPoint(3,x,y);
    cout << "x = " << x << " y = " << y << endl;


    TGraphErrors *gr_res = new TGraphErrors();
    gr_res->SetMarkerStyle(20);
    gr_res->SetMarkerColor(kBlack);
    gr_res->SetMarkerSize(0.5);
    gr_res->SetLineColor(kBlack);
    gr_res->SetLineWidth(3);

    //fill the tgrapherrors
    for (int i=0; i<gr_data2->GetN(); i++) {
        double x,y;
        gr_data2->GetPoint(i,x,y);
        gr_res->SetPoint(i,x,y);
        gr_res->SetPointError(i,0,gr_data2->GetErrorY(i));
    }

    
    gr_res->Draw("AP");
//gr_data2->Draw("AP");
//save tgraph in a root file
    TFile *f2 = new TFile("residuals.root","RECREATE");
    gr_res->SetName("residuals");
    gr_res->Write();
    f2->Close();
}