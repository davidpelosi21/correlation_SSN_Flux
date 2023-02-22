void residuals3() {

    TFile *f = new TFile("/Users/davidpelosi/Desktop/AMS-PAMELA_SOHO_BESS/PAMELA_LxPgModel/OUT/residuals.root");
    TGraphErrors *gr = (TGraphErrors*)f->Get("residuals");


    //TFile *f = new TFile("/Users/davidpelosi/Desktop/AMS-PAMELA_SOHO_BESS/PAMELA_LxPgModel/OUT/TGraph_Flux_data_LIS4_Z1_A1_BR30_DCStandard.root");
    //TGraphErrors *gr = (TGraphErrors*)f->Get("graph30");
       
    gr->Draw("AP");

}

