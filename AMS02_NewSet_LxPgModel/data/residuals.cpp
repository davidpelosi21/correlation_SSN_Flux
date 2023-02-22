
void residuals() {

//create canvas and pads
   TCanvas *c1 = new TCanvas("c1"," Residuals Graphs",700,900);

   auto *p2 = new TPad("p2","p3",0.,0.,1.,0.3); p2->Draw();
   //pad title
   p2->SetTopMargin(0.001);
   p2->SetBottomMargin(0.3);
   
   p2->SetLogx ();
   p2->SetLogy();

   auto *p1 = new TPad("p1","p1",0.,0.3,1.,1.);  p1->Draw();

   p1->SetBottomMargin(0.001);
   p1->cd();
   p1->SetLogx ();
   p1->SetLogy();

//extract tgrapherrors from root file

//graph model
    TFile *f1 = new TFile("AMS02.root");
    TGraph *gr_ams_new_set = (TGraph*)f1->Get("grProton_vs_Ekn_T0");


//graph data with errors 
    TFile *f2 = new TFile("grAMS02_ProtonHelium_vs_Ekn_vs_Time_DEC2017.root");
    TGraphErrors *gr_ams = (TGraphErrors*)f2->Get("grProton_vs_Ekn_T0");

//draw options
    gr_ams->SetMarkerStyle(20);
    gr_ams->SetMarkerColor(kBlack);
    gr_ams->SetMarkerSize(0.8);

    gr_ams_new_set->SetLineColor(kAzure+1);
    gr_ams_new_set->SetLineWidth(3);

//title and axis  
   gr_ams->GetYaxis()->SetTitle("J(E) [MeV^{ -1} m^{ -2} s^{ -1} sr^{ -1} ]");

//drawing  in canva
   gr_ams->Draw("AP");
   gr_ams_new_set->Draw("P same");

//adding legend
   TLegend *leg = new TLegend(0.15,0.75,0.5,0.85);
   leg->AddEntry(gr_ams_new_set,"ams new ser","lp");
   leg->AddEntry(gr_ams,"ams old set (dec 2017)","lp");
   leg->Draw();

// ratio plot
//get x and y values from TGraphErrors

    int n = gr_ams->GetN();
    double x[n], y1[n], y2[n];

    for (int i=0; i<n; i++) {
        gr_ams->GetPoint(i,x[i],y1[i]);
        y2[i] = gr_ams_new_set->Eval(x[i]);
    }

//draw ratio plot
   p2->cd();
   p2->SetGridy();
   TGraph*r = new TGraph(n); r->SetTitle("");
   for (int i=0; i<n; i++) {
   r->SetPoint(i, x[i], y2[i]/y1[i]);
   }


      r->GetXaxis()->SetLabelSize(0.08);
      r->GetYaxis()->SetLabelSize(0.075);
      r->GetXaxis()->SetTitleSize(0.08);
      r->GetYaxis()->SetTitleSize(0.0875);
      
      r->GetXaxis()->SetTitle("kinetic energy [MeV/n]");
      r->GetYaxis()->SetTitleOffset(0.5);
      r->GetYaxis()->SetTitle("new /old ");
      r->SetMarkerColor(kRed);
      r->SetMarkerStyle(20);
      r->SetMarkerSize(0.8);
      r->Draw("APL");


// write canvas in root file
   TFile *f = new TFile("residuals.root","RECREATE");
   r->Write();
   
   f->Close();

}




