
void residuals() {
TString root_file_path = "TGraph_Flux_data_LIS4_Z1_A1_BR30_DCStandard.root";
TString root_file_path2 = "residuals.root";


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
    TFile *f1 = new TFile(root_file_path);
    TGraph *gr_model = (TGraph*)f1->Get("Graph");


//graph data with errors 
    TFile *f2 = new TFile(root_file_path2);
    TGraphErrors *gr_data = (TGraphErrors*)f2->Get("residuals");

//draw options
    gr_model->SetLineColor(kAzure+1);
    gr_model->SetLineWidth(3);

//title and axis  

//drawing  in canva

    TGraphErrors *gr_data2 = new TGraphErrors();
    //clone
    for (int i=0; i<gr_data->GetN(); i++) {
        double x,y;
        gr_data->GetPoint(i,x,y);
        gr_data2->SetPoint(i,x,y);
        cout<<"x = "<<x<<" y = "<<y<<endl;
        gr_data2->SetPointError(i,0,0);
    }
    
    gr_data2->SetMarkerStyle(20);
    gr_data2->SetMarkerColor(kBlack);
    gr_data2->SetMarkerSize(0.8);
    gr_data2->GetYaxis()->SetTitle("J(E) [MeV^{ -1} m^{ -2} s^{ -1} sr^{ -1} ]");

   gr_data2->Draw("AP");
   gr_model->Draw("L same");

//adding legend
   TLegend *leg = new TLegend(0.15,0.75,0.5,0.85);
   leg->AddEntry(gr_model,"model","lp");
   leg->AddEntry(gr_data,"exp data","lp");
   leg->Draw();

// ratio plot
//get x and y values from TGraphErrors

    int n = gr_data->GetN();
    double x[n], y1[n], y2[n];

    for (int i=0; i<n; i++) {
        gr_data->GetPoint(i,x[i],y1[i]);
        y2[i] = gr_model->Eval(x[i]);
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
      r->GetYaxis()->SetTitle("model/data");
      r->SetMarkerColor(kRed);
      r->SetMarkerStyle(20);
      r->SetMarkerSize(0.8);
      r->Draw("APL");


// write canvas in root file
   TFile *f = new TFile("residuals_final.root","RECREATE");
   r->Write();
   gr_data2->Write();
   f->Close();


}




