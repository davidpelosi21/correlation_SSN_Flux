
void residuals_N0_exp() {
double NPAMELA;
double NAMS02 = 1.;
double residuo = 0.0;

//create an array of integers
//int br_pamela[26] = {57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82};
//int br_ams02[26] = {0,1,4,6,7,8,9,10,12,14,17,18,19,20,21,22,23,24,28,29,30,31,32,34,36};

int br_pamela[10] = {57,58,59,60,61,62,63,64,65,66};

int br_ams02[10] = {0,1,4,6,7,8,9,10,12,14};

int i = 6;
//TString root_file_path_AMS02 = Form("/Users/davidpelosi/Desktop/Correlation/AMS02_NewSet_LxPgModel/OUT/TGraph_Flux_data_LIS4_Z1_A1_BR%d_DCPotgeiter.root",br_ams02[i]);
TString root_file_path_AMS02 = Form("/Users/davidpelosi/Desktop/Correlation/AMS02_LxPgModel/OUT/TGraph_Flux_data_LIS4_Z1_A1_BR%d_DCPotgeiter.root",br_ams02[i]);
TString root_file_path_PAMELA = Form("/Users/davidpelosi/Desktop/Correlation/PAMELA_LxPgModel/OUT/TGraph_Flux_data_LIS4_Z1_A1_BR%d_DCPotgeiter.root",br_pamela[i]);


//graph AMS02 with errors 
    TFile *f1 = new TFile(root_file_path_AMS02);
    TGraphErrors *gr_AMS02 = (TGraphErrors*)f1->Get(Form("grProton_vs_Ekn_T%d",br_ams02[i]));


//graph PAMELA data with errors 
    TFile *f2 = new TFile(root_file_path_PAMELA);
    TGraphErrors *gr_PAMELA = (TGraphErrors*)f2->Get(Form("graph%d",br_pamela[i]));

    
//create canvas and pads
   TCanvas *c1 = new TCanvas("c1"," Residuals Graphs",700,900);

   auto *p2 = new TPad("p2","p3",0.,0.,1.,0.3); p2->Draw();
   //pad title
   p2->SetTopMargin(0.001);
   p2->SetBottomMargin(0.3);
   
   p2->SetLogx ();
   //p2->SetLogy();

   auto *p1 = new TPad("p1","p1",0.,0.3,1.,1.);  p1->Draw();

   p1->SetBottomMargin(0.001);

    int n = gr_AMS02->GetN();
    double x[n], y1[n], y2[n];


    for (int i=0; i<n; i++) {
        gr_AMS02->GetPoint(i,x[i],y1[i]);
        y2[i] = gr_PAMELA->Eval(x[i]);
    }

//draw ratio plot
   p2->cd();
   p2->SetGridy();
   TGraph*r = new TGraph(n); r->SetTitle("");


NPAMELA = 1;
double step= 0.0001;
double bestNPAMELA;
double bestresiduo = 100000;
   while (NPAMELA <= 1.15) {
    residuo = 0.0; 
    
for (int i=0; i<n; i++) {
    if (x[i]>490 && x[i]<9644) {
     residuo = residuo + abs((NAMS02*y1[i])-(NPAMELA*y2[i]));
     //residuo = residuo + (NAMS02*y1[i])/(NPAMELA*y2[i]);
    }
}
   
if (residuo < bestresiduo) {
    bestresiduo = residuo;
    bestNPAMELA = NPAMELA;
}

     NPAMELA = NPAMELA + step; 
}
   
    cout << "best NPAMELA: " << bestNPAMELA << endl;
    cout<< "best residuo: " << bestresiduo << endl;


int ii=0;
    for (int i=0; i<n; i++) {
   // if (x[i]>4.e3)
    if (x[i]>490 && x[i]<9644)

    {
     //r->SetPoint(i, x[i], (NPAMELA*y2[i])/(NAMS02*y1[i]));
     //r->SetPoint(ii, x[i], abs(y1[i]-y2[i]));
      r->SetPoint(ii, x[i], abs(y1[i]-(NPAMELA*y2[i])));
     ii++;

    }
   }

      r->GetXaxis()->SetLabelSize(0.08);
      r->GetYaxis()->SetLabelSize(0.075);
      r->GetXaxis()->SetTitleSize(0.08);
      r->GetYaxis()->SetTitleSize(0.0875);
      
      r->GetXaxis()->SetTitle("kinetic energy [MeV/n]");
      r->GetYaxis()->SetTitleOffset(0.5);
      //r->GetYaxis()->SetTitle("ams02/pamela");
      r->GetYaxis()->SetTitle("|Flux(ams) - Flux(pamela) |");
      r->SetMarkerColor(kRed);
      r->SetMarkerStyle(20);
      r->SetMarkerSize(0.8);



      r->Draw("AP");



//pad 1
   p1->cd();
   p1->SetLogx ();
   p1->SetLogy();

//extract tgrapherrors from root file

     int j=0; int k=0;
//copy the content of a tgrapherrors in a tgrapherrors point by point
    TGraphErrors *gr_PAMELA_copy = new TGraphErrors();
    for (int i=0; i<gr_PAMELA->GetN(); i++) {
        double x,y;
        gr_PAMELA->GetPoint(i,x,y);
            if (x>490 && x<9644){
        gr_PAMELA_copy->SetPoint(j,x, 0.96*y);
        //gr_PAMELA_copy->SetPoint(j,x,y);
             j++;
            }
    }

    TGraphErrors *gr_AMS02_copy = new TGraphErrors();
    for (int i=0; i<gr_AMS02->GetN(); i++) {
        double x,y;
        gr_AMS02->GetPoint(i,x,y);
        if (x>490 && x<9644){
        gr_AMS02_copy->SetPoint(k,x,y);
        k++;
        }

    }

//draw options
    gr_PAMELA_copy->SetMarkerStyle(20);
    gr_PAMELA_copy->SetMarkerColor(kBlack);
    gr_PAMELA_copy->SetMarkerSize(0.8);

    gr_AMS02_copy->SetMarkerStyle(20);
    gr_AMS02_copy->SetMarkerColor(kAzure+2);
    gr_AMS02_copy->SetMarkerSize(0.8);



//title and axis  
   gr_PAMELA_copy->GetYaxis()->SetTitle("J(E) [MeV^{ -1} m^{ -2} s^{ -1} sr^{ -1} ]");

//drawing  in canva
   gr_PAMELA_copy->Draw("AP");
   gr_AMS02_copy->Draw("P same");

//adding legend
   TLegend *leg = new TLegend(0.15,0.75,0.5,0.85);
   leg->AddEntry(gr_AMS02_copy,"ams02","lp");
   leg->AddEntry(gr_PAMELA_copy,"pamela","lp");
   leg->Draw();

// ratio plot
//get x and y values from TGraphErrors


// write canvas in root file
   TFile *f = new TFile("residuals.root","RECREATE");
   r->Write();
   
   
   f->Close();
}



