void residuals_N0_exp() {
double NPAMELA;
double residuo = 0.0;


NPAMELA = 1.01;

//create an array of integers
int br_pamela[26] = {57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82};
int br_ams02[26] = {0,1,4,6,7,8,9,10,12,14,17,18,19,20,21,22,23,24,28,29,30,31,32,34,36};

//int br_pamela[10] = {57,58,59,60,61,62,63,64,65,66};

//int br_ams02[10] = {0,1,4,6,7,8,9,10,12,14};

//create canvas and pads
//graph PAMELA data with errors 

for(int j=0; j <26; j++){

    TString root_file_path_AMS02 = Form("/Users/davidpelosi/Desktop/Correlation/AMS02_NewSet_LxPgModel/OUT/TGraph_Flux_data_LIS4_Z1_A1_BR%d_DCPotgeiter.root",br_ams02[j]);
//TString root_file_path_AMS02 = Form("/Users/davidpelosi/Desktop/Correlation/AMS02_LxPgModel/OUT/TGraph_Flux_data_LIS4_Z1_A1_BR%d_DCPotgeiter.root",br_ams02[i]);
TString root_file_path_PAMELA = Form("/Users/davidpelosi/Desktop/Correlation/PAMELA_LxPgModel/OUT/TGraph_Flux_data_LIS4_Z1_A1_BR%d_DCPotgeiter.root",br_pamela[j]);

//create canvas and pads
    TFile *f1 = new TFile(root_file_path_AMS02);
    TFile *f2 = new TFile(root_file_path_PAMELA);
  
//graph AMS02 with errors 
    TGraphErrors *gr_AMS02 = (TGraphErrors*)f1->Get(Form("grProton_vs_Ekn_T%d",br_ams02[j]));
//graph PAMELA data with errors 
    TGraphErrors *gr_PAMELA = (TGraphErrors*)f2->Get(Form("graph%d",br_pamela[j]));




    int n = gr_AMS02->GetN();
    double x[n], y1[n], y2[n];

    for (int i=0; i<n; i++) {
        gr_AMS02->GetPoint(i,x[i],y1[i]);
        y2[i] = gr_PAMELA->Eval(x[i]);
    }




for (int i=0; i<n; i++) {
    if (x[i]>490 && x[i]<9644) {
     residuo = residuo + abs((y1[i])-(NPAMELA*y2[i]));
     //residuo = residuo + (NAMS02*y1[i])/(NPAMELA*y2[i]);
    }
}

}
   
    cout << " NPAMELA: " << NPAMELA << endl;
    cout<< "tot residuo: " << residuo/26 << endl;

}



