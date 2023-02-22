void N0_PAMELA() {
double NPAMELA = 1;
double NAMS02 = 1.;
double residuo = 0.0;

//create an array of integers
int br_pamela[26] = {57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82};
int br_ams02[26] = {0,1,4,6,7,8,9,10,12,14,17,18,19,20,21,22,23,24,28,29,30,31,32,34,36};

//int br_pamela[10] = {57,58,59,60,61,62,63,64,65,66};

//int br_ams02[10] = {0,1,4,6,7,8,9,10,12,14};

double bestNPAMELA_tot[25];

for(int i = 0; i<25; i++){
    
TString root_file_path_AMS02 = Form("/Users/davidpelosi/Desktop/Correlation/AMS02_NewSet_LxPgModel/OUT/TGraph_Flux_data_LIS4_Z1_A1_BR%d_DCPotgeiter.root",br_ams02[i]);
//TString root_file_path_AMS02 = Form("/Users/davidpelosi/Desktop/Correlation/AMS02_LxPgModel/OUT/TGraph_Flux_data_LIS4_Z1_A1_BR%d_DCPotgeiter.root",br_ams02[i]);
TString root_file_path_PAMELA = Form("/Users/davidpelosi/Desktop/Correlation/PAMELA_LxPgModel/OUT/TGraph_Flux_data_LIS4_Z1_A1_BR%d_DCPotgeiter.root",br_pamela[i]);


//extract tgrapherrors from root file
//graph AMS02 with errors 
    TFile *f1 = new TFile(root_file_path_AMS02);
    TGraph *gr_AMS02 = (TGraph*)f1->Get(Form("grProton_vs_Ekn_T%d",br_ams02[i]));


//graph PAMELA data with errors 
    TFile *f2 = new TFile(root_file_path_PAMELA);
    TGraphErrors *gr_PAMELA = (TGraphErrors*)f2->Get(Form("graph%d",br_pamela[i]));


// ratio plot
//get x and y values from TGraphErrors

    int n = gr_AMS02->GetN();
    double x[n], y1[n], y2[n];


    for (int i=0; i<n; i++) {
        gr_AMS02->GetPoint(i,x[i],y1[i]);
        y2[i] = gr_PAMELA->Eval(x[i]);
    }



NPAMELA = 0.85;
double step= 0.0001;
double bestNPAMELA;
double bestresiduo = 100000;


   while (NPAMELA <= 4) {
    residuo = 0.0; 
    
for (int i=0; i<n; i++) {
    //if (x[i]>490 && x[i]<9644) {
     residuo = residuo + abs((NAMS02*y1[i])-(NPAMELA*y2[i]));
   // }
}
   
if (residuo < bestresiduo) {
    bestresiduo = residuo;
    bestNPAMELA = NPAMELA;
}

     NPAMELA = NPAMELA + step; 
}
   
cout << "best NPAMELA: " << bestNPAMELA << endl;

bestNPAMELA_tot[i] = bestNPAMELA;

}

//print mean of bestNPAMELA_tot
double sum = 0;
for (int i = 0; i<25; i++){
    sum = sum + bestNPAMELA_tot[i];
}

cout << "mean of bestNPAMELA_tot: " << sum/25 << endl;


}


