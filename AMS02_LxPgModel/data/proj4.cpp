/*

void ConvertToTGraphErrors() {

int nTimePAMELA_Proton  = 83; // 79+2 with UTTPS off: 81

TGraphAsymmErrors* grPAMELA_ProtonFluxVSEkn_Asymm[nTimePAMELA_Proton];
TGraphErrors* grPAMELA_ProtonFluxVSEkn[nTimePAMELA_Proton];

TFile* inFilePAMELA= new TFile("ssdc_canvas.root","READ");
  
for(int tt=0;tt<nTimePAMELA_Proton;tt++){
      grPAMELA_ProtonFluxVSEkn_Asymm[tt]=(TGraphAsymmErrors*)inFilePAMELA->Get(Form("graph%d",tt+1));
     
  double* x=grPAMELA_ProtonFluxVSEkn_Asymm[tt]->GetX();
  double* y=grPAMELA_ProtonFluxVSEkn_Asymm[tt]->GetY();
  double* exl=grPAMELA_ProtonFluxVSEkn_Asymm[tt]->GetEXlow();
  double* exh=grPAMELA_ProtonFluxVSEkn_Asymm[tt]->GetEXhigh();
  double* eyl=grPAMELA_ProtonFluxVSEkn_Asymm[tt]->GetEYlow();
  double* eyh=grPAMELA_ProtonFluxVSEkn_Asymm[tt]->GetEYhigh();
  
  for(int ii=0;ii<78;ii++){
    double ex=0.5*(exl[ii]+exh[ii]);
    double ey=0.5*(eyl[ii]+eyh[ii]);
    grPAMELA_ProtonFluxVSEkn[tt]->SetPoint(ii, x[ii], y[ii]);
    grPAMELA_ProtonFluxVSEkn[tt]->SetPointError(ii, ex, ey);
    }
  }


//save tgrapherrors in a root file
TFile* outFile= new TFile("PAMELA.root","RECREATE");
for(int tt=0;tt<nTimePAMELA_Proton;tt++){
  grPAMELA_ProtonFluxVSEkn[tt]->SetName(Form(" grPAMELA_ProtonFluxVSEkn_T%d",tt));
  grPAMELA_ProtonFluxVSEkn[tt]->Write();
  }
}

*/

void ConvertToTGraphErrors() {
int nTimePAMELA_Proton  = 83; // 79+2 with UTTPS off: 81

double ex,ey;
double *x,*y,*exl,*exh,*eyl,*eyh;


TGraphAsymmErrors* grPAMELA_ProtonFluxVSEkn_Asymm[nTimePAMELA_Proton];
TGraphErrors* grPAMELA_ProtonFluxVSEkn[nTimePAMELA_Proton];

TFile* inFilePAMELA= new TFile("ssdc_canvas.root","READ");
  
for(int tt=0;tt<nTimePAMELA_Proton;tt++){

      grPAMELA_ProtonFluxVSEkn_Asymm[tt]=(TGraphAsymmErrors*)inFilePAMELA->Get(Form("graph%d",tt+1));
}

for(int tt=0;tt<nTimePAMELA_Proton;tt++){
  /*double* x=grPAMELA_ProtonFluxVSEkn_Asymm[tt]->GetX();
  double* y=grPAMELA_ProtonFluxVSEkn_Asymm[tt]->GetY();
  double* exl=grPAMELA_ProtonFluxVSEkn_Asymm[tt]->GetEXlow();
  double* exh=grPAMELA_ProtonFluxVSEkn_Asymm[tt]->GetEXhigh();
  double* eyl=grPAMELA_ProtonFluxVSEkn_Asymm[tt]->GetEYlow();
  double* eyh=grPAMELA_ProtonFluxVSEkn_Asymm[tt]->GetEYhigh();
  */

  for(int ii=0;ii<78;ii++){
     ex=0.5*(exl[ii]+exh[ii]);
     ey=0.5*(eyl[ii]+eyh[ii]);
      grPAMELA_ProtonFluxVSEkn[0]->SetPoint(ii, 0, 1);
   // grPAMELA_ProtonFluxVSEkn[tt]->SetPoint(ii, x[ii], y[ii]);
    //grPAMELA_ProtonFluxVSEkn[tt]->SetPointError(ii, ex, ey);
    }
  }



}