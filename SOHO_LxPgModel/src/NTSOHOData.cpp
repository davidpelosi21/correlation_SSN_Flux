#include "NTSOHOData.h"

NTSOHOData::NTSOHOData(){
  kTimeUnits= 0; // [ UT | FYR ]
}


NTSOHOData::NTSOHOData(int TimeUnits){
  kTimeUnits= TimeUnits; // [0:UT | 1:FYR]
}


void NTSOHOData::SetAllData(){
  SetProtonData();
}

void NTSOHOData::SetProtonData() {
  // ---- set SOHO data ----
  double xEkn_SOHO[13]={ 0.292, 0.336, 0.387, 0.446, 0.513, 0.591, 0.681, 0.784, 0.903, 1.04, 1.198, 1.38, 1.589 }; 
  double xErrEkn_SOHO[13]={ 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.};
  
  double xTime_SOHO[20]={ 8.04683e+08, 8.36264e+08, 8.67845e+08, 8.99381e+08, 9.30917e+08, 9.62496e+08, 9.94075e+08, 1.02561e+09, 1.05715e+09, 1.08873e+09, 1.12031e+09, 1.15184e+09, 1.18338e+09, 1.21496e+09, 1.24654e+09, 1.27807e+09, 1.30961e+09, 1.34119e+09, 1.37277e+09, 1.4043e+09 }; 

  // double xErrTime_SOHO[20];
  // for(int tt=0;tt<20;tt++) xErrTime_SOHO[tt]= (3.154e+7)/2.0001; // 6-month

  double ErrTime_SOHO = (3.154e+7)/2.0001; // time error: 6 months
  for(int tt=0;tt<20;tt++) eTimeSOHO_ProtonFlux[tt]= ErrTime_SOHO;

  
  double valFluxVSTimeVSEkn_SOHO[20][13]= {
    { 1960, 1830, 1830, 1690, 1700, 1400, 1350, 1190, 1130, 1010, 940, 760, 600,  }, 
    { 2000, 1900, 1900, 1710, 1740, 1410, 1360, 1210, 1140, 1010, 930, 760, 590,  }, 
    { 2120, 2020, 2020, 1800, 1830, 1490, 1430, 1260, 1180, 1050, 960, 780, 610,  }, 
    { 1840, 1710, 1710, 1540, 1580, 1300, 1250, 1110, 1030, 960, 870, 720, 570,  }, 
    { 1220, 1150, 1180, 1070, 1110, 920, 900, 820, 800, 730, 680, 570, 460,  }, 
    { 580, 560, 580, 510, 550, 490, 500, 460, 450, 460, 430, 380, 310,  }, 
    { 530, 510, 510, 490, 540, 470, 460, 460, 480, 480, 450, 410, 340,  }, 
    { 610, 590, 560, 540, 590, 520, 550, 510, 500, 490, 460, 420, 330,  }, 
    { 560, 520, 530, 500, 550, 480, 490, 470, 450, 460, 430, 370, 310,  }, 
    { 680, 680, 690, 640, 670, 590, 600, 560, 550, 540, 500, 430, 350,  }, 
    { 940, 940, 950, 900, 910, 750, 780, 730, 720, 660, 600, 520, 410,  }, 
    { 1270, 1220, 1240, 1130, 1170, 980, 970, 860, 840, 760, 700, 580, 450,  }, 
    { 1670, 1590, 1600, 1440, 1470, 1220, 1170, 1040, 990, 880, 800, 660, 510,  }, 
    { 1860, 1760, 1760, 1590, 1600, 1320, 1260, 1120, 1050, 930, 850, 680, 530,  }, 
    { 2270, 2150, 2130, 1900, 1910, 1540, 1470, 1290, 1200, 1060, 940, 760, 570,  }, 
    { 1960, 1850, 1830, 1630, 1660, 1350, 1290, 1140, 1060, 970, 850, 680, 520,  }, 
    { 1250, 1260, 1240, 1090, 1160, 940, 920, 830, 800, 750, 650, 520, 410,  }, 
    { 820, 770, 820, 710, 780, 650, 650, 590, 560, 530, 480, 390, 300,  }, 
    { 790, 760, 760, 710, 750, 640, 650, 620, 590, 580, 490, 410, 320,  }, 
    { 740, 720, 720, 670, 730, 640, 640, 610, 570, 580, 460, 380, 300,  }, 
  };
    
    double errFluxVSTimeVSEkn_SOHO[20][13]= {
      { 294.653, 275.109, 275.109, 254.063, 255.566, 210.466, 202.95, 178.896, 169.876, 151.836, 141.313, 114.253, 90.1998,  }, 
      { 300.666, 285.633, 285.633, 257.069, 261.579, 211.969, 204.453, 181.903, 171.38, 151.836, 139.81, 114.253, 88.6964,  }, 
      { 318.706, 303.673, 303.673, 270.599, 275.109, 223.996, 214.976, 189.42, 177.393, 157.85, 144.32, 117.26, 91.7031,  }, 
      { 276.613, 257.069, 257.069, 231.513, 237.526, 195.433, 187.916, 166.87, 154.843, 144.32, 130.79, 108.24, 85.6898,  }, 
      { 183.406, 172.883, 177.393, 160.856, 166.87, 138.306, 135.3, 123.273, 120.266, 109.743, 102.226, 85.6898, 69.1532,  }, 
      { 87.1931, 84.1865, 87.1931, 76.6698, 82.6831, 73.6632, 75.1665, 69.1532, 67.6498, 69.1532, 64.6432, 57.1265, 46.6032,  }, 
      { 79.6765, 76.6698, 76.6698, 73.6632, 81.1798, 70.6565, 69.1532, 69.1532, 72.1598, 72.1598, 67.6498, 61.6365, 51.1132,  }, 
      { 91.7031, 88.6964, 84.1865, 81.1798, 88.6964, 78.1731, 82.6831, 76.6698, 75.1665, 73.6632, 69.1532, 63.1398, 49.6099,  }, 
      { 84.1865, 78.1731, 79.6765, 75.1665, 82.6831, 72.1598, 73.6632, 70.6565, 67.6498, 69.1532, 64.6432, 55.6232, 46.6032,  }, 
      { 102.226, 102.226, 103.73, 96.2131, 100.723, 88.6964, 90.1998, 84.1865, 82.6831, 81.1798, 75.1665, 64.6432, 52.6165,  }, 
      { 141.313, 141.313, 142.816, 135.3, 136.803, 112.75, 117.26, 109.743, 108.24, 99.2198, 90.1998, 78.1731, 61.6365,  }, 
      { 190.923, 183.406, 186.413, 169.876, 175.89, 147.326, 145.823, 129.286, 126.28, 114.253, 105.233, 87.1931, 67.6498,  }, 
      { 251.056, 239.029, 240.533, 216.479, 220.989, 183.406, 175.89, 156.346, 148.83, 132.293, 120.266, 99.2198, 76.6698,  }, 
      { 279.619, 264.586, 264.586, 239.029, 240.533, 198.44, 189.42, 168.373, 157.85, 139.81, 127.783, 102.226, 79.6765,  }, 
      { 341.256, 323.216, 320.209, 285.633, 287.136, 231.513, 220.989, 193.93, 180.4, 159.353, 141.313, 114.253, 85.6898,  }, 
      { 294.653, 278.116, 275.109, 245.043, 249.553, 202.95, 193.93, 171.38, 159.353, 145.823, 127.783, 102.226, 78.1731,  }, 
      { 187.916, 189.42, 186.413, 163.863, 174.386, 141.313, 138.306, 124.776, 120.266, 112.75, 97.7164, 78.1731, 61.6365,  }, 
      { 123.273, 115.756, 123.273, 106.736, 117.26, 97.7164, 97.7164, 88.6964, 84.1865, 79.6765, 72.1598, 58.6299, 45.0999,  }, 
      { 118.763, 114.253, 114.253, 106.736, 112.75, 96.2131, 97.7164, 93.2064, 88.6964, 87.1931, 73.6632, 61.6365, 48.1065,  }, 
      { 111.246, 108.24, 108.24, 100.723, 109.743, 96.2131, 96.2131, 91.7031, 85.6898, 87.1931, 69.1532, 57.1265, 45.0999,  }, 
    };


    // ---- CORRECT SOHO DATA / RENORM ---- PORK
    // for(int tt=0;tt<20;tt++){
    //   for(int ee=0;ee<13;ee++){
    // 	valFluxVSTimeVSEkn_SOHO[tt][ee] *= 1.10;
    //   }
    // }

    
    for(int tt=0;tt<20;tt++){
      xTimeSOHO_ProtonFlux[tt]= xTime_SOHO[tt]; // TIME GRID
      grSOHO_ProtonFluxVSEkn[tt]= new TGraphErrors(13, xEkn_SOHO, valFluxVSTimeVSEkn_SOHO[tt], xErrEkn_SOHO, errFluxVSTimeVSEkn_SOHO[tt]);
      grSOHO_ProtonFluxVSEkn[tt]->SetName(Form("grSOHO_ProtonFluxVSEkn_T%d",tt));
      grSOHO_ProtonFluxVSEkn[tt]->SetMarkerStyle(25);
      grSOHO_ProtonFluxVSEkn[tt]->SetMarkerSize(1.0);
      grSOHO_ProtonFluxVSEkn[tt]->SetMarkerColor(kGreen+3);
      grSOHO_ProtonFluxVSEkn[tt]->SetLineColor(kGreen+3);
    }


    // ---- SOHO flux VS time ----
    double valFluxVSEknVSTime_SOHO[13][20];
    double errFluxVSEknVSTime_SOHO[13][20];

    for(int tt=0;tt<20;tt++){
      for(int ee=0;ee<13;ee++){
	valFluxVSEknVSTime_SOHO[ee][tt]=valFluxVSTimeVSEkn_SOHO[tt][ee];
	errFluxVSEknVSTime_SOHO[ee][tt]=errFluxVSTimeVSEkn_SOHO[tt][ee];
      }
    }



    for(int ee=0;ee<13;ee++){
      xEknSOHO_ProtonFlux[ee]=xEkn_SOHO[ee]; // EKN grid
      grSOHO_ProtonFluxVSTime[ee]= new TGraphErrors(20, xTime_SOHO, valFluxVSEknVSTime_SOHO[ee], eTimeSOHO_ProtonFlux, errFluxVSEknVSTime_SOHO[ee]);
      grSOHO_ProtonFluxVSTime[ee]->SetName(Form("grSOHO_ProtonFluxVSTime_E%d",ee));
      grSOHO_ProtonFluxVSTime[ee]->SetMarkerStyle(25);
      grSOHO_ProtonFluxVSTime[ee]->SetMarkerSize(1.0);
      grSOHO_ProtonFluxVSTime[ee]->SetMarkerColor(kGreen+3);
      grSOHO_ProtonFluxVSTime[ee]->SetLineColor(kGreen+3);
    }


    // ---- convert time units UT - FYR ----

    if(kTimeUnits==1){
      
      // convert graph
      for(int ee=0;ee<nEknSOHO_Proton;ee++){ 
	ConvertUnixTime2FractYear( grSOHO_ProtonFluxVSTime[ee] );
      }
      
      // update time array
      double* xTimeSOHO= (double*)grSOHO_ProtonFluxVSTime[0]->GetX();
      for(int tt=0;tt<nTimeSOHO_Proton;tt++) xTimeSOHO_ProtonFlux[tt]=xTimeSOHO[tt];
      for(int tt=0;tt<nTimeSOHO_Proton;tt++){ // tutti uguali
	eTimeSOHO_ProtonFlux[tt]= 0.45*(xTimeSOHO_ProtonFlux[1]-xTimeSOHO_ProtonFlux[0]);
      }
    }
    
}









// ---- auxiliary functions ----

void NTSOHOData::ScaleGraphErrors(TGraphErrors* gr, double ScaleValue){
  int NP= gr->GetN();
  for(int ii=0;ii<NP;ii++){
    double Ekn, Val, ErrX, ErrY;
    gr->GetPoint(ii, Ekn, Val);
    ErrX=gr->GetErrorX(ii);
    ErrY=gr->GetErrorY(ii);
    Val  *= ScaleValue;
    ErrY *= ScaleValue;
    gr->SetPoint(ii, Ekn, Val);  
    gr->SetPointError(ii, ErrX, ErrY);  
  }
}



void NTSOHOData::CleanXErrors(TGraphErrors* gr){
  int NP=gr->GetN();
  for(int ii=0;ii<NP;ii++){
    double e_ekn=0.;
    double e_val=gr->GetErrorY(ii);
    gr->SetPointError(ii, e_ekn,e_val); // SOMETHING WRONG HERE...?
  }
}



void NTSOHOData::RemoveZero(TGraphErrors* gr){  
  for(int ii=gr->GetN()-1;ii>=0;ii--){
    if( !(gr->GetY()[ii] > 0) ){
      gr->RemovePoint(ii);
    }
  }
}


void NTSOHOData::ConvertUnixTime2FractYear(TGraphErrors* gr){

  for(int tt=0;tt<gr->GetN();tt++){
    UInt_t UT = (UInt_t)(gr->GetX()[tt]);

    // --- extract date ---
    TDatime date(UT);
    int YY = date.GetYear();
    int MM = date.GetMonth();
    int DD = date.GetDay();

    // --- day of the year, approximately ---
    // double DOY = (MM-1.0)*30.416667 + DD*1.0 + 0.4; 

    // --- extract day of the year ---
    double DOY= 0;
    if(MM>=1) DOY = (double)DD;
    if(MM>=2) DOY+= 31; // add Jan
    if(MM>=3) DOY+= 28; // add Feb
    if(MM>=4) DOY+= 31; // add Mar
    if(MM>=5) DOY+= 30; // add Apr
    if(MM>=6) DOY+= 31; // add May
    if(MM>=7) DOY+= 30; // add Jun
    if(MM>=8) DOY+= 31; // add Jul
    if(MM>=9) DOY+= 31; // add Aug
    if(MM>=10)DOY+= 30; // add Sep
    if(MM>=11)DOY+= 31; // add Oct
    if(MM>=12)DOY+= 30; // add Nov
    
    // --- compute fractional year ---
    double FD = DOY/365.; // fraction of year
    if(FD>=1.)FD=0.999999999;
    double FYR = YY*1. + FD;
    //cout<<"Y: "<<YY<<"   M: "<<MM<<"    D: "<<DD<<"       DOY: "<<DOY<<"    FYR: "<<FYR<<endl;

    gr->GetX()[tt]= FYR;
    gr->GetEX()[tt]= 0.;
  }
}
