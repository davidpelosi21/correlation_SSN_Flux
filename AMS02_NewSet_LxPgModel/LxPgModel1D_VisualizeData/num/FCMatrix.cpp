#include "FCMatrix.h"
#include <cstdio>


FCmatrix::FCmatrix(double** fm, int lin, int col){
    if(fm){
        for(int i =0; i<lin; i++)
            M.push_back(Vec(col,fm[i]));
    }
    
    else{
        for(int i =0; i<lin; i++)
            M.push_back(Vec(col,0.));
    }
}

FCmatrix::FCmatrix(double* fm, int lin, int col){
    double aux[col];
    
    for(int i =0; i<lin; i++){
        for(int j=0; j<col; j++){
            aux[j]=fm[j+col*i];
        }
        M.push_back(Vec(col,aux));
    }
}

FCmatrix::FCmatrix(vector<Vec> aux){
    for(int i=0; i<aux.size();i++)
        M.push_back(aux[i]);
}

FCmatrix::FCmatrix(const FCmatrix& mbase){
    M.clear();
    
    for(int i=0; i<mbase.M.size(); i++)
        M.push_back(mbase.M[i]);
}

FCmatrix::~FCmatrix(){;}

