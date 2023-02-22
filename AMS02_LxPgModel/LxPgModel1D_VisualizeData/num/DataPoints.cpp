#include "DataPoints.h"
#include <fstream>
#include <cstdio>
#include <iomanip>

DataPoints::DataPoints(int n, double* xref, double* yref): N(n){
	x = new double[N];
	y = new double[N];


	for(int j=0; j<N; j++){
		x[j]=xref[j];
		y[j]=yref[j];
	}
}

DataPoints::DataPoints(const DataPoints& D){
	x = new double[D.N];
	y = new double[D.N];


	for(int j=0; j<D.N; j++){
		x[j]=D.x[j];
		y[j]=D.y[j];
	}
}

DataPoints::~DataPoints(){
	if(x && y){
		delete[] x;
		delete[] y;
	} 

	if(g)
		delete g;
}

void DataPoints::Print(string namefile){
	ofstream F(namefile.c_str());

	for(int i=0; i<N; i++){
		F << setw(3) << x[i] << setw(14) << y[i] << endl;
	}

	F.close();
}

