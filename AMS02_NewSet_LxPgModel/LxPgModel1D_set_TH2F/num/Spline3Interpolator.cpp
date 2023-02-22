#include "Spline3Interpolator.h"
#include "EqSolver.h"
#include "TLegend.h"
#include <cmath>
#include "Vec.h"
#include <vector>
#include <cstdio>
using namespace std;

Spline3Interpolator::Spline3Interpolator(int fN, double *fx, double *fy, TF1* fF0) : DataPoints(fN,fx,fy) {
	F0=fF0;
	K = new double[N];
	SetCurvatureLines(); //define segment interpolators
	FInterpolator = new TF1("FInterpolator", this, &Spline3Interpolator::fInterpolator, x[0]-(x[N-1]-x[0])/20,x[N-1]+(x[N-1]-x[0])/20, 0, "Spline3Interpolator", "fInterpolator");
	D1FInterpolator = new TF1("D1FInterpolator", this, &Spline3Interpolator::fDerivator, x[0]-0.1,x[N-1]+0.1, 0, "Spline3Interpolator", "fDerivator");
	D2FInterpolator = new TF1("D2FInterpolator", this, &Spline3Interpolator::fSecDerivator, x[0],x[N-1], 0, "Spline3Interpolator", "fSecDerivator");
    FInterpolatorInverse = new TF1("FInterpolatorInverse", this, &Spline3Interpolator::fInterpolatorInverse, x[0]-0.1,x[N-1]+0.1, 0, "Spline3Interpolator", "fInterpolatorInverse");
	
}

Spline3Interpolator::Spline3Interpolator(const Spline3Interpolator& S){
	N=S.GetN();
	x = new double[N];
	y = new double[N];


	for(int j=0; j<N; j++){
		x[j]=S.x[j];
		y[j]=S.y[j];
	}

	if(S.F0)
	F0= new TF1(*(S.F0));
	FInterpolator = new TF1("FInterpolator", this, &Spline3Interpolator::fInterpolator,
	x[0]-1,x[N-1]+1, 0, "Spline3Interpolator", "fInterpolator");
	D1FInterpolator = new TF1("D1FInterpolator", this, &Spline3Interpolator::fDerivator,
	x[0]-0.1,x[N-1]+0.1, 0, "Spline3Interpolator", "fDerivator");
	D2FInterpolator = new TF1("D2FInterpolator", this, &Spline3Interpolator::fSecDerivator,
	x[0],x[N-1], 0, "Spline3Interpolator", "fSecDerivator");
	K = new double[N];
	SetCurvatureLines(); 
}

Spline3Interpolator::Spline3Interpolator(const DataPoints& D){
	N=D.GetN();
	x = new double[N];
	y = new double[N];


	for(int j=0; j<N; j++){
		x[j]=D.GetX()[j];
		y[j]=D.GetY()[j];
	}

	FInterpolator = new TF1("FInterpolator", this, &Spline3Interpolator::fInterpolator,
	x[0]-1,x[N-1]+1, 0, "Spline3Interpolator", "fInterpolator");
	D1FInterpolator = new TF1("D1FInterpolator", this, &Spline3Interpolator::fDerivator,
	x[0]-1,x[N-1]+1, 0, "Spline3Interpolator", "fDerivator");
	D2FInterpolator = new TF1("D2FInterpolator", this, &Spline3Interpolator::fSecDerivator,
	x[0],x[N-1], 0, "Spline3Interpolator", "fSecDerivator");
	K = new double[N];
	SetCurvatureLines(); 
}

Spline3Interpolator::~Spline3Interpolator(){

	if(F0)
		delete F0;
	if(K){
	delete[] K;
		delete D1FInterpolator;
		delete D2FInterpolator;
		delete FInterpolator;
	}
	
} 

void Spline3Interpolator::SetCurvatureLines() {
	// define tri-diagonal matrix and array of constants
	double* v = new double[N];
	vector<Vec> aux;
	double* b = new double[N];

	for(int i=0; i<N-2; i++){

		for(int j=0; j<N-2; j++){
			if(j==i)
				v[j]=2*(x[i]-x[i+2]);
			else if(j==i-1)
				v[j]=x[i]-x[i+1];
			else if(j==i+1)
				v[j]=x[i+1]-x[i+2];
			else
				v[j]=0;
		}
		aux.push_back(Vec(N,v));

		b[i]=6*((y[i]-y[i+1])/(x[i]-x[i+1])-(y[i+1]-y[i+2])/(x[i+1]-x[i+2]));
	}
	// solve system and get the 2nd derivative coefficients
	Vec solution(N);

	EqSolver S(FCmatrixFull(aux),Vec(N,b));

	solution = S.LUdecompositionSolver(); 

	// store coeffs on internal array K
	K[0]=0;
	for(int i=1; i<N-1; i++)
		K[i]=solution[i-1];

	K[N-1]=0;
	delete[] v;
	delete[] b;
}

double Spline3Interpolator::Interpolate(double fx) {
	// detect in which segment is x
	int i=0;
	for (i=0; i<N; i++) {
		if ((fx-x[i])<=0.) break;
		} //upper bound returned
	if (i==0) // out of range
		i++;
	if(i==N)
		i--;

	//retrieve segment interpolator and return function value
	 return ( (K[i-1]/6)*((pow(fx-x[i],3))/(x[i-1]-x[i])-(fx-x[i])*(x[i-1]-x[i]))-(K[i]/6)*((pow(fx-x[i-1], 3))/(x[i-1]-x[i])-(fx-x[i-1])*(x[i-1]-x[i]))+(y[i-1]*(fx-x[i])-y[i]*(fx-x[i-1]))/(x[i-1]-x[i]));
}

	double Spline3Interpolator::Derivate(double fx) {
	// detect in which segment is x
	int i=0;
	for (i=0; i<N; i++) {
		if ((fx-x[i])<=0.) break;
		} //upper bound returned
	if (i==0) // out of range
		i++;
	if(i==N)
		i--;

	//retrieve segment interpolator and return function value
	 return ((K[i-1]/6)*(3*(fx-x[i])*(fx-x[i])/(x[i-1]-x[i])-(x[i-1]-x[i]))-(K[i]/6)*(3*(fx-x[i-1])*(fx-x[i-1])/(x[i-1]-x[i])-(x[i-1]-x[i]))+(y[i-1]-y[i])/(x[i-1]-x[i]));
}

double  Spline3Interpolator::InterpolateInverse(double fx){ return (1/FInterpolator->Eval(fx));}

double Spline3Interpolator::SecDerivate(double fx) {
	// detect in which segment is x
	int i=0;
	for (i=0; i<N; i++) {
		if ((fx-x[i])<=0.) break;
		} //upper bound returned
	if (i==0) // out of range
		i++;
	if(i==N)
		i--;
	//retrieve segment interpolator and return function value
	return ((K[i-1]*(fx-x[i])-K[i]*(fx-x[i-1]))/(x[i-1]-x[i]));
}


TF1* Spline3Interpolator::GetF0() const { return F0;}
TF1* Spline3Interpolator::GetInterpolationFunction() const { return FInterpolator;}
TF1* Spline3Interpolator::GetInterpolationInverseFunction() const { return FInterpolatorInverse;}
TF1* Spline3Interpolator::GetInterpolationFunctionDerivative() const{return D1FInterpolator;}
TF1* Spline3Interpolator::GetInterpolationFunctionSecDerivative() const{ return D2FInterpolator;}

