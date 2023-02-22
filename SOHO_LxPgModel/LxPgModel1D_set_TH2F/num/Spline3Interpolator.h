#ifndef __Spline3Interpolator__
#define __Spline3Interpolator__
#include "DataPoints.h"
#include <stdio.h>
#include <fstream>
#include "TGraph.h"
#include "TF1.h"

class Spline3Interpolator : public DataPoints {
public:
    Spline3Interpolator(int N=0, double *x=NULL, double *y=NULL, TF1* fF0=NULL);
    Spline3Interpolator(const Spline3Interpolator&);
    Spline3Interpolator(const DataPoints&);
    ~Spline3Interpolator();
    double Interpolate(double x);
    double Derivate(double x);
    double InterpolateInverse(double);
    
    double SecDerivate(double x);
    
    double* GetX() const {return x;}
    double* GetY() const {return y;}
    int GetN() const {return N;}
    
    TF1* GetF0() const; //draw everything (points and interpolation function)
    TF1* GetInterpolationFunction() const;
    TF1* GetInterpolationFunctionDerivative() const ;
    TF1* GetInterpolationFunctionSecDerivative() const;
    TF1* GetInterpolationInverseFunction() const;
    
private:
    void SetCurvatureLines();
    double fInterpolator(double *fx, double *par) {
        return Interpolate(fx[0]);
    }
    double fDerivator(double *fx, double *par) {
        return Derivate(fx[0]);
    }
    double fSecDerivator(double *fx, double *par) {
        return SecDerivate(fx[0]);
    }
    double fInterpolatorInverse(double *fx, double *par) {
        return InterpolateInverse(fx[0]);
    }
    TF1* FInterpolator; //interpolation function
    TF1* D1FInterpolator = NULL; //interpolation derivative function
    TF1* D2FInterpolator = NULL;//interpolation second derivative function
    TF1* FInterpolatorInverse = NULL;
    TF1* F0 = NULL; //eventual underlying function
    double* K = NULL; //2nd derivatives
};

#endif
