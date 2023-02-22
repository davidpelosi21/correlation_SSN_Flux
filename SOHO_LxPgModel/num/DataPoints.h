#ifndef __DataPoints__
#define __DataPoints__

#include "TGraph.h"
#include "TPad.h"
#include <string>
using namespace std;


class DataPoints {
public:
	DataPoints(int = 0, double* = NULL, double* = NULL);
	DataPoints(const DataPoints&);
	virtual ~DataPoints();
	virtual void Print(string namefile = "DataPoints.txt");
    virtual double* GetX() const{return x;}
    virtual double* GetY() const {return y;}
    virtual int GetN() const {return N;}

protected:
	int N; // number of data points
	double *x, *y; // arrays
	TGraph *g=NULL;
};

#endif
