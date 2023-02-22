#ifndef __FCmatrix__
#define __FCmatrix__
#include <string>
#include <vector>
#include "Vec.h"
using namespace std;

class FCmatrix {
public:
    FCmatrix(double**fm=NULL,int lin=0, int col=0); //Matrix fM x fn
    FCmatrix(double*,int, int);
    FCmatrix(vector<Vec>);
    FCmatrix(const FCmatrix&);
    virtual ~FCmatrix();
    
    virtual Vec& operator[](int) = 0;
    virtual Vec operator[](int) const = 0;
    virtual FCmatrix& operator=(const FCmatrix&)=0;
    
    virtual Vec GetRow(int i) const = 0;
    virtual Vec GetCol(int i) const = 0;
    virtual double Determinant() = 0;
    virtual void Print() = 0;
    virtual int GetRowMax(int i=0, int a=0, double s = 1) = 0;
    virtual int GetColMax(int i=0, int a=0, double s = 1) = 0;
    virtual void swapRows(int i, int j) = 0;
    virtual void swapCols(int, int)=0;
    virtual string GetClassname() const =0;
    virtual vector <Vec> GetM() const = 0;
    virtual int GetBanda() const = 0;
    virtual int RowSize() = 0;
    virtual int ColSize() = 0;
    virtual void GaussElimination() = 0;
    virtual double Trace() const = 0;
    virtual int Getrowindices(int) const = 0;
    virtual int Getcolindices(int) const = 0;
    
    
protected:
    vector <Vec> M;
    string classname;
};

#endif
