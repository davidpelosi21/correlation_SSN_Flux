#ifndef	__FCmatrixFull__
#define __FCmatrixFull__

#include "FCMatrix.h"


class FCmatrixFull : public FCmatrix {
public:
    FCmatrixFull(double**fm = NULL, int lin=0 , int col=0);
    FCmatrixFull(double*, int, int);
    FCmatrixFull(vector<Vec>);
    ~FCmatrixFull();
    
    FCmatrixFull(const FCmatrix&);
    
    FCmatrixFull operator+(const FCmatrix&) const;
    FCmatrixFull operator-(const FCmatrix&) const;
    FCmatrixFull operator*(const FCmatrix&) const;
    FCmatrixFull operator*(double) const;
    FCmatrixFull operator*(const Vec&) const;
    FCmatrixFull& operator=(const FCmatrix&);
    friend FCmatrixFull operator* (double, FCmatrix&);
    
    
    Vec& operator[](int);
    Vec operator[](int) const;
    int GetBanda() const;
    Vec GetRow(int) const;
    Vec GetCol(int) const;
    FCmatrixFull Delete(int, int);
    
    int GetRowMax(int i=0, int a=0, double s = 1);
    int GetColMax(int i=0, int a=0, double s = 1);
    void GaussElimination ();
    double Determinant();
    FCmatrixFull Cofatores();
    vector <Vec> GetM() const;
    int RowSize();
    int ColSize();
    double Trace() const;
    FCmatrixFull Transposta() const;
    FCmatrixFull Inversa();
    
    int Getrowindices(int i) const {return rowindices[i];}
    int Getcolindices(int i) const {return colindices[i];}
    
    string GetClassname() const;
    void Print();
    void swapRows(int,int);
    void swapCols (int, int);
    static FCmatrixFull ReadFile(string);
    
private:
    int* rowindices=NULL;
    int* colindices=NULL;
    int banda;
};

#endif 
