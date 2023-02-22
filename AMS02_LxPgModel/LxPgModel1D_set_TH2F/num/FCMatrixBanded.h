#ifndef __FCmatrixBanded__
#define __FCmatrixBanded__

#include "FCMatrix.h"
#include "FCMatrixFull.h"

class FCmatrixBanded : public FCmatrix{
public:
    FCmatrixBanded(double** fm=NULL,int lin=0, int col=0);
    FCmatrixBanded(vector<Vec>);
    FCmatrixBanded(const FCmatrix&);
    ~FCmatrixBanded();
    
    FCmatrixBanded operator+(const FCmatrix&) const;
    FCmatrixBanded operator-(const FCmatrix&) const;
    FCmatrixFull operator*(const FCmatrix&) const;
    FCmatrixBanded operator*(double) const;
    FCmatrixFull operator*(const Vec&) const;
    FCmatrixBanded& operator=(const FCmatrix&);
    friend FCmatrixBanded operator* (double, FCmatrixBanded&);
    FCmatrixBanded Delete(int, int);
    
    Vec& operator[](int);
    Vec operator[](int) const;
    int GetBanda() const;
    Vec GetRow(int) const;
    Vec GetCol(int) const;
    int GetRowMax(int i=0, int a=0, double s = 1);
    int GetColMax(int i=0, int a=0, double s = 1);
    void GaussElimination ();
    double Determinant();
    FCmatrixBanded Cofatores();
    vector <Vec> GetM() const;
    int RowSize();
    int ColSize();
    double Trace() const;
    FCmatrixBanded Transposta() const;
    FCmatrixBanded Inversa();
    
    int Getrowindices(int i) const{return rowindices[i];}
    int Getcolindices(int i) const{return colindices[i];}
    
    string GetClassname() const;
    void Print();
    void swapRows(int,int);
    void swapCols (int, int);
    
private:
    int* rowindices=NULL;
    int* colindices=NULL;
    int banda;
};
#endif
