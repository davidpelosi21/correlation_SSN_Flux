#ifndef __EqSolver__
#define __EqSolver__
#include "Vec.h"
#include <stdio.h>
#include "FCMatrix.h"
#include "FCMatrixFull.h"
#include "FCMatrixBanded.h"
#include <cmath>

class EqSolver {
public:
    EqSolver();
    EqSolver(const FCmatrix&, const Vec&);
    void SetConstants(const Vec&);
    void SetMatrix(const FCmatrixFull&);
    Vec GaussEliminationSolver();
    Vec LUdecompositionSolver();
    Vec RecursiveTriDiagonal(); // only for BANDED Matrices
    Vec GSIterator(double tol=1.E-4);
    Vec JIterator(double tol=1.E-4);
    void Print();
    
private:
    void LUdecomposition(FCmatrix* A){ //retorna A rescrito com L na subdiagonal inferior, com diagonal de L toda com 1's, e com U na diagonal e subdiagonal inferior
        int n = A->RowSize();
        if(A->GetClassname() == "FULL"){
            for(int i=0; i<n; ++i)
            {
                for(int j=0; j<i; ++j)
                {
                    double a = (*A)[i][j];
                    for(int p=0; p<j; ++p) a = a - (*A)[i][p]*(*A)[p][j];
                    (*A)[i][j] = a/(*A)[j][j];
                }
                for(int u=i; u<n; ++u)
                {
                    double a  = (*A)[i][u];
                    for(int s=0; s<i; ++s) a = a - (*A)[i][s]*(*A)[s][u];
                    (*A)[i][u] = a;
                }
            }
        }
        if(A->GetClassname() == "BANDED") {;}
    }
    
    
    void GaussElimination(FCmatrix* A, Vec& index){
        double**x = new double*[A->RowSize()];
        for(int i=0; i<A->RowSize(); ++i) x[i] = new double[A->RowSize()+1];
        for(int i=0; i<A->RowSize(); ++i){
            x[i][A->RowSize()] = index[i];
            for (int j=0; j<A->RowSize(); ++j) {
                x[i][j] = (*A)[i][j];
            }
        }
        FCmatrixFull B(x, A->RowSize(), A->RowSize()+1);
        B.GaussElimination();
        A->GaussElimination();
        for(int i=0; i<A->RowSize(); ++i) index[i] = B[i][A->RowSize()];
    }
    
    FCmatrix *M;
    Vec b; 
};

#endif
