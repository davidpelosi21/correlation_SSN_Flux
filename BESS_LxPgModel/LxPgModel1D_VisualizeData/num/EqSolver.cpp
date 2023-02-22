#include "EqSolver.h"

EqSolver::EqSolver(){
    double** x = new double*[2];
    for(int i=0; i<2; ++i)
    {
        x[i]  =new double[2];
    }
    x[0][0] = 1;
    x[0][1] = 2;
    x[1][0] = 2;
    x[1][1] = 1;
    FCmatrixFull A(x, 2, 2);
    Vec c(2, 1);
    *M = A;
    b=c;
}

EqSolver::EqSolver(const FCmatrix& A, const Vec& v)
{
    if(A.GetClassname() == "FULL"){
        M  = new FCmatrixFull(A.GetM());
        b = v;
    }
    if(A.GetClassname() == "BANDED"){
        M = new FCmatrixBanded(A.GetM());
        b = v;
    }
}

void EqSolver::SetConstants(const Vec& v){
    b = v;
}

void EqSolver::SetMatrix(const FCmatrixFull& A){
    *M = A;
}

Vec EqSolver::LUdecompositionSolver(){
    this->LUdecomposition(M);
    int n = M->RowSize();
    FCmatrixFull L;
    FCmatrixFull U;
    L = *M;
    U = *M;
    for(int i=0; i<n; ++i) {
        for(int j=0; j<n; ++j){
            if(i>j) {U[i][j] = 0;}
            if(i<j) {L[i][j] = 0;}
            if(i==j) {L[i][j] = 1;}
        }
    }
    
    Vec x(n, 0.);
    Vec y(n, 0.);
    for (int k=0; k<n; k++) {
        double sumC = 0;
        for (int i=0; i<k; i++) {
            sumC += y[i]*L[k][i];
        }
        y[k] = b[k] - sumC;
    }
    for (int k=n-1; k>=0; k--) {
        double sumC = 0;
        for (int i=k+1; i<n; i++) {
            sumC += x[i]*U[k][i];}
        x[k] = (y[k] - sumC)/U[k][k];
    }
    return x;
    
}

Vec EqSolver::GaussEliminationSolver(){
    this->GaussElimination(M, b);
    Vec sol(M->RowSize(),0.);
    for(int i=M->RowSize()-1; i>=0; --i){
        double sum=0;
        for(int j=i+1; j<M->RowSize(); j++) sum+=(*M)[i][j]*sol[j];
        sol[i] = (b[i] - sum)/(*M)[i][i];
    }
    
    return sol;
}

Vec EqSolver::GSIterator(double eps){
    double  m = M->RowSize();
    Vec x(m,0.);
    Vec x_aux(m,0.);
    double aux;
    bool btol = false;
    int it = 0.;
    while (!btol && (it++ < 1000)) {
        x_aux = x;
        for (int i=0; i<m; i++) {
            aux = 0.;
            for (int j=0; j<m; j++){
                if (i!= j){
                    aux -= (*M)[i][j]*x[j];
                }
            }
            x[i] = (1/((*M)[i][i]))*(aux+b[i]);
            if (fabs(x[i]-x_aux[i]) < eps)
                btol = true;
            else
                btol = false;
        }
    }
    return x;
    
}
Vec EqSolver::JIterator(double eps){
    double  m = M->RowSize();
    Vec x(m,0.);
    Vec x_aux(m,0.);
    double aux;
    bool btol = false;
    int it = 0.;
    while (!btol && (it++ < 1000)) {
        x_aux = x;
        for (int i=0; i<m; i++) {
            aux = 0.;
            for (int j=0; j<m; j++){
                if (i!= j){
                    aux -= (*M)[i][j]*x_aux[j];
                }
            }
            x[i] = (1/((*M)[i][i]))*(aux+b[i]);
            if (fabs(x[i]-x_aux[i]) < eps)
                btol = true;
            else
                btol = false;
        }
    }
    return x;
    
}

void EqSolver::Print(){
    int lin = b.size();
    //int col = M->[0].size();
    printf("\n");
    
    for(int i=0; i<lin; i++){
        printf("|");
        for(int j=0; j<lin; j++){
            printf("%8.3lf", (*M)[i][j]);
        }
        printf("  |%8.3lf|", b[i]);
        printf("\n");
    }
    printf("\n");
}

Vec EqSolver::RecursiveTriDiagonal(){
    //n -> row index , starting at 0
    int n = M->RowSize()-1;
    double* alpha = new double[n];//Lower diagonal
    double* beta = new double[n+1]; //main diagonal
    double* gama = new double[n]; //upper diagonal

    for (int i=0; i<n; ++i) alpha[i] = (*M)[i+1][i];
    for (int i=0; i<n+1; ++i) beta[i] = (*M)[i][i];
    for (int i=0; i<n; ++i) gama[i] = (*M)[i][i+1];

    //introducing variables g and h
    double* g = new double[n];
    double* h = new double[n];
    
    g[n-1] = -(alpha[n-1]/beta[n]);
    h[n-1] = b[n]/beta[n];
    
    //backwards loop to obtain all values of parameter g[n-1] and h[n-1]
    for(int i=n-1; i>0; --i)  g[i-1] = -(alpha[i-1]/(beta[i]+gama[i]*g[i]));
    for(int i=n-1; i>0; --i)  h[i-1] = (b[i]-gama[i]*h[i])/(beta[i]+gama[i]*g[i]);

    //recursive solution  
    Vec solve(n+1,0.);
    solve[0] = (b[0]-gama[0]*h[0])/(beta[0]+gama[0]*g[0]);
    for (int i=1;i<n+1; ++i) solve[i] = g[i-1]*solve[i-1] + h[i-1];
           
    delete[] alpha;
    delete[] beta;
    delete[] gama;
    delete[] g;
    delete[] h;
    return solve;
}
