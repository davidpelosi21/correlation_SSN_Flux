#include "FCMatrixBanded.h"
#include <cmath>
#include <cstdio>


FCmatrixBanded::FCmatrixBanded(double** fm, int lin, int col): FCmatrix(fm,lin,col){
    classname="BANDED";
    
    if(fm){
        rowindices = new int[lin];
        colindices = new int[col];
        
        for(int i = 0; i<lin; i++)
            rowindices[i]=i;
        
        for(int i = 0; i<col; i++)
            colindices[i]=i;
        
        for(int i=0; fm[0][col-i-1]==0; i++)
            banda=col-i-1;
    }
}

FCmatrixBanded::FCmatrixBanded(vector<Vec> aux):  FCmatrix(aux){
    
    classname="BANDED";
    
    rowindices = new int[aux.size()];
    colindices = new int[aux[0].size()];
    
    for(int i = 0; i<aux.size(); i++)
        rowindices[i]=i;
    
    for(int i = 0; i<aux[0].size(); i++)
        colindices[i]=i;
}

FCmatrixBanded::FCmatrixBanded(const FCmatrix& mbase){
    classname = "BANDED";
    
    M.clear();
    
    for(int i=0; i<mbase.GetM().size(); i++)
        M.push_back((mbase.GetM())[i]);
    
    delete [] rowindices;
    delete [] colindices;
    
    rowindices = new int[M.size()];
    colindices = new int[M[0].size()];
    
    for(int i = 0; i<M.size(); i++)
        rowindices[i]=mbase.Getrowindices(i);
    
    for(int i = 0; i<M[0].size(); i++)
        colindices[i]=mbase.Getcolindices(i);
}

FCmatrixBanded::~FCmatrixBanded(){
    if(rowindices && colindices){
        delete[] rowindices;
        delete[] colindices;
    }
}

Vec& FCmatrixBanded::operator[](int i){
    return M[rowindices[i]];
}


Vec FCmatrixBanded::operator[](int i) const{
    return M[rowindices[i]];
}

int FCmatrixBanded::GetBanda() const{
    return banda;
}

FCmatrixBanded FCmatrixBanded::operator+(const FCmatrix& m1) const{
    
    if(M[0].size() != m1[0].size() || M.size() != m1.GetM().size()){
        printf("Operação Abortada. Matrizes Somadas possuem Diferentes Tamanhos.\n");
        return FCmatrixBanded(M);
    }
    
    vector<Vec> A;
    
    for(int i=0; i<M.size(); ++i){
        A.push_back(Vec(M[rowindices[0]].size(), 0.));
        A[i][i]=M[i][i]+m1[i][i];
        
        
        for(int j=1; (j<=banda || j<=m1.GetBanda()) && (i-j)>=0 ; j++){
            A[i][i-j]=M[i][i-j]+m1[i][i-j];
        }
        
        
        for(int j=1; (j<=banda || j<=m1.GetBanda()) && (i+j)<=M.size()-1; j++){
            A[i][i+j]=M[i][i+j]+m1[i][i+j];
        }
    }
    
    return FCmatrixBanded(A);
}

FCmatrixBanded FCmatrixBanded::operator-(const FCmatrix& m1) const{
    if(M[0].size() != m1[0].size() || M.size() != m1.GetM().size()){
        printf("Operação Abortada. Matrizes Somadas Subtraídas possuem Diferentes Tamanhos.\n");
        return FCmatrixBanded(M);
    }
    
    vector<Vec> A;
    
    for(int i=0; i<M.size(); ++i){
        A.push_back(Vec(M[rowindices[0]].size(), 0.));
        A[i][i]=M[i][i]-m1[i][i];
        
        
        for(int j=1; (j<=banda || j<=m1.GetBanda()) && (i-j)>=0 ; j++){
            A[i][i-j]=M[i][i-j]-m1[i][i-j];
        }
        
        
        for(int j=1; (j<=banda || j<=m1.GetBanda()) && (i+j)<=M.size()-1; j++){
            A[i][i+j]=M[i][i+j]-m1[i][i+j];
        }
    }
    
    return FCmatrixBanded(A);
}

FCmatrixFull FCmatrixBanded::operator*(const FCmatrix& m1) const{
    
    if(M[0].size() != m1.GetM().size()){
        printf("Operação Abortada. Numero de Linhas da Primeira Matriz diferente do Numero de Linhas da Segunda Matriz.\n");
        return FCmatrixFull(M);
    }
    
    vector<Vec> A;
    
    for(int i=0; i<M.size() ; ++i){
        A.push_back(Vec(m1[0].size(),0.));
        for(int j=0; j<m1[0].size(); ++j){
            A[i][j] = M[rowindices[i]].dot(m1.GetCol(j));
        }
    }
    return FCmatrixFull(A);
}

FCmatrixBanded FCmatrixBanded::operator*(double l) const{
    vector<Vec> A;
    
    
    for(int i=0; i<M.size(); ++i){
        A.push_back(Vec(M[rowindices[0]].size(), 0.));
        A[i][i]=M[i][i]*l;
        
        
        for(int j=1; j<=banda && (i-j)>=0 ; j++){
            A[i][i-j]=M[i][i-j]*l;
        }
        
        
        for(int j=1; j<=banda && (i+j)<=M.size()-1; j++){
            A[i][i+j]=M[i][i+j]*l;
        }
    }
    
    return FCmatrixBanded(A);
}

FCmatrixFull FCmatrixBanded::operator*(const Vec& v) const {
    
    if(M[0].size() != v.size()){
        printf("Operação Abortada. Numero de Linhas da Primeira Matriz diferente do Numero de Entradas do Vetor\n");
        return FCmatrixFull(M);
    }
    
    
    vector<Vec> A;
    
    
    for(int i=0; i<M.size(); ++i){
        A.push_back(Vec(1,0.));
        A[i][0] = M[i].dot(v);
    }
    
    return FCmatrixFull(A);
}

FCmatrixBanded& FCmatrixBanded::operator=(const FCmatrix& mtx){
    
    if(this!=&mtx){
        delete[] rowindices;
        delete[] colindices;
        
        classname="BANDED";
        
        M.clear();
        
        for(int i=0; i<mtx.GetCol(0).size() ;i++)
            M.push_back(mtx[i]);
        
        
        rowindices = new int[M.size()];
        colindices = new int[M[0].size()];
        
        
        for(int i = 0; i<M.size(); i++)
            rowindices[i]=i;
        
        for(int i = 0; i<M[0].size(); i++)
            colindices[i]=i;
    }
    
    return *this;
}


FCmatrixBanded operator* (double m,FCmatrixBanded& mtx){
    return (FCmatrixBanded(mtx*m));
}


Vec FCmatrixBanded::GetRow(int i) const{
    return M[rowindices[i]];
}

Vec FCmatrixBanded::GetCol(int i) const{
    double column[M.size()];
    
    for(int j = 0; j<M.size(); j++)
        column[j]=M[rowindices[j]][colindices[i]];
    
    return Vec(M.size(),column);
}

int FCmatrixBanded::GetRowMax(int i, int a, double s){
    double max=M[i][a];
    double pos=a;
    
    for(int j=a; j<M.size(); j++){
        if(fabs(max) < fabs(M[rowindices[i]][colindices[j]])){
            max=M[rowindices[i]][colindices[j]];
            pos=j;
        }
    }
    
    return (s*pos);
}


int FCmatrixBanded::GetColMax(int i, int a,double s){
    double max=M[a][i];
    double pos=a;
    
    for(int j=a; j<M.size(); j++){
        if(fabs(max) < fabs(M[rowindices[j]][colindices[i]])){
            max=M[rowindices[j]][colindices[i]];
            pos = j;
        }
    }
    
    return (s*pos);
}

void FCmatrixBanded::GaussElimination(){
    
    for(int i=0; i<M.size() && i<M[0].size(); i++){
        if(i!=this->GetColMax(i,i)){
            this->swapRows(i,this->GetColMax(i,i));
        }
    }
    
    for(int j=0; j<M[0].size() && j<M.size(); j++){
        for(int i=j+1; i<M.size(); i++){
            M[rowindices[i]]=M[rowindices[i]]-(M[rowindices[j]]*(M[rowindices[i]][colindices[j]]/(M[rowindices[j]][colindices[j]])));
        }
    }
}

double FCmatrixBanded::Determinant(){
    FCmatrixBanded aux(M);
    
    aux.GaussElimination();
    
    double determinant=1;
    
    for(int i=0; i<aux.M.size();i++){
        determinant*=aux[i][i];
    }
    
    return determinant;
}

void FCmatrixBanded::swapRows(int i,int j){
    for(int k=0; k<M[0].size(); ++k){
        double aux = M[k][i];
        M[k][i] = M[k][j];
        M[k][j] = aux;
    }
}

void FCmatrixBanded::swapCols(int i,int j)
{
    int aux = colindices[i];
    colindices[i] = colindices[j];
    colindices[j] = aux;
}

string FCmatrixBanded::GetClassname() const{
    return "BANDED";
}

void FCmatrixBanded::Print(){
    int lin = M.size();
    int col = M[0].size();
    
    printf("\n");
    
    for(int i=0; i<lin; i++){
        for(int j=0; j<col; j++){
            printf(" %8.3lf ", M[rowindices[i]][colindices[j]]);
        }
        printf("\n");
    }
    printf("\n");
}

vector <Vec> FCmatrixBanded::GetM() const{
    return M;
}

int FCmatrixBanded::RowSize(){
    return M.size();
}

int FCmatrixBanded::ColSize(){
    return M[0].size();
}

double FCmatrixBanded::Trace() const{
    double trace=0;
    
    for(int j=0; j<M.size(); j++)
        trace+=M[j][j];
    
    return trace;
}

FCmatrixBanded FCmatrixBanded::Delete(int linha, int coluna){
    
    int lin = M.size();
    int col = M[0].size();
    int in1=0;
    int in2=0;
    double** xaux = new double*[lin-1];
    for(int m=0; m<lin-1; ++m)
        xaux[m] = new double[col-1];
    
    for(int k=0; k<lin; ++k)
    {
        if(k==linha) {++k;in1=1;}
        if(k>=lin) break;
        for(int l=0; l<col; ++l)
        {
            if(l==coluna) {++l;in2=1; if(l>=col)break;}
            if(in1==0 && in2==0)
                xaux[k][l] = M[rowindices[k]][colindices[l]];
            if(in1==1 && in2==0)
                xaux[k-1][l] = M[rowindices[k]][colindices[l]];
            if(in1==0 && in2==1)
                xaux[k][l-1] = M[rowindices[k]][colindices[l]];
            if(in1==1 && in2==1)
                xaux[k-1][l-1] = M[rowindices[k]][colindices[l]];
        }
        in2=0;
    }
    
    FCmatrixBanded A(xaux, lin-1, col-1);
    
    for(int m=0; m<lin-1; ++m)
        delete[] xaux[m];
    
    delete[] xaux;
    return A;
}

FCmatrixBanded FCmatrixBanded::Cofatores(){
    double** x = new double*[M.size()];
    
    for(int m=0; m<M.size(); ++m)
        x[m] = new double[M[0].size()];
    
    for(int i=0; i<M.size(); ++i){
        for(int j=0; j<M[0].size(); ++j){
            FCmatrixBanded A = Delete(i,j);
            double d = A.Determinant();
            if(d==0)
                x[i][j] = 0;
            else
                x[i][j] = pow(-1, i+j)*d;
        }
    }
    
    FCmatrixBanded B(x, M.size(), M[0].size());
    
    for(int m=0; m<M.size(); ++m)
        delete[] x[m];
    
    delete[] x;
    
    return B;
}

FCmatrixBanded FCmatrixBanded::Transposta() const {
    double** x = new double*[M[0].size()];
    
    for(int m=0; m<M[0].size(); ++m)
        x[m] = new double[M.size()];
    
    for(int i=0; i<M[0].size(); ++i){
        for(int j=0; j<M.size(); ++j){
            x[i][j] = M[rowindices[j]][colindices[i]];
        }
    }
    
    FCmatrixBanded B(x, M.size(), M[0].size());
    
    for(int m=0; m<M[0].size(); ++m)
        delete[] x[m];
    
    delete[] x;
    return B;
}

FCmatrixBanded FCmatrixBanded::Inversa(){
    if(M.size()!=M[0].size()){
        printf("Operação Abortada.Matriz Não Quadrada\n");
        return FCmatrixFull(M);
    }
    if(Determinant()!=0)
        return (Cofatores().Transposta())*(1/Determinant());
    else{
        printf("Operação Abortada.Matriz é singular\n");
        return FCmatrixBanded(M);
    }
}
